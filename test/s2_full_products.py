"""Whole Sentinel-2 L2A products, EOPF Zarr -> HEALPix Zarr: bytes on disk and quality.

For every tile of the configuration (one product per tile, see ``s2_ten_tiles.json``):

1. ``mirror``: the EOPF product is copied object by object (metadata and chunk
   objects exactly as served) into a local store, so its size per array is
   the size ESA distributes;
2. ``convert``: ``healpix_convert.create_healpix_dataset`` on the whole product
   with the legacy-converters L2A preset (all groups: reflectance, masks,
   classification, atmosphere, probabilities, detector footprints, geometry,
   meteorology); the reflectance groups are stored with
   ``hcc-sentinel2-pyramid`` (Q and O from the bands' ``scale_factor`` and
   ``add_offset``), the other groups with the preset's codecs. Every
   codestream written is decoded again and compared with the cells it was
   made from: the codec error per array, in digital numbers;
3. ``sizes``: bytes per array of both stores, and the UTM reflectance bands
   under lossless JPEG 2000 and JPEG XL per source chunk (the SAFE reference);
4. ``quality``: on a random sample of interior chunks of every reflectance
   resolution, the native pixels, the PSF transfer of those pixels without
   codec and the stored cells, both read back onto the native pixels with the
   transfer's own operator: transfer RMSE, RMSE added by the codec, total;
5. ``upload``: the results (not the mirror) to ``dest/<site>/``.

Every step leaves a file in ``work_dir/<site>/`` and is skipped when that file
exists, so a tile resumes at its first unfinished step. Tiles run in parallel
processes (``--parallel``), each with ``--threads`` torch threads.

Usage
-----
    python s2_full_products.py s2_ten_tiles.json --parallel 4 --threads 4
    python s2_full_products.py s2_ten_tiles.json --site paris_urban --extent 2.33 48.84 2.37 48.87   # quick check
    python s2_full_products.py s2_ten_tiles.json --summary                                      # the table
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eopf_s2_to_healpix as conv  # noqa: E402  (EopfSource.chunk_keys, lossless_sizes)

REFLECTANCE = {"r10m": (20, 10.0), "r20m": (19, 20.0), "r60m": (17, 60.0)}
CHUNK_LEVEL = 11
META_V2 = (".zgroup", ".zattrs", ".zarray", ".zmetadata")


def log(site_dir: Path, msg: str) -> None:
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with open(site_dir / "progress.log", "a") as f:
        f.write(line + "\n")


def dir_bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*")
               if p.is_file() and p.name not in META_V2 + ("zarr.json",))


# --- 1. the EOPF product, object by object ------------------------------------------
def mirror(href: str, dest: Path, storage_options: dict | None = None) -> dict:
    import fsspec
    import zarr

    href = href.rstrip("/")
    fs, _ = fsspec.core.url_to_fs(href, **(storage_options or {}))
    root = zarr.open_group(href, mode="r", storage_options=storage_options or None)
    v3 = root.metadata.zarr_format == 3
    meta_names = ("zarr.json",) if v3 else (".zgroup", ".zattrs")
    keys = list(meta_names) + ([] if v3 else [".zmetadata"])
    arrays = {}
    for path, node in root.members(max_depth=None):
        if isinstance(node, zarr.Array):
            keys += [f"{path}/zarr.json"] if v3 else [f"{path}/.zarray", f"{path}/.zattrs"]
            ck = [f"{path}/{k}" for k in conv.EopfSource.chunk_keys(node)]
            arrays[path] = {"shape": list(node.shape), "chunks": list(node.chunks), "dtype": str(node.dtype),
                            "keys": ck}
            keys += ck
        else:
            keys += [f"{path}/{m}" for m in meta_names]
    if dest.exists():
        shutil.rmtree(dest)
    got = {}
    for i in range(0, len(keys), 32):
        batch = keys[i:i + 32]
        data = fs.cat([f"{href}/{k}" for k in batch], on_error="omit")
        for k in batch:
            blob = data.get(f"{href}/{k}")
            if blob is None or isinstance(blob, Exception):
                continue
            f = dest / k
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_bytes(blob)
            got[k] = len(blob)
    out = {}
    for path, a in arrays.items():
        out[path] = {k: v for k, v in a.items() if k != "keys"}
        out[path]["stored_bytes"] = int(sum(got.get(k, 0) for k in a["keys"]))
        out[path]["chunk_objects"] = int(sum(k in got for k in a["keys"]))
    meta_bytes = sum(v for k, v in got.items() if k.rsplit("/", 1)[-1] in META_V2 + ("zarr.json",))
    return {"href": href, "zarr_format": 3 if v3 else 2, "arrays": out, "metadata_bytes": int(meta_bytes),
            "total_bytes": int(sum(got.values()))}


def rechunk_copy(local: Path, dest: Path, max_chunk: int = 2048) -> dict:
    """A copy of the mirror for the converter, with the single-chunk 10980 x 10980
    arrays (masks, atmosphere, footprints) cut into ``max_chunk`` blocks: every
    level-11 chunk of the conversion reads its window of those arrays, and with
    one chunk per array that window costs the decompression of the whole array.
    Same codec, same values; the sizes reported are the mirror's, not this copy's."""
    import numcodecs

    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(local, dest)
    zmeta_path = dest / ".zmetadata"
    zmeta = json.loads(zmeta_path.read_text())
    changed = []
    for zarray in sorted(dest.rglob(".zarray")):
        meta = json.loads(zarray.read_text())
        shape, chunks = meta["shape"], meta["chunks"]
        if len(shape) != 2 or all(c <= max_chunk for c in chunks) or meta.get("filters"):
            continue
        adir = zarray.parent
        rel = adir.relative_to(dest).as_posix()
        sep = meta.get("dimension_separator", ".")
        comp = numcodecs.get_codec(meta["compressor"]) if meta.get("compressor") else None
        dtype = np.dtype(meta["dtype"])
        fill = meta.get("fill_value")
        grid = [math.ceil(s / c) for s, c in zip(shape, chunks)]
        data = np.full(shape, 0 if fill is None else fill, dtype=dtype)
        for idx in np.ndindex(*grid):
            f = adir / sep.join(map(str, idx))
            if not f.exists():
                continue
            raw = f.read_bytes()
            block = np.frombuffer(comp.decode(raw) if comp else raw, dtype=dtype).reshape(chunks)
            sl = tuple(slice(i * c, min((i + 1) * c, s)) for i, c, s in zip(idx, chunks, shape))
            data[sl] = block[tuple(slice(0, x.stop - x.start) for x in sl)]
            f.unlink()
        new = [min(c, max_chunk) for c in chunks]
        for idx in np.ndindex(*[math.ceil(s / c) for s, c in zip(shape, new)]):
            block = np.full(new, 0 if fill is None else fill, dtype=dtype)
            sl = tuple(slice(i * c, min((i + 1) * c, s)) for i, c, s in zip(idx, new, shape))
            block[tuple(slice(0, x.stop - x.start) for x in sl)] = data[sl]
            block = np.ascontiguousarray(block)
            (adir / sep.join(map(str, idx))).write_bytes(bytes(comp.encode(block)) if comp else block.tobytes())
        meta["chunks"] = new
        zarray.write_text(json.dumps(meta, indent=4))
        zmeta["metadata"][f"{rel}/.zarray"] = meta
        changed.append(rel)
    zmeta_path.write_text(json.dumps(zmeta, indent=4))
    return {"rechunked": changed, "max_chunk": max_chunk}


# --- 2. the conversion ---------------------------------------------------------------
def make_settings(local: Path, cfg: dict) -> dict:
    """The legacy-converters L2A preset, reflectance groups under the pyramid codec."""
    import zarr
    from legacy_converters.settings import sentinel2 as presets

    settings = copy.deepcopy(presets.SENTINEL2_L2A_CONVERT_SETTINGS)
    root = zarr.open_group(str(local), mode="r")
    for res in REFLECTANCE:
        group = f"measurements/reflectance/{res}"
        g = root[group]
        band = next(n for n in sorted(g.array_keys()) if n not in ("x", "y"))
        attrs = g[band].attrs.asdict()
        sf, ao = float(attrs.get("scale_factor", 1e-4)), float(attrs.get("add_offset", -0.1))
        codec = copy.deepcopy(cfg["codec"])
        conf = codec.setdefault("configuration", {})
        conf.setdefault("quantification_value", 1.0 / sf)
        conf.setdefault("add_offset", ao / sf)
        settings["group_settings"][group]["codecs"] = [codec]
    # The other quantised fields come out of the converter decoded (float64 for aot,
    # wvp and scl, float32 for the probabilities) and the preset stores them under
    # zstd; put them back on their integer lattice with the lossless serializer
    # (digital numbers under JPEG XL): DN = round(Q v - O), Q = 1 / scale_factor.
    # Where 0 is a valid value (classes, probabilities) the lattice is shifted by
    # one, since DN 0 is the serializer's no-data.
    lossless = {"quality/atmosphere": 0, "conditions/mask/l2a_classification": 1, "quality/probability": 1}
    for group, gs in settings["group_settings"].items():
        head = group.rsplit("/", 1)[0]
        if head not in lossless or group not in root:
            continue
        g = root[group]
        band = next(n for n in sorted(g.array_keys()) if n not in ("x", "y"))
        attrs = g[band].attrs.asdict()
        sf, ao = float(attrs.get("scale_factor", 1.0)), float(attrs.get("add_offset", 0.0))
        gs["codecs"] = [{"name": "hcc-sentinel2",
                         "configuration": {"quantification_value": 1.0 / sf, "add_offset": ao / sf - lossless[head]}}]
    return settings


class CodecErrors:
    """Decode every pyramid codestream written and compare it with its cells."""

    def __init__(self):
        self.current = None
        self.stats = {}

    def install(self):
        import zarr
        from hcc.compression import sentinel2_pyramid as sp

        original_set, original_encode = zarr.Array.__setitem__, sp.encode
        me = self

        def setitem(arr, key, value):
            me.current = arr.path
            try:
                return original_set(arr, key, value)
            finally:
                me.current = None

        def encode(cells, **kwargs):
            stream = original_encode(cells, **kwargs)
            q = float(kwargs.get("quantification_value", 1e4))
            cells = np.asarray(cells, dtype=np.float64)
            dec = sp.decode(stream, quantification_value=q, add_offset=float(kwargs.get("add_offset", -1000.0)),
                            dtype="float64")
            ok = np.isfinite(cells) & np.isfinite(dec)
            err = (dec[ok] - cells[ok]) * q
            info = sp.stream_info(stream)
            s = me.stats.setdefault(me.current or "?", dict(chunks=0, cells=0, valid=0, nan_mismatch=0,
                                                             sum_sq=0.0, max_abs=0.0, sum_step=0.0,
                                                             sum_sq_signal=0.0, bytes=0))
            s["chunks"] += 1
            s["cells"] += int(cells.size)
            s["valid"] += int(ok.sum())
            s["nan_mismatch"] += int((np.isfinite(cells) != np.isfinite(dec)).sum())
            s["sum_sq"] += float((err ** 2).sum())
            s["max_abs"] = max(s["max_abs"], float(np.abs(err).max()) if err.size else 0.0)
            s["sum_step"] += float(info.get("step", 0.0)) if isinstance(info, dict) else 0.0
            v = cells[ok] * q
            s["sum_sq_signal"] += float(((v - v.mean()) ** 2).sum()) if v.size else 0.0
            s["bytes"] += len(stream)
            return stream

        zarr.Array.__setitem__ = setitem
        sp.encode = encode

    def summary(self):
        out = {}
        for path, s in self.stats.items():
            n = max(1, s["valid"])
            out[path] = dict(chunks=s["chunks"], valid_cells=s["valid"], nan_mismatch=s["nan_mismatch"],
                             rmse_dn=math.sqrt(s["sum_sq"] / n), max_abs_dn=s["max_abs"],
                             mean_step_dn=s["sum_step"] / max(1, s["chunks"]),
                             within_chunk_std_dn=math.sqrt(s["sum_sq_signal"] / n), bytes=s["bytes"])
        return out


def _string_coordinates_as_zarr_strings() -> None:
    """Workaround: the string coordinates of ``conditions/geometry`` (``angle``,
    ``band``, ``detector``) reach ``zarr.create`` as numpy ``object`` arrays,
    which zarr 3 refuses as ambiguous; create them as variable-length strings."""
    import healpix_convert.healpix_converters as hc
    import numcodecs.zarr3  # noqa: F401  -- registers numcodecs.delta with zarr

    original = hc.zarr.create
    if getattr(original, "_strings_patched", False):
        return

    def create(*args, **kwargs):
        if "dtype" in kwargs and np.dtype(kwargs["dtype"]) == object:
            kwargs["dtype"] = str
            kwargs["fill_value"] = ""
        # The cell ids of a dense chunk are one contiguous range, which zstd alone keeps
        # at about 1 byte per cell (as much as a quantised band): delta-code them first.
        if str(kwargs.get("path", "")).endswith("/cell_ids") and kwargs.get("codecs") is None:
            kwargs["codecs"] = [{"name": "numcodecs.delta", "configuration": {"dtype": "<u8"}},
                                {"name": "bytes", "configuration": {"endian": "little"}},
                                {"name": "zstd", "configuration": {"level": 3, "checksum": False}}]
        return original(*args, **kwargs)

    create._strings_patched = True
    hc.zarr.create = create

    # Workaround: on float32 coordinates (CAMS and ECMWF latitude/longitude) the nearest
    # neighbour search of healpix_resample 2026.6.0 mixes float32 and float64 tensors
    # ("Index put requires the source and destination dtypes match"); give it float64.
    import healpix_resample.knn as knn
    import torch

    original_knn = knn.healpix_weighted_nearest

    def as64(v):
        if isinstance(v, np.ndarray) and v.dtype == np.float32:
            return v.astype(np.float64)
        if isinstance(v, torch.Tensor) and v.dtype == torch.float32:
            return v.double()
        return v

    def nearest(*args, **kwargs):
        try:
            return original_knn(*args, **kwargs)
        except RuntimeError as e:
            if "dtypes match" not in str(e):
                raise
            return original_knn(*[as64(a) for a in args], **{k: as64(v) for k, v in kwargs.items()})

    knn.healpix_weighted_nearest = nearest


SKIPPED_CHUNKS: dict[str, int] = {}


def _skip_chunks_without_support() -> None:
    """Workaround: on a chunk at the edge of the tile with only a handful of input
    points, the nearest-neighbour resampler of healpix_resample 2026.6.0 raises
    "No HEALPix cell passed the threshold" and stops the whole conversion. Such a
    chunk has no cell to fill: leave it at the fill value and count it."""
    import healpix_convert.healpix_converters as hc

    for cls in (hc.DenseChunkConverter, getattr(hc, "SparseChunkConverter", None)):
        if cls is None or getattr(cls.convert, "_skip_patched", False):
            continue
        original = cls.convert

        def convert(self, *args, _original=original, **kwargs):
            try:
                return _original(self, *args, **kwargs)
            except RuntimeError as e:
                if "No HEALPix cell passed the threshold" not in str(e):
                    raise
                key = str(self.output_path)
                SKIPPED_CHUNKS[key] = SKIPPED_CHUNKS.get(key, 0) + 1

        convert._skip_patched = True
        cls.convert = convert


def _share_neighbour_search(cache_dir: Path) -> None:
    """The masks, atmosphere, probability, footprint and classification groups of a
    resolution sit on the same UTM grid and go to the same HEALPix cells with the
    same resampler: the neighbour search of a chunk (most of the time of a
    nearest-neighbour group) is the same for all of them. Keep each search on disk,
    keyed by the exact input coordinates and the search parameters, and reuse it.
    Only small searches (``N * Npt <= 2e6``, i.e. not the 9-point PSF ones, which
    are made once per chunk anyway) are kept."""
    import hashlib

    import healpix_resample.knn as knn
    import torch

    cache_dir.mkdir(parents=True, exist_ok=True)
    original = knn.healpix_weighted_nearest

    def search(longitude1, latitude1, level, Npt, **kwargs):
        n = int(longitude1.numel()) if isinstance(longitude1, torch.Tensor) else len(longitude1)
        # the cache lives for one conversion: after it (quality step) search directly
        if n * Npt > 2_000_000 or not cache_dir.exists():
            return original(longitude1, latitude1, level, Npt, **kwargs)
        h = hashlib.blake2b(digest_size=20)
        for v in (longitude1, latitude1, kwargs.get("out_cell_ids")):
            if v is not None:
                a = v.cpu().numpy() if isinstance(v, torch.Tensor) else np.asarray(v)
                h.update(str(a.dtype).encode())
                h.update(np.ascontiguousarray(a).tobytes())
        h.update(repr((level, Npt, sorted((k, repr(v)) for k, v in kwargs.items() if k != "out_cell_ids"))).encode())
        f = cache_dir / f"{h.hexdigest()}.pt"
        if f.exists():
            search.hits += 1
            return torch.load(f)
        out = original(longitude1, latitude1, level, Npt, **kwargs)
        torch.save(out, f)
        search.misses += 1
        return out

    search.hits = search.misses = 0
    knn.healpix_weighted_nearest = search
    return search


def _sync(out: Path, remote: str | None, storage_options: dict, since: float) -> int:
    """Upload the files of the output store written since ``since``."""
    if remote is None:
        return 0
    import fsspec

    fs, root = fsspec.core.url_to_fs(remote.rstrip("/"), **storage_options)
    n = 0
    for f in out.rglob("*"):
        if f.is_file() and f.stat().st_mtime >= since:
            fs.put_file(str(f), f"{root}/{f.relative_to(out).as_posix()}")
            n += 1
    return n


def convert(local: Path, out: Path, cfg: dict, extent=None, state_path: Path | None = None,
            remote: str | None = None, batch: int = 64) -> dict:
    """The conversion, group by group and ``batch`` chunks at a time.

    With ``state_path`` (outside the container's /tmp) and ``remote`` (the store's
    copy in the bucket), every finished batch is uploaded and recorded: after a
    restart of the server the store is downloaded again and the conversion goes
    on from the next batch, with the codec statistics gathered so far."""
    import math as _math

    import hcc  # noqa: F401  -- registers the codecs with zarr
    import shapely
    import zarr
    from healpix_convert import convert_group_to_healpix, create_staging_cache, prepare_healpix_dataset
    from healpix_convert.healpix_converters import init_converter
    from healpix_convert.settings.common import validate_convert_settings

    so = cfg.get("storage_options", {})
    errors = CodecErrors()
    errors.install()
    _string_coordinates_as_zarr_strings()
    _skip_chunks_without_support()
    cache_dir = out.parent / f"{out.name}.knn_cache"
    search = _share_neighbour_search(cache_dir)
    settings = validate_convert_settings(make_settings(local, cfg))
    cache = create_staging_cache([str(local)], settings,
                                 output_extent=None if extent is None else shapely.box(*extent))
    state = json.loads(state_path.read_text()) if state_path is not None and state_path.exists() else {}

    def save():
        if state_path is not None:
            state["codec_stats"] = errors.stats
            state["skipped_chunks"] = dict(SKIPPED_CHUNKS)
            state_path.write_text(json.dumps(state))

    if state.get("prepared") and not out.exists():
        import fsspec

        # Only the metadata comes back (a few hundred small files): the chunks already
        # written stay in the bucket, the conversion only adds new ones, and the sizes
        # and the quality are measured on the bucket's copy.
        fs, root = fsspec.core.url_to_fs(remote.rstrip("/"), **so)
        meta = [p for p in fs.find(root) if p.rsplit("/", 1)[-1] == "zarr.json"]
        fs.get(meta, [str(out / Path(p).relative_to(root)) for p in meta])
        state["resumed"] = True
        print(f"resumed: {len(meta)} metadata files from {remote}", flush=True)
    if not state.get("prepared"):
        if out.exists():
            shutil.rmtree(out)
        since = time.time()
        prepare_healpix_dataset(cache, settings, str(out))
        _sync(out, remote, so, since)
        state.update(prepared=True, groups={}, seconds=0.0)
        save()
    errors.stats = state.get("codec_stats", {})
    SKIPPED_CHUNKS.update(state.get("skipped_chunks", {}))
    # Groups on the same UTM grid converted with the same nearest-neighbour resampler
    # (masks, atmosphere, probabilities, footprints, classification of one resolution)
    # go through the chunks together, batch by batch: the neighbour search of a chunk
    # is reused across them while it is fresh, and the cache only ever holds one batch.
    families: dict[tuple, list] = {}
    for path in sorted(cache.input_spatial_groups, key=str):
        g = str(path)
        converter = init_converter(path, cache=cache, settings=settings,
                                   output_store=zarr.open_group(str(out), mode="a").store)
        no_chunk = converter.settings.chunk.method == "no_chunk"
        n = 1 if no_chunk else _math.ceil(converter.cell_dim_size / converter.cell_dim_chunk_size)
        res = g.rsplit("/", 1)[-1]
        shared = converter.settings.resampler.name == "nearest" and res in REFLECTANCE and not no_chunk
        key = ("nearest", res, n) if shared else (g,)
        families.setdefault(key, []).append((g, n, no_chunk))
    for members in families.values():
        n = members[0][1]
        start = min(state["groups"].get(g, 0) for g, _, _ in members)
        while start < n:
            stop = n if members[0][2] else min(start + batch, n)
            for g, _, no_chunk in members:
                first = max(start, state["groups"].get(g, 0))
                if first >= stop:
                    continue
                since, t0 = time.time(), time.perf_counter()
                convert_group_to_healpix(f"/{g}", cache=cache, settings=settings, output_path=str(out),
                                         chunk_indices=None if no_chunk else list(range(first, stop)))
                state["seconds"] = state.get("seconds", 0.0) + time.perf_counter() - t0
                _sync(out, remote, so, since)
                state["groups"][g] = stop
                save()
            print(f"  {', '.join(g for g, _, _ in members)}: {stop}/{n} chunks, "
                  f"{state['seconds'] / 60:.0f} min of conversion", flush=True)
            for f in cache_dir.glob("*.pt"):
                f.unlink()
            start = stop
    since = time.time()
    zarr.consolidate_metadata(str(out))
    _sync(out, remote, so, since)
    shutil.rmtree(cache_dir, ignore_errors=True)
    return {"seconds": round(state.get("seconds", 0.0), 1), "codec_errors": errors.summary(),
            "neighbour_searches": {"computed": search.misses, "reused": search.hits},
            "skipped_chunks_without_support": dict(SKIPPED_CHUNKS),
            "resumed": bool(state.get("resumed"))}


# --- 3. the sizes --------------------------------------------------------------------
def healpix_name(utm_group: str) -> str:
    head, _, res = utm_group.rpartition("/")
    levels = {"r10m": "20", "r20m": "19", "r60m": "17"}
    return f"{head}/{levels[res]}" if res in levels else utm_group


def _store_files(out: Path, remote: str | None, storage_options: dict) -> dict[str, int]:
    """Size of every file of the HEALPix store, relative path -> bytes: the bucket's
    copy when there is one (the local store may hold only the chunks written since
    the last restart), else the local store."""
    if remote is None:
        return {f.relative_to(out).as_posix(): f.stat().st_size for f in out.rglob("*") if f.is_file()}
    import fsspec

    fs, root = fsspec.core.url_to_fs(remote.rstrip("/"), **storage_options)
    return {p[len(root) + 1:]: int(info["size"]) for p, info in fs.find(root, detail=True).items()}


def _open_store(out: Path, remote: str | None, storage_options: dict):
    import zarr

    if remote is None:
        return zarr.open_group(str(out), mode="r")
    return zarr.open_group(remote, mode="r", storage_options=storage_options)


def sizes(local: Path, out: Path, src: dict, remote: str | None = None, storage_options: dict | None = None) -> dict:
    import zarr

    import hcc  # noqa: F401

    so = storage_options or {}
    files = _store_files(out, remote, so)
    meta_names = META_V2 + ("zarr.json",)

    def stored(prefix: str) -> int:
        return sum(v for k, v in files.items() if k.startswith(prefix + "/") and k.rsplit("/", 1)[-1] not in meta_names)

    # the store that holds it all (the bucket's copy when there is one)
    root = data_root = _open_store(out, remote, so)
    hpx = {}
    for path, node in root.members(max_depth=None):
        if isinstance(node, zarr.Array):
            hpx[path] = {"shape": list(node.shape), "chunks": list(node.chunks), "dtype": str(node.dtype),
                         "codecs": [c.get("name") for c in node.metadata.to_dict().get("codecs", [])],
                         "stored_bytes": stored(path)}
    root = data_root
    # what the arrays whose codec this study changed would weigh under the preset's
    # default (bytes + zstd): estimated on a random sample of chunks
    import numcodecs

    zstd = numcodecs.Zstd(level=0)
    rng = np.random.default_rng(0)
    changed = ("quality/atmosphere/", "conditions/mask/l2a_classification/", "quality/probability/")
    for path, a in hpx.items():
        if not (path.endswith("/cell_ids") or path.startswith(changed)) or len(a["shape"]) != 1:
            continue
        arr = root[path]
        n_chunks = math.ceil(a["shape"][0] / a["chunks"][0])
        sample = rng.choice(n_chunks, size=min(48, n_chunks), replace=False)
        c = a["chunks"][0]
        nb = [len(zstd.encode(np.ascontiguousarray(arr[i * c:(i + 1) * c]))) for i in sample]
        a["preset_default_bytes_estimate"] = int(np.mean(nb) * n_chunks)
    by_group = {}
    for path, a in src["arrays"].items():
        g = path.rsplit("/", 1)[0]
        by_group.setdefault(g, {"utm_bytes": 0, "healpix_bytes": 0})["utm_bytes"] += a["stored_bytes"]
    for path, a in hpx.items():
        g = path.rsplit("/", 1)[0]
        utm_g = next((u for u in by_group if healpix_name(u) == g), g)
        row = by_group.setdefault(utm_g, {"utm_bytes": 0, "healpix_bytes": 0})
        row["healpix_bytes"] += a["stored_bytes"]
        row["healpix_preset_default_bytes"] = (row.get("healpix_preset_default_bytes", 0)
                                               + a.get("preset_default_bytes_estimate", a["stored_bytes"]))
    # the SAFE reference: the UTM reflectance under lossless JPEG 2000 (and JPEG XL), per source chunk
    lossless = {}
    src_root = zarr.open_group(str(local), mode="r")
    for res in REFLECTANCE:
        g = f"measurements/reflectance/{res}"
        for band in sorted(src_root[g].array_keys()):
            if band in ("x", "y"):
                continue
            arr = src_root[g][band]
            dn = arr[:]
            lossless[f"{g}/{band}"] = dict(conv.lossless_sizes(dn, tuple(arr.chunks)),
                                           valid_pixels=int((dn != 0).sum()))
    return {"utm_total_bytes": src["total_bytes"],
            "healpix_total_bytes": sum(v for k, v in files.items() if k.rsplit("/", 1)[-1] not in meta_names),
            "healpix_metadata_bytes": sum(v for k, v in files.items() if k.rsplit("/", 1)[-1] in meta_names),
            "groups": by_group, "healpix_arrays": hpx, "utm_reflectance_lossless": lossless}


# --- 4. the quality ------------------------------------------------------------------
def _rss_gb() -> float:
    return int(open("/proc/self/statm").read().split()[1]) * os.sysconf("SC_PAGE_SIZE") / 1e9


def _phase(msg: str) -> None:
    print(time.strftime("%H:%M:%S ") + f"[quality] {msg} (rss {_rss_gb():.1f} GB)", flush=True)


def _rss_guard(limit_gb: float) -> None:
    """Stop this process before it takes the server down: the quality step went from 2
    to 48 GB in 20 s once, and the out-of-memory kill restarted the whole container.
    (RLIMIT_AS cannot be used: the CUDA build of torch reserves a huge address space.)"""
    import threading

    def watch():
        while True:
            rss = _rss_gb()
            if rss > limit_gb:
                print(f"[quality] RSS GUARD: {rss:.1f} GB > {limit_gb} GB, abort", flush=True)
                os._exit(3)
            time.sleep(0.02)

    threading.Thread(target=watch, daemon=True).start()



def quality(local: Path, out: Path, cfg: dict, n_chunks: int, seed: int = 0, remote: str | None = None) -> dict:
    import healpix_geo.nested as hn
    import healpix_resample as hr
    import pyproj
    import shapely
    import torch
    import zarr

    import hcc  # noqa: F401

    src_root = zarr.open_group(str(local), mode="r")
    props = src_root.attrs["stac_discovery"]["properties"]
    crs = pyproj.CRS.from_user_input(props.get("proj:code") or f"EPSG:{props['proj:epsg']}")
    to_utm = pyproj.Transformer.from_crs(4326, crs, always_xy=True)
    to_ll = pyproj.Transformer.from_crs(crs, 4326, always_xy=True)
    preset = make_settings(local, cfg)["group_settings"]
    hpx_root = _open_store(out, remote, cfg.get("storage_options", {}))
    result = {}
    for res, (level, pixel) in REFLECTANCE.items():
        group = f"measurements/reflectance/{res}"
        gs = preset[group]
        q = gs["codecs"][0]["configuration"]["quantification_value"]
        o = gs["codecs"][0]["configuration"]["add_offset"]
        src = src_root[group]
        hpx = hpx_root[healpix_name(group)]
        bands = [b for b in sorted(src.array_keys()) if b not in ("x", "y")]
        n = 4 ** (level - CHUNK_LEVEL)
        x, y = np.asarray(src["x"][:], dtype=float), np.asarray(src["y"][:], dtype=float)
        _phase(f"{res}: reading the chunk cell ids")
        cell_ids = np.asarray(hpx["cell_ids"][::n], dtype=np.int64)
        _phase(f"{res}: {cell_ids.size} chunks")
        chunk_ids = cell_ids >> (2 * (level - CHUNK_LEVEL))
        clon, clat = hn.healpix_to_lonlat(chunk_ids.astype(np.uint64), CHUNK_LEVEL, ellipsoid="WGS84")
        cx, cy = to_utm.transform(np.asarray(clon), np.asarray(clat))
        margin = 3000.0
        interior = np.nonzero((cx > x.min() + margin) & (cx < x.max() - margin)
                              & (cy > y.min() + margin) & (cy < y.max() - margin))[0]
        order = np.random.default_rng(seed).permutation(interior)
        rows, used = [], 0
        for i in order:
            if used >= n_chunks:
                break
            cid = int(chunk_ids[i])
            cells = cid * n + np.arange(n, dtype=np.int64)
            vlon, vlat = hn.vertices(np.array([cid], dtype=np.uint64), CHUNK_LEVEL, ellipsoid="WGS84")
            poly = shapely.Polygon(list(zip(*to_utm.transform(np.asarray(vlon[0]), np.asarray(vlat[0])))))
            outer, inner = shapely.buffer(poly, 6 * pixel), shapely.buffer(poly, -6 * pixel)
            xmin, ymin, xmax, ymax = shapely.bounds(outer)
            ci, ri = np.nonzero((x >= xmin) & (x <= xmax))[0], np.nonzero((y >= ymin) & (y <= ymax))[0]
            c0, c1, r0, r1 = ci.min(), ci.max() + 1, ri.min(), ri.max() + 1
            xx, yy = np.meshgrid(x[c0:c1], y[r0:r1])
            inside = shapely.contains_xy(outer, xx, yy)
            dn = {b: np.asarray(src[b][r0:r1, c0:c1], dtype=np.float64) for b in bands}
            if np.mean(dn[bands[0]][inside] > 0) < 0.99:
                continue  # chunks with no-data pixels: the read-back is not defined everywhere
            lon, lat = to_ll.transform(xx[inside], yy[inside])
            op = hr.PSFResampler(lon_deg=np.asarray(lon), lat_deg=np.asarray(lat), level=level, out_cell_ids=cells,
                                 verbose=False, nest=True, ellipsoid="WGS84", **gs["resampler"]["init_params"])
            op_ids = np.asarray(op.get_cell_ids() if hasattr(op, "get_cell_ids") else op.cell_ids, dtype=np.int64)
            pos = op_ids - int(cells[0])
            # read-back onto the pixels: v @ MT, computed as sparse @ dense with MT transposed
            mtt = op.MT.to_sparse_coo().t().coalesce()                    # (N, K)
            _phase(f"{res} chunk {i}: operator K={op_ids.size} N={len(lon)}")
            norm = torch.sparse.mm(mtt, torch.ones(op_ids.size, 1, dtype=op.MT.dtype)).squeeze(1)
            good = (norm > 0.5).numpy() & shapely.contains_xy(inner, xx[inside], yy[inside])
            # DN 0 is the L2A no-data (e.g. lakes the atmospheric correction took below
            # rho = -0.1): the converter reads it as NaN through the CF decoding and the
            # transfer leaves it out, so the reference transfer and the scores must too
            valid = np.stack([dn[b][inside] > 0 for b in bands])                # (B, P)
            native = np.stack([np.where(dn[b][inside] > 0, (dn[b][inside] + o) / q, np.nan) for b in bands])
            params = gs["resampler"]["resample_params"]
            # one band at a time, with its own parameters and its own no-data, as the converter does
            raw = np.stack([
                np.asarray(op.resample(native[k].astype(np.float32),
                                       **{p: (v[b] if isinstance(v, dict) else v) for p, v in params.items()}).cell_data,
                           dtype=np.float64).reshape(-1)
                for k, b in enumerate(bands)])
            _phase(f"{res} chunk {i}: transfer done")
            stored = np.stack([np.asarray(hpx[b][i * n:(i + 1) * n], dtype=np.float64)[pos] for b in bands])
            _phase(f"{res} chunk {i}: stored cells read")

            def back(v):
                # normalised by the weight of the finite cells only, so that a cell left
                # NaN (no valid pixel under it) does not read back as a zero reflectance
                t = torch.as_tensor(np.nan_to_num(v), dtype=op.MT.dtype)
                w = torch.as_tensor(np.isfinite(v), dtype=op.MT.dtype)
                den = torch.sparse.mm(mtt, w.T).T
                return (torch.sparse.mm(mtt, t.T).T / torch.where(den > 0, den, norm)).numpy()

            nat, raw_b, sto_b = native[:, good] * q, back(raw)[:, good] * q, back(stored)[:, good] * q
            for k, b in enumerate(bands):
                v = valid[k, good]
                nk, rk, sk = nat[k][v], raw_b[k][v], sto_b[k][v]
                rows.append(dict(chunk=int(i), chunk_id=cid, band=b, pixels=int(v.sum()),
                                 nodata_pixels=int((~v).sum()),
                                 native_std_dn=float(nk.std()),
                                 transfer_rmse_dn=float(np.sqrt(np.mean((rk - nk) ** 2))),
                                 added_rmse_dn=float(np.sqrt(np.mean((sk - rk) ** 2))),
                                 total_rmse_dn=float(np.sqrt(np.mean((sk - nk) ** 2)))))
            used += 1
        per_band = {}
        for b in bands:
            r = [row for row in rows if row["band"] == b]
            if not r:
                continue
            w = np.array([row["pixels"] for row in r], dtype=float)

            def pooled(key):
                return float(np.sqrt(np.sum(w * np.array([row[key] for row in r]) ** 2) / w.sum()))

            per_band[b] = dict(chunks=len(r), native_std_dn=pooled("native_std_dn"),
                               transfer_rmse_dn=pooled("transfer_rmse_dn"), added_rmse_dn=pooled("added_rmse_dn"),
                               total_rmse_dn=pooled("total_rmse_dn"))
        result[group] = {"level": level, "pixel_m": pixel, "chunks_sampled": used, "per_band": per_band,
                         "per_chunk": rows}
    return result


# --- 5. upload ------------------------------------------------------------------------
def upload(site_dir: Path, dest: str, storage_options: dict) -> int:
    import fsspec

    fs, root = fsspec.core.url_to_fs(dest.rstrip("/"), **storage_options)
    n = 0
    for f in sorted(site_dir.rglob("*")):
        parts = f.relative_to(site_dir).parts
        # the mirror stays out; the HEALPix store is already there (uploaded batch by batch)
        if f.is_file() and not {"eopf.zarr", "eopf_rechunked.zarr", "healpix.zarr"} & set(parts):
            fs.put_file(str(f), f"{root}/{f.relative_to(site_dir).as_posix()}")
            n += 1
    return n


# --- one tile -------------------------------------------------------------------------
def run_site(cfg: dict, site: str, threads: int, extent=None, quality_now: bool = False) -> None:
    import torch

    torch.set_num_threads(threads)
    tile = cfg["tiles"][site]
    if quality_now and _remote(cfg, site, "quality.json") is not None:
        print(f"{site}: quality already measured", flush=True)
        return
    if quality_now:
        _rss_guard(cfg.get("quality_rss_limit_gb", 24.0))
    # The hub image ships healpix-resample v2026.06.0. There, one NaN input point makes
    # the PSF result of a whole band NaN (the chunk is then not written at all), and the
    # fallback for weakly supported out_cell_ids can grow past the container memory
    # (GRID4EARTH/g4e-desp-argo#28). Both are fixed in v2026.07.0: use the copy
    # installed next to the run when there is one.
    hr_path = cfg.get("healpix_resample_path")
    if hr_path and Path(hr_path).expanduser().exists():
        sys.path.insert(0, str(Path(hr_path).expanduser()))
        sys.modules.pop("healpix_resample", None)
        import healpix_resample

        print(f"{site}: healpix_resample from {healpix_resample.__file__}", flush=True)
    if not quality_now and extent is None and cfg.get("dest") and _remote(cfg, site, "DONE") is not None:
        print(f"{site}: already in {cfg['dest']}", flush=True)
        return
    work = Path(cfg["work_dir"]).expanduser() / site
    work.mkdir(parents=True, exist_ok=True)
    local, out = work / "eopf.zarr", work / "healpix.zarr"
    steps_path = work / "steps.json"
    steps = json.loads(steps_path.read_text()) if steps_path.exists() else {}
    # The results of the finished steps are also kept in the bucket, so that a
    # restart of the server (which empties /tmp) does not redo them.
    kept = ("mirror", "convert", "sizes", "quality")
    remote_steps = (json.loads(_remote(cfg, site, "results_steps.json") or b"{}")
                    if extent is None and cfg.get("dest") else {})
    for name in kept:
        if name in remote_steps and name not in steps:
            (work / f"{name}.json").write_bytes(_remote(cfg, site, f"{name}.json"))
            steps[name] = remote_steps[name]
    if not local.exists():
        steps.pop("rechunk", None)
        if not {"sizes", "quality"} <= set(steps):
            steps.pop("mirror", None)  # the next steps read the product again

    def step(name, fn):
        if name in steps:
            return
        log(work, f"{site}: {name}...")
        t0 = time.perf_counter()
        res = fn()
        if res is not None:
            (work / f"{name}.json").write_text(json.dumps(res, indent=1))
        steps[name] = round(time.perf_counter() - t0, 1)
        steps_path.write_text(json.dumps(steps, indent=1))
        if name in kept and extent is None and cfg.get("dest"):
            import fsspec

            fs, root = fsspec.core.url_to_fs(cfg["dest"].rstrip("/"), **cfg.get("storage_options", {}))
            if res is not None:
                fs.pipe_file(f"{root}/{site}/{name}.json", json.dumps(res, indent=1).encode())
            fs.pipe_file(f"{root}/{site}/results_steps.json",
                         json.dumps({k: v for k, v in steps.items() if k in kept}).encode())
        log(work, f"{site}: {name} done in {steps[name] / 60:.1f} min")

    step("mirror", lambda: dict(mirror(tile["href"], local, tile.get("storage_options")), product=tile["product"]))
    conv_input = work / "eopf_rechunked.zarr"
    if "convert" not in steps:
        step("rechunk", lambda: rechunk_copy(local, conv_input))
    if extent is not None:
        # quick check of the whole chain on a small box, outside the tile's steps
        check = work.parent / f"{site}_extent"
        check.mkdir(exist_ok=True)
        res = convert(conv_input, check / "healpix.zarr", cfg, extent)
        (check / "convert.json").write_text(json.dumps(res, indent=1))
        log(work, f"{site}: extent check converted in {res['seconds']} s")
        return
    state_dir = Path(cfg.get("state_dir", "~/hcc-test/s2_ten_state")).expanduser()
    state_dir.mkdir(parents=True, exist_ok=True)
    remote = f"{cfg['dest'].rstrip('/')}/{site}/healpix.zarr" if cfg.get("dest") else None
    step("convert", lambda: convert(conv_input, out, cfg, extent, state_path=state_dir / f"{site}.json",
                                    remote=remote, batch=cfg.get("batch_chunks", 64)))
    if "convert" in steps and conv_input.exists():
        shutil.rmtree(conv_input)
    src = json.loads((work / "mirror.json").read_text())
    step("sizes", lambda: sizes(local, out, src, remote, cfg.get("storage_options", {})))
    if cfg.get("defer_quality") and not quality_now:
        # measured afterwards, one tile at a time (``--quality``): the step has taken
        # the whole server down several times while five conversions were running
        log(work, f"{site}: quality deferred")
    else:
        step("quality", lambda: quality(local, out, cfg, cfg.get("quality_chunks", 12), remote=remote))
        if quality_now:
            return
    if cfg.get("dest"):
        step("upload", lambda: {"files": upload(work, f"{cfg['dest'].rstrip('/')}/{site}",
                                                cfg.get("storage_options", {}))})
        # the marker last: the work directory (on the container's /tmp) does not survive
        # a restart of the server, the bucket does
        import fsspec

        fs, root = fsspec.core.url_to_fs(cfg["dest"].rstrip("/"), **cfg.get("storage_options", {}))
        fs.pipe_file(f"{root}/{site}/steps.json", json.dumps(steps, indent=1).encode())
        fs.pipe_file(f"{root}/{site}/DONE", b"")
        log(work, f"{site}: in {cfg['dest']}/{site}")
    if cfg.get("delete_mirror_after_upload") and "upload" in steps and local.exists():
        shutil.rmtree(work)


# --- the table ------------------------------------------------------------------------
def _remote(cfg: dict, site: str, name: str) -> bytes | None:
    """A file of a finished tile in the bucket, or None."""
    import fsspec

    fs, root = fsspec.core.url_to_fs(cfg["dest"].rstrip("/"), **cfg.get("storage_options", {}))
    try:
        return fs.cat_file(f"{root}/{site}/{name}")
    except FileNotFoundError:
        return None


def summary(cfg: dict) -> dict:
    work = Path(cfg["work_dir"]).expanduser()
    rows = {}
    for site, tile in cfg["tiles"].items():
        if cfg.get("dest") and _remote(cfg, site, "quality.json") is not None:
            s, c, qd = (json.loads(_remote(cfg, site, f"{n}.json")) for n in ("sizes", "convert", "quality"))
        elif (work / site / "quality.json").exists():
            s, c, qd = (json.loads((work / site / f"{n}.json").read_text()) for n in ("sizes", "convert", "quality"))
        else:
            continue
        refl_utm = sum(v["utm_bytes"] for g, v in s["groups"].items() if g.startswith("measurements/reflectance"))
        refl_hpx = sum(v["healpix_bytes"] for g, v in s["groups"].items() if g.startswith("measurements/reflectance"))
        jp2 = sum(v["jpeg2000_lossless_bytes"] for v in s["utm_reflectance_lossless"].values())
        rows[site] = dict(
            product=tile["product"], utm_total_mb=s["utm_total_bytes"] / 1e6, healpix_total_mb=s["healpix_total_bytes"] / 1e6,
            reflectance_utm_mb=refl_utm / 1e6, reflectance_jp2_mb=jp2 / 1e6, reflectance_healpix_mb=refl_hpx / 1e6,
            reflectance_x_jp2=refl_hpx / jp2, total_x_utm=s["healpix_total_bytes"] / s["utm_total_bytes"],
            convert_min=c["seconds"] / 60,
            quality={g: v["per_band"] for g, v in qd.items()},
            codec_errors=c["codec_errors"])
    (work / "summary.json").write_text(json.dumps(rows, indent=1))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", type=Path)
    ap.add_argument("--site", help="run one tile in this process")
    ap.add_argument("--parallel", type=int, default=1, help="tiles converted at the same time")
    ap.add_argument("--threads", type=int, default=None, help="torch threads per tile")
    ap.add_argument("--extent", type=float, nargs=4, metavar=("W", "S", "E", "N"),
                    help="with --site: convert only this lon/lat box (quick check, separate directory)")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--quality", action="store_true", help="with --site: measure the deferred quality step")
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    if args.summary:
        print(json.dumps(summary(cfg), indent=1))
        return 0
    threads = args.threads or max(1, len(os.sched_getaffinity(0)) // args.parallel)
    if args.site:
        run_site(cfg, args.site, threads, args.extent, quality_now=args.quality)
        return 0
    work = Path(cfg["work_dir"]).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    pending, running, failed = [s for s in cfg["tiles"] if not (cfg.get("dest") and _remote(cfg, s, "DONE") is not None)], {}, []
    print("tiles to do:", pending, flush=True)
    while pending or running:
        while pending and len(running) < args.parallel:
            site = pending.pop(0)
            (work / site).mkdir(exist_ok=True)
            logf = open(work / site / "run.log", "a")
            running[site] = subprocess.Popen([sys.executable, "-u", str(Path(__file__).resolve()), str(args.config),
                                              "--site", site, "--threads", str(threads)],
                                             stdout=logf, stderr=subprocess.STDOUT)
            print(time.strftime("%H:%M:%S"), "started", site, flush=True)
        time.sleep(30)
        for site, p in list(running.items()):
            if p.poll() is not None:
                print(time.strftime("%H:%M:%S"), site, "exit", p.returncode, flush=True)
                if p.returncode:
                    failed.append(site)
                del running[site]
    rows = summary(cfg)
    if cfg.get("dest"):
        import fsspec

        fs, root = fsspec.core.url_to_fs(cfg["dest"].rstrip("/"), **cfg.get("storage_options", {}))
        fs.pipe_file(f"{root}/summary.json", json.dumps(rows, indent=1).encode())
    print("failed:", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
