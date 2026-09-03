"""Codec adapters and parameter sweeps -- single source of truth.

These adapters were duplicated verbatim in hpx_codec_bench.ipynb,
utm_codec_bench.ipynb and _adapters_extract.py; they now live here so the two
grids can never drift apart (a silent difference between the HEALPix-side and
UTM-side encoder settings would directly bias every ratio in the paper).

Every adapter returns the encoded byte count actually produced and a decoded
array, so RMSE is always measured against a real round trip.
"""

from __future__ import annotations

import io
import time
from dataclasses import dataclass, field

import numpy as np


@dataclass
class CodecRunResult:
    codec: str
    config: dict
    total_bytes: int
    encode_seconds: float
    decode_seconds: float
    reconstruction: np.ndarray = field(repr=False)


def rmse(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    mask = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2))) if mask.any() else np.nan


# --- float codecs ------------------------------------------------------------
def encode_zfp(tile, *, tolerance):
    import zfpy
    tile = np.asarray(tile, dtype=np.float32)
    t0 = time.perf_counter()
    enc = zfpy.compress_numpy(tile, tolerance=float(tolerance))
    t1 = time.perf_counter()
    dec = zfpy.decompress_numpy(enc)
    t2 = time.perf_counter()
    return CodecRunResult("zfp", dict(tolerance=float(tolerance)), len(enc),
                          t1 - t0, t2 - t1, dec)


def encode_sz3(tile, *, abs_error):
    from pysz import sz as _sz, szConfig, szErrorBoundMode
    tile = np.ascontiguousarray(np.asarray(tile, dtype=np.float32))
    config = szConfig()
    config.errorBoundMode = szErrorBoundMode.ABS
    config.absErrorBound = float(abs_error)
    t0 = time.perf_counter()
    compressed, _ = _sz.compress(tile, config)
    t1 = time.perf_counter()
    dec, _ = _sz.decompress(compressed, np.float32, tile.shape)
    t2 = time.perf_counter()
    total_bytes = int(getattr(compressed, "nbytes", None) or len(compressed))
    return CodecRunResult("sz3", dict(abs_error=float(abs_error)), total_bytes,
                          t1 - t0, t2 - t1, np.asarray(dec, dtype=np.float32))


# --- integer image codecs ----------------------------------------------------
def _to_u16(tile):
    tile = np.asarray(tile, dtype=np.float32)
    t_min, t_max = float(tile.min()), float(tile.max())
    span = (t_max - t_min) or 1.0
    return ((tile - t_min) / span * 65535).astype(np.uint16), t_min, span


def _from_u16(u16, t_min, span):
    return u16.astype(np.float32) / 65535.0 * span + t_min


def encode_jpeg2000_lossy(tile, *, quality):
    from PIL import Image
    u16, t_min, span = _to_u16(tile)
    factor = max(1, int(quality))     # never round to 0 (past bug)
    buf = io.BytesIO()
    t0 = time.perf_counter()
    Image.fromarray(u16).save(buf, format="jpeg2000", quality_mode="rates",
                              quality_layers=[factor])
    t1 = time.perf_counter()
    total_bytes = buf.tell()
    buf.seek(0)
    dec_u16 = np.array(Image.open(buf), dtype=np.uint16)
    t2 = time.perf_counter()
    return CodecRunResult("jpeg2000_lossy", dict(quality=factor), total_bytes,
                          t1 - t0, t2 - t1, _from_u16(dec_u16, t_min, span))


def encode_jpeg2000_lossless(tile, **_ignored):
    import imagecodecs
    u16, t_min, span = _to_u16(tile)
    t0 = time.perf_counter()
    enc = imagecodecs.jpeg2k_encode(u16, level=0)
    t1 = time.perf_counter()
    dec_u16 = imagecodecs.jpeg2k_decode(enc)
    t2 = time.perf_counter()
    return CodecRunResult("jpeg2000_lossless", {}, len(enc), t1 - t0, t2 - t1,
                          _from_u16(dec_u16, t_min, span))


_jpegxl_warned = False


def encode_jpegxl_lossy(tile, *, level=None, distance=None):
    """JPEG XL lossy.

    `level` is imagecodecs' shared 0-100 quality convention; `distance` is
    JPEG XL's own Butteraugli target, where smaller is better and 0 is
    lossless. The level ladder saturates around RMSE 1e-3 on this imagery, so
    the tight end of the sweep is reached with `distance` instead -- without
    it JPEG XL is simply *absent* from the high-fidelity targets and would win
    or lose those by forfeit rather than on merit.
    """
    global _jpegxl_warned
    import imagecodecs
    if (level is None) == (distance is None):
        raise ValueError("pass exactly one of level= or distance=")
    u16, t_min, span = _to_u16(tile)
    t0 = time.perf_counter()
    try:
        enc = (imagecodecs.jpegxl_encode(u16, level=int(level)) if level is not None
               else imagecodecs.jpegxl_encode(u16, distance=float(distance)))
    except TypeError as e:
        if not _jpegxl_warned:
            print(f"[jpegxl] encode signature mismatch ({e}) -- check your "
                  "imagecodecs version's jpegxl_encode kwargs. Skipping jpegxl.")
            _jpegxl_warned = True
        raise
    t1 = time.perf_counter()
    dec_u16 = imagecodecs.jpegxl_decode(enc)
    t2 = time.perf_counter()
    conf = dict(level=int(level)) if level is not None else dict(distance=float(distance))
    return CodecRunResult("jpegxl_lossy", conf, len(enc),
                          t1 - t0, t2 - t1, _from_u16(dec_u16, t_min, span))


def encode_jpegxl_lossless(tile, **_ignored):
    import imagecodecs
    u16, t_min, span = _to_u16(tile)
    t0 = time.perf_counter()
    enc = imagecodecs.jpegxl_encode(u16, lossless=True)
    t1 = time.perf_counter()
    dec_u16 = imagecodecs.jpegxl_decode(enc)
    t2 = time.perf_counter()
    return CodecRunResult("jpegxl_lossless", {}, len(enc), t1 - t0, t2 - t1,
                          _from_u16(dec_u16, t_min, span))


def encode_dwt_zstd(tile, *, wavelet="db4", levels=4, quantize_bits):
    import pywt
    import zstandard as zstd
    tile = np.asarray(tile, dtype=np.float32)
    t0 = time.perf_counter()
    coeffs = pywt.wavedec2(tile, wavelet=wavelet, level=levels, mode="periodization")
    arr, slices = pywt.coeffs_to_array(coeffs)
    max_level = (1 << (quantize_bits - 1)) - 1
    scale = float(np.max(np.abs(arr))) or 1.0
    dtype_q = np.int8 if quantize_bits <= 8 else np.int16
    q = np.round(np.clip(arr / scale * max_level, -max_level, max_level)).astype(dtype_q)
    payload = zstd.ZstdCompressor(level=19).compress(q.tobytes())
    t1 = time.perf_counter()
    total_bytes = len(payload) + 4                      # + float32 scale sidecar
    q2 = np.frombuffer(zstd.ZstdDecompressor().decompress(payload),
                       dtype=dtype_q).reshape(arr.shape)
    arr2 = (q2.astype(np.float32) / max_level * scale).astype(np.float64)
    coeffs2 = pywt.array_to_coeffs(arr2, slices, output_format="wavedec2")
    dec = pywt.waverec2(coeffs2, wavelet=wavelet, mode="periodization").astype(np.float32)
    t2 = time.perf_counter()
    return CodecRunResult("dwt_zstd",
                          dict(wavelet=wavelet, levels=levels,
                               quantize_bits=quantize_bits),
                          total_bytes, t1 - t0, t2 - t1, dec)


def encode_uncompressed(tile, **_ignored):
    tile = np.asarray(tile, dtype=np.float32)
    return CodecRunResult("uncompressed", {}, tile.nbytes, 0.0, 0.0, tile.copy())


# --- the package's own zarr serializer path ----------------------------------
def encode_hcc_zarr(flat_cells, *, name="hcc-dct", delta=2, top_k=64):
    """Round-trip a flat NESTED chunk through hcc's zarr codec.

    This is the software artifact the paper describes (standard codecs made
    usable as zarr serializers on HEALPix chunks); it is reported separately
    from the five established codecs and never mixed into their comparison.
    """
    import zarr
    import hcc  # noqa: F401  (registers the entry points)
    from hcc.codecs.dct import HccDctCodec
    from hcc.codecs.grad_dct import HccGradDctCodec
    from hcc.codecs.grad_dwt import HccGradDwtCodec

    serializers = {"hcc-dct": HccDctCodec, "hcc-grad-dct": HccGradDctCodec,
                   "hcc-grad-dwt": HccGradDwtCodec}
    if name not in serializers:
        raise KeyError(f"unknown hcc codec {name!r} (have {sorted(serializers)})")
    data = np.asarray(flat_cells, dtype=np.float32).ravel()

    store = zarr.storage.MemoryStore()
    root = zarr.api.synchronous.create_group(store, zarr_format=3)
    t0 = time.perf_counter()
    z = root.create_array("cells", data=data, fill_value=float("nan"),
                          chunks=(data.size,),
                          serializer=serializers[name](delta=delta, top_k=top_k),
                          filters=None, compressors=None,
                          dimension_names=["cells"])
    t1 = time.perf_counter()
    dec = np.asarray(z[:], dtype=np.float32)
    t2 = time.perf_counter()
    total_bytes = sum(len(v) for k, v in store._store_dict.items()
                      if not str(k).endswith("zarr.json"))
    return CodecRunResult(name.replace("-", "_"),
                          dict(delta=delta, top_k=top_k), int(total_bytes),
                          t1 - t0, t2 - t1, dec)


# --- sweeps ------------------------------------------------------------------
ZFP_SWEEP = [dict(tolerance=float(v)) for v in np.logspace(-5, -1, 8)]
SZ3_SWEEP = [dict(abs_error=float(v)) for v in np.logspace(-5, -1, 8)]
# jpeg2000 "quality" is a compression-RATIO target: smaller = higher quality.
JPEG2000_LOSSY_SWEEP = [dict(quality=q) for q in (1, 2, 3, 5, 8, 12, 20, 35, 60, 100)]
JPEGXL_LOSSY_SWEEP = (
    [dict(level=lv) for lv in (10, 25, 40, 55, 68, 78, 86, 92, 96, 99)]
    # distance points extend the sweep past the level ladder's RMSE floor
    # (~1e-3 on Sentinel-2 B04), so JPEG XL is actually evaluable at the
    # high-fidelity targets instead of dropping out of them.
    + [dict(distance=d) for d in (2.0, 1.0, 0.5, 0.25, 0.1, 0.05, 0.01)]
)
DWT_ZSTD_SWEEP = [dict(quantize_bits=b) for b in (4, 6, 8, 10, 12, 14, 16)]

CODEC_RUNS = {
    "zfp": (encode_zfp, ZFP_SWEEP),
    "sz3": (encode_sz3, SZ3_SWEEP),
    "jpeg2000_lossy": (encode_jpeg2000_lossy, JPEG2000_LOSSY_SWEEP),
    "jpeg2000_lossless": (encode_jpeg2000_lossless, [dict()]),
    "jpegxl_lossy": (encode_jpegxl_lossy, JPEGXL_LOSSY_SWEEP),
    "jpegxl_lossless": (encode_jpegxl_lossless, [dict()]),
    "dwt_zstd": (encode_dwt_zstd, DWT_ZSTD_SWEEP),
    "uncompressed": (encode_uncompressed, [dict()]),
}

HCC_SWEEP = [dict(name=n, delta=d, top_k=k)
             for n in ("hcc-dct", "hcc-grad-dwt")
             for d in (1, 2, 3)
             for k in (16, 32, 64, 128, 256)]


def match_shape(tile, run):
    """Force a decoded array back to the input shape, or reject the point.

    Transform codecs pad: a multi-level wavelet decomposition of an odd-sized
    image (the 637x637 equal-area UTM crop) reconstructs one or more pixels
    larger, which is why this never shows up on the 1024x1024 HEALPix side.
    Cropping the reconstruction is the standard fix, and it is the honest one
    here: the padding was really encoded, so its bytes stay charged to the
    codec; only the comparison is trimmed.
    """
    want = np.shape(tile)
    got = np.shape(run.reconstruction)
    if got == want:
        return run
    if len(got) != len(want) or any(g < w for g, w in zip(got, want)):
        raise ValueError(f"{run.codec} returned {got}, smaller than the input "
                         f"{want}; cannot compare")
    run.reconstruction = run.reconstruction[tuple(slice(0, w) for w in want)]
    run.config = dict(run.config, decoded_shape=list(got), cropped_to=list(want))
    return run


def iter_runs(tile, codecs=None, on_error="warn", stats=None):
    """Yield CodecRunResult for every (codec, config) in the sweep.

    `stats`, if given, is filled with per-codec attempted/ok/skipped counts and
    the reason of every skip -- a codec silently dropping out of the sweep is
    how a comparison ends up decided by which codec was evaluable rather than
    by which one was better.
    """
    codecs = list(CODEC_RUNS) if codecs is None else list(codecs)
    for codec in codecs:
        if codec not in CODEC_RUNS:
            raise KeyError(f"unknown codec {codec!r} (have {sorted(CODEC_RUNS)})")
        fn, sweep = CODEC_RUNS[codec]
        for conf in sweep:
            if stats is not None:
                st = stats.setdefault(codec, {"attempted": 0, "ok": 0, "skipped": []})
                st["attempted"] += 1
            try:
                run = match_shape(tile, fn(tile, **conf))
            except Exception as exc:                      # noqa: BLE001
                if on_error == "raise":
                    raise
                msg = f"{type(exc).__name__}: {exc}"
                print(f"    [skip] {codec} {conf}: {msg}")
                if stats is not None:
                    stats[codec]["skipped"].append({"config": dict(conf),
                                                     "error": msg})
                continue
            if stats is not None:
                stats[codec]["ok"] += 1
            yield run
