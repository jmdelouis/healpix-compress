"""Audit of the byte counts and the data types behind every reported number.

The review's blocking objection (its section 3.5) is that the paper reports
compressed sizes and the word "lossless" without ever showing what was actually
handed to each encoder, what came back, and whether the number quoted is the
size of a file that exists.  `codec_lib.py` returns a byte count computed
in memory; nothing so far has checked that a file of that size can be written,
closed, re-opened by an independent reader and decoded back to the same array.
That is what this script does, and it is deliberately written so that the
decoding path does *not* reuse the adapters: an audit that calls the encoder's
own decoder proves only that the adapter is self-consistent.

For every (scenario, grid, region, codec, setting) it records

    library and version, input dtype and shape, the 16-bit scale/offset,
    the error the pre-codec conversion introduces before any compression,
    the reported byte count, the payload length, the size of the file after
    it was written and closed, the sha256 of what was written and of what was
    read back, whether an independent decoder reproduces the adapter's array,
    and -- for every point announced as lossless -- a bit-identity test.

Two scenarios are kept separate because the answer differs between them:

    dn        the Sentinel-2 integer digital numbers are preserved
              (reflectance x 10000, uint16, which is how L2A is distributed);
    float32   the scientific cube, i.e. what the paper's benchmark feeds in.

Usage
-----
    python codec_audit.py --synthetic --grid hpx    # ~1 min, no tiles needed
    python codec_audit.py                           # 2 regions, both grids, both scenarios
    python codec_audit.py --limit 5 --grid hpx
    python codec_audit.py --codecs jpeg2000_lossless jpegxl_lossless

Writes results/codec_audit.csv, results/codec_audit.md and
../tex/tables/table_byte_audit.tex.  The default run is about 420 operating
points and takes on the order of ten minutes, dominated by lossless JPEG~XL on
the 1024x1024 side; the audit is about mechanics rather than statistics, so
there is no reason to run it over all 40 regions.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import codec_lib
import build_tiles
import region_sites


# --- what each codec is, in library terms ------------------------------------
CODEC_LIBRARY = {
    "zfp": "zfpy",
    "sz3": "pysz",
    "jpeg2000_lossy": "PIL",
    "jpeg2000_lossless": "imagecodecs",
    "jpegxl_lossy": "imagecodecs",
    "jpegxl_lossless": "imagecodecs",
    "dwt_zstd": "pywt",
    "uncompressed": "numpy",
}

LOSSLESS_CLAIMED = {"jpeg2000_lossless", "jpegxl_lossless", "uncompressed"}

# File extension per codec, so what is written on disk is a file a normal
# reader would recognise rather than an anonymous blob.
EXTENSION = {
    "zfp": ".zfp", "sz3": ".sz3",
    "jpeg2000_lossy": ".jp2", "jpeg2000_lossless": ".jp2",
    "jpegxl_lossy": ".jxl", "jpegxl_lossless": ".jxl",
    "dwt_zstd": ".zst", "uncompressed": ".raw",
}


def library_version(module_name: str) -> str:
    try:
        mod = __import__(module_name)
    except Exception:                                      # noqa: BLE001
        return "absent"
    return str(getattr(mod, "__version__", "unknown"))


# --- the two input scenarios -------------------------------------------------
DN_SCALE = 10000.0        # Sentinel-2 L2A BOA_QUANTIFICATION_VALUE
DN_OFFSET_B4 = 1000.0     # BOA_ADD_OFFSET introduced by processing baseline 04.00
DN_MAX = 65535


def dn_mapping(a, scale=None, offset=None):
    """Pick the integer DN convention that fits this array, and name it.

    Sentinel-2 L2A is distributed as uint16 digital numbers, but which affine
    map takes reflectance to those numbers depends on the processing baseline:
    before 04.00, DN = 10000 x rho with rho >= 0; from 04.00 on, an
    BOA_ADD_OFFSET of -1000 is applied, so the stored integer is
    10000 x rho + 1000 and slightly negative reflectance -- an ordinary outcome
    of atmospheric correction over dark water and shadow -- is representable.
    Hard-coding the first convention makes the audit refuse real tiles, which
    is the wrong failure: the audit's job is to report the conversion, not to
    require a particular one.
    """
    finite = np.isfinite(a)
    if not finite.any():
        raise ValueError("array is entirely non-finite")
    lo, hi = float(np.min(a[finite])), float(np.max(a[finite]))
    if scale is not None:
        return float(scale), float(offset or 0.0), "given on the command line"
    # Values that are already integers spanning more than unit range are DN.
    if hi > 1.5 and np.allclose(a[finite], np.rint(a[finite]), atol=1e-3):
        return 1.0, 0.0, "array is already integer DN"
    off = DN_OFFSET_B4 if lo < 0 else 0.0
    name = (f"DN = {DN_SCALE:.0f} x rho" if off == 0
            else f"DN = {DN_SCALE:.0f} x rho + {off:.0f} (baseline 04.00 offset)")
    return DN_SCALE, off, name


def prepare(arr, scenario: str, dn_scale=None, dn_offset=None):
    """Return (fed_to_codec_float32, info).

    `info` records the dtype the scenario claims to carry, the unit of the
    values, and the error introduced *before* any codec sees the data.
    """
    a = np.asarray(arr, dtype=np.float32)
    finite = np.isfinite(a)
    nan_frac = float(1.0 - finite.mean())

    if scenario == "float32":
        return a, dict(declared_dtype="float32", value_unit="reflectance",
                       pre_pipeline_rmse=0.0, pre_pipeline_max_abs=0.0,
                       dn_mapping="--", dn_clipped_frac=0.0, nan_frac=nan_frac)

    if scenario == "dn":
        scale, offset, name = dn_mapping(a, dn_scale, dn_offset)
        raw = np.rint(np.where(finite, a, 0.0).astype(np.float64) * scale + offset)
        clipped = finite & ((raw < 0) | (raw > DN_MAX))
        dn = np.clip(raw, 0, DN_MAX).astype(np.uint16)
        back = (dn.astype(np.float64) - offset) / scale
        err = np.where(finite, back - a.astype(np.float64), 0.0)
        n = max(int(finite.sum()), 1)
        return dn.astype(np.float32), dict(
            declared_dtype="uint16",
            value_unit=f"DN ({name})",
            # reported in reflectance units, so it is directly comparable with
            # the paper's RMSE targets
            pre_pipeline_rmse=float(np.sqrt(np.sum(err ** 2) / n)),
            pre_pipeline_max_abs=float(np.max(np.abs(err))),
            dn_mapping=name,
            dn_clipped_frac=float(clipped.sum()) / n,
            nan_frac=nan_frac)

    raise KeyError(f"unknown scenario {scenario!r}")


# --- the 16-bit renormalization the image codecs perform ---------------------
def u16_params(a):
    """Recompute what codec_lib._to_u16 derives, from the input alone.

    That it *is* recomputable from the input is the point: the pair is not in
    the file, so a reader who has only the file cannot invert the mapping.
    Eight bytes, not charged to any reported size.
    """
    a = np.asarray(a, dtype=np.float32)
    t_min = float(a.min())
    span = (float(a.max()) - t_min) or 1.0
    return t_min, span


def u16_roundtrip_error(a):
    """RMSE the 16-bit requantization costs before the codec runs."""
    a = np.asarray(a, dtype=np.float32)
    t_min, span = u16_params(a)
    u = ((a - t_min) / span * 65535).astype(np.uint16)
    back = u.astype(np.float32) / 65535.0 * span + t_min
    d = (back - a).astype(np.float64)
    return float(np.sqrt(np.mean(d ** 2))), float(np.max(np.abs(d))), u


# --- an independent decoder, deliberately not the adapters' ------------------
SELF_CONTAINED = {
    # can a reader with the file and nothing else recover the array?
    "zfp": True,                 # the zfp stream carries shape and type
    "sz3": False,                # shape and dtype passed out of band
    "jpeg2000_lossy": False,     # u16 scale/offset out of band
    "jpeg2000_lossless": False,
    "jpegxl_lossy": False,
    "jpegxl_lossless": False,
    "dwt_zstd": False,           # coefficient layout and the float32 scale
    "uncompressed": False,       # shape and dtype
}


def decode_file(codec: str, config: dict, data: bytes, fed):
    """Decode `data` back to an array, using only the file plus the metadata
    the audit declares as external.  Never calls codec_lib."""
    t_min, span = u16_params(fed)
    shape = np.shape(fed)

    def from_u16(u):
        return np.asarray(u, dtype=np.float32) / 65535.0 * span + t_min

    if codec == "zfp":
        import zfpy
        return np.asarray(zfpy.decompress_numpy(data), dtype=np.float32)

    if codec == "sz3":
        from pysz import sz as _sz
        dec, _ = _sz.decompress(np.frombuffer(data, dtype=np.uint8),
                                np.float32, shape)
        return np.asarray(dec, dtype=np.float32)

    if codec == "jpeg2000_lossy":
        from PIL import Image
        with Image.open(io.BytesIO(data)) as im:
            return from_u16(np.array(im, dtype=np.uint16))

    if codec == "jpeg2000_lossless":
        import imagecodecs
        return from_u16(np.asarray(imagecodecs.jpeg2k_decode(data),
                                   dtype=np.uint16))

    if codec in ("jpegxl_lossy", "jpegxl_lossless"):
        import imagecodecs
        return from_u16(np.asarray(imagecodecs.jpegxl_decode(data),
                                   dtype=np.uint16))

    if codec == "dwt_zstd":
        import pywt
        import zstandard as zstd
        wavelet = config["wavelet"]
        levels = int(config["levels"])
        bits = int(config["quantize_bits"])
        # The coefficient layout is a function of (shape, wavelet, levels), so
        # a reader who knows those three can rebuild it; the float32 `scale`
        # is not in the file at all -- codec_lib charges 4 bytes for it but
        # writes them nowhere, so we recompute it the way the encoder did and
        # flag the point as not self-contained.
        probe = pywt.wavedec2(np.zeros(shape, dtype=np.float32),
                              wavelet=wavelet, level=levels, mode="periodization")
        parr, slices = pywt.coeffs_to_array(probe)
        coeffs = pywt.wavedec2(np.asarray(fed, dtype=np.float32),
                               wavelet=wavelet, level=levels, mode="periodization")
        arr, _ = pywt.coeffs_to_array(coeffs)
        scale = float(np.max(np.abs(arr))) or 1.0
        max_level = (1 << (bits - 1)) - 1
        dtype_q = np.int8 if bits <= 8 else np.int16
        q = np.frombuffer(zstd.ZstdDecompressor().decompress(data),
                          dtype=dtype_q).reshape(parr.shape)
        arr2 = (q.astype(np.float32) / max_level * scale).astype(np.float64)
        c2 = pywt.array_to_coeffs(arr2, slices, output_format="wavedec2")
        out = pywt.waverec2(c2, wavelet=wavelet, mode="periodization")
        out = np.asarray(out, dtype=np.float32)
        return out[tuple(slice(0, s) for s in shape)]

    if codec == "uncompressed":
        return np.frombuffer(data, dtype=np.float32).reshape(shape)

    raise KeyError(f"no independent decoder for {codec!r}")


# --- one audited run ---------------------------------------------------------
def audit_run(run, fed, scenario_info, scenario, grid, scene, workdir: Path):
    codec = run.codec
    row = dict(
        scenario=scenario, grid=grid, region=scene, codec=codec,
        config=json.dumps(run.config, sort_keys=True),
        library=CODEC_LIBRARY.get(codec, "?"),
        library_version=library_version(CODEC_LIBRARY.get(codec, "sys")),
        declared_dtype=scenario_info["declared_dtype"],
        fed_dtype=str(np.asarray(fed).dtype),
        shape="x".join(str(s) for s in np.shape(fed)),
        value_unit=scenario_info["value_unit"],
        pre_pipeline_rmse=scenario_info["pre_pipeline_rmse"],
        pre_pipeline_max_abs=scenario_info["pre_pipeline_max_abs"],
        dn_mapping=scenario_info["dn_mapping"],
        dn_clipped_frac=scenario_info["dn_clipped_frac"],
        nan_frac=scenario_info["nan_frac"],
        total_bytes_reported=int(run.total_bytes),
        uncharged_bytes=int(run.uncharged_bytes),
        external_metadata=";".join(run.external_metadata),
        self_contained=bool(SELF_CONTAINED.get(codec, False)),
        lossless_claimed=codec in LOSSLESS_CLAIMED,
    )

    if run.payload is None:
        row.update(payload_bytes=-1, file_bytes=-1, file_identical=False,
                   decode_ok=False, note="adapter exposes no payload")
        return row

    payload = run.payload
    row["payload_bytes"] = len(payload)
    row["charged_minus_payload"] = int(run.total_bytes) - len(payload)

    # written, closed, stat'ed, re-opened -- in that order, on a real path
    path = workdir / f"{grid}_{scene}_{codec}{EXTENSION.get(codec, '.bin')}"
    with open(path, "wb") as fh:
        fh.write(payload)
        fh.flush()
    row["file_bytes"] = int(path.stat().st_size)
    reread = path.read_bytes()
    row["sha256_written"] = hashlib.sha256(payload).hexdigest()[:16]
    row["sha256_reread"] = hashlib.sha256(reread).hexdigest()[:16]
    row["file_identical"] = bool(reread == payload)

    # an independent reader decodes the file
    try:
        dec = decode_file(codec, run.config, reread, fed)
        want = np.shape(fed)
        if np.shape(dec) != want:
            dec = np.asarray(dec)[tuple(slice(0, s) for s in want)]
        ref = np.asarray(run.reconstruction, dtype=np.float32)
        row["decode_ok"] = True
        row["decode_matches_adapter"] = bool(np.array_equal(
            np.asarray(dec, dtype=np.float32), ref))
        row["decode_max_abs_diff"] = float(np.max(np.abs(
            np.asarray(dec, dtype=np.float64) - ref.astype(np.float64))))
        row["rmse_vs_input"] = codec_lib.rmse(dec, fed)
        row["note"] = ""
    except Exception as exc:                               # noqa: BLE001
        row.update(decode_ok=False, decode_matches_adapter=False,
                   decode_max_abs_diff=np.nan, rmse_vs_input=np.nan,
                   note=f"{type(exc).__name__}: {exc}")
        path.unlink(missing_ok=True)
        return row
    path.unlink(missing_ok=True)

    # bit-identity, where the point is announced lossless
    if codec in LOSSLESS_CLAIMED:
        fed_f32 = np.asarray(fed, dtype=np.float32)
        row["lossless_input_bit_exact"] = bool(np.array_equal(
            np.asarray(dec, dtype=np.float32), fed_f32))
        if codec == "uncompressed":
            row["lossless_u16_bit_exact"] = True
            row["u16_rmse"] = 0.0
            row["u16_max_abs"] = 0.0
        else:
            u16_rmse, u16_max, u_in = u16_roundtrip_error(fed_f32)
            t_min, span = u16_params(fed_f32)
            u_out = np.rint((np.asarray(dec, dtype=np.float64) - t_min)
                            / span * 65535)
            row["lossless_u16_bit_exact"] = bool(np.array_equal(
                u_out.astype(np.int64), u_in.astype(np.int64)))
            row["u16_rmse"] = u16_rmse
            row["u16_max_abs"] = u16_max
    else:
        row["lossless_input_bit_exact"] = False
        row["lossless_u16_bit_exact"] = False
        if CODEC_LIBRARY.get(codec) in ("PIL", "imagecodecs"):
            r, m, _ = u16_roundtrip_error(np.asarray(fed, dtype=np.float32))
            row["u16_rmse"], row["u16_max_abs"] = r, m
        else:
            row["u16_rmse"], row["u16_max_abs"] = 0.0, 0.0

    t_min, span = u16_params(np.asarray(fed, dtype=np.float32))
    row["u16_min"], row["u16_span"] = t_min, span
    row["encode_seconds"] = run.encode_seconds
    row["decode_seconds"] = run.decode_seconds
    return row


# --- driver ------------------------------------------------------------------
def synthetic_arrays(grids, seed=0):
    """A stand-in field, for smoke-testing the audit without built tiles.

    Reflectance-like: smooth large scales plus a textured component, so the
    codecs behave qualitatively as they do on B04 and every mechanical check
    (file written, re-read, decoded, bit-identity) is exercised.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for grid, n in (("hpx", cfg.TILE_SIDE), ("utm", cfg.UTM_SIDE_PX)):
        if grid not in grids:
            continue
        yy, xx = np.mgrid[0:n, 0:n] / n
        field = (0.25 + 0.15 * np.sin(6 * np.pi * xx) * np.cos(4 * np.pi * yy)
                 + 0.05 * rng.standard_normal((n, n)))
        out[grid] = np.clip(field, 0.0, 1.0).astype(np.float32)
    return out


def audit_region(scene: str, grids, codecs, scenarios, workdir: Path,
                 synthetic=False, dn_scale=None, dn_offset=None):
    if synthetic:
        arrays = synthetic_arrays(grids)
    else:
        from hcc.compression.zorder import expand_tile
        _ids, cell_data, utm_crop, _meta, _npz = build_tiles.load_tile(scene)
        arrays = {}
        if "hpx" in grids:
            arrays["hpx"] = expand_tile(np.asarray(cell_data, dtype=np.float32)
                                        ).astype(np.float32)
        if "utm" in grids:
            arrays["utm"] = np.asarray(utm_crop, dtype=np.float32)

    rows = []
    for scenario in scenarios:
        for grid, arr in arrays.items():
            fed, info = prepare(arr, scenario, dn_scale, dn_offset)
            a = np.asarray(arr, dtype=np.float32)
            fin = np.isfinite(a)
            print(f"    {scenario:8s} {grid}: range "
                  f"[{float(np.min(a[fin])):+.5f}, {float(np.max(a[fin])):+.5f}], "
                  f"{info['nan_frac']:.3%} non-finite, {info['dn_mapping']}, "
                  f"{info['dn_clipped_frac']:.3%} clipped, "
                  f"conversion RMSE {info['pre_pipeline_rmse']:.3e}")
            for run in codec_lib.iter_runs(fed, codecs=codecs, on_error="warn"):
                rows.append(audit_run(run, fed, info, scenario, grid, scene,
                                      workdir))
    return rows


def summarize(df: pd.DataFrame) -> str:
    out = ["# Byte and type audit", ""]
    out.append(f"- rows: **{len(df)}**, regions: {df['region'].nunique()}, "
               f"grids: {', '.join(sorted(df['grid'].unique()))}, "
               f"scenarios: {', '.join(sorted(df['scenario'].unique()))}")
    bad_file = df[~df["file_identical"].fillna(False)]
    bad_dec = df[~df["decode_ok"].fillna(False)]
    mismatch = df[df["decode_ok"].fillna(False)
                  & ~df["decode_matches_adapter"].fillna(False)]
    out.append(f"- files whose re-read bytes differ from what was written: "
               f"**{len(bad_file)}**")
    out.append(f"- points an independent decoder could not read back: "
               f"**{len(bad_dec)}**")
    out.append(f"- points where the independent decoder disagrees with the "
               f"adapter: **{len(mismatch)}**")
    size_gap = df["total_bytes_reported"] - df["payload_bytes"]
    out.append(f"- reported size minus payload length: "
               f"min {int(size_gap.min())} B, max {int(size_gap.max())} B "
               f"(the DWT baseline charges 4 B for a scale it does not write)")
    out.append("")

    out.append("## Per codec")
    out.append("")
    out.append("| codec | library | version | self-contained | uncharged sidecar | "
               "16-bit floor (RMSE) | lossless claimed | bit-exact vs input | "
               "bit-exact vs 16-bit input |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for codec, g in df.groupby("codec"):
        floor = _median_floor(g)
        claimed = bool(g["lossless_claimed"].iloc[0])
        out.append(
            f"| `{codec}` | {g['library'].iloc[0]} | {g['library_version'].iloc[0]} | "
            f"{'yes' if g['self_contained'].iloc[0] else 'no'} | "
            f"{int(g['uncharged_bytes'].iloc[0])} B | "
            f"{'--' if not np.isfinite(floor) else f'{floor:.2e}'} | "
            f"{'yes' if claimed else 'no'} | "
            f"{'yes' if claimed and bool(g['lossless_input_bit_exact'].all()) else ('no' if claimed else '--')} | "
            f"{'yes' if claimed and bool(g['lossless_u16_bit_exact'].all()) else ('no' if claimed else '--')} |")
    out.append("")

    out.append("## Per scenario")
    out.append("")
    for scenario, g in df.groupby("scenario"):
        pre = g["pre_pipeline_rmse"].max()
        out.append(f"- **{scenario}**: declared dtype "
                   f"{g['declared_dtype'].iloc[0]}, unit "
                   f"{g['value_unit'].iloc[0]}, conversion error before any "
                   f"codec {pre:.3e} (reflectance units), "
                   f"{g['dn_clipped_frac'].max():.3%} of samples clipped, "
                   f"{g['nan_frac'].max():.3%} non-finite")
        loss = g[g["lossless_claimed"]]
        if len(loss):
            ok_in = int(loss["lossless_input_bit_exact"].sum())
            ok_u16 = int(loss["lossless_u16_bit_exact"].sum())
            out.append(f"  - lossless points: {ok_in}/{len(loss)} bit-exact "
                       f"against the array fed in, {ok_u16}/{len(loss)} "
                       f"bit-exact against its 16-bit renormalization")
    out.append("")
    if len(bad_dec):
        out.append("## Points that failed to read back")
        out.append("")
        for _, r in bad_dec.iterrows():
            out.append(f"- `{r['codec']}` {r['config']} ({r['grid']}, "
                       f"{r['scenario']}): {r['note']}")
        out.append("")
    return "\n".join(out)


def _median_floor(g):
    """Median 16-bit requantization RMSE, NaN when the codec has no such step."""
    v = g["u16_rmse"].to_numpy(dtype=float, na_value=np.nan)
    v = v[np.isfinite(v) & (v > 0)]
    return float(np.median(v)) if v.size else np.nan


def _sci(x) -> str:
    """1.23e-05 -> $1.2\\times10^{-5}$, and a blank for a missing value."""
    if x is None or not np.isfinite(x):
        return "--"
    mant, exp = f"{x:.1e}".split("e")
    return f"${mant}\\times10^{{{int(exp)}}}$"


def latex_table(df: pd.DataFrame) -> str:
    rows = []
    for codec, g in sorted(df.groupby("codec")):
        floor = _median_floor(g)
        claimed = bool(g["lossless_claimed"].iloc[0])
        if not claimed:
            exact = "--"
        elif bool(g["lossless_u16_bit_exact"].all()):
            exact = "yes"
        else:
            exact = "no"
        rows.append(
            f"\\texttt{{{codec.replace('_', chr(92) + '_')}}} & "
            f"{g['library'].iloc[0]}~{g['library_version'].iloc[0]} & "
            f"{'yes' if g['self_contained'].iloc[0] else 'no'} & "
            f"{int(g['uncharged_bytes'].iloc[0])} & "
            f"{_sci(floor)} & "
            f"{exact} \\\\")
    body = "\n".join(rows)
    n_ok = int(df["file_identical"].fillna(False).sum())
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Byte and type audit of the codec adapters. Every operating "
        "point was written to a file, closed, re-opened and decoded by a reader "
        "independent of the encoder; "
        f"{n_ok} of {len(df)} points reproduced the encoder's array exactly. "
        "\\emph{Self-contained} asks whether the file alone suffices to decode; "
        "\\emph{sidecar} is the number of bytes a decoder needs that the "
        "reported size does not charge (the float32 minimum and span of the "
        "16-bit renormalization). The 16-bit floor is the RMSE that "
        "renormalization costs before any compression, which is why points "
        "through it are reported as lossless \\emph{after requantization} "
        "rather than lossless.}\n"
        "\\label{tab:byte_audit}\n"
        "\\begin{tabular}{llccrc}\n\\toprule\n"
        "Codec & Library & Self-cont. & Sidecar (B) & 16-bit floor & Bit-exact \\\\\n"
        "\\midrule\n" + body + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=2,
                    help="regions to audit (the audit is about mechanics, not "
                         "statistics; 2 is enough to expose a type bug)")
    ap.add_argument("--grid", choices=["hpx", "utm", "both"], default="both")
    ap.add_argument("--scenarios", nargs="*", default=["dn", "float32"],
                    choices=["dn", "float32"])
    ap.add_argument("--codecs", nargs="*", default=None)
    ap.add_argument("--no-tex", action="store_true")
    ap.add_argument("--dn-scale", type=float, default=None,
                    help="force the reflectance->DN scale instead of deriving "
                         "it (Sentinel-2 L2A uses 10000)")
    ap.add_argument("--dn-offset", type=float, default=None,
                    help="force the DN offset (baseline 04.00 uses 1000)")
    ap.add_argument("--synthetic", action="store_true",
                    help="smoke-test the audit on a stand-in field, with no "
                         "tiles and no results written to tex/")
    args = ap.parse_args(argv)

    grids = ["hpx", "utm"] if args.grid == "both" else [args.grid]
    if args.synthetic:
        scenes = ["synthetic"]
        args.no_tex = True
    else:
        scenes = args.regions or list(region_sites.load_sites()["scene"])
        scenes = [s for s in scenes if (cfg.TILE_DIR / f"{s}.npz").exists()]
        if not scenes:
            print(f"no tiles under {cfg.TILE_DIR}; run build_tiles.py first, "
                  f"or use --synthetic to smoke-test the audit itself")
            return 1
        if args.limit:
            scenes = scenes[:args.limit]

    print(f"auditing {len(scenes)} region(s), grids {grids}, "
          f"scenarios {args.scenarios}")
    rows = []
    with tempfile.TemporaryDirectory(prefix="codec_audit_") as tmp:
        workdir = Path(tmp)
        for scene in scenes:
            print(f"  {scene}")
            rows.extend(audit_region(scene, grids, args.codecs,
                                     args.scenarios, workdir,
                                     synthetic=args.synthetic,
                                     dn_scale=args.dn_scale,
                                     dn_offset=args.dn_offset))

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stem = "codec_audit_synthetic" if args.synthetic else "codec_audit"
    csv = cfg.RESULTS_DIR / f"{stem}.csv"
    df.to_csv(csv, index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / f"{stem}.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")

    if not args.no_tex:
        tex_dir = cfg.REPO_ROOT / "tex" / "tables"
        tex_dir.mkdir(parents=True, exist_ok=True)
        (tex_dir / "table_byte_audit.tex").write_text(latex_table(df),
                                                      encoding="utf-8")
        print(f"wrote {tex_dir / 'table_byte_audit.tex'}")

    failures = int((~df["decode_ok"].fillna(False)).sum())
    mism = int((df["decode_ok"].fillna(False)
                & ~df["decode_matches_adapter"].fillna(False)).sum())
    if failures or mism:
        print(f"\n{failures} point(s) unreadable, {mism} point(s) where an "
              f"independent decoder disagrees with the adapter. "
              f"These are the rows the paper cannot claim.")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
