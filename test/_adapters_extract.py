# --- codec adapters ----------------------------------------------------------
from dataclasses import dataclass


@dataclass
class CodecRunResult:
    codec: str
    config: dict
    total_bytes: int
    encode_seconds: float
    decode_seconds: float
    reconstruction: np.ndarray


def _rmse(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2))) if mask.any() else np.nan


def encode_zfp(tile, *, tolerance):
    import zfpy
    tile = np.asarray(tile, dtype=np.float32)
    t0 = time.perf_counter()
    enc = zfpy.compress_numpy(tile, tolerance=float(tolerance))
    t1 = time.perf_counter()
    dec = zfpy.decompress_numpy(enc)
    t2 = time.perf_counter()
    return CodecRunResult("zfp", dict(tolerance=float(tolerance)), len(enc), t1 - t0, t2 - t1, dec)


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
    return CodecRunResult("sz3", dict(abs_error=float(abs_error)), total_bytes, t1 - t0, t2 - t1,
                           np.asarray(dec, dtype=np.float32))


def _to_u16(tile):
    tile = np.asarray(tile, dtype=np.float32)
    t_min, t_max = float(tile.min()), float(tile.max())
    span = (t_max - t_min) or 1.0
    return ((tile - t_min) / span * 65535).astype(np.uint16), t_min, span


def _from_u16(u16, t_min, span):
    return u16.astype(np.float32) / 65535.0 * span + t_min


def encode_jpeg2000_lossy(tile, *, quality):
    import io
    from PIL import Image
    u16, t_min, span = _to_u16(tile)
    factor = max(1, int(quality))  # never let it round to 0 (past bug, see main notebook)
    buf = io.BytesIO()
    t0 = time.perf_counter()
    Image.fromarray(u16).save(buf, format="jpeg2000", quality_mode="rates", quality_layers=[factor])
    t1 = time.perf_counter()
    total_bytes = buf.tell()
    buf.seek(0)
    dec_u16 = np.array(Image.open(buf), dtype=np.uint16)
    t2 = time.perf_counter()
    return CodecRunResult("jpeg2000_lossy", dict(quality=factor), total_bytes, t1 - t0, t2 - t1,
                           _from_u16(dec_u16, t_min, span))


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


def encode_jpegxl_lossy(tile, *, level):
    # imagecodecs' `level` for jpegxl follows its shared 0-100 convention
    # across codecs (100 ~ visually lossless / near-lossless, not literally
    # lossless -- see encode_jpegxl_lossless for that). If the installed
    # imagecodecs version uses a different kwarg name, this fails loudly
    # once with a clear message instead of silently corrupting the sweep.
    global _jpegxl_warned
    import imagecodecs
    u16, t_min, span = _to_u16(tile)
    t0 = time.perf_counter()
    try:
        enc = imagecodecs.jpegxl_encode(u16, level=int(level))
    except TypeError as e:
        if not _jpegxl_warned:
            print(f"[jpegxl] encode signature mismatch ({e}) -- check your imagecodecs "
                  "version's jpegxl_encode kwargs and adjust this adapter. Skipping jpegxl.")
            _jpegxl_warned = True
        raise
    t1 = time.perf_counter()
    dec_u16 = imagecodecs.jpegxl_decode(enc)
    t2 = time.perf_counter()
    return CodecRunResult("jpegxl_lossy", dict(level=int(level)), len(enc), t1 - t0, t2 - t1,
                           _from_u16(dec_u16, t_min, span))


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
    # Deliberately simple, well-known technique (PyWavelets + uniform
    # quantisation + zstd) -- a generic scale-correlation-exploiting
    # baseline, NOT the custom hcc SPIHT/zero-tree machinery this notebook
    # is moving away from.
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
    total_bytes = len(payload) + 4  # + float32 scale sidecar
    q2 = np.frombuffer(zstd.ZstdDecompressor().decompress(payload), dtype=dtype_q).reshape(arr.shape)
    arr2 = (q2.astype(np.float32) / max_level * scale).astype(np.float64)
    coeffs2 = pywt.array_to_coeffs(arr2, slices, output_format="wavedec2")
    dec = pywt.waverec2(coeffs2, wavelet=wavelet, mode="periodization").astype(np.float32)
    t2 = time.perf_counter()
    return CodecRunResult("dwt_zstd", dict(wavelet=wavelet, levels=levels, quantize_bits=quantize_bits),
                           total_bytes, t1 - t0, t2 - t1, dec)


def encode_uncompressed(tile, **_ignored):
    tile = np.asarray(tile, dtype=np.float32)
    return CodecRunResult("uncompressed", {}, tile.nbytes, 0.0, 0.0, tile.copy())