import numpy as np
import pytest
from zarr.buffer.cpu import numpy_buffer_prototype
from zarr.core.array_spec import ArrayConfig, ArraySpec
from zarr.dtype import Float32, parse_dtype

from hcc.codecs import grad_dwt
from hcc.codecs.decompose import ArrayConfigWithScalars
from hcc.tests import is_noise


@pytest.mark.parametrize("levels", [3, 4])
@pytest.mark.parametrize("quantize_bits", [8, 16])
@pytest.mark.parametrize("clip_percentile", [0.5, 0.1])
def test_codec_serialization_roundtrip(levels, quantize_bits, clip_percentile):
    codec = grad_dwt.HccGradDwtCodec(
        wavelet="db4",
        levels=levels,
        quantize_bits=quantize_bits,
        clip_percentile=clip_percentile,
    )

    serialized = codec.to_dict()
    reconstructed = grad_dwt.HccGradDwtCodec.from_dict(serialized)

    assert codec == reconstructed


@pytest.mark.parametrize("quantize_bits", [8, 16])
def test_resolve_metadata(quantize_bits):
    prototype = numpy_buffer_prototype()
    chunk_spec = ArraySpec(
        shape=(64, 64),
        dtype=Float32(endianness="little"),
        prototype=prototype,
        fill_value=np.nan,
        config=ArrayConfig.from_dict({}),
    )
    codec = grad_dwt.HccGradDwtCodec(
        wavelet="db4", levels=3, quantize_bits=quantize_bits, clip_percentile=0.5
    )

    actual = codec.resolve_metadata(chunk_spec)
    dtype = f"int{quantize_bits}"

    assert actual.prototype == prototype
    assert actual.shape == chunk_spec.shape
    assert actual.dtype == parse_dtype(dtype, zarr_format=3)
    assert actual.fill_value == np.iinfo(dtype).min


def _make_tile(edge_length=64):
    x = np.linspace(0, 4, edge_length, dtype="float32")
    y = np.linspace(0, 4, edge_length, dtype="float32")
    xx, yy = np.meshgrid(x, y)
    return (np.sin(xx) * np.cos(yy)).astype("float32")


@pytest.mark.asyncio
@pytest.mark.parametrize("levels", [2, 3])
@pytest.mark.parametrize("quantize_bits", [8, 16])
@pytest.mark.parametrize("clip_percentile", [0.5, 0.1])
async def test_codec_standalone_fulltile_roundtrip(levels, quantize_bits, clip_percentile):
    """Recommended usage (see module docstring): the codec applied directly
    to a tile, with no HccDecomposeCodec wrapper and no
    `chunk_spec.config.context["prediction"]` populated."""
    data = _make_tile(64)

    codec = grad_dwt.HccGradDwtCodec(
        wavelet="db4",
        levels=levels,
        quantize_bits=quantize_bits,
        clip_percentile=clip_percentile,
    )

    prototype = numpy_buffer_prototype()
    array_spec = ArraySpec(
        shape=data.shape,
        dtype=Float32(endianness="little"),
        prototype=prototype,
        fill_value=np.nan,
        config=ArrayConfigWithScalars.from_dict({}),
    )
    # Deliberately NOT setting array_spec.config.context["prediction"] --
    # standalone usage must not require it.

    input_buffer = prototype.nd_buffer.from_ndarray_like(data)
    encoded = await codec._encode_single(input_buffer, array_spec)
    decoded = await codec._decode_single(encoded, array_spec)
    roundtripped = decoded.as_numpy_array()

    assert roundtripped.shape == data.shape
    assert roundtripped.dtype == data.dtype

    ravelled_data = np.ravel(data)
    ravelled_roundtripped = np.ravel(roundtripped)
    assert (
        is_noise(ravelled_data) and is_noise(ravelled_roundtripped)
    ) or np.corrcoef(ravelled_roundtripped, ravelled_data)[0, 1] >= 0.5

    # both sidecar scalars must have been produced
    assert set(array_spec.config.scalars.keys()) == {"q_scale", "ll_values"}


@pytest.mark.asyncio
async def test_codec_legacy_residual_path_still_works():
    """Backward compatibility: still usable inside a residual pipeline with
    context["prediction"] populated (the value is simply ignored)."""
    data = _make_tile(64)
    block_shape = 32
    downsampled = np.mean(
        np.reshape(data, (block_shape, 2, block_shape, 2)), axis=(1, 3)
    )
    prediction = np.reshape(
        np.broadcast_to(downsampled[:, None, :, None], (block_shape, 2, block_shape, 2)),
        data.shape,
    ).astype("float32")
    residuals = (data - prediction).astype("float32")

    codec = grad_dwt.HccGradDwtCodec(
        wavelet="db4", levels=3, quantize_bits=8, clip_percentile=0.5
    )

    prototype = numpy_buffer_prototype()
    array_spec = ArraySpec(
        shape=residuals.shape,
        dtype=Float32(endianness="little"),
        prototype=prototype,
        fill_value=np.nan,
        config=ArrayConfigWithScalars.from_dict({}),
    )
    array_spec.config.context["prediction"] = prediction

    input_buffer = prototype.nd_buffer.from_ndarray_like(residuals)
    encoded = await codec._encode_single(input_buffer, array_spec)
    decoded = await codec._decode_single(encoded, array_spec)
    roundtripped = decoded.as_numpy_array()

    assert roundtripped.shape == residuals.shape
    assert np.corrcoef(np.ravel(roundtripped), np.ravel(residuals))[0, 1] >= 0.3


@pytest.mark.parametrize("name", ["low", "medium", "high", "ultra"])
def test_lookup_preset_legacy(name):
    from hcc.codecs.decompose import HccDecomposeCodec
    from hcc.compression.grad_dwt import lookup_preset as lookup_raw_preset

    presets = lookup_raw_preset(name, codec=True)
    codec = grad_dwt.lookup_preset(name)

    assert isinstance(codec, HccDecomposeCodec)
    assert codec.factor == presets["downsample_factor"]

    grad_dwt_codec = codec.codecs["residuals"][0]
    assert grad_dwt_codec.wavelet == presets["wavelet"]
    assert grad_dwt_codec.levels == presets["levels"]
    assert grad_dwt_codec.quantize_bits == presets["quantize_bits"]
    assert grad_dwt_codec.clip_percentile == presets["clip_percentile"]


@pytest.mark.parametrize("name", ["low", "medium", "high", "ultra"])
def test_build_fulltile_pipeline(name):
    from hcc.compression.grad_dwt import lookup_preset as lookup_raw_preset

    presets = lookup_raw_preset(name, codec=False)
    codecs = grad_dwt.build_fulltile_pipeline(name)

    assert isinstance(codecs, list)
    dwt_codec = codecs[0]
    assert isinstance(dwt_codec, grad_dwt.HccGradDwtCodec)
    assert dwt_codec.wavelet == presets["wavelet"]
    assert dwt_codec.levels == presets["levels"]
    assert dwt_codec.quantize_bits == presets["quantize_bits"]
    assert dwt_codec.clip_percentile == presets["clip_percentile"]
    # last codec in the chain should be the gzip compressor carrying the
    # preset's zlib_level
    assert codecs[-1].level == presets["zlib_level"]


@pytest.mark.asyncio
@pytest.mark.parametrize("name", ["low", "medium"])
async def test_build_fulltile_pipeline_roundtrip(name):
    """End-to-end sanity check of the recommended standalone pipeline: run
    just the HccGradDwtCodec stage (the ArrayArrayCodec) from the pipeline
    build_fulltile_pipeline() returns, standalone on a full tile."""
    data = _make_tile(64)
    codecs = grad_dwt.build_fulltile_pipeline(name)
    dwt_codec = codecs[0]

    prototype = numpy_buffer_prototype()
    array_spec = ArraySpec(
        shape=data.shape,
        dtype=Float32(endianness="little"),
        prototype=prototype,
        fill_value=np.nan,
        config=ArrayConfigWithScalars.from_dict({}),
    )

    input_buffer = prototype.nd_buffer.from_ndarray_like(data)
    encoded = await dwt_codec._encode_single(input_buffer, array_spec)
    decoded = await dwt_codec._decode_single(encoded, array_spec)
    roundtripped = decoded.as_numpy_array()

    assert roundtripped.shape == data.shape
    assert np.corrcoef(np.ravel(roundtripped), np.ravel(data))[0, 1] >= 0.5
