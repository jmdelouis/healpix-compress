"""Zarr ArrayArrayCodec that applies DWT compression to a tile or residual.

This module exposes :class:`HccGradDwtCodec`, a Zarr v3 ``ArrayArrayCodec``
implementing the wavelet-based per-subband quantisation scheme defined in
:mod:`hcc.compression.grad_dwt`.

Pipeline position
------------------
As of the full-tile-DWT architecture change (see the project's paper
Discussion and ``hcc.compression.grad_dwt``'s module docstring for the
empirical motivation: a controlled A/B comparison found that a prior
"downsample → bilinear-upsample prediction, then DWT the residual" design
was consistently worse -- in *both* compressed size and reconstruction
RMSE -- than applying the DWT directly to the un-decimated tile, because
the bilinear prediction step throws away exactly the low-frequency energy
the DWT's own LL subband is best at compacting, leaving the DWT to work on
an already-decorrelated, harder-to-compress signal), **the recommended
production pipeline runs this codec directly on the tile**, with no
``HccDecomposeCodec`` wrapper::

    float32 tile
        │
        ▼
    HccGradDwtCodec   ← this codec (recommended: standalone)
        │
        ▼
    int8 / int16 DWT *detail* coefficients
    (+ sidecar scalars: q_scale, ll_values)
        │
        ▼
    BytesCodec → GzipCodec / ZstdCodec

See :func:`build_fulltile_pipeline` for a ready-made ``codecs=[...]`` list.

This codec no longer depends on being nested inside
:class:`~hcc.codecs.decompose.HccDecomposeCodec` (it previously required
``chunk_spec.config.context["prediction"]`` to be pre-populated by that
wrapper, even though the DWT path never actually used the prediction
array). It can still be used inside ``HccDecomposeCodec``'s residual
sub-pipeline for callers who specifically want the two-stage design (e.g.
for its cheap coarse progressive/multi-resolution "low" layer) --
:func:`lookup_preset` still builds that pipeline shape for backward
compatibility -- but that is no longer the recommended default.

The LL approximation subband is now stored **losslessly** (as a raw
float32 sidecar, ``ll_values``) rather than being quantised like the detail
subbands -- see :func:`hcc.compression.grad_dwt.compress_residuals` for the
rationale (LL is a tiny fraction of the coefficient count but concentrates
most of the signal's energy, so quantising it is a false economy).

Compared to :class:`~hcc.codecs.grad_dct.HccGradDctCodec`:

* Uses a **multi-level 2-D DWT** instead of a block DCT, providing better
  energy compaction and no block-boundary artefacts.
* Two extra parameters: ``wavelet`` (wavelet family) and ``levels``
  (decomposition depth), replacing the ``window`` (block size) of the DCT
  codec.

Public API
----------
HccGradDwtCodec : Zarr ArrayArrayCodec implementing DWT tile/residual coding.
build_fulltile_pipeline : Recommended standalone (no-HccDecomposeCodec)
                  codec pipeline for a named quality preset.
lookup_preset   : Construct a legacy :class:`~hcc.codecs.decompose.HccDecomposeCodec`
                  (downsample + residual-DWT) pipeline from a named preset.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from zarr.abc.codec import ArrayArrayCodec
from zarr.core.common import parse_named_configuration
from zarr.dtype import Float32, Int8, Int16

from hcc.codecs.scalars import decode_scalar, encode_scalar
from hcc.compression import grad_dwt

if TYPE_CHECKING:
    from typing import Self

    from zarr.abc.buffer import Buffer, NDBuffer
    from zarr.core.array_spec import ArraySpec
    from zarr.core.common import JSON
    from zarr.dtype import ZDType

    from hcc.codecs.decompose import HccDecomposeCodec


@dataclass(frozen=True)
class HccGradDwtCodec(ArrayArrayCodec):
    """Zarr ArrayArrayCodec — DWT per-subband quantisation of a tile.

    This codec transforms a float32 array (recommended: the un-decimated
    tile itself, not a residual -- see the module docstring) into quantised
    DWT detail coefficients (int8 or int16), plus two sidecars, by:

    1. Applying a multi-level 2-D DWT (``'periodization'`` boundary mode)
       to the input.  The coefficient array has the same shape as the
       input because periodization preserves the total element count.
    2. Storing the LL approximation subband **losslessly** as a float32
       sidecar (``ll_values``).
    3. Computing an independent quantisation scale per *detail* subband
       (sidecar ``q_scale``, one value per subband) from each subband's own
       ``(100 − clip_percentile)``-th percentile of absolute coefficients.
    4. Uniformly quantising each detail subband to int8 or int16 using its
       own scale.

    Both sidecar scalars (``q_scale``, ``ll_values``) are stored in the
    packed chunk buffer and recovered automatically during decoding.

    Parameters
    ----------
    wavelet : str, default ``'db4'``
        PyWavelets wavelet name.  Recommended choices:

        * ``'db2'`` — fast, moderate quality.
        * ``'db4'`` — good default (Daubechies 4).
        * ``'db6'`` — high quality, slightly slower.
        * ``'bior4.4'`` — biorthogonal, symmetric; close to JPEG 2000's CDF
          9/7 wavelet.
    levels : int, default 4
        Number of DWT decomposition levels.  Higher values push more energy
        into the LL approximation subband (which is stored losslessly
        regardless of *quantize_bits*, so more levels also means fewer
        losslessly-stored LL samples), and give detail subbands more
        opportunity for cross-scale energy compaction.  Must not exceed
        ``floor(log2(n / filter_length))``.
    quantize_bits : {8, 16}, default 8
        Quantisation depth for detail subbands (does not affect LL, which
        is always lossless).

        * ``8``  → int8, 254 usable levels (±127), 1 byte per coefficient.
        * ``16`` → int16, 65 534 usable levels (±32 767), 2 bytes per
          coefficient.
    clip_percentile : float, default 0.5
        Tail percentile for each detail subband's ``q_scale`` derivation. A
        larger value allows more coefficient saturation in exchange for
        finer resolution of the majority.  Typical range: 0.02–5.0.

    Notes
    -----
    This codec can be used **standalone**, directly in a Zarr array's own
    codec chain (``codecs=[HccGradDwtCodec(...), BytesCodec(),
    GzipCodec()]`` -- see :func:`build_fulltile_pipeline`), which is now the
    recommended production configuration (see module docstring). It no
    longer requires ``chunk_spec.config.context["prediction"]`` to be
    populated: the *prediction* parameter threaded through
    :func:`~hcc.compression.grad_dwt.compress_residuals` /
    :func:`~hcc.compression.grad_dwt.decompress_residuals` is unused by the
    DWT path (kept only for those functions' API symmetry), so this codec
    passes the chunk array itself where a prediction would otherwise be
    expected. It can still be nested inside
    :class:`~hcc.codecs.decompose.HccDecomposeCodec`'s residual
    sub-pipeline for callers who want that legacy two-stage design.

    See Also
    --------
    hcc.compression.grad_dwt.compress_residuals : Core compression function.
    hcc.compression.grad_dwt.decompress_residuals : Core decompression function.
    build_fulltile_pipeline : Recommended standalone pipeline factory.
    hcc.codecs.decompose.HccDecomposeCodec : Legacy outer codec managing a
        separate low-frequency layer.
    hcc.codecs.grad_dct.HccGradDctCodec : Sibling codec using block DCT.
    lookup_preset : Legacy factory building a full HccDecomposeCodec-wrapped
        pipeline from a named preset.

    Examples
    --------
    Build a codec directly:

    >>> codec = HccGradDwtCodec(wavelet='db4', levels=4,
    ...                         quantize_bits=16, clip_percentile=0.1)
    >>> codec.to_dict()   # doctest: +NORMALIZE_WHITESPACE
    {'name': 'hcc-grad-dwt', 'configuration': {'wavelet': 'db4', 'levels': 4, 'quantize_bits': 16, 'clip_percentile': 0.1}}

    Build the recommended standalone (full-tile) pipeline from a preset:

    >>> codecs = build_fulltile_pipeline("high")
    >>> codecs[0].wavelet
    'bior4.4'
    """

    name: ClassVar[str] = "hcc-grad-dwt"

    wavelet: str = "db4"
    levels: int = 4
    quantize_bits: int = 8
    clip_percentile: float = 0.5

    #: Sidecar scalars produced/consumed by this codec.
    #: ``q_scale`` is the per-detail-subband quantisation step size.
    #: ``ll_values`` is the losslessly-stored LL approximation subband (see
    #: :func:`hcc.compression.grad_dwt.compress_residuals`). Both are
    #: per-chunk data-dependent values that cannot be recomputed without the
    #: original input.
    dtypes: ClassVar[dict[str, ZDType]] = {
        "q_scale": Float32(endianness="little"),
        "ll_values": Float32(endianness="little"),
    }

    @property
    def scalars(self) -> list[str]:
        """Ordered list of sidecar scalar names produced/consumed by this codec."""
        return list(self.dtypes.keys())

    @classmethod
    def from_dict(cls, data: dict[str, JSON]) -> Self:
        """Deserialise from a Zarr codec JSON object.

        Parameters
        ----------
        data : dict
            Codec JSON as stored in ``zarr.json``, e.g.::

                {
                    "name": "hcc-grad-dwt",
                    "configuration": {
                        "wavelet": "db4",
                        "levels": 4,
                        "quantize_bits": 16,
                        "clip_percentile": 0.1
                    }
                }

        Returns
        -------
        HccGradDwtCodec
        """
        _, configuration_parsed = parse_named_configuration(
            data, cls.name, require_configuration=True
        )
        configuration_parsed = configuration_parsed or {}
        return cls(**configuration_parsed)

    def to_dict(self) -> dict[str, JSON]:
        """Serialise to a Zarr codec JSON object.

        Returns
        -------
        dict
            JSON-serialisable representation for storage in ``zarr.json``.
        """
        return {
            "name": self.name,
            "configuration": {
                "wavelet": self.wavelet,
                "levels": self.levels,
                "quantize_bits": self.quantize_bits,
                "clip_percentile": self.clip_percentile,
            },
        }

    def resolve_metadata(self, chunk_spec: ArraySpec) -> ArraySpec:
        """Compute the output ``ArraySpec`` produced by this codec.

        The output dtype is int8 (for ``quantize_bits=8``) or int16 (for
        ``quantize_bits=16``); the fill value is the minimum of the integer
        dtype.  Shape is unchanged because the DWT in ``'periodization'``
        mode preserves the total element count.

        Parameters
        ----------
        chunk_spec : ArraySpec
            Input spec (float32 residuals).

        Returns
        -------
        ArraySpec
            Output spec (int8 or int16 DWT coefficients).

        Raises
        ------
        ValueError
            If :attr:`quantize_bits` is not 8 or 16.
        """
        target_dtypes: dict[int, ZDType] = {
            8: Int8(),
            16: Int16(endianness="little"),
        }
        target_zdtype = target_dtypes.get(self.quantize_bits)
        if target_zdtype is None:
            raise ValueError(
                f"quantize_bits must be 8 or 16; got {self.quantize_bits!r}."
            )
        new_fill = np.iinfo(target_zdtype.to_native_dtype()).min
        return replace(chunk_spec, dtype=target_zdtype, fill_value=new_fill)

    async def _decode_single(
        self, chunk_data: NDBuffer, chunk_spec: ArraySpec
    ) -> NDBuffer:
        """Decode a single chunk: dequantise DWT coefficients → tile.

        Reads ``q_scale`` and ``ll_values`` from ``chunk_spec.config.scalars``,
        dequantises the detail DWT coefficients, restores the losslessly
        stored LL subband, and applies a multi-level IDWT to recover the
        float32 array.

        Parameters
        ----------
        chunk_data : NDBuffer
            Quantised DWT detail coefficients (int8 or int16).
        chunk_spec : ArraySpec
            Must carry ``scalars["q_scale"]`` and ``scalars["ll_values"]``
            populated by the surrounding pipeline.

        Returns
        -------
        NDBuffer
            Reconstructed float32 array.
        """
        prototype = chunk_spec.prototype

        scalars = {
            name: decode_scalar(
                prototype, buffer, dtype=self.dtypes[name]
            ).as_numpy_array()
            for name, buffer in chunk_spec.config.scalars.items()
        }

        # `prediction` is accepted by decompress_residuals only for API
        # symmetry with the (legacy) residual-based pipeline -- the DWT
        # path never reads it. Standalone use of this codec (recommended;
        # see module docstring) has no prediction array to supply, so pass
        # the chunk's own (about-to-be-overwritten) data as a harmless
        # placeholder instead of requiring
        # chunk_spec.config.context["prediction"] to be pre-populated.
        chunk_array = grad_dwt.decompress_residuals(
            chunk_data.as_numpy_array(),
            prediction=chunk_data.as_numpy_array(),
            wavelet=self.wavelet,
            levels=self.levels,
            quantize_bits=self.quantize_bits,
            **scalars,
        )
        return prototype.nd_buffer.from_ndarray_like(chunk_array)

    async def _encode_single(
        self, chunk_array: NDBuffer, chunk_spec: ArraySpec
    ) -> Buffer | None:
        """Encode a single chunk: tile → quantised DWT detail coefficients.

        Applies a multi-level DWT and per-subband quantisation to the
        float32 input, writes ``q_scale`` and ``ll_values`` into
        ``chunk_spec.config.scalars``, and returns the quantised detail
        coefficient array (LL region zeroed -- it is carried separately in
        the ``ll_values`` sidecar).

        Parameters
        ----------
        chunk_array : NDBuffer
            Float32 input array. Recommended: the tile itself (see module
            docstring); also accepts a residual (``original −
            bilinear_prediction``) for callers using the legacy
            :class:`~hcc.codecs.decompose.HccDecomposeCodec`-wrapped
            pipeline.
        chunk_spec : ArraySpec
            Output scalars are written to ``chunk_spec.config.scalars``.

        Returns
        -------
        Buffer or None
            Quantised DWT detail coefficients (int8 or int16), or ``None``
            to skip.
        """
        prototype = chunk_spec.prototype

        # See _decode_single: `prediction` is unused by the DWT path, so no
        # context["prediction"] lookup is required for standalone use.
        coeffs_q_arr, q_scale_arr, ll_values_arr = grad_dwt.compress_residuals(
            chunk_array.as_numpy_array(),
            prediction=chunk_array.as_numpy_array(),
            wavelet=self.wavelet,
            levels=self.levels,
            quantize_bits=self.quantize_bits,
            clip_percentile=self.clip_percentile,
        )

        encoded_data = prototype.nd_buffer.from_ndarray_like(coeffs_q_arr)
        q_scale_buf = prototype.nd_buffer.from_ndarray_like(q_scale_arr)
        ll_values_buf = prototype.nd_buffer.from_ndarray_like(ll_values_arr)

        chunk_spec.config.scalars.update(
            {
                "q_scale": encode_scalar(prototype, q_scale_buf),
                "ll_values": encode_scalar(prototype, ll_values_buf),
            }
        )
        return encoded_data


def lookup_preset(name: str) -> HccDecomposeCodec:
    """Construct a full :class:`~hcc.codecs.decompose.HccDecomposeCodec`
    pipeline from a named quality preset.

    The returned pipeline uses:

    * A ``BytesCodec → GzipCodec`` sub-pipeline for the lossless
      low-frequency layer.
    * An ``HccGradDwtCodec → BytesCodec → GzipCodec`` sub-pipeline for the
      DWT-compressed residual.

    .. important::

       Store the explicit parameter dict (from
       :func:`hcc.compression.grad_dwt.lookup_preset`) in your array
       metadata rather than the preset name.  This ensures future changes to
       preset defaults cannot make existing data unreadable.

    Parameters
    ----------
    name : {'low', 'medium', 'high', 'ultra'}
        Quality preset name.

        +----------+---------+--------+---------------+-----------------+------------+-------------------+
        | name     | wavelet | levels | quantize_bits | clip_percentile | zlib_level | downsample_factor |
        +==========+=========+========+===============+=================+============+===================+
        | low      | db2     | 3      | 8             | 3.0             | 9          | 4                 |
        +----------+---------+--------+---------------+-----------------+------------+-------------------+
        | medium   | db4     | 4      | 8             | 0.5             | 9          | 2                 |
        +----------+---------+--------+---------------+-----------------+------------+-------------------+
        | high     | db4     | 5      | 16            | 0.1             | 6          | 2                 |
        +----------+---------+--------+---------------+-----------------+------------+-------------------+
        | ultra    | db6     | 6      | 16            | 0.02            | 1          | 1                 |
        +----------+---------+--------+---------------+-----------------+------------+-------------------+

    Returns
    -------
    HccDecomposeCodec
        Fully assembled pipeline.

    Raises
    ------
    ValueError
        If *name* is not a known preset.

    See Also
    --------
    hcc.compression.grad_dwt.lookup_preset : Returns the raw parameter dict.

    Examples
    --------
    >>> pipeline = lookup_preset("medium")
    >>> pipeline.factor
    2
    >>> pipeline.codecs["residuals"][0].wavelet
    'db4'
    >>> pipeline.codecs["residuals"][0].levels
    4
    """
    from zarr.codecs import BytesCodec, GzipCodec

    from hcc.codecs.decompose import HccDecomposeCodec
    from hcc.compression.grad_dwt import lookup_preset as lookup_raw_preset

    presets = lookup_raw_preset(name, codec=True)
    factor = presets.pop("downsample_factor")
    level = presets.pop("zlib_level")

    return HccDecomposeCodec(
        factor=factor,
        codecs={
            "downsampled": [BytesCodec(), GzipCodec(level=level)],
            "residuals": [
                HccGradDwtCodec(**presets),
                BytesCodec(),
                GzipCodec(level=level),
            ],
        },
    )


def build_fulltile_pipeline(name: str) -> list:
    """Build the **recommended** standalone codec pipeline for a named
    quality preset: :class:`HccGradDwtCodec` applied directly to the tile,
    with no :class:`~hcc.codecs.decompose.HccDecomposeCodec` wrapper.

    Use this (rather than :func:`lookup_preset`) for new arrays. A
    controlled A/B comparison (see the project's paper Discussion and the
    module docstring above) found this full-tile configuration consistently
    matches or beats the legacy downsample-then-residual-DWT pipeline built
    by :func:`lookup_preset`, on both compressed size and reconstruction
    RMSE, at every matched low-frequency resolution tested -- while also
    being architecturally simpler (one codec, one sidecar pair, no separate
    low-frequency sub-pipeline).

    .. important::

       As with :func:`lookup_preset`, store the explicit parameter dict
       (from :func:`hcc.compression.grad_dwt.lookup_preset` with
       ``codec=False``) in your array metadata rather than the preset name,
       so future preset-default changes cannot make existing data
       unreadable.

    Parameters
    ----------
    name : {'low', 'medium', 'high', 'ultra'}
        Quality preset name -- the *standalone* preset table (see
        :func:`hcc.compression.grad_dwt.lookup_preset`, ``codec=False``),
        which has no ``downsample_factor`` (not applicable here: there is
        no separate downsampling stage).

    Returns
    -------
    list[Codec]
        ``[HccGradDwtCodec(...), BytesCodec(), GzipCodec(level=...)]`` --
        pass directly as the ``codecs=`` argument when creating a Zarr
        array.

    Raises
    ------
    ValueError
        If *name* is not a known preset.

    See Also
    --------
    lookup_preset : Legacy HccDecomposeCodec-wrapped (downsample + residual
        DWT) pipeline factory, kept for backward compatibility.
    hcc.compression.grad_dwt.lookup_preset : Returns the raw parameter dict.

    Examples
    --------
    >>> codecs = build_fulltile_pipeline("medium")
    >>> codecs[0].wavelet
    'bior4.4'
    >>> codecs[0].levels
    5
    """
    from zarr.codecs import BytesCodec, GzipCodec

    from hcc.compression.grad_dwt import lookup_preset as lookup_raw_preset

    presets = lookup_raw_preset(name, codec=False)
    level = presets.pop("zlib_level")

    return [
        HccGradDwtCodec(**presets),
        BytesCodec(),
        GzipCodec(level=level),
    ]
