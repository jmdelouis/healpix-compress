import math
from typing import TypeVar

import numpy as np
import numpy.typing as npt

T = TypeVar("T", bound=np.floating)


def _decode_zorder(delta: int) -> tuple[npt.NDArray[np.uint64], npt.NDArray[np.uint64]]:
    size = 4**delta
    x = np.zeros((size,), dtype="uint64")
    y = np.zeros((size,), dtype="uint64")

    code = np.arange(size, dtype="uint64")
    for bit in range(delta):
        x |= ((code >> (bit << 1)) & 1) << bit
        y |= ((code >> ((bit << 1) + 1)) & 1) << bit

    return x, y


def encode_zorder(
    x: npt.NDArray[np.integer], y: npt.NDArray[np.integer], delta: int
) -> npt.NDArray[np.uint64]:
    """Exact inverse of :func:`_decode_zorder`: interleave (x, y) integer
    coordinates into a Morton/Z-order code.

    Unlike :func:`_decode_zorder` (which always decodes the *canonical*
    range ``0 .. 4**delta - 1``), this accepts arbitrary integer
    coordinate arrays -- in particular coordinates that have been shifted
    (e.g. with periodic wraparound via ``% side``) to look up a
    neighbouring pixel's NESTED index. This is the building block used by
    :mod:`hcc.compression.quadtree_wavelet` to do neighbour lookups
    directly in NESTED-index space, without ever forming a 2-D raster.

    Bit convention (matching :func:`_decode_zorder` exactly): bit ``b`` of
    ``x`` becomes bit ``2*b`` of the code, bit ``b`` of ``y`` becomes bit
    ``2*b + 1``.

    Parameters
    ----------
    x, y : numpy.ndarray of integer dtype, broadcastable shapes
        Coordinates in ``[0, 2**delta)`` (values outside this range should
        be reduced first, e.g. via periodic wraparound).
    delta : int
        Number of bits per coordinate (i.e. side length is ``2**delta``).

    Returns
    -------
    numpy.ndarray of dtype uint64
        NESTED (Morton) indices, same broadcast shape as ``x``/``y``.

    See Also
    --------
    _decode_zorder : The inverse operation, restricted to the canonical
        index range.
    """
    x = np.asarray(x, dtype="uint64")
    y = np.asarray(y, dtype="uint64")
    code = np.zeros(np.broadcast_shapes(x.shape, y.shape), dtype="uint64")
    for bit in range(delta):
        code |= ((x >> np.uint64(bit)) & np.uint64(1)) << np.uint64(bit << 1)
        code |= ((y >> np.uint64(bit)) & np.uint64(1)) << np.uint64((bit << 1) + 1)
    return code


def flatten_tile(data: npt.NDArray[T]) -> npt.NDArray[T]:
    delta = int(np.log(data.size) / np.log(4))
    x, y = _decode_zorder(delta)

    return np.squeeze(data[..., x, y])


def expand_tile(data: npt.NDArray[T]) -> npt.NDArray[T]:
    nside = int(math.sqrt(data.shape[-1]))
    delta = int(math.log2(nside))

    x, y = _decode_zorder(delta)

    shape = (*data.shape[:-1], nside, nside)

    # TODO: figure out how to invert x, y to become a 2D array of indexes into data
    out = np.empty(shape, dtype=data.dtype)
    out[..., x, y] = data

    return out
