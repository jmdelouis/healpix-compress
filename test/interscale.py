"""Does the link between scales pay? (idea 2)

Equal-area cells make the HEALPix quadtree an exact Haar pyramid. A parent's
value *is* the mean of its four children, so the transform

    a = (c0+c1+c2+c3)/2,   h = (c0-c1+c2-c3)/2,
    v = (c0+c1-c2-c3)/2,   d = (c0-c1-c2+c3)/2

is orthonormal, critically sampled, exactly invertible, and free of the two
things that make a planar wavelet awkward on the sphere: there is no filter to
design and no tile boundary to extend. Applied recursively it turns a level-20
tile into an approximation at a coarse level plus three detail subbands per
level, and every approximation channel is itself a *bit-exact* HEALPix product
at its own level rather than a decoded approximation of one.

`level_choice.py` already measured what those subbands carry. Per tile, the
median shares of the variance are

    level 20 detail   3.07%   in 75.00% of the coefficients
    level 19 detail   7.15%   in 18.75%
    level 18 detail   9.77%   in  4.69%
    level 17 approx  79.60%   in  1.56%

-- about a factor of nine in variance per coefficient at every step, which is
the decay a subband coder lives on. That much is settled and is not the
question. The question is the one the design rests on: the structures are
*linked across scales*, so a coefficient should cost less once its parent is
known. This measures whether it does.

The measurement is a conditional entropy on quantised coefficients under four
context models:

    H(q)                 no context; what an order-0 coder pays
    H(q | parent)        the inter-scale link alone -- the claim under test
    H(q | neighbour)     the intra-band link alone -- what EBCOT already uses
    H(q | parent, nbr)   both

The third is the control, and it is what lets this test fail. JPEG 2000 already
exploits intra-band context; if the parent merely reproduces what the neighbour
gives and adds nothing on top of it, then the link across scales carries no
information deployed coders are not already charging for, and the design has no
edge to trade for Haar's poor compaction.

Rates are ideal entropy-coder bounds -- no coding overhead, no signalling,
empirical per-tile distributions. That flatters this coder against the measured
bytes of real codecs on the same tiles, deliberately: a bound that loses is a
dead idea, and a bound that wins is necessary but not sufficient.

Two transforms, so that the context result is not confounded with the filter
-------------------------------------------------------------------------------
`--transform haar` is the exact quadtree Haar above. `--transform cdf53` is the
quadtree-native CDF 5/3 of `hcc.compression.quadtree_wavelet`, computed by
lifting with neighbour lookups in NESTED index space. The distinction matters
because `expand_tile` is a Morton de-interleave: a level-10 parent's
descendants *are* a square raster, so JPEG 2000 already runs a 9/7 wavelet on
this very quadtree. Haar's one vanishing moment against its four is a
compaction deficit no context model can repay, and running both separates "the
link between scales pays" from "Haar is the wrong filter".

The CDF 5/3 lifting is not normalised, so each subband is quantised with
step / gain, the gain measured by pushing a unit impulse through the inverse.
That is JPEG 2000's own L2-norm weighting, and it is what makes the two
transforms comparable at one allocation instead of two.

Idea 3, on the same coefficients
--------------------------------
The same run prices the other half of the design. Section V-F says the reader
dominates the error, which is another way of saying its kernel is large: the
carrier is oversampled relative to the operator that reads it, so a gather read
never touches most cells and anything a codec spends on them is spent on
nothing. What the coder should minimise is ||R(x_hat - x)||, not ||x_hat - x||.

Two different things can make the product see less than the codec is charged
for, and they call for opposite designs, so the measurement separates them:

    transmission ||Re|| / ||e||    near 1 with a large kernel means the error
                                   per read pixel is unchanged and the gain is
                                   in *allocation*: rate spent on unreachable
                                   cells.
                                   well below 1 means the reader *attenuates* a
                                   white coding error, and a coarser step
                                   reaches the same product distortion.

Both are measured against the two readers that matter: the nearest read the
paper's floor uses, and the operator read of the companion. The gain belongs to
the reader, not to HEALPix, and improving the reader removes it -- which is why
this is reported next to the rate, not instead of it.

Usage
-----
    python interscale.py --limit 4
    python interscale.py
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import build_tiles                        # noqa: E402
import common_reference                 # noqa: E402
import region_sites                       # noqa: E402

DEPTH_MAX = 6                  # detail at levels 20..15, approximation at 14
TRANSFORM = "haar"             # or "cdf53": the quadtree-native CDF 5/3
N_CTX = 6                      # magnitude buckets in a context symbol
BENCH = "results/hpx_codec_bench.csv"
ORIENT = ("h", "v", "d")
MODELS = (("order0", "no context"),
          ("parent", "parent (inter-scale)"),
          ("neighbour", "neighbour (intra-band)"),
          ("both", "parent + neighbour"))
READERS = ("nearest", "psf")   # the paper's floor reader, and the
                               # companion's forward operator
CLASS_COLOR = {"urban": "#D55E00", "water": "#0072B2",
               "forest": "#009E73", "agriculture": "#E69F00"}


# --------------------------------------------------------------------------
# the transform
# --------------------------------------------------------------------------
def haar_step(a):
    """One level of the exact HEALPix Haar transform.

    `a` is a contiguous NESTED block, so `a.reshape(-1, 4)` are sibling quads.
    Orthonormal, so a coefficient error of e contributes e to the RMSE of the
    reconstructed cells and rate--distortion can be read off directly.
    """
    q = np.asarray(a, dtype=np.float64).reshape(-1, 4)
    c0, c1, c2, c3 = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    return ((c0 + c1 + c2 + c3) / 2.0,
            {"h": (c0 - c1 + c2 - c3) / 2.0,
             "v": (c0 + c1 - c2 - c3) / 2.0,
             "d": (c0 - c1 - c2 + c3) / 2.0})


def haar_inverse_step(a, det):
    h, v, d = det["h"], det["v"], det["d"]
    out = np.empty((a.size, 4), dtype=np.float64)
    out[:, 0] = (a + h + v + d) / 2.0
    out[:, 1] = (a - h + v - d) / 2.0
    out[:, 2] = (a + h - v - d) / 2.0
    out[:, 3] = (a - h - v + d) / 2.0
    return out.ravel()


def pyramid(x, depth_max=DEPTH_MAX):
    """Return (approximation, [details per level, finest first])."""
    a = np.asarray(x, dtype=np.float64).ravel()
    levels = []
    for _ in range(depth_max):
        a, det = haar_step(a)
        levels.append(det)
    return a, levels


def rebuild(a, levels):
    for det in reversed(levels):
        a = haar_inverse_step(a, det)
    return a


class Pyr:
    """A pyramid on the HEALPix quadtree, with the gains its inverse applies.

    The Haar pyramid is orthonormal, so a coefficient error of e becomes an
    error of e in the tile and every subband takes the same quantiser step.
    The CDF 5/3 lifting is *not* normalised: a unit coefficient in subband s
    reconstructs to a tile perturbation of norm `gain[s]`, which varies here
    by a factor of seven between subbands. Quantising it with a single step
    would therefore not be the uniform allocation it is for Haar, and any
    comparison between the two transforms would be a comparison of
    allocations instead. Each subband is given step / gain[s], which is the
    L2-norm weighting JPEG 2000 itself uses, so the two transforms are
    compared at the same distortion allocation and only the transform
    differs.
    """

    def __init__(self, approx, levels, gains, inverse, name, source):
        self.source = np.asarray(source, dtype=np.float64).ravel()
        self.approx = approx
        self.levels = levels
        self.gains = gains
        self.inverse = inverse
        self.name = name
        self.n_leaf = approx.size * (4 ** len(levels))

    def step_of(self, key, step):
        return step / self.gains[key]


def _cdf53_gains(n, depth, qw):
    """Synthesis gain of each subband, by pushing a unit impulse through."""
    zero = np.zeros(n, dtype=np.float64)
    ll, det = qw.qt_dwt53_forward(zero, levels=depth)
    ll = np.asarray(ll, dtype=np.float64)
    det = [tuple(np.asarray(a, dtype=np.float64) for a in d) for d in det]

    def norm(fill):
        a = np.zeros_like(ll)
        d = [tuple(np.zeros_like(v) for v in dd) for dd in det]
        fill(a, d)
        return float(np.linalg.norm(np.asarray(qw.qt_dwt53_inverse(a, d))))

    gains = {"a": norm(lambda a, d: a.__setitem__(a.size // 2, 1.0))}
    n_lv = len(det)
    for j, o in enumerate(ORIENT):
        for c in range(n_lv):                    # c: coarsest-first index
            k = n_lv - 1 - c                     # k: finest-first index
            gains[(k, o)] = norm(
                lambda a, d, c=c, j=j: d[c][j].__setitem__(d[c][j].size // 2, 1.0))
    return gains


_GAIN_CACHE = {}


def build_pyramid(x, depth_max=DEPTH_MAX, transform=TRANSFORM):
    """The pyramid this run codes, with its per-subband synthesis gains."""
    x = np.asarray(x, dtype=np.float64).ravel()
    if transform == "haar":
        a, levels = pyramid(x, depth_max)
        gains = {"a": 1.0}
        gains.update({(k, o): 1.0 for k in range(len(levels)) for o in ORIENT})
        return Pyr(a, levels, gains, rebuild, "haar", x)
    if transform != "cdf53":
        raise ValueError(f"unknown transform {transform!r}")

    from hcc.compression import quadtree_wavelet as qw
    ll, det = qw.qt_dwt53_forward(x, levels=depth_max)
    ll = np.asarray(ll, dtype=np.float64)
    # qt_dwt53 returns details coarsest-first; this module works finest-first
    levels = [{o: np.asarray(d[j], dtype=np.float64)
               for j, o in enumerate(ORIENT)} for d in reversed(det)]

    def inverse(a, lv):
        back = [tuple(lv[k][o] for o in ORIENT) for k in range(len(lv) - 1, -1, -1)]
        return np.asarray(qw.qt_dwt53_inverse(a, back), dtype=np.float64).ravel()

    key = (x.size, depth_max)
    if key not in _GAIN_CACHE:
        _GAIN_CACHE[key] = _cdf53_gains(x.size, depth_max, qw)
    return Pyr(ll, levels, _GAIN_CACHE[key], inverse, "cdf53", x)


# --------------------------------------------------------------------------
# quantisation, rate, distortion
# --------------------------------------------------------------------------
def quantise(c, step):
    """Uniform deadzone quantiser; returns (indices, reconstruction)."""
    q = np.floor(np.abs(c) / step) * np.sign(c)
    r = np.where(q == 0, 0.0, np.sign(q) * (np.abs(q) + 0.5) * step)
    return q.astype(np.int64), r


def rmse_at(pyr, step):
    """Distortion only -- cheap, so the step can be bisected.

    Exact for an orthonormal transform. For a biorthogonal one the sum of
    weighted coefficient errors is the standard L2-norm estimate, and the
    published RMSE is always the one measured after an actual inverse in
    `coder_point`, never this.
    """
    _, r = quantise(pyr.approx, pyr.step_of("a", step))
    err2 = float((((pyr.approx - r) * pyr.gains["a"]) ** 2).sum())
    for k, det in enumerate(pyr.levels):
        for o in ORIENT:
            g = pyr.gains[(k, o)]
            _, rr = quantise(det[o], pyr.step_of((k, o), step))
            err2 += float((((det[o] - rr) * g) ** 2).sum())
    return float(np.sqrt(err2 / pyr.n_leaf))


def rmse_exact(pyr, step):
    """Distortion measured by inverting the quantised pyramid."""
    e = np.asarray(reconstruct(pyr, step), dtype=np.float64).ravel() - pyr.source
    return float(np.sqrt(np.mean(e ** 2)))


def solve_step(pyr, target_rmse, sd, iters=30):
    """Bisect the quantiser step onto a given distortion.

    The bisection is on the *exact* distortion, not on the L2-norm estimate.
    For an orthonormal transform the two agree; for a biorthogonal one the
    estimate lands several percent off target, and since the whole point of
    running two transforms is to compare their rates at the same distortion,
    an off-target biorthogonal point would credit it with a saving it did not
    earn.
    """
    lo, hi = 1e-8 * sd, 8.0 * sd
    for _ in range(iters):
        mid = float(np.sqrt(lo * hi))
        if rmse_exact(pyr, mid) > target_rmse:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.0005:
            break
    return float(np.sqrt(lo * hi))


def _entropy(counts):
    n = float(counts.sum())
    if n <= 0:
        return 0.0
    p = counts[counts > 0] / n
    return float(-(p * np.log2(p)).sum())


def entropy0(sym):
    _, counts = np.unique(np.asarray(sym).ravel(), return_counts=True)
    return _entropy(counts.astype(np.float64))


def entropy_cond(sym, ctx):
    """H(sym | ctx) in bits per symbol, empirical over one tile."""
    sym = np.asarray(sym).ravel()
    ctx = np.asarray(ctx).ravel()
    order = np.argsort(ctx, kind="stable")
    s, c = sym[order], ctx[order]
    edges = np.flatnonzero(np.diff(c)) + 1
    total = 0.0
    for part in np.split(s, edges):
        _, counts = np.unique(part, return_counts=True)
        total += part.size * _entropy(counts.astype(np.float64))
    return total / s.size if s.size else 0.0


def bucket(q, n_ctx=N_CTX):
    """Magnitude class of a quantised coefficient, on the usual log ladder."""
    m = np.abs(np.asarray(q, dtype=np.int64))
    return np.minimum(np.floor(np.log2(m + 1)).astype(np.int64), n_ctx - 1)


def neighbour_bucket(q, n_ctx=N_CTX):
    """Class of the preceding coefficient of the same subband.

    Coefficients are in NESTED order, so index i-1 is a spatially adjacent cell
    for three indices out of four. This deliberately understates the intra-band
    context a real coder builds from a full causal neighbourhood: the control
    is meant to be conservative, i.e. to favour the hypothesis under test.
    """
    q = np.asarray(q, dtype=np.int64)
    prev = np.empty_like(q)
    prev[0] = 0
    prev[1:] = q[:-1]
    return bucket(prev, n_ctx)


def coder_point(pyr, step):
    """Quantise every subband; return (rmse, bits per model).

    The RMSE is measured by actually inverting the quantised pyramid, so it
    is exact for a biorthogonal transform too rather than an L2-norm estimate.
    """
    qa, ra = quantise(pyr.approx, pyr.step_of("a", step))
    dec = []
    h_a = pyr.approx.size * entropy0(qa)      # no parent above the top
    bits = {k: h_a for k, _ in MODELS}
    prev = {o: None for o in ORIENT}
    per_band = []
    levels = pyr.levels
    rec = [dict() for _ in levels]
    for k in range(len(levels) - 1, -1, -1):  # coarsest detail first
        for o in ORIENT:
            c = levels[k][o]
            q, r = quantise(c, pyr.step_of((k, o), step))
            rec[k][o] = r
            nb = neighbour_bucket(q)
            pb = (np.zeros(q.size, dtype=np.int64) if prev[o] is None
                  else np.repeat(bucket(prev[o]), 4))
            b = {"order0": entropy0(q),
                 "parent": entropy_cond(q, pb),
                 "neighbour": entropy_cond(q, nb),
                 "both": entropy_cond(q, pb * N_CTX + nb)}
            for m in b:
                bits[m] += q.size * b[m]
            per_band.append(dict(level=cfg.LEVEL - k, orient=o, n_coef=q.size,
                                 nonzero=float((q != 0).mean()),
                                 **{f"bpc_{m}": v for m, v in b.items()}))
            prev[o] = q
    # exact distortion: invert the quantised pyramid rather than model it
    dec = pyr.inverse(ra, rec)
    err = np.asarray(dec, dtype=np.float64).ravel() - pyr.source
    return (float(np.sqrt(np.mean(err ** 2))), bits, pyr.n_leaf, per_band)


# --------------------------------------------------------------------------
# idea 3: what does the reader never see?
# --------------------------------------------------------------------------
def reconstruct(pyr, step):
    """The decoded tile at one quantiser step."""
    _, ra = quantise(pyr.approx, pyr.step_of("a", step))
    det = [{o: quantise(lv[o], pyr.step_of((k, o), step))[1] for o in ORIENT}
           for k, lv in enumerate(pyr.levels)]
    return pyr.inverse(ra, det)


def rms_masked(field, mask):
    v = np.asarray(field, dtype=np.float64)[np.asarray(mask)]
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v ** 2))) if v.size else np.nan


def transmission(cr, e_cells):
    """||R e|| on the mask over ||e|| in cell space.

    Every reader the paper defines is linear, so R(x_hat - x) is exactly the
    product-level error a coding error e produces. Two very different things
    can make this small, and telling them apart is the point:

    * a *kernel*  -- cells the reader never touches. A gather read has one,
      and it is large precisely because the carrier is oversampled. It does
      not shrink the error per read pixel, so transmission stays near 1; what
      it wastes is rate, spent controlling cells nothing depends on.
    * *attenuation* -- a smoothing read suppresses the high-frequency part of
      a white coding error, so transmission falls well below 1 and the coder
      may quantise coarser for the same product distortion.
    """
    num = rms_masked(cr.predict(e_cells), cr.mask)
    den = float(np.sqrt(np.mean(np.asarray(e_cells, dtype=np.float64) ** 2)))
    return num / den if den > 0 else np.nan


def kernel_energy(cr, e_cells):
    """(share of the error energy in ker(R), share of cells R can reach)."""
    sup = cr.read_support()
    e2 = np.asarray(e_cells, dtype=np.float64).ravel() ** 2
    tot = float(e2.sum())
    return (float(e2[~sup].sum()) / tot if tot > 0 else np.nan,
            float(sup.mean()))


def solve_step_product(cr, pyr, x, target_rmse, seed, iters=20):
    """Bisect the step so the *product-level* error hits the target.

    The bits saved between this step and the cell-matched one are what a coder
    would gain from being told the target is a product target -- available
    with no new codec, only a corrected distortion norm. Any further gain
    needs a coder that actively steers its error, which this does not measure.
    """
    lo, hi = seed / 64.0, seed * 64.0
    for _ in range(iters):
        mid = float(np.sqrt(lo * hi))
        e = reconstruct(pyr, mid) - x
        if rms_masked(cr.predict(e), cr.mask) > target_rmse:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.002:
            break
    return float(np.sqrt(lo * hi))


def zerotree_stats(pyr, step):
    """P(child = 0 | parent = 0) against P(child = 0): the zerotree premise."""
    levels = pyr.levels
    out = []
    for k in range(len(levels) - 1, 0, -1):
        for o in ORIENT:
            qp, _ = quantise(levels[k][o], pyr.step_of((k, o), step))
            qc, _ = quantise(levels[k - 1][o], pyr.step_of((k - 1, o), step))
            pz = np.repeat(qp == 0, 4)
            if pz.all() or not pz.any():
                continue
            out.append(dict(level_parent=cfg.LEVEL - k, orient=o,
                            p_child_zero=float((qc == 0).mean()),
                            p_zero_given_parent_zero=float((qc[pz] == 0).mean()),
                            p_zero_given_parent_nonzero=float(
                                (qc[~pz] == 0).mean())))
    return out


# --------------------------------------------------------------------------
def region_rows(scene, front, verbose=True, transform=TRANSFORM,
                kernel=None):
    kernel = kernel or build_tiles.DEFAULT_KERNEL
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene, kernel)
    x = np.asarray(cell_data, dtype=np.float64).ravel()
    if x.size != cfg.N_CHILDREN:
        raise ValueError(f"{x.size} cells, expected {cfg.N_CHILDREN}")
    sd = float(np.std(x))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("degenerate tile")
    pyr = build_pyramid(x, DEPTH_MAX, transform)

    # The transform must be exact before anything measured from it means
    # anything.
    rt = float(np.max(np.abs(np.asarray(
        pyr.inverse(pyr.approx, pyr.levels)).ravel() - x)))
    scale = max(1.0, float(np.max(np.abs(x))))
    if rt > 1e-9 * scale:
        raise RuntimeError(f"Haar round-trip error {rt:.3e}")

    cls = scene.split("__")[1]
    energies = [dict(scene=scene, cls=cls, level=cfg.LEVEL - k, orient=o,
                     transform=pyr.name, n_coef=pyr.levels[k][o].size,
                     var_frac=float(((pyr.levels[k][o] * pyr.gains[(k, o)]) ** 2
                                     ).sum() / (x.size * sd ** 2)))
                for k in range(len(pyr.levels)) for o in ORIENT]
    energies.append(dict(scene=scene, cls=cls, level=cfg.LEVEL - DEPTH_MAX,
                         orient="a", transform=pyr.name, n_coef=pyr.approx.size,
                         var_frac=float((((pyr.approx - pyr.approx.mean())
                                          * pyr.gains["a"]) ** 2).sum()
                                        / (x.size * sd ** 2))))

    obs = np.asarray(utm_crop, dtype=np.float64)
    crs, failures = {}, []
    for name in READERS:
        try:
            cr = common_reference.build_for_tile(npz, mode=name)
            # A misbuilt operator is not detectable from a transmission ratio,
            # so each reader publishes its own representation floor next to
            # the ratio: if this is worse than the nearest read's, the reader
            # is wrong and idea 3's numbers for it mean nothing.
            floor = cr.rmse_on_mask(cr.predict(x), obs)
            if not np.isfinite(floor):
                raise RuntimeError("non-finite representation floor")
            crs[name] = (cr, float(floor))
        except Exception as exc:            # noqa: BLE001
            print(f"    reader {name!r} FAILED: {type(exc).__name__}: {exc}")
            failures.append(dict(scene=scene, cls=scene.split("__")[1],
                                 reader=name, error=f"{type(exc).__name__}: {exc}"))

    targets = front.get(scene) or [(t, np.nan) for t in
                                   np.geomspace(1e-4 * sd, 0.5 * sd, 6)]
    rows, ztree, bands, reads = [], [], [], list(failures)
    for target_rmse, ref_bytes in targets:
        step = solve_step(pyr, target_rmse, sd)
        rmse, bits, n_leaf, per_band = coder_point(pyr, step)
        row = dict(scene=scene, cls=cls, transform=pyr.name, kernel=kernel,
                   step=step, target_rmse=target_rmse,
                   rmse=rmse, tile_std=sd, n_coef=n_leaf, ref_bytes=ref_bytes)
        for m, _ in MODELS:
            row[f"bits_{m}"] = bits[m]
            row[f"bpc_{m}"] = bits[m] / n_leaf
            row[f"bytes_{m}"] = bits[m] / 8.0
        b0 = bits["order0"]
        row["gain_parent"] = 1.0 - bits["parent"] / b0
        row["gain_neighbour"] = 1.0 - bits["neighbour"] / b0
        row["gain_both"] = 1.0 - bits["both"] / b0
        # The number the design lives on: what the parent adds once a
        # within-subband context is already being used.
        row["gain_parent_over_neighbour"] = (
            1.0 - bits["both"] / bits["neighbour"])
        # --- idea 3: the reader's kernel and its attenuation -----------
        e = reconstruct(pyr, step) - x
        for name, (cr, floor) in crs.items():
            t = transmission(cr, e)
            ker, sup = kernel_energy(cr, e)
            rec = dict(scene=scene, cls=cls, transform=pyr.name, kernel=kernel,
                       reader=name,
                       target_rmse=target_rmse, rmse_cell=rmse,
                       rmse_product=t * rmse if t == t else np.nan,
                       transmission=t, kernel_energy_frac=ker,
                       support_fraction=sup, reader_floor=floor,
                       bits_cell=bits["both"], error="")
            if t == t and t > 0 and t < 0.999:
                s2 = solve_step_product(cr, pyr, x, target_rmse,
                                        seed=step / t)
                _, b2, _, _ = coder_point(pyr, s2)
                rec["step_product"] = s2
                rec["bits_product"] = b2["both"]
                rec["gain_product_target"] = 1.0 - b2["both"] / bits["both"]
            reads.append(rec)

        rows.append(row)
        ztree += [dict(scene=scene, target_rmse=target_rmse, **z)
                  for z in zerotree_stats(pyr, step)]
        bands += [dict(scene=scene, cls=cls, target_rmse=target_rmse, **b)
                  for b in per_band]
        if verbose:
            tail = (f"  | bound/codec {row['bytes_both']/ref_bytes:.2f}x"
                    if ref_bytes == ref_bytes else "")
            print(f"    rmse {rmse:.3e}  order0 {row['bpc_order0']:.3f}"
                  f"  parent {row['bpc_parent']:.3f}"
                  f"  nbr {row['bpc_neighbour']:.3f}"
                  f"  both {row['bpc_both']:.3f} b/coef"
                  f"  (+{100*row['gain_parent_over_neighbour']:.1f}% over nbr)"
                  f"{tail}")
    return rows, energies, ztree, bands, reads


def pareto_front(path=BENCH, n_points=6):
    """Per scene, the cheapest measured bytes reaching each distortion."""
    try:
        d = pd.read_csv(path)
    except FileNotFoundError:
        print(f"  (no {path}; falling back to a geometric step grid)")
        return {}
    d = d.dropna(subset=["rmse", "total_bytes"])
    d = d[(d.rmse > 0) & (d.total_bytes > 0)]
    front = {}
    for scene, g in d.groupby("scene"):
        g = g.sort_values("rmse")
        best, keep = np.inf, []
        for r, b in zip(g.rmse.values, g.total_bytes.values):
            if b < best:                     # not dominated by a tighter point
                best = b
                keep.append((float(r), float(b)))
        if len(keep) > n_points:
            idx = np.linspace(0, len(keep) - 1, n_points).round().astype(int)
            keep = [keep[i] for i in sorted(set(idx.tolist()))]
        front[scene] = keep
    return front


# --------------------------------------------------------------------------
def summarize(df, en, zt, rd=None):
    n = en.scene.nunique()
    out = ["# Does the link between scales pay?", "",
           f"- regions: {n}; rates are ideal entropy bounds, so they flatter "
           f"this coder against the measured bytes of real codecs",
           f"- transform {df['transform'].iloc[0] if 'transform' in df else 'haar'}, "
           f"tiles {df['kernel'].iloc[0] if 'kernel' in df else 'clough_tocher'}",
           ""]

    tr = df["transform"].iloc[0] if "transform" in df else "haar"
    name = {"haar": "exact HEALPix Haar", "cdf53": "quadtree CDF 5/3"}.get(tr, tr)
    out += [f"## Subband energy of the {name} pyramid", "",
            "| level | coefficients | % of cells | % of variance |",
            "|---|---|---|---|"]
    g = en.groupby(["level", "orient"]).agg(nc=("n_coef", "mean"),
                                            v=("var_frac", "median"))
    g = g.groupby("level").agg(nc=("nc", "sum"), v=("v", "sum"))
    for lvl, r in g.sort_index(ascending=False).iterrows():
        out.append(f"| {lvl} | {r.nc:.0f} | "
                   f"{100*r.nc/cfg.N_CHILDREN:.2f}% | {100*r.v:.2f}% |")
    out.append("")

    out += ["## Bits per coefficient by context model", "",
            "| model | b/coef (median) | saving vs order-0 |", "|---|---|---|"]
    for m, lab in MODELS:
        gain = 0.0 if m == "order0" else df[f"gain_{m}"].median()
        out.append(f"| {lab} | {df[f'bpc_{m}'].median():.3f} | "
                   f"{100*gain:.1f}% |")
    out.append("")

    gp = float(df.gain_parent.median())
    gn = float(df.gain_neighbour.median())
    gpn = float(df.gain_parent_over_neighbour.median())
    out += [f"- the inter-scale link alone saves {100*gp:.1f}% of the order-0 "
            f"rate; the intra-band link alone saves {100*gn:.1f}%",
            f"- **the parent adds {100*gpn:.1f}% on top of the intra-band "
            f"context**, and that is the number the design lives on, because "
            f"deployed coders already charge for the intra-band link", ""]
    if gpn < 0.03:
        out.append("- On this evidence the link across scales carries almost "
                   "nothing a within-subband context does not already carry. "
                   "A HEALPix Haar coder would have to win on the exact "
                   "multiresolution property alone, against a transform with "
                   "one vanishing moment -- a hard trade, and the honest "
                   "conclusion is that the gain claimed for coding the "
                   "connection between levels is not there.")
    else:
        out.append(f"- A context-modelled inter-scale coder has {100*gpn:.1f}% "
                   f"of rate to trade against Haar's compaction deficit "
                   f"relative to a 9/7 wavelet. The head-to-head below says "
                   f"whether that is enough.")
    out.append("")

    z = zt.groupby("orient")[["p_child_zero", "p_zero_given_parent_zero",
                              "p_zero_given_parent_nonzero"]].median()
    out += ["## The zerotree premise", "",
            "| orientation | P(child = 0) | given parent = 0 | "
            "given parent != 0 |", "|---|---|---|---|"]
    for o, r in z.iterrows():
        out.append(f"| {o} | {r.p_child_zero:.3f} | "
                   f"{r.p_zero_given_parent_zero:.3f} | "
                   f"{r.p_zero_given_parent_nonzero:.3f} |")
    lift = float((z.p_zero_given_parent_zero - z.p_child_zero).median())
    out += ["", f"- a zerotree is worth signalling only if an insignificant "
                f"parent really predicts insignificant children; here it lifts "
                f"P(child = 0) by {100*lift:.1f} points", ""]

    q = pd.qcut(np.log10(df.rmse), 4,
                labels=["tightest", "tight", "loose", "loosest"])
    out += ["## By distortion regime", "",
            "| quartile | RMSE | parent adds over neighbour | "
            "bound / codec bytes |", "|---|---|---|---|"]
    tmp = df.assign(_q=q, _r=df.bytes_both / df.ref_bytes)
    for lab, g in tmp.groupby("_q", observed=True):
        out.append(f"| {lab} | {g.rmse.median():.2e} | "
                   f"{100*g.gain_parent_over_neighbour.median():.1f}% | "
                   f"{g._r.median():.2f} |")
    out += ["", "- the parent earns its keep where most coefficients are "
                "zero, which is the low-rate end; that is the regime a "
                "zerotree coder is built for", ""]

    h = df.dropna(subset=["ref_bytes"])
    if len(h):
        out += ["## Head to head against the measured codec front", "",
                "| model | entropy bytes / codec bytes (median) |", "|---|---|"]
        for m, lab in MODELS:
            out.append(f"| {lab} | {(h[f'bytes_{m}']/h.ref_bytes).median():.2f} |")
        best = float((h.bytes_both / h.ref_bytes).median())
        out.append("")
        if best >= 1.0:
            out.append(f"- At equal distortion the pyramid's *entropy bound* "
                       f"is {best:.2f} times the bytes real codecs actually "
                       f"spend on the same tiles. A bound above one is "
                       f"decisive against the design: no implementation can "
                       f"reach its own bound, and this one is already losing.")
        else:
            out.append(f"- At equal distortion the pyramid's entropy bound is "
                       f"{best:.2f} of the bytes real codecs spend on the same "
                       f"tiles. That is headroom, and it is necessary but not "
                       f"sufficient -- the bound charges nothing for coding "
                       f"overhead, signalling, or the loss from a per-tile "
                       f"adaptive model.")
    if rd is not None and len(rd):
        if "error" in rd:
            # an empty marker round-trips through CSV as NaN, so a bare
            # .eq("") would silently drop every good row
            rd = rd[rd["error"].fillna("").astype(str).str.len() == 0]
        out += ["", "## Idea 3: what the reader never sees", "",
                "| reader | floor | transmission ||Re||/||e|| | error energy in "
                "ker(R) | cells R can reach | rate saved by a product "
                "target |", "|---|---|---|---|---|---|"]
        for name, g in rd.groupby("reader"):
            gain = (g.gain_product_target.median()
                    if "gain_product_target" in g else float("nan"))
            gtxt = "---" if gain != gain else f"{100*gain:.1f}%"
            out.append(f"| {name} | {g.reader_floor.median():.3e} | "
                       f"{g.transmission.median():.3f} | "
                       f"{100*g.kernel_energy_frac.median():.1f}% | "
                       f"{100*g.support_fraction.median():.1f}% | {gtxt} |")
        out.append("")
        # Column-wise medians hide a paired comparison. A reader whose own
        # floor is worse than the nearest read's is reading a carrier it was
        # not matched to, and its transmission then prices noise against the
        # wrong estimand.
        piv = rd.pivot_table(index="scene", columns="reader",
                             values="reader_floor", aggfunc="median")
        if "nearest" in piv:
            for name in [c for c in piv.columns if c != "nearest"]:
                ratio = (piv[name] / piv["nearest"]).dropna()
                if not len(ratio):
                    continue
                worse = int((ratio > 1).sum())
                out.append(
                    f"- paired floors, {name} over nearest: median "
                    f"{ratio.median():.3f}, worse on {worse}/{len(ratio)} "
                    f"regions" + ("  **-- this reader is not matched to the "
                                  "carrier these tiles were built with, so its "
                                  "numbers below price noise against the wrong "
                                  "estimand**" if ratio.median() > 1.02 else ""))
            out.append("")
        for name, g in rd.groupby("reader"):
            t = float(g.transmission.median())
            ker = float(g.kernel_energy_frac.median())
            sup = float(g.support_fraction.median())
            if t > 0.9 and ker > 0.2:
                out.append(
                    f"- **{name}**: transmission {t:.2f} with {100*ker:.0f}% of "
                    f"the error energy in the kernel. The kernel does not make "
                    f"the product error smaller -- it reveals wasted rate. The "
                    f"product depends on {100*sup:.1f}% of the cells, so at "
                    f"high rate a reader-aware coder has a bound of about "
                    f"{1/sup:.1f}x on rate, and it is a bound on *allocation*, "
                    f"not something a coarser step buys.")
            elif t <= 0.9:
                out.append(
                    f"- **{name}**: transmission {t:.2f}. This reader "
                    f"*attenuates* a white coding error rather than ignoring "
                    f"particular cells, so the same product distortion is "
                    f"reached with a coarser step -- a saving that needs no "
                    f"new codec, only the target expressed on the product.")
            else:
                out.append(f"- **{name}**: transmission {t:.2f}, kernel "
                           f"{100*ker:.1f}%: nothing free here.")
        out += ["", "- The gain is a property of the *reader*, not of HEALPix. "
                    "A gather read has a large kernel exactly because the "
                    "carrier is oversampled relative to it, which is what "
                    "Section V-F calls the dominant error; fixing the reader "
                    "closes the kernel and removes this gain with it. Shaping "
                    "coding noise into a smoothing reader's stopband is the "
                    "safer half, and it is safe only for a product consumed "
                    "through that operator: the companion's point is that the "
                    "small scales there are real, so a viewer reading the "
                    "carrier as an image would see the noise put in.", ""]
    return "\n".join(out)


def figure(df, en, rd=None, path="figures/fig_interscale.pdf"):
    n_ax = 3 if rd is None or not len(rd) else 4
    fig, axes = plt.subplots(1, n_ax, figsize=(3.7 * n_ax, 3.3))

    ax = axes[0]
    d = en[en.orient != "a"].groupby(["scene", "level"]).var_frac.sum()
    med = d.groupby("level").median()
    ax.semilogy(med.index.values, 100 * med.values, "o-", color="0.2")
    ax.set_xlabel("detail level")
    ax.set_ylabel("% of tile variance")
    ax.set_title("subband energy")
    ax.invert_xaxis()
    ax.grid(alpha=.3)

    ax = axes[1]
    for m, lab in MODELS:
        s = df.groupby(pd.cut(np.log10(df.rmse), 10),
                       observed=True)[f"bpc_{m}"].median()
        ax.plot([iv.mid for iv in s.index], s.values, "o-", ms=3, label=lab)
    ax.set_xlabel("log10 RMSE")
    ax.set_ylabel("bits per coefficient")
    ax.set_title("rate by context model")
    ax.legend(fontsize=6, frameon=False)
    ax.grid(alpha=.3)

    ax = axes[2]
    for c, g in df.groupby("cls"):
        ax.scatter(g.rmse, 100 * g.gain_parent_over_neighbour, s=14,
                   color=CLASS_COLOR.get(c, "0.4"), label=c)
    ax.set_xscale("log")
    ax.axhline(0, color="0.6", lw=.8)
    ax.set_xlabel("RMSE")
    ax.set_ylabel("% rate the parent adds\nover the neighbour context")
    ax.set_title("what the link across scales adds")
    ax.legend(fontsize=6, frameon=False)
    ax.grid(alpha=.3)

    if n_ax == 4:
        ax = axes[3]
        for name, g in rd.groupby("reader"):
            ax.scatter(g.rmse_cell, g.transmission, s=14, label=name)
        ax.set_xscale("log")
        ax.axhline(1.0, color="0.6", lw=.8)
        ax.set_xlabel("cell-space RMSE")
        ax.set_ylabel(r"$\|Re\|/\|e\|$")
        ax.set_title("what the reader transmits")
        ax.legend(fontsize=6, frameon=False)
        ax.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--transform", choices=("haar", "cdf53"), default=TRANSFORM,
                    help="haar is the exact quadtree Haar; cdf53 is the "
                         "quadtree-native CDF 5/3 of hcc.compression."
                         "quadtree_wavelet, i.e. JPEG 2000's own filter on "
                         "this quadtree")
    ap.add_argument("--kernel", default=None,
                    choices=list(build_tiles.KERNELS),
                    help="which transfer's tiles to read (default: the "
                         "paper's Clough-Tocher tiles). Idea 3's PSF reader "
                         "is matched only to the PSF carrier, so its numbers "
                         "are meaningful only with --kernel psf")
    ap.add_argument("--tag", default=None,
                    help="suffix for the result files; defaults to the "
                         "transform when it is not the default")
    args = ap.parse_args(argv)
    if args.tag is not None:
        tag = args.tag
    else:
        tag = "" if args.transform == "haar" else f"_{args.transform}"
        if args.kernel and args.kernel != build_tiles.DEFAULT_KERNEL:
            tag += f"_{args.kernel}"

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if build_tiles.tile_path(s).exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print(f"no tile under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1

    bench = BENCH
    if args.kernel and args.kernel != build_tiles.DEFAULT_KERNEL:
        bench = BENCH.replace(".csv", f"__{args.kernel}.csv")
    front = pareto_front(bench)
    print(f"{len(scenes)} region(s), pyramid depth {DEPTH_MAX} "
          f"(levels {cfg.LEVEL}..{cfg.LEVEL - DEPTH_MAX + 1} detail)")

    rows, ens, zts, bands, reads = [], [], [], [], []
    for s in scenes:
        print(f"  {s}")
        try:
            r, e, z, b, k = region_rows(s, front, transform=args.transform,
                                        kernel=args.kernel)
        except Exception as exc:            # noqa: BLE001
            print(f"    skipped: {exc}")
            continue
        rows += r
        ens += e
        zts += z
        bands += b
        reads += k

    if not rows:
        print("no region produced a measurement")
        return 1

    ok = [r for r in reads if not r.get("error")]
    if READERS and not ok:
        errs = sorted({r.get("error", "?") for r in reads}) or ["no rows"]
        raise SystemExit(
            "idea 3 measured nothing: no reader in "
            f"{list(READERS)} could be built on any of the {len(scenes)} "
            "tiles. Distinct failures:\n  " + "\n  ".join(errs) +
            "\nFix that and rerun rather than publishing a run without it.")

    df = pd.DataFrame(rows)
    en = pd.DataFrame(ens)
    zt = pd.DataFrame(zts)
    bd = pd.DataFrame(bands)
    rd = pd.DataFrame(reads)
    df.to_csv(f"results/interscale{tag}.csv", index=False)
    en.to_csv(f"results/interscale_energy{tag}.csv", index=False)
    zt.to_csv(f"results/interscale_zerotree{tag}.csv", index=False)
    bd.to_csv(f"results/interscale_bands{tag}.csv", index=False)
    rd.to_csv(f"results/interscale_reader{tag}.csv", index=False)
    md = summarize(df, en, zt, rd)
    with open(f"results/interscale{tag}.md", "w", encoding="utf-8") as fh:
        fh.write(md)
    figure(df, en, rd,
           path=f"figures/fig_interscale{tag}.pdf")
    print("\n" + md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
