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
dead idea.

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
import region_sites                       # noqa: E402

DEPTH_MAX = 6                  # detail at levels 20..15, approximation at 14
N_CTX = 6                      # magnitude buckets in a context symbol
BENCH = "results/hpx_codec_bench.csv"
ORIENT = ("h", "v", "d")
MODELS = (("order0", "no context"),
          ("parent", "parent (inter-scale)"),
          ("neighbour", "neighbour (intra-band)"),
          ("both", "parent + neighbour"))
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


# --------------------------------------------------------------------------
# quantisation, rate, distortion
# --------------------------------------------------------------------------
def quantise(c, step):
    """Uniform deadzone quantiser; returns (indices, reconstruction)."""
    q = np.floor(np.abs(c) / step) * np.sign(c)
    r = np.where(q == 0, 0.0, np.sign(q) * (np.abs(q) + 0.5) * step)
    return q.astype(np.int64), r


def rmse_at(approx, levels, step):
    """Distortion only -- cheap, so the step can be bisected."""
    n_leaf = approx.size * (4 ** len(levels))
    _, r = quantise(approx, step)
    err2 = float(((approx - r) ** 2).sum())
    for det in levels:
        for o in ORIENT:
            _, rr = quantise(det[o], step)
            err2 += float(((det[o] - rr) ** 2).sum())
    return float(np.sqrt(err2 / n_leaf))


def solve_step(approx, levels, target_rmse, sd, iters=34):
    """Bisect the quantiser step onto a given distortion."""
    lo, hi = 1e-8 * sd, 8.0 * sd
    for _ in range(iters):
        mid = float(np.sqrt(lo * hi))
        if rmse_at(approx, levels, mid) > target_rmse:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.0002:
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


def coder_point(approx, levels, step):
    """Quantise every subband at one step; return (rmse, bits per model)."""
    n_leaf = approx.size * (4 ** len(levels))
    qa, ra = quantise(approx, step)
    err2 = float(((approx - ra) ** 2).sum())
    h_a = approx.size * entropy0(qa)          # no parent above the top
    bits = {k: h_a for k, _ in MODELS}
    prev = {o: None for o in ORIENT}
    per_band = []
    for k in range(len(levels) - 1, -1, -1):  # coarsest detail first
        for o in ORIENT:
            c = levels[k][o]
            q, r = quantise(c, step)
            err2 += float(((c - r) ** 2).sum())
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
    return float(np.sqrt(err2 / n_leaf)), bits, n_leaf, per_band


def zerotree_stats(levels, step):
    """P(child = 0 | parent = 0) against P(child = 0): the zerotree premise."""
    out = []
    for k in range(len(levels) - 1, 0, -1):
        for o in ORIENT:
            qp, _ = quantise(levels[k][o], step)
            qc, _ = quantise(levels[k - 1][o], step)
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
def region_rows(scene, front, verbose=True):
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    x = np.asarray(cell_data, dtype=np.float64).ravel()
    if x.size != cfg.N_CHILDREN:
        raise ValueError(f"{x.size} cells, expected {cfg.N_CHILDREN}")
    sd = float(np.std(x))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("degenerate tile")
    approx, levels = pyramid(x)

    # The transform must be exact before anything measured from it means
    # anything.
    rt = float(np.max(np.abs(rebuild(approx, levels) - x)))
    scale = max(1.0, float(np.max(np.abs(x))))
    if rt > 1e-9 * scale:
        raise RuntimeError(f"Haar round-trip error {rt:.3e}")

    cls = scene.split("__")[1]
    energies = [dict(scene=scene, cls=cls, level=cfg.LEVEL - k, orient=o,
                     n_coef=levels[k][o].size,
                     var_frac=float((levels[k][o] ** 2).sum()
                                    / (x.size * sd ** 2)))
                for k in range(len(levels)) for o in ORIENT]
    energies.append(dict(scene=scene, cls=cls, level=cfg.LEVEL - DEPTH_MAX,
                         orient="a", n_coef=approx.size,
                         var_frac=float(((approx - approx.mean()) ** 2).sum()
                                        / (x.size * sd ** 2))))

    targets = front.get(scene) or [(t, np.nan) for t in
                                   np.geomspace(1e-4 * sd, 0.5 * sd, 6)]
    rows, ztree, bands = [], [], []
    for target_rmse, ref_bytes in targets:
        step = solve_step(approx, levels, target_rmse, sd)
        rmse, bits, n_leaf, per_band = coder_point(approx, levels, step)
        row = dict(scene=scene, cls=cls, step=step, target_rmse=target_rmse,
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
        rows.append(row)
        ztree += [dict(scene=scene, target_rmse=target_rmse, **z)
                  for z in zerotree_stats(levels, step)]
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
    return rows, energies, ztree, bands


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
def summarize(df, en, zt):
    n = en.scene.nunique()
    out = ["# Does the link between scales pay?", "",
           f"- regions: {n}; rates are ideal entropy bounds, so they flatter "
           f"this coder against the measured bytes of real codecs", ""]

    out += ["## Subband energy of the exact HEALPix Haar pyramid", "",
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
    return "\n".join(out)


def figure(df, en, path="figures/fig_interscale.pdf"):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.3))

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
        s = df.groupby(pd.cut(np.log10(df.rmse), 10))[f"bpc_{m}"].median()
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
    args = ap.parse_args(argv)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if build_tiles.tile_path(s).exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print(f"no tile under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1

    front = pareto_front()
    print(f"{len(scenes)} region(s), pyramid depth {DEPTH_MAX} "
          f"(levels {cfg.LEVEL}..{cfg.LEVEL - DEPTH_MAX + 1} detail)")

    rows, ens, zts, bands = [], [], [], []
    for s in scenes:
        print(f"  {s}")
        try:
            r, e, z, b = region_rows(s, front)
        except Exception as exc:            # noqa: BLE001
            print(f"    skipped: {exc}")
            continue
        rows += r
        ens += e
        zts += z
        bands += b

    if not rows:
        print("no region produced a measurement")
        return 1

    df = pd.DataFrame(rows)
    en = pd.DataFrame(ens)
    zt = pd.DataFrame(zts)
    bd = pd.DataFrame(bands)
    df.to_csv("results/interscale.csv", index=False)
    en.to_csv("results/interscale_energy.csv", index=False)
    zt.to_csv("results/interscale_zerotree.csv", index=False)
    bd.to_csv("results/interscale_bands.csv", index=False)
    md = summarize(df, en, zt)
    with open("results/interscale.md", "w", encoding="utf-8") as fh:
        fh.write(md)
    figure(df, en)
    print("\n" + md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
