"""What a deconvolved HEALPix product costs to store (V3, step 3).

The question this answers is not "does the PSF transfer compress worse", which
is meaningless without saying at what quality, but:

    at the same quality *relative to what each product already carries*, how
    many more bytes does the deconvolved product need?

Two products, two references. `clough_tocher` represents the observed field;
`psf` represents the latent one. Comparing them at the same absolute RMSE
compares errors on two different quantities (Section V-H of the paper), so the
target is set per product:

    kappa = e_comp / e_resample        ->    target_p = kappa * e_resample(p)

and the surcharge is S_psf(target_psf) / S_ct(target_ct), paired by region.

The measurement splits into two halves that do NOT have the same standing, and
this module reports them separately rather than blending them:

    absolute   Both products swept on real Sentinel-2 tiles and compared at the
               same absolute RMSE. Nothing synthetic enters. This half is a
               fact about how well a deconvolved raster compresses, and it is
               unaffected by the recovery discrepancy of Section V-H.

    kappa      The same sweeps compared at equal quality relative to each
               product's own reconstruction error, which is measured on the
               synthetic harness. This half inherits that harness's caveats.

As it happens the two nearly coincide today: because the PSF reconstruction
recovers almost none of the blur, its error against the latent field is within
a few percent of Clough--Tocher's error against the same field, so equal kappa
is almost equal absolute target. That coincidence is itself reported, since a
reader should be able to see that the kappa result is not manufactured by the
choice of denominator.

A second transfer is swept as a control. Without it, "deconvolution costs
bytes" cannot be told apart from "any transfer other than Clough--Tocher costs
bytes"; `--control linear` is the cheap way to separate them.

Usage
-----
    python run_codec_bench.py --grid hpx --resampler psf     # produces the sweep
    python superres_cost.py                                 # then this
    python superres_cost.py --control linear
"""

from __future__ import annotations

import argparse
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import paper_stats as ps                  # noqa: E402
import build_tiles                        # noqa: E402

KAPPA_GRID = (0.1, 0.2, 0.5, 1.0)
CLASS_COLOR = {"urban": "#D55E00", "water": "#0072B2",
               "forest": "#009E73", "agriculture": "#E69F00"}
# Which reference each transfer answers to, as in resampling_error.py.
ESTIMAND = {"clough_tocher": "observed", "linear": "observed",
            "nearest": "observed", "bicubic": "observed", "psf": "latent"}


def inflation(kappa):
    return np.sqrt(1.0 + np.asarray(kappa, dtype=float) ** 2) - 1.0


def sweep(kernel: str) -> pd.DataFrame:
    """The HEALPix rate-distortion sweep for one transfer."""
    tag = "" if kernel == build_tiles.DEFAULT_KERNEL else f"__{kernel}"
    path = cfg.RESULTS_DIR / f"hpx_codec_bench{tag}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing; run\n"
            f"    python build_tiles.py --resampler {kernel}\n"
            f"    python run_codec_bench.py --grid hpx --resampler {kernel}")
    return pd.read_csv(path)


def utm_sweep():
    """The native-UTM sweep, which is the same whichever transfer is studied.

    The native product carries no transfer error, so kappa is undefined for it:
    there is no reconstruction error to divide by. The comparison is therefore
    made at the HEALPix product's own absolute target -- how many bytes each
    representation needs to hold this ground to that accuracy -- which is the
    convention the rest of the paper already uses.
    """
    path = cfg.RESULTS_DIR / "utm_codec_bench.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def denominators(kernel: str, denominator: str) -> dict:
    """e_resample per scene class for one transfer.

    `own` scores each transfer against the reference it claims; `latent` scores
    both against the unblurred scene, which is the only denominator under which
    the two products are answering the same question and therefore the one the
    headline surcharge uses.
    """
    path = cfg.RESULTS_DIR / "resampling_error.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run resampling_error.py")
    d = pd.read_csv(path).dropna(subset=["e_resample"])
    d = d[d.kernel == kernel]
    if not len(d):
        raise ValueError(f"no rows for kernel {kernel!r} in {path}; re-run "
                         f"resampling_error.py --kernels clough_tocher psf "
                         f"{kernel}")
    d = d.assign(cls=d.scene.str.split("__").str[1])
    est = ESTIMAND.get(kernel, "observed")
    if denominator == "latent" and est != "latent":
        col = "e_against_other_estimand"      # this transfer vs the latent scene
    else:
        col = "e_resample"
    return d.groupby("cls")[col].median().to_dict()


def sizes(df, scenes, codecs, targets_by_scene):
    """Interpolated compressed size at each scene's own target."""
    out = {}
    for s in scenes:
        t = targets_by_scene.get(s)
        if t is None or not np.isfinite(t):
            continue
        out[s] = dict(
            size=ps.best_of_codecs_size(df, s, codecs, t),
            bracketed=ps.any_bracketed(df, s, codecs, t, "rmse"),
            target=float(t))
    return out


def _vs_utm(scene, A, B, Ub, Ut, gamma=cfg.GAMMA_UTM):
    """Each HEALPix product against the operational native comparison.

    Composed per region rather than as a product of medians: the surcharge and
    the HPX/UTM ratio are both skewed, so multiplying their medians is not the
    median of the product.
    """
    def r(num, den):
        if den is None or not np.isfinite(den) or den == 0:
            return np.nan
        return num / (gamma * den)
    ub = Ub.get(scene, {}).get("size")
    ut = Ut.get(scene, {}).get("size")
    return dict(size_utm_base=ub, size_utm_test=ut,
                ratio_utm_base=r(A[scene]["size"], ub),
                ratio_utm_test=r(B[scene]["size"], ut))


def build(kernel="psf", baseline="clough_tocher", denominator="latent",
          codecs=None):
    codecs = codecs or cfg.COMPRESSOR_CODECS
    a, b, u = sweep(baseline), sweep(kernel), utm_sweep()
    shared = sorted(set(a.scene.unique()) & set(b.scene.unique()))
    if u is not None:
        shared = sorted(set(shared) & set(u.scene.unique()))
    if not shared:
        raise ValueError("the two sweeps share no region")
    e_base = denominators(baseline, denominator)
    e_test = denominators(kernel, denominator)

    cls = {s: s.split("__")[1] for s in shared}
    rows = []
    for kappa in KAPPA_GRID:
        tb = {s: kappa * e_base[cls[s]] for s in shared if cls[s] in e_base}
        tt = {s: kappa * e_test[cls[s]] for s in shared if cls[s] in e_test}
        A, B = sizes(a, shared, codecs, tb), sizes(b, shared, codecs, tt)
        Ub = sizes(u, shared, codecs, tb) if u is not None else {}
        Ut = sizes(u, shared, codecs, tt) if u is not None else {}
        for s in shared:
            if s not in A or s not in B:
                continue
            rows.append(dict(
                mode="kappa", kappa=kappa, scene=s, cls=cls[s],
                **_vs_utm(s, A, B, Ub, Ut),
                target_base=A[s]["target"], target_test=B[s]["target"],
                size_base=A[s]["size"], size_test=B[s]["size"],
                ratio=B[s]["size"] / A[s]["size"]
                if A[s]["size"] and np.isfinite(A[s]["size"]) else np.nan,
                bracketed_base=bool(A[s]["bracketed"]),
                bracketed_test=bool(B[s]["bracketed"]),
                bracketed=bool(A[s]["bracketed"] and B[s]["bracketed"])))
    # The same comparison at a shared absolute target: no synthetic denominator
    # enters, so this half stands on the real tiles alone.
    for t in cfg.TARGET_RMSE:
        tb = {s: t for s in shared}
        A, B = sizes(a, shared, codecs, tb), sizes(b, shared, codecs, tb)
        Ub = sizes(u, shared, codecs, tb) if u is not None else {}
        for s in shared:
            if s not in A or s not in B:
                continue
            rows.append(dict(
                mode="absolute", kappa=np.nan, scene=s, cls=cls[s],
                **_vs_utm(s, A, B, Ub, Ub),
                target_base=t, target_test=t,
                size_base=A[s]["size"], size_test=B[s]["size"],
                ratio=B[s]["size"] / A[s]["size"]
                if A[s]["size"] and np.isfinite(A[s]["size"]) else np.nan,
                bracketed_base=bool(A[s]["bracketed"]),
                bracketed_test=bool(B[s]["bracketed"]),
                bracketed=bool(A[s]["bracketed"] and B[s]["bracketed"])))
    df = pd.DataFrame(rows)
    df["kernel"], df["baseline"] = kernel, baseline
    df["denominator"] = denominator
    return df, e_base, e_test


def summarize(df, e_base, e_test, kernel, baseline, denominator):
    out = [f"# Storage surcharge of `{kernel}` over `{baseline}`", "",
           f"- denominator: {denominator}; regions: {df.scene.nunique()}", ""]
    out += ["| class | e_resample base | e_resample test | ratio |",
            "|---|---|---|---|"]
    for c in sorted(set(e_base) | set(e_test)):
        b, t = e_base.get(c, np.nan), e_test.get(c, np.nan)
        out.append(f"| {c} | {b:.4e} | {t:.4e} | {t/b:.3f} |")
    out += ["", "The two denominators being within a few percent of each other "
            "is not an assumption: it is what a recovery fraction near zero "
            "implies, and it means the kappa comparison below is close to the "
            "absolute one.", ""]

    for mode, label in (("kappa", "At equal kappa"),
                        ("absolute", "At equal absolute RMSE")):
        sub = df[(df["mode"] == mode)]
        if not len(sub):
            continue
        out += [f"## {label}", "",
                "| target | n | bracketed | cens. base | cens. test | "
                "median size ratio | IQR | p |",
                "|---|---|---|---|---|---|---|---|"]
        key = "kappa" if mode == "kappa" else "target_base"
        meds = {}
        for t, g in sub.groupby(key):
            ok = g[g.bracketed]
            frac = len(ok) / len(g) if len(g) else 0.0
            cb = 1.0 - g.bracketed_base.mean()
            ctst = 1.0 - g.bracketed_test.mean()
            st = ps.paired_ratio_test(g.ratio)
            meds[float(t)] = (st["median"], frac)
            flag = "" if frac >= cfg.MIN_BRACKETED else " *(withheld)*"
            extra = (f" (+{100*float(inflation(t)):.1f}% total error)"
                     if mode == "kappa" else "")
            out.append(f"| {t:g}{extra}{flag} | {st['n']} | {100*frac:.0f}% | "
                       f"{100*cb:.0f}% | {100*ctst:.0f}% | "
                       f"{st['median']:.3f} | [{st['q25']:.3f}, "
                       f"{st['q75']:.3f}] | {st['p']:.1e} |")
        out.append("")
        # A censored region contributes the size of the coarsest swept setting,
        # which is larger than the size a genuinely coarser setting would need.
        # The smoother product reaches a loose target first, so censoring
        # inflates the baseline more than the deconvolved side and pulls the
        # ratio *down*: where it bites, the reported surcharge is a floor.
        keys = sorted(meds)
        drops = [(a, b) for a, b in zip(keys, keys[1:])
                 if meds[b][0] < meds[a][0] and meds[b][1] < 1.0]
        if drops:
            a, b = drops[-1]
            out.append(
                f"- The ratio falls from {meds[a][0]:.2f} at {a:g} to "
                f"{meds[b][0]:.2f} at {b:g}, where only "
                f"{100*meds[b][1]:.0f}% of regions are bracketed on both "
                f"sides. Censoring inflates the smoother product's size more "
                f"than the deconvolved one's, so it biases the ratio towards "
                f"1: the fall is a property of the sweep's reach, not a turn "
                f"in the surcharge, and the affected values are floors.")
            out.append("")

    # The decision-relevant number: does the deconvolved product still beat
    # the operational native comparison, or does the surcharge spend the whole
    # advantage the non-deconvolving product had?
    if "ratio_utm_test" in df and df.ratio_utm_test.notna().any():
        out += [f"## Against native UTM (gamma = {cfg.GAMMA_UTM})", "",
                "Each HEALPix product's compressed size over the equal-area "
                "native crop's, scaled by the assumed overhead, at that "
                "product's own target. Composed per region, not as a product "
                "of medians. Values below 1 favour HEALPix. Restricted to the "
                "regions bracketed on both HEALPix sweeps, so these differ "
                "slightly from the paper's own HPX/UTM table, which admits "
                "censored bounds over all 40 regions; the comparison that "
                "matters here is between the two columns, which are computed "
                "on the same regions.", "",
                "| target | n | " + f"`{baseline}`/UTM | `{kernel}`/UTM | "
                "still favours HEALPix |", "|---|---|---|---|---|"]
        for mode, key in (("kappa", "kappa"), ("absolute", "target_base")):
            sub = df[(df["mode"] == mode) & df.bracketed]
            if not len(sub):
                continue
            for t, g in sub.groupby(key):
                if len(g) / max(len(df[(df["mode"] == mode)
                                       & (df[key] == t)]), 1) < cfg.MIN_BRACKETED:
                    continue
                rb = float(np.nanmedian(g.ratio_utm_base))
                rt = float(np.nanmedian(g.ratio_utm_test))
                if not np.isfinite(rb) or not np.isfinite(rt):
                    continue
                lab = (f"kappa {t:g}" if mode == "kappa" else f"RMSE {t:g}")
                verdict = ("both" if rt < 1 and rb < 1 else
                           "baseline only" if rb < 1 else "neither")
                out.append(f"| {lab} | {len(g)} | {rb:.2f} | {rt:.2f} | "
                           f"{verdict} |")
        out.append("")
        won = df[(df["mode"] == "absolute") & df.bracketed]
        if len(won):
            per_t = won.groupby("target_base")[["ratio_utm_base",
                                                "ratio_utm_test"]].median()
            flipped = per_t[(per_t.ratio_utm_base < 1)
                            & (per_t.ratio_utm_test >= 1)]
            if len(flipped):
                out.append(
                    f"- **At {', '.join(f'{t:g}' for t in flipped.index)} the "
                    f"non-deconvolving product is the cheaper representation "
                    f"and the deconvolved one is not**: the surcharge spends "
                    f"the whole margin. This is the storage price of "
                    f"super-resolution, stated as the thing a product designer "
                    f"has to trade away.")
                out.append("")

    k = df[(df["mode"] == "kappa")]
    rep = []
    for t, g in k.groupby("kappa"):
        if len(g) and (g.bracketed.sum() / len(g)) >= cfg.MIN_BRACKETED:
            rep.append((t, float(np.nanmedian(g.ratio))))
    if rep:
        t, r = rep[0]
        out.append(f"- **At kappa = {t:g}, the deconvolved product costs "
                   f"{r:.2f} times the bytes of the non-deconvolving one** "
                   f"for the same relative quality "
                   f"(+{100*float(inflation(t)):.1f}% total error each).")
        out.append(f"- Across the reportable kappa grid the surcharge runs "
                   f"{min(r for _, r in rep):.2f} to "
                   f"{max(r for _, r in rep):.2f}.")
    else:
        out.append("- no kappa target is bracketed on both sweeps for at "
                   f"least {100*cfg.MIN_BRACKETED:.0f}% of regions; nothing "
                   "is reportable and the sweep range needs extending.")
    return "\n".join(out)


def figure(df, kernel, baseline, publish=True):
    k = df[(df["mode"] == "kappa") & df.bracketed]
    if not len(k):
        return
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    for cls, g in k.groupby("cls"):
        m = g.groupby("kappa").ratio.median()
        ax.plot(m.index, m.values, "o-", ms=4, lw=1.6,
                color=CLASS_COLOR.get(cls, "0.4"), label=cls)
    ax.axhline(1.0, color="0.6", ls="--", lw=1.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\kappa=e_{\rm comp}/e_{\rm resample}$", fontsize=9)
    ax.set_ylabel(f"{kernel} / {baseline} compressed size", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    name = f"fig_superres_cost_{kernel}"
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight",
                    dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", default="psf",
                    help="the transfer whose storage cost is in question")
    ap.add_argument("--baseline", default="clough_tocher")
    ap.add_argument("--control", default=None,
                    help="a second, non-deconvolving transfer swept the same "
                         "way. Without it a surcharge cannot be attributed to "
                         "deconvolution rather than to changing kernel at all")
    ap.add_argument("--denominator", choices=["latent", "own"],
                    default="latent",
                    help="'latent' scores both products against the unblurred "
                         "scene, the only reference under which they answer "
                         "the same question; 'own' scores each against its own "
                         "estimand, which makes the two kappa grids "
                         "incomparable and is offered for inspection only")
    ap.add_argument("--codecs", nargs="*", default=None)
    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args(argv)

    todo = [args.kernel] + ([args.control] if args.control else [])
    # Check every transfer's denominator up front. Discovering a missing one
    # after the first analysis has already written its outputs leaves half a
    # comparison on disk and no control.
    missing = []
    for kern in [args.baseline] + todo:
        try:
            denominators(kern, args.denominator)
        except (FileNotFoundError, ValueError) as exc:
            missing.append(f"  {kern}: {exc}")
    if missing:
        print("no reconstruction error available for every transfer, so "
              "nothing is comparable yet:")
        print("\n".join(missing))
        print(f"\nrun\n    python resampling_error.py --kernels "
              f"{' '.join(dict.fromkeys([args.baseline] + todo))}")
        return 1
    for kern in todo:
        df, e_base, e_test = build(kern, args.baseline, args.denominator,
                                   args.codecs)
        tag = f"_{kern}"
        cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        csv = cfg.RESULTS_DIR / f"superres_cost{tag}.csv"
        df.to_csv(csv, index=False)
        md = summarize(df, e_base, e_test, kern, args.baseline,
                       args.denominator)
        (cfg.RESULTS_DIR / f"superres_cost{tag}.md").write_text(md,
                                                                encoding="utf-8")
        print("\n" + md)
        print(f"wrote {csv}")
        figure(df, kern, args.baseline, publish=not args.no_publish)
    if args.control:
        print("\nRead the two together: a surcharge that the control also "
              "shows is a cost of changing transfer, not of deconvolving.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
