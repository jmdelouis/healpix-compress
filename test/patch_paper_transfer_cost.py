# -*- coding: utf-8 -*-
"""Adds Section V-I (what the transfer costs to store) and rewrites the two
limitation items the new measurement settles."""
import sys
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else "tgrs_healpix_compress.tex")
s = p.read_text(encoding="utf-8")
E = []

SECTION = r"""
\subsection{What the Transfer Costs to Store}
\label{sec:transfer_cost}
Section~\ref{sec:error_budget} priced the compression step against the error
each product already carries. This section prices the transfer itself. The same
rate--distortion protocol is run on carrier tiles rebuilt from the same native
windows by two further transfers: the PSF-aware one, which restores small-scale
structure by construction, and plain bilinear interpolation, which restores
none and serves as the control. That control is what makes the question
falsifiable. A deconvolved tile costing more bytes than a Clough--Tocher tile
would tell us nothing on its own, because every property of a deconvolving
transfer is confounded with the mere fact of its not being Clough--Tocher; only
a non-deconvolving transfer swept the same way separates the two. The PSF tiles
are built at the same damping the reconstruction error was measured at
($\lambda=10^{-3}$) and at the nominal \SI{12.5}{m} response width, both
recorded in each tile. Coverage is 40 regions for the PSF transfer and 39 for
bilinear, one open-water region being absent from the bilinear sweep.

\input{tables/table_transfer_cost}

The control fires. At the same absolute RMSE, the deconvolved tile costs $1.27$
to $2.11$ times the Clough--Tocher tile of the same region, rising with the
target; bilinear costs $1.21$ to $2.01$, which is very nearly the same. Taken
against each other, the deconvolving transfer exceeds the non-deconvolving one
by $1.04$ to $1.09$ --- separable from parity at the three tightest targets
($p=4\times10^{-10}$ to $1.7\times10^{-4}$) and not at the two loosest.
Whatever the roughly two-fold cost at loose targets is, it is not the price of
restoring spatial detail, because the transfer that restores none pays it too.
The isolated cost of deconvolving, over the closest non-deconvolving transfer
available, is a few percent; the rest is the cost of leaving Clough--Tocher.

\input{tables/table_compactability}

What separates the transfers is not how much each raster carries but how
compactly. Their lossless sizes are within $9\%$ of each other
(Table~\ref{tab:compactability}): \SI{1.31}{MB} for Clough--Tocher,
\SI{1.37}{MB} for the PSF transfer, \SI{1.41}{MB} for bilinear, so none of them
holds appreciably more information than the others. Their rate--distortion
curves, however, have very different slopes. At the $10^{-2}$ target the
Clough--Tocher tile compresses to $1.1\%$ of its own lossless size, the
bilinear tile to $2.3\%$ and the deconvolved tile to $2.5\%$; the native UTM
raster, for comparison, reaches $3.4\%$ of its own, which is the oversampling
redundancy this paper is about and which every transfer inherits. The
explanation is energy compaction. Clough--Tocher is a $C^1$ piecewise-cubic
interpolant that approximately minimises curvature, so its energy sits in
coarse wavelet coefficients and the fine ones are near zero and free to
discard; bilinear interpolation is only $C^0$ and leaves gradient
discontinuities at cell boundaries, while the damped deconvolution restores
fine structure and amplifies the noise along with it, and both put energy where
a coder targeting an RMSE has to keep it. The lossless size does not predict
the byte ratio across regions --- Spearman $\rho=-0.47$ for the PSF transfer
and $+0.06$ for bilinear --- which is what rules out total information content
as the mechanism and leaves compaction.

The consequence for this paper's headline is direct, and it is the largest
single sensitivity we have measured. Under Clough--Tocher the HEALPix carrier
is the cheaper representation at the two loosest targets ($0.94$ and $0.78$
against the operational native comparison). Under either of the other two
transfers it is cheaper at none: $2.00$ and $1.74$ for the PSF transfer,
$1.90$ and $1.77$ for bilinear. A factor of two in the transfer's storage cost
is more than the whole margin, so the regime-dependent verdict of
Section~\ref{sec:results} is a statement about a carrier populated by this
transfer, not about the grid.

\input{tables/table_kernel_ablation}

That is not, however, a reason to read the choice as opportunistic, and the
accuracy half of the ablation is why. Scored against an analytic field known at
every position --- the only measure in Table~\ref{tab:kernel_ablation} that
gives a transfer no credit for copying its input samples --- Clough--Tocher is
at $1.19\%$ of the field's standard deviation against $8.78\%$ for the PSF
transfer, $9.84\%$ for nearest and $12.07\%$ for bilinear, on three regions.
Clough--Tocher is therefore simultaneously the most accurate transfer tested
and the most compressible one, and the two properties have a common cause: a
transfer that represents a band-limited field faithfully does not invent
fine-scale content, and content it does not invent is content the codec is not
asked to pay for. Bilinear is dominated rather than an alternative --- twice
the bytes and ten times the transfer error --- so the comparison that matters
is with the PSF transfer, which buys its extra bytes with an estimand this
paper's targets are not written against.

The honest scope is therefore narrower than a grid comparison and more useful
than a kernel preference. What is measured is that a HEALPix level-20 carrier
becomes competitive at aggressive compression when it is populated by a
transfer in Clough--Tocher's smoothness class, that leaving that class costs
about a factor of two at the loosest targets and removes the advantage, and
that this cost is a property of how compactly the transfer's output can be
coded rather than of how much detail it restores.

"""

E.append(("\n\\section{Discussion}\n", SECTION + "\n\\section{Discussion}\n"))

# --- the kernel limitation is now measured, not pending ---------------------
OLD_KERNEL = (
    "    \\item \\textbf{The resampling kernel used to populate the HEALPix grid "
    "plausibly affects the result, but this is not yet quantified.}")
NEW_KERNEL = (
    "    \\item \\textbf{The verdict is conditional on the transfer, by about a "
    "factor of two, and the transfer this paper uses is the most favourable of "
    "those tested.} Section~\\ref{sec:transfer_cost} measures this rather than "
    "leaving it open: rebuilding the carrier with a PSF-aware or a bilinear "
    "transfer costs $1.2$ to $2.1$ times the bytes at the same accuracy, which "
    "is more than the margin by which HEALPix leads at the loosest targets, so "
    "under either of them it leads at none. Clough--Tocher is also the most "
    "accurate of the transfers tested against an analytic field, and we argue "
    "the two properties share a cause, but the argument does not extend the "
    "result to transfers we have not swept --- the original inverse-distance "
    "weighting among them --- and it does not make the reported ratios "
    "properties of the grid. Three residual gaps remain in this ablation: the "
    "accuracy column of Table~\\ref{tab:kernel_ablation} is measured on three "
    "regions rather than forty, we have not swept the IDW resampler that "
    "Section~\\ref{sec:resampler} replaced, and we have not tested a "
    "non-deconvolving transfer of Clough--Tocher's smoothness class, which is "
    "what would isolate deconvolution from smoothness properly rather than "
    "bounding it as we do;")
assert s.count(OLD_KERNEL) == 1, "kernel limitation anchor"
i = s.index(OLD_KERNEL)
j = s.index("\n    \\item \\textbf{The equal-area UTM square", i)
s = s[:i] + NEW_KERNEL + s[j:]

# --- the super-resolution item is now measured ------------------------------
OLD_SR = "    \\item \\textbf{The super-resolution surcharge is not measured, and"
assert s.count(OLD_SR) == 1, "super-resolution limitation anchor"
i = s.index(OLD_SR)
j = s.index("\n    \\item", i + 10)
NEW_SR = (
    "    \\item \\textbf{The super-resolution surcharge is measured, but our "
    "deconvolution control does not reproduce the companion's recovery, so the "
    "product it prices is not the one that study validates.} "
    "Section~\\ref{sec:transfer_cost} prices the deconvolved carrier at $1.27$ "
    "to $2.11$ times the Clough--Tocher one, of which only a few percent is "
    "attributable to deconvolving rather than to leaving Clough--Tocher's "
    "smoothness class. That is a measurement on real imagery and does not "
    "depend on the synthetic harness. What does depend on it is the claim that "
    "the extra bytes buy a better estimand: in our harness the PSF-aware "
    "reconstruction recovers $+3.0\\%$ to $+6.5\\%$ of the blur without noise "
    "and nothing at the measured noise level, against the $15.4\\%$ to "
    "$19.6\\%$ the companion study reports noise-free. We identify the added "
    "noise and a latent spectrum extrapolated to the level-20 Nyquist as "
    "contributors, exclude the depth of the cell average by measurement, and "
    "are left with a residual gap we cannot attribute. The storage price is "
    "therefore reportable and the quality it purchases is not, which is the "
    "asymmetry a reader should carry away rather than a net verdict on "
    "PSF-aware carriers;")
s = s[:i] + NEW_SR + s[j:]

# --- abstract: the conditionality belongs there -----------------------------
E.append((
    """We summarise this as a critical overhead
falling from $2.07$ to $1.08$ with the distortion target, and report both
measures rather than electing one.""",
    """We summarise this as a critical overhead
falling from $2.07$ to $1.08$ with the distortion target, and report both
measures rather than electing one. Both are conditional on the transfer that
populates the carrier: rebuilding the tiles with a PSF-aware or a bilinear
transfer costs $1.2$ to $2.1$ times the bytes at the same accuracy and leaves
HEALPix ahead at no target, and we show the cause to be energy compaction
rather than added content."""))

# --- conclusion -------------------------------------------------------------
E.append((
    """Separating the transfer error from the read, and from the resampling operator that produces the first, is the natural continuation of this work, and the reason this study shares its sites and its acquisitions with a companion study of HEALPix resampling~\\cite{delouis2026resample}.""",
    """Separating the transfer error from the read, and from the resampling operator that produces the first, is the natural continuation of this work, and the reason this study shares its sites and its acquisitions with a companion study of HEALPix resampling~\\cite{delouis2026resample}. How much the operator matters is now measured rather than assumed: rebuilding the carrier with a PSF-aware or a bilinear transfer costs $1.2$ to $2.1$ times the bytes at the same accuracy (Section~\\ref{sec:transfer_cost}), more than the margin by which HEALPix leads at the loosest targets. Only a few percent of that is the price of deconvolving; the rest is the price of leaving a $C^1$ curvature-minimising interpolant, whose output a lossy coder can discard far more of at equal error. The competitive regime this paper reports therefore belongs to a carrier populated by a transfer in that class, and Clough--Tocher being also the most accurate transfer tested against an analytic field is, we argue, the same fact seen from the other side."""))

for old, new in E:
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)

p.write_text(s, encoding="utf-8")
print("paper patched")
