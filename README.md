# ⚛️ An Upper Bound on the Compactification Scale of SU(7) Grand Gauge-Higgs Unification

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22087251-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.22087251)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Gates](https://img.shields.io/badge/gates-15_green-1B6F8C)](paper/)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](.)
[![Reproduction](https://img.shields.io/badge/scripts-70%2F70_byte--for--byte-1B6F8C)](outputs/)

**…and the dijet angular distribution that tests it.**

**📄 Paper (EN + ES) and every verification script on Zenodo → https://doi.org/10.5281/zenodo.22087251**

That is the *concept* DOI: it always resolves to the current version. The paper, its sources and
every archived run are also in [`paper/`](paper/), [`scripts/`](scripts/) and [`outputs/`](outputs/).

> ### 📚 Part **VII** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU***
>   → [github.com/karlesmarin/ghu-su4-completion](https://github.com/karlesmarin/ghu-su4-completion) · [Zenodo 10.5281/zenodo.21432625](https://doi.org/10.5281/zenodo.21432625)
> - **Part II — *Three Gates to a Quark Generation***
>   → [github.com/karlesmarin/su4-sm-cell-criterion](https://github.com/karlesmarin/su4-sm-cell-criterion) · [Zenodo 10.5281/zenodo.21432627](https://doi.org/10.5281/zenodo.21432627)
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential***
>   → [github.com/karlesmarin/centre-parity-selection](https://github.com/karlesmarin/centre-parity-selection) · [Zenodo 10.5281/zenodo.21438226](https://doi.org/10.5281/zenodo.21438226)
> - **Part IV — *Schur Functions at (1,−1,t,t⁻¹)***
>   → [github.com/karlesmarin/schur-nonidentity-o4](https://github.com/karlesmarin/schur-nonidentity-o4) · [Zenodo 10.5281/zenodo.21463000](https://doi.org/10.5281/zenodo.21463000)
> - **Part V — *What the Higgs Potential Cannot See***
>   → [github.com/karlesmarin/higgs-blind-class](https://github.com/karlesmarin/higgs-blind-class) · [Zenodo 10.5281/zenodo.21727094](https://doi.org/10.5281/zenodo.21727094)
> - **Part VI — *Proton Decay in SU(7) Grand Gauge-Higgs Unification***
>   → [github.com/karlesmarin/su7-proton-row](https://github.com/karlesmarin/su7-proton-row) · [Zenodo 10.5281/zenodo.22033302](https://doi.org/10.5281/zenodo.22033302)
> - **Part VII — *An Upper Bound on the Compactification Scale of SU(7) GHU*** (this repo)
> - **Part VIII — *A Certified 2.68 TeV Gap in the Closed-Form Map of the Compactification Scale***
>   → [github.com/karlesmarin/su7-certified-gap](https://github.com/karlesmarin/su7-certified-gap) · [Zenodo 10.5281/zenodo.22159036](https://doi.org/10.5281/zenodo.22159036)
> - **Part IX-A — *The Alphabet of Orbifold Boundary Conditions***
>   → [github.com/karlesmarin/orbifold-alphabet](https://github.com/karlesmarin/orbifold-alphabet) · [Zenodo 10.5281/zenodo.22254861](https://doi.org/10.5281/zenodo.22254861)
> - **Part IX-B — *An Affine Semigroup from Orbifold Boundary Conditions***
>   → [github.com/karlesmarin/orbifold-semigroup](https://github.com/karlesmarin/orbifold-semigroup) · [Zenodo 10.5281/zenodo.22254863](https://doi.org/10.5281/zenodo.22254863)

In gauge-Higgs unification the Higgs **is** the Wilson-line phase, so the electroweak hierarchy is a
minimum of a potential rather than a ratio of scales. Part VI computed the curvature of that vacuum
for the SU(7) model of Komori and Maru (arXiv:2503.04090) and reached a sign. This paper computes
the vacuum itself, and then asks how large the compactification scale can be over *arbitrary*
multiplicities of the representations they put in the bulk.

The answer is finite, and the second half of the title is why that matters: a collider recast of
CMS's dijet angular distribution approaches the same scale **from below**, so the two calculations
cross.

## 🙏 With thanks

This work exists because Komori and Maru wrote their model out fully — every parity matrix, every
mode expansion, the whole one-loop potential — so that a reader outside their group can rebuild it
from its own equations rather than transcribe it. That is not the norm and it is what made this
paper possible.

Fifty-nine quoted equations from eleven sources carry the load, and it is not spread evenly:
twenty-five from Komori–Maru, fourteen from Cacciapaglia, Cossu and Deandrea on orbifold vacuum
stability, seven from Haba, Hosotani and Kawamura and five from Haba and Yamashita — the two papers
that made Wilson-line potentials on orbifolds something one computes rather than models. §12 of the
paper says which piece is whose, one by one, and a gate keeps that table honest.

## 🧭 What is in it

| | |
|---|---|
| **The ladder** | Expanding the one-loop potential about the symmetric point gives one identity read at three rungs, `p = 5, 3, 1`. The value, the curvature and the fourth moment are the *same* formula. At `p = 1` the antiperiodic half cancels its pole into a finite `ln 2` while the periodic pole survives as the `x⁴ln x` branch point — which is why no polynomial exists. |
| **The closed form** | `α_min` in Lambert-`W`, `x² = −d/W(−d e^{−b})`, in three moments. Not an algebraic root — it carries `ln x` — but *solved* rather than iterated, with `W₋₁` the minimum and `W₀` the maximum, proved rather than assigned. Median error `0.13 %` over 272 contents. |
| **Both columns for free** | The logarithm and `G` cancel identically at the stationary point, so `α_min` and `m_h` both follow from two moments with **no minimisation anywhere**. |
| **The arithmetic** | `8D ≡ 2A₄ + 3 (mod 6)`, because `c⁴ ≡ c² (mod 3)`. The matter half holds whatever the gauge sector does; the residue is the gauge base point's, and both seeds considered give the same `9`. |
| **The ceiling** | An integer program over a rational cone, closed by an exact two-variable dual: `1/R₅ ≤ 10.03 TeV` for arbitrary bulk content, sharpened to `10.01` attained once the lattice is lifted to four dimensions, and to **`9.22 TeV`** once the vacuum is required to be the *true* one — on the three-rung algebraic surface. At the measured Higgs mass, `9.09 TeV`. |
| **Which state a collider bounds** | Read off their three parity matrices and nothing else. Every colour-changing gauge boson is odd under the small-circle parity, so its tower starts at `1/R₆` and its wavefunction **vanishes** on the quark brane — which answers, for the gauge channel, the rapid proton decay their conclusions flag and leave to future work. |
| **The recast** | The surviving coloured tower acts through a colour-**octet** operator in the dijet *angular* distribution, so published contact-interaction limits, stated for a singlet, do not directly apply. Carried through CMS's unfolded distribution, 77×77 covariance and NNLO baseline: a nominal sensitivity of `1/R₅ > 10.86 TeV`. **Deliberately not called an exclusion** — a Gaussian χ² over unfolded data is not a detector-level `CLₛ`, and the paper measures that gap rather than waving at it. |
| **Where it stops** | One inherited coefficient this paper does not compute: the gauge weight of their eq. (68). Everything needing it to be *odd* is marked conditional, in the abstract, on the front page, in the theorem's own signature and in the ledger. |

## ⚠️ The one condition, stated plainly

Theorem 1 carries its hypothesis **in its own statement**: *if* the gauge sector's contribution to
`8D` is odd, then `8D` is an odd integer for every bulk content. On the coefficients Komori and Maru
print, it is. On a parity-resolved covariant count it would not be — and this paper does part of
that count, finding the split `ξ`-independent and equal to `3 + 1`, which is the candidate.

So the fork is **not symmetric**, and the paper says so: its own machinery points twice, by
different routes, at the branch where that hypothesis fails. What is *not* settled is whether the
ghosts' boundary conditions are the assumed ones, and whether their spectrum construction hides a
cancellation. Both branches are costed in a table, and no number is quoted without saying which
branch it lives on.

## 🔍 Reproducing it

Every displayed number regenerates from the ancillary scripts, whose standard output is archived
beside them in `outputs/`. The ledger in the paper cites the **script** for each entry rather than
the number alone, so a reader who disagrees with a claim can re-run the thing that made it.

```bash
python check_reproduces.py     # runs every cited script and diffs it against its archived output
```

Fifteen gates guard the two editions — attribution, branch conditionality, figures, formulas,
inline maths, framed boxes, layout, narrative, numbers, parity, prose, references, scripts,
EN/ES synchrony, and reproduction. `comprobar.py` decides which of them a given change can break.

Three scripts need CIJET's output, whose licence permits academic use and forbids redistribution;
they are named, with the reason, rather than silently skipped. Five `.sage` files need the
container.

## 📁 Layout

```
paper/     the two editions, their sources, figures, and the fifteen gates
scripts/   every cited computation
outputs/   the archived standard output of each, byte for byte
```

---

*Carles Marín · independent researcher · [ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688)*
*Written with Claude (Anthropic) as a research instrument.*
