# Sampling-Budget Audit of Dimensionality Reduction for Conformational Landscapes

Auditing whether Gaussian-mixture component counts in PCA, TICA, and VAE
projections are properties of the free-energy landscape or partly artifacts of
how much data was collected. These are apparent basins, not kinetically
validated metastable states.

## The question

A very common claim in the molecular simulation literature has the form:

> Projecting our trajectory with method M reveals N metastable basins.

That sentence presents N as a property of the molecule. This repository tests
whether it is also, substantially, a property of the sampling budget — the
number of frames the analyst happened to save.

The test is only meaningful if the true answer is known independently of the
method being audited. So the primary systems are potentials where the number,
location, and population of basins are fixed **by construction** (Müller-Brown,
3 basins; Prinz, 4 basins), with the low-dimensional dynamics lifted through a
fixed random nonlinear map into 30 observed dimensions so that dimensionality
reduction is a non-trivial task. Ground-truth labels always come from the
latent coordinate, never from a fitted model.

## The control that makes it defensible

Two sampling modes are run over the identical landscape:

| mode | what varies | what it isolates |
|---|---|---|
| `short` | a genuinely short trajectory | exploration deficit **and** sample size, confounded |
| `subsample` | uniform thinning of one long trajectory | sample size alone, coverage held fixed |

If the effect survives `subsample`, "short simulations explore less" does not
explain it, and the inflation lives in the estimator. If it vanishes, it is an
exploration deficit. Running only `short` and asserting either one would not be
a result.

A second analysis distinguishes criterion-specific behavior. Five diagnostics
are computed on every condition (BIC, AIC, ICL, silhouette, and a simple elbow
rule), but they do not behave alike and are not presented as five independent
confirmations. BIC supplies the monotone component-count result; the others are
sensitivity checks with their own limitations.

A third control targets a constant number of VAE gradient updates rather than
a constant number of epochs. Full-epoch rounding adds at most 32 updates (1.1%)
and is disclosed in the manuscript; it does not scale monotonically with data.

An IID equilibrium control samples the exact synthetic latent coordinates from
gridded Boltzmann distributions. It removes trajectory autocorrelation and
projection learning; the BIC trend remains (+0.92 on Müller–Brown and +0.86 on
Prinz), directly isolating an estimator-splitting contribution.

## Layout

```
src/
  potentials.py   Müller-Brown and Prinz potentials; ground-truth basin assignment
  simulate.py     overdamped Langevin integrator; nonlinear lift to R^30
  embed.py        PCA, TICA (explicit generalised eigenproblem), VAE, t-SNE
  selection.py    five model-selection criteria for "how many states?"
  cluster.py      clustering + chance-corrected recovery metrics
  sweep.py        the experiment driver
  analyze.py      bootstrap CIs over seeds; Spearman budget-trend tests
  figures.py         all publication figures
  figures_alanine.py alanine-specific panels (Ramachandran reference;
                     TICA effective-lag association)
  alanine.py         real-system validation (alanine dipeptide, optional)
paper/
  manuscript.md   the full manuscript (markdown source)
  main.tex        submission LaTeX, generated from manuscript.md
  refs.bib        29 references, built from publisher metadata
  tools/          one-way markdown -> LaTeX pipeline and its checks
  outline.md      section-by-section plan mapped to figures
  threats.md      threats to validity and how each is addressed
results/
  NOTE.md         inventory: which sweep file is authoritative, and why
  iid_equilibrium_control.csv  exact-latent IID Boltzmann control
```

## Reproducing

```bash
pip install -r requirements-lock.txt  # recorded analysis environment
# or: pip install -r requirements.txt # supported minimum versions

# fresh outputs; does not silently resume from committed result files
./run_all.sh

# main sweep: ~3-6 hours on CPU at these settings
python src/sweep.py \
  --methods pca tica vae \
  --budgets 250 500 1000 2000 4000 8000 16000 32000 \
  --seeds 20 --modes short subsample \
  --out results/sweep.csv

# second landscape, for generality
python src/sweep.py --potential prinz1d --seeds 20 \
  --out results/sweep_prinz.csv

# figures
python src/figures.py --results results/sweep.csv
python src/figures.py --results results/sweep_prinz.csv \
  --outdir figures_prinz --potential prinz1d
```

`--resume` skips conditions already written, so a long run can be interrupted
and restarted safely.

### Optional real-system validation

```bash
pip install -r requirements-md-lock.txt  # recorded MD environment
# or: pip install -r requirements-md.txt # supported minimum versions

# 100 ns; ~1 h on a GPU, overnight on CPU
python src/alanine.py simulate --ns 100

python src/alanine.py sweep --dcd data/alanine/traj_seed0.dcd \
                            --top data/alanine_dipeptide.pdb \
                            --budgets 100 250 500 1000 2000 4000 \
                                      8000 16000 32000 64000 \
                            --out results/alanine_sweep_100ns.csv

# figures. --dcd is optional and only adds the Ramachandran ground-truth panel
python src/figures_alanine.py --results results/alanine_sweep_100ns.csv \
                              --dcd data/alanine/traj_seed0.dcd \
                              --top data/alanine_dipeptide.pdb
```

The 100 ns length is not arbitrary. An initial 5 ns pilot visited the rarest
basin (alpha_L/C7ax) inconsistently at the sampling budgets under study, which
confounded the alanine result with a basin-visitation deficit instead of
isolating the sample-size effect. The pilot sweeps are kept in `results/` for
provenance only; see `results/NOTE.md`.

## Statistical conventions

- The **seed is the replicate.** Frames within a trajectory are heavily
  autocorrelated, so all confidence intervals bootstrap over seeds, never over
  frames.
- The headline statistic is a **Spearman rank correlation** between
  `log(n_frames)` and reported state count, because the claim is monotone
  inflation rather than any particular functional form.
- **ARI** is the recovery metric because it is chance-corrected: it does not
  reward a method merely for producing more clusters, which is precisely the
  failure mode under investigation.
- Ceiling hits (`k_selected == kmax`) are reported explicitly. Where the rate is
  non-negligible, the inflation is right-censored and therefore understated.

## Status

Both synthetic sweeps are complete, along with the TICA lag and
embedding-dimension ablations, and the alanine dipeptide validation now runs on
two independent 100 ns trajectories rather than one. The manuscript is in
`paper/manuscript.md`, and `paper/main.pdf` is the built submission version.

The cross-seed comparison largely replicates: no budget-inflation correlation
moves by more than 0.02, and the method ranking is unchanged. The exception is
TICA, whose high-recovery crossover occurs one budget doubling later on seed 1 and which
therefore shows no decline within the tested range. Section 7 of the manuscript
reports this rather than averaging it away.

## Building the paper

```bash
python paper/tools/build.py                  # manuscript.md -> main.tex, then check
python paper/tools/build.py --pdf            # ...and build the PDF
python paper/tools/build.py --audit --pdf    # ...and re-verify all 29 references
python paper/tools/build.py --bib            # rebuild refs.bib from Crossref, then audit
```

`--pdf` runs pdflatex, bibtex and two more pdflatex passes, then reports
undefined references and overfull boxes. It locates the TeX binaries itself,
including MiKTeX's per-user install directory, which is often not on PATH.

`build.py` converts, then verifies: every numeric token in the markdown must
survive into the LaTeX, every `\cite` key must exist in `refs.bib`, every
`\ref` must resolve, and every `\includegraphics` path must exist.

`--audit` additionally checks all 29 references against OpenAlex and Semantic
Scholar. Those are independent of Crossref, which is what `refs.bib` is built
from, so the check is not circular. Five known disagreements are adjudicated in
`audit_refs.py`, each with its reason and the source consulted: four are
Semantic Scholar dating a preprint rather than the published paper, and one is
an OpenAlex page-range error contradicted by the publisher of record. Anything
else is reported as a regression. The
conversion is one-way -- once you begin editing `main.tex` directly, retire
`manuscript.md` rather than editing both.

## Citing

See `CITATION.cff`. Preprint: [https://doi.org/10.26434/chemrxiv.15008954/v1](https://doi.org/10.26434/chemrxiv.15008954/v1)

## License

Code is released under the MIT License (`LICENSE`). The manuscript text and
figures in `paper/` and `figures*/` are released under CC BY 4.0.
