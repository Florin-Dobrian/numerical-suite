# CLAUDE.md

Conventions for `numerical-suite`. Read this before writing anything in this repository.

---

## What this repository is

A working collection on **numerical algorithms and methods** — implementations, experiments, and reference documents. The unifying theme is not a single project but a subject: how numerical methods actually behave, established by running them rather than by citing them.

Two kinds of content, and they are held to different standards:

| | what it is | standard |
|---|---|---|
| **experiments** | scripts and notebooks exploring something | works, is readable, no further obligation |
| **reference documents** | prose meant to be re-read and relied on | everything in "Writing standard" below |

Most of the repository is currently the first kind. `methods/` is the second.

---

## Layout

```
numerical-suite/
├── methods/              reference documents on numerical methods
│   ├── COMPUTATIONAL.md      equations and solvers by field
│   ├── DISCRETIZATION.md     finite difference / volume / element
│   ├── TIME_INTEGRATION.md   explicit against implicit
│   ├── code/                 scripts that produce every number and figure
│   └── figs/                 their output, committed
├── experiments/          scratch work
├── *.py, *.ipynb         loose experiments at the root, by topic:
│                           bt*            Bradley–Terry ranking
│                           doubleml_*     causal inference
│                           paraboloid, rosenbrock, himmelblau,
│                           saddlepoint    optimization test functions
├── pyproject.toml        uv-managed; numpy, scipy, matplotlib, pandas,
│                         pyspark, doubleml, xgboost, lightgbm
└── uv.lock
```

**The root is unorganized and known to be.** When a topic there accumulates enough to be worth documenting, it should move into its own directory following the `methods/` pattern — document, `code/`, `figs/`. Optimization is the obvious next candidate: `paraboloid`, `rosenbrock`, `himmelblau` and `saddlepoint` are already a coherent set, and `COMPUTATIONAL.md` Part IV is the document they belong under.

---

## Environment and mechanics

- **Run with `uv`**: `uv run script.py`. Not `uv run python script.py` — uv dispatches `.py` files directly.
- **Dependencies live in `pyproject.toml`.** Check before adding: numpy, scipy, matplotlib, pandas and pyspark are already there. Do not add a dependency a script does not need — JAX was considered for the time-integration work and rejected, because the implicit side needs `scipy.linalg.solve_banded` and nothing needed autodiff.
- **All paths inside scripts must be relative to the script's own location**, via `os.path.dirname(os.path.abspath(__file__))`. Never absolute. A script must run correctly from any working directory.
- **Figures are committed.** The root `.gitignore` has a global `*.png`; `methods/figs/*.png` is re-included by negation because those figures are documentation, not output. Any new documented directory needs the same negation.
- **Intermediates are not committed.** JSON dumps passed between scripts are ignored.

---

## Writing standard

This is the part that matters, and it was arrived at the hard way. **A reference document is judged by whether someone who does not already know the material can follow it**, not by whether the statements in it are true.

**Define before use.** Every symbol gets a definition at or before its first appearance. No symbol arrives mid-argument. If a letter is already taken elsewhere in the document, pick another and say why.

**Motivate before asserting.** State the question a section answers before answering it. A reader who does not know why a thing is being computed cannot evaluate the result.

**No forward references in the argument.** Pointing ahead for *more detail* is fine; relying on something not yet established is not. If §2 needs a fact, §1 supplies it.

**Derive rather than jump.** Show the intermediate steps. If two derivatives are being discretized, discretize them one at a time and show both.

**Every number traceable.** No figure appears in prose that is not produced by a script in `code/` and listed in the document's run table. Numbers recalled rather than computed are marked `[unverified]` and collected in a dedicated section.

**State the consequence, not just the fact.** A measurement without its implication is a fact thrown at the reader. "The fastest mode has a lifetime of $7\times10^{-6}$" is incomplete; what it costs is the point.

**Repeat context across parallel sections.** If four sections are meant to be compared, they must answer the same questions in the same order with the same vocabulary. Introducing new terminology per section makes comparison impossible.

**Tables carry the argument; prose says what to read from them.** Do not restate a table in sentences. Do say which row matters and why.

**Mark what is not established.** `[unverified]` inline, collected in a closing section, with a note on what would settle it. Caveats belong in their own section, stated plainly, including the ones that weaken the document's own claims.

---

## Notation used in `methods/`

Consistent across the three documents. Reuse these; do not invent alternatives.

| symbol | meaning |
|---|---|
| $\mathbf{L}$ | the discrete spatial operator, after the method of lines |
| $\mathbf{L}_0$ | the same with physical factors stripped out, e.g. $\mathrm{tridiag}(1,-2,1)$ |
| $\lambda_m$ | eigenvalue of mode $m$; $M$ is the mode count, $m$ the only index |
| $\phi_m = k_m\Delta x$ | phase shift per cell. **Not $\theta$** — that is the method parameter |
| $z = \lambda\Delta t$ | the scaled eigenvalue; the only form in which $\lambda$ enters a stability statement |
| $g(z)$ | amplification factor: what one step does to one mode |
| $\theta$ | the $\theta$-method parameter: 0 forward Euler, ½ Crank–Nicolson, 1 backward Euler |
| $E_n$, $\ell_n$ | global and local error. **Subscripts**, so they cannot read as powers of $e$ |
| $r = \alpha\Delta t/\Delta x^2$ | Fourier number (diffusion) |
| $\nu = c\Delta t/\Delta x$ | Courant number (advection, waves) |
| $\rho$ | cost of one implicit solve divided by one explicit step |
| $S$ | safety factor on a stability limit, 0.8–0.9 |
| $p$, $q$ | order in time, order in space |
| $i$, $j$ | structured-grid indices: $i$ with $x$, $j$ with $y$. In a mesh-shaped display $i$ therefore runs across the **columns** — the transpose of matrix convention. Say so wherever an array is laid out on the mesh |
| $\mathbf{A}$ | the discrete operator of a steady problem, $\mathbf{A}\mathbf{T} \approx -\nabla^2 T$. Related to $\mathbf{L}$ but not the same: $\mathbf{L}$ is the method-of-lines right-hand side, so it carries the opposite sign |
| $A_e$ | area of one element. **Subscripted**, so it cannot be read as the system matrix $\mathbf{A}$ set in italic |
| $I$, $J$ | global degree-of-freedom index: the flattened node or cell number, as in $K_{IJ}$, $F_I$. **Not $p$** — that is the order of accuracy. Capital, so it cannot be read as a grid index |

---

## Document structure in `methods/`

Each of the three documents follows the same shape, and a new one should:

1. **Status block** — what belongs here, scope, what is borrowed, how claims are marked.
2. **A plan section** — what is covered and why *these* examples, before any machinery.
3. **Machinery**, built in dependency order.
4. **Worked sections with identical skeletons.** `TIME_INTEGRATION.md` uses five questions per section: where are the modes, where is the stability ceiling, what step does accuracy need, how big is the gap, what happens to the unresolved modes.
5. **A runs table** — script, sections, what it computes.
6. **Caveats**, **unverified claims**, **sources**.

Cross-references between documents are by section number and are expected to be maintained.

---

## Verification discipline

- **Check before stating.** If a claim can be computed, compute it before writing it down.
- **Report predicted against measured** wherever theory gives a prediction. Agreement is the result; disagreement is more interesting still.
- **Timings are indicative, step counts are results.** Everything here is Python, and wall-clock numbers are dominated by interpreter overhead on small problems. Conclusions should rest on operation counts.
- **State what the measurement does not cover.** One-dimensional experiments do not establish multidimensional behaviour; a single manufactured solution does not establish general convergence.
