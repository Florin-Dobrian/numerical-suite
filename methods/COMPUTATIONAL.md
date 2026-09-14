# Computational methods — an inventory of equations and the solvers they force

**Field by field: what is being solved, what mathematical character it has, and why that character leaves you with the solver you end up with.**

---

> **What belongs here.** One section per computational field. Each section carries the governing equations, the reductions in common use, and a methodology discussion that runs from the PDE to the linear algebra. The organising claim is that **the solver is largely not a choice** — it is forced, in stages, by the character of the equations and by the question being asked of them.
>
> **Scope.** Two fields: **structural mechanics** (Part I) and **fluid dynamics** (Part II). Part III sets them side by side, which is where most of the value is — the contrasts are sharper than either section alone. **Part IV is optimization**, which is not a third field but a **layer** that wraps one: the simulation becomes a constraint, and everything above it becomes the inner problem. §12 says what is still missing and applies that distinction to it.
>
> **Borrowed material.** The fluid equations in §3.2–§3.3 and the PDE-character framing in §0.3 are condensed from the `fusion` repository — `notes/theory/mhd_model.md` §1 and `notes/theory/pde_character.md`. Those documents are the fuller treatment and should be read for the derivations; what is here is the part that generalises beyond plasma.
>
> **Marking.** Claims that could not be checked against a source or a run are marked **[unverified]** and collected in §13. Order-of-magnitude figures carry their formula so they can be rechecked.

---

# Part 0 — The chain, and the axes

## 0.1 From physics to solver, in five steps

Nothing in this document is a single decision. It is a chain, and each link constrains the next without determining it:

| | Step | What is chosen | What constrains the choice |
|---|---|---|---|
| 1 | **Continuum statement** | balance laws, kinematics, constitutive law | the physics, and what you are willing to neglect |
| 2 | **PDE system** | the closed set of equations | step 1, plus the closure |
| 3 | **Character** | hyperbolic / parabolic / elliptic | **not chosen** — read off the system |
| 4 | **Discretization** | FEM, FVM, FD, spectral; and where the time level sits | step 3, plus geometry and conservation requirements |
| 5 | **Algebraic system** | linear or nonlinear, symmetric or not, its size and sparsity | steps 2 and 4 |
| 6 | **Solver** | direct or iterative; Newton or Picard; which preconditioner | step 5 |

**Step 3 is the pivot.** It is the only one that is not a choice, and it is the one that removes options from every step below it. A great deal of confusion about numerical methods comes from arguing about step 6 while the constraint actually lives at step 3.

## 0.2 The four axes, and what is independent of what

Four binary distinctions get used interchangeably in conversation and they are not the same distinction.

| Axis | What it is about | Decided by |
|---|---|---|
| **Hyperbolic / parabolic / elliptic** | how information propagates | the PDE system, at a given state |
| **Linear / nonlinear** | whether the operator depends on the solution | the constitutive law, the kinematics, or the convective term |
| **Explicit / implicit** | where the right-hand side is evaluated in time | a discretization choice — where one exists |
| **Direct / iterative** | how the resulting linear system is solved | size, sparsity, symmetry, and how accurately it must be solved |

**They are close to orthogonal, and the four most common confusions are these:**

- **Implicit does not mean nonlinear.** Backward Euler on the heat equation is implicit and perfectly linear.
- **Nonlinear does not mean implicit.** Explicit crash codes handle plasticity, contact and fracture — the most violently nonlinear problems in engineering — without ever forming a tangent or solving a system.
- **Coupling does not force implicit.** Tightly coupled systems can be fully explicit whenever the coupling is algebraic rather than global. What forces implicit treatment is an *elliptic constraint* or *stiffness*, and neither is the same as coupling.
- **Elliptic problems are not "implicit because implicit is better."** There is no time level to be explicit about. The word does not apply.

The last of these is the one worth internalising: **explicit and implicit are only meaningful where a time derivative exists.**

## 0.3 How character is determined

### The textbook version

For a second-order scalar PDE in two variables,

$$A\,\phi_{xx} + B\,\phi_{xy} + C\,\phi_{yy} + \text{lower order} = 0$$

the classification is by the discriminant $B^2 - 4AC$: negative is elliptic, zero parabolic, positive hyperbolic. This is correct, it is what every textbook opens with, and it is close to useless in practice — most systems of interest are first-order systems of several equations, not one second-order scalar equation.

### The version that gets used

For a first-order system $\partial_t \mathbf{q} + \mathbf{A}(\mathbf{q})\,\partial_x \mathbf{q} = 0$, the character is read off the eigenvalues of the flux Jacobian $\mathbf{A}$:

| Eigenvalues of $\mathbf{A}$ | Character |
|---|---|
| real and a complete set of eigenvectors | **hyperbolic** — the eigenvalues are the wave speeds |
| real and distinct | strictly hyperbolic |
| complex | **elliptic** |
| some zero, second-order terms present | **parabolic** |

The eigenvalues *are* the characteristic speeds, which is why this version is useful: it tells you not only the type but the numbers, and those numbers are what set the timestep and the boundary conditions.

### What each type means

| | information travels | domain of dependence | can it be marched? |
|---|---|---|---|
| **Hyperbolic** | at finite speed, along characteristics | a cone in space–time | **yes** — explicit, $\Delta t \lesssim \Delta x/c$ |
| **Parabolic** | instantly, but damped with distance | everything, with rapid decay | **yes** — but $\Delta t \lesssim \Delta x^2/2D$ |
| **Elliptic** | no time direction at all | everything, undamped | **no** — global solve required |

**The rule of thumb that follows:** first-order derivatives in space with a time derivative give propagation; second-order spatial derivatives give diffusion; no time derivative at all gives an equilibrium problem. Odd-order spatial derivatives are directional and produce real frequencies; even-order ones are symmetric and produce imaginary frequencies, hence decay.

### Boundary and initial conditions follow from the type

This is the practical payoff of the classification, and it is often the first place a badly posed problem announces itself.

| | What is well posed |
|---|---|
| **Elliptic** | boundary conditions on the **entire closed boundary**, and nothing else. No initial condition; there is no time. |
| **Parabolic** | an initial condition, plus boundary conditions on the whole spatial boundary for all time |
| **Hyperbolic** | an initial condition, plus **as many boundary conditions at each boundary point as there are characteristics entering the domain there** |

The hyperbolic rule is the one that bites. At a supersonic inflow every characteristic enters, so every variable is prescribed. At a subsonic outflow one characteristic still runs upstream, so exactly one condition — usually pressure — may be imposed and no more. Imposing two overdetermines the problem and produces reflections; imposing none underdetermines it.

### The type is a property of the state, not of the equation's name

A system can change type as the solution changes. The steady small-disturbance potential equation,

$$(1 - M^2)\,\phi_{xx} + \phi_{yy} = 0$$

is elliptic for $M < 1$ and hyperbolic for $M > 1$, so a transonic flow field is **mixed** — elliptic in the subsonic pocket, hyperbolic outside it, with the boundary between them part of the unknown. That is why transonic flow was historically hard and why it needed schemes that switch behaviour locally.

**A nuance worth keeping straight:** it is the *steady* equations that change type. The unsteady compressible Euler equations are hyperbolic at every Mach number, because the time derivative supplies a characteristic direction regardless. Steady-state solvers that march in pseudo-time are exploiting exactly this.

### A fourth category, in passing

Some systems are neither hyperbolic nor parabolic: **dispersive** ones, where $\omega$ is real but nonlinear in $k$, so waves propagate and small scales oscillate ever faster without damping. The whistler wave of Hall physics is the standard example, with $\omega \propto k^2$ giving a $\Delta t \propto \Delta x^2$ restriction for the *opposite* reason diffusion does. Neither structural mechanics nor ordinary fluid dynamics produces this, so it is noted and set aside; `fusion`'s `notes/theory/pde_character.md` §4 treats it.

## 0.4 What actually forces an implicit method

Two things, and the list is short.

**1. An elliptic constraint.** There is no time direction, so there is nothing to march. The pressure in incompressible flow is the canonical case (§3.3). A global solve is unavoidable, and it is unavoidable per timestep, not once.

**2. Stiffness.** A time direction exists, but the fastest mode the grid supports is far faster than the physics of interest, and an explicit method must still resolve it or blow up. Formally, for $\mathbf{u}' = \mathbf{A}\mathbf{u}$, stiffness is the spread of eigenvalues $S = |\lambda_{\max}|/|\lambda_{\min}|$; practically, it is the ratio of the fastest resolvable timescale to the timescale you actually care about.

**Stiffness is not a property of the equations alone.** It is a property of the equations, the grid, **and the question**. Refine the mesh and a well-behaved problem becomes stiff. Change the question — decide you want to see the fast mode after all — and a stiff problem stops being stiff, because nothing is being wasted any more.

**This is the whole answer to "why do structures go implicit and fluids often not."** It is worked out in §2.4 and §4.4 and summarised in §5.2.

## 0.5 The trade, in one table

| | matrix operation | cost per step | step restriction |
|---|---|---|---|
| **Explicit** | **multiply** by $\mathbf{A}$ | cheap, local, parallel-friendly | $\Delta t \lesssim \Delta x/c$ or $\Delta x^2/2D$ |
| **Implicit** | **solve** with $\mathbf{A}$ | expensive, global | none from stability; accuracy still binds |
| **Elliptic** | **solve** with $\mathbf{A}$ | expensive, global | no time to restrict |

Explicit multiplies; implicit solves. Elliptic problems leave no choice.

**One caveat that costs people real time.** Implicit does not mean unconditionally stable — backward Euler on diffusion happens to be, but implicit schemes in general are not. And unconditionally stable does not mean accurate. A very large implicit step is stable and *wrong*, and it fails quietly rather than blowing up, which is the more dangerous failure of the two.

---

# Part I — Structural mechanics

# 1. The equations

## 1.1 The three-part structure

Every problem in solid mechanics is built from exactly three ingredients, and knowing which is which tells you where the difficulty will come from.

| | What it says | Universal? |
|---|---|---|
| **Balance** | Newton's second law for a continuum: $\nabla\cdot\boldsymbol\sigma + \mathbf{b} = \rho\,\ddot{\mathbf{u}}$ | **yes** — holds for any material whatsoever |
| **Kinematics** | how strain is defined from displacement | **yes**, but the definition is a modelling choice (small vs finite) |
| **Constitutive** | how stress responds to strain: $\boldsymbol\sigma = \boldsymbol\sigma(\boldsymbol\varepsilon, \dots)$ | **no** — this is the material |

**The three nonlinearities of structural mechanics map one-to-one onto these three slots:**

| Nonlinearity | Enters through | Example |
|---|---|---|
| **Geometric** | kinematics — strain is quadratic in displacement | a fishing rod, a buckling column, a membrane |
| **Material** | constitutive — stress is not linear in strain | plasticity, hyperelasticity, creep, damage |
| **Boundary** | the boundary conditions themselves change | contact, friction, follower loads |

That taxonomy is standard and it is worth stating early, because it is also a taxonomy of *how the solver will misbehave*. Geometric nonlinearity gives limit points and needs arc-length control. Material nonlinearity gives path dependence and needs history variables and a consistent tangent. Boundary nonlinearity gives a non-smooth problem where Newton's method has no derivative to work with, and it is the reason explicit codes exist in this field at all.

## 1.2 Balance — and the same equation as fluids

$$\nabla\cdot\boldsymbol\sigma + \mathbf{b} = \rho\,\ddot{\mathbf{u}}$$

with $\boldsymbol\sigma$ the Cauchy stress tensor, $\mathbf{b}$ the body force per unit volume, and $\mathbf{u}$ the displacement.

**This is the Cauchy momentum equation, and it is identical to the one the Navier–Stokes equations are built on** (§3.1). The difference between a solid and a fluid is not in the balance law at all: it is entirely in the constitutive relation supplied for $\boldsymbol\sigma$, and in whether the frame is attached to the material or to space.

That single observation organises most of Part III.

### The three regimes

Dropping terms from the right-hand side gives the three problems structural analysis actually solves:

| | Equation | Character | When |
|---|---|---|---|
| **Static** | $\nabla\cdot\boldsymbol\sigma + \mathbf{b} = 0$ | **elliptic** | inertia negligible; the answer is an equilibrium |
| **Quasi-static** | $\nabla\cdot\boldsymbol\sigma + \mathbf{b} = 0$, with load a function of a parameter | elliptic at each load step | slow loading, but path-dependent material or geometry |
| **Dynamic** | $\nabla\cdot\boldsymbol\sigma + \mathbf{b} = \rho\ddot{\mathbf{u}}$ | **hyperbolic** | inertia matters — impact, vibration, wave propagation |

**Quasi-static is not a fourth kind of physics.** It is a sequence of elliptic problems marched along a load parameter, and the parameter exists only because the material or the geometry remembers the path taken. If neither does, the intermediate steps are wasted work and one solve at the final load is enough.

## 1.3 Kinematics — small strain and finite strain

**Small strain** (infinitesimal, engineering):

$$\boldsymbol\varepsilon = \frac{1}{2}\left(\nabla\mathbf{u} + (\nabla\mathbf{u})^{\mathsf{T}}\right)$$

Linear in $\mathbf{u}$. This is the assumption that makes linear elasticity linear, and it requires both small strains **and** small rotations. A structure that rotates rigidly through a large angle while straining by $10^{-4}$ violates it.

**Finite strain** keeps the quadratic term. With the deformation gradient $\mathbf{F} = \mathbf{I} + \partial\mathbf{u}/\partial\mathbf{X}$, the Green–Lagrange strain is

$$\mathbf{E} = \frac{1}{2}\left(\mathbf{F}^{\mathsf{T}}\mathbf{F} - \mathbf{I}\right) = \frac{1}{2}\left(\nabla_0\mathbf{u} + (\nabla_0\mathbf{u})^{\mathsf{T}} + (\nabla_0\mathbf{u})^{\mathsf{T}}\nabla_0\mathbf{u}\right)$$

The final term is the geometric nonlinearity, and it is the entire difference. $\mathbf{E}$ is objective — it vanishes under rigid rotation, which $\boldsymbol\varepsilon$ does not.

### Reference configuration versus current configuration

Solid mechanics is naturally **Lagrangian**: the mesh is attached to the material and moves with it, and $\mathbf{u}$ is the displacement of a material point from where it started. Fluid mechanics is naturally **Eulerian**: the mesh is fixed in space and material flows through it.

This choice has consequences far beyond bookkeeping, and it explains two things at once:

- **Why fluids have a convective term and solids do not.** $(\mathbf{u}\cdot\nabla)\mathbf{u}$ exists because the frame is fixed and material moves past it. In a Lagrangian frame, the material derivative *is* the partial derivative, and the term is absent.
- **Why solids have geometric nonlinearity and fluids mostly do not.** The domain itself is an unknown when the mesh follows the material.

Two stress measures follow from the two frames, and mixing them is a standard source of error: **Cauchy stress** $\boldsymbol\sigma$ is force per unit *current* area (what a gauge measures), while the **second Piola–Kirchhoff stress** $\mathbf{S}$ is referred to the undeformed configuration and is the work conjugate of $\mathbf{E}$. They coincide only in the small-strain limit.

## 1.4 Constitutive laws — the ladder

This is where the material lives, and where nearly all the modelling effort goes.

### Linear elastic

$$\boldsymbol\sigma = \mathbb{C} : \boldsymbol\varepsilon$$

with $\mathbb{C}$ the fourth-order stiffness tensor. In general anisotropy $\mathbb{C}$ has **21 independent constants**; orthotropy reduces this to 9, transverse isotropy to 5, and **isotropy to 2**. For an isotropic material,

$$\boldsymbol\sigma = \lambda\,\mathrm{tr}(\boldsymbol\varepsilon)\,\mathbf{I} + 2\mu\,\boldsymbol\varepsilon$$

with $\lambda$ and $\mu$ the Lamé parameters, $\mu$ the shear modulus, and the engineering pair $(E, \nu)$ related by $\mu = E/2(1+\nu)$ and $\lambda = E\nu/[(1+\nu)(1-2\nu)]$.

**Note the $(1-2\nu)$ in the denominator of $\lambda$.** As $\nu \to 1/2$ the material becomes incompressible and $\lambda \to \infty$. That singularity is not cosmetic: it is the source of volumetric locking, it makes the stiffness matrix catastrophically ill-conditioned, and the standard fix — a mixed displacement–pressure formulation — produces exactly the same saddle-point system as incompressible flow. **Rubber and incompressible Navier–Stokes are the same numerical problem in this respect** (§2.6, §4.3).

### Hyperelastic

For large elastic strains, stress is derived from a **strain energy density** $W(\mathbf{F})$:

$$\mathbf{S} = \frac{\partial W}{\partial \mathbf{E}}$$

Common forms: **Neo-Hookean**, $W = C_1(\bar{I}_1 - 3) + \frac{1}{D_1}(J-1)^2$; **Mooney–Rivlin**, adding an $\bar{I}_2$ term; **Ogden**, in principal stretches. All are nonlinear but **path-independent and reversible** — the material has no memory, which keeps them well behaved. Energy is a function of state, so a potential exists and the tangent stays symmetric.

### Elastoplastic

The first genuinely difficult case, because it is **path dependent**. The ingredients:

| | |
|---|---|
| **Yield function** | $f(\boldsymbol\sigma, q) \le 0$ — elastic inside, plastic on the surface. Von Mises: $f = \sqrt{3J_2} - \sigma_y(q)$ |
| **Flow rule** | $\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\partial f/\partial\boldsymbol\sigma$ — associated; non-associated when the plastic flow direction differs (soils, some metals under pressure) |
| **Hardening** | how $\sigma_y$ evolves: isotropic (surface grows), kinematic (surface translates — needed for cyclic loading and the Bauschinger effect), or both |
| **Consistency** | $\dot\gamma \ge 0$, $f \le 0$, $\dot\gamma f = 0$ — a complementarity condition, not an equation |

**Three consequences for the solver:**

1. **History variables at every integration point.** Plastic strain and hardening variables are stored per Gauss point and updated each step. The state of the model is far larger than the displacement vector.
2. **The constitutive update is itself a local nonlinear solve.** The standard algorithm is **return mapping**: take an elastic trial step, and if it lands outside the yield surface, project back onto it. For von Mises with isotropic hardening this reduces to a scalar nonlinear equation per point — the radial return.
3. **Non-associated flow destroys symmetry.** The tangent is no longer symmetric, and a Cholesky factorization no longer applies. This is a real cost in geomechanics.

### Viscoelastic and creep

Rate dependence: stress depends on strain *history*, not just strain. Represented by spring–dashpot networks (Maxwell, Kelvin–Voigt, generalized Maxwell) and implemented as a **Prony series** with internal variables per relaxation time. Adds a real timescale to a problem that might otherwise have none, which converts a static problem into a parabolic one.

### Damage and fracture

Stiffness degrades with a damage variable $d$: $\boldsymbol\sigma = (1-d)\,\mathbb{C}:\boldsymbol\varepsilon$. **Softening** — the tangent becomes negative definite — which makes the problem locally ill-posed and mesh-dependent unless regularised (nonlocal, gradient, or phase-field formulations). The most numerically hostile constitutive class in routine use.

### Summary — what each costs

| Constitutive class | Nonlinear? | Path dependent? | Tangent symmetric? | Solver consequence |
|---|---|---|---|---|
| Linear elastic | no | no | **yes**, SPD | one linear solve |
| Hyperelastic | yes | no | yes | Newton, well behaved |
| Plasticity, associated | yes | **yes** | yes (consistent tangent) | Newton + local return mapping + history storage |
| Plasticity, non-associated | yes | yes | **no** | Newton with a nonsymmetric solver |
| Viscoelastic | mildly | yes | yes | internal variables, real timescale introduced |
| Damage / softening | yes | yes | **indefinite** | needs regularisation and arc-length control |

## 1.5 The canonical problems and their character

Substituting Hooke's law and small strain into the balance equation gives the **Navier–Cauchy** (or Lamé) equation:

$$(\lambda + \mu)\,\nabla(\nabla\cdot\mathbf{u}) + \mu\,\nabla^2\mathbf{u} + \mathbf{b} = \rho\,\ddot{\mathbf{u}}$$

Everything below is a specialisation of this.

| Problem | Equation | Character | What the solver does |
|---|---|---|---|
| **Static linear elasticity** | drop $\rho\ddot{\mathbf{u}}$ | **elliptic** | one linear solve, $\mathbf{K}\mathbf{u} = \mathbf{f}$ |
| **Elastodynamics** | keep it | **hyperbolic** | march, explicitly or implicitly |
| **Modal analysis** | $(\mathbf{K} - \omega^2\mathbf{M})\boldsymbol\phi = 0$ | eigenproblem | Lanczos / subspace iteration |
| **Linear buckling** | $(\mathbf{K} + \lambda\mathbf{K}_\sigma)\boldsymbol\phi = 0$ | eigenproblem | same machinery, different matrices |
| **Heat conduction / thermal** | $\rho c\,\dot{T} = \nabla\cdot(k\nabla T) + Q$ | **parabolic** | implicit almost always |
| **Harmonic / frequency response** | $(\mathbf{K} - \omega^2\mathbf{M} + i\omega\mathbf{C})\hat{\mathbf{u}} = \hat{\mathbf{f}}$ | complex linear system per frequency | direct solve per frequency, or modal superposition |

**Two of these are not PDE solves at all.** Modal analysis and buckling are eigenvalue problems, and the numerical methods they need — Lanczos, shift-and-invert, subspace iteration — have nothing to do with time integration or with the explicit/implicit question. This is a large fraction of what structural analysis software actually spends its time on, and it is invisible if you think of the field only as "solve a PDE."

**Elastodynamics carries two wave speeds**, which is what sets the explicit timestep in §2.3:

$$c_p = \sqrt{\frac{\lambda + 2\mu}{\rho}}, \qquad c_s = \sqrt{\frac{\mu}{\rho}}$$

the dilatational (pressure) and shear speeds. For structural steel with $E = 210$ GPa, $\nu = 0.3$, $\rho = 7850$ kg/m³: $\lambda + 2\mu = 283$ GPa, so $c_p \approx 6000$ m/s and $c_s \approx 3200$ m/s. **The larger one is what binds.**

## 1.6 Reduced models — beams, plates, shells

Structures are usually thin, and exploiting that is where structural mechanics diverges most from fluid mechanics. A thin body's through-thickness behaviour can be integrated out analytically, leaving a lower-dimensional problem with a **higher-order** operator.

| Model | Equation | Order | Assumption dropped |
|---|---|---|---|
| **Euler–Bernoulli beam** | $EI\,w'''' = q$ | 4th | shear deformation, rotary inertia |
| **Timoshenko beam** | coupled $w$, $\theta$ system | 2nd, coupled | — restores shear and rotary inertia |
| **Kirchhoff–Love plate** | $D\,\nabla^4 w = q$, $D = Eh^3/12(1-\nu^2)$ | 4th, biharmonic | transverse shear |
| **Reissner–Mindlin plate** | coupled $w$, $\boldsymbol\theta$ | 2nd, coupled | — restores transverse shear |
| **Shell** | membrane + bending, coupled | mixed | — |

All of these are still **elliptic** in statics. What changes is the numerical difficulty, in three specific ways:

**1. Fourth-order operators require $C^1$ continuity.** A conforming finite element for the Kirchhoff plate must have continuous first derivatives across element boundaries, which is genuinely awkward to construct. This is the main practical reason Reissner–Mindlin and Timoshenko formulations dominate in software despite being the more elaborate physics: their weak forms need only $C^0$ continuity.

**2. Conditioning degrades as $h^{-4}$, not $h^{-2}$.** A second-order elliptic operator gives a stiffness matrix with condition number $O(h^{-2})$; a fourth-order one gives $O(h^{-4})$. Refining a plate mesh by a factor of ten costs four orders of magnitude in conditioning, which is why iterative solvers struggle on shell models and why direct solvers persist there.

**3. Locking.** As thickness $h \to 0$, low-order elements become artificially stiff — **shear locking** in Timoshenko and Mindlin formulations, **membrane locking** in curved shells, **volumetric locking** at $\nu \to 1/2$. All three are the same phenomenon: a constraint that should be satisfied in a weak sense is enforced pointwise by an element that cannot represent the constrained solution. The fixes — reduced and selective integration, assumed strain, enhanced assumed strain, mixed formulations — are among the most-cited results in the FEM literature.

## 1.7 Contact — why it is a different kind of problem

Contact is not a hard nonlinearity. It is a **non-smooth** one, which is worse.

The Signorini conditions on a candidate contact surface, with $g$ the gap and $p_N$ the normal pressure:

$$g \ge 0, \qquad p_N \ge 0, \qquad g\,p_N = 0$$

This is a **complementarity condition**, and the resulting problem is a variational inequality rather than an equation. The consequences:

- **The residual is not differentiable** at the moment of contact. Newton's method is built on a Taylor expansion that does not exist there.
- **The active set changes discontinuously.** One node coming into contact changes the sparsity structure of the tangent.
- **Friction makes the tangent nonsymmetric**, and Coulomb friction adds a second complementarity condition on the tangential direction.

Standard treatments: **penalty** (allow small penetration, add a stiff spring — simple, but introduces a user-chosen parameter and worsens conditioning), **Lagrange multipliers** (exact, but adds unknowns and produces a saddle point), **augmented Lagrangian** (a compromise, and the usual default), and **mortar methods** for non-matching meshes.

**This is the single largest source of non-convergence in implicit structural analysis**, and it is the strongest argument for explicit methods in the impact and forming industries — an explicit code never asks Newton to converge, so a problem that is merely non-smooth costs it nothing.

---

# 2. Methodology — how a structural problem becomes algebra

## 2.1 Why the finite element method

FEM is close to universal in structural mechanics, and the reasons are specific rather than historical.

**1. The problem is elliptic and has a variational structure.** Static elasticity minimises potential energy. The weak form is the principle of virtual work:

$$\int_\Omega \boldsymbol\sigma : \delta\boldsymbol\varepsilon\,\mathrm{d}V = \int_\Omega \mathbf{b}\cdot\delta\mathbf{u}\,\mathrm{d}V + \int_{\Gamma} \mathbf{t}\cdot\delta\mathbf{u}\,\mathrm{d}S$$

For a self-adjoint positive-definite operator, the Galerkin solution is the **best approximation in the energy norm** — not merely a good one. Céa's lemma makes this precise. That optimality is a genuine mathematical advantage and it does *not* transfer to the advection-dominated operators of fluid dynamics (§4.1), which is the deepest reason the two fields discretize differently.

**2. Natural boundary conditions are tractions.** Integrating by parts moves one derivative onto the test function and produces a boundary term that *is* the applied traction. Loads — the thing an engineer actually specifies — enter the formulation for free, without needing to be imposed. In a strong-form method they have to be imposed explicitly and awkwardly.

**3. The weak form lowers the continuity requirement.** The strong form of elasticity needs $\mathbf{u}$ twice differentiable; the weak form needs it once. That is what permits simple piecewise-polynomial elements at all.

**4. Complex geometry, and the thin-structure problem.** Unstructured meshes of arbitrary shape, and dimensionally-reduced beam, plate and shell elements that carry the analytic through-thickness integration of §1.6 inside them. Neither finite differences nor finite volumes offer anything comparable.

**What FEM is not chosen for:** conservation. Finite volumes exist because conservation of a flux across a cell face must hold discretely, which matters enormously for shocks and not at all for a static stress field. Structures rarely need it.

## 2.2 What comes out

**Statics:**

$$\mathbf{K}\mathbf{u} = \mathbf{f}$$

with $\mathbf{K}$ the assembled stiffness matrix. Its properties, and where each comes from:

| Property | Origin | Lost when |
|---|---|---|
| **Sparse** | elements couple only neighbouring nodes | never |
| **Symmetric** | the operator is self-adjoint; a potential energy exists | non-associated plasticity, friction, follower loads, some contact |
| **Positive definite** | after essential BCs remove rigid-body modes | buckled or softening states — an indefinite $\mathbf{K}_T$ is the signature of a limit point |
| **Banded / profile** | node numbering, improvable by reordering | — |

**Before boundary conditions, $\mathbf{K}$ is singular**, with a null space of dimension 6 in 3D — three translations and three rotations. A model that will not solve because of a "zero pivot" is almost always one where a part is insufficiently restrained, and that diagnosis follows directly from the null space.

**Dynamics**, after semi-discretization in space only (the method of lines):

$$\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{K}\mathbf{u} = \mathbf{f}(t)$$

A system of second-order ODEs. Two mass matrix forms matter:

- **Consistent mass**, from the same shape functions as $\mathbf{K}$ — more accurate, but not diagonal.
- **Lumped mass**, mass concentrated at nodes — **diagonal**, hence trivially invertible.

That distinction is not a detail. **Lumped mass is what makes explicit dynamics matrix-free**, and the entire explicit branch of §2.3 depends on it.

$\mathbf{C}$ is the least physical object in the system. **Rayleigh damping**, $\mathbf{C} = \alpha\mathbf{M} + \beta\mathbf{K}$, is used overwhelmingly, and its justification is that it keeps $\mathbf{C}$ diagonalisable in the same basis as $\mathbf{M}$ and $\mathbf{K}$, not that any material behaves that way. **[unverified]** — that a mass-and-stiffness-proportional form is chosen for algebraic convenience rather than physical fidelity is standard commentary, but this document has not checked a primary source.

## 2.3 Time integration

### Implicit — the Newmark family

Newmark's scheme parameterises the update by $\beta$ and $\gamma$:

$$\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\dot{\mathbf{u}}_n + \frac{\Delta t^2}{2}\left[(1-2\beta)\ddot{\mathbf{u}}_n + 2\beta\,\ddot{\mathbf{u}}_{n+1}\right]$$

$$\dot{\mathbf{u}}_{n+1} = \dot{\mathbf{u}}_n + \Delta t\left[(1-\gamma)\ddot{\mathbf{u}}_n + \gamma\,\ddot{\mathbf{u}}_{n+1}\right]$$

Substituting into the equation of motion gives an **effective stiffness** system to be solved each step:

$$\left(\frac{1}{\beta\Delta t^2}\mathbf{M} + \frac{\gamma}{\beta\Delta t}\mathbf{C} + \mathbf{K}\right)\mathbf{u}_{n+1} = \mathbf{f}_{\text{eff}}$$

| $\gamma$, $\beta$ | Scheme | Properties |
|---|---|---|
| $\gamma = 1/2$, $\beta = 1/4$ | average acceleration / trapezoidal | **unconditionally stable**, second order, **no numerical damping** |
| $\gamma = 1/2$, $\beta = 0$ | central difference | explicit, conditionally stable |
| $\gamma > 1/2$ | — | adds numerical damping, drops to first order |

Unconditional stability requires $\gamma \ge 1/2$ and $\beta \ge (\gamma + 1/2)^2/4$.

**Why you often want numerical damping.** The highest modes of a finite element mesh are discretization artifacts, not physics — they are the modes the mesh cannot represent accurately anyway. Left undamped, they contribute noise. But raising $\gamma$ in plain Newmark damps them only at the cost of second-order accuracy everywhere. **HHT-$\alpha$** (Hilber–Hughes–Taylor) and the **generalized-$\alpha$** method solve exactly this: they achieve controllable high-frequency dissipation while remaining second-order accurate and unconditionally stable. Generalized-$\alpha$ is the default in most modern implicit structural codes, and the reason is this one property.

### Explicit — central difference

$$\ddot{\mathbf{u}}_n = \mathbf{M}^{-1}\left(\mathbf{f}_n - \mathbf{f}^{\text{int}}(\mathbf{u}_n)\right)$$

With lumped $\mathbf{M}$, the inverse is elementwise division. **No stiffness matrix is ever assembled and no system is ever solved.** Internal forces are computed element by element from the current stresses, which are obtained from the constitutive law at each integration point. Memory scales with the mesh; work per step is a fixed small multiple of the element count.

Stability requires

$$\Delta t \le \frac{2}{\omega_{\max}} \approx \frac{L_{\min}}{c_p}$$

the transit time of a dilatational wave across the smallest element. **The smallest element in the entire model sets the step for the entire model**, which is why mesh quality control in explicit codes is largely about eliminating small elements, and why mass scaling — artificially inflating the density of the offending elements — is a standard and slightly disreputable practice.

**The numbers.** For steel at $c_p = 6000$ m/s and a 1 mm element:

$$\Delta t \approx \frac{10^{-3}\ \text{m}}{6000\ \text{m/s}} = 1.7\times10^{-7}\ \text{s}$$

## 2.4 Why structural analysis so often goes implicit

The observation is correct, and it has three separate causes that are worth keeping apart.

### Cause 1 — most structural problems have no time in them at all

Statics, modal analysis, buckling, frequency response. **There is no time derivative to place at level $n$ or $n+1$**, so implicit is not chosen over explicit; the question does not arise. This is a large share of all structural analysis performed, and it is the least discussed reason.

A sharper version of the same point: **explicit integration requires a mass matrix to divide by.** A static problem has none. Explicit statics is not a slow method — it is not a method.

### Cause 2 — the timescale separation is enormous

Where time does exist, compare the event duration against the explicit step from §2.3:

| Application | Event duration | Explicit steps required at $\Delta t = 1.7\times10^{-7}$ s |
|---|---|---|
| Ballistic impact | 100 µs | ~600 |
| Vehicle crash | 100 ms | ~6 × 10⁵ |
| Earthquake response | 30 s | ~2 × 10⁸ |
| Quasi-static load | minutes to hours | absurd |
| Creep | months | absurd |

**The dilatational wave crosses a 1 mm element in 170 nanoseconds; a building sways in seconds.** That ratio — around $10^8$ — is the stiffness of the problem in the numerical sense of §0.4, and it is precisely what implicit integration exists to escape. Implicit steps are set by **accuracy**, meaning by the highest mode you want to *represent* — perhaps 20 steps per period of the highest mode of interest — not by the highest mode the mesh happens to support.

### Cause 3 — the linear algebra is favourable

$\mathbf{K}$ is sparse, symmetric and positive definite, and a sparse Cholesky factorization of it is fast and utterly reliable. The cost per implicit step is therefore lower relative to explicit than it would be in a field with nonsymmetric indefinite matrices. Structures can afford implicit in a way that many fields cannot (§2.6).

### And yet structures use explicit methods too

Crash, blast, ballistic impact, drop test, metal forming, machining. Two reasons, and the second is the one usually underweighted:

**1. The event is fast, so nothing is wasted.** At 100 µs duration, 600 explicit steps is not a penalty. **The criterion is universal: is the fastest resolvable mode close to the physics you care about? If yes, explicit is efficient.** That is exactly why fluids go explicit (§4.4) — the same test, a different answer.

**2. Explicit never asks Newton to converge.** Contact, fragmentation, self-contact, buckling collapse, material failure — these make implicit Newton iterations fail, and a failed Newton iteration means the analysis stops. An explicit code has no Newton iteration to fail. **It trades a guaranteed small step for a guaranteed step**, and for problems where robustness is the binding constraint rather than cost, that is the better trade.

That trade is the honest summary of the whole implicit/explicit question in structural mechanics: **implicit is chosen for efficiency, explicit for robustness**, and the timescale decides which one you can afford.

## 2.5 Solving the nonlinear problem

For a nonlinear static or quasi-static problem, define the residual

$$\mathbf{r}(\mathbf{u}) = \mathbf{f}^{\text{ext}} - \mathbf{f}^{\text{int}}(\mathbf{u}) = \mathbf{0}$$

**Newton–Raphson** linearises about the current iterate:

$$\mathbf{K}_T(\mathbf{u}^k)\,\Delta\mathbf{u} = \mathbf{r}(\mathbf{u}^k), \qquad \mathbf{u}^{k+1} = \mathbf{u}^k + \Delta\mathbf{u}$$

with $\mathbf{K}_T = \partial\mathbf{f}^{\text{int}}/\partial\mathbf{u}$ the tangent stiffness. Quadratic convergence near the solution, and no convergence at all far from it.

| Variant | What it changes | When |
|---|---|---|
| **Full Newton** | re-form and re-factor $\mathbf{K}_T$ every iteration | default; most robust per iteration |
| **Modified Newton** | reuse a factorization for several iterations | when factorization dominates cost; linear convergence |
| **Quasi-Newton (BFGS)** | low-rank updates to the inverse | between the two |
| **Line search** | scale $\Delta\mathbf{u}$ by $s$ minimising the residual | when full steps overshoot |

**Load incrementation is continuation.** Applying the load in steps is not about capturing the path for its own sake when the material is elastic — it is about starting each Newton solve from a good initial guess. For a path-dependent material it is both.

### Arc-length methods, and why load control fails

At a **limit point** — snap-through in a shallow arch, buckling collapse — the load–displacement curve turns back on itself. The structure passes through states where displacement increases while load *decreases*. A load-controlled Newton solve cannot find these: at the limit point $\mathbf{K}_T$ becomes singular, and past it there is no equilibrium at the prescribed load.

**Arc-length methods (Riks, Crisfield)** fix this by making the load factor an unknown and adding a constraint on the arc length travelled in load–displacement space. The path is then parameterised by something that increases monotonically even where load does not. This is the standard tool for post-buckling analysis, and it exists purely because the naive parameterisation of the problem is wrong.

### The consistent tangent

For plasticity, there are two candidate tangents: the **continuum tangent**, derived from the rate form of the constitutive law, and the **algorithmic (consistent) tangent**, obtained by differentiating the *discrete return-mapping algorithm* actually used. They differ.

**Only the consistent tangent preserves quadratic convergence.** Using the continuum tangent degrades Newton to linear convergence, sometimes severely. This is a well-known result (Simo and Taylor, 1985) and it is a good illustration of a general principle: **the Jacobian must be the derivative of the discrete operator you are actually applying, not of the continuous one you have in mind.** The same principle reappears in fluids as the argument for exact or matrix-free Jacobians.

## 2.6 The linear algebra

Structural mechanics is unusual in that **direct solvers remain the default.** The reasons are specific and they mostly do not apply elsewhere.

### Why direct

| Reason | |
|---|---|
| **SPD matrices** | sparse Cholesky needs no pivoting for stability, so the fill-in pattern can be computed symbolically in advance and the numerical factorization is a fixed, highly optimisable sequence |
| **Many right-hand sides** | multiple load cases, and Newton iterations under modified Newton, reuse one factorization at the cost of a forward-back substitution each |
| **Robustness** | a direct solve either succeeds or reports a singularity. An iterative solve can stall, and a commercial code cannot ask its user to tune a preconditioner |
| **Bad conditioning** | shells at $O(h^{-4})$, high stiffness contrasts, near-incompressibility. Direct methods are insensitive to conditioning in a way iterative methods are not |

The machinery: **fill-reducing ordering** (approximate minimum degree, nested dissection) to limit the fill created during elimination, then a **supernodal** or **multifrontal** factorization that groups columns with identical sparsity structure so the inner loops become dense BLAS-3 operations. Nearly all the performance of a modern sparse direct solver is in that last step.

### When iterative

Large 3D solid models, where fill-in makes direct factorization exceed memory. The tool is **preconditioned conjugate gradient** — CG applies exactly because $\mathbf{K}$ is SPD — with either incomplete Cholesky or, better, **algebraic multigrid**.

**AMG is unusually well suited to elasticity**, and for a structural reason: multigrid convergence depends on the coarse grids being able to represent the operator's near-null space, which for elasticity is the six **rigid-body modes**. Smoothed aggregation AMG supplied with those modes explicitly converges at a rate nearly independent of mesh size. Supplied without them, it does not. This is a case where knowing the physics changes the solver's asymptotic behaviour, not just its constant.

### Domain decomposition

**FETI** and **BDDC** were developed in structural mechanics and remain most natural there. The idea: partition the mesh, solve each subdomain independently with a direct solver, and enforce continuity across interfaces with Lagrange multipliers, iterating on the interface problem alone. This combines the robustness of direct methods with the memory scaling of iterative ones, and it parallelises naturally. **[unverified]** — attribution of FETI's origin to structural mechanics (Farhat and Roux, around 1991) is from memory and has not been checked against the paper.

### Where the saddle points come from

Three distinct situations in structural mechanics produce an indefinite system of the form

$$\begin{pmatrix} \mathbf{A} & \mathbf{B}^{\mathsf T} \\ \mathbf{B} & \mathbf{0} \end{pmatrix}\begin{pmatrix} \mathbf{u} \\ \boldsymbol\lambda \end{pmatrix} = \begin{pmatrix} \mathbf{f} \\ \mathbf{g} \end{pmatrix}$$

- **Near-incompressibility** ($\nu \to 1/2$), where pressure becomes an independent unknown
- **Lagrange-multiplier contact**, where $\boldsymbol\lambda$ is the contact pressure
- **Multi-point constraints** and rigid links

**All three are structurally identical to incompressible Navier–Stokes** (§4.3), and all three require the same **inf-sup (LBB) condition** on the choice of discrete spaces. CG no longer applies — the matrix is indefinite — so the solver moves to MINRES, or to a block-preconditioned Krylov method with a Schur-complement approximation. This is the clearest single point of contact between the two fields, and it arrives from completely different physics on each side.

---

# Part II — Fluid dynamics

# 3. The equations

## 3.1 The same balance law, a different closure

The starting point is the equation already written in §1.2:

$$\rho\frac{D\mathbf{u}}{Dt} = \nabla\cdot\boldsymbol\sigma + \rho\mathbf{f}$$

**This is Cauchy's momentum equation, and it is valid for any material whatsoever** — solid, fluid, or anything else. It is Newton's second law for a continuum and nothing has yet been said about how the material responds.

Splitting the stress into an isotropic pressure and a deviatoric part, $\boldsymbol\sigma = -p\,\mathbf{I} + \boldsymbol\tau$, and supplying the **Newtonian closure**

$$\boldsymbol\tau = \mu\left(\nabla\mathbf{u} + (\nabla\mathbf{u})^{\mathsf T}\right) + \left(\zeta - \frac{2}{3}\mu\right)(\nabla\cdot\mathbf{u})\,\mathbf{I}$$

turns Cauchy's equation into the Navier–Stokes equation. **The closure, not the conservation law, carries the name.**

"Newtonian" means the stress is assumed **linear in the velocity gradient**. Imposing isotropy and symmetry then forces this exact expression: a linear isotropic relation between two symmetric rank-2 tensors admits exactly two independent scalars, and $\mu$ (shear viscosity) and $\zeta$ (bulk viscosity) are they. The $-\frac{2}{3}\mu$ exists only to remove the trace from the first term so the two pieces do not overlap.

**The parallel with §1.4 is exact.** Hooke's law and the Newtonian closure occupy the same slot in the same universal balance law. The difference is what stress responds to:

| | Stress responds to | Frame |
|---|---|---|
| **Elastic solid** | **strain** — accumulated deformation | Lagrangian |
| **Newtonian fluid** | **strain rate** — rate of deformation | Eulerian |

A solid remembers where it started; a fluid does not. Everything else in Part III follows from those two rows.

## 3.2 The compressible system

Five scalar equations, closed by two constitutive relations.

**Continuity:**

$$\frac{\partial\rho}{\partial t} + \nabla\cdot(\rho\mathbf{u}) = 0$$

**Momentum:**

$$\rho\frac{D\mathbf{u}}{Dt} = -\nabla p + \nabla\cdot\boldsymbol\tau + \rho\mathbf{f}$$

**Energy**, in internal-energy form:

$$\rho\frac{D\varepsilon}{Dt} = -p\,\nabla\cdot\mathbf{u} + \boldsymbol\tau : \nabla\mathbf{u} + \nabla\cdot(k\nabla T) + Q$$

The double contraction $\boldsymbol\tau : \nabla\mathbf{u} = \sum_{i,j}\tau_{ij}\partial u_i/\partial x_j$ is the **viscous dissipation rate**, always non-negative for a Newtonian fluid — the second law appearing in the equations. It is the term that heats a fluid when you stir it.

**Closure — two relations, not one**, and sources routinely supply both silently under the single heading "equation of state":

$$p = p(\rho, T) \quad\text{(thermal EOS)}, \qquad \varepsilon = \varepsilon(T) \quad\text{(caloric relation)}$$

For an ideal gas these are $p = \rho R T$ and $\varepsilon = c_v T$, which combine into $p = (\gamma - 1)\rho\varepsilon$ — which is why they are often not distinguished. They must be distinguished for a **tabulated EOS**, where ionization and excitation change the two relations independently.

**The two time derivatives** are different operators, not notational variants:

$$\frac{\partial}{\partial t} \quad\text{— at a fixed point in space}, \qquad \frac{D}{Dt} = \frac{\partial}{\partial t} + (\mathbf{u}\cdot\nabla) \quad\text{— following a fluid parcel}$$

They differ by the advective term, which accounts for a parcel moving into a region where the quantity has a different value. **This term is the origin of nearly every difficulty in computational fluid dynamics**, and it exists only because the frame is Eulerian.

### Conservative versus primitive form

$$\frac{\partial\mathbf{q}}{\partial t} + \nabla\cdot\mathbf{F}(\mathbf{q}) = \mathbf{S}, \qquad \mathbf{q} = (\rho,\ \rho\mathbf{u},\ E)^{\mathsf T}$$

The **conservative** form evolves the conserved quantities themselves. The **primitive** form evolves $(\rho, \mathbf{u}, p)$ and is easier to read. They are algebraically equivalent for smooth solutions and **not equivalent across a shock**, where the solution is not differentiable and only the integral form has meaning. Finite-volume codes discretize the conservative form for exactly this reason (§4.1).

## 3.3 The incompressible system

Four scalar equations, and a structurally different problem.

**Continuity becomes a constraint:**

$$\nabla\cdot\mathbf{u} = 0$$

**The time derivative is gone.** This is no longer an evolution equation but a condition the velocity field must satisfy at every instant.

**Momentum:**

$$\frac{D\mathbf{u}}{Dt} = -\frac{1}{\rho}\nabla p + \nu\,\nabla^2\mathbf{u} + \mathbf{f}$$

The viscous term collapses to a plain Laplacian only because $\mu$ is constant *and* the flow is divergence-free: incompressibility kills both the bulk-viscosity term and the transpose term, since $[\nabla\cdot(\nabla\mathbf{u})^{\mathsf T}]_i = \partial_i(\nabla\cdot\mathbf{u}) = 0$.

**Energy decouples.** With $\rho$ fixed and $\mu$ independent of $T$, nothing in the momentum equation depends on temperature, so $T$ can be advected afterwards through a known velocity field as a passive scalar. The coupling is one-way, $\mathbf{u}\to T$. It returns two-way through buoyancy, which is what the **Boussinesq approximation** exists to represent: keep $\nabla\cdot\mathbf{u}=0$ everywhere except in a body-force term $-\rho_0\beta(T-T_0)\mathbf{g}$.

### Where the pressure comes from — the central fact of incompressible CFD

**There is no evolution equation for $p$, and no equation of state for it either.** Pressure here is not a thermodynamic variable at all: it is a **Lagrange multiplier** enforcing $\nabla\cdot\mathbf{u} = 0$. Only $\nabla p$ appears, so its absolute value is undetermined.

Taking the divergence of the momentum equation and imposing the constraint gives

$$\nabla^2 p = -\rho\,\nabla\cdot\left[(\mathbf{u}\cdot\nabla)\mathbf{u}\right]$$

an **elliptic Poisson equation, to be solved at every timestep.** This single fact dominates the cost, the parallel scaling and the solver design of every incompressible flow code in existence.

### How a variable is obtained — the axis that matters

| How the variable is obtained | What you do | Reach |
|---|---|---|
| **evolved** — has $\partial/\partial t$ | march it | local |
| **constrained, algebraic** — pointwise or a finite stencil | evaluate it | local |
| **constrained, elliptic** — depends on the whole domain at once | **solve** for it | nonlocal |

| System | Evolved | Constrained, algebraic | Constrained, elliptic |
|---|---|---|---|
| **Compressible** | $\rho$, $\mathbf{u}$, $\varepsilon$ | $p$ — from the EOS | — |
| **Incompressible** | $\mathbf{u}$ | — | **$p$** |
| **Static elasticity** (§1.5) | — | — | **$\mathbf{u}$** |

**The distinction is not compressible against incompressible. It is whether the closure is local.** An equation of state is local: given $\rho$ and $T$ *here*, $p$ follows *here*. The incompressible constraint is not — it must hold everywhere simultaneously, and enforcing it couples the entire domain.

Systems mixing evolved and algebraically-constrained variables are **differential-algebraic equations**. Incompressible Navier–Stokes is a DAE of **index 2**: the index counts how many times the constraint must be differentiated before it becomes an evolution equation, which is a formal measure of how badly it resists being marched.

### The irony of incompressibility

Incompressible flow is the $M \to 0$ limit — the limit in which the **sound speed becomes infinite**. Compressible flow is hyperbolic with signals travelling at $|u| \pm c_s$, and explicit stepping needs $\Delta t \lesssim \Delta x/(|u| + c_s)$, which is entirely workable when $c_s \sim u$. Send $c_s\to\infty$ and that timestep goes to zero.

**So incompressibility does not avoid the fast wave. It makes the wave infinitely fast, removes it from the equations by fiat, and replaces it with a constraint.** The elliptic solve is the price paid for a timestep set by $u$ rather than by $c_s$. It is a trade, not a defect.

The converse holds too, and it is the reason low-Mach compressible solvers are hard: **at $M \ll 1$, an explicit compressible solver is also bad**, because every step resolves acoustic waves nobody cares about. That is stiffness, not ellipticity, and the fixes are different (§4.4).

## 3.4 The reduced-model ladder

Almost all practical CFD is done on a reduction, not on the full compressible Navier–Stokes system. Each reduction drops a term and changes the character.

| Model | What is dropped | Character | Typical method |
|---|---|---|---|
| **Compressible Navier–Stokes** | — | hyperbolic + parabolic | FVM, explicit or dual-time implicit |
| **Euler** | viscosity, conduction | **purely hyperbolic** | Godunov-type FVM, Riemann solvers |
| **Incompressible Navier–Stokes** | compressibility | parabolic + **elliptic constraint** | projection / SIMPLE / PISO |
| **Stokes** | the convective term entirely | **elliptic saddle point** | direct or block-preconditioned Krylov |
| **Potential flow** | viscosity and vorticity | $\nabla^2\phi = 0$, **elliptic** | boundary element, panel methods |
| **Steady small-disturbance potential** | as above, linearised | elliptic if $M<1$, **hyperbolic if $M>1$** | type-dependent switching |
| **Boundary layer (Prandtl)** | streamwise diffusion | **parabolic** in the streamwise direction | space-marching — no global solve at all |
| **Shallow water** | vertical structure | hyperbolic | FVM, same machinery as Euler |
| **Lubrication (Reynolds eq.)** | inertia, in a thin gap | elliptic | direct solve |

**Two entries deserve comment.**

**Stokes flow is the fluid problem that looks most like structural mechanics.** Dropping the convective term removes the only nonlinearity, leaving a linear, self-adjoint, elliptic saddle-point problem. Its discrete form is essentially the same object as near-incompressible elasticity (§2.6), and the same inf-sup condition governs both.

**The boundary-layer equations are parabolic, and that is their whole point.** Dropping streamwise diffusion removes upstream influence, so the equations can be **marched in space** from the leading edge downstream, solving a small system at each station rather than a global one. This is what made boundary-layer theory computationally tractable decades before anyone could solve Navier–Stokes, and it fails exactly where upstream influence returns — at separation.

## 3.5 Turbulence — closures again

Turbulence is not a separate equation. The Navier–Stokes equations already contain it; the problem is that resolving it is unaffordable.

| Approach | What is resolved | What is modelled | Cost |
|---|---|---|---|
| **DNS** | everything down to the Kolmogorov scale | nothing | grid points scale roughly as $\mathrm{Re}^{9/4}$; total cost as $\mathrm{Re}^{3}$ **[unverified]** |
| **LES** | the energy-containing eddies | subgrid stresses (Smagorinsky, WALE, dynamic) | far cheaper, but wall-resolved LES still scales steeply with $\mathrm{Re}$ |
| **RANS** | the mean flow only | the entire Reynolds stress tensor | independent of $\mathrm{Re}$ in cost; carries the modelling burden |

**RANS is the industrial default.** Averaging the momentum equation produces an unclosed term, the **Reynolds stress** $-\rho\overline{u_i'u_j'}$, and the model supplies it. The most common route is an **eddy-viscosity** hypothesis — assume the Reynolds stress is proportional to the mean strain rate, exactly as the Newtonian closure assumes for the molecular stress — which reduces the problem to specifying a scalar $\nu_t$.

| Model | Transport equations added |
|---|---|
| **Spalart–Allmaras** | one, for a variable closely related to $\nu_t$ |
| **$k$–$\varepsilon$** | two: turbulent kinetic energy and its dissipation rate |
| **$k$–$\omega$ SST** | two, blending $k$–$\omega$ near walls with $k$–$\varepsilon$ outside |
| **Reynolds stress models** | seven — each stress component plus a scale |

**A turbulence model is the same kind of object as a constitutive law.** It is not a conservation statement; it is an assertion about how the medium responds, supplied to make the equation count work. The parallel with §1.4 is exact, and so is the consequence: this is where models diverge while the conservation laws above them stay fixed.

**What it costs the solver**, concretely:

- Extra advection–diffusion–**reaction** transport equations, and the reaction (source) terms are often **stiff** and can be much stiffer than the flow itself.
- **Positivity requirements** — $k$, $\varepsilon$, $\omega$ and $\nu_t$ must remain non-negative, and a scheme that produces a negative value produces a crash, not an inaccuracy. Limiters and clipping are ubiquitous and are a real source of unphysical behaviour.
- **Near-wall resolution.** Resolving to $y^+ \sim 1$ puts the first cell height in the micrometre range for typical industrial Reynolds numbers, with cell aspect ratios of $10^3$–$10^4$. That drives both the linear-solver conditioning and the case for implicit time integration (§4.4).

---

# 4. Methodology — how a fluid problem becomes algebra

## 4.1 Why the finite volume method

The dominance of FVM in fluid dynamics is as specific as the dominance of FEM in structures, and for a different reason.

**1. Conservation is the point.** FVM discretizes the *integral* form over each cell:

$$\frac{\mathrm{d}}{\mathrm{d}t}\int_{V_i}\mathbf{q}\,\mathrm{d}V + \oint_{\partial V_i}\mathbf{F}\cdot\mathbf{n}\,\mathrm{d}S = \int_{V_i}\mathbf{S}\,\mathrm{d}V$$

Each interior face is shared by exactly two cells, and the flux leaving one is the flux entering the other. **The fluxes telescope, so mass, momentum and energy are conserved to machine precision on any mesh, at any resolution.** That is a discrete statement of the physics, not an approximation that improves with refinement.

**2. Shocks come out right for free.** The Lax–Wendroff theorem states that if a consistent **conservative** scheme converges, it converges to a weak solution of the conservation law — with the correct shock speed. A non-conservative scheme can converge smoothly to a solution with a shock in the wrong place, which is the worst possible failure mode: plausible, stable, and wrong.

**3. Arbitrary cell shapes and finite-volume geometry.** Only face areas and normals are needed, so unstructured polyhedral meshes are natural.

### Why not Galerkin FEM, given §2.1

The optimality result that justifies FEM in structures **does not transfer**. Céa's lemma requires the operator to be coercive; for a self-adjoint elliptic operator the Galerkin solution is the energy-norm best approximation. The advection operator is **not self-adjoint**, and standard Galerkin applied to advection-dominated problems is equivalent to central differencing — oscillatory and, at high cell Péclet number, useless.

FEM is used in fluids, but only with **stabilisation** added:

| Method | What it adds |
|---|---|
| **SUPG** (streamline-upwind Petrov–Galerkin) | test functions weighted along the streamline — upwinding, in a variational setting |
| **PSPG** | a pressure-stabilising term that circumvents the inf-sup condition, permitting equal-order velocity and pressure |
| **GLS**, **variational multiscale** | generalisations of the same idea |

The other discretizations in use: **finite difference** on structured grids (simple, high order, geometrically limited); **spectral and spectral element** methods, which dominate DNS because they resolve a given accuracy with far fewer points; and **lattice Boltzmann**, which solves a kinetic equation whose macroscopic limit is Navier–Stokes rather than solving Navier–Stokes at all.

## 4.2 Advection is the hard part

Nearly every difficulty specific to CFD traces to the convective term, and the structure of that difficulty is worth stating precisely.

**The physics is upwind.** For $\partial_t\phi + c\,\partial_x\phi = 0$ with $c > 0$, information arrives from the left. A stencil that looks symmetrically both ways is asking about a region the exact solution never consults. Central differencing does this, and it is unconditionally unstable on pure advection — a useful reminder that the CFL condition is necessary and not sufficient.

**But upwinding is only first-order accurate**, and first-order upwinding introduces numerical diffusion so large that on typical meshes it can exceed the physical viscosity. In a turbulence calculation that means the scheme is doing the modelling.

**Godunov's theorem is why this cannot simply be fixed:** any *linear* scheme that preserves monotonicity is at most first-order accurate. There is no linear escape. Every high-resolution scheme is therefore **nonlinear by construction**, even when applied to a linear equation:

| Family | Idea |
|---|---|
| **Flux limiters / TVD** | blend a high-order and a low-order flux, with the blend depending on the local solution smoothness |
| **MUSCL** | reconstruct a linear profile in each cell, limit its slope, solve a Riemann problem at each face |
| **ENO / WENO** | choose or weight stencils to avoid interpolating across a discontinuity; high order in smooth regions |
| **Discontinuous Galerkin** | high order within cells, with a Riemann-solver coupling between them |

**Riemann solvers** are the machinery at each face: given two constant states, what is the exact evolution of the discontinuity between them? Exact solvers are iterative and expensive; **Roe**, **HLL**, **HLLC** and **Rusanov** are approximations trading accuracy for cost. This entire apparatus exists because the equations are hyperbolic and admit discontinuous solutions, and it has no counterpart at all in structural mechanics.

## 4.3 Pressure–velocity coupling

For incompressible flow this is the central algorithmic problem, and it exists because §3.3 leaves no equation to advance $p$ with.

### The families

| Approach | How it works |
|---|---|
| **Projection / fractional step** (Chorin) | advance momentum ignoring the constraint to get an intermediate $\mathbf{u}^*$; solve a Poisson equation for pressure; subtract the pressure gradient to project $\mathbf{u}^*$ onto the divergence-free subspace |
| **SIMPLE** and relatives (SIMPLEC, SIMPLER) | a segregated iteration: guess $p$, solve momentum, derive a pressure *correction* from the continuity residual, correct both, repeat. Under-relaxation is essential |
| **PISO** | SIMPLE with additional corrector steps, aimed at transient rather than steady problems |
| **Artificial compressibility** | add $\partial p/\partial\tau$ with a pseudo-time, making the system hyperbolic and marchable; the constraint is satisfied only at pseudo-steady state |
| **Coupled / monolithic** | solve for $\mathbf{u}$ and $p$ simultaneously as one saddle-point system, with a block preconditioner |

The projection viewpoint is the most illuminating: **the Helmholtz decomposition splits any vector field into a divergence-free part and a gradient, and the pressure solve is exactly the extraction of the gradient part.** Pressure is whatever gradient field must be removed to make the velocity divergence-free.

### Checkerboarding and inf-sup

On a collocated grid with central differences, the pressure at a node couples to $p_{i-2}$ and $p_{i+2}$ but not to $p_{i\pm1}$. The odd and even nodes decouple, and an oscillatory pressure field of alternating high and low values is invisible to the discrete gradient operator — a spurious mode. This is **checkerboarding**, and it is the finite-difference face of the same object as the **inf-sup (LBB) condition** in finite elements: velocity and pressure spaces must be compatible, and equal-order interpolation generally is not.

Three standard fixes:

- **Staggered grids** (MAC): store pressure at cell centres and velocities at faces, so each pressure difference drives an adjacent velocity. Structurally correct, awkward on unstructured meshes.
- **Rhie–Chow interpolation**: on a collocated grid, reconstruct face velocities in a way that reintroduces the missing coupling. This is what most modern unstructured codes do.
- **Inf-sup-stable element pairs** (Taylor–Hood: quadratic velocity, linear pressure), or **pressure stabilisation** (PSPG) to sidestep the requirement.

**All of this reappears verbatim in near-incompressible elasticity and Lagrange-multiplier contact** (§2.6). It is the same mathematics arriving from unrelated physics.

## 4.4 Time integration in fluids

### Why explicit is viable here and not in structures

The convective CFL condition,

$$\Delta t \le \frac{\Delta x}{|u| + c_s} \quad\text{(compressible)}, \qquad \Delta t \le \frac{\Delta x}{|u|} \quad\text{(incompressible)}$$

**is not a large penalty, because in fluids the timescale of interest is usually the advective timescale itself.** An eddy turnover, a vortex shedding cycle, a shock transit — these happen at exactly the rate at which information crosses the mesh. The fastest resolvable mode and the physics of interest are close together, so an explicit method wastes very little.

**That is the whole asymmetry with structural mechanics**, where the fastest resolvable mode is a dilatational wave crossing a millimetre and the physics of interest is a building swaying over seconds. Same criterion, opposite answer (§5.2).

### Where implicit returns in fluids

It returns whenever that alignment of timescales breaks, and there are four common ways for it to break:

| Situation | Why the timescales separate | Standard treatment |
|---|---|---|
| **Low Mach number** | acoustic waves at $c_s \gg u$ are resolved and uninteresting | low-Mach preconditioning, all-speed schemes, or fully implicit |
| **Steady RANS** | the transient is not wanted at all, only its endpoint | pseudo-time marching, **local time stepping** (a different $\Delta t$ in every cell — accuracy in time deliberately abandoned), implicit residual smoothing, multigrid on the steady residual |
| **Near-wall viscous cells** | at $y^+\sim1$ the wall-normal cell is micrometres, and the viscous limit $\Delta t \le \Delta y^2/2\nu$ is brutal | implicit in the wall-normal direction, or dual-time stepping |
| **Stiff chemistry (combustion)** | reaction timescales orders of magnitude below flow timescales | operator splitting with an implicit ODE integrator per cell |

**Note what happens in row 2.** In a steady calculation, time accuracy has no value, so every trick that trades it for convergence rate becomes available. Local time stepping — where neighbouring cells are at different times and the field is not a physical state at any instant — is the clearest example, and it would be nonsense in a transient calculation.

### The labels are per-term, not per-code

**Even a fully explicit incompressible solver contains an elliptic solve at every step.** The advective and viscous terms are marched; the pressure is solved. Calling such a code "explicit" describes the momentum treatment only.

This is the general situation, and the honest framing of the whole explicit/implicit question: **each term gets the treatment its character demands.** The formalisation is **IMEX** — implicit for the stiff terms, explicit for the rest, paying for a solve only where stiffness actually lives — and **operator splitting**, which sequences the terms and applies a suitable integrator to each. Fluid dynamics has been doing this in practice, under other names, for longer than the terminology has existed.

## 4.5 Nonlinearity in fluids

**The nonlinearity is in the operator itself**, not in a constitutive law: $(\mathbf{u}\cdot\nabla)\mathbf{u}$ is quadratic in the unknown. Three consequences that distinguish it from the structural case:

**1. It is mild in form and severe in effect.** A quadratic term is the gentlest possible nonlinearity, and the reason turbulence is difficult is not that the nonlinearity is complicated but that it couples every scale to every other.

**2. Picard is a genuine competitor to Newton.** Lagging one factor,

$$(\mathbf{u}^k\cdot\nabla)\mathbf{u}^{k+1}$$

gives a *linear* problem in $\mathbf{u}^{k+1}$ — the **Oseen** problem. Convergence is linear rather than quadratic, but the basin of attraction is much larger. Newton adds the missing term $(\mathbf{u}^{k+1}\cdot\nabla)\mathbf{u}^{k}$ and converges quadratically, near the solution and only there. **The standard practice is Picard first, then switch to Newton when the residual is small enough.** Nothing analogous is standard in structural mechanics, where Newton with a consistent tangent is the default and the difficulty lies in the constitutive update.

**3. Discontinuous solutions are a nonlinearity of a different kind.** A shock is not a nonlinearity of the operator — the operator is perfectly smooth — but of the *solution*, which ceases to be differentiable. That is why §4.2's machinery is nonlinear even for the linear advection equation, and it has no structural counterpart except in the softening and fracture regimes.

Additional sources of nonlinearity: a non-ideal or tabulated EOS, temperature-dependent transport properties, free surfaces (where the domain is unknown), and turbulence models with their own coupled nonlinear source terms.

## 4.6 The linear algebra

Where structures default to direct solvers, **fluids default to iterative ones**, and the reasons invert one by one.

| | Structures | Fluids |
|---|---|---|
| Matrix symmetry | SPD | **nonsymmetric** (advection) or **indefinite** (saddle point) |
| Applicable Krylov method | CG | GMRES, BiCGSTAB, or MINRES |
| Right-hand sides per factorization | many (load cases, modified Newton) | **one** — the operator changes every step |
| Required solve tolerance | tight | **loose** — it sits inside an outer nonlinear iteration |
| Typical size | large | **larger** |

The fourth row is underrated. Inside an outer Newton or SIMPLE iteration, an over-solved linear system is wasted work: the right-hand side is about to change anyway. **Inexact Newton** methods exploit this deliberately, with Eisenstat–Walker forcing terms tightening the linear tolerance as the nonlinear residual falls. There is nothing to exploit for a static structural analysis with one load case, which is part of why direct methods survive there.

### The pressure Poisson equation is the dominant cost

In incompressible flow, the elliptic pressure solve typically dominates both runtime and parallel communication, because it is the only genuinely global operation in the algorithm. It is also, conveniently, **SPD** — the one piece of an incompressible flow solver that looks like a structural problem — which makes it the natural target for **conjugate gradient with algebraic or geometric multigrid**. Getting a scalable pressure solver is most of getting a scalable incompressible code.

### Preconditioning the coupled system

For monolithic solvers the system is the saddle point of §2.6 with a nonsymmetric $\mathbf{A}$. The preconditioners are built around approximating the **Schur complement** $\mathbf{S} = -\mathbf{B}\mathbf{A}^{-1}\mathbf{B}^{\mathsf T}$, which is dense and never formed:

- **PCD** (pressure convection–diffusion) and **LSC** (least-squares commutator) — spectrally-motivated approximations
- **SIMPLE used as a preconditioner** rather than as an outer iteration, which is a common and effective repurposing
- **Physics-based block preconditioners** generally, where the blocks correspond to the terms of the PDE

### Matrix-free and JFNK

For large problems the Jacobian is often never assembled. **Jacobian-free Newton–Krylov** needs only the action $\mathbf{J}\mathbf{v}$, which can be approximated by a directional finite difference of the residual,

$$\mathbf{J}\mathbf{v} \approx \frac{\mathbf{r}(\mathbf{u} + \epsilon\mathbf{v}) - \mathbf{r}(\mathbf{u})}{\epsilon}$$

or computed exactly by automatic differentiation. **The preconditioner then becomes the entire difficulty**, since a Krylov method without one converges too slowly to be useful, and the usual construction is a cheap approximate Jacobian assembled from a simplified physics model.

---

# Part III — The two side by side

# 5. Comparison

## 5.1 The map, in one table

| | **Structural mechanics** | **Fluid dynamics** |
|---|---|---|
| **Balance law** | Cauchy momentum | **the same Cauchy momentum** |
| **Closure** | stress responds to **strain** (Hooke, plasticity, hyperelasticity) | stress responds to **strain rate** (Newtonian) |
| **Frame** | **Lagrangian** — mesh moves with material | **Eulerian** — material moves through mesh |
| **Convective term** | absent | $(\mathbf{u}\cdot\nabla)\mathbf{u}$ — the central difficulty |
| **Canonical static problem** | $\mathbf{K}\mathbf{u} = \mathbf{f}$, **elliptic** | Stokes flow, **elliptic saddle point** |
| **Canonical transient problem** | elastodynamics, **hyperbolic** | compressible flow, **hyperbolic + parabolic** |
| **Where an elliptic solve hides** | the problem *is* elliptic | the **pressure** in incompressible flow |
| **Discretization** | **FEM** — variational, optimal for elliptic operators | **FVM** — conservative, correct shock speeds |
| **Time integration** | mostly **implicit** (generalized-$\alpha$); explicit for impact | mostly **explicit**; implicit for low Mach, steady, near-wall, chemistry |
| **Nonlinearity from** | constitutive law, finite strain, contact | the operator itself; plus EOS and turbulence closure |
| **Nonlinear solver** | **Newton** with consistent tangent; arc-length at limit points | **Picard then Newton**; SIMPLE-type outer iterations |
| **Matrix** | sparse **SPD**, banded | sparse **nonsymmetric** or **indefinite** |
| **Linear solver** | **sparse direct** (Cholesky, multifrontal); PCG+AMG when large | **iterative** (GMRES, BiCGSTAB); multigrid on the pressure Poisson |
| **Right-hand sides per operator** | **many** — factorizations are reused | **one** — the operator changes every step |
| **Characteristic difficulty** | contact, locking, softening, shell conditioning | shocks, turbulence, pressure coupling, wall resolution |
| **Eigenvalue problems** | **central** — modal analysis, buckling | peripheral — stability analysis only |

## 5.2 The four asymmetries, and what causes each

Everything in that table reduces to four differences. It is worth being precise about which is a cause and which is a consequence.

### 1. The frame — Lagrangian against Eulerian

**This is the root cause, and nearly everything else follows.**

A solid remembers its undeformed configuration, so attaching the mesh to the material is natural and $\mathbf{u}$ is a displacement from a definite reference. A fluid does not remember, and a Lagrangian mesh in a fluid tangles within a few turnover times, so the mesh stays fixed and material flows through it.

Two consequences, and they run in opposite directions:

- **Fluids get the convective term and solids do not.** $(\mathbf{u}\cdot\nabla)$ exists only because the frame is fixed. In a Lagrangian frame the material derivative *is* the partial derivative.
- **Solids get geometric nonlinearity and fluids mostly do not.** When the mesh follows the material, the domain itself is an unknown.

So each field pays for its frame, in different coin. Fluids pay with a nonlinear, non-self-adjoint operator that needs upwinding, limiters and Riemann solvers. Solids pay with a moving reference configuration, stress measures that must be transformed between configurations, and a tangent that changes because the geometry changes.

### 2. Where the nonlinearity lives

| | Nonlinearity is in | Consequence |
|---|---|---|
| **Structures** | the **constitutive law** and the kinematics | it is local, and evaluated per Gauss point; history variables; return mapping; a consistent tangent recovers quadratic convergence |
| **Fluids** | the **operator** | it is global and quadratic; Picard is competitive; the difficulty is scale coupling, not material response |

This is why the two fields' nonlinear solvers look different despite both being Newton. In structures, most of the engineering effort goes into the material model and its consistent linearisation. In fluids, most of it goes into the discretization of a single quadratic term.

### 3. Timescale separation — the answer to "why do structures go implicit"

This is the question the document was written around, and it has a clean answer.

| | Fastest resolvable mode | Physics of interest | Ratio |
|---|---|---|---|
| **Structure, seismic** | $c_p/\Delta x$ — wave across a 1 mm element, ~170 ns | building sway, ~1 s | $\sim10^{7}$ |
| **Structure, crash** | the same, ~170 ns | 100 ms event | $\sim10^{6}$, but only $6\times10^5$ steps — affordable |
| **Structure, ballistic** | the same | 100 µs event | $\sim10^{3}$ — explicit is the natural choice |
| **Fluid, external aero** | $\Delta x/(|u|+c_s)$ | convective transit of the body | $\sim1$–$10$ |
| **Fluid, low Mach** | acoustic transit | convective transit | $1/M$ — large, and this is why low-Mach is hard |

**The criterion is the same everywhere: is the fastest resolvable mode close to the physics you care about?** If yes, explicit wastes nothing. If not, the system is stiff and implicit becomes attractive.

Structures usually answer no, because a stiff material supports very fast waves while structures are loaded slowly. Fluids usually answer yes, because the mesh is sized to the flow and the flow *is* the physics. **The apparent field-level difference is not a field-level difference at all — it is the same test giving different answers on different problems**, and the counterexamples in both directions (crash on one side, low-Mach on the other) confirm it.

Two further reasons stand behind the timescale one and are easy to miss:

- **Many structural problems have no time.** Statics, modal, buckling. Implicit is not preferred there; it is the only thing that exists, because explicit integration requires a mass matrix to divide by.
- **Explicit is chosen for robustness, not only for speed.** Contact and material failure make Newton fail, and a failed Newton iteration stops the analysis. Explicit never asks it to converge.

### 4. Matrix character — and hence direct against iterative

| | | |
|---|---|---|
| **Structures** | SPD from a symmetric positive-definite variational problem | Cholesky needs no pivoting; fill-in is predictable; CG applies; AMG works well given the rigid-body near-null space |
| **Fluids** | nonsymmetric from advection, indefinite from the pressure constraint | needs GMRES or MINRES; preconditioning is the whole problem; multigrid is applied to the pressure Poisson, which is the one SPD piece |

The reuse asymmetry compounds this. A structural factorization is amortised over load cases and Newton iterations; a fluid operator changes every step and never is.

## 5.3 Where they meet

Three genuine points of contact, and they are worth naming because each is a place where a technique transfers.

**Saddle points.** Near-incompressible elasticity, Lagrange-multiplier contact, Stokes flow and incompressible Navier–Stokes all produce the same block system and all require the same inf-sup condition. A block preconditioner developed for one applies to the others. **This is the most transferable single piece of machinery between the two fields.**

**Fluid–structure interaction.** The coupling itself, in two flavours:

- **Monolithic** — assemble one system for fluid and structure together. Robust, expensive, and requires a solver that copes with the union of both matrix characters.
- **Partitioned** — solve each field with its own specialised code and iterate on the interface. Reuses existing solvers, which is decisive in practice, but suffers the **added-mass instability**: when the fluid and structural densities are comparable, a weakly coupled (single-pass) scheme is unconditionally unstable regardless of timestep. Strong coupling with sub-iterations, or a monolithic formulation, is then required. **[unverified]** — the unconditional character of the added-mass instability for weakly coupled schemes is recalled from the literature and has not been checked against a source here.

The mesh treatment is **ALE** (arbitrary Lagrangian–Eulerian), which lets the mesh move with the material at the interface and independently of it in the interior — an explicit interpolation between the two frames of §5.2.

**The Lagrangian fluid methods.** SPH, material point methods and particle methods are fluid solvers that adopt the solid frame, and they inherit both sides of the trade: no convective term, but a moving reference and difficulty enforcing conservation and boundary conditions. They are also where the two fields' impact codes actually merge, in high-velocity impact and fragmentation.

## 5.4 Misconceptions, stated flatly

Collected because each one costs time when believed.

- **Implicit does not mean nonlinear**, and nonlinear does not mean implicit.
- **Coupling does not force implicit.** An elliptic constraint or stiffness does. A tightly coupled system whose closure is algebraic can be fully explicit.
- **Elliptic problems are not implicit by preference.** There is no time level to be explicit about.
- **Unconditionally stable is not accurate.** A large implicit step gives a stable, smooth, wrong answer, and it does so quietly.
- **"Stiffness" has two unrelated meanings in structural mechanics.** The stiffness matrix $\mathbf{K}$ is a physical rigidity; a stiff ODE system has widely separated timescales. They are related — a rigid structure has high natural frequencies, $\omega \sim \sqrt{K/M}$ — but they are not synonyms, and a system can be numerically stiff without being physically rigid.
- **Stiffness is not a property of the equations alone.** It depends on the equations, the grid, and the question. Refining the mesh creates it; changing the question destroys it.
- **A converged run is not a correct run.** In both fields, tightening the numerics can convert a loud failure into a smooth wrong answer. A diagnostic with a known exact answer — a conserved quantity, a component that must remain identically zero, an analytic solution — is the only thing that distinguishes them.

---

# 6. Quick reference

**Given an equation, what do I do with it?**

| If the system is | Then | And the solver is |
|---|---|---|
| **elliptic** | there is nothing to march | a global linear or nonlinear solve. SPD → CG/AMG or Cholesky; indefinite → MINRES or a block preconditioner |
| **parabolic** | explicit costs $\Delta t \propto \Delta x^2$ | usually implicit; a tridiagonal-like solve per step, cheap and unconditionally stable |
| **hyperbolic** | explicit costs $\Delta t \propto \Delta x$ | usually explicit, with upwinding and a limiter. Implicit only if the wave is uninteresting |
| **hyperbolic + an elliptic constraint** | the constraint must hold at every instant | split: march the hyperbolic part, solve the constraint. Projection, SIMPLE, PISO |
| **stiff but not elliptic** | the fast mode is stable and uninteresting | IMEX, operator splitting, or fully implicit with a Newton–Krylov solve per step |
| **an eigenproblem** | no marching, no boundary-value solve | Lanczos, shift-and-invert, subspace iteration |
| **non-smooth (contact, complementarity)** | Newton has no derivative | active set, penalty, augmented Lagrangian — or go explicit and avoid the question |

---

# Part IV — Optimization

> **This part is a different kind of section, and the difference matters.** Parts I and II are **fields**: a physical domain, its governing equations, and the numerical methods those equations force. Optimization is not a field. It is a **layer** that wraps a field — the simulation of Parts I and II becomes a *constraint*, and everything in §0.1's chain becomes the inner problem, solved repeatedly inside an outer loop.
>
> That is why it belongs in this document rather than beside it. Adding an objective does not sit on top of the chain; it changes every step of it, and §11 works through how. It is also why it is not "Part III" alongside two fields: filing a layer as a field would be a category error, and the next section added should be checked against this distinction before it is placed.

---

# 7. The shape of the problem

## 7.1 The canonical statement

$$\min_{x\in\mathbb{R}^n} f(x) \quad\text{subject to}\quad c_E(x) = 0, \quad c_I(x) \le 0$$

$f$ is the **objective**, $x$ the **design** or **decision variables**, and the two constraint sets the feasible region. Everything in optimization is a variation on this, and the variations that matter are listed in §7.3.

## 7.2 The PDE-constrained form, and the two ways to arrange it

In the setting of this document the constraint set contains a simulation. Separate the variables accordingly:

| | |
|---|---|
| $x$ | **design variables** — thicknesses, boundary shape, material distribution, an inflow condition |
| $u$ | **state variables** — the displacement, velocity or pressure field. What Parts I and II solve for |
| $R(u, x) = 0$ | the **state equation** — the discretized PDE. This *is* $\mathbf{K}\mathbf{u} - \mathbf{f}$, or the fluid residual |
| $J(u, x)$ | the objective — compliance, drag, mass, a misfit against data |

Two formulations follow, and the choice between them is the first architectural decision in any optimization code.

| | **Reduced space** (nested, NAND, "black box") | **Full space** (all-at-once, SAND, one-shot) |
|---|---|---|
| Unknowns | $x$ only; $u = u(x)$ implicitly | $x$ **and** $u$, treated as independent |
| $R = 0$ | satisfied exactly at every iterate | satisfied only at convergence |
| Each outer iteration costs | **a full simulation**, plus a gradient | one Newton-type step on a larger system |
| Optimizer sees | a smooth unconstrained (or lightly constrained) problem in $n$ variables | a large equality-constrained problem |
| Intermediate iterates | are physically meaningful | are **not** — the "flow field" does not satisfy the flow equations |
| Reuses an existing solver | **yes** — this is decisive in practice | no; the solver must be rewritten |
| Total cost | more outer iterations, each expensive | potentially far cheaper; harder to make robust |

**The reduced form dominates industrial practice, and the reason is organisational rather than mathematical**: it lets an existing, trusted, validated simulation code be used unmodified as a function evaluator. The full-space form is the better idea asymptotically — solving the state equation to tight tolerance at a design point that is about to be discarded is obviously wasteful — and it is where the research effort goes.

**Note the parallel with §5.3.** Reduced-space against full-space optimization is the same trade as partitioned against monolithic fluid–structure interaction: reuse existing solvers and iterate, or build one system and solve it once. Both fields reach for the partitioned form first for the same reason, and both find its convergence is the weak point.

## 7.3 What makes an optimization problem hard

Six axes. Most of the field's taxonomy is a product of them.

| Axis | Consequence |
|---|---|
| **Convex or not** | This is the real dividing line, not linear against nonlinear. A convex problem has one minimum and local methods find it globally. A nonconvex one has no such guarantee, and *every* practical method returns a local minimum |
| **Smooth or not** | Non-smooth objectives kill Newton the same way contact does in §1.7. Sources: shocks, limiters in a turbulence model, `if` statements in the residual code, $\ell_1$ terms |
| **Constrained or not** | Constraints add multipliers, complementarity, and the KKT system of §7.5 |
| **Gradients available?** | The single most consequential axis. §8 is about it |
| **Dimension $n$** | $n \sim 10$ admits derivative-free and surrogate methods; $n\sim10^6$ admits only gradient-based ones with limited-memory Hessian approximations |
| **Continuous or discrete** | Integer variables destroy the derivative machinery entirely and move the problem to branch-and-bound |

**On convexity.** It is worth stating plainly that essentially every problem in Parts I and II is nonconvex once wrapped in an objective, and that the standard practice is to accept a local minimum, start from a sensible design, and check that the answer is not an artifact of the starting point. Claims of global optimality in engineering shape or topology optimization should be read as claims about a local minimum unless the problem is provably convex — compliance minimization in topology optimization with a fixed material interpolation being one of the few cases where something can actually be said.

## 7.4 Optimality conditions

**Unconstrained**, the familiar case: $\nabla f(x^\ast) = 0$ and $\nabla^2 f(x^\ast) \succeq 0$.

**Constrained.** Form the Lagrangian

$$\mathcal{L}(x,\lambda,\mu) = f(x) + \lambda^{\mathsf T}c_E(x) + \mu^{\mathsf T}c_I(x)$$

The **Karush–Kuhn–Tucker conditions** at a solution are

$$\nabla_x\mathcal{L} = 0, \qquad c_E = 0, \qquad c_I \le 0, \qquad \mu \ge 0, \qquad \mu_i\,c_{I,i} = 0$$

The first is stationarity — the gradient of the objective is a combination of constraint gradients, so no feasible direction improves $f$. The last is **complementary slackness**: each inequality is either active ($c_i = 0$) or has zero multiplier. Either the constraint is pushing back, or it is not there.

**That last pair of conditions is exactly the Signorini contact condition of §1.7**, with the gap playing the role of $-c_I$ and the contact pressure the role of $\mu$. This is not an analogy. **Frictionless contact is an optimization problem** — minimise potential energy subject to non-penetration — and the active-set, penalty and augmented-Lagrangian methods used for contact in §1.7 are the constrained-optimization methods of §9.2, arriving in structural mechanics under different names.

**The caveat.** KKT are necessary conditions only under a **constraint qualification** (LICQ: the active constraint gradients are linearly independent). Redundant or degenerate constraints — common when geometric constraints are generated automatically — violate it, and the symptom is unbounded or wildly oscillating multipliers rather than an error message.

## 7.5 The KKT system is the fourth saddle point

For the equality-constrained problem, a Newton step on the KKT conditions solves

$$\begin{pmatrix} \mathbf{H} & \mathbf{A}^{\mathsf T} \\ \mathbf{A} & \mathbf{0}\end{pmatrix}\begin{pmatrix}\Delta x \\ \Delta\lambda\end{pmatrix} = \begin{pmatrix}-\nabla_x\mathcal{L} \\ -c_E\end{pmatrix}$$

with $\mathbf{H} = \nabla^2_{xx}\mathcal{L}$ and $\mathbf{A} = \nabla c_E$ the constraint Jacobian.

**This is the same block structure as §2.6 and §4.3**, and it is the fourth arrival at it in this document:

| Arrival | $\mathbf{A}$-block | $\mathbf{B}$-block | The multiplier is |
|---|---|---|---|
| Near-incompressible elasticity (§2.6) | elastic stiffness | volumetric constraint | **pressure** |
| Lagrange-multiplier contact (§2.6, §1.7) | stiffness | non-penetration | **contact pressure** |
| Incompressible Navier–Stokes (§4.3) | momentum operator | divergence | **pressure** |
| **Constrained optimization (here)** | Hessian of the Lagrangian | constraint Jacobian | **the Lagrange multiplier itself** |

**And the same conditions govern all four.** Non-singularity of the KKT matrix requires $\mathbf{A}$ to have full rank (LICQ — the optimization spelling of the **inf-sup condition**) and $\mathbf{H}$ to be positive definite *on the null space of* $\mathbf{A}$ (second-order sufficiency — the optimization spelling of coercivity on the constrained subspace). The block preconditioners of §4.6 apply directly, and so does MINRES, since the matrix is symmetric indefinite whenever $\mathbf{H}$ is symmetric.

**This is the strongest single argument for optimization living in this document rather than beside it.** Four unrelated physical situations produce one algebraic object, and the machinery built for any of them transfers to the rest.

---

# 8. Gradients — the decisive question

Which optimization algorithm is viable is determined almost entirely by how expensive a gradient is. This chapter is the crux of the part.

## 8.1 The four routes, and their costs

Let $n$ be the number of design variables and $m$ the number of scalar outputs whose derivatives are wanted (the objective plus any constraints that depend on the state).

| Method | Cost, in units of one linear solve | Accuracy | Intrusiveness |
|---|---|---|---|
| **Finite differences** | $n$ *nonlinear* solves | poor — half the available digits | none; treats the code as a black box |
| **Complex step** | $n$ solves, in complex arithmetic | machine precision | moderate — the code must be complex-capable throughout |
| **Forward sensitivity** (tangent linear) | $n$ linear solves, **one operator, $n$ right-hand sides** | exact | high — requires $\partial R/\partial u$ and $\partial R/\partial x$ |
| **Adjoint** | $m$ linear solves, **independent of $n$** | exact | high — requires the transposed operator |
| **Algorithmic differentiation** | forward mode $\sim n$; reverse mode $\sim m$ | exact | applied at the source level rather than the equation level |

**The forward/adjoint choice is decided by $n$ against $m$, and nothing else.** Aerodynamic shape optimization has $n \sim 10^3$–$10^6$ design variables and $m \sim 1$–$10$ functionals, so the adjoint wins by orders of magnitude. A model with three design parameters and two hundred monitored outputs is the reverse case, and forward sensitivities win.

## 8.2 Why finite differences are worse than they look

$$\frac{\partial f}{\partial x_i} \approx \frac{f(x + h e_i) - f(x)}{h}$$

Two errors move in opposite directions: truncation error $O(h)$, and roundoff error $O(\epsilon_{\text{mach}}/h)$ from cancelling two nearly-equal numbers. The optimum is around $h \sim \sqrt{\epsilon_{\text{mach}}}$, which delivers roughly **half the digits of the underlying computation** — about eight in double precision, and far fewer if the objective is itself computed by an iterative solver converged to a loose tolerance.

That last point is the one that bites in practice. **A simulation converged to $10^{-6}$ cannot support a finite-difference gradient at all**, because the "noise floor" of the objective exceeds the difference being measured. The symptom is an optimizer that stalls or wanders for no visible reason, and the diagnosis is to plot the objective along one design variable at fine spacing and look for the staircase.

**The complex-step derivative** removes half the problem:

$$\frac{\partial f}{\partial x_i} = \frac{\operatorname{Im}\left[f(x + i h e_i)\right]}{h} + O(h^2)$$

There is no subtraction, so there is no cancellation, and $h$ can be taken as small as $10^{-30}$. It gives an exact derivative at the cost of complex arithmetic through the entire code, and it requires care with any function that branches on comparison or uses `abs`.

## 8.3 The adjoint, derived

This derivation is short and the whole of the adjoint method is in it.

The reduced objective is $\hat J(x) = J(u(x), x)$, with $u(x)$ defined implicitly by $R(u, x) = 0$. Differentiating,

$$\frac{\mathrm{d}\hat J}{\mathrm{d}x} = \frac{\partial J}{\partial x} + \frac{\partial J}{\partial u}\frac{\mathrm{d}u}{\mathrm{d}x}$$

and differentiating the state equation gives the sensitivity of the state,

$$\frac{\partial R}{\partial u}\frac{\mathrm{d}u}{\mathrm{d}x} = -\frac{\partial R}{\partial x} \qquad\Longrightarrow\qquad \frac{\mathrm{d}u}{\mathrm{d}x} = -\left(\frac{\partial R}{\partial u}\right)^{-1}\frac{\partial R}{\partial x}$$

Substituting:

$$\frac{\mathrm{d}\hat J}{\mathrm{d}x} = \frac{\partial J}{\partial x} - \underbrace{\frac{\partial J}{\partial u}\left(\frac{\partial R}{\partial u}\right)^{-1}}_{\text{bracket here}}\ \underbrace{\frac{\partial R}{\partial x}}_{\text{or here}}$$

**Everything turns on which product is formed first.** Bracketing to the right computes $\mathrm{d}u/\mathrm{d}x$, an $n_u \times n$ object requiring $n$ linear solves — that is forward sensitivity. Bracketing to the left computes a single row vector, requiring **one** solve. Define the **adjoint variable** $\lambda$ by

$$\left(\frac{\partial R}{\partial u}\right)^{\mathsf T}\lambda = \left(\frac{\partial J}{\partial u}\right)^{\mathsf T}$$

and the gradient is

$$\frac{\mathrm{d}\hat J}{\mathrm{d}x} = \frac{\partial J}{\partial x} - \lambda^{\mathsf T}\frac{\partial R}{\partial x}$$

**One linear solve, with the transpose of the Jacobian, for the full gradient with respect to every design variable at once.** The entire content of the method is the associativity of matrix multiplication. Everything else — checkpointing, discrete against continuous, turbulence-model adjoints — is engineering around that one observation.

**Two things worth noticing.** The operator is $(\partial R/\partial u)^{\mathsf T}$ — the transpose of the Jacobian that a Newton solve of the forward problem already needs, so a code with a Newton solver is most of the way to an adjoint. And the same $\lambda$ appears in §7.5: **the adjoint variable is the Lagrange multiplier of the state constraint.** The reduced-space gradient and the full-space KKT system are two views of one thing.

## 8.4 Forward sensitivity is the structural-mechanics situation

Worth flagging because it connects two chapters. The forward sensitivity equation

$$\frac{\partial R}{\partial u}\frac{\mathrm{d}u}{\mathrm{d}x_i} = -\frac{\partial R}{\partial x_i}, \qquad i = 1\dots n$$

is **one operator with $n$ right-hand sides**. That is precisely the case §2.6 identified as favouring a **direct** solver: factor once, back-substitute $n$ times. For a linear structural problem with a modest number of design variables and many outputs, forward sensitivities on a stored Cholesky factorization are close to free, and reaching for an adjoint is unnecessary sophistication.

## 8.5 The practical difficulties

**Unsteady adjoints run backward in time.** The adjoint of a time-dependent problem is a terminal-value problem integrated from $T$ to $0$, and at each backward step it needs the *forward* state at that time. Storing the entire forward trajectory is usually impossible. The standard answer is **checkpointing** — store the state at selected times and recompute the intervening segments — with Griewank's `revolve` algorithm giving the optimal schedule for a given memory budget. The cost is a logarithmic factor in recomputation.

**The adjoint of a chaotic system diverges.** Adjoint sensitivities grow like $e^{\lambda_{\max} T}$ with the leading Lyapunov exponent, so for a long-time-averaged quantity in a chaotic flow — which is what a turbulent simulation produces — the computed gradient is enormous and meaningless. This is a genuine obstruction, not an implementation problem, and it is why unsteady adjoint optimization of turbulent flows remains hard. Least-squares shadowing and related methods attack it. **[unverified]** — stated from general knowledge of the literature; no source checked here.

**Discrete against continuous adjoint** — and this is the same choice as the consistent tangent of §2.5.

| | |
|---|---|
| **Continuous adjoint** | derive the adjoint PDE analytically, then discretize it. Cleaner derivation, insight into boundary conditions, discretization independent of the forward code |
| **Discrete adjoint** | transpose the discrete Jacobian of the code as written. Harder to implement, but gives the **exact** gradient of the discrete objective the optimizer is actually minimising |

**The optimizer minimises the discrete objective, so it wants the discrete gradient.** A continuous adjoint gives a gradient that is consistent with the continuous problem and slightly inconsistent with the discrete one, and the mismatch shows up as an optimizer that converges to a tolerance and then stalls, or that fails a line search near the optimum. The principle is the one already stated in §2.5: **the derivative must be the derivative of the discrete operator you are actually applying, not of the continuous one you have in mind.**

**Non-differentiable objectives.** Shock position, limiter switches, turbulence-model clipping, remeshing, and any branch in the residual code all make $\hat J$ non-smooth. The gradient still evaluates to a number; it is simply the wrong number, and there is no warning. The standard check is a **dot-product test** — verify $\lambda^{\mathsf T}(\partial R/\partial x)\,v$ against a directional finite difference for a random $v$ — which catches both implementation errors and genuine non-smoothness.

**The frozen-turbulence approximation** deserves naming because it is so common. Differentiating a RANS solver requires differentiating the turbulence model; many implementations skip it, holding $\nu_t$ fixed. The resulting gradient is wrong, sometimes in sign, and the practice persists because the full adjoint of a turbulence model is substantial work. **[unverified]** — the characterisation of how wrong is recalled, not verified.

## 8.6 Algorithmic differentiation

The same two modes, applied to source code rather than to equations.

| | Cost | Memory | Corresponds to |
|---|---|---|---|
| **Forward mode** | $O(n)$ passes | small | forward sensitivity |
| **Reverse mode** | $O(m)$ passes | **large** — the tape | the adjoint |

Reverse mode records every intermediate operation on a tape during the forward pass and replays it backward. **The memory cost is the tape**, and checkpointing is the same remedy as for the unsteady adjoint above — indeed it is the same problem, seen at a lower level.

**AD applied to a whole simulation code is usually a mistake.** Differentiating through an iterative solver's convergence history is wasteful and can be inaccurate; the better construction differentiates the *converged* residual, using the implicit-function argument of §8.3 and applying AD only to the residual assembly. That is the standard structure of a modern discrete adjoint: AD for the local pieces, the implicit-function theorem for the solve.

---

# 9. The algorithms

## 9.1 Unconstrained

| Method | Uses | Cost per iteration | When |
|---|---|---|---|
| **Gradient descent** | $\nabla f$ | one gradient | almost never on its own — convergence rate depends on the condition number of the Hessian, and engineering problems are badly conditioned |
| **Nonlinear CG** | $\nabla f$ | one gradient | large $n$, cheap gradients, no Hessian storage |
| **Newton** | $\nabla f$, $\nabla^2 f$ | Hessian + solve | quadratic convergence near the solution; the Hessian is $n^2$ and indefinite far from it |
| **Quasi-Newton (BFGS)** | $\nabla f$ | rank-2 update | the standard for moderate $n$; builds curvature from gradient history |
| **L-BFGS** | $\nabla f$ | stores $\sim10$ vector pairs | **the workhorse for large $n$**; the default for PDE-constrained problems |
| **Gauss–Newton / Levenberg–Marquardt** | residual Jacobian | one solve | least-squares objectives specifically — approximates the Hessian as $\mathbf{J}^{\mathsf T}\mathbf{J}$, which is free once $\mathbf{J}$ is available |

**Two globalization strategies** turn a locally-convergent method into a reliable one, and they are alternatives rather than complements:

- **Line search** — pick a direction, then a step length satisfying sufficient-decrease and curvature conditions (Wolfe). This is the same construction as §2.5's line search in structural Newton iterations.
- **Trust region** — pick a radius within which the local model is believed, minimise the model inside it, then expand or shrink the radius based on how well the model predicted the actual decrease. More robust with indefinite Hessians, because it does not require a descent direction to exist.

**A connection worth making explicit.** For a quadratic $f(x) = \frac{1}{2}x^{\mathsf T}\mathbf{A}x - b^{\mathsf T}x$ with $\mathbf{A}$ symmetric positive definite, minimising $f$ and solving $\mathbf{A}x = b$ are the same problem, and **the conjugate gradient method of §2.6 is an optimization algorithm.** The linear solver and the optimizer are one object; they only look different because in Parts I and II the quadratic is never written down.

## 9.2 Constrained

| Method | Idea | Where it is strong |
|---|---|---|
| **Penalty** | add $\frac{1}{2\rho}\|c\|^2$ to the objective; drive $\rho\to0$ | simple; the subproblem becomes ill-conditioned as $\rho\to0$, which is the method's fatal flaw |
| **Augmented Lagrangian** | penalty **plus** an explicit multiplier estimate, updated each outer iteration | exact convergence without $\rho\to0$. The same construction as augmented-Lagrangian contact in §1.7 |
| **Sequential quadratic programming (SQP)** | at each iterate solve a QP with a quadratic model of $\mathcal{L}$ and linearised constraints | **few, expensive function evaluations** — the PDE-constrained case. Each QP's optimality system is the KKT block of §7.5 |
| **Interior point / barrier** | replace $c_I\le0$ with $-\mu\sum\log(-c_i)$ in the objective; solve a sequence with $\mu\to0$ | **very large $n$ and many inequality constraints**. IPOPT and its relatives |
| **Active set** | guess which inequalities are active, solve the equality-constrained problem, correct the guess | excellent **warm starts** — successive similar problems. Identical in spirit to active-set contact |

**Two observations that tie this back.**

**The barrier parameter is a continuation.** Interior-point methods solve an easy problem (large $\mu$, solution far inside the feasible region) and follow the solution as $\mu$ decreases — the central path. That is structurally the same device as **load incrementation** in §2.5 and **pseudo-time marching to steady state** in §4.4: make the problem easy, then deform it to the one you want, using each solution as the initial guess for the next. Continuation is one of the small number of ideas that recurs across every chapter of this document.

**Three of the five methods above appear in §1.7 under other names.** Penalty, augmented Lagrangian and active set are the standard contact treatments. Structural mechanics has been doing constrained optimization for decades without labelling it that way, because contact *is* a constrained minimisation of potential energy.

## 9.3 Derivative-free and global methods

Legitimate when the gradient is genuinely unavailable — a non-differentiable objective, a legacy binary, a physical experiment in the loop, discrete variables — and a poor choice otherwise.

| Method | Practical dimension limit |
|---|---|
| Nelder–Mead simplex | $n \lesssim 10$, and no convergence theory |
| Pattern / direct search | $n \lesssim 20$, with convergence guarantees under mild conditions |
| **Bayesian optimization / surrogate-based** | $n \lesssim 20$, but tolerates **very expensive** objectives — tens to hundreds of evaluations total |
| Genetic and evolutionary algorithms, particle swarm | popular, easy to implement, weak guarantees, and cost that scales badly with $n$ |

**On evolutionary methods, stated fairly.** They handle discrete variables, discontinuous objectives and multiple objectives naturally, and they parallelise trivially, which is a real advantage when a function evaluation is an hour of compute on its own node. What they do not do is compete with a gradient-based method on a smooth problem with an available adjoint — there the difference is not a constant factor but a scaling in $n$. The failure mode to watch for is a genetic algorithm used on a smooth differentiable problem because building the adjoint looked like more work.

**Multi-objective** problems deserve one note: with competing objectives there is no single optimum but a **Pareto front**, and the deliverable is the front rather than a point. Weighted-sum scalarisation is the cheap route and misses non-convex portions of the front; $\epsilon$-constraint and Pareto-tracking methods do not.

## 9.4 Discrete and mixed-integer

Integer design variables — how many stiffeners, which material from a catalogue, ply angles from a discrete set, on/off decisions — remove the derivative machinery entirely. The tools are **branch and bound** over a relaxation, **mixed-integer linear** or **conic programming** when the model can be cast that way, and heuristics when it cannot. This is a genuinely different field from everything else in this document and is noted rather than developed.

---

# 10. Where the layer meets Parts I and II

## 10.1 Structural optimization — three levels

The standard hierarchy, in increasing order of design freedom and difficulty.

| Level | Design variables | Difficulty |
|---|---|---|
| **Sizing** | thicknesses, cross-sectional areas, ply counts | modest $n$, fixed mesh, gradients straightforward |
| **Shape** | boundary geometry, control points of a spline | needs a **mesh deformation** scheme and design velocity fields; mesh quality degrades as the shape moves |
| **Topology** | material presence or absence at every point | largest $n$ (one variable per element), and the answer is not a shape but a field to be interpreted |

### SIMP, and a self-adjoint gradient

The dominant topology formulation is **SIMP** — solid isotropic material with penalisation. Introduce a density field $\rho_e \in [0,1]$ per element and interpolate the stiffness as

$$E(\rho_e) = \rho_e^{\,p}\,E_0, \qquad p \approx 3$$

The exponent penalises intermediate densities: at $p=3$ a half-dense element gives an eighth of the stiffness for half the material, so intermediate values are uneconomical and the optimizer drives the field toward $0$ or $1$. Minimise compliance $c = \mathbf{f}^{\mathsf T}\mathbf{u}$ subject to a volume constraint.

**Compliance minimisation is self-adjoint**, and this is worth seeing. With $J = \mathbf{f}^{\mathsf T}\mathbf{u}$ and $R = \mathbf{K}\mathbf{u} - \mathbf{f}$, the adjoint equation of §8.3 is $\mathbf{K}^{\mathsf T}\lambda = \mathbf{f}$, and since $\mathbf{K}$ is symmetric this is the state equation itself — so $\lambda = \mathbf{u}$ and **no additional solve is needed at all.** The sensitivity reduces to

$$\frac{\partial c}{\partial \rho_e} = -p\,\rho_e^{\,p-1}\,\mathbf{u}_e^{\mathsf T}\mathbf{k}_0\,\mathbf{u}_e$$

an element-local quantity computed from the displacement field already in hand. **The gradient with respect to a million design variables costs nothing beyond the forward solve**, which is why topology optimization became practical early and why its textbook implementations are famously short.

**Checkerboarding, again — but a different one.** Unfiltered SIMP produces alternating solid and void elements, because the discretization overestimates the stiffness of that pattern. The remedy is a **sensitivity or density filter** over a radius, which also imposes a minimum length scale and removes mesh dependence. This resembles the pressure checkerboarding of §4.3 in appearance and in being a spurious mode the discrete operator cannot see, but the mechanism differs — there it is an inf-sup failure of the mixed formulation, here it is an artificial stiffness of the element pattern. **Worth not conflating.**

**Level-set methods** are the main alternative: represent the boundary implicitly and evolve it with a Hamilton–Jacobi equation driven by the shape derivative. Crisper boundaries, no intermediate densities, and a harder time creating new holes.

## 10.2 Aerodynamic shape optimization

**This is the application the adjoint method was made for**, and the numbers show why: $10^3$–$10^6$ shape design variables against a single objective (drag), with lift, pitching moment, volume and thickness as constraints. Forward sensitivities would need one linear solve per design variable; the adjoint needs one per functional.

The practical difficulties, all of which are §8.5 made concrete:

- **Mesh deformation** must be smooth, robust to large shape changes, and itself differentiable — the sensitivity of the mesh to the design variables is part of $\partial R/\partial x$
- **The turbulence model adjoint**, and the frozen-turbulence temptation
- **Shocks**, where the discrete objective is non-smooth in the design variables
- **Geometric constraints** from manufacturing, which are frequently degenerate and violate LICQ (§7.4)

## 10.3 Inverse problems and parameter estimation

The same machinery pointed at a different question. Instead of designing a device, infer a field or parameter set from observations:

$$\min_x\ \frac{1}{2}\|\,\mathcal{O}(u(x)) - d\,\|^2 + \alpha\,\mathcal{R}(x)$$

with $\mathcal{O}$ an observation operator, $d$ the data, and $\mathcal{R}$ a **regularisation** term. The regularisation is not optional: inverse problems are typically ill-posed, the data constrain only a low-dimensional subspace of $x$, and without $\mathcal{R}$ the optimizer fits noise. Tikhonov ($\|x\|^2$ or $\|\nabla x\|^2$) is the standard choice, with $\alpha$ selected by the L-curve or by discrepancy against the known noise level.

**The Gauss–Newton Hessian is diagnostic here.** Its eigenvalue spectrum shows how many directions in design space the data actually resolve; the rest are determined by the regulariser alone. Reporting a fitted parameter without checking that it lies in the resolved subspace is a standard error.

**Variational data assimilation (4D-Var)** is exactly this construction with time in it, and its adjoint is the backward-in-time adjoint of §8.5 — which is where checkpointing was largely developed.

## 10.4 The relation to training neural networks

Worth one paragraph because it is the version of this material most readers have met first, and because the differences are instructive.

The differentiation machinery is identical: backpropagation **is** reverse-mode algorithmic differentiation, and the adjoint of §8.3 is the same computation. The optimization regime is completely different. Objectives are stochastic (a minibatch gives a noisy gradient), $n$ is $10^9$ or more, the problem is wildly nonconvex, and the objective is cheap to evaluate but evaluated an enormous number of times. That regime favours **SGD, Adam and their relatives** — first-order methods with per-parameter step sizes — where PDE-constrained optimization, with expensive and accurate gradients and $n$ merely large, favours **L-BFGS and SQP**. The methods diverge not because the mathematics differs but because the cost model does.

---

# 11. What the optimization layer does to the chain

Returning to §0.1. Adding an objective does not append a step; it adds a layer above every step, and each acquires a dual.

| Step | The forward problem | What optimization adds |
|---|---|---|
| **1. Continuum statement** | balance, kinematics, constitutive | an **objective functional** and constraints; the design parameterisation |
| **2. PDE system** | the closed forward system | the **adjoint PDE**, whose operator is the formal adjoint of the linearised forward operator |
| **3. Character** | hyperbolic / parabolic / elliptic | the adjoint inherits the type but **reverses the direction of information flow** |
| **4. Discretization** | FEM, FVM, explicit or implicit | the choice between the **discretization of the adjoint** and the **adjoint of the discretization** |
| **5. Algebraic system** | $\mathbf{K}\mathbf{u}=\mathbf{f}$, or a nonlinear residual | the **KKT saddle point** of §7.5 |
| **6. Solver** | direct or iterative, Newton or Picard | an **outer optimizer** wrapping the inner solve, with inexactness now permissible at both levels |

**Step 3 is the one that surprises people.** The adjoint of an advection equation advects in the *opposite* direction; the adjoint of a parabolic problem is parabolic **backward in time**, which is why an unsteady adjoint is a terminal-value problem and why checkpointing exists. Elliptic problems are self-adjoint or nearly so, which is why elliptic adjoints are the easy case — and why compliance minimisation in §10.1 turned out to need no extra solve at all. **The difficulty of the adjoint is predictable from the character of the forward problem**, which is the classification of §0.3 paying off one more time.

**Step 4 restates a principle already met twice.** The consistent tangent of §2.5 and the discrete adjoint of §8.5 are the same demand: **the derivative must be the derivative of the discrete operator actually applied.** Anything else converges to a different problem, or fails to converge at all near the solution.

**Step 6 introduces a genuinely new degree of freedom: inexactness.** In Parts I and II the linear solve sits inside a Newton iteration, and §4.6 noted that over-solving it is waste. Optimization adds another level of the same idea — the *state* equation need not be solved tightly at a design point that will be discarded, nor the gradient computed to full precision when the design step is large. **Inexact and multilevel optimization methods** exploit this deliberately, and the full-space formulation of §7.2 is the limiting case in which the state equation is never solved exactly until the end.

---

# 12. What is not here

Named so the gaps are deliberate rather than accidental.

**Optimization is no longer absent** — it is Part IV, added after the two fields were written. It is a **layer** rather than a field, and the distinction is worth applying to anything added next: a field brings its own governing equations, a layer wraps a field and changes what the existing equations are used for. Uncertainty quantification and reduced-order modelling are layers in the same sense and would belong in the same position; the entries below are fields.

**Other fields**, in rough order of how much they would add:

- **Heat transfer and multiphysics coupling** — parabolic, and the simplest place to discuss coupling strategies without FSI's added-mass difficulty
- **Electromagnetics** — Maxwell's equations, and the curl-conforming (edge) elements that a naive nodal discretization gets wrong
- **Porous media and geomechanics** — Darcy and Biot poroelasticity; another saddle point, and a fourth arrival at the same structure
- **Molecular and particle methods** — MD, DEM, SPH; the discrete alternative to a continuum entirely
- **Plasma physics** — kinetic, hybrid and MHD models. The `fusion` repository is the reference and this document borrows from it

**Cross-cutting topics deliberately omitted**, each of which would apply to every section: verification and validation, uncertainty quantification, mesh generation and adaptivity, parallel decomposition, and reduced-order modelling.

---

# 13. Unverified claims

Collected per the convention in `fusion`'s `CLAUDE.md`. Each is a statement made above that could not be checked against a primary source while writing.

1. **§2.2** — that Rayleigh damping's mass-and-stiffness-proportional form is chosen for algebraic convenience rather than physical fidelity. Standard commentary, not checked against a source.
2. **§2.6** — attribution of FETI to structural mechanics (Farhat and Roux, c. 1991). Recalled, not verified.
3. **§3.5** — the DNS cost scalings, $\mathrm{Re}^{9/4}$ in grid points and $\mathrm{Re}^{3}$ overall. These are the figures usually quoted and the exponents depend on the definition of $\mathrm{Re}$ and on whether the time step is included. Treat as orders of magnitude.
4. **§3.5** — the claim that wall-resolved LES at $y^+\sim1$ puts the first cell in the micrometre range. Correct in order of magnitude for industrial Reynolds numbers, but no specific case was worked.
5. **§5.3** — that weakly coupled partitioned FSI schemes are *unconditionally* unstable at comparable densities, rather than merely conditionally so. The strong form of the claim is recalled from the literature and should be checked.
6. **Throughout** — the characterisations of what commercial codes do by default (generalized-$\alpha$ in implicit structural codes, Rhie–Chow in unstructured finite-volume codes, direct solvers as the structural default). These reflect general practice as understood here and have not been checked against documentation.
7. **§8.5** — that adjoint sensitivities of a chaotic system diverge as $e^{\lambda_{\max}T}$ with the leading Lyapunov exponent, and that this is an obstruction in principle rather than an implementation defect. Stated from general familiarity with the literature; no source checked.
8. **§8.5** — the claim that the frozen-turbulence approximation can produce a gradient wrong in sign, not merely in magnitude. Recalled, not verified, and the answer plausibly depends on the case.
9. **§9.3, §10.4** — the comparative statements about method families (evolutionary methods scaling badly in $n$ against gradient-based ones; first-order stochastic methods dominating in neural-network training for cost-model reasons). These are consensus positions rather than results checked here.

**What would settle any of these:** a single primary source each. None is load-bearing for the argument of the document; they are supporting detail, and are marked so a later reader does not quote them as established.

---

# 14. Sources

**Structural mechanics and FEM**

- Hughes T.J.R., *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis* — the reference for the variational formulation, element technology, and Newmark-family integration
- Belytschko T., Liu W.K., Moran B., *Nonlinear Finite Elements for Continua and Structures* — finite strain, explicit dynamics, contact
- Simo J.C., Hughes T.J.R., *Computational Inelasticity* — return mapping and the consistent tangent
- Zienkiewicz O.C., Taylor R.L., *The Finite Element Method* — the standard multi-volume treatment
- Bathe K.-J., *Finite Element Procedures* — time integration and eigensolution in particular
- Crisfield M.A., *Non-linear Finite Element Analysis of Solids and Structures* — arc-length methods

**Fluid dynamics and CFD**

- Ferziger J.H., Perić M., Street R.L., *Computational Methods for Fluid Dynamics* — the general reference, and the standard treatment of SIMPLE-family coupling
- Versteeg H.K., Malalasekera W., *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*
- LeVeque R.J., *Finite Volume Methods for Hyperbolic Problems* — conservation, Godunov's theorem, limiters
- Toro E.F., *Riemann Solvers and Numerical Methods for Fluid Dynamics*
- Pope S.B., *Turbulent Flows* — the closure problem and the model hierarchy
- Wesseling P., *Principles of Computational Fluid Dynamics*

**Numerical linear algebra and solvers**

- Saad Y., *Iterative Methods for Sparse Linear Systems*
- Elman H., Silvester D., Wathen A., *Finite Elements and Fast Iterative Solvers* — the definitive treatment of saddle-point preconditioning for both Stokes and elasticity
- Davis T.A., *Direct Methods for Sparse Linear Systems* — ordering, supernodal and multifrontal factorization
- Kelley C.T., *Iterative Methods for Linear and Nonlinear Equations* — Newton–Krylov, inexact Newton
- Trottenberg U., Oosterlee C., Schüller A., *Multigrid*

**Classification and time integration**

- Hairer E., Wanner G., *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems* — the definition of stiffness, and the DAE index
- Ascher U.M., Petzold L.R., *Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations*

**Optimization**

- Nocedal J., Wright S.J., *Numerical Optimization* — the standard reference for everything in §7 and §9: KKT conditions, line search and trust region, quasi-Newton, SQP, interior point
- Boyd S., Vandenberghe L., *Convex Optimization* — convexity, duality, and why the convex/nonconvex line is the one that matters
- Biegler L.T., Ghattas O., Heinkenschloss M., van Bloemen Waanders B. (eds.), *Large-Scale PDE-Constrained Optimization* — the reduced-space against full-space question of §7.2
- Griewank A., Walther A., *Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation* — forward and reverse mode, and the `revolve` checkpointing schedule
- Giles M.B., Pierce N.A., "An Introduction to the Adjoint Approach to Design," *Flow, Turbulence and Combustion* **65**, 393 (2000) — the clearest short treatment of the discrete/continuous adjoint question
- Jameson A., "Aerodynamic Design via Control Theory," *Journal of Scientific Computing* **3**, 233 (1988) — the origin of adjoint aerodynamic shape optimization
- Bendsøe M.P., Sigmund O., *Topology Optimization: Theory, Methods and Applications* — SIMP, filtering, and the self-adjoint compliance result of §10.1
- Martins J.R.R.A., Ning A., *Engineering Design Optimization* — a modern engineering-facing treatment covering §8 and §10 together
- Tarantola A., *Inverse Problem Theory* — regularisation and resolution for §10.3

**In the `fusion` repository**

- `notes/theory/pde_character.md` — the three classes worked in full, stiffness, and what forces an implicit method. §0.3–§0.5 above condense it
- `notes/theory/mhd_model.md` §1 — the Navier–Stokes equations, the compressible and incompressible systems, the equation count, and the evolved/constrained axis. §3.2–§3.3 above condense it
- `notes/theory/hierarchy.md` — the model ladder from kinetic to fluid, and the closure problem in its general form
