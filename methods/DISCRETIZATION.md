# Discretization — three ways to turn a PDE into a matrix

**One 2D domain, one equation, three methods, worked end to end. What each one actually asserts, where they coincide, and the four tests that separate them.**

---

> **Companion to `COMPUTATIONAL.md`.** That document says *which* method each field uses and why — §2.1 for finite elements in structural mechanics, §4.1 for finite volumes in fluid dynamics. It does not say what each method *does*, because a map should not contain a derivation. This document is the derivation, and it exists for the reason `fusion`'s `pde_character.md` exists: both sections need it and neither owns it.
>
> **Everything numerical here was computed, not recalled.** The scripts are in `code/` and §9 is the run table. Any number in the text can be traced to a script and re-run. Where a claim could not be tested, it is marked **[unverified]** and collected in §11.
>
> **Scope.** Finite difference, finite volume, and finite element, on a scalar second-order problem in two dimensions. Spectral, discontinuous Galerkin and meshfree methods are named in §8 and not developed.

---

# 0. The problem

The same domain throughout: the **unit square** $\Omega = (0,1)^2$, and steady heat conduction,

$$-\nabla\cdot\left(k\,\nabla T\right) = f \quad\text{in }\Omega$$

with $T$ prescribed on the boundary. It is **elliptic** — no time derivative, and the second-order coefficients $\operatorname{diag}(k,k)$ have eigenvalues of one sign, which `COMPUTATIONAL.md` §0.3 works out against the parabolic and hyperbolic cases — and it is the simplest problem on which all three methods apply without qualification or special pleading.

**Every symbol in it.**

| | meaning | units |
|---|---|---|
| $T(x,y)$ | **temperature** — the unknown, the field being solved for | K |
| $k$ | **thermal conductivity** — how readily the material carries heat | W m$^{-1}$K$^{-1}$ |
| $f(x,y)$ | **volumetric source** — heat generated per unit volume per unit time, given | W m$^{-3}$ |
| $\mathbf{q} = -k\nabla T$ | **heat flux** — energy crossing unit area per unit time, not solved for but implied | W m$^{-2}$ |

$\mathbf{q} = -k\nabla T$ is **Fourier's law**: heat flows down the temperature gradient, at a rate proportional to it, and the minus sign is what makes it flow from hot to cold. $k$ is in general a positive function of position, and in an anisotropic material a tensor rather than a scalar. Ellipticity requires $k \ge k_{\min} > 0$ — if $k$ could vanish somewhere, the equation would stop constraining $T$ there. In this document $k$ is constant and equal to 1 everywhere except §4.3, which makes it discontinuous, and §5, which uses the interface that creates.

**$k$ is a function of position here, never of $T$.** Real conductivities do depend on temperature, and allowing $k = k(T)$ makes the equation **nonlinear** — the operator applied to $T$ now depends on $T$, so there is no fixed matrix, and everything §3 establishes by comparing matrices entry by entry lapses. Every method in this document would then be applied inside an outer iteration rather than once. `TIME_INTEGRATION.md` §0.5 raises the same assumption from the other side, counting a temperature-dependent coefficient among the modelling errors that mesh refinement cannot touch.

**What "steady" means, and what it does not.** The unsteady equation is

$$\underbrace{\rho c_p\frac{\partial T}{\partial t}}_{\text{accumulation}} \;-\; \underbrace{\nabla\cdot(k\,\nabla T)}_{\text{conduction}} \;=\; \underbrace{f}_{\text{generation}}$$

with two further symbols and one further independent variable:

| | meaning | units |
|---|---|---|
| $t$ | **time**, the independent variable that §2–§5 do not have | s |
| $\rho$ | **density** of the material | kg m$^{-3}$ |
| $c_p$ | **specific heat capacity** — energy to raise unit *mass* by one degree | J kg$^{-1}$K$^{-1}$ |
| $\rho c_p$ | **volumetric heat capacity** — energy to raise unit *volume* by one degree | J m$^{-3}$K$^{-1}$ |

**Where $f$ acts: throughout the interior, at every point.** It is *volumetric*, W m$^{-3}$, a property of the material at each location rather than something imposed at the edge. Heat supplied at the boundary is a different object entirely — a **Neumann condition** $-k\,\partial T/\partial n = q_{\text{wall}}$, which §4.2 covers. This document's $f$ happens to be positive everywhere inside and to vanish on all four sides, but that is a property of the particular choice made below, not a requirement on $f$.

Read term by term it is an energy budget for a point. **Accumulation** is the rate at which thermal energy is being stored there, fast if the material is dense and hard to heat. **Conduction** is the net energy arriving by conduction — the divergence of the flux, with the sign such that $-\nabla\cdot\mathbf{q}$ counts energy flowing *in*. **Generation** is energy appearing from some other process entirely: resistive heating, a chemical reaction, absorbed radiation. All three terms carry the same units, W m$^{-3}$, each being a volumetric power density; the equation says stored equals received plus produced.

**Dividing through by $\rho c_p$ collects the entire material into one number.** With $k$ constant,

$$\frac{\partial T}{\partial t} = \alpha\nabla^2 T + \frac{f}{\rho c_p}, \qquad \boxed{\;\alpha = \frac{k}{\rho c_p}\;}$$

where $\alpha$ is the **thermal diffusivity**, units m$^2$s$^{-1}$. It is a ratio, and reading it as one is the point: conduction over storage. A material that carries heat well but also stores a great deal of it spreads a disturbance slowly — copper and water differ far more in $\alpha$ than in $k$ alone. Length² over time is the signature of a diffusion coefficient, and it means $\alpha$ sets a time scale $L^2/\alpha$ for a domain of size $L$ rather than a speed.

**That equation with $f = 0$ is what is normally called the heat equation**, and it is `TIME_INTEGRATION.md` §2's subject, written there as $u_t = \alpha u_{xx}$ in one dimension. It is also where `CLAUDE.md`'s Fourier number $r = \alpha\Delta t/\Delta x^2$ comes from. **Six equations, one family**, each obtained from the one above by assuming something away:

| | name | what is assumed | used in |
|---|---|---|---|
| $\rho c_p T_t - \nabla\cdot(k\nabla T) = f$ | general heat conduction | nothing | §0, as the parent |
| $T_t = \alpha\nabla^2T + f/\rho c_p$ | heat equation with a source | $k$ constant in space | — |
| $T_t = \alpha\nabla^2T$ | **the heat equation** | and $f = 0$ | `TIME_INTEGRATION.md` §2 |
| $-\nabla\cdot(k\nabla T) = f$ | steady conduction | $\partial_t = 0$, $k$ free again | **this document** |
| $-\nabla^2T = f/k$ | **Poisson** | and $k$ constant | §2 onward, with $k = 1$ |
| $-\nabla^2T = 0$ | **Laplace** | and $f = 0$ | — |

**What conduction is, physically.** Energy transfer through matter that does not itself move. Molecules at a hot spot vibrate or travel faster than their neighbours, collisions pass that kinetic energy along, and statistically it flows from fast to slow. In a solid the carriers are lattice vibrations — phonons — together with free electrons in a metal, which is why metals conduct heat and electricity well for the same underlying reason. What defines conduction is the absence of bulk transport: **convection** carries energy by physically moving the fluid, **radiation** carries it electromagnetically and needs no medium at all, and conduction moves energy while the material stays put. This equation models conduction alone, which is why §4.4's advection test is a genuine departure rather than a variation — it adds the transport mechanism the equation was built without.

**What conduction does with the energy: nothing. It only moves it.** The three terms have three different jobs and it is worth keeping them apart. Generation *creates* thermal energy from something else — electrical, chemical, nuclear. Accumulation is the only term that *converts* energy into temperature, with $\rho c_p$ as the exchange rate. Conduction is pure *transport*: it changes where the energy is, never how much of it there is, and never directly what temperature anything is. A point warms through conduction only by way of accumulation — energy arrives, accumulation converts it, temperature rises.

**Over the whole domain conduction is zero-sum.** By the divergence theorem,

$$\int_\Omega \nabla\cdot\mathbf{q}\;\mathrm{d}\Omega = \oint_{\partial\Omega}\mathbf{q}\cdot\mathbf{n}\;\mathrm{d}s$$

everything conduction does internally cancels — what leaves one point enters its neighbour — and the only net effect is what crosses the boundary. **Total energy inside changes solely through generation and boundary flux, never through interior conduction.** That is exactly the statement §5 tests, and why finite volume's cell balances telescope. What conduction changes is the *distribution*: it flattens gradients, moving energy from hot to cold and never the reverse, which is the Second Law appearing in this equation and the reason diffusion has an arrow of time.

**How it is observed depends on whether the problem is steady.** Unsteady, it shows up as temperature changing in place — a thermocouple in a cooling bar, the reading dropping. But the thermocouple is reading accumulation; conduction is inferred from the drop being fastest where the curvature is greatest. Steady is the sharper case: **nothing in the domain changes at all.** A wall's temperature field is fixed, every point has zero accumulation, nothing is stored or released anywhere inside it — and energy flows through continuously. The observable is not in the wall. It is the power the heater must supply, exactly $\oint\mathbf{q}\cdot\mathbf{n}$, forever, to hold the room steady. The wall is a conduit, and what one measures is the flux passing through it, with a heat-flux sensor or as $k\,\Delta T/L$.

**The conduction term is an imbalance, not a rate of flow.** Fourier's law gives the flux at a point; the equation contains its *divergence*, which measures how much more leaves than arrives. Heat can pour through a location at enormous rate and contribute nothing, if inflow and outflow match — the steady wall again. So $-\nabla\cdot\mathbf{q}$ is the net gain, and a point warms only when it receives more than it passes on.

**With $k$ constant, that imbalance is a comparison with the neighbourhood.** Then $-\nabla\cdot(k\nabla T) = -k\nabla^2T$, and the Laplacian compares a point to the average around it: $\nabla^2T > 0$ means $T$ sits *below* its local average. The budget reads: **a point warms in proportion to how much colder it is than its surroundings.** Conduction is an averaging operator and diffusion is everything being dragged toward its local mean. That is also the physical content of §2.1's matrix, before any of it is derived — in two dimensions

$$\nabla^2T \;\approx\; \frac{T_E + T_W + T_N + T_S - 4T_P}{h^2} \;=\; \frac{4}{h^2}\left(\overline{T}_{\text{neighbours}} - T_P\right)$$

so the row $4T_P - T_E - T_W - T_N - T_S$ is not an arbitrary stencil. It is how far a point falls below the average of its four neighbours, and the $h^2$ is what turns that comparison into a curvature.

**What flows is energy, not temperature.** $\mathbf{q}$ is a *heat* flux, W m$^{-2}$ — joules per second across unit area. Temperature does not flow and cannot: it is an intensive property, the **potential** that drives the transport, not a transportable quantity. Energy is what is conserved and what moves. Holding the two apart is what explains why $\rho c_p$ has to be in the equation at all: **the conservation law is about energy, the unknown is temperature, and $\rho c_p$ is the exchange rate between them.**

**And heat conduction is one instance of a structure that is everywhere in physics.** Every conservation law has the same shape,

$$\frac{\partial\,(\text{density})}{\partial t} + \nabla\cdot(\text{flux}) = \text{source}$$

and in each case a **constitutive law** — an empirical statement about the material, not a conservation principle — supplies the flux as a coefficient times the gradient of a potential:

| conserved quantity | potential | constitutive law | coefficient |
|---|---|---|---|
| energy | temperature $T$ | **Fourier** $\mathbf{q} = -k\nabla T$ | conductivity $k$ |
| species mass | concentration $c$ | **Fick** $\mathbf{J} = -D\nabla c$ | diffusivity $D$ |
| charge | voltage $\phi$ | **Ohm** $\mathbf{J} = -\sigma\nabla\phi$ | conductivity $\sigma$ |
| fluid in a porous medium | pressure $p$ | **Darcy** $\mathbf{u} = -(\kappa/\mu)\nabla p$ | permeability over viscosity |
| momentum | velocity $u$ | **Newton** $\tau = \mu\,\partial u/\partial y$ | viscosity $\mu$ |

Every row is minus a coefficient times a gradient, and the minus sign always means downhill. Substitute any of them into the conservation law and the same second-order equation comes out.

**So this document is not really about heat.** Steady, source-free and constant-coefficient is $-\nabla^2\phi = 0$: electrostatics in a charge-free region, groundwater through uniform soil, steady diffusion of a solute. Add a source and it is Poisson — $-\nabla^2\phi = \rho/\varepsilon_0$ is §0's equation with different letters. **The three methods, §3's equivalence, §4's four tests and §5's conservation result all carry over unchanged**, because none of them ever uses the fact that the potential is a temperature. Heat is the vocabulary, not the restriction.

Momentum is the row that does not transfer cleanly: its flux is a tensor rather than a vector, and the transported quantity is also what does the transporting, which is where the nonlinearity of Navier–Stokes enters. §4.4's advection test is the first step in that direction.

**Steady state is $\partial T/\partial t = 0$**, which deletes the accumulation term and leaves the equation at the top. **The source survives.**

**It is worth being exact about which term goes away.** Conduction is not "the steady term" and accumulation "the unsteady one" — both are active at every instant, and heat conducts vigorously while a body is warming up. The only asymmetry is that accumulation is the term carrying $\partial_t$. What happens as a system settles is that **accumulation decays to zero** while conduction and generation keep working and come into balance with each other. Steady state is not conduction taking over from accumulation; it is accumulation running out. And it is approached asymptotically rather than reached: the decay below never formally arrives at zero. Steady does not mean nothing is happening; it means nothing changes in time. The physical picture is heat generated at rate $f$ and conducted away exactly as fast, at every point, indefinitely — a resistor held at constant power against a heat sink, source running, temperature field fixed.

Dropping $f$ as well would give $-\nabla\cdot(k\nabla T) = 0$, the **Laplace equation**: steady *and* source-free. Here $f = 2\pi^2\sin\pi x\sin\pi y$ is strictly positive inside, the boundary values are zero, and the solution peaks at 1 in the middle — the bump exists precisely because the source is on.

**So with no source at all, where does any temperature come from?** There are only two other places it can come from, and removing $f$ leaves exactly those two.

**From the boundary.** Steady and sourceless is the case of a wall with a warm room on one side and cold air on the other. Nothing generates heat inside the wall and nothing accumulates, yet there is a full temperature field and a steady flux passing through it. That is most of building insulation, heat exchangers, and any thermal-resistance calculation — Laplace is not a degenerate case but the commonest one in engineering practice. Its solutions are the **harmonic** functions, and they obey the **maximum principle**: the interior temperature lies strictly between the boundary extremes, with no interior peak, because sustaining a peak would require a source. Our solution has an interior peak of exactly 1 with zero boundary data, which is the maximum principle failing in the only way it can — the source.

**From the past.** The unsteady equation with zero boundary values and no source keeps everything it was given in its initial condition and spreads it out until it is gone. A bar heated once and then left alone.

**And the two combine exactly, because the problem is linear.** Any steady problem splits as

$$T = T_{\text{source}} + T_{\text{boundary}}$$

where $T_{\text{source}}$ solves $-\nabla^2T = f$ with *zero* boundary values and $T_{\text{boundary}}$ solves $-\nabla^2T = 0$ with the real ones. §0 picks the first extreme deliberately: homogeneous Dirichlet means every part of the solution is attributable to $f$, with nothing entering from the edge. A Laplace problem is the opposite extreme. Real problems are a sum of the two, and studying Laplace alone is worthwhile precisely because it is the half carrying all the geometry and boundary behaviour with nothing else in the way — which is what §4.1 and §4.2 go on to disturb.

**And steady state is a limit rather than a separate problem.** Start the unsteady equation from any initial field with this document's $f$ and boundary values, and $T$ approaches the steady solution; the difference decays like $e^{-2\pi^2\alpha t}$, because the slowest mode that fits in the square with zero edges is $\sin\pi x\sin\pi y$ and $-\nabla^2$ multiplies it by $2\pi^2$. With $\alpha = 1$ the time constant is $1/2\pi^2 \approx 0.0507$, so the square is within a thousandth of steady by $t \approx 0.35$. **The $2\pi^2$ is the same number that appears in $f$** — the manufactured source of §0 and the decay rate of the unsteady problem are one eigenvalue seen twice, once as a multiplier and once as a rate.

Two consequences run through the rest of the document. **The balance §5 measures is local and nonzero:** for every cell, flux out minus flux in equals $\int_V f$, which is why §2.2's right-hand side is what it is. Were steady the same as source-free, every cell would balance against zero and there would be nothing to conserve. And **there is no time step anywhere here:** no $\partial_t$ means no stability limit, no CFL condition, no explicit-against-implicit question. That axis is the subject of `TIME_INTEGRATION.md`, and removing it is what allows §2 through §5 to be about the spatial discretization alone.

**Why the divergence form and not $-k\nabla^2 T = f$.** The two are not interchangeable. Expanding the divergence,

$$-\nabla\cdot(k\,\nabla T) = -k\,\nabla^2 T - \nabla k\cdot\nabla T$$

so they agree only when $\nabla k = 0$. The Laplacian form is the **constant-coefficient special case**, and with $k$ constant it reduces further to the Poisson equation $-\nabla^2 T = f/k$. Two reasons to keep the general form as the starting point. It is the one that survives §4.3, where $k$ jumps across an interface and $\nabla k$ does not exist as a function at all — there the Laplacian form is not merely inconvenient but meaningless. And it is the one finite volume needs: §2.2 integrates the equation over a cell and applies the divergence theorem, which requires a divergence to apply it to. The flux is the physical object, $-\nabla\cdot\mathbf{q} = f$ is the energy budget above with the accumulation term deleted — what is produced must be what leaves, and everything in §5 is about which methods preserve that statement discretely.

**Why this problem and not something more interesting.** Because the differences between the methods are invisible on it, and that is the point. §3 shows that on a uniform mesh with constant $k$ the three produce *literally the same matrix*. Everything that distinguishes them appears only when one assumption is removed at a time, which is what §4 does — one test per assumption, each isolating one difference.

**Manufacturing a solution, and manufacturing is the operative word.** No solution is found here. One is *chosen*, and the problem is then built backwards around it. Take

$$T(x,y) = \sin(\pi x)\,\sin(\pi y), \qquad k \equiv 1$$

and define $f$ by applying the operator to $T$. With $k$ constant the operator is the Laplacian, so

$$f = -\nabla\cdot(k\nabla T) = -\nabla^2 T = \pi^2\sin(\pi x)\sin(\pi y) + \pi^2\sin(\pi x)\sin(\pi y) = 2\pi^2\sin(\pi x)\sin(\pi y)$$

one $\pi^2$ from each second derivative. $T$ vanishes on all four sides, so the boundary condition is homogeneous Dirichlet, $T = 0$ on $\partial\Omega$. This pair is used for every convergence measurement in §3.

**The point of manufacturing rather than solving** is that the exact answer is then known at every point by construction, so a computed error is the true error and an observed order of convergence is a measurement rather than an estimate. The cost is that $T$ is as smooth as a function can be, which flatters every method equally — §4 exists because of that.

---

## 0.1 The same operator in fluid dynamics, electromagnetism and plasma physics

Four families where §0's equation is not an analogy but the thing itself.

| | the equation | which part is §0's operator |
|---|---|---|
| **Navier–Stokes**, viscous term | $\rho(\mathbf{u}_t + \mathbf{u}\cdot\nabla\mathbf{u}) = -\nabla p + \mu\nabla^2\mathbf{u}$ | $\mu\nabla^2\mathbf{u}$, one Laplacian per velocity component |
| **Navier–Stokes**, pressure | $\nabla^2p = -\rho\,\nabla\cdot(\mathbf{u}\cdot\nabla\mathbf{u})$ | the whole equation: Poisson |
| **Maxwell**, electrostatics | $-\nabla\cdot(\varepsilon\nabla\phi) = \rho_{\text{charge}}$ | the whole equation, $\varepsilon$ where $k$ was |
| **Maxwell**, magnetostatics | $\nabla^2\mathbf{A} = -\mu_0\mathbf{J}$ | one Poisson per component of $\mathbf{A}$ |
| **MHD**, induction | $\mathbf{B}_t = \nabla\times(\mathbf{u}\times\mathbf{B}) + \eta\nabla^2\mathbf{B}$ | $\eta\nabla^2\mathbf{B}$, resistive diffusion |
| **MHD**, equilibrium | $\Delta^*\psi = -\mu_0R^2p'(\psi) - FF'(\psi)$ | $\Delta^*$, a modified Laplacian |

**Viscosity is momentum diffusion.** The kinematic viscosity $\nu = \mu/\rho$ plays exactly the part $\alpha$ plays above, same units m$^2$s$^{-1}$, and treating the term implicitly gives $(\mathbf{I} - \Delta t\,\nu\nabla^2)\mathbf{u} = \text{rhs}$ — a **screened** Poisson problem rather than a pure one, better conditioned than the pure operator because the identity dominates at small $\Delta t$.

**The pressure Poisson equation is the larger appearance.** Take the divergence of the momentum equation and impose $\nabla\cdot\mathbf{u} = 0$; pressure is then determined by an elliptic problem with no time derivative in it. Every projection or fractional-step incompressible solver solves it at least once per time step, and it is usually the dominant cost of the entire simulation — an elliptic solve embedded inside a problem that is otherwise hyperbolic and parabolic. Stokes flow, at vanishing Reynolds number, is elliptic outright, and the streamfunction–vorticity formulation $\nabla^2\psi = -\omega$ is Poisson again.

**Electrostatics is §0 with the letters changed**, permittivity $\varepsilon$ where conductivity $k$ was and charge density where the heat source was. That is the field solve in an electrostatic particle-in-cell code, and it is why a block-structured AMR framework ships a multigrid solver at all. Particle-in-cell codes also run a Poisson solve for **divergence cleaning**, projecting $\nabla\cdot\mathbf{E}$ back onto Gauss's law after a field update has drifted off it.

**In MHD the diffusion is the heat equation outright.** With $\eta = 1/\mu_0\sigma$ the **magnetic diffusivity**, units m$^2$s$^{-1}$, a magnetic field diffuses through a conductor exactly as temperature diffuses through a solid. The ratio of the two terms in the induction equation is the **magnetic Reynolds number** $R_m = uL/\eta$ — which is §4.4's Péclet number with different letters: advection against diffusion, the same competition, and the same numerical trouble when it is large. Equilibrium is elliptic too: the **Grad–Shafranov** equation is §0's operator with a source that depends on the solution, which is the $k = k(T)$ nonlinearity of §0 arrived at from the source side rather than the coefficient side.

**And magnetized plasma is where the tensor coefficient of §0 stops being a footnote.** Thermal conduction there is

$$\nabla\cdot\Big(\kappa_\parallel\,\hat{\mathbf{b}}\hat{\mathbf{b}}\cdot\nabla T + \kappa_\perp(\mathbf{I} - \hat{\mathbf{b}}\hat{\mathbf{b}})\cdot\nabla T\Big)$$

with $\hat{\mathbf{b}}$ the local field direction — §0's divergence form with $k$ a tensor built from the magnetic field. Heat runs almost freely along field lines and barely crosses them: the anisotropy ratio $\kappa_\parallel/\kappa_\perp$ reaches $10^7$–$10^{10}$ in tokamak conditions, and estimates at the plasma edge run higher still. **Three consequences follow, and all three are §0's material.** The divergence form is not a stylistic preference but the only available one, since a tensor coefficient admits no Laplacian form at all. The condition number of the discrete operator scales with the anisotropy ratio, so the linear solve is near-singular by construction. And a small discretization error in the parallel direction is magnified by that ratio into the perpendicular dynamics, which are the ones of physical interest — the reason this is a live numerical-methods problem rather than a solved one.

---

# 1. What each method actually asserts

**This is the whole conceptual content of the comparison, and everything else is consequence.** The three methods do not differ in accuracy, or in which equation they solve. They differ in **which form of the same conservation statement they take as primitive**, and therefore in what the discrete unknowns mean and what is satisfied exactly.

| | **Finite difference** | **Finite volume** | **Finite element** |
|---|---|---|---|
| **Starts from** | the **differential** form, at points | the **integral** form, over cells | the **weak** form, against test functions |
| **The statement made** | $-\nabla^2T = f$ holds at each node, with derivatives replaced by difference quotients | $\oint_{\partial V}(-k\nabla T)\cdot\mathbf{n}\,\mathrm{d}S = \int_V f$ holds **exactly** on each cell | $\int_\Omega k\nabla T\cdot\nabla v = \int_\Omega fv$ holds for every $v$ in a finite space |
| **The unknown is** | a **point value** $T_{ij} \approx T(x_i,y_j)$ | a **cell average** $\bar T_I \approx \frac{1}{\|V_I\|}\int_{V_I}T$ | a **basis coefficient**; $T_h = \sum_J T_J\varphi_J$ is a function defined everywhere |
| **Exactly satisfied** | nothing — it is a truncation | **the balance on every cell**, and hence on every union of cells | orthogonality of the residual to the trial space (Galerkin) |
| **The approximation is in** | the difference quotient | the **face flux** — and nowhere else | the finiteness of the space |
| **Needs** | a structured grid, or a mapping to one | a partition into cells, any shape | a mesh and a basis, any shape |

**Read the fourth row again.** Finite volume does not approximate the divergence theorem — it *uses* it. The volume integral of a divergence becomes a sum of surface integrals identically, so the only thing left to approximate is the flux through each face. That single structural fact is the origin of every property finite volumes are prized for, and §5 measures it.

**And the third row.** The three unknowns are different objects. A point value and a cell average of the same function differ by $O(h^2)$, so on a smooth problem this is invisible in the convergence rate — but it is not invisible in what "the error" means, and §3 has to define it separately for each method before it can be measured.

---

# 2. The three, worked on a $4\times4$ mesh

![The same 4x4 discretization, three ways](figs/fig1_three_meshes.png)

Same domain, same resolution, three different sets of unknowns. The red markings show the support of one equation in each case — the five nodes of a stencil, the four faces of a control volume, the six triangles surrounding a node. **The numbers are the unknown ordering**, $i$ running fastest, and every matrix printed below is indexed by them: row 4 of the finite-difference system is the equation at the centre node, row 5 of the finite-volume system is the shaded cell.

Note the degree-of-freedom counts already differ: **9 for finite difference and finite element** (the interior nodes of a $5\times5$ grid; the 16 boundary nodes are prescribed) against **16 for finite volume** (every cell is an unknown; there are no boundary cells to prescribe, only boundary *faces*). That asymmetry never goes away and is a recurring source of confusion when comparing "the same mesh" between codes.

## 2.1 Finite difference

Take Taylor expansions about the node $(x_i, y_j)$ with uniform spacing $h$:

$$T_{i+1,j} = T_{ij} + h\,\partial_xT + \frac{h^2}{2}\partial_x^2T + \frac{h^3}{6}\partial_x^3T + O(h^4)$$

Adding the expansion for $T_{i-1,j}$ cancels every odd derivative:

$$T_{i+1,j} + T_{i-1,j} = 2T_{ij} + h^2\,\partial_x^2T + O(h^4) \quad\Longrightarrow\quad \partial_x^2 T = \frac{T_{i+1,j} - 2T_{ij} + T_{i-1,j}}{h^2} + O(h^2)$$

**The second-order accuracy comes from the cancellation of the odd terms**, which requires the stencil to be symmetric and the spacing uniform. Both of those will be revisited: symmetry fails for advection (§4.4), and uniformity fails on a stretched grid (§5).

Doing the same in $y$ and substituting gives the **5-point stencil**:

$$\frac{4T_{ij} - T_{i+1,j} - T_{i-1,j} - T_{i,j+1} - T_{i,j-1}}{h^2} = f(x_i, y_j)$$

One equation per interior node, 9 of them on the $4\times4$ mesh. Nodes adjacent to the boundary reference a prescribed value, which moves to the right-hand side.

**The right-hand side is the point value of $f$.** Remember that; it is one of the three things that still differ after §3 proves the matrices identical.

**Assembled.** Nine equations, nine unknowns. Write $\mathbf{A}$ for the discrete operator — the matrix with $\mathbf{A}\mathbf{T} \approx -\nabla^2 T$, carrying $4/h^2$ on the diagonal and $-1/h^2$ off it. Every matrix in this section is printed as $h^2\mathbf{A}$, which clears the $h^2$ and leaves integers, and it is the form in which the three methods can be set side by side. Written out in full — unknown vector, right-hand side, and the numbers it evaluates to. Both vectors are indexed the way the figure numbers the interior nodes, $i$ running fastest:

```
 4 -1  0 | -1  0  0 |  0  0  0       T0         h^2 f0      0.616850
-1  4 -1 |  0 -1  0 |  0  0  0       T1         h^2 f1      0.872358
 0 -1  4 |  0  0 -1 |  0  0  0       T2         h^2 f2      0.616850
---------+----------+---------       --         ------      --------
-1  0  0 |  4 -1  0 | -1  0  0       T3         h^2 f3      0.872358
 0 -1  0 | -1  4 -1 |  0 -1  0   x   T4    =    h^2 f4  =   1.233701
 0  0 -1 |  0 -1  4 |  0  0 -1       T5         h^2 f5      0.872358
---------+----------+---------       --         ------      --------
 0  0  0 | -1  0  0 |  4 -1  0       T6         h^2 f6      0.616850
 0  0  0 |  0 -1  0 | -1  4 -1       T7         h^2 f7      0.872358
 0  0  0 |  0  0 -1 |  0 -1  4       T8         h^2 f8      0.616850
```

$T_I$ is the unknown at node $I$ and $f_I = f(x_i,y_j)$ is the source sampled at the same node — one index set, used for the rows of the matrix, the entries of the unknown vector, and the entries of the right-hand side.

**Where that last column comes from.** The nodes are $x_i = i/4$ and $y_j = j/4$ for $i,j \in \{1,2,3\}$, and the unknown at that node is component

$$I = 3(j-1) + (i-1)$$

of the vector. Every entry of the source comes from one formula — substitute the node coordinates into $f$:

$$f_{ij} = f(x_i, y_j) = 2\pi^2\,\sin\frac{i\pi}{4}\,\sin\frac{j\pi}{4}, \qquad i, j \in \{1,2,3\}$$

Laid out on the mesh this is the array below — **rows are $j$ and columns are $i$**, because $x$ is horizontal in the figure and $i$ indexes $x$. That is the structured-grid convention and it is the transpose of matrix indexing, where $i$ would be the row. Reading left to right then top to bottom gives $I$ in order:

$$\mathbf{f} = 2\pi^2\begin{pmatrix}
\sin\frac{\pi}{4}\sin\frac{\pi}{4} & \sin\frac{\pi}{2}\sin\frac{\pi}{4} & \sin\frac{3\pi}{4}\sin\frac{\pi}{4}\\[3pt]
\sin\frac{\pi}{4}\sin\frac{\pi}{2} & \sin\frac{\pi}{2}\sin\frac{\pi}{2} & \sin\frac{3\pi}{4}\sin\frac{\pi}{2}\\[3pt]
\sin\frac{\pi}{4}\sin\frac{3\pi}{4} & \sin\frac{\pi}{2}\sin\frac{3\pi}{4} & \sin\frac{3\pi}{4}\sin\frac{3\pi}{4}
\end{pmatrix}$$

The $x$ argument advances across a row and the $y$ argument down a column. Only two sine values occur, $\sin\tfrac{\pi}{4} = \sin\tfrac{3\pi}{4} = \tfrac{\sqrt2}{2}$ and $\sin\tfrac{\pi}{2} = 1$, so the nine entries take three values: $\pi^2$ at the four corner nodes, $\sqrt2\,\pi^2$ at the four edge nodes, $2\pi^2$ at the centre. Multiplying by $h^2 = 1/16$ gives the last column of the block above.

**Nothing in that construction asked what $f$ does between the nodes.** A spike halfway between nodes 3 and 4 would change the true solution and leave every entry of $\mathbf{f}$ untouched. §2.2 and §2.3 both integrate instead, and that is the difference.

**Only the middle row carries the full five-point stencil.** The other eight are missing one or two $-1$'s, and that is exactly where a prescribed boundary value was eliminated. The $4\times4$ mesh is the coarsest on which an uncontaminated interior row exists at all.

**The missing entries and the right-hand side are the same bookkeeping seen from two sides.** The equation at an interior node is

$$4T_{ij} - T_{i+1,j} - T_{i-1,j} - T_{i,j+1} - T_{i,j-1} = h^2f_{ij} + \!\!\!\sum_{\text{neighbours on }\partial\Omega}\!\!\! g$$

so each $-1$ dropped from a row adds the prescribed value $g$ at that neighbour to the same row's right-hand side. Here $g = 0$ everywhere, the sum vanishes, and the right-hand side is nine values of $f$ and nothing else — a property of this boundary data, not of the method. Row 4 is the only row with neither a missing entry nor a boundary term; row 0 lost two entries and with nonzero data would gain two.

## 2.2 Finite volume

Integrate the equation over one cell $V_{ij} = [x_i, x_{i+1}]\times[y_j,y_{j+1}]$ and apply the divergence theorem:

$$\int_{V_{ij}} -\nabla\cdot(k\nabla T)\,\mathrm{d}V = \oint_{\partial V_{ij}} (-k\nabla T)\cdot\mathbf{n}\,\mathrm{d}S = \int_{V_{ij}} f\,\mathrm{d}V$$

**Nothing has been approximated yet.** This is an exact statement about the exact solution, and it is an exact statement about *any* cell, of any shape. The four faces give

$$\underbrace{q_E + q_W + q_N + q_S}_{\text{outward fluxes}} = \int_{V_{ij}} f\,\mathrm{d}V$$

Now the single approximation: the flux through the east face, of area $h$ (in 2D, a length), evaluated from the two adjacent cell averages,

$$q_E = -k\left.\frac{\partial T}{\partial x}\right|_{E}\cdot h \;\approx\; -k\,\frac{\bar T_{i+1,j} - \bar T_{ij}}{h}\cdot h = -k\left(\bar T_{i+1,j} - \bar T_{ij}\right)$$

Summing the four and dividing by $h^2$ recovers the same five coefficients as §2.1. But **the right-hand side is now $\int_V f$, the cell integral**, not a point value.

**Assembled.** Sixteen cells, sixteen equations. No scaling is needed — the face flux $-k(\bar T_E - \bar T_P)$ carries no $h$, so the assembled balance matrix already *is* $h^2\mathbf{A}$ and is directly comparable with §2.1. Written out as §2.1 was, the unknown column written $T_I$ for the cell average $\bar T_I$, and $b_I = \int_{V_I} f$ the source integrated over that cell:

```
 6 -1  .  . | -1  .  .  . |  .  .  .  . |  .  .  .  .     T0      b0      0.171573
-1  5 -1  . |  . -1  .  . |  .  .  .  . |  .  .  .  .     T1      b1      0.414214
 . -1  5 -1 |  .  . -1  . |  .  .  .  . |  .  .  .  .     T2      b2      0.414214
 .  . -1  6 |  .  .  . -1 |  .  .  .  . |  .  .  .  .     T3      b3      0.171573
------------+-------------+-------------+------------     ---     ---     --------
-1  .  .  . |  5 -1  .  . | -1  .  .  . |  .  .  .  .     T4      b4      0.414214
 . -1  .  . | -1  4 -1  . |  . -1  .  . |  .  .  .  .     T5      b5      1.000000
 .  . -1  . |  . -1  4 -1 |  .  . -1  . |  .  .  .  .     T6      b6      1.000000
 .  .  . -1 |  .  . -1  5 |  .  .  . -1 |  .  .  .  .     T7      b7      0.414214
------------+-------------+-------------+------------     ---     ---     --------
 .  .  .  . | -1  .  .  . |  5 -1  .  . | -1  .  .  .  x  T8   =  b8   =  0.414214
 .  .  .  . |  . -1  .  . | -1  4 -1  . |  . -1  .  .     T9      b9      1.000000
 .  .  .  . |  .  . -1  . |  . -1  4 -1 |  .  . -1  .     T10     b10     1.000000
 .  .  .  . |  .  .  . -1 |  .  . -1  5 |  .  .  . -1     T11     b11     0.414214
------------+-------------+-------------+------------     ---     ---     --------
 .  .  .  . |  .  .  .  . | -1  .  .  . |  6 -1  .  .     T12     b12     0.171573
 .  .  .  . |  .  .  .  . |  . -1  .  . | -1  5 -1  .     T13     b13     0.414214
 .  .  .  . |  .  .  .  . |  .  . -1  . |  . -1  5 -1     T14     b14     0.414214
 .  .  .  . |  .  .  .  . |  .  .  . -1 |  .  . -1  6     T15     b15     0.171573
```

**The diagonal takes three values and the reason is geometric.** A face shared with another cell contributes $k h/h = 1$; a face on the wall contributes $2k$, because the prescribed value sits at the half-cell distance $h/2$ rather than $h$. So corner cells carry 6, edge cells 5, and only the four cells of the middle $2\times2$ block carry 4. **Twelve of the sixteen rows are boundary-modified**, which is the asymmetry §4.2 describes, here as a count.

**Where that last column comes from.** $f$ is separable, so the cell integral factors into two one-dimensional integrals and the $2\pi^2$ collapses — each integral carries a $1/\pi$:

$$b_I = \int_{V_I} f = 2\pi^2\int_{x_i}^{x_{i+1}}\!\!\sin\pi x\,\mathrm{d}x \int_{y_j}^{y_{j+1}}\!\!\sin\pi y\,\mathrm{d}y = 2\,C_i\,C_j, \qquad C_i = \cos\pi x_i - \cos\pi x_{i+1}$$

Cell $(i,j)$ spans $[x_i, x_{i+1}]\times[y_j, y_{j+1}]$ for $i,j \in \{0,1,2,3\}$, and its average is component

$$I = 4j + i$$

of the vector — no shift, because no cell was eliminated. Writing it out in the index, with $x_i = i/4$:

$$b_{ij} = 2\left[\cos\frac{i\pi}{4} - \cos\frac{(i{+}1)\pi}{4}\right]\left[\cos\frac{j\pi}{4} - \cos\frac{(j{+}1)\pi}{4}\right], \qquad i, j \in \{0,1,2,3\}$$

The four intervals give

$$\mathbf{C} = \left(1-\tfrac{\sqrt2}{2},\;\; \tfrac{\sqrt2}{2},\;\; \tfrac{\sqrt2}{2},\;\; 1-\tfrac{\sqrt2}{2}\right)$$

so $\mathbf{b}$ is twice the outer product $\mathbf{C}\otimes\mathbf{C}$, written here on the mesh rather than as a 16-row column because the rank-one structure is the content — rows are $j$, columns are $i$, reading order gives $I$:

$$\mathbf{b} = 2\begin{pmatrix}
(1-\tfrac{\sqrt2}{2})^2 & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & (1-\tfrac{\sqrt2}{2})^2\\[3pt]
(1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & \tfrac12 & \tfrac12 & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2}\\[3pt]
(1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & \tfrac12 & \tfrac12 & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2}\\[3pt]
(1-\tfrac{\sqrt2}{2})^2 & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & (1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} & (1-\tfrac{\sqrt2}{2})^2
\end{pmatrix}$$

Three distinct entries: $2(1-\tfrac{\sqrt2}{2})^2 = 3-2\sqrt2$ on a corner cell, $2(1-\tfrac{\sqrt2}{2})\tfrac{\sqrt2}{2} = \sqrt2-1$ on an edge cell, and $2(\tfrac{\sqrt2}{2})^2 = 1$ exactly on an interior cell — $b_0$ through $b_{15}$ of the block above.

**Sixteen unknowns and no boundary elimination.** Every cell is an unknown, so nothing was removed from the vector and no prescribed value moved to the right-hand side. The wall enters through the diagonal instead: where §2.1 lost a $-1$ and gained a boundary term, §2.2 keeps every row full-width and pays with a larger diagonal. Same information, opposite bookkeeping.

**Two properties follow immediately and neither depends on the flux approximation being good.** The flux leaving cell $i$ through a face is, by construction, **the same number** that enters its neighbour, with the opposite sign. So summing the discrete balance over any set of cells makes every interior face cancel, leaving only the fluxes through the outer boundary of that set. Conservation holds on every cell, on every patch, and on the whole domain — exactly, at any resolution, on any mesh. §5 measures it.

## 2.3 Finite element

Start from the **strong form** — the equation itself, $-\nabla\cdot(k\nabla T) = f$, required to hold at every point of $\Omega$. Multiply by a test function $v$ vanishing where $T$ is prescribed, integrate over $\Omega$, and integrate by parts:

$$\int_\Omega -\nabla\cdot(k\nabla T)\,v\,\mathrm{d}\Omega = \int_\Omega k\,\nabla T\cdot\nabla v\,\mathrm{d}\Omega - \oint_{\partial\Omega}\left(k\frac{\partial T}{\partial n}\right)v\,\mathrm{d}S$$

giving the **weak form**: find $T$ such that

$$\int_\Omega k\,\nabla T\cdot\nabla v\,\mathrm{d}\Omega = \int_\Omega f\,v\,\mathrm{d}\Omega + \oint_{\partial\Omega}\left(k\frac{\partial T}{\partial n}\right)v\,\mathrm{d}S \qquad\text{for all admissible } v$$

**What $v$ is.** It is not part of the answer and it is never solved for. Write the **residual** $r = -\nabla\cdot(k\nabla T) - f$, the amount by which a candidate $T$ fails to satisfy the equation. The strong form demands $r = 0$ at every point. The weak form demands only

$$\int_\Omega r\,v\,\mathrm{d}\Omega = 0$$

for every admissible $v$ — that $r$ average to zero against every probe. Allowed *all* $v$, the two demands are equivalent. **Restricting $v$ to finitely many probes is the discretization**, and it is why §1's table says finite elements exactly satisfy orthogonality of the residual to the trial space rather than the equation itself.

**Here $v$ is a pyramid.** The probe attached to node $I$ is the function $\varphi_I$ that equals 1 at node $I$, 0 at every other node, and is linear on each triangle. On this mesh it has a closed form. In local coordinates $\xi = (x - x_I)/h$ and $\eta = (y - y_I)/h$,

$$\varphi_I = \max\big(0,\; 1 - \max(|\xi|,\; |\eta|,\; |\xi - \eta|)\big)$$

which is a pyramid of height 1 standing on the hexagonal patch, with six flat triangular faces and zero outside. Face by face around the node the pieces are $1-\xi$, $1-\eta$, $1+\xi-\eta$, $1+\xi$, $1+\eta$, $1-\xi+\eta$ — each one linear, each one 1 at the centre and 0 along the outer edge of its triangle.

So the $I$th equation reads: the residual, weighted by that pyramid, averages to zero over the patch. Nine nodes, nine pyramids, nine equations. It is a local weighted average, not a point check.

**The third term in that maximum is the mesh diagonal.** $|\xi|$ and $|\eta|$ survive every sign flip; $|\xi-\eta|$ survives $\xi\leftrightarrow\eta$ and not $\xi\leftrightarrow-\eta$. One term in the formula for a single basis function is the whole of §3.3.

**The pyramid does two jobs.** As a probe it weights the residual, which is the paragraph above. It is also the building block of the answer: $T_h = \sum_J T_J\varphi_J$, so summing the pyramids with the computed coefficients gives a function defined at *every* point of $\Omega$, not only at nodes. That is why §1's third row calls the unknown a basis coefficient rather than a point value — $T_I$ is not a sample of $T$, it is how much of pyramid $I$ goes into the mix. Using one family for both jobs is exactly what **Galerkin** names.

**It is not a spike.** $\varphi_I$ spans two mesh widths, rising from zero at each neighbour node, so it overlaps its neighbours' pyramids. That overlap is why $\mathbf{K}$ has off-diagonal entries at all: $K_{IJ} \ne 0$ precisely when nodes $I$ and $J$ share a triangle. A function confined to one cell would give a diagonal $\mathbf{K}$ and no coupling between nodes.

**Three requirements on any such family, and they are what rule choices in or out.** Local support, or $\mathbf{K}$ is dense. One derivative, or the integral $\int\nabla\varphi_I\cdot\nabla\varphi_J$ does not exist. And $\sum_J\varphi_J \equiv 1$ — the **partition of unity** — or a constant cannot be represented and the method cannot be consistent. The pyramids satisfy the third to $1.1\times10^{-16}$, which is just the statement that they interpolate the value 1 at every node.

**The choice used here has a name: the P1 element**, piecewise linear on triangles, also called the Courant element after its 1943 appearance. **$P$ is for polynomial** — $\mathbb{P}_k$ denotes polynomials of total degree at most $k$ well outside finite elements, and the finite-element usage inherits it through the French school, *polynômes*. **$Q$ is simply the next letter**, taken for the tensor-product space because $P$ was already spoken for; the gloss "quadrilateral" is a mnemonic attached afterwards and an imperfect one, since $Q_k$ is used on hexahedra in 3D where nothing is quadrilateral, and $Q_1$ is bilinear rather than quadratic. **[unverified]** — the $P$ half is standard, the origin of $Q$ is inference.

**In both, the subscript is the degree the space runs up to**, and the difference is how degree is counted: $P_k$ is total degree $\le k$ on simplices — triangles in 2D — while $Q_k$ is degree $\le k$ *in each variable separately* on boxes. So $P_1$ spans $\{1, x, y\}$, three coefficients matching a triangle's three vertices; $Q_1$ spans $\{1, x, y, xy\}$, four matching a square's four corners; $P_2$ spans $\{1,x,y,x^2,xy,y^2\}$, six matching three vertices plus three edge midpoints.

**The name belongs to the local space, the pyramid to the global one.** Restricted to any single triangle, $\varphi_I$ is a member of $P_1$ — that is the whole content of the label. Glue the six such pieces together so they agree across shared edges and the pyramid is what comes out. P1 is the rule applied element by element; $\varphi_I$ is what the rule produces on the mesh.

**Other choices, and the two that appear later in this document:**

| element | shape functions | global basis on a uniform mesh | support | appears in |
|---|---|---|---|---|
| **P1** on triangles | the three barycentric coordinates $\lambda_a$ | $\max\big(0,\,1-\max(\lvert\xi\rvert,\lvert\eta\rvert,\lvert\xi-\eta\rvert)\big)$ | 6 triangles, a hexagon | §2.3, §3 |
| **Q1** on squares | $\tfrac14(1+\xi_a\xi)(1+\eta_a\eta)$ on the reference square | $(1-\lvert\xi\rvert)(1-\lvert\eta\rvert)$ for $\lvert\xi\rvert,\lvert\eta\rvert \le 1$ | 4 squares, a $2h\times2h$ box | §3.2 |
| **P2** on triangles | $\lambda_a(2\lambda_a-1)$ at vertices, $4\lambda_a\lambda_b$ at edge midpoints | no single closed form; 6 nodes per triangle | 6 triangles plus edge nodes | not used here |

**The Q1 basis is the tensor product of two one-dimensional hats**, a bilinear tent on a $2h\times2h$ square rather than a hexagon. It reaches the four diagonal neighbours that P1 misses, and its stiffness matrix is the nine-point stencil

$$\mathbf{K}_{Q1} = \frac{1}{3}\begin{pmatrix}-1 & -1 & -1\\ -1 & 8 & -1\\ -1 & -1 & -1\end{pmatrix}$$

whose diagonal $8/3$ against P1's 4 is the 1.3333 that §3 measures. Both integrate to $h^2$ and both satisfy the partition of unity, checked to $2.2\times10^{-16}$. P2 is the standard route to third-order accuracy in $L^2$ and is included for the naming only — it is not assembled anywhere in this document. **[unverified]**

**One more family is already in this document under another name.** The piecewise constant $P_0$, the indicator of a cell, fails the one-derivative requirement and so cannot be a trial space for this weak form — but as a *test* weight it is exactly the finite-volume choice of §3.1, and the Dirac spike there is the collocation limit. The three methods of §2 are three points in this one catalogue.

**And $v$ vanishes where $T$ is prescribed** because there is nothing to test there — the value is already known, so no equation is wanted — which is also what removes the boundary integral over that part of $\partial\Omega$.

In structural terms $v$ is a **virtual displacement** and the weak form is the statement of virtual work, which is where the vocabulary of §2.3 comes from. §3.1 makes the other half of the point: finite difference and finite volume are the same construction with a spike and a cell indicator in place of the hat.

**Using a test function does not by itself mean Galerkin.** The general construction — require $\int_\Omega r\,w_I = 0$ for a finite family of weights $w_I$ — is the **method of weighted residuals**, and Galerkin is one choice within it: the choice to weight with the very functions used to build the solution. The alternatives are equally legitimate:

| weight | name | where it appears |
|---|---|---|
| $w_I$ = the trial basis itself | **Galerkin**, strictly Bubnov–Galerkin | §2.3, this section |
| $w_I \ne$ the trial basis | **Petrov–Galerkin** | SUPG, §4.4 |
| $w_I = \delta(\mathbf{x} - \mathbf{x}_I)$ | **collocation** | finite difference, §3.1 |
| $w_I = \mathbb{1}_{V_I}$ | **subdomain collocation** | finite volume, §3.1 |
| $w_I = \partial r/\partial T_I$ | **least squares** | least-squares FEM, not used here |

§3.1's table of weights is this family, which is why all three methods of §2 fall inside it.

**And "Galerkin" and "finite element" are two independent choices, not a framework and an instance of it.** One answers *how the approximation space is built*, the other *how the equation is tested*, and neither implies the other:

| | Galerkin testing | collocation testing |
|---|---|---|
| **finite-element basis** — piecewise polynomials on a mesh | standard FEM, this document | orthogonal collocation FEM |
| **global basis** — Fourier, Chebyshev | spectral Galerkin | pseudospectral, spectral collocation |

Finite elements say nothing about testing: the term is a statement about local support and piecewise polynomials. Galerkin says nothing about the basis: it applies to Fourier modes or radial basis functions just as well. **The full name "Galerkin finite element method" is used precisely because both halves are decisions.**

This document contains the evidence that they separate. §4.4's SUPG is finite elements *without* Galerkin — the same mesh and the same hats for the trial space, with test functions deliberately skewed upstream to stabilise advection. §8's spectral methods are Galerkin *without* finite elements. So the layering is: **weighted residuals is the framework**, Galerkin is one instance of its testing choice, finite elements one instance of its basis choice, and what §2.3 builds is the pairing of the two.

**Two things happened in that one line, and both are load-bearing.**

**One derivative moved from $T$ to $v$.** The strong form needs $T$ twice differentiable; the weak form needs it once. That is what allows piecewise-linear basis functions, whose second derivative does not exist.

**The boundary term is the flux.** A prescribed heat flux — a Neumann condition — enters as a known term on the right-hand side and requires nothing further. It is a **natural** boundary condition. In finite differences the same condition requires ghost nodes or one-sided differences; here it falls out of the integration by parts. §4.2 returns to this.

**Discretizing.** Choose $T_h = \sum_J T_J\varphi_J$, the same pyramids, and take $v = \varphi_I$ in turn (**Galerkin**: the test space equals the trial space). The result is $\mathbf{K}\mathbf{T} = \mathbf{F}$ with

$$K_{IJ} = \int_\Omega k\,\nabla\varphi_I\cdot\nabla\varphi_J\,\mathrm{d}\Omega, \qquad F_I = \int_\Omega f\,\varphi_I\,\mathrm{d}\Omega$$

**$\mathbf{K}$ is the stiffness matrix and $\mathbf{F}$ the load vector**, and the names are worth stating once because every later section uses them. $K_{IJ}$ measures how much the $I$th equation resists a unit of the $J$th coefficient: it is nonzero only where $\varphi_I$ and $\varphi_J$ overlap, and its size is set by how steeply they vary there. The vocabulary is structural — for a spring, force $=$ stiffness $\times$ displacement, and $\mathbf{K}$ is the discrete form of exactly that — and it survives into problems like this one where the physical quantity is a conductance and nothing is stiff.

Computed element by element and **assembled**: each triangle contributes a $3\times3$ block to the rows and columns of its own three nodes. For a linear (P1) triangle of area $A_e$ with vertices $(x_a,y_a)$, the gradients of the shape functions are constant, and

$$K^e_{ab} = \frac{k}{4A_e}\left(\beta_a\beta_b + \gamma_a\gamma_b\right), \qquad \beta_a = y_b - y_c,\quad \gamma_a = x_c - x_b$$

with $(a,b,c)$ cyclic. **$\mathbf{K}$ is symmetric by inspection** — it is symmetric in $I$ and $J$ because the integrand is — which is the algebraic trace of the operator being self-adjoint, and the reason `COMPUTATIONAL.md` §2.6 can reach for Cholesky.

**The right-hand side is $\int f\varphi_I$**, a weighted integral. Three methods, three different right-hand sides: a point value, a cell integral, a weighted integral.

**Assembled.** The $5\times5$ grid gives 25 nodes; the 16 on the boundary are constrained, leaving the same nine unknowns as §2.1, in the same order. $\mathbf{K}$ restricted to those nine, with no scaling — in two dimensions $\int\nabla\varphi_I\cdot\nabla\varphi_J$ is dimensionless, so $\mathbf{K}$ carries no $h$. The load is evaluated with the three-point midside rule `code/disc.py` uses, exact for quadratics and in error by about 0.3% on this source:

```
 4 -1  0 | -1  0  0 |  0  0  0     T0     F0     0.585547
-1  4 -1 |  0 -1  0 |  0  0  0     T1     F1     0.785504
 0 -1  4 |  0  0 -1 |  0  0  0     T2     F2     0.525323
---------+----------+---------     --     --     --------
-1  0  0 |  4 -1  0 | -1  0  0     T3     F3     0.785504
 0 -1  0 | -1  4 -1 |  0 -1  0  x  T4  =  F4  =  1.110870
 0  0 -1 |  0 -1  4 |  0  0 -1     T5     F5     0.785504
---------+----------+---------     --     --     --------
 0  0  0 | -1  0  0 |  4 -1  0     T6     F6     0.525323
 0  0  0 |  0 -1  0 | -1  4 -1     T7     F7     0.785504
 0  0  0 |  0  0 -1 |  0 -1  4     T8     F8     0.585547
```

**This is §2.1's matrix, entry for entry.** §3 measures how far that goes.

**Where that last column comes from, and it does not factor.** The nodes are those of §2.1, indexed the same way — node $(i,j)$ with $i,j \in \{1,2,3\}$ is component

$$I = 3(j-1) + (i-1)$$

of $\mathbf{F}$, which is why $\mathbf{K}$ and $h^2\mathbf{A}_{\text{FD}}$ can be compared entry by entry at all. Each triangle has area $A_e = h^2/2$, and the midside rule gives a vertex $A_e/6$ times the sum of $f$ at the midpoints of the two edges meeting there. Each of the six edges at node $I$ is shared by two of the six triangles in the patch, so the two contributions add and every edge midpoint carries the same weight $A_e/3 = h^2/6$:

$$F_I = \frac{h^2}{6}\sum_{\text{6 edge midpoints}} f$$

The midpoints sit at offsets $(\pm\tfrac{h}{2},0)$, $(0,\pm\tfrac{h}{2})$, $(+\tfrac{h}{2},+\tfrac{h}{2})$ and $(-\tfrac{h}{2},-\tfrac{h}{2})$ from the node. Written out at node $I$, at position $(x,y)$, with $\tfrac{h^2}{6}\cdot 2\pi^2 = \tfrac{\pi^2}{48}$ at $h = \tfrac14$, so every argument lands on an eighth:

$$F_I = \frac{\pi^2}{48}\Big[\sin\pi(x{+}\tfrac{h}{2})\sin\pi y + \sin\pi(x{+}\tfrac{h}{2})\sin\pi(y{+}\tfrac{h}{2}) + \sin\pi x\sin\pi(y{+}\tfrac{h}{2})$$

$$+\; \sin\pi(x{-}\tfrac{h}{2})\sin\pi y + \sin\pi(x{-}\tfrac{h}{2})\sin\pi(y{-}\tfrac{h}{2}) + \sin\pi x\sin\pi(y{-}\tfrac{h}{2})\Big]$$

**In the index this is entirely a formula in eighths.** The node sits at $(2i, 2j)$ in units of $h/2 = \tfrac18$, and the six midpoints are one eighth away in the six directions, so with $s_m = \sin\frac{m\pi}{8}$

$$F_{ij} = \frac{\pi^2}{48}\Big[s_{2i+1}s_{2j} + s_{2i+1}s_{2j+1} + s_{2i}s_{2j+1} + s_{2i-1}s_{2j} + s_{2i-1}s_{2j-1} + s_{2i}s_{2j-1}\Big]$$

for $i, j \in \{1,2,3\}$. Unlike §2.1's $f_{ij}$ and §2.2's $b_{ij}$, **this does not factor into an $i$ part times a $j$ part** — six products, no common factor — which is the algebraic form of the same fact §2.3 keeps running into.

**Read the six offsets and the asymmetry is there before any arithmetic.** They include $(+\tfrac{h}{2},+\tfrac{h}{2})$ and $(-\tfrac{h}{2},-\tfrac{h}{2})$ and neither $(+\tfrac{h}{2},-\tfrac{h}{2})$ nor $(-\tfrac{h}{2},+\tfrac{h}{2})$, because the mesh diagonal runs SW–NE and there are no anti-diagonal edges to contribute. The stencil is symmetric under $x\leftrightarrow y$ and not under $x\leftrightarrow 1-x$. So although $f$ takes the identical value at $(\tfrac14,\tfrac14)$ and $(\tfrac34,\tfrac14)$, and $\mathbf{K}$ is perfectly symmetric, $F_0$ and $F_2$ differ by a tenth of their own size. §3.3 is that discrepancy.

This expression is the quadrature, not $\int f\varphi_I$: it reproduces the assembled load to $1.1\times10^{-16}$ and differs from the exact integral by about 0.3%. The asymmetry is a property of the patch, not of the rule, and survives an exact integration.

**That substitution has a name: a variational crime.** The term covers any method that approximates the bilinear or linear form rather than evaluating it, so what §2.3 assembles is not quite the Galerkin system for this problem — it is the Galerkin system for a nearby one. It is an entirely conventional crime, committed by essentially every production code, and it is bounded here at 0.3%. Two others belong to the same category and are worth naming together: **lumping the mass matrix**, which §3.1 reaches in any transient problem, and approximating a **curved boundary by straight element edges**. Neither arises here — §0's domain is a square and a steady problem has no mass matrix — but they are the same kind of departure, and the theory that bounds their effect (Strang's lemmas) is the theory that bounds this one.

**Everything else in §2.3 is the default on both axes** of the table above: P1 Lagrange elements on a uniform right-triangle mesh, continuous across element boundaries, conforming, with Bubnov–Galerkin testing. That is deliberate. §3's equality only means something because nothing clever was done — had the finite-element side used a stabilised formulation or an exotic element, "the matrices are identical" would be a statement about that choice rather than about the methods.

## 2.4 The three systems side by side

| | unknowns | matrix | right-hand side | value at the centre, $\div h^2$ |
|---|---|---|---|---|
| **Finite difference** | 9, point values | $9\times9$, diagonal 4 throughout | $f(x_i,y_j)$ | 19.739209 |
| **Finite volume** | 16, cell averages | $16\times16$, diagonal 4, 5 or 6 | $\int_V f$ | 16.000000 |
| **Finite element** | 9, basis coefficients | $9\times9$, diagonal 4 throughout | $\int f\varphi_I$ | 17.773924 |

with $f(\tfrac12,\tfrac12) = 19.739209$ for comparison.

**The matrices are two-thirds of a coincidence and the right-hand sides are not close.** Finite difference and finite element produced the same nine equations; finite volume produced sixteen, of which four match and twelve carry the wall term. The three sources disagree by up to 19% at the same physical location — and since two of the three methods are solving with the *same* matrix, on this mesh the right-hand side is the entire difference between them. §3 measures the matrices properly and §3.1 takes the right-hand sides apart.

---

# 3. First result: on a uniform mesh they are the same operator

Assemble all three on the unit square with $k$ constant and a uniform mesh, and extract one interior row of each, scaled so that each approximates $-\nabla^2$:

![Interior stencils](figs/fig2_stencils.png)

Measured directly from the assembled matrices at $n=8$:

| Comparison | max difference over the interior row |
|---|---|
| finite difference against finite volume | **0.0** |
| finite difference against FEM P1 (right triangles) | **0.0** |
| finite difference against FEM Q1 (bilinear quads) | 1.333 |

**Not "similar" — identical, to the last bit.** On a uniform mesh with constant coefficients, the 5-point finite-difference Laplacian, the cell-centred finite-volume Laplacian, and the P1 finite-element stiffness matrix on the diagonal triangulation are the same matrix. Three derivations from three different starting points converge on one operator.

**How far "the same matrix" extends is not the same for the two pairs**, and §2's assembled systems show where the line falls. Finite difference against P1 is a whole-matrix statement: on the $4\times4$ mesh $\max|h^2\mathbf{A}_{\text{FD}} - \mathbf{K}_{\text{P1}}| = 0.0$ over all 81 entries, boundary-adjacent rows included. Finite volume against either is a strictly interior statement: its matrix is a different size, and only the four cells of the middle $2\times2$ block carry the 5-point row — the other twelve carry the wall term of §2.2. The row extracted above is one of the four. **The equivalence is a statement about the interior operator, not about the discrete problem**, and the boundary treatment is where the three methods have already parted company before any test in §4.

**Bilinear quadrilaterals do not**, and this is worth seeing because it shows the equivalence is a coincidence of the *element*, not a law. Q1 gives the 9-point stencil

$$\frac{1}{3}\begin{pmatrix} -1 & -1 & -1 \\ -1 & 8 & -1 \\ -1 & -1 & -1\end{pmatrix}$$

which is also a second-order Laplacian, with a wider support, more nonzeros per row, and a different constant in its error.

## 3.1 What still differs, even here

Three things, and each matters somewhere.

**1. The right-hand side — three weights on one construction.** Point value, cell integral, weighted integral look like three unrelated recipes. They are one:

$$b_I = \int_\Omega f\,w_I, \qquad w_I = \begin{cases} h^2\,\delta(\mathbf{x} - \mathbf{x}_I) & \text{finite difference}\\ \mathbb{1}_{V_I} & \text{finite volume}\\ \varphi_I & \text{finite element}\end{cases}$$

A spike, an indicator, a tent. The finite-difference case does not look like an integral only because the delta collapses it to a sample — which is precisely why $f$ must be defined pointwise there and need only be integrable in the other two. That is the practical content of "not identical when $f$ has a peak or a discontinuity": the finite-volume and finite-element right-hand sides cannot be changed by moving a spike between sample points, and the finite-difference one can.

**All three weights carry the same total mass $h^2$.** Immediate for the first two. For the third, a P1 hat restricted to one triangle is a barycentric coordinate and integrates to $A_e/3$, so over the six triangles of its patch $\int\varphi_I = 6\cdot\tfrac{h^2/2}{3} = h^2$. So each right-hand side divided by $h^2$ is a weighted average of $f$ with weights summing to one, which is what makes §2.4's three centre values comparable: 19.739209, 16.000000 and 17.773924 against $f(\tfrac12,\tfrac12) = 19.739209$. **They order by how spread out the weight is** — the spike returns the peak exactly, the cell average pulls furthest down, the tent lands between because it is still concentrated at the node.

**The weights are not arbitrary: each method tests the equation with the same kind of object it uses to represent the solution.** Point values tested at points, cell averages tested over cells, basis coefficients tested against the basis — collocation, Petrov–Galerkin with piecewise-constant test functions, and Galerkin proper. That is §1's table read along its other axis.

**On the $4\times4$ mesh that $O(h^2)$ is 5%, and it is the whole error.** The manufactured solution is an exact discrete eigenvector of the 5-point operator, with eigenvalue $(4-4\cos\pi h)/h^2$, so each method's answer is its source divided by that number and the comparison is clean:

| | what the solution does | $n=4$ | $n=8$ | $n=16$ | $n=32$ |
|---|---|---|---|---|---|
| **Finite difference** | exact point value $\times\;2\pi^2h^2/(4-4\cos\pi h)$ | 1.053029288 | 1.012950747 | 1.003218964 | 1.000803578 |
| **Finite volume** | exact point value, $L^\infty$ discrepancy | 3.3e-16 | 4.4e-16 | 8.9e-16 | 5.3e-15 |

**The finite-difference row is its $O(h^2)$ error in closed form**, predicted and measured agreeing to nine digits. **The finite-volume row is zero** — those entries are roundoff, growing with problem size as roundoff does, and they will differ in the last digit from machine to machine. It is not luck: integrating the source exactly over a cell contributes a factor $\left[\sin(\pi h/2)/(\pi h/2)\right]^2$, and $2\pi^2 h^2\left[\sin(\pi h/2)/(\pi h/2)\right]^2 = 8\sin^2(\pi h/2) = 4 - 4\cos\pi h$ identically. The cell-integrated source and the discrete eigenvalue carry the same factor and it cancels.

**Do not over-read the second row.** It says that for a source which is a single eigenmode integrated exactly, the cell-centred scheme is pointwise exact at any resolution — not that finite volume is generally more accurate. Measured against what it actually approximates, the cell average, the same solution is in error by $4.3\times10^{-2}$ at $n=4$; the entire gap is the $O(h^2)$ difference between a point value and a cell average, which is §3.1 item 3 arriving as a number.

**2. The mass matrix — which is where they part company in any transient problem.** For $\partial_t T = \nabla^2 T$, the semi-discrete system is $\mathbf{M}\dot{\mathbf{T}} + \mathbf{K}\mathbf{T} = \mathbf{F}$. Finite difference and cell-centred finite volume give $\mathbf{M} = h^2\mathbf{I}$, diagonal. Galerkin finite elements give the **consistent mass matrix** $M_{IJ} = \int\varphi_I\varphi_J$, which is not diagonal — so an "explicit" FEM time step still requires a solve unless the mass is **lumped** — the second of §2.3's variational crimes, and the one with a practical payoff rather than a cost. `COMPUTATIONAL.md` §2.2 records that lumping is what makes explicit structural dynamics matrix-free; this is where that requirement comes from. **Equal stiffness matrices do not imply equal transient behaviour.**

**3. What the unknown means.** A point value, a cell average, a basis coefficient. To measure an error one must first decide what to compare against, and the honest comparison for a cell-centred method is against the exact **cell average**, not the exact point value at the centre — those differ by $O(h^2)$, the same order as the error being measured.

## 3.2 Convergence

With the error defined per method as above:

![Convergence](figs/fig3_convergence.png)

$L^2$ errors on the manufactured solution, from `code/disc.py`:

| $n$ | $h$ | FD | FVM | FEM P1 | FEM Q1 |
|---|---|---|---|---|---|
| 8 | 0.1250 | 6.4754e-03 | 6.3926e-03 | 6.5378e-03 | 6.4749e-03 |
| 16 | 0.0625 | 1.6095e-03 | 1.6043e-03 | 1.6340e-03 | 1.6095e-03 |
| 32 | 0.0312 | 4.0179e-04 | 4.0147e-04 | 4.0846e-04 | 4.0179e-04 |
| 64 | 0.0156 | 1.0041e-04 | 1.0039e-04 | 1.0211e-04 | 1.0041e-04 |
| 128 | 0.0078 | 2.5100e-05 | 2.5099e-05 | 2.5528e-05 | 2.5100e-05 |
| **observed order** | | **2.00** | **2.00** | **2.00** | **2.00** |

**All four are second order, and the constants agree to within 2%.** On a smooth problem, on a uniform mesh, with constant coefficients, **the choice of method does not matter**. Anyone claiming a large accuracy advantage for one of these three on this class of problem is comparing implementations, not methods.

Note also that Q1 nodal errors match FD to five digits (6.4749e-03 against 6.4754e-03) despite a completely different stencil. That is **nodal superconvergence**: the finite-element solution is more accurate *at the nodes* than its own global order would suggest.

## 3.3 The triangulation's diagonal, which enters through $\mathbf{F}$ and not $\mathbf{K}$

§2.3's load vector is not symmetric and the problem is. $f$ is invariant under reflection in either diagonal of the square; the node set is too; $\mathbf{K}$ is §2.1's matrix exactly. Yet $F$ at $(\tfrac14,\tfrac14)$ is 0.585547 and at $(\tfrac34,\tfrac14)$ is 0.525323.

**The cause is the one asymmetric object in the construction.** Splitting every square along its SW–NE diagonal makes each interior node's support a hexagon whose long axis lies along that diagonal. The hexagon is symmetric under reflection in the mesh diagonal and not under reflection in the anti-diagonal, so two nodes that the *problem* cannot distinguish sit differently inside their own basis functions.

**It does not reach $\mathbf{K}$.** Splitting the squares in alternating orientation, chessboard fashion, leaves the stiffness matrix bit-identical — $\max|\mathbf{K}_{\text{fixed}} - \mathbf{K}_{\text{alt}}| = 0.0$ at $n=4$ and again at $n=8$. Both orientations are right-triangle meshes on a uniform grid and both assemble to the 5-point operator. **The diagonal is invisible in the operator and visible in the source.**

| | corner-to-corner spread in $F$ | same in the solution | $L^\infty$ error |
|---|---|---|---|
| **one diagonal throughout** | $1.03\times10^{-1}$ | $3.13\times10^{-2}$ | 5.1813e-02 |
| **alternating diagonals** | $1.5\times10^{-16}$ | $1.1\times10^{-16}$ | 7.2986e-02 |

**Alternating restores the symmetry exactly and makes the answer worse**, by 41% in $L^\infty$. The two effects are unrelated: symmetry of the load is a statement about how the mesh treats equivalent points, accuracy is a statement about how well $\sum F_I\varphi_I$ represents $f$, and the chessboard mesh is better at the first and worse at the second. A symmetry defect is not an error bound, and repairing one is not the same as reducing the other.

**What this costs the reader of §3.** Finite difference and finite element solve with the same matrix on this mesh and still return different answers, one of them symmetric and one not. Everything separating them is in $\mathbf{F}$.

**So the interesting question is not which is more accurate. It is which assumptions each one needs**, and what happens when they fail.

---

# 4. Where they diverge — four tests

Each test removes one assumption from §3 and nothing else.

## 4.1 Geometry — remove the structured grid

| | On an unstructured or non-rectangular mesh |
|---|---|
| **Finite difference** | **The method largely does not apply.** It needs a logically rectangular grid, so the options are a body-fitted curvilinear mapping (which introduces metric terms and coordinate singularities), an immersed or cut-cell boundary treatment, or overset grids. Each is a substantial complication and each degrades accuracy at the boundary |
| **Finite volume** | **Applies unchanged.** Cells may be arbitrary polygons or polyhedra. The only requirement is that each face has an area and a normal. Skewness and non-orthogonality degrade the *flux approximation* and require correction terms, but the conservation property is untouched |
| **Finite element** | **Applies unchanged, and is the most comfortable.** Elements are mapped from a reference element by an isoparametric map, so curved boundaries are represented by curved elements; the entire apparatus of mesh generation was built for it |

**This alone decides the matter for most engineering geometry**, and it is the reason finite differences survive mainly where the geometry is a box: direct numerical simulation of turbulence in periodic domains, seismic wave propagation, and academic test problems. `COMPUTATIONAL.md` §4.1 lists that usage.

## 4.2 Boundary conditions — remove pure Dirichlet

| Condition | Finite difference | Finite volume | Finite element |
|---|---|---|---|
| **Dirichlet**, $T = g$ | direct: eliminate the node, or set the row to identity | face value is $g$, at a half-cell distance $h/2$ — a slight asymmetry in the stencil | **essential**: constrain the coefficient, and lift the known value to the right-hand side |
| **Neumann**, $-k\partial_nT = q$ | **awkward**: needs a ghost node and a mirrored expansion, or a one-sided difference that loses an order | **the most natural of the three**: the face flux is simply set to $q$ and no gradient is ever needed | **natural**: it *is* the boundary term of the weak form, added to $\mathbf{F}$, requiring nothing else |
| **Robin**, $-k\partial_nT = h_c(T-T_\infty)$ | ghost node plus an algebraic elimination | face flux written in terms of the cell value; contributes to both the diagonal and the right-hand side | boundary term contributes to both $\mathbf{K}$ and $\mathbf{F}$ |

**The pattern is worth naming.** In finite volume and finite element, a flux condition is *the thing the method already works with* — a face flux in one case, a boundary integral in the other. In finite difference it is a foreign object, because the method is built on point values of the solution and a flux is a derivative of it.

Conversely, Dirichlet is cleanest in finite difference and finite element (both carry nodal values, and a node can sit exactly on the boundary) and slightly awkward in cell-centred finite volume, where nothing sits on the boundary and the half-cell distance breaks the stencil's symmetry.

## 4.3 Discontinuous coefficients — remove constant $k$

Two materials, $k_1 = 1$ on the left half and $k_2 = 1000$ on the right, with $T=0$ at $x=0$, $T=1$ at $x=1$, and insulated top and bottom. The exact solution is piecewise linear with a **uniform flux**

$$q = \frac{\Delta T}{L_1/k_1 + L_2/k_2} = \frac{1}{0.5/1 + 0.5/1000} = 1.998001998\ldots$$

![Two-material slab](figs/fig5_interface.png)

**The question is what conductivity to use on the face between two cells of different material.** Two candidates, and the difference is not cosmetic. From `code/disc2.py`, with the interface landing exactly on a cell face:

| $n$ | harmonic mean face $k$ | arithmetic mean face $k$ |
|---|---|---|
| 8 | rel. flux error **8.9e-16** | rel. flux error **14.2%** |
| 16 | **1.3e-15** | **6.6%** |
| 32 | **4.4e-16** | **3.2%** |

**The harmonic mean is exact at every resolution**, including a mesh of eight cells. The arithmetic mean is first-order convergent and badly wrong on a coarse mesh — and note that it converges, so a mesh-refinement study would eventually hide the error while never revealing that a better choice existed.

**Why the harmonic mean.** The face sits between a half-cell of $k_P$ and a half-cell of $k_E$, and heat must pass through both. Thermal resistances in series **add**:

$$R = \frac{h/2}{k_P} + \frac{h/2}{k_E} \quad\Longrightarrow\quad q = \frac{T_P - T_E}{R} = \frac{T_P - T_E}{h}\cdot\underbrace{\frac{2k_Pk_E}{k_P + k_E}}_{k_{\text{face}}}$$

which is the harmonic mean. **The arithmetic mean averages the wrong quantity** — it averages conductances where the physics adds resistances, and with $k_2/k_1 = 1000$ it is dominated by the conductor and effectively ignores the insulator that is actually controlling the flux.

**Finite element on an interface-aligned mesh is also exact**, and for the same reason expressed differently: the stiffness integral $\int k\nabla\varphi_I\cdot\nabla\varphi_J$ is evaluated element by element with each element's own $k$, which reproduces the series-resistance arithmetic automatically. Measured, via the nodal reaction at the $x=1$ boundary:

| $n$ | 8 | 16 | 32 |
|---|---|---|---|
| FEM P1 flux, relative error | 2.8e-13 | 2.0e-12 | 2.4e-12 |

**Naive finite difference has no good answer here at all.** The interface passes through a node, at which $k$ is undefined; sampling $k$ at nodes and differencing $k\,\partial_x T$ gives the arithmetic-mean answer or worse. The standard repair is to write the scheme in flux form with a harmonic face conductivity — **which is to say, to stop doing finite differences and start doing finite volumes.**

**The general lesson generalises past this problem.** A discontinuity in a coefficient is not a small perturbation to be resolved by refinement. It is a place where the *form* in which the method is written decides whether the answer is right, and the integral form is the one that gets it right.

## 4.4 Advection — remove self-adjointness

Add a uniform velocity: $U\,\partial_x T = \alpha\,\partial_x^2 T$ with $T(0)=0$, $T(1)=1$, at global Péclet number $UL/\alpha = 50$. The exact solution is a boundary layer at $x=1$:

$$T(x) = \frac{e^{\mathrm{Pe}\,x} - 1}{e^{\mathrm{Pe}} - 1}$$

The operator is no longer self-adjoint, and every guarantee of §3 is void.

![Advection at high cell Peclet number](figs/fig4_advection.png)

With 10 cells, so **cell Péclet number** $\mathrm{Pe}_h = Uh/\alpha = 5$:

| Scheme | min $T$ | behaviour |
|---|---|---|
| FVM, central face value | $-0.1224$ | **oscillates** |
| FVM, upwind face value | $0.0$ | monotone, but heavily over-diffused |
| FEM, Galerkin P1 | $-0.8577$ | **oscillates violently** |
| FEM, SUPG | $0.0$ | monotone, and nodally near-exact |

**Galerkin finite elements and central finite volumes fail in the same way**, and they fail together because they are the same scheme in this respect: linear Galerkin on pure advection reduces to a central difference. The optimality result that makes Galerkin the best possible choice for §3's elliptic problem (`COMPUTATIONAL.md` §2.1, Céa's lemma) requires the operator to be self-adjoint. Advection is not, and the guarantee evaporates precisely where it would be most useful.

### The threshold is exactly 2, and it can be derived

For the central scheme the coefficient multiplying the downstream neighbour is

$$a_E = \underbrace{\frac{\alpha}{h}}_{\text{diffusion}} - \underbrace{\frac{U}{2}}_{\text{advection}}, \qquad a_E < 0 \iff \frac{Uh}{\alpha} > 2$$

A positive off-diagonal (after moving to the standard form) destroys the **M-matrix** property, and with it the discrete maximum principle that guarantees the solution stays between its boundary values. Measured, by bisection over $n$:

| $n$ | $\mathrm{Pe}_h$ | $\min T$ | |
|---|---|---|---|
| 24 | 2.083 | $-9.996\times10^{-3}$ | oscillates |
| **25** | **2.000** | $0.0$ | **monotone** |
| 26 | 1.923 | $5\times10^{-44}$ | monotone |

**The predicted threshold and the measured one agree exactly.** This is the sharpest single result in the document: a stability limit derived from the sign of one matrix entry, confirmed to the resolution of the bisection.

### What upwinding costs

Upwind takes the face value from the upstream cell, giving $a_E = \alpha/h$ — never negative, monotone at any $\mathrm{Pe}_h$. The price is visible in the figure and is quantifiable exactly:

$$a_E^{\text{upwind}} - a_E^{\text{central}} = \frac{U}{2} \quad\Longleftrightarrow\quad \alpha_{\text{effective}} = \alpha + \frac{Uh}{2}$$

**First-order upwinding is central differencing plus an artificial diffusivity of $Uh/2$.** At $\mathrm{Pe}_h = 5$ that is 2.5 times the physical diffusivity — the scheme contributes more diffusion than the physics does. The measured error at the last cell centre is 2.04e-01 at $n=10$, 1.58e-01 at $n=20$, 6.0e-02 at $n=50$: converging, slowly, and by brute force.

**This is exactly the situation `COMPUTATIONAL.md` §4.2 describes**, and Godunov's theorem is why there is no linear escape from it: a linear monotone scheme is at most first order. Flux limiters, MUSCL and WENO are the nonlinear escape.

### SUPG is the finite-element answer

Streamline-upwind Petrov–Galerkin modifies the *test* function rather than the flux, weighting it along the streamline:

$$v \;\longrightarrow\; v + \tau\,(\mathbf{u}\cdot\nabla v), \qquad \tau = \frac{h}{2U}\left(\coth\frac{\mathrm{Pe}_h}{2} - \frac{2}{\mathrm{Pe}_h}\right)$$

Test space no longer equals trial space, so it is Petrov–Galerkin rather than Galerkin. In one dimension with this $\tau$ the scheme is **nodally exact at any $\mathrm{Pe}_h$**, which the figure shows. That exactness is a one-dimensional accident and does not survive in 2D — but the stabilisation does, and it is why finite elements are usable in fluid dynamics at all.

**The deeper point.** Upwind finite volume and SUPG finite element are the same idea reached from different directions: **respect the direction information travels**. One does it by choosing which cell supplies the face value, the other by tilting the test function upstream. Both are the discrete form of the domain-of-dependence argument in `COMPUTATIONAL.md` §0.3.

---

# 5. Conservation, precisely

"Finite volume is conservative" is repeated so often that it has stopped carrying information. The precise statements are worth separating, because two of the three are true and one of them is commonly overstated.

## 5.1 What is exactly true of finite volume

Summing the discrete balance over any set of cells makes every interior face appear twice, with opposite signs, and cancel. What remains is the flux through the outer boundary of that set, equated to the source integrated over it. **This holds on any mesh, at any resolution, to machine precision, whether or not the solution is accurate.**

Measured on the two-material slab of §4.3 at $n=32$, where the exact answer is a flux that is the same through every vertical plane:

| | |
|---|---|
| flux through each of the 31 interior planes | min $-1.998001998008$, max $-1.998001997987$ |
| **spread** | $2.1\times10^{-11}$ |
| exact | $1.998001998\ldots$ |

**The spread is roundoff.** The conservation statement is not converging to true — it is true, and the residual is the condition number of the solve times machine epsilon.

**Why this matters and where it does not.** For a shock, it decides whether the shock travels at the right speed (`COMPUTATIONAL.md` §4.1, the Lax–Wendroff theorem). For a long transient, it decides whether mass slowly disappears. For a steady linear elliptic problem solved to convergence, it buys much less than its reputation suggests, and the equivalence of §3 is the proof — on a uniform mesh the finite-difference scheme *is* the finite-volume scheme and is therefore equally conservative.

## 5.2 Where finite difference stops being conservative

The equivalence of §3 required a **uniform** mesh. On a stretched or non-uniform grid, a finite-difference scheme written by Taylor expansion about each node does not generally telescope: the coefficient with which a node's equation references its neighbour is not the negative of the coefficient with which the neighbour references it, so summing over a patch leaves interior residue. There is no flux to cancel because the method never constructed one.

**The repair is to write the scheme in flux form** — compute a face quantity, then difference it — at which point one is doing finite volumes and may as well say so. This is the same conclusion §4.3 reached from the interface problem, by a different route.

## 5.3 What is true of finite element

More subtle, and usually stated too weakly or too strongly.

**Too strong:** "the finite-element solution is locally conservative." Taking $-k\nabla T_h$ from the computed solution and integrating it around a patch does **not** generally give the source integral exactly.

**Too weak:** "finite elements are not conservative." They are, in the sense that matters, but the conservative flux has to be extracted rather than read off. The **consistent nodal flux** — the reaction $\mathbf{r} = \mathbf{K}\mathbf{T} - \mathbf{F}$ at constrained nodes — satisfies the balance exactly, and this follows from the weak form: taking $v \equiv 1$ on a patch makes the weak statement *be* the conservation statement on that patch.

**And it is measurable.** The §4.3 finite-element flux was computed exactly this way, as the summed reaction along the $x=1$ boundary, and it reproduced the exact flux to $2.8\times10^{-13}$ on an 8-element mesh with a 1000:1 conductivity jump. **The information is there; it is in the reactions, not in the gradient of the solution.** This result is due to Hughes and co-workers and it is the correct answer to the question. **[unverified]** — the attribution is recalled, not checked.

---

# 6. What each method costs and gives

| | Finite difference | Finite volume | Finite element |
|---|---|---|---|
| **Implementation effort, simplest case** | **lowest** — a stencil loop | moderate — face loops, geometry | **highest** — shape functions, quadrature, assembly, mapping |
| **Arbitrary geometry** | poor | good | **excellent** |
| **Local conservation** | on uniform grids only | **exact, always** | exact in the consistent-flux sense |
| **Natural handling of flux BCs** | poor | **excellent** | **excellent** |
| **Discontinuous coefficients** | poor | **good** (harmonic face values) | **good** (interface-aligned mesh) |
| **High order** | **easy** — widen the stencil | hard — reconstruction over cell averages, and hardest on unstructured meshes | **easy** — raise the polynomial degree, $p$-refinement |
| **Advection without extra work** | no | no (needs upwinding) | no (needs SUPG) |
| **Rigorous error theory** | Taylor / maximum principle | weaker in general | **strongest** — Céa, best approximation, a posteriori estimators |
| **Adaptivity** | poor | good | **excellent** — $h$-, $p$-, and $hp$-refinement, with computable error indicators |
| **Matrix** | banded, structured | sparse, one row per cell | sparse, symmetric when the operator is |
| **Where it dominates** | DNS in boxes, seismic, finance | industrial CFD, compressible flow | structural mechanics, electromagnetics, anything with geometry |

**Two of these rows deserve a sentence.**

**High order.** Finite differences and finite elements raise order by a mechanism they already have — a wider stencil, a higher polynomial. Finite volume must *reconstruct* a high-order profile from cell averages, which on an unstructured mesh requires a least-squares stencil that is expensive, and is why high-order unstructured CFD drifted toward discontinuous Galerkin instead (§8).

**Error theory.** The finite-element framework yields computable **a posteriori** error estimators — quantities assembled from the solution itself that bound the error — and those estimators are what drive automatic mesh adaptation. Nothing of comparable strength exists for the other two, and it is a large part of why finite elements dominate where a certified answer is wanted.

---

# 7. Choosing

| If the problem has | Reach for |
|---|---|
| a rectangular domain and a smooth solution | **finite difference** — simplest thing that works, and easiest to make high order |
| complex geometry and a stress or displacement field | **finite element** |
| shocks, or a long transient where mass must not drift | **finite volume** |
| discontinuous material properties | **finite volume** with harmonic face values, or **finite element** with an interface-aligned mesh |
| flux boundary conditions as the primary data | **finite volume** or **finite element** |
| advection dominance | either, **with the appropriate stabilisation** — this is not a differentiator |
| an energy or variational principle | **finite element** — it is the discretization of the principle itself |
| a requirement for a certified error bound and adaptive refinement | **finite element** |

**This is the substance behind `COMPUTATIONAL.md` §2.1 and §4.1.** Structural mechanics reaches for finite elements because it has geometry, a variational principle, and no shocks. Fluid dynamics reaches for finite volumes because it has conservation laws, discontinuous solutions, and a modest need for geometric fidelity in the interior. **Neither choice is about accuracy**, and §3 is the evidence: on the problem where both apply cleanly, they give the same answer to within 2%.

---

# 8. The fourth family, in passing

Named because each one dissolves part of the comparison above rather than sitting inside it.

**Discontinuous Galerkin** is the interesting case. It carries a polynomial basis within each cell (finite element) and couples cells only through numerical fluxes computed by a Riemann solver (finite volume). It is **locally conservative like finite volume and arbitrarily high order like finite element**, on unstructured meshes, which is precisely the combination §6 said was hard. The cost is many more degrees of freedom for a given mesh — the solution is duplicated at every interface — and a restrictive explicit time step that scales as $h/p^2$. **[unverified]** — the $p^2$ scaling is recalled from the literature and not checked.

**Spectral and spectral-element methods** replace piecewise polynomials with global or high-degree local expansions, converging exponentially for analytic solutions. This is why direct numerical simulation of turbulence uses them: at a fixed error, they need far fewer points, which is decisive when the point count scales as $\mathrm{Re}^{9/4}$.

**Meshfree and particle methods** — SPH, material point, radial basis functions — abandon the mesh. They are natural where the domain fragments or deforms beyond what a mesh survives, and they pay for it with difficulty imposing boundary conditions and enforcing conservation.

**Lattice Boltzmann** does not discretize the Navier–Stokes equations at all. It evolves a discrete kinetic equation whose macroscopic limit is Navier–Stokes, on a fixed lattice, with an entirely local update.

---

# 9. Runs

Provenance for every number quoted above, following the convention in `fusion`'s `CLAUDE.md`: a number that cannot be traced to the invocation that produced it is not a result.

| Script | Section | What it computes | Key output |
|---|---|---|---|
| `code/disc.py` | §2, §3 | **E0**: assembles all three at $n=4$ and prints the full systems, the three right-hand sides, the eigenvector amplification factor, and the fixed-against-alternating diagonal comparison. **E1**: interior stencils at $n=8$. **E2**: refinement study at $n = 8,16,32,64,128$ | $h^2\mathbf{A}_{\text{FD}} = \mathbf{K}_{\text{P1}}$ to 0.0 over the whole $9\times9$; FVM pointwise exact to 3e-16; FD $=$ FVM $=$ P1 interior row to 0.0; all four second order |
| `code/disc2.py` | §4.3, §4.4, §5.1 | two-material slab with harmonic and arithmetic face conductivity; FEM P1 reaction flux; advection–diffusion with central, upwind, Galerkin and SUPG; bisection on the monotonicity threshold | harmonic exact to 9e-16; arithmetic 14.2% at $n=8$; threshold at $\mathrm{Pe}_h = 2$ exactly |
| `code/figs.py` | all | the five figures in `figs/`, from the JSON dumped by the two scripts above | — |

```bash
cd code && python3 disc.py && python3 disc2.py && python3 figs.py
```

Dependencies are NumPy, SciPy and Matplotlib only, all already pinned in the repository's `pyproject.toml`; run with `uv run python disc.py` from `code/` to use the project environment. The numbers above were produced with NumPy 2.4.4, SciPy 1.17.1 and Matplotlib 3.10.8. Every problem here is small enough for a sparse direct factorization, so nothing depends on the solver.

---

# 10. Caveats

Stated while fresh, per `COMPUTATIONAL.md`'s own standard. None of these invalidates a result above; each limits what it may be quoted for.

**The §4.3 and §4.4 tests are one-dimensional in $x$.** Both problems are posed on the 2D square with insulated top and bottom, so the exact solution is $y$-invariant and the discrete problem reduces exactly to its $x$-line. That is a genuine reduction, not an approximation — but it means neither test exercises mesh skewness, cross-diffusion, or a flow oblique to the grid, which is where multidimensional upwinding actually becomes difficult. **A 2D test with the velocity at 45° to the mesh would show far worse behaviour from every scheme here**, and would be the natural next experiment.

**The FD/FVM/FEM-P1 equivalence of §3 is specific.** It requires a uniform mesh and constant $k$; a non-uniform mesh or a variable coefficient breaks it. It does *not* require a fixed diagonal orientation — §3.3 measures the alternating triangulation and finds the same stiffness matrix — but it does require that every element be a right triangle with legs on the grid. A criss-cross triangulation, splitting each square into four triangles about an added centre node, introduces a degree of freedom the finite-difference grid does not have and cannot give the same matrix. **[unverified]** — that last case was not assembled here. The result should be read as "these methods coincide in the simplest case", not as a general identity.

**Second-order convergence was measured on one smooth manufactured solution.** Solutions with corner singularities — the re-entrant corner of an L-shaped domain is the standard example — converge more slowly for every method, and the *relative* ranking on such problems has not been tested here.

**The SUPG $\tau$ is the one-dimensional nodally-exact formula.** Its exactness is a 1D accident. Multidimensional $\tau$ definitions are numerous, are not equivalent, and the choice matters; nothing here bears on that.

**The advection tests use a Dirichlet outflow.** The face there takes the upwind value convectively and the wall value diffusively. Other consistent treatments exist and would shift the last cell's value; the oscillation threshold of §4.4 is an interior property and is unaffected.

**Nothing here is a performance measurement.** All timings would be dominated by the fact that these are Python assembly loops on very small problems. The cost rows of §6 are qualitative and are not backed by measurement in this document.

---

# 11. Unverified claims

1. **§5.3** — attribution of the consistent-nodal-flux result to Hughes and co-workers. The *result* is demonstrated numerically here to 2.8e-13; only the attribution is unverified.
2. **§8** — the $h/p^2$ explicit time-step restriction for discontinuous Galerkin. Recalled, not checked.
3. **§4.1** — the characterisation of where finite differences remain dominant in practice (DNS, seismic, finance). General impression, not surveyed.
4. **§10** — that a four-triangle criss-cross triangulation cannot reproduce the 5-point matrix. Argued from the degree-of-freedom count, not assembled. **What would settle it:** add a `criss` option alongside `p1_solve`'s `diagonal` argument and compare.
5. **§2.3** — that P2 gives third-order $L^2$ accuracy. Standard textbook result, stated for completeness of the naming; no P2 element is assembled here. **What would settle it:** add a `p2_solve` to `code/disc.py` and extend §3.2's table.
6. **§2.3** — that $Q$ was chosen as the letter following $P$. The reading is consistent with the notation's history and with $Q$ meaning neither quadratic nor, in 3D, quadrilateral, but no source was consulted. **What would settle it:** Ciarlet's *The Finite Element Method for Elliptic Problems*, where the notation is set.
7. **§6** — the row on a posteriori error estimation being materially stronger for finite elements than for the other two. This is the standard position and is believed correct, but no comparison of the finite-volume a posteriori literature was made.

**What would settle these:** one primary source each for 1, 2, 3, 6 and 7; items 4 and 5 are computations and are noted above. None is load-bearing — every quantitative claim in this document rests on §9's scripts rather than on recall.

---

# 12. Sources

**General and comparative**

- Strikwerda J.C., *Finite Difference Schemes and Partial Differential Equations* — the Taylor/stability apparatus of §2.1
- LeVeque R.J., *Finite Difference Methods for Ordinary and Partial Differential Equations*
- Patankar S.V., *Numerical Heat Transfer and Fluid Flow* — the origin of the harmonic-mean face conductivity of §4.3 and of the cell-Péclet analysis of §4.4
- Versteeg H.K., Malalasekera W., *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* — the clearest elementary treatment of §2.2 and §4.4
- Ferziger J.H., Perić M., Street R.L., *Computational Methods for Fluid Dynamics*
- Hughes T.J.R., *The Finite Element Method* — §2.3, the consistent flux of §5.3, and SUPG
- Brenner S., Scott L.R., *The Mathematical Theory of Finite Element Methods* — Céa's lemma and the approximation theory behind §3.2
- Eymard R., Gallouët T., Herbin R., "Finite Volume Methods," *Handbook of Numerical Analysis* VII — the rigorous finite-volume analysis that §6's "weaker in general" row refers to

**Specific results**

- Godunov S.K. (1959) — the barrier theorem invoked in §4.4
- Brooks A.N., Hughes T.J.R., "Streamline upwind/Petrov-Galerkin formulations…," *CMAME* **32**, 199 (1982) — SUPG and the $\tau$ of §4.4
- Cockburn B., Shu C.-W. — the Runge–Kutta discontinuous Galerkin series behind §8
- Chacón L. et al., *A fully implicit, asymptotic-preserving, semi-Lagrangian algorithm for the time dependent anisotropic heat transport equation*, arXiv:2404.08771 — the $10^7$–$10^{10}$ anisotropy ratio of §0.1, and the observation that the operator's condition number scales with it
- Green D.L. et al., *Mesh refinement for anisotropic diffusion in magnetized plasmas*, arXiv:2210.16442 — edge-region anisotropy estimates and the boundary layers it produces

**Companion documents**

- `COMPUTATIONAL.md` §0.3 — hyperbolic, parabolic and elliptic, and the domain-of-dependence argument §4.4 appeals to
- `COMPUTATIONAL.md` §2.1 — why structural mechanics reaches for finite elements; §7 above is the evidence
- `COMPUTATIONAL.md` §4.1–§4.2 — why fluid dynamics reaches for finite volumes, and Godunov's theorem
- `COMPUTATIONAL.md` §2.2 — lumped against consistent mass, whose origin is §3.1 above
