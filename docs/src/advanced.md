# Advanced

The purpose of this article is to provide additional information about magnon-holon representation.
We will cover here details of the transformation to magnon-holon basis, including:
- sublattice-wise transformation of standard operators,
- discussion of magnon-magnon interaction scaling on commutation of the Hamiltonian with translation operator
- definition and numerical representation of momentum states and representative states in mangon-holon basis,
- derivations for operators action on magnon-holon representative basis states

---

Let us start by introducing bosonic magnon ``\hat{a}`` and fermionic holon ``\hat{h}`` operators.
What we call the magnon-holon transformation is the procedure of rewriting operators in terms of ``\hat{a}`` and ``\hat{h}`` operators.
Accordingly, we say that operator ``\hat{O}`` expressed in terms of ``\hat{a}`` and ``\hat{h}`` is in magnon-holon representation (or in magnon-holon basis).
The magnon-holon transformation is a combination of Holstein-Primakoff transformation and Slave-fermion transformation (see e.g. [G. Martinez and P. Horsch, PRB 44 (1991)](https://doi.org/10.1103/PhysRevB.44.317)).
Considering magnons as excitations around perfect antiferromagnetic spin background, below set of equations define the magnon-holon representation 
of the standard operators in the [tJMagnonHolon](@ref). We assume the system is one-dimensional with ``L = 2N`` sites, 0-based indexed, and ``N`` - positive integer.

On sublattice A (sites with even indices ``\{0, 2, 4, ..., L - 2\}``):
```math
\begin{aligned}
	\hat{\tilde{c}}_{i\uparrow}^\dag &= \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\uparrow} &= \hat{h}_i^\dag \hat{P}_i, \\
	\hat{\tilde{c}}_{i\downarrow}^\dag &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\downarrow} &= \hat{h}_i^\dag \hat{P}_i \hat{a}_i, \\
    \hat{S}_i^+ &= \hat{h}_i \hat{h}_i^\dag \hat{P}_i \hat{a}_i, &\quad \hat{S}_i^z &= \left(\frac{1}{2} - \hat{a}_i^\dag \hat{a}_i \right) \hat{h}_i \hat{h}_i^\dag, \\
    \hat{S}_i^- &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i \hat{h}_i^\dag, &\quad \hat{\tilde{n}}_i &= 1 - \hat{h}_i^\dag \hat{h}_i = \hat{h}_i \hat{h}_i^\dag.
\end{aligned}
```

On sublattice B (sites with odd indices ``\{1, 2, 5, ..., L - 1\}``):
```math
\begin{aligned}
	\hat{\tilde{c}}_{i\uparrow}^\dag &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\uparrow} &= \hat{h}_i^\dag \hat{P}_i \hat{a}_i, \\
	\hat{\tilde{c}}_{i\downarrow}^\dag &= \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\downarrow} &= \hat{h}_i^\dag \hat{P}_i, \\
    \hat{S}_i^+ &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i \hat{h}_i^\dag, &\quad \hat{S}_i^z &= \left(\hat{a}_i^\dag \hat{a}_i - \frac{1}{2}\right) \hat{h}_i \hat{h}_i^\dag, \\
    \hat{S}_i^- &= \hat{h}_i \hat{h}_i^\dag \hat{P}_i \hat{a}_i, &\quad \hat{\tilde{n}}_i &= 1 - \hat{h}_i^\dag \hat{h}_i = \hat{h}_i \hat{h}_i^\dag.
\end{aligned}
```

Above, ``\hat{P}_i = \sqrt{1 - \hat{a}_i^\dag \hat{a}_i}`` is a projection operator ensuring that we cannot leave the physically meaningful subspace of ``\{0, 1\}`` magnons at each site.
If one consideres ``\hat{a}`` to be a hardcore boson instead, then projection operators can be replaced with identity operators. Then, one can skip them in the notation.
Both approaches are equally valid and consistent with [tJMagnonHolon](@ref) code base.

