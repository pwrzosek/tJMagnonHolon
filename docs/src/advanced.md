# Advanced

The purpose of this section is to provide additional information about magnon-holon representation.
We will cover here details of the transformation to magnon-holon basis, including:
- sublattice-wise transformation of standard operators,
- standard ``t``-``J`` Hamiltonian, its magnon-holon representation and generalization to ``t``-``J_{xy}``-``J_{z}(\lambda)`` model,
- quantities preserved by the generalized Hamiltonian and commutation relation with translation operator,
- definition and numerical representation of momentum states and representative states in magnon-holon basis,
- example derivations for operators action on magnon-holon representative basis states.

---

Let us start by denoting annihilation operators for magnon (hard-core boson) as ``\hat{a}`` and holon (spinless fermion) as ``\hat{h}``. 
The *hard-core* nature of the boson comes from ``s = 1/2`` spin in the considered model.
What we call the magnon-holon transformation is the procedure of rewriting operators in terms of ``\hat{a}`` and ``\hat{h}`` operators.
Accordingly, we say that operator ``\hat{O}`` expressed in terms of ``\hat{a}`` and ``\hat{h}`` is in magnon-holon representation (or in magnon-holon basis).
Considering magnons as excitations around perfect antiferromagnetic spin background, below set of equations define the magnon-holon representation 
of the standard operators in the [tJMagnonHolon](@ref). We assume the system is one-dimensional with ``L = 2N`` sites, 0-based indexed, and ``N`` - positive integer.

On sublattice A (sites with even indices ``\{0, 2, 4, ..., L - 2\}``):
```math
\begin{aligned}
	\hat{\tilde{c}}_{i\uparrow}^\dag &= \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\uparrow} &= \hat{h}_i^\dag \hat{P}_i, \\
	\hat{\tilde{c}}_{i\downarrow}^\dag &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\downarrow} &= \hat{h}_i^\dag \hat{P}_i \hat{a}_i, \\
    \hat{S}_i^+ &= \hat{h}_i \hat{h}_i^\dag \hat{a}_i, &\quad \hat{S}_i^z &= \left(\frac{1}{2} - \hat{a}_i^\dag \hat{a}_i \right) \hat{h}_i \hat{h}_i^\dag, \\
    \hat{S}_i^- &= \hat{a}_i^\dag \hat{h}_i \hat{h}_i^\dag, &\quad \hat{\tilde{n}}_i &= 1 - \hat{h}_i^\dag \hat{h}_i = \hat{h}_i \hat{h}_i^\dag.
\end{aligned}
```

On sublattice B (sites with odd indices ``\{1, 2, 5, ..., L - 1\}``):
```math
\begin{aligned}
	\hat{\tilde{c}}_{i\uparrow}^\dag &= \hat{a}_i^\dag \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\uparrow} &= \hat{h}_i^\dag \hat{P}_i \hat{a}_i, \\
	\hat{\tilde{c}}_{i\downarrow}^\dag &= \hat{P}_i \hat{h}_i, &\quad \hat{\tilde{c}}_{i\downarrow} &= \hat{h}_i^\dag \hat{P}_i, \\
    \hat{S}_i^+ &= \hat{a}_i^\dag \hat{h}_i \hat{h}_i^\dag, &\quad \hat{S}_i^z &= \left(\hat{a}_i^\dag \hat{a}_i - \frac{1}{2}\right) \hat{h}_i \hat{h}_i^\dag, \\
    \hat{S}_i^- &= \hat{h}_i \hat{h}_i^\dag \hat{a}_i, &\quad \hat{\tilde{n}}_i &= 1 - \hat{h}_i^\dag \hat{h}_i = \hat{h}_i \hat{h}_i^\dag.
\end{aligned}
```

On both sublattices, ``\hat{P}_i = 1 - \hat{a}_i^\dag \hat{a}_i``.

---

With the above definition we can take the standard ``t``-``J`` Hamiltonian formulation,
```math
\hat{H} = -t \sum_{\langle i,j \rangle, \sigma} \left( \hat{\tilde{c}}_{i\sigma}^{\dag} \hat{\tilde{c}}_{j\sigma} + \mathrm{H.c.} \right)
        +  J \sum_{\langle i,j \rangle} \left( \hat{S}^x_i \hat{S}^x_j + \hat{S}^y_i \hat{S}^y_j + \hat{S}^z_i \hat{S}^z_j - \frac{1}{4} \hat{\tilde{n}}_i \hat{\tilde{n}}_j \right),
```
and rewrite it to magnon-holon basis,
```math
\hat{H} = \hat{H}_t + \hat{H}_{J},
```
```math
\hat{H}_t = t \sum_{\langle i,j \rangle} \hat{P}_i \left( \hat{h}_i^\dag \hat{h}_j \hat{a}_i + \hat{h}_i^\dag \hat{h}_j \hat{a}_j^\dag \right) \hat{P}_j + \textrm{H.c.},
```
```math
\hat{H}_{J} = \frac{J}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left( \hat{a}_i \hat{a}_j + \hat{a}_i^\dag \hat{a}_j^\dag 
            + \hat{a}_i^\dag \hat{a}_i + \hat{a}_j^\dag \hat{a}_j - 2 \hat{a}_i^\dag \hat{a}_i \hat{a}_j^\dag \hat{a}_j - 1 \right) \hat{h}_j \hat{h}_j^\dag.
```

We extend the above model by introducing two additional components:
- ``\textcolor{orange}{\alpha}`` - anisotropy between quatization axis and perpendicular axes (XXZ anisotropy),
- ``\textcolor{orange}{\lambda}`` - scaling parameter for interactions between neighboring magnons (magnon-magnon interaction).

```math
\hat{H} = \hat{H}_t + \hat{H}_{xy} + \hat{H}_z
```
```math
\hat{H}_t = t \sum_{\langle i,j \rangle} \hat{P}_i \left( \hat{h}_i^\dag \hat{h}_j \hat{a}_i + \hat{h}_i^\dag \hat{h}_j \hat{a}_j^\dag \right) \hat{P}_j + \textrm{H.c.},
```
```math
\hat{H}_{xy} = \frac{\textcolor{orange}{\alpha} J}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left( \hat{a}_i \hat{a}_j + \hat{a}_i^\dag \hat{a}_j^\dag \right) \hat{h}_j \hat{h}_j^\dag
```
```math
\hat{H}_{z} = \frac{J}{2} \sum_{\langle i,j \rangle} \hat{h}_i \hat{h}_i^\dag \left(\hat{a}_i^\dag \hat{a}_i + \hat{a}_j^\dag \hat{a}_j - 2 \textcolor{orange}{\lambda} \hat{a}_i^\dag \hat{a}_i \hat{a}_j^\dag \hat{a}_j - 1 \right) \hat{h}_j \hat{h}_j^\dag
```

---

For numerical diagnalization of an abstract model like the one above, it is important to determine its symmetries. Those symmetries generate quatities that are conserved by the model. 
When writing the matrix of the model Hamiltonian in a basis, the basis may incorporate those symmetries, i.e. the basis states can be taken as eigen-states of the symmetry operators.
Eigen-states that correspond to different eigen-values of the symmetry operator are orthogonal. Thus Hamiltonian written in such a basis will consist of orthogonal blocks that can be
enumerated by the eigen-values of the symmetry operators used for the basis construction. Each such block can be diagonalized separately, decreasing the amount of memory needed
to store the Hamiltonian matrix and resulting eigen-states. But what constitutes a symmetry of our model? The word symmetry may be misleading in this case, since the symmetry does not have
to reflect only geometrical properties. What we really look for when asking about the symmetries of the Hamiltonian, are operators that commute with the Hamiltonian,
```math
[\hat{H}, \hat{O}] = 0.
```

Some of the conserved quantities are pretty easy to notice. For example, size ``L`` of the system is one of them. What is the operator for the system size? We don't have to define it explicitely.
It is enough for us that ``L`` enumerates its eigen-values.

Next, we can observe that our Hamiltonian does not change the number of electrons (holes). The corresponding operator is ``\hat{N} = \sum_i \hat{\tilde{n}}_i``. One can check that,
```math
[\hat{H}, \hat{N}] = 0.
```
The eigen-values of ``\hat{N}`` are simply ``N_e \in \{0, 1, 2, ..., L\}``. This operator does not influence the system size either.
Thus subspace for given ``L`` can be split into futher subspaces with different electron occupancy ``N_e``.

In similar way, once we set a system size ``L`` and number of electrons ``N_e``, such subspace can be split into subspaces with different magnetization ``S^z_{\mathrm{tot}}``.
The corresponding operator is ``\hat{S}^z_{\mathrm{tot}} = \sum_i \hat{S}^z_i``. It indeed commutes with the Hamiltonian and does not break considered so far symmetries 
(it commutes with their operators too),
```math
[\hat{H}, \hat{S}^z_{\mathrm{tot}}] = 0,
```
```math
[\hat{N}, \hat{S}^z_{\mathrm{tot}}] = 0.
```

The standard Hamiltonian (and other previously mentioned operators) commute also with translation operator ``\hat{T}``. But it is not true for the generalized Hamiltonian.
Generalized model commutes with ``\hat{T}`` only for ``\lambda = 1``.
This is related to distinguished sublattices in the transformation to magnon-holon basis.
At the same time, the square of the translation operator, ``\hat{D} = \hat{T}^2`` commutes with the generalized
Hamiltonian for all values of ``\lambda``,
```math
[\hat{H}, \hat{D}] = 0.
```
We will call ``\hat{D}`` an *even translation* operator. It of course commute with all operatators that ``\hat{T}`` commutes with.

All the above symmetries are incorporated in the [tJMagnonHolon](@ref). But there is still a room for improvements. For example, parity symmetry is another symmetry that can be taken into account.
It is important to remember that additional symmetries will make the basis in itself more complex. This can make the obtained results difficult to analyze or interpret.
Thus stopping at the translation operator is optimal. Allowed system sizes are only few sites smaller than maximum for available computers. At the same time, 
basis construction is still easy to understand.

---

Let us now construct the basis for our problem. Our requirements are that all the basis states are simultaneously eigen-states of all the mentioned symmetry operators:
system size, number of electrons ``\hat{N}``, magnetization ``\hat{S}^z_{\mathrm{tot}}`` and even translation ``\hat{D}``. Let us start with representation that will
automatically include the first three of them. We will use the following notation:
- system size: ``L = 2N, N \in \mathbb{N}``,
- number of electrons: ``N_e \in \{0, 1, ..., L\}``,
- number of spins up: ``n_\uparrow \in \{0, 1, ..., N_e\}``.
Notice that for given ``N_e`` the magnetization ``S^z_{\mathrm{tot}}`` is uniquely defiend by ``n_\uparrow``,
```math
    S^z_{\mathrm{tot}} \propto n_\uparrow - n_\downarrow = n_\uparrow - (N_e - n_\uparrow) = 2n_\uparrow - N_e,
```
where ``n_\downarrow`` stands for number of spins down in the system. 

We will use two 64-bit numbers to represent configuration of electrons and spins. This will set a limit of 64
to available system sizes, but apart from trivial cases we won't be able to calculate more than around 20-30 sites anyway. 
Let us name the two numbers as `charges` and `spins` respectively.
We will use ``L`` first bits of each of the number to represent presence (absence) of electrons and spins up.
Our rules are simple. We connect bit values of ``\{0, 1\}`` with occupancy.
- for `charges`:
    - 1 at ``i``-th bit ``\to`` electron at ``i``-th site,
    - 0 at ``i``-th bit ``\to`` hole at ``i``-th site,
- for `spins`:
    - 1 at ``i``-th bit ``\to`` spin up at ``i``-th site,
    - 0 at ``i``-th bit ``\to`` spin down ``i``-th site.

Since bits are indexed from 0, we use the same indexing for sites. It may not be the most intuitive choice for physicists, but it is the most practical one (and natural for most programmers).
It is also important to specify that only electron can carry a spin (holons are spinless fermions). To disambiguate the definition of the site occupied by a hole, we always set
the bit of the `spin` variable to ``0`` if corresponding bit of `charges` variable is ``0``.
Therefore only 3 out of 4 combinations are meaningful for us (``i``[`charges`], ``i``[`spins`]) ``\in\{(0,0), (1,0), (1,1)\}``.

Before we follow with translation operator, let us transform the above representation to magnon-holon basis.
We look at the details of the transformation to magnon-holon basis (operators defined at the begining), and see that the transformation does not alter the electron/hole configuration. 
We therefore use the above defined `charges` without changes.
(Of course one can swap the meaning of ``0`` and ``1`` bit values, but using ``0`` for hole seems more intuitive).
On the other hand, we have no spins in magnon-holon basis, but magnons instead. Let us define `magnons` - number that represents configuration of magnons in the system.
On odd sites spins up directly transform to magnons, so there's nothing to do there. On even sites on the other hand,
spins down transform to magnons, so we have to transform all zeros to ones and all ones to zeros on even sites. 
We call this last step a rotation of sublattice A (sublattice of even-indexed sites assuming 0-based indexing).
There is of course a deeper meaning behind calling it a rotation but it won't be important here.

After the rotation of sublattice A is performed, we land at the following representation,
- for `charges`:
    - 1 at ``i``-th bit ``\to`` electron at ``i``-th site,
    - 0 at ``i``-th bit ``\to`` hole at ``i``-th site,
- for `magnons`:
    - 1 at ``i``-th bit ``\to`` magnon at ``i``-th site,
    - 0 at ``i``-th bit ``\to`` no magnon ``i``-th site.
We disambiguate the hole representation in the same way as before (there cannot be both the hole and the magnon at the same site).
The above corresponds to the [`State`](@ref Main.tJmodel1D.State) structure in the [tJMagnonHolon](@ref) code.

Let us now follow with introduction of *momentum states* - eigen-states of the even translation operator ``\hat{D}``.
Let ``\vert s, L, N_e, n_\uparrow \rangle`` stand for a `state`  with certain configuration of `charges` and `magnons` denoted by `s` 
for a system with ``L = 2N, N \in \mathbb{N}`` sites, ``N_e`` electrons, and ``n_\uparrow`` spins up. Let us also shortly write ``\vert s \rangle``
when exact values of ``L, N_e, n_\uparrow`` are not important. Notice that set of all `s` can be well-ordered since `s` is represented by a pair of integers.
We define the magnon-holon momentum states in a standard way (we just use ``\hat{D}`` instead of ``\hat{T}``),
```math
\vert s(k) \rangle = \frac{1}{\sqrt{N_s}} \sum_{r = 0}^{N - 1} \exp\left(-\frac{2\pi i}{N}kr\right) \hat{D}^r \vert s \rangle.
```
One can see that indeed ``\vert s(k) \rangle`` are eigen-states of ``\hat{D}``,
```math
\hat{D} \vert s(k) \rangle = \exp\left(\frac{2\pi i}{N}k\right) \vert s(k) \rangle,
```
with eigen-values ``\exp\left(\frac{2\pi i}{N}k\right)``, where ``k \in \{0, 1, ..., N-1\}``. This ``k`` variable enumerates the subspaces generated by even translation operator.
When `system` variable of type [`System`](@ref Main.tJmodel1D.System) is created, ``k`` stands for `system.momentum` field. 

Two more steps are required to form a basis out ot the above defined momentum states. 

- First: we need to make sure that they have unique representation (for linear independence of basis states).
For this we need to pick a *representative* state ``\vert \tilde{s} \rangle`` out of the set of even translations ``\hat{D}^r \vert s \rangle`` of particular configuration ``s``.
Now the fact that set of all ``s`` can be well-ordered comes in handy. We define our ordering relation between ``s_1`` and ``s_2``. It will allow us to answer e.g. if ``s_1`` is greater than ``s_2``.
For example we can compare `charges` of the two configurations, and if `charges` are the same we compare `magnons`. With this we can pick the representative state ``\vert \tilde{s} \rangle``
such that ``\tilde{s}`` it the *smallest* a out of all translations of particular configuration ``s``.

- Second: we need to make sure that state ``\vert s \rangle`` is compatible with momentum subspace ``k \in \{0, ..., N-1\}``.
This is connected to normalizability of momentum states. Let us define periodicity of state ``\vert s \rangle``,
as the minimal number ``R_s > 0`` of even translations that transform ``\vert s \rangle`` onto itself,
```math
\min \{ R_s : R_s > 0 \land \hat{D}^{R_s} \vert s \rangle = \vert s \rangle \}.
```
The compatibility of ``\vert s \rangle`` with ``k`` requires that,
```math
kR_s\mod N = 0.
```
If the above is fulfilled, the normalization factor ``N_s = \frac{N^2}{R_s}``. If the above is *not* fulfilled, then ``\vert s(k) \rangle \equiv 0`` is not normalizable.

Set of representative states compatible with ``k`` form a momentum basis for subspace with given ``L, N_e, S^z_{\mathrm{tot}}, k``.

---

Equiped with the introduced notation, we are ready to derive some expressions for operators action on magnon-holon (momentum) basis states.
To shorten the notation, from now on, we write ``\exp(-ikr)`` instead of ``\exp\left(-\frac{2\pi i}{N}kr\right)``.
It should be clear to which momentum subspace we refer even if we use redefined ``k`` to enumerate it.
Just remember that in the code `momentum` field of [`System`](@ref Main.tJmodel1D.System) is an integer from set ``\{0, ..., N-1\}``.

We start with the simplest case of ``\hat{S}^z_k`` operator,
```math
\hat{S}^z_k = \frac{1}{\sqrt{L}} \sum_{r=0}^{L-1} \exp(-ikr) \hat{S}^z_r.
```
We evaluate its action on a magnon-holon momentum basis state ``\vert \tilde{s}(p) \rangle = \vert \tilde{s}(p), L, N_e, n_\uparrow \rangle``,
```math
\begin{aligned}
\hat{S}^z_k \vert \tilde{s}(p) \rangle &= \frac{1}{\sqrt{N_s}} \sum_{r=0}^{N-1} \exp(-ipr) \hat{S}^z_k \hat{D}^r \vert \tilde{s} \rangle = \frac{1}{\sqrt{L N_s}} \sum_{r=0}^{N-1} \sum_{r'=0}^{L-1} \exp(-ipr - ikr') \hat{S}^z_{r'} \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') (-1)^{r'} \left( \frac{1}{2} - \hat{a}^{\dag}_{r'} \hat{a}_{r'} \right) \hat{h}_{r'} \hat{h}^{\dag}_{r'} \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') \hat{D}^r (-1)^{r'} \left( \frac{1}{2} - \hat{a}^{\dag}_{2r + r'} \hat{a}_{2r + r'} \right) \hat{h}_{2r + r'} \hat{h}^{\dag}_{2r + r'} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp[-i(p - 2k)r -ik(2r + r')] \hat{D}^r (-1)^{r'} \left( \frac{1}{2} - \hat{a}^{\dag}_{2r + r'} \hat{a}_{2r + r'} \right) \hat{h}_{2r + r'} \hat{h}^{\dag}_{2r + r'} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p - 2k)r -ikR] \hat{D}^r (-1)^{R - 2r} \left( \frac{1}{2} - \hat{a}^{\dag}_{R} \hat{a}_{R} \right) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p - 2k)r -ikR] \hat{D}^r (-1)^{R} \left( \frac{1}{2} - \hat{a}^{\dag}_{R} \hat{a}_{R} \right) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L}} \sum_{R=0}^{L-1} \exp(-ikR) \alpha_R(\tilde{s}) \vert \tilde{s}(p-2k) \rangle,
\end{aligned}
```
where ``(-1)^{R} \left( \frac{1}{2} - \hat{a}_{R} \hat{a}^{\dag}_{R} \right) \hat{h}^{\dag}_{R} \hat{h}_{R} \vert \tilde{s} \rangle = \alpha_R(\tilde{s}) \vert \tilde{s} \rangle``.
We intentionally provide all the transformations step by step to showcase the derivation procedure. We use the definition of the momentum representation (Fourier transform) for both the operator and the state.
We follow with the transformation of the operators to magnon-holon basis. Then we commute the even translation operator ``\hat{D}`` with other operators. Next 3 lines are probably least transparent.
We change the summation order, rendering the action of the operators (other than even translation operator) independent from the Fourier transform to momentum space.
Finally, we evaluate the action of operators on the state ``\vert \tilde{s} \rangle`` and we write the state back in the momentum space.
It is immediately visible that ``\alpha_R(\tilde{s})`` evaluates to ``0`` when site ``R`` in state ``\vert \tilde{s} \rangle`` is occupied by a hole. When site ``R`` is not occupied by the hole,
possible values of ``\alpha_R(\tilde{s})`` are given according to the table below,

| ``R`` in sublattice | magnons at ``R`` in ``\vert \tilde{s} \rangle`` | ``\alpha_R(\tilde{s})`` |
| :------------------ | :---------------------------------------------: | ----------------------: |
| ``R \in A``         | ``0``                                           | ``\frac{1}{2}``         |
| ``R \in B``         | ``0``                                           | ``-\frac{1}{2}``        |
| ``R \in A``         | ``1``                                           | ``-\frac{1}{2}``        |
| ``R \in B``         | ``1``                                           | ``\frac{1}{2}``         |

Using the result for ``\hat{S}^z_k \vert \tilde{s}(p) \rangle``, we easily obtain,
```math
\hat{S}^z_r \vert \tilde{s}(p) \rangle = \frac{1}{\sqrt{L}} \sum_{q} \exp(iqr) \hat{S}^z_q \vert \tilde{s}(p) \rangle = \frac{1}{L} \sum_{q,R} \exp[-iq(R-r)]\alpha_R(\tilde{s}) \vert \tilde{s} (p-2q) \rangle.
```
In the above sum, we assert the periodicity-momentum match between state ``\vert \tilde{s} \rangle`` and  momentum ``p-2q`` by excluding from the sum terms that do not fulfill periodicity-momentum relation.

Let us now work out the ladder oparator ``\hat{S}^+_k``. Compared to the previous example, it will require few more steps, since it affects the number of spins up in the system.
```math
\begin{aligned}
\hat{S}^+_k \vert \tilde{s}(p) \rangle &= \frac{1}{\sqrt{N_s}} \sum_{r=0}^{N-1} \exp(-ipr) \hat{S}^+_k \hat{D}^r \vert \tilde{s} \rangle = \frac{1}{\sqrt{L N_s}} \sum_{r=0}^{N-1} \sum_{r'=0}^{L-1} \exp(-ipr - ikr') \hat{S}^+_{r'} \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') (\delta_{r' \in A} \hat{a}_{r'} +\delta_{r' \in B} \hat{a}^{\dag}_{r'} ) \hat{h}_{r'} \hat{h}^{\dag}_{r'} \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') \hat{D}^r (\delta_{r' \in A} \hat{a}_{2r + r'} + \delta_{r' \in B} \hat{a}^{\dag}_{2r + r'} ) \hat{h}_{2r + r'} \hat{h}^{\dag}_{2r + r'} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r (\delta_{R-2r \in A} \hat{a}_{R} + \delta_{R-2r \in B} \hat{a}^{\dag}_{R} ) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r \left( \frac{1+(-1)^R}{2} \hat{a}_{R} + \frac{1-(-1)^R}{2} \hat{a}^{\dag}_{R} \right) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r \frac{1}{2}(\hat{a}_{R} + \hat{a}^{\dag}_{R}) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r \frac{(-1)^R}{2}(\hat{a}_{R} - \hat{a}^{\dag}_{R}) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r [\alpha^+_R(\tilde{s}) \vert s^+_R \rangle + \alpha^-_R(\tilde{s}) \vert s^-_R \rangle] \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r [\beta^+_R(\tilde{s}) \vert s^+_R \rangle + \beta^-_R(\tilde{s}) \vert s^-_R \rangle] \\
&= \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^+_R}} \exp(-ikR) [\alpha^+_R(\tilde{s}) + \beta^+_R(\tilde{s})] \vert s^+_R (p - 2k) \rangle \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^-_R}} \exp(-ikR) [\alpha^-_R(\tilde{s}) + \beta^-_R(\tilde{s})] \vert s^-_R (p - 2k) \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^+_R}} \exp[-ikR -i(p-2k)d_{s^+_R}] [\alpha^+_R(\tilde{s}) + \beta^+_R(\tilde{s})] \vert \tilde{s}^+_R (p - 2k) \rangle \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^-_R}} \exp[-ikR -i(p-2k)d_{s^-_R}] [\alpha^-_R(\tilde{s}) + \beta^-_R(\tilde{s})] \vert \tilde{s}^-_R (p - 2k) \rangle,
\end{aligned}
```
with ``\frac{1}{2}(\hat{a}_{R} + \hat{a}^{\dag}_{R}) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle = [\alpha^+_R(\tilde{s}) \vert s^+_R \rangle + \alpha^-_R(\tilde{s}) \vert s^-_R \rangle]``
and ``\frac{(-1)^R}{2}(\hat{a}_{R} - \hat{a}^{\dag}_{R}) \hat{h}_{R} \hat{h}^{\dag}_{R} \vert \tilde{s} \rangle = [\beta^+_R(\tilde{s}) \vert s^+_R \rangle + \beta^-_R(\tilde{s}) \vert s^-_R \rangle]``,
where ``\vert s^{\pm}_R \rangle \equiv \vert s^{\pm}_R, L, N_e, n_\uparrow \pm 1 \rangle`` comes from flipping a spin at site ``R`` in state ``\vert \tilde{s} \rangle``. 
The normalization factors ``\sqrt{N_{s^{\pm}_R}}`` appear since new states ``\vert s^{\pm}_R \rangle`` may in general have periodicity different than ``\vert \tilde{s} \rangle``. 
These states may also not be representative states, so we need to find their representatives ``\vert \tilde{s}^{\pm}_R \rangle`` and include proper phase shifts,
as determined by the distance ``d_{s^{\pm}_R}`` to the representative.
Of course, ``\beta^-_R(\tilde{s}) = -\alpha^-_R(\tilde{s})`` meaning that we can only increase the number of spins up (as expected from ``\hat{S}^+_k``). On the other hand,
``\beta^+_R(\tilde{s}) = \alpha^+_R(\tilde{s})`` with ``\alpha^+_R(\tilde{s})`` equal to ``0`` when ``R`` is occupied by a hole or otherwise it is given according to the table below.

| ``R`` in sublattice | magnons at ``R`` in ``\vert \tilde{s} \rangle`` | ``\alpha^+_R(\tilde{s})`` |
| :------------------ | :---------------------------------------------: | ------------------------: |
| ``R \in A``         | ``0``                                           | ``0``                     |
| ``R \in B``         | ``0``                                           | ``\frac{1}{2}``           |
| ``R \in A``         | ``1``                                           | ``\frac{1}{2}``           |
| ``R \in B``         | ``1``                                           | ``0``                     |

From the above considerations, we see that the code for this operator shall first check if the value of ``\alpha^+_R(\tilde{s})`` is non-zero and only proceed if that's the case.
Then the operators action on state ``\vert \tilde{s} \rangle`` simply alternates the ``R``-th bit of `magnons` variable, which equates to XORing `magnons` with ``2^R``.

As one can expect, derivation for ``\hat{S}^-_k`` is analogical, with the same result form and following relations between coefficients,
```math
\alpha^-_R(\tilde{s}) = \beta^-_R(\tilde{s}) = \frac{1}{2} - \alpha^+_R(\tilde{s}), \\
\beta^+_R(\tilde{s}) = -\alpha^+_R(\tilde{s}).
```
Also, ``\hat{S}^{\pm}_r`` can be obtained the same way as ``\hat{S}^z_r``.

Creation and annihilation operators for electrons follow the same scheme as spin ladder operators. For example,
```math
\begin{aligned}
\hat{\tilde{c}}_{k\downarrow} \vert \tilde{s}(p) \rangle &= \frac{1}{\sqrt{N_s}} \sum_{r=0}^{N-1} \exp(-ipr) \hat{\tilde{c}}_{k\downarrow} \hat{D}^r \vert \tilde{s} \rangle = \frac{1}{\sqrt{L N_s}} \sum_{r=0}^{N-1} \sum_{r'=0}^{L-1} \exp(-ipr - ikr') \hat{\tilde{c}}_{r'\downarrow} \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') (\delta_{r' \in A} \hat{h}^{\dag}_{r'} \hat{P}_{r'} \hat{a}_{r'} + \delta_{r' \in B} \hat{h}^{\dag}_{r'} \hat{P}_{r'}) \hat{D}^r \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,r'} \exp(-ipr -ikr') \hat{D}^r (\delta_{r' \in A} \hat{h}^{\dag}_{2r + r'} \hat{P}_{2r + r'} \hat{a}_{2r + r'} + \delta_{r' \in B} \hat{h}^{\dag}_{2r + r'} \hat{P}_{2r + r'}) \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p - 2k)r -ikR] \hat{D}^r (\delta_{R-2r \in A} \hat{h}^{\dag}_{R} \hat{P}_{R} \hat{a}_{R} + \delta_{R-2r \in B} \hat{h}^{\dag}_{R} \hat{P}_{R}) \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p - 2k)r -ikR] \hat{D}^r \frac{1}{2} \hat{h}^{\dag}_{R} \hat{P}_{R} [\hat{a}_{R} + 1 + (-1)^R (\hat{a}_{R} - 1)] \vert \tilde{s} \rangle \\
&= \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r [\alpha^+_R(\tilde{s}) \vert s^+_R \rangle + \alpha^-_R(\tilde{s}) \vert s^-_R \rangle] \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{r,R} \exp[-i(p-2k)r -ikR] \hat{D}^r [\beta^+_R(\tilde{s}) \vert s^+_R \rangle + \beta^-_R(\tilde{s}) \vert s^-_R \rangle] \\
&= \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^+_R}} \exp[-ikR -i(p-2k)d_{s^+_R}] [\alpha^+_R(\tilde{s}) + \beta^+_R(\tilde{s})] \vert \tilde{s}^+_R (p - 2k) \rangle \\
&+ \frac{1}{\sqrt{L N_s}} \sum_{R} \sqrt{N_{s^-_R}} \exp[-ikR -i(p-2k)d_{s^-_R}] [\alpha^-_R(\tilde{s}) + \beta^-_R(\tilde{s})] \vert \tilde{s}^-_R (p - 2k) \rangle.
\end{aligned}
```
where the hole is added at site R in state ``\vert s^{\pm}_R \rangle`` resulting in the increased(+)/decreased(-) magnetization with respect to ``\vert \tilde{s} \rangle``. 
The only non-zero contributions come from sites not occupied by the holes. Of course, even then the terms for decreased magnetization state must vanish, and indeed we have
``\beta^-_R(\tilde{s}) = -\alpha^-_R(\tilde{s})``. In general, the coefficients in the above result follow the same relations as for ``\hat{S}^+_k``. On the other hand,
for creation operator ``\hat{\tilde{c}}^{\dag}_{k\downarrow}``, only sites occupied by a hole give contribution. Then we are guaranteed that there are no magnons and in such case,
``\beta^+_R(\tilde{s}) = -\alpha^+_R(\tilde{s})`` and simply ``\beta^-_R(\tilde{s}) = \alpha^-_R(\tilde{s}) = \frac{1}{2}``.

In the end, let us discuss operators with more than one index. Such operators may appear when calculating n-particle correlation functions.
For example, let ``\hat{O}_{k,q} = \hat{\tilde{c}}_{k\uparrow}\hat{\tilde{c}}_{q\downarrow}``. Applying previous results, we obtain,
```math
\begin{aligned}
\hat{O}_{k,q} \vert \tilde{s}(p) \rangle &= \hat{\tilde{c}}_{k\uparrow} \sum_{R} \sqrt{\frac{N_{s^+_R}}{L N_s}} \exp[-iqR -i(p-2q)d_{s^+_R}] [\alpha^+_R(\tilde{s}) + \beta^+_R(\tilde{s})] \vert \tilde{s}^+_R (p - 2q) \rangle \\
&= \sum_{R,R'} \sqrt{\frac{N_{s^{+-}_{R,R'}}}{L^2 N_s}} \exp[-iqR -i(p-2q)d_{s^+_R}] \exp[-ikR' -i(p-2q-2k)d_{s^{+-}_{R,R'}}] \times \\
&\times [\alpha^+_R(\tilde{s}) + \beta^+_R(\tilde{s})] [\alpha^-_{R'}(\tilde{s}^+_R) + \beta^-_{R'}(\tilde{s}^+_R)] \vert \tilde{s}^{+-}_{R,R'} (p - 2q - 2k) \rangle.
\end{aligned}
```
It is important to mention that periodicity-momentum correspondence has to be checked for all the intermediate states, not only for the final state.
This means that the sum over ``R'`` only makes sense if state ``\vert \tilde{s}^+_R \rangle`` is compatible with momentum ``p-2q``. Otherwise, one may introduce
false non-zero contributions to the final result. Other than that, there are no suprises compared to previous examples. 

This summarizes our derivations. If you plan to introduce your own custom operators, read through [Custom Operators](@ref) subsection of the [Guide](@ref). 
You will learn there how the opreators functions are desined in [tJMagnonHolon](@ref). Then you may want to compare presented here derivations,
with corresponding operators functions in the `operators.jl` file in `.../tJMagnonHolon/src/modules/`.

