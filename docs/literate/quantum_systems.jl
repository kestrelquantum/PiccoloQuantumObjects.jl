# ```@meta
# CollapsedDocStrings = true
# ```

# # Abstract Quantum Systems

using PiccoloQuantumObjects
using SparseArrays # for visualization
⊗ = kron;

#=

## Quantum Systems

The [`QuantumSystem`](@ref) type is used to represent a quantum system with a drift 
Hamiltonian and a set of drive Hamiltonians,

```math
H = H_{\text{drift}} + \sum_i a_i H_{\text{drives}}^{(i)}
```

```@docs; canonical = false
QuantumSystem
```

`QuantumSystem`'s are containers for quantum dynamics. Internally, they compute the
necessary isomorphisms to perform the dynamics in a real vector space.

=#

H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
system = QuantumSystem(H_drift, H_drives)

a_drives = [1, 0]
system.H(a_drives)

#=
To extract the drift and drive Hamiltonians from a `QuantumSystem`, use the 
[`get_drift`](@ref) and [`get_drives`](@ref) functions. 

=#

get_drift(system) |> sparse

#
drives = get_drives(system)
drives[1] |> sparse

# 
drives[2] |> sparse

#=
!!! note
    We can also construct a `QuantumSystem` directly from a Hamiltonian function. Internally,
    `ForwardDiff.jl` is used to compute the drives.
=#

H(a) = PAULIS[:Z] + a[1] * PAULIS[:X] + a[2] * PAULIS[:Y]
system = QuantumSystem(H, 2)
get_drives(system)[1] |> sparse

#=
## Open quantum systems

We can also construct an [`OpenQuantumSystem`](@ref) with Lindblad dynamics, enabling
a user to pass a list of dissipation operators.

```@docs; canonical = false
OpenQuantumSystem
```
=#

H_drives = [PAULIS[:X]]
dissipation_operators = [PAULIS[:Z], PAULIS[:X]]
system = OpenQuantumSystem(H_drives, dissipation_operators=dissipation_operators)
system.dissipation_operators[1] |> sparse

#=
!!! warning
    The Hamiltonian part `system.H` excludes the Lindblad operators. This is also true
    for functions that report properties of `system.H`, such as [`get_drift`](@ref), 
    [`get_drives`](@ref), and [`is_reachable`](@ref).
=#

get_drift(system) |> sparse


#=
## Composite quantum systems

A [`CompositeQuantumSystem`](@ref) is constructed from a list of subsystems and their 
interactions. The interaction, in the form of drift or drive Hamiltonian, acts on the full
Hilbert space. The subsystems, with their own drift and drive Hamiltonians, are internally
lifted to the full Hilbert space.

=#

system_1 = QuantumSystem([PAULIS[:X]])
system_2 = QuantumSystem([PAULIS[:Y]])
H_drift = PAULIS[:Z] ⊗ PAULIS[:Z]
system = CompositeQuantumSystem(H_drift, [system_1, system_2]);

# _The drift Hamiltonian is the ZZ coupling._
get_drift(system) |> sparse

# _The drives are the X and Y operators on the first and second subsystems._
drives = get_drives(system)
drives[1] |> sparse

#
drives[2] |> sparse

#=
### The `lift` operation

To lift operators acting on a subsystem into the full Hilbert space, use [`lift`](@ref).
=#

# _Create an `a + a'` operator acting on the 1st subsystem of a qutrit and qubit system._
subspace_levels = [3, 2]
lift(create(3) + annihilate(3), 1, subspace_levels) .|> real |> sparse

# _Create IXI operator on the 2nd qubit in a 3-qubit system._
lift(PAULIS[:X], 2, 3) .|> real |> sparse

# _Create an XX operator acting on qubits 3 and 4 in a 4-qubit system._
lift([PAULIS[:X], PAULIS[:X]], [3, 4], 4) .|> real |> sparse


#=
# Reachability tests

Whether a quantum system can be used to reach a target state or operator can be tested
by computing the dynamical Lie algebra, which is provided by the [`is_reachable`](@ref)
function.
```@docs; canonical = false
is_reachable
```
=#

# _Y can be reached by commuting Z and X._
system = QuantumSystem(PAULIS[:Z], [PAULIS[:X]])
is_reachable(PAULIS[:Y], system)

# _Y cannot be reached by X alone._
system = QuantumSystem([PAULIS[:X]])
is_reachable(PAULIS[:Y], system)
