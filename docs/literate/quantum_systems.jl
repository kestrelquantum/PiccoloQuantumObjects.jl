# ```@meta
# CollapsedDocStrings = true
# ```
using PiccoloQuantumObjects
using SparseArrays # for visualization

#=
# Quantum Systems

The [`QuantumSystem`](@ref) type is used to represent a quantum system with a drift 
Hamiltonian and a set of drive Hamiltonians,

```math
H = H_{\text{drift}} + \sum_i a_i H_{\text{drives}}^{(i)}
```

They are the containers for the dynamics.

```@docs
QuantumSystem
```

=#
H_drift = GATES[:Z]
H_drives = [GATES[:X], GATES[:Y]]
system = QuantumSystem(H_drift, H_drives)


#=
## Open quantum systems
Construct a `QuantumSystem` with Lindblad operators.
=#
dissipation_operators = [GATES[:Z], annihilate(2)]
system = QuantumSystem(H_drift, H_drives, dissipation_operators)

# TODO: put warning box here
# The Hamiltonian part `system.H` excludes the Lindblad operators.

#=
Construct a `QuantumSystem` from a Hamiltonian function and ForwardDiff.jl
=# 
H(a) = GATES[:Z] + a[1] * GATES[:X] + a[2] * GATES[:Y]
system = QuantumSystem(H, 2)

#=
# Composite systems

A composite quantum system is constructed from a set of subsystems and their interactions.
A [`CompositeQuantumSystem`](@ref) can contain interaction in the form of drift and drive
Hamiltonians acting on the full Hilbert space. It also contains subsystems and their drift
and drive Hamiltonians, which are internally lifted to the full Hilbert space.

## The `lift` operation`

To lift operators acting on a subsystem into the full Hilbert space, use [`lift`](@ref).
=#

# create the a + a' operator acting on the 1st subsystem
subspace_levels = [3, 2]
lift(create(3) + annihilate(3), 1, subspace_levels)

# create IXI operator 
lift(PAULIS[:X], 2, 3) |> abs |> sparse

# create an XX operator acting on qubits 3 and 4 in a 4-qubit system
lift([PAULIS[:X], PAULIS[:X]], [3, 4], 4) |> abs |> sparse


#=
# Reachability tests

Whether a quantum system can be used to reach a target state or operator can be tested
by computing the dynamical Lie algebra, which is provided by the [`is_reachable`](@ref)
function.
```@docs
is_reachable
```
=#

# Y can be reached by commuting Z and X.
system = QuantumSystem(PAULIS[:Z], [PAULIS[:X]])
is_reachable(GATES[:Y], sys)

# Y cannot be reached by X alone.
system = QuantumSystem([PAULIS[:X]])
is_reachable(GATES[:Y], sys)
