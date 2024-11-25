```@meta
EditURL = "../../literate/quantum_systems.jl"
```

```@meta
CollapsedDocStrings = true
```

# `Abstract Quantum Systems`

````@example quantum_systems
using PiccoloQuantumObjects
using SparseArrays # for visualization
⊗ = kron;
nothing #hide
````

```@docs; canonical = false
AbstractQuantumSystem
```

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

````@example quantum_systems
H_drift = PAULIS[:Z]
H_drives = [PAULIS[:X], PAULIS[:Y]]
system = QuantumSystem(H_drift, H_drives)

a_drives = [1, 0]
system.H(a_drives)
````

To extract the drift and drive Hamiltonians from a `QuantumSystem`, use the
[`get_drift`](@ref) and [`get_drives`](@ref) functions.

````@example quantum_systems
get_drift(system) |> sparse
````

_Get the X drive._

````@example quantum_systems
drives = get_drives(system)
drives[1] |> sparse
````

_And the Y drive._

````@example quantum_systems
drives[2] |> sparse
````

!!! note
    We can also construct a `QuantumSystem` directly from a Hamiltonian function. Internally,
    `ForwardDiff.jl` is used to compute the drives.

````@example quantum_systems
H(a) = PAULIS[:Z] + a[1] * PAULIS[:X] + a[2] * PAULIS[:Y]
system = QuantumSystem(H, 2)
get_drives(system)[1] |> sparse
````

## Open quantum systems

We can also construct an [`OpenQuantumSystem`](@ref) with Lindblad dynamics, enabling
a user to pass a list of dissipation operators.

```@docs; canonical = false
OpenQuantumSystem
```

_Add a dephasing and annihilation error channel._

````@example quantum_systems
H_drives = [PAULIS[:X]]
a = annihilate(2)
dissipation_operators = [a'a, a]
system = OpenQuantumSystem(H_drives, dissipation_operators=dissipation_operators)
system.dissipation_operators[1] |> sparse
````

````@example quantum_systems
system.dissipation_operators[2] |> sparse
````

!!! warning
    The Hamiltonian part `system.H` excludes the Lindblad operators. This is also true
    for functions that report properties of `system.H`, such as [`get_drift`](@ref),
    [`get_drives`](@ref), and [`is_reachable`](@ref).

````@example quantum_systems
get_drift(system) |> sparse
````

## Composite quantum systems

A [`CompositeQuantumSystem`](@ref) is constructed from a list of subsystems and their
interactions. The interaction, in the form of drift or drive Hamiltonian, acts on the full
Hilbert space. The subsystems, with their own drift and drive Hamiltonians, are internally
lifted to the full Hilbert space.

````@example quantum_systems
system_1 = QuantumSystem([PAULIS[:X]])
system_2 = QuantumSystem([PAULIS[:Y]])
H_drift = PAULIS[:Z] ⊗ PAULIS[:Z]
system = CompositeQuantumSystem(H_drift, [system_1, system_2]);
nothing #hide
````

_The drift Hamiltonian is the ZZ coupling._

````@example quantum_systems
get_drift(system) |> sparse
````

_The drives are the X and Y operators on the first and second subsystems._

````@example quantum_systems
drives = get_drives(system)
drives[1] |> sparse
````

````@example quantum_systems
drives[2] |> sparse
````

### The `lift` function

To lift operators acting on a subsystem into the full Hilbert space, use [`lift`](@ref).
```@docs; canonical = false
lift
```

_Create an `a + a'` operator acting on the 1st subsystem of a qutrit and qubit system._

````@example quantum_systems
subspace_levels = [3, 2]
lift(create(3) + annihilate(3), 1, subspace_levels) .|> real |> sparse
````

_Create IXI operator on the 2nd qubit in a 3-qubit system._

````@example quantum_systems
lift(PAULIS[:X], 2, 3) .|> real |> sparse
````

_Create an XX operator acting on qubits 3 and 4 in a 4-qubit system._

````@example quantum_systems
lift([PAULIS[:X], PAULIS[:X]], [3, 4], 4) .|> real |> sparse
````

We can also lift an operator that entangles different subspaces by passing the indices
of the entangled subsystems.

````@example quantum_systems
#_Here's another way to create an XX operator acting on qubits 3 and 4 in a 4-qubit system._
lift(kron(PAULIS[:X], PAULIS[:X]), [3, 4], 4) .|> real |> sparse
````

_Lift a CX gate acting on the 1st and 3rd qubits in a 3-qubit system._
_The result is independent of the state of the second qubit._

````@example quantum_systems
lift(GATES[:CX], [1, 3], 3) .|> real |> sparse
````

# Reachability tests

Whether a quantum system can be used to reach a target state or operator can be tested
by computing the dynamical Lie algebra. Access to this calculation is provided by the
[`is_reachable`](@ref) function.
```@docs; canonical = false
is_reachable
```

_Y can be reached by commuting Z and X._

````@example quantum_systems
system = QuantumSystem(PAULIS[:Z], [PAULIS[:X]])
is_reachable(PAULIS[:Y], system)
````

_Y cannot be reached by X alone._

````@example quantum_systems
system = QuantumSystem([PAULIS[:X]])
is_reachable(PAULIS[:Y], system)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

