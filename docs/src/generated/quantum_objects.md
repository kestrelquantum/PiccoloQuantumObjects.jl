```@meta
EditURL = "../../literate/quantum_objects.jl"
```

```@meta
CollapsedDocStrings = true
```

````@example quantum_objects
using PiccoloQuantumObjects
using SparseArrays # for visualization

#=
````

Quantum Objects

````@example quantum_objects
PiccoloQuantumObjects.jl provides a set of tools for working with quantum objects,
such as quantum states and operators.

# Quantum states

Construct quantum states from bitstrings or string representations using atomic notation
(ground state `g`, excited state `e`, etc.).
=#

ket_from_string("g", [2])

ket_from_string("gg", [2,2])

ket_from_string("(g+e)g", [2,2])

ket_from_bitstring("0")

ket_from_bitstring("01")
````

## Quantum Operators

Common operators are provided in [`PAULIS`](@ref) and [`GATES`](@ref).
```@docs
GATES
```

Quantum operators can also be constructed from strings.

````@example quantum_objects
operator_from_string("X")

operator_from_string("XZ")
````

Annihilation and creation operators are provided for oscillator systems.

````@example quantum_objects
annihilate(2)

annihilate(3)

create(2) + annihilate(2)

create(3)'annihilate(3)
````

Random unitary operators can be generated using the Haar measure.

````@example quantum_objects
haar_random(3)
````

If we want to generate random operations that are close to the identity:

````@example quantum_objects
haar_identity(2, 0.1)

haar_identity(2, 0.01)
````

## Embedded operators
Sometimes we want to embed a quantum operator into a larger Hilbert space, $\mathcal{H}$,
which we decompose into subspace and leakage components:
```math
    \mathcal{H} = \mathcal{H}_{\text{subspace}} \oplus \mathcal{H}_{\text{leakage}},
```
In quantum computing, the computation is encoded in a `subspace`, while the remaining
`leakage` states should be avoided.

### The `embed` and `unembed` functions

The [`embed`](@ref) function allows to embed a quantum operator in a larger Hilbert space.
```@docs
embed
```

The [`unembed`](@ref) function allows to unembed a quantum operator from a larger Hilbert space.
```@docs
unembed
```

For a single qubit X gate embedded in a multilevel system:

````@example quantum_objects
levels = 3
X = GATES[:X]
subspace_indices = 1:2
````

Embed the two-level X gate in the full system

````@example quantum_objects
X_embedded = embed(X, subspace_indices, levels)
````

We can retrieve the original operator:

````@example quantum_objects
X_original = unembed(X_embedded, subspace_indices)
````

### The `EmbeddedOperator` type
The [`EmbeddedOperator`](@ref) type stores information about an operator embedded in the subspace
of a larger quantum system.
```@docs
EmbeddedOperator
```

We can construct an embedded operator in the same manner as the `embed` function:
```@docs
EmbeddedOperator(subspace_operator::Matrix{<:Number}, subspace::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})
```

For an X gate on the first qubit of two qubit, 3-level system:

````@example quantum_objects
gate = GATES[:X] ⊗ GATES[:I]
subsystem_levels = [3, 3]
subspace_indices = get_subspace_indices([1:2, 1:2], subsystem_levels)
op = EmbeddedOperator(gate, subspace_indices, subsystem_levels)
````

Show the full operator.

````@example quantum_objects
op.operator .|> abs |> sparse
````

We can get the original operator back:

````@example quantum_objects
gate_unembeded = unembed(op)
gate_unembeded .|> abs |> sparse
````

## Subspace and leakage indices

### The `get_subspace_indices` function
The [`get_subspace_indices`](@ref) function is a convenient way to get the indices of a subspace in
a larger quantum system.
# ```@docs
# get_subspace_indices(subspace::AbstractVector{Int}, levels::Int)
# ```

Its dual function is [`get_leakage_indices`](@ref).

get the indices of the lowest two levels of a 3-level system

````@example quantum_objects
get_subspace_indices([1:2, 1:2], [3, 3])
````

qubits are assumed if the indices are not provided

````@example quantum_objects
get_subspace_indices([3, 3])
````

equivalently, get the indices of the leakage states

````@example quantum_objects
get_leakage_indices([3, 3])
````

#### Excitation number restrictions
Choose only the ground state and single excitation states of two 3-level systems:

````@example quantum_objects
get_enr_subspace_indices(1, [3, 3])
````

### The `get_iso_vec_subspace_indices` function
For isomorphic operators, the [`get_iso_vec_subspace_indices`](@ref) function can be used
to find the appropriate vector indices of the equivalent operator subspace.
```@docs
get_iso_vec_subspace_indices
```

Its dual function is [`get_iso_vec_leakage_indices`](@ref), which by default only returns
the leakage indices of the blocks:
```math
\mathcal{H}_{\text{subspace}} \otimes \mathcal{H}_{\text{subspace}}
\mathcal{H}_{\text{subspace}} \otimes \mathcal{H}_{\text{leakage}},
\mathcal{H}_{\text{leakage}} \otimes \mathcal{H}_{\text{leakage}},
```
allowing for leakage-suppressing code to disregard the uncoupled pure-leakage space.

````@example quantum_objects
get_iso_vec_subspace_indices(1:2, 3)

get_iso_vec_leakage_indices(1:2, 3)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

