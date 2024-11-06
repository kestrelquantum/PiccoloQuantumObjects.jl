module Gates

export GATES
export PAULIS

@doc raw"""
A constant dictionary `GATES` containing common quantum gate matrices as complex-valued matrices. Each gate is represented by its unitary matrix.

- `GATES[:I]` - Identity gate: Leaves the state unchanged.
- `GATES[:X]` - Pauli-X (NOT) gate: Flips the qubit state.
- `GATES[:Y]` - Pauli-Y gate: Rotates the qubit state around the Y-axis of the Bloch sphere.
- `GATES[:Z]` - Pauli-Z gate: Flips the phase of the qubit state.
- `GATES[:H]` - Hadamard gate: Creates superposition by transforming basis states.
- `GATES[:CX]` - Controlled-X (CNOT) gate: Flips the second qubit (target) if the first qubit (control) is |1⟩.
- `GATES[:CZ]` - Controlled-Z (CZ) gate: Flips the phase of the second qubit (target) if the first qubit (control) is |1⟩.
- `GATES[:XI]` - Complex gate: A specific gate used for complex operations.
- `GATES[:sqrtiSWAP]` - Square root of iSWAP gate: Partially swaps two qubits with a phase.

```julia
julia> GATES[:Z]
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im

julia> get_gate(:CX)
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
```
"""
const GATES = Dict{Symbol, Matrix{ComplexF64}}(
    :I => [1 0;
           0 1],

    :X => [0 1;
           1 0],

    :Y => [0 -im;
           im 0],

    :Z => [1 0;
           0 -1],

    :H => [1 1;
           1 -1]/√2,

    :CX => [1 0 0 0;
            0 1 0 0;
            0 0 0 1;
            0 0 1 0],

    :CZ => [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 -1],

    :XI => [0 0 -im 0;
            0 0 0 -im;
            -im 0 0 0;
            0 -im 0 0],

    :sqrtiSWAP => [1 0 0 0;
                   0 1/sqrt(2) 1im/sqrt(2) 0;
                   0 1im/sqrt(2) 1/sqrt(2) 0;
                   0 0 0 1]
)

const PAULIS = Dict{Symbol, Matrix{ComplexF64}}(
    :I => GATES[:I],
    :X => GATES[:X],
    :Y => GATES[:Y],
    :Z => GATES[:Z]
)

end
