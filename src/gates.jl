module Gates

export PAULIS
export GATES

using TestItemRunner


const PAULIS = (
    I = ComplexF64[1 0;
                    0 1],
    X = ComplexF64[0 1;
                    1 0],
    Y = ComplexF64[0 -im;
                    im 0],
    Z = ComplexF64[1 0;
                    0 -1],
)

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
GATES = (
    I = PAULIS.I,

    X = PAULIS.X,

    Y = PAULIS.Y,

    Z = PAULIS.Z,

    H = ComplexF64[1 1;
                    1 -1]/√2,

    CX = ComplexF64[1 0 0 0;
                    0 1 0 0;
                    0 0 0 1;
                    0 0 1 0],

    CZ = ComplexF64[1 0 0 0;
                    0 1 0 0;
                    0 0 1 0;
                    0 0 0 -1],

    XI = ComplexF64[0 0 -im 0;
                    0 0 0 -im;
                    -im 0 0 0;
                    0 -im 0 0],

    sqrtiSWAP = ComplexF64[1 0 0 0;
                            0 1/√2 1im/√2 0;
                            0 1im/√2 1/√2 0;
                            0 0 0 1],
)

# ******************************************************************************* #

@testitem "Gates are complex matrices" begin
    # test call
    @test GATES.X == [0 1; 1 0]
    @test GATES[:X] == [0 1; 1 0]

    for k in keys(GATES)
        @test typeof(GATES[k]) == Matrix{ComplexF64}
    end
end

@testitem "Paulis are complex matrices" begin
    # test call
    @test PAULIS.X == [0 1; 1 0]
    @test PAULIS[:X] == [0 1; 1 0]
    
    for k in keys(PAULIS)
        @test typeof(PAULIS[k]) == Matrix{ComplexF64}
    end
end

end
