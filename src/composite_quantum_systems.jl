module CompositeQuantumSystems

export lift
export CompositeQuantumSystem

using ..QuantumSystems
using ..Isomorphisms

using LinearAlgebra
using SparseArrays
using TestItemRunner

# ----------------------------------------------------------------------------- #
# Lift operators
# ----------------------------------------------------------------------------- #

@doc raw"""
    lift(operator::AbstractMatrix{<:Number}, i::Int, subsystem_levels::Vector{Int})
    lift(operator::AbstractMatrix{<:Number}, i::Int, n_qubits::Int; kwargs...)
    lift(operators::AbstractVector{<:AbstractMatrix{T}}, indices::AbstractVector{Int}, subsystem_levels::Vector{Int})
    lift(operators::AbstractVector{<:AbstractMatrix{T}}, indices::AbstractVector{Int}, n_qubits::Int; kwargs...)

Lift an `operator` acting on the `i`-th subsystem within `subsystem_levels` to an operator
acting on the entire system spanning `subsystem_levels`.
"""
function lift end

function lift(operator::AbstractMatrix{T}, i::Int, subsystem_levels::Vector{Int}
) where T <: Number
    @assert size(operator, 1) == subsystem_levels[i] "Operator must match subsystem level."
    Is = [Matrix{T}(I(l)) for l ∈ subsystem_levels]
    Is[i] = operator
    return reduce(kron, Is)
end

function lift(
    operator::AbstractMatrix{T}, i::Int, n_qubits::Int; 
    levels::Int=size(operator, 1)
) where T <: Number
    return lift(operator, i, fill(levels, n_qubits))
end

function lift(
    operators::AbstractVector{<:AbstractMatrix{T}},
    indices::AbstractVector{Int},
    subsystem_levels::Vector{Int}
) where T <: Number
    @assert length(operators) == length(indices)
    return prod([lift(op, i, subsystem_levels) for (op, i) ∈ zip(operators, indices)])
end

function lift(
    operators::AbstractVector{<:AbstractMatrix{T}},
    indices::AbstractVector{Int},
    n_qubits::Int;
    levels::Int=size(operators[1], 1)
) where T <: Number
    return prod(
        [lift(op, i, n_qubits, levels=levels) for (op, i) ∈ zip(operators, indices)]
    )
end

# ----------------------------------------------------------------------------- #
# Composite Quantum Systems
# ----------------------------------------------------------------------------- #

"""
    CompositeQuantumSystem <: AbstractQuantumSystem

A composite quantum system consisting of `subsystems`. Couplings between subsystems can
be additionally defined. Subsystem drives are always appended to any new coupling drives.
"""
struct CompositeQuantumSystem <: AbstractQuantumSystem
    H::Function
    G::Function
    ∂G::Function
    n_drives::Int
    levels::Int
    params::Dict{Symbol, Any}
    subsystem_levels::Vector{Int}
    subsystems::Vector{QuantumSystem}
end

function CompositeQuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    H_drives::AbstractVector{<:AbstractMatrix{<:Number}},
    subsystems::AbstractVector{QuantumSystem};
    params::Dict{Symbol, Any}=Dict{Symbol, Any}()
)
    subsystem_levels = [sys.levels for sys ∈ subsystems]
    levels = prod(subsystem_levels)

    H_drift = sparse(H_drift)
    for (i, sys) ∈ enumerate(subsystems)
        H_drift += lift(get_H_drift(sys), i, subsystem_levels)
    end

    H_drives = sparse.(H_drives)
    for (i, sys) ∈ enumerate(subsystems)
        for H_drive ∈ get_H_drives(sys)
            push!(H_drives, lift(H_drive, i, subsystem_levels))
        end
    end

    n_drives = length(H_drives)
    H_drives = sparse.(H_drives)
    G_drives = sparse.(Isomorphisms.G.(H_drives))

    if n_drives == 0
        H = a -> H_drift
        G = a -> Isomorphisms.G(H_drift)
        ∂G = a -> 0
    else
        H = a -> H_drift + sum(a .* H_drives)
        G = a -> G_drift + sum(a .* G_drives)
        ∂G = a -> G_drives
    end

    return CompositeQuantumSystem(
        H,
        G,
        ∂G,
        n_drives,
        levels,
        params,
        subsystem_levels,
        subsystems
    )
end

function CompositeQuantumSystem(
    H_drives::AbstractVector{<:AbstractMatrix{T}},
    subsystems::AbstractVector{QuantumSystem};
    kwargs...
) where T <: Number
    @assert !isempty(H_drives) "At least one drive is required"
    return CompositeQuantumSystem(
        spzeros(T, size(H_drives[1])),
        H_drives,
        subsystems;
        kwargs...
    )
end

function CompositeQuantumSystem(
    H_drift::AbstractMatrix{T},
    subsystems::AbstractVector{QuantumSystem};
    kwargs...
) where T <: Number
    return CompositeQuantumSystem(H_drift, Matrix{T}[], subsystems; kwargs...)
end

function CompositeQuantumSystem(
    subsystems::AbstractVector{QuantumSystem};
    kwargs...
)
    @assert !isempty(subsystems) "At least one subsystem is required"
    T = eltype(get_H_drift(subsystems[1]))
    levels = prod([sys.levels for sys ∈ subsystems])
    return CompositeQuantumSystem(
        spzeros(T, (levels, levels)), Matrix{T}[], subsystems; kwargs...
    )
end

# ****************************************************************************** #

@testitem "Lift subsystems" begin
    using LinearAlgebra
    @test lift(PAULIS[:X], 1, [2, 3]) ≈ kron(PAULIS[:X], I(3))
    @test lift(PAULIS[:Y], 2, [4, 2]) ≈ kron(I(4), PAULIS[:Y])
    @test lift(PAULIS[:X], 2, [3, 2, 4]) ≈ reduce(kron, [I(3), PAULIS[:X], I(4)])
end

@testitem "Lift qubits" begin
    using LinearAlgebra
    @test lift(PAULIS[:X], 1, 2) ≈ kron(PAULIS[:X], I(2))
    @test lift(PAULIS[:Y], 2, 2) ≈ kron(I(2), PAULIS[:Y])
    @test lift(PAULIS[:X], 2, 3) ≈ reduce(kron, [I(2), PAULIS[:X], I(2)])
end

@testitem "Lift multiple operators" begin
    using LinearAlgebra
    pair = [PAULIS[:X], PAULIS[:Y]]
    @test lift(pair, [1, 2], [2, 2]) ≈ kron(PAULIS[:X], PAULIS[:Y])
    @test lift(pair, [2, 1], [2, 2]) ≈ kron(PAULIS[:Y], PAULIS[:X])
    @test lift(pair, [1, 2], [2, 2, 3]) ≈ kron(PAULIS[:X], PAULIS[:Y], I(3))
    @test lift(pair, [2, 3], [4, 2, 2]) ≈ kron(I(4), PAULIS[:X], PAULIS[:Y])

    # Pass number of qubits
    @test lift(pair, [1, 2], 3) ≈ kron(PAULIS[:X], PAULIS[:Y], I(2))
end


@testitem "Composite system" begin
    subsystem_levels = [4, 2, 2]
    sys1 = QuantumSystem(kron(PAULIS[:Z], PAULIS[:Z]), [kron(PAULIS[:X], PAULIS[:Y])])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]])
    sys3 = QuantumSystem(zeros(2, 2))
    subsystems = [sys1, sys2, sys3]
    g12 = 0.1 * lift([kron(PAULIS[:X], PAULIS[:X]), PAULIS[:X]], [1, 2], subsystem_levels)
    g23 = 0.2 * lift([PAULIS[:Y], PAULIS[:Y]], [2, 3], subsystem_levels)

    # Construct composite system
    csys = CompositeQuantumSystem(g12, [g23], [sys1, sys2, sys3])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == 1 + sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_H_drift(csys) ≈ g12 + lift(kron(PAULIS[:Z], PAULIS[:Z]), 1, subsystem_levels)end

@testitem "Composite system from drift" begin
    using LinearAlgebra

    subsystem_levels = [2, 2]
    sys1 = QuantumSystem([PAULIS[:X], PAULIS[:Y]])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]])
    subsystems = [sys1, sys2]
    g12 = 0.1 * kron(PAULIS[:X], PAULIS[:X])

    # Construct composite system from drift
    csys = CompositeQuantumSystem(g12, [sys1, sys2])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_H_drift(csys) ≈ g12
end

@testitem "Composite system from drives" begin
    subsystem_levels = [2, 2, 2]
    sys1 = QuantumSystem(PAULIS[:Z], [PAULIS[:X], PAULIS[:Y]])
    sys2 = QuantumSystem([PAULIS[:Y], PAULIS[:Z]])
    sys3 = QuantumSystem(zeros(2, 2))
    subsystems = [sys1, sys2, sys3]
    g12 = 0.1 * lift([PAULIS[:X], PAULIS[:X]], [1, 2], subsystem_levels)
    g23 = 0.2 * lift([PAULIS[:Y], PAULIS[:Y]], [2, 3], subsystem_levels)

    csys = CompositeQuantumSystem([g12, g23], [sys1, sys2, sys3])
    @test csys.levels == prod(subsystem_levels)
    @test csys.n_drives == 2 + sum([sys.n_drives for sys ∈ subsystems])
    @test csys.subsystems == subsystems
    @test csys.subsystem_levels == subsystem_levels
    @test get_H_drift(csys) ≈ lift(PAULIS[:Z], 1, subsystem_levels)
end

end
