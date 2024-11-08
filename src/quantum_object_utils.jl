module QuantumObjectUtils

export operator_from_string
export ket_from_string
export ket_from_bitstring

export haar_random
export haar_identity

export create
export annihilate

using ..Gates

using LinearAlgebra
using TestItemRunner

# TODO:
# [ ] Remove need for oscillator operators (used by tests)
# [ ] Allow multi-character symbols for operators_from_string
# [ ] Remove need for otimes symbol or avoid import conflicts with other packages


@doc raw"""
operator_from_string(operator::String; lookup::Dict{Symbol, AbstractMatrix}=PAULIS)

    Reduce the string (each character is one key) via operators from a dictionary.

"""
function operator_from_string(
    operator::String;
    lookup::Dict{Symbol, <:AbstractMatrix}=PAULIS
)::Matrix{ComplexF64}
    # TODO: allow multi-character keys, ['(', ')']

    # split string into keys and replace with operators
    characters = [Symbol(c) for c ∈ operator]
    operators = replace(characters, lookup...)

    return foldr(kron, operators)
end

function cavity_state(state::Int, levels::Int)::Vector{ComplexF64}
    @assert state ≤ levels - 1 "Level $state is not allowed for $levels levels"
    ket = zeros(levels)
    ket[state + 1] = 1
    return ket
end

@doc raw"""
    ket_from_string(
        ket::String,
        levels::Vector{Int};
        level_dict=Dict(:g => 0, :e => 1, :f => 2, :h => 2),
        return_states=false
    )

Construct a quantum state from a string ket representation.

# Example

# TODO: add example
"""
function ket_from_string(
    ket::String,
    levels::Vector{Int};
    level_dict=Dict(:g => 0, :e => 1, :f => 2, :h => 2),
    return_states=false
)::Vector{ComplexF64}
    kets = []

    for x ∈ split(ket, ['(', ')'])
        if x == ""
            continue
        elseif all(Symbol(xᵢ) ∈ keys(level_dict) for xᵢ ∈ x)
            append!(kets, x)
        elseif occursin("+", x)
            superposition = split(x, '+')
            @assert all(all(Symbol(xᵢ) ∈ keys(level_dict) for xᵢ ∈ x) for x ∈ superposition) "Invalid ket: $x"
            @assert length(superposition) == 2 "Only two states can be superposed for now"
            push!(kets, x)
        else
            error("Invalid ket: $x")
        end
    end

    states = []

    for (ψᵢ, l) ∈ zip(kets, levels)
        if ψᵢ isa AbstractString && occursin("+", ψᵢ)
            superposition = split(ψᵢ, '+')
            superposition_states = [level_dict[Symbol(x)] for x ∈ superposition]
            @assert all(state ≤ l - 1 for state ∈ superposition_states) "Level $ψᵢ is not allowed for $l levels"
            superposition_state = sum([
                cavity_state(state, l) for state ∈ superposition_states
            ])
            normalize!(superposition_state)
            push!(states, superposition_state)
        else
            state = level_dict[Symbol(ψᵢ)]
            @assert state ≤ l - 1 "Level $ψᵢ is not allowed for $l levels"
            push!(states, cavity_state(state, l))
        end
    end

    if return_states
        return states
    else
        return kron([1.0], states...)
    end
end

@doc raw"""
    ket_from_bitstring(ket::String)

Get the state vector for a qubit system given a ket string `ket` of 0s and 1s.
"""
function ket_from_bitstring(ket::String)::Vector{ComplexF64}
    cs = [c for c ∈ ket]
    @assert all(c ∈ "01" for c ∈ cs)
    states = [c == '0' ? [1, 0] : [0, 1] for c ∈ cs]
    return foldr(kron, states)
end

###
### Random operators
###

@doc raw"""
    haar_random(n::Int)

Generate a random unitary matrix using the Haar measure for an `n`-dimensional system.
"""
function haar_random(n::Int)
    # Ginibre matrix
    Z = (randn(n, n) + im * randn(n, n)) / √2
    F = qr(Z)
    # QR correction (R main diagonal is real, strictly positive)
    Λ = diagm(diag(F.R) ./ abs.(diag(F.R)))
    return F.Q * Λ
end

@doc raw"""
    haar_identity(n::Int, radius::Number)

Generate a random unitary matrix close to the identity matrix using the Haar measure for an `n`-dimensional system with a given `radius`.
"""
function haar_identity(n::Int, radius::Number)
    # Ginibre matrix
    Z = (I + radius * (randn(n, n) + im * randn(n, n)) / √2) / (1 + radius)
    F = qr(Z)
    # QR correction (R main diagonal is real, strictly positive)
    Λ = diagm(diag(F.R) ./ abs.(diag(F.R)))
    return F.Q * Λ
end

###
### Oscillator operators
###

@doc raw"""
    annihilate(levels::Int)

Get the annihilation operator for a system with `levels` levels.
"""
function annihilate(levels::Int)::Matrix{ComplexF64}
    return diagm(1 => map(sqrt, 1:levels - 1))
end

@doc raw"""
    create(levels::Int)

Get the creation operator for a system with `levels` levels.
"""
function create(levels::Int)
    return collect(annihilate(levels)')
end

# ============================================================================= #

@testitem "Test ket_from_bitstring function" begin
    using LinearAlgebra
    @test ket_from_bitstring("0") == [1, 0]
    @test ket_from_bitstring("1") == [0, 1]
    @test ket_from_bitstring("00") == [1, 0, 0, 0]
    @test ket_from_bitstring("01") == [0, 1, 0, 0]
    @test ket_from_bitstring("10") == [0, 0, 1, 0]
    @test ket_from_bitstring("11") == [0, 0, 0, 1]
end


end
