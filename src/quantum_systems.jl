module QuantumSystems

export AbstractQuantumSystem
export QuantumSystem
export CompositeQuantumSystem
export QuantumSystemCoupling

export lift

using ..Isomorphisms

using LinearAlgebra
using SparseArrays
using TestItemRunner
using ForwardDiff


# TODO:
# [ ] subtypes? SingleQubitSystem, TwoQubitSystem, TransmonSystem, MultimodeSystem, etc.
# [ ] add frame info to type
# [ ] add methods to combine composite quantum systems

# ----------------------------------------------------------------------------- #
# AbstractQuantumSystem
# ----------------------------------------------------------------------------- #

"""
    AbstractQuantumSystem

Abstract type for defining systems.
"""
abstract type AbstractQuantumSystem end

# ----------------------------------------------------------------------------- #
# QuantumSystem
# ----------------------------------------------------------------------------- #

"""
    QuantumSystem <: AbstractQuantumSystem

A struct for storing the isomorphisms of the system's drift and drive Hamiltonians,
as well as the system's parameters.
"""
struct QuantumSystem <: AbstractQuantumSystem
    H::Function
    G::Function
    ∂G::Function
    levels::Int
    n_drives::Int
end

"""
    QuantumSystem(
        H_drift::Matrix{<:Number},
        H_drives::Vector{Matrix{<:Number}};
        params=Dict{Symbol, Any}(),
        kwargs...
    )::QuantumSystem

Constructs a `QuantumSystem` object from the drift and drive Hamiltonian terms.
"""
function QuantumSystem(
    H_drift::AbstractMatrix{<:Number},
    H_drives::Vector{<:AbstractMatrix{<:Number}}
)
    H_drift = sparse(H_drift)
    H_drives = sparse.(H_drives)
    G_drift = sparse(Isomorphisms.G(H_drift))
    G_drives = sparse.(Isomorphisms.G.(H_drives))
    H = a -> H_drift + sum(a .* H_drives)
    G = a -> G_drift + sum(a .* G_drives)
    ∂G = a -> G_drives
    levels = size(H_drift, 1)
    return QuantumSystem(
        H,
        G,
        ∂G,
        levels,
        length(H_drives)
    )
end

function QuantumSystem(H_drives::Vector{<:AbstractMatrix{<:Number}}; kwargs...)
    return QuantumSystem(
        spzeros(eltype(H_drives[1]), size(H_drives[1])),
        H_drives;
        kwargs...
    )
end

function QuantumSystem(H_drift::AbstractMatrix{<:Number}; kwargs...)
    return QuantumSystem(
        H_drift,
        Matrix{ComplexF64}[];
        kwargs...
    )
end

function generator_jacobian(H::Function)
    return function ∂G(a::Vector{Float64})
        ∂G⃗ = ForwardDiff.jacobian(a_ -> vec(sparse(H(a_))), a)
        dim = Int(sqrt(size(∂G⃗, 1)))
        return [reshape(∂G⃗ⱼ, dim, dim) for ∂G⃗ⱼ ∈ eachcol(∂G⃗)]
    end
end

function QuantumSystem(H::Function, n_drives::Int)
    G = a -> Isomorphisms.G(sparse(H(a)))
    ∂G = generator_jacobian(H)
    levels = size(H(zeros(n_drives)), 1)
    return QuantumSystem(H, G, ∂G, levels, n_drives)
end



function Base.copy(sys::QuantumSystem)
    return QuantumSystem(
        copy(sys.H_drift),
        copy.(sys.H_drives)
    )
end



# ============================================================================= #

@testitem "System creation" begin
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    system = QuantumSystem(H_drift, H_drives)
end


end
