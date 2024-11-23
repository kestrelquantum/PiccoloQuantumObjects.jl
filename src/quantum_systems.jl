module QuantumSystems

export AbstractQuantumSystem
export QuantumSystem

export get_drift
export get_drives

using ..Isomorphisms
using ..QuantumObjectUtils

using LinearAlgebra
using SparseArrays
using TestItemRunner
using ForwardDiff


# ----------------------------------------------------------------------------- #
# AbstractQuantumSystem
# ----------------------------------------------------------------------------- #

"""
    AbstractQuantumSystem

Abstract type for defining systems.
"""
abstract type AbstractQuantumSystem end

# ----------------------------------------------------------------------------- #
# AbstractQuantumSystem methods
# ----------------------------------------------------------------------------- #

"""
    get_drift(sys::AbstractQuantumSystem)
    
Returns the drift Hamiltonian of the system.
"""
get_drift(sys::AbstractQuantumSystem) = sys.H(zeros(sys.n_drives))

"""
    get_drives(sys::AbstractQuantumSystem)

Returns the drive Hamiltonians of the system.
"""
function get_drives(sys::AbstractQuantumSystem)
    H_drift = get_drift(sys)
    # Basis vectors for controls will extract drive operators
    return [sys.H(I[1:sys.n_drives, i]) - H_drift for i ∈ 1:sys.n_drives]
end


# ----------------------------------------------------------------------------- #
# QuantumSystem
# ----------------------------------------------------------------------------- #

"""
    QuantumSystem <: AbstractQuantumSystem

A struct for storing quantum dynamics and the appropriate gradients.

# Fields
- `H::Function`: The Hamiltonian function, excluding dissipation: a -> H(a).
- `G::Function`: The isomorphic generator function, including dissipation, a -> G(a).
- `∂G::Function`: The generator jacobian function, a -> ∂G(a).
- `levels::Int`: The number of levels in the system.
- `n_drives::Int`: The number of drives in the system.
"""
struct QuantumSystem <: AbstractQuantumSystem
    H::Function
    G::Function
    ∂G::Function
    n_drives::Int
    levels::Int
    params::Dict{Symbol, Any}
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
    H_drives::Vector{<:AbstractMatrix{<:Number}}; 
    params::Dict{Symbol, Any}=Dict{Symbol, Any}(),
)
    levels = size(H_drift, 1)
    H_drift = sparse(H_drift)
    G_drift = sparse(Isomorphisms.G(H_drift))
    
    n_drives = length(H_drives)
    H_drives = sparse.(H_drives)
    G_drives = sparse.(Isomorphisms.G.(H_drives))
    
    if n_drives == 0
        H = a -> H_drift
        G = a -> G_drift
        ∂G = a -> 0
    else
        H = a -> H_drift + sum(a .* H_drives)
        G = a -> G_drift + sum(a .* G_drives)
        ∂G = a -> G_drives
    end

    return QuantumSystem(
        H,
        G,
        ∂G,
        n_drives,
        levels,
        params
    )
end

function QuantumSystem(H_drives::Vector{<:AbstractMatrix{T}}; kwargs...) where T <: Number
    @assert !isempty(H_drives) "At least one drive is required"
    return QuantumSystem(spzeros(T, size(H_drives[1])), H_drives; kwargs...)
end

QuantumSystem(H_drift::AbstractMatrix{T}; kwargs...) where T <: Number = 
    QuantumSystem(H_drift, Matrix{T}[]; kwargs...)

function generator_jacobian(G::Function)
    return function ∂G(a::Vector{Float64})
        ∂G⃗ = ForwardDiff.jacobian(a_ -> vec(G(a_)), a)
        dim = Int(sqrt(size(∂G⃗, 1)))
        return [reshape(∂G⃗ⱼ, dim, dim) for ∂G⃗ⱼ ∈ eachcol(∂G⃗)]
    end
end

function QuantumSystem(H::Function, n_drives::Int; params=Dict{Symbol, Any}())
    G = a -> Isomorphisms.G(sparse(H(a)))
    ∂G = generator_jacobian(G)
    levels = size(H(zeros(n_drives)), 1)
    return QuantumSystem(H, G, ∂G, levels, n_drives, params)
end

function QuantumSystem(
    H_drift::AbstractMatrix,
    H_drives::Vector{<:AbstractMatrix},
    dissipation_operators::Vector{<:AbstractMatrix};
    params::Dict{Symbol, Any}=Dict{Symbol, Any}()
)
    levels = size(H_drift, 1)
    H_drift = sparse(H_drift)
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift))

    n_drives = length(H_drives)
    H_drives = sparse.(H_drives)
    𝒢_drives = Isomorphisms.G.(Isomorphisms.ad_vec.(H_drives))

    if isempty(dissipation_operators)
        𝒟 = zeros(size(𝒢_drift))
    else
        𝒟 = Isomorphisms.iso(sum(
            kron(conj(L), L) - 1 / 2 * Isomorphisms.ad_vec(L'L, anti=true) 
            for L ∈ sparse.(dissipation_operators)
        ))
    end

    if n_drives == 0
        H = a -> H_drift
        𝒢 = a -> 𝒢_drift + 𝒟
        ∂𝒢 = a -> 0
    else
        H = a -> H_drift + sum(a .* H_drives)
        𝒢 = a -> 𝒢_drift + sum(a .* 𝒢_drives) + 𝒟
        ∂𝒢 = a -> 𝒢_drives
    end

    return QuantumSystem(
        H,
        𝒢,
        ∂𝒢,
        levels,
        n_drives,
        params
    )

end

# ****************************************************************************** #

@testitem "System creation" begin
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    system = QuantumSystem(H_drift, H_drives)
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives

    # test jacobians
    a = randn(n_drives)
    ∂G = system.∂G(a)
    @test length(∂G) == system.n_drives
    @test all(∂G .≈ QuantumSystems.generator_jacobian(system.G)(a))

    # repeat with a bigger system
    H_drift = kron(GATES[:Z], GATES[:Z])
    H_drives = [kron(GATES[:X], GATES[:I]), kron(GATES[:I], GATES[:X]),
                kron(GATES[:Y], GATES[:I]), kron(GATES[:I], GATES[:Y])]
    n_drives = length(H_drives)

    system = QuantumSystem(H_drift, H_drives)
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives

    # test jacobians
    a = randn(n_drives)
    ∂G = system.∂G(a)
    @test length(∂G) == system.n_drives
    @test all(∂G .≈ QuantumSystems.generator_jacobian(system.G)(a))
end

@testitem "No drift system creation" begin
    H_drift = zeros(2, 2)
    H_drives = [GATES[:X], GATES[:Y]]

    sys1 = QuantumSystem(H_drift, H_drives)
    sys2 = QuantumSystem(H_drives)

    @test get_drift(sys1) == get_drift(sys2) == H_drift
    @test get_drives(sys1) == get_drives(sys2) == H_drives
end

@testitem "No drive system creation" begin
    H_drift = GATES[:Z]
    H_drives = Matrix{ComplexF64}[]

    sys1 = QuantumSystem(H_drift, H_drives)
    sys2 = QuantumSystem(H_drift)

    @test get_drift(sys1) == get_drift(sys2) == H_drift
    @test get_drives(sys1) == get_drives(sys2) == H_drives
end

@testitem "System creation with Hamiltonian function" begin
    H(a) = GATES[:Z] + a[1] * GATES[:X] + a[2] * GATES[:Y]
    system = QuantumSystem(H, 2)
    @test system isa QuantumSystem
    @test get_drift(system) == GATES[:Z]
    @test get_drives(system) == [GATES[:X], GATES[:Y]]

    # test jacobians
    compare = QuantumSystem(GATES[:Z], [GATES[:X], GATES[:Y]])
    a = randn(system.n_drives)
    @test system.∂G(a) == compare.∂G(a)
end

@testitem "System creation with dissipation" begin
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    dissipation_operators = [GATES[:Z], GATES[:X]]

    system = QuantumSystem(H_drift, H_drives, dissipation_operators)
    @test system isa QuantumSystem
    @test get_drift(system) == H_drift
    @test get_drives(system) == H_drives

    # test dissipation
    𝒢_drift = Isomorphisms.G(Isomorphisms.ad_vec(H_drift))
    @test system.G(zeros(system.n_drives)) != 𝒢_drift

    # test jacobians (disspiation is constant)
    a = randn(system.n_drives)
    ∂G = system.∂G(a)
    @test length(∂G) == system.n_drives
    @test all(∂G .≈ QuantumSystems.generator_jacobian(system.G)(a))
    
end

end
