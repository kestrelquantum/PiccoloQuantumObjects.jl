module QuantumSystems

export AbstractQuantumSystem
export QuantumSystem

export get_H_drift
export get_H_drives

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
# QuantumSystem
# ----------------------------------------------------------------------------- #

"""
    QuantumSystem <: AbstractQuantumSystem

A struct for storing the isomorphisms of the system's drift and drive Hamiltonians,
as well as the system's parameters.

# Fields
- `H::Function`: The Hamiltonian function, excluding dissipation: a -> H(a).
- `G::Function`: The isomorphic generator function, including dissipation, a -> G(a).
- `∂G::Function`: The generator jacobian function, a -> ∂G(a).
- `levels::Int`: The number of levels in the system.
- `n_drives::Int`: The number of drives in the system.

# Examples

Construct a `QuantumSystem` object from the drift and drive Hamiltonian terms:

```jldoctest
julia> H_drift = GATES[:Z]
julia> H_drives = [GATES[:X], GATES[:Y]]
julia> system = QuantumSystem(H_drift, H_drives)

julia> a = randn(system.n_drives)
```

Construct a `QuantumSystem` with Lindblad operators:

```jldoctest
julia> H_drift = GATES[:Z]
julia> H_drives = [GATES[:X], GATES[:Y]]
julia> dissipation_operators = [GATES[:Z], annihilate(2)]
julia> system = QuantumSystem(H_drift, H_drives, dissipation_operators)
```

Construct a `QuantumSystem` from a Hamiltonian function and ForwardDiff.jl:
```jldoctest

julia> H(a) = GATES[:Z] + a[1] * GATES[:X] + a[2] * GATES[:Y]
julia> system = QuantumSystem(H, 2)
```
"""
struct QuantumSystem <: AbstractQuantumSystem
    H::Function
    G::Function
    ∂G::Function
    levels::Int
    n_drives::Int
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
    params=Dict{Symbol, Any}(),
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
        H_drives = sparse.(H_drives)
        G_drives = sparse.(Isomorphisms.G.(H_drives))

        H = a -> H_drift + sum(a .* H_drives)
        G = a -> G_drift + sum(a .* G_drives)
        ∂G = a -> G_drives
    end

    return QuantumSystem(
        H,
        G,
        ∂G,
        levels,
        n_drives,
        params
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
        Matrix{eltype(H_drift)}[];
        kwargs...
    )
end

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
    params=Dict{Symbol, Any}()
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

# ----------------------------------------------------------------------------- #
# QuantumSystem methods
# ----------------------------------------------------------------------------- #

get_H_drift(sys::QuantumSystem) = sys.H(zeros(sys.n_drives))

function get_H_drives(sys::QuantumSystem)
    # Basis vectors for controls will extract drive operators
    H_drift = get_H_drift(sys)
    return [sys.H(I[1:sys.n_drives, i]) - H_drift for i ∈ 1:sys.n_drives]
end

# ****************************************************************************** #

@testitem "System creation" begin
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    system = QuantumSystem(H_drift, H_drives)
    @test system isa QuantumSystem
    @test get_H_drift(system) == H_drift
    @test get_H_drives(system) == H_drives

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

    @test get_H_drift(sys1) == get_H_drift(sys2) == H_drift
    @test get_H_drives(sys1) == get_H_drives(sys2) == H_drives
end

@testitem "No drive system creation" begin
    H_drift = GATES[:Z]
    H_drives = Matrix{ComplexF64}[]

    sys1 = QuantumSystem(H_drift, H_drives)
    sys2 = QuantumSystem(H_drift)

    @test get_H_drift(sys1) == get_H_drift(sys2) == H_drift
    @test get_H_drives(sys1) == get_H_drives(sys2) == H_drives
end

@testitem "System creation with Hamiltonian function" begin
    H(a) = GATES[:Z] + a[1] * GATES[:X] + a[2] * GATES[:Y]
    system = QuantumSystem(H, 2)
    @test system isa QuantumSystem
    @test get_H_drift(system) == GATES[:Z]
    @test get_H_drives(system) == [GATES[:X], GATES[:Y]]

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
    @test get_H_drift(system) == H_drift
    @test get_H_drives(system) == H_drives

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
