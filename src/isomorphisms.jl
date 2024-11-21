module Isomorphisms

# Export state isomorphisms
export mat
export ket_to_iso
export iso_to_ket
export dm_to_iso
export iso_to_dm
export iso_vec_to_operator
export iso_vec_to_iso_operator
export operator_to_iso_vec
export iso_operator_to_iso_vec
export iso_operator_to_operator
export operator_to_iso_operator

# Do not export Hamiltonian isomorphisms

using LinearAlgebra
using SparseArrays
using TestItemRunner


@doc raw"""
    mat(x::AbstractVector)

Convert a vector `x` into a square matrix. The length of `x` must be a perfect square.
"""
function mat(x::AbstractVector)
    n = isqrt(length(x))
    @assert n^2 ≈ length(x) "Vector length must be a perfect square"
    return reshape(x, n, n)
end


# ----------------------------------------------------------------------------- #
#                                Kets                                           #
# ----------------------------------------------------------------------------- #

@doc raw"""
    ket_to_iso(ψ)

Convert a ket vector `ψ` into a complex vector with real and imaginary parts.
"""
ket_to_iso(ψ) = [real(ψ); imag(ψ)]

@doc raw"""
    iso_to_ket(ψ̃)

Convert a complex vector `ψ̃` with real and imaginary parts into a ket vector.
"""
iso_to_ket(ψ̃) = ψ̃[1:div(length(ψ̃), 2)] + im * ψ̃[(div(length(ψ̃), 2) + 1):end]

# ----------------------------------------------------------------------------- #
#                             Unitaries                                         #
# ----------------------------------------------------------------------------- #

@doc raw"""
    iso_vec_to_operator(Ũ⃗::AbstractVector)

Convert a real vector `Ũ⃗` into a complex matrix representing an operator.

Must be differentiable.
"""
function iso_vec_to_operator(Ũ⃗::AbstractVector{R}) where R <: Real
    Ũ⃗_dim = div(length(Ũ⃗), 2)
    N = Int(sqrt(Ũ⃗_dim))
    U = Matrix{complex(R)}(undef, N, N)
    for i=0:N-1
        U[:, i+1] .= @view(Ũ⃗[i * 2N .+ (1:N)]) + one(R) * im * @view(Ũ⃗[i * 2N .+ (N+1:2N)])
    end
    return U
end

@doc raw"""
    iso_vec_to_iso_operator(Ũ⃗::AbstractVector)

Convert a real vector `Ũ⃗` into a real matrix representing an isomorphism operator.

Must be differentiable.
"""
function iso_vec_to_iso_operator(Ũ⃗::AbstractVector{R}) where R <: Real
    N = Int(sqrt(length(Ũ⃗) ÷ 2))
    Ũ = Matrix{R}(undef, 2N, 2N)
    U_real = Matrix{R}(undef, N, N)
    U_imag = Matrix{R}(undef, N, N)
    for i=0:N-1
        U_real[:, i+1] .= @view(Ũ⃗[i*2N .+ (1:N)])
        U_imag[:, i+1] .= @view(Ũ⃗[i*2N .+ (N+1:2N)])
    end
    Ũ[1:N, 1:N] .= U_real
    Ũ[1:N, (N + 1):end] .= -U_imag
    Ũ[(N + 1):end, 1:N] .= U_imag
    Ũ[(N + 1):end, (N + 1):end] .= U_real
    return Ũ
end

@doc raw"""
    operator_to_iso_vec(U::AbstractMatrix{<:Complex})

Convert a complex matrix `U` representing an operator into a real vector.

Must be differentiable.

# Example 

```jldoctest
julia> U = [1 5; 2 6] + im * [3 7; 4 8]
2×2 Matrix{Complex{Int64}}:
 1+3im  5+7im
 2+4im  6+8im

julia> operator_to_iso_vec(U)
8-element Vector{Int64}:
    1
    2
    3
    4
    5
    6
    7
    8

```
"""
function operator_to_iso_vec(U::AbstractMatrix{R}) where R <: Number
    N = size(U,1)
    Ũ⃗ = Vector{real(R)}(undef, N^2 * 2)
    for i=0:N-1
        Ũ⃗[i*2N .+ (1:N)] .= real(@view(U[:, i+1]))
        Ũ⃗[i*2N .+ (N+1:2N)] .= imag(@view(U[:, i+1]))
    end
    return Ũ⃗
end

@doc raw"""
    iso_operator_to_iso_vec(Ũ::AbstractMatrix)

Convert a real matrix `Ũ` representing an isomorphism operator into a real vector.

Must be differentiable.
"""
function iso_operator_to_iso_vec(Ũ::AbstractMatrix{R}) where R <: Real
    N = size(Ũ, 1) ÷ 2
    Ũ⃗ = Vector{R}(undef, N^2 * 2)
    for i=0:N-1
        Ũ⃗[i*2N .+ (1:2N)] .= @view Ũ[:, i+1]
    end
    return Ũ⃗
end

iso_operator_to_operator(Ũ) = iso_vec_to_operator(iso_operator_to_iso_vec(Ũ))

operator_to_iso_operator(U) = iso_vec_to_iso_operator(operator_to_iso_vec(U))

# ----------------------------------------------------------------------------- #
#                             Density matrix                                    #
# ----------------------------------------------------------------------------- #

@doc raw"""
    dm_to_iso(ρ::AbstractMatrix)

Returns the isomorphism `ρ⃗̃ = ket_to_iso(vec(ρ))` of a density matrix `ρ`
"""
dm_to_iso(ρ::AbstractMatrix{<:Number}) = ket_to_iso(vec(ρ))

@doc raw"""
    iso_to_dm(ρ⃗̃::AbstractVector)

Returns the density matrix `ρ` from its isomorphism `ρ⃗̃`
"""
iso_to_dm(ρ⃗̃::AbstractVector) = mat(iso_to_ket(ρ⃗̃))

# ----------------------------------------------------------------------------- #
#                             Hamiltonians                                      #
# ----------------------------------------------------------------------------- #

const Im2 = [
    0 -1;
    1  0
]

@doc raw"""
   iso(H::AbstractMatrix{<:Number})

Returns the isomorphism of ``H``:

```math
iso(H) = \widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

where ``\Im(H)`` and ``\Re(H)`` are the imaginary and real parts of ``H`` and the tilde 
indicates the standard isomorphism of a complex valued matrix:

```math
\widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

Note that ``iso(-iH) = G(H)``.
"""
iso(H::AbstractMatrix{<:Number}) = kron(I(2), real(H)) + kron(Im2, imag(H))

@doc raw"""
    G(H::AbstractMatrix)::Matrix{Float64}

Returns the isomorphism of ``-iH``:

```math
G(H) = \widetilde{- i H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Im(H) - \mqty(0 & -1 \\ 1 & 0) \otimes \Re(H)
```

where ``\Im(H)`` and ``\Re(H)`` are the imaginary and real parts of ``H`` and the tilde 
indicates the standard isomorphism of a complex valued matrix:

```math
\widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

Note that ``G(H) = iso(-iH)``.
"""
G(H::AbstractMatrix{<:Number}) = iso(-im * H)


@doc raw"""
    H(G::AbstractMatrix{<:Number})::Matrix{ComplexF64}

Returns the inverse of ``G(H) = iso(-iH)``, i.e. returns H

"""
function H(G::AbstractMatrix{<:Number})
    dim = size(G, 1) ÷ 2
    H_imag = G[1:dim, 1:dim]
    H_real = -G[dim+1:end, 1:dim]
    return H_real + 1.0im * H_imag
end


@doc raw"""
    ad_vec(H::AbstractMatrix{R}; anti::Bool=false) where R <: Number

Returns the vectorized adjoint action of a matrix `H`:

```math
\text{ad_vec}(H) = \mqty(1 & 0 \\ 0 & 1) \otimes H - (-1)^{\text{anti}} \mqty(0 & 1 \\ 1 & 0) \otimes H^*
```
"""
function ad_vec(H::AbstractMatrix{R}; anti::Bool=false) where R <: Number
    Id = sparse(R, I, size(H)...)
    return kron(Id, H) - (-1)^anti * kron(conj(H)', Id)
end


# *************************************************************************** #

@testitem "Test ket isomorphisms" begin
    @test ket_to_iso([1.0, 2.0]) ≈ [1.0, 2.0, 0.0, 0.0]
    @test ket_to_iso([-im, 2.0 + 3.0im]) ≈ [0.0, 2.0, -1.0, 3.0]
    @test iso_to_ket([1.0, 2.0, 0.0, 0.0]) ≈ [1.0, 2.0]
    @test iso_to_ket([0.0, 2.0, -1.0, 3.0]) ≈ [-im, 2.0 + 3.0im]
end

@testitem "Test operator isomorphisms" begin
    iso_vec_I = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    @test mat([1.0, 2.0, 3.0, 4.0]) ≈ [1.0 3.0; 
                                        2.0 4.0]
    @test iso_vec_to_operator(iso_vec_I) ≈ [1.0 0.0;
                                           0.0 1.0]
    @test iso_vec_to_iso_operator(iso_vec_I) ≈ [1.0 0.0 0.0 0.0; 
                                               0.0 1.0 0.0 0.0; 
                                               0.0 0.0 1.0 0.0; 
                                               0.0 0.0 0.0 1.0]
    @test operator_to_iso_vec(Complex[1.0 0.0; 0.0 1.0]) ≈ iso_vec_I
    @test iso_operator_to_iso_vec(iso_vec_to_iso_operator(iso_vec_I)) ≈ iso_vec_I

    iso_vec_XY = [0, 1, 0, 1, 1, 0, -1, 0]
    @test iso_vec_to_operator(iso_vec_XY) ≈ [0 1-im; 
                                             1+im 0]
    @test iso_vec_to_iso_operator(iso_vec_XY) ≈ [0 1 0 1; 
                                                 1 0 -1 0; 
                                                 0 -1 0 1; 
                                                 1 0 1 0]
    @test operator_to_iso_vec(Complex[0.0 1-im; 1+im 0.0]) ≈ iso_vec_XY
    @test iso_operator_to_iso_vec(iso_vec_to_iso_operator(iso_vec_XY)) ≈ iso_vec_XY
end

@testitem "Test density matrix isomorphisms" begin
    # Totally mixed state
    ρ = [1.0 0.0; 0.0 1.0]
    ρ_iso = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    @test dm_to_iso(ρ) ≈ ρ_iso
    @test iso_to_dm(ρ_iso) ≈ ρ

    # Density matrix of a Bell state
    ρ = [1.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 1.0] / 2
    @test iso_to_dm(dm_to_iso(ρ)) ≈ ρ

    # Random
    ρ1 = [1.0 1.0; 1.0 1.0] / 2
    U1 = [-0.831976-0.101652im  -0.422559-0.344857im;
          -0.527557+0.138444im   0.799158+0.252713im]
    ρ2 = [1.0 0.0; 0.0 0.0]
    U2 = [-0.784966-0.163279im   -0.597246-0.0215881im
           0.597536+0.0109124im  -0.792681+0.120364im]
    ρ = (U1*ρ1*U1' + U2*ρ2*U2') / 2
    @test iso_to_dm(dm_to_iso(ρ)) ≈ ρ
    @test iso_to_dm(dm_to_iso(ρ)) ≈ ρ
end

@testitem "Test Hamiltonian isomorphisms" begin
    using PiccoloQuantumObjects: Isomorphisms.G, Isomorphisms.H, Isomorphisms.iso, Isomorphisms.ad_vec

    H_real = [1.0 2.0; 3.0 4.0]
    H_imag = [0.0 1.0; 1.0 0.0]
    H_complex = H_real + 1.0im * H_imag
    G_H = G(H_complex)

    @test H(G_H) ≈ H_complex

    @test G_H ≈ [0 1 1 2; 1 0 3 4; -1 -2 0 1; -3 -4 1 0]

    @test iso(H_complex) ≈ [1 2 0 -1; 3 4 -1 0; 0 1 1 2; 1 0 3 4]

    @test iso(-im * H_complex) ≈ G_H

    H = [0 1; 1 0]
    ad_H = ad_vec(H)
    @test ad_H ≈ [0 1 -1 0; 1 0 0 -1; -1 0 0 1; 0 -1 1 0]

    H = [0 -im; im 0]
    ad_H = ad_vec(H)
    @test ad_H ≈ [0 -im -im 0; im 0 0 -im; im 0 0 -im; 0 im im 0]
end

end
