
# ```@meta
# CollapsedDocStrings = true
# ```
using PiccoloQuantumObjects
using SparseArrays # for visualization

#=
# Isomorphisms for quantum objects

Linear algebra operations on quantum objects are often performed on real vectors and
matrices. We provide isomorphisms to convert between complex and real representations of
quantum objects. These isomorphisms are used internally by the [`QuantumSystem`](@ref) type
to perform quantum dynamics.

=#

#=
## Quantum states

- [`ket_to_iso`](@ref) is the real isomorphism of a quantum state `ψ ∈ ℂⁿ`
- [`iso_to_ket`](@ref) is the inverse isomorphism of a real vector `ψ̃ ∈ ℝ²ⁿ`

=#

ψ = [1; 2] + im * [3; 4]
ψ̃ = ket_to_iso(ψ)

#
iso_to_ket(ψ̃)

#=
## Quantum operators

We often need to convert a complex matrix `U` to a real vector `Ũ⃗`. We provoide the
following isomorphisms to convert between the two representations.
- [`iso_vec_to_operator(Ũ⃗::AbstractVector{ℝ})`](@ref)
- [`operator_to_iso_vec(U::AbstractVector{ℂ})`](@ref)
- [`iso_vec_to_iso_operator(Ũ⃗::AbstractVector{ℝ})`](@ref)
- [`iso_operator_to_iso_vec(Ũ::AbstractMatrix{ℝ})`](@ref)
- [`iso_operator_to_operator(Ũ::AbstractMatrix{ℝ})`](@ref)
- [`operator_to_iso_operator(U::AbstractMatrix{ℂ})`](@ref)

In additon, we provide [`mat(x::AbstractVector)`](@ref) to convert a vector `x` into a 
square matrix, as the inverse to Base's `vec`.

Notice that Julia uses column-major order.
=#

U = [1 5; 2 6] + im * [3 7; 4 8]
Ũ⃗ = operator_to_iso_vec(U)

# 
iso_vec_to_operator(Ũ⃗)


#=
## Density matrices

The isomorphisms for density matrices are:
- [`density_to_iso_vec(ρ::AbstractMatrix{ℂ})`](@ref)
- [`iso_vec_to_density(ρ̃::AbstractVector{ℝ})`](@ref)

!!! warning
    The isomorphism `density_to_iso_vec` is not the same as `operator_to_iso_vec`.

=#

ρ = [1 2; 3 4] + im * [5 6; 7 8]
ρ̃⃗ = density_to_iso_vec(ρ)


#=
# Quantum dynamics

The quantum dynamics isomorphisms, which correspond to these state transformations, are 
handled internally by the [`QuantumSystem`](@ref) type.

The [`Isomorphisms.iso`](@ref) isomorphism of a Hamiltonian ``H`` is:
```math
\text{iso}(H) := \widetilde{H} = \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```
where ``\Im(H)`` and ``\Re(H)`` are the imaginary and real parts of ``H`` and the tilde 
indicates the standard isomorphism of a complex valued matrix:
```math
\widetilde{H} := \mqty(1 & 0 \\ 0 & 1) \otimes \Re(H) + \mqty(0 & -1 \\ 1 & 0) \otimes \Im(H)
```

Hence, the generator [`Isomorphisms.G`](@ref) associated to a Hamiltonian ``H`` is:
```math
G(H) := \text{iso}(- i \widetilde{H}) = \mqty(1 & 0 \\ 0 & 1) \otimes \Im(H) - \mqty(0 & -1 \\ 1 & 0) \otimes \Re(H)
```

=#
