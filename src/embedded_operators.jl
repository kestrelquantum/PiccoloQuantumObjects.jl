module EmbeddedOperators

export AbstractPiccoloOperator
export EmbeddedOperator

export embed
export unembed
export get_subspace_indices
export get_subspace_enr_indices
export get_subspace_leakage_indices
export get_iso_vec_leakage_indices
export get_iso_vec_subspace_indices

using ..Gates
using ..Isomorphisms
using ..QuantumObjectUtils
using ..QuantumSystems
using ..CompositeQuantumSystems

using LinearAlgebra
using TestItemRunner

# ----------------------------------------------------------------------------- #
#                             Embedding operations                              #
# ----------------------------------------------------------------------------- #

"""
    embed(operator::AbstractMatrix{<:Number}, subspace_indices::AbstractVector{Int}, levels::Int)

Embed an `operator` at the `subspace_indices` in a larger matrix of size `levels x levels`.
"""
function embed(
    operator::AbstractMatrix{R}, subspace_indices::AbstractVector{Int}, levels::Int
) where R <: Number
    @assert size(operator, 1) == size(operator, 2) "Operator must be square."
    op_embedded = zeros(R, levels, levels)
    op_embedded[subspace_indices, subspace_indices] = operator
    return op_embedded
end

"""
    unembed(matrix::AbstractMatrix{<:Number}, subspace_indices::AbstractVector{Int})

Unembed a subspace operator from the `matrix`.

This is equivalent to calling `matrix[subspace_indices, subspace_indices]`.
"""
function unembed(matrix::AbstractMatrix{<:Number}, subspace_indices::AbstractVector{Int})
    return matrix[subspace_indices, subspace_indices]
end

# ----------------------------------------------------------------------------- #
#                             Embedded Operator                                 #
# ----------------------------------------------------------------------------- #

@doc raw"""
    EmbeddedOperator

Embedded operator type to represent an operator embedded in a subspace of a larger 
quantum system.

The larger system $\mathcal{X}$ can be decomposed into its subspace and leakage components:

```math
    \mathcal{X} = \mathcal{X}_{\text{subspace}} \oplus \mathcal{X}_{\text{leakage}},
```

# Fields
- `operator::Matrix{ComplexF64}`: Embedded operator of size `prod(subsystem_levels) x prod(subsystem_levels)`.
- `subspace_indices::Vector{Int}`: Indices of the subspace the operator is embedded in.
- `subsystem_levels::Vector{Int}`: Levels of the subsystems in the composite system.
"""
struct EmbeddedOperator
    operator::Matrix{ComplexF64}
    subspace_indices::Vector{Int}
    subsystem_levels::Vector{Int}

    @doc raw"""
        EmbeddedOperator(op::Matrix{<:Number}, subspace_indices::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})

    Create an embedded operator. The operator `op` is embedded in the subspace defined by 
    `subspace_indices` in `subsystem_levels`.

    # Example

    ```jldoctest
    julia> operator = kron([0 1; 1 0], [0 1; 1 0])
    4×4 Matrix{Int64}:
        0  0  0  1
        0  0  1  0
        0  1  0  0
        1  0  0  0
    julia> subspace_indices = get_subspace_indices([1:2, 1:2], [3, 3])
    4-element Vector{Int64}:
        1
        2
        4
        5
    julia> subsystem_levels = [3, 3]
    julia> EmbeddedOperator(operator, subspace_indices, subsystem_levels)
    9×9 Matrix{ComplexF64}:
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
    ```
    """
    function EmbeddedOperator(
        op::Matrix{<:Number},
        subspace_indices::AbstractVector{Int},
        subsystem_levels::AbstractVector{Int}
    )

        op_embedded = embed(Matrix{ComplexF64}(op), subspace_indices, prod(subsystem_levels))
        return new(op_embedded, subspace_indices, subsystem_levels)
    end
end

# ----------------------------------------------------------------------------- #
#                             AbstractPiccoloOperator                           #
# ----------------------------------------------------------------------------- #

"""
    AbstractPiccoloOperator

Union type for operators.
"""
const AbstractPiccoloOperator = Union{AbstractMatrix{<:Number}, EmbeddedOperator}

# ----------------------------------------------------------------------------- #
#                  Additional EmbeddedOperator constructors                     #
# ----------------------------------------------------------------------------- #

EmbeddedOperator(op::Matrix{<:Number}, subspace_indices::AbstractVector{Int}, levels::Int) =
    EmbeddedOperator(op, subspace_indices, [levels])

function EmbeddedOperator(
    op::AbstractMatrix{<:Number},
    system::QuantumSystem;
    subspace_indices=1:size(op, 1),
    levels=system.levels
)
    return EmbeddedOperator(
        op,
        get_subspace_indices(subspace_indices, levels),
        [system.levels]
    )
end

function EmbeddedOperator(
    op::AbstractMatrix{<:Number},
    csystem::CompositeQuantumSystem,
    op_subsystem_indices::AbstractVector{Int};
    subspaces=fill(1:2, length(csystem.subsystems)),
)
    @assert all(diff(op_subsystem_indices) .== 1) "op_subsystem_indices must be consecutive (for now)."

    if size(op, 1) == prod(length.(subspaces[op_subsystem_indices]))
        Is = Matrix{ComplexF64}.(I.(length.(subspaces)))
        Is[op_subsystem_indices[1]] = op
        deleteat!(Is, op_subsystem_indices[2:end])
        op = kron(Is...)
    else
        @assert(
            size(op, 1) == prod(length.(subspaces)),
            """\n
                Operator size ($(size(op, 1))) must match product of subsystem subspaces ($(prod(length.(subspaces)))).
            """
        )
    end

    subspace_indices = get_subspace_indices(subspaces, csystem.subsystem_levels)

    return EmbeddedOperator(
        op,
        subspace_indices,
        csystem.subsystem_levels
    )
end

function EmbeddedOperator(
    op::AbstractMatrix{<:Number},
    csystem::CompositeQuantumSystem,
    op_subsystem_index::Int;
    kwargs...
)
    return EmbeddedOperator(
        op,
        csystem,
        [op_subsystem_index];
        kwargs...
    )
end

function EmbeddedOperator(op::Symbol, args...; kwargs...)
    if op ∉ keys(gates)
        throw(ArgumentError("Operator must be a valid gate. "
            *"See PiccoloQuantumObjects.gates.GATES dict for available gates."))
    end
    return EmbeddedOperator(GATES[op], args...; kwargs...)
end

function EmbeddedOperator(
    ops::AbstractVector{Symbol},
    sys::CompositeQuantumSystem,
    op_indices::AbstractVector{Int}
)
    ops_embedded = [
        EmbeddedOperator(op, sys, op_indices[i])
            for (op, i) ∈ zip(ops, op_indices)
    ]
    return *(ops_embedded...)
end

# ----------------------------------------------------------------------------- #
#                           EmbeddedOperator operations                         #
# ----------------------------------------------------------------------------- #

"""
    embed(matrix::AbstractMatrix{<:Number}, op::EmbeddedOperator)

Embed an operator `matrix` in the subspace of a larger system defined by `op`.
"""
function embed(matrix::AbstractMatrix{<:Number}, op::EmbeddedOperator)
    return embed(matrix, op.subspace_indices, prod(op.subsystem_levels))
end

"""
    unembed(op::EmbeddedOperator)::Matrix{ComplexF64}

Unembed an embedded operator, returning the original operator.
"""
function unembed(op::EmbeddedOperator)::Matrix{ComplexF64}
    return op.operator[op.subspace_indices, op.subspace_indices]
end

"""
    unembed(matrix::AbstractMatrix, op::EmbeddedOperator)

Unembed an operator `matrix` from the subspace defined by `op`.
"""
function unembed(matrix::AbstractMatrix{<:Number}, op::EmbeddedOperator)
    return matrix[op.subspace_indices, op.subspace_indices]
end

Base.size(op::EmbeddedOperator) = size(op.operator)
Base.size(op::EmbeddedOperator, dim::Union{Int, Nothing}) = size(op.operator, dim)

function Base.:*(op1::EmbeddedOperator, op2::EmbeddedOperator)
    @assert size(op1) == size(op2) "Operators must be of the same size."
    @assert op1.subspace_indices == op2.subspace_indices "Operators must have the same subspace."
    @assert op1.subsystem_levels == op2.subsystem_levels "Operators must have the same subsystem levels."
    return EmbeddedOperator(
        unembed(op1) * unembed(op2),
        op1.subspace_indices,
        op1.subsystem_levels
    )
end

function Base.kron(op1::EmbeddedOperator, op2::EmbeddedOperator)
    levels = [size(op1, 1), size(op2, 2)]
    indices = get_subspace_indices(
        [op1.subspace_indices, op2.subspace_indices], levels
    )
    return EmbeddedOperator(kron(unembed(op1), unembed(op2)), indices, levels)
end

# ----------------------------------------------------------------------------- #
#                            Subspace Indices                                   #
# ----------------------------------------------------------------------------- #

basis_labels(subsystem_levels::AbstractVector{Int}; baseline::Int=1) =
    kron([""], [string.(baseline:levels - 1 + baseline) for levels ∈ subsystem_levels]...)

basis_labels(subsystem_level::Int; kwargs...) = basis_labels([subsystem_level]; kwargs...)

"""
    get_subspace_indices(subspaces::Vector{<:AbstractVector{Int}}, subsystem_levels::AbstractVector{Int})
    get_subspace_indices(subspace::AbstractVector{Int}, levels::Int)
    get_subspace_indices(levels::AbstractVector{Int}; subspace=1:2)

Get the indices for the provided subspace(s) of a composite quantum system.

# Example

A two-qubit subspace of two 3-level systems:
```jldoctest
julia> subspaces = [1:2, 1:2]
julia> subsystem_levels = [3, 3]
julia> get_subspace_indices(subspaces, subsystem_levels)
4-element Vector{Int64}:
    1
    2
    4
    5
```

A two qubit subspace of a single 3-level system:
```jldoctest
julia> get_subspace_indices([1, 2], 3)
[1, 2]

julia> get_subspace_indices([1:2, 1:2], [3, 3])
[1, 2, 4, 5]
````

Two 3-level systems with a default (qubit) subspace:
```jldoctest
julia> get_subspace_indices([3, 3])
4-element Vector{Int64}:
    1
    2
    4
    5
```

"""
function get_subspace_indices end

function get_subspace_indices(
    subspaces::AbstractVector{<:AbstractVector{Int}},
    subsystem_levels::AbstractVector{Int}
)
    @assert length(subspaces) == length(subsystem_levels)
    return findall(
        b -> all(l ∈ subspaces[i] for (i, l) ∈ enumerate([parse(Int, bᵢ) for bᵢ ∈ b])),
        basis_labels(subsystem_levels, baseline=1)
    )
end

get_subspace_indices(subspace::AbstractVector{Int}, levels::Int) =
    get_subspace_indices([subspace], [levels])

get_subspace_indices(levels::AbstractVector{Int}; subspace=1:2) =
    get_subspace_indices(fill(subspace, length(levels)), levels)

"""
    get_subspace_enr_indices(excitation_restriction::Int, subsystem_levels::AbstractVector{Int})

Get the indices for the subspace of composite quantum system with an excitation restriction.

# Example

Choose only the ground state and single excitation states of two 3-level systems:
```jldoctest

julia> get_subspace_enr_indices(1, [3, 3])
3-element Vector{Int64}:
    1
    2
    4

```

"""
function get_subspace_enr_indices(excitation_restriction::Int, subsystem_levels::AbstractVector{Int})
    # excitation_number uses baseline of zero
    return findall(
        b -> sum([parse(Int, bᵢ) for bᵢ ∈ b]) ≤ excitation_restriction,
        basis_labels(subsystem_levels, baseline=0)
    )
end

"""
    get_subspace_leakage_indices(subspaces::AbstractVector{<:AbstractVector{Int}}, subsystem_levels::AbstractVector{Int})
    get_subspace_leakage_indices(subspace::AbstractVector{Int}, levels::Int)
    get_subspace_leakage_indices(op::EmbeddedOperator)

Get the indices for the states that are outside of the provided subspace(s).

# Example

```jldoctest

julia> subspaces = [1:2, 1:2]
julia> subsystem_levels = [3, 3]
julia> get_subspace_leakage_indices(subspaces, subsystem_levels)
5-element Vector{Int64}:
    3
    6
    7
    8
    9

julia> subspace = 1:2
julia> levels = 3
julia> get_subspace_leakage_indices(subspace, levels)
1-element Vector{Int64}:
    3

```

"""
function get_subspace_leakage_indices end

function get_subspace_leakage_indices(
    subspaces::AbstractVector{<:AbstractVector{Int}},
    subsystem_levels::AbstractVector{<:Int}
)
    return get_subspace_leakage_indices(
        get_subspace_indices(subspaces, subsystem_levels),
        prod(subsystem_levels)
    )
end

get_subspace_leakage_indices(subspace_indices::AbstractVector{Int}, levels::Int) =
    setdiff(1:levels, subspace_indices)

get_subspace_leakage_indices(op::EmbeddedOperator) =
    get_subspace_leakage_indices(op.subspace_indices, size(op)[1])

"""
    get_iso_vec_subspace_indices(subspace_indices::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})
    get_iso_vec_subspace_indices(op::EmbeddedOperator)

Get the indices for the subspace in the isomorphic vector space for operators.

# TODO: Example

"""
function get_iso_vec_subspace_indices end

function get_iso_vec_subspace_indices(
    subspace_indices::AbstractVector{Int},
    subsystem_levels::AbstractVector{Int}
)
    N = prod(subsystem_levels)
    iso_subspace_indices = Int[]
    for sⱼ ∈ subspace_indices
        for sᵢ ∈ subspace_indices
            push!(iso_subspace_indices, index(sⱼ, sᵢ, 2N))
        end
        for sᵢ ∈ subspace_indices
            push!(iso_subspace_indices, index(sⱼ, sᵢ + N, 2N))
        end
    end
    return iso_subspace_indices
end

get_iso_vec_subspace_indices(op::EmbeddedOperator) =
    get_iso_vec_subspace_indices(op.subspace_indices, op.subsystem_levels)


"""
    get_iso_vec_leakage_indices(subspace_indices::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})
    get_iso_vec_leakage_indices(op::EmbeddedOperator)

Get the indices for the leakage in the isomorphic vector space for operators.

# TODO: Example

"""
function get_iso_vec_leakage_indices end

function get_iso_vec_leakage_indices(
    subspace_indices::AbstractVector{Int},
    subsystem_levels::AbstractVector{Int}
)
    N = prod(subsystem_levels)
    leakage_indices = get_subspace_leakage_indices(subspace_indices, N)
    iso_leakage_indices = Int[]
    for sⱼ ∈ subspace_indices
        for lᵢ ∈ leakage_indices
            push!(iso_leakage_indices, index(sⱼ, lᵢ, 2N))
        end
        for lᵢ ∈ leakage_indices
            push!(iso_leakage_indices, index(sⱼ, lᵢ + N, 2N))
        end
    end
    return iso_leakage_indices
end

get_iso_vec_leakage_indices(op::EmbeddedOperator) =
    get_iso_vec_leakage_indices(op.subspace_indices, op.subsystem_levels)

# =========================================================================== #

@testitem "Basis labels" begin
    levels = [3, 3]
    labels = ["11", "12", "13", "21", "22", "23", "31", "32", "33"]
    @test EmbeddedOperators.basis_labels(levels, baseline=1) == labels

    labels = ["1", "2", "3"]
    @test EmbeddedOperators.basis_labels(3, baseline=1) == labels
    @test EmbeddedOperators.basis_labels([3], baseline=1) == labels

    labels = ["0", "1", "2"]
    @test EmbeddedOperators.basis_labels(3, baseline=0) == labels
    @test EmbeddedOperators.basis_labels([3], baseline=0) == labels

    levels = [2, 2]
    labels = ["00", "01", "10", "11"]
    @test EmbeddedOperators.basis_labels(levels, baseline=0) == labels
end

@testitem "Subspace Indices" begin
    @test get_subspace_indices([1, 2], 3) == [1, 2]
    # 2 * 2 = 4 elements
    @test get_subspace_indices([1:2, 1:2], [3, 3]) == [1, 2, 4, 5]
    # 1 * 1 = 1 element
    @test get_subspace_indices([[2], [2]], [3, 3]) == [5]
    # 1 * 2 = 2 elements
    @test get_subspace_indices([[2], 1:2], [3, 3]) == [4, 5]
end

@testitem "Subspace ENR Indices" begin
    # 00, 01, 02x, 10, 11x, 12x, 20x, 21x, 22x
    @test get_subspace_enr_indices(1, [3, 3]) == [1, 2, 4]
    # 00, 01, 02, 10, 11, 12x, 20, 21x, 22x
    @test get_subspace_enr_indices(2, [3, 3]) == [1, 2, 3, 4, 5, 7]
    # 00, 01, 02, 10, 11, 12, 20, 21, 22x
    @test get_subspace_enr_indices(3, [3, 3]) == [1, 2, 3, 4, 5, 6, 7, 8]
    # 00, 01, 02, 10, 11, 12, 20, 21, 22
    @test get_subspace_enr_indices(4, [3, 3]) == [1, 2, 3, 4, 5, 6, 7, 8, 9]
end

@testitem "Subspace Leakage Indices" begin
    leakage = get_subspace_leakage_indices([1:2, 1:2], [3, 3])
    @test leakage == [3, 6, 7, 8, 9]

    leakage = get_subspace_leakage_indices([1, 2], 3)
    @test leakage == [3]
end

@testitem "Embedded operator" begin
    using LinearAlgebra: I

    # Embed X
    op = Matrix{ComplexF64}([0 1; 1 0])
    embedded_op = Matrix{ComplexF64}([0 1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0])
    @test embed(op, 1:2, 4) == embedded_op
    embedded_op_struct = EmbeddedOperator(op, 1:2, 4)
    @test embedded_op_struct.operator == embedded_op
    @test embedded_op_struct.subspace_indices == 1:2
    @test embedded_op_struct.subsystem_levels == [4]

    # Properties
    @test size(embedded_op_struct) == size(embedded_op)
    @test size(embedded_op_struct, 1) == size(embedded_op, 1)

    # X^2 = I
    x2 = (embedded_op_struct * embedded_op_struct).operator
    id = embed(I(2), embedded_op_struct)
    @test x2 == id

    # Embed X twice
    op2 = kron(op, op)
    embedded_op2 = [
        0  0  0  0  1  0  0  0  0;
        0  0  0  1  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0;
        0  1  0  0  0  0  0  0  0;
        1  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0
    ]
    subspace_indices = get_subspace_indices([1:2, 1:2], [3, 3])
    @test embed(op2, subspace_indices, 9) == embedded_op2
    embedded_op2_struct = EmbeddedOperator(op2, subspace_indices, [3, 3])
    @test embedded_op2_struct.operator == embedded_op2
    @test embedded_op2_struct.subspace_indices == subspace_indices
    @test embedded_op2_struct.subsystem_levels == [3, 3]
end

@testitem "Embedded operator from system" begin
    CZ = GATES[:CZ]
    a = [0 1 0; 0 0 1; 0 0 0]
    σ_x = a + a'
    σ_y = -1im*(a - a')
    system = QuantumSystem([kron(σ_x, σ_x), kron(σ_y, σ_y)])

    op_explicit_qubit = EmbeddedOperator(
        CZ,
        system,
        subspace_indices=get_subspace_indices([1:2, 1:2], [3, 3])
    )
    op_implicit_qubit = EmbeddedOperator(CZ, system)
    # This does not work (implicit puts indicies in 1:4)
    @test op_implicit_qubit.operator != op_explicit_qubit.operator
    # But the ops are the same
    @test unembed(op_explicit_qubit) == unembed(op_implicit_qubit)
    @test unembed(op_implicit_qubit) == CZ
end

@testitem "Embedded operator from composite system" begin
    @test_skip nothing
end

@testitem "Embedded operator kron" begin
    Z = GATES[:Z]
    Ẑ = EmbeddedOperator(Z, 1:2, [4])
    @test unembed(kron(Ẑ, Ẑ)) == kron(Z, Z)
end


end
