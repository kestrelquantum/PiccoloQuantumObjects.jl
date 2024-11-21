module EmbeddedOperators

export AbstractPiccoloOperator
export EmbeddedOperator

export embed
export unembed
export get_subspace_indices
export get_subspace_enr_indices
export get_leakage_indices
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
    embed(operator::AbstractMatrix{<:Number}, subspace::AbstractVector{Int}, levels::Int)

Embed an `operator` in the `subspace` of a larger matrix of size `levels x levels`.
"""
function embed(
    operator::AbstractMatrix{R}, subspace::AbstractVector{Int}, levels::Int
) where R <: Number
    @assert size(operator, 1) == size(operator, 2) "Operator must be square."
    op_embedded = zeros(R, levels, levels)
    op_embedded[subspace, subspace] = operator
    return op_embedded
end

"""
    unembed(matrix::AbstractMatrix{<:Number}, subspace::AbstractVector{Int})

Unembed a subspace operator from the `matrix`.

This is equivalent to calling `matrix[subspace, subspace]`.
"""
function unembed(matrix::AbstractMatrix{<:Number}, subspace::AbstractVector{Int})
    return matrix[subspace, subspace]
end

# ----------------------------------------------------------------------------- #
#                             Embedded Operator                                 #
# ----------------------------------------------------------------------------- #

@doc raw"""
    EmbeddedOperator

Embedded operator type to represent an operator embedded in a subspace of a larger 
quantum system.

The larger system, $\mathcal{X}$, can be decomposed into subspace and leakage components:

```math
    \mathcal{X} = \mathcal{X}_{\text{subspace}} \oplus \mathcal{X}_{\text{leakage}},
```

# Fields
- `operator::Matrix{ComplexF64}`: Embedded operator of size `prod(subsystem_levels) x prod(subsystem_levels)`.
- `subspace::Vector{Int}`: Indices of the subspace the operator is embedded in.
- `subsystem_levels::Vector{Int}`: Levels of the subsystems in the composite system.
"""
struct EmbeddedOperator
    operator::Matrix{ComplexF64}
    subspace::Vector{Int}
    subsystem_levels::Vector{Int}

    @doc raw"""
        EmbeddedOperator(op::Matrix{<:Number}, subspace::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})

    Create an embedded operator. The `operator` is embedded at the `subspace` of the
    system spanned by the `subsystem_levels`.

    # Example

    ```jldoctest
    julia> operator = kron([0 1; 1 0], [0 1; 1 0])
    4×4 Matrix{Int64}:
        0  0  0  1
        0  0  1  0
        0  1  0  0
        1  0  0  0
    julia> subspace = get_subspace_indices([1:2, 1:2], [3, 3])
    4-element Vector{Int64}:
        1
        2
        4
        5
    julia> subsystem_levels = [3, 3]
    julia> EmbeddedOperator(operator, subspace, subsystem_levels)
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
        subspace_operator::Matrix{<:Number},
        subspace::AbstractVector{Int},
        subsystem_levels::AbstractVector{Int}
    )
        embedded_operator = embed(
            Matrix{ComplexF64}(subspace_operator), subspace, prod(subsystem_levels)
        )
        return new(embedded_operator, subspace, subsystem_levels)
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

EmbeddedOperator(
    subspace_operator::Matrix{<:Number}, subspace::AbstractVector{Int}, levels::Int
) = EmbeddedOperator(subspace_operator, subspace, [levels])

function EmbeddedOperator(
    subspace_operator::AbstractMatrix{<:Number},
    system::QuantumSystem;
    subspace=1:size(subspace_operator, 1),
    levels=system.levels
)
    return EmbeddedOperator(
        subspace_operator, get_subspace_indices(subspace, levels), [system.levels]
    )
end

function EmbeddedOperator(op::Symbol, args...; kwargs...)
    if op ∉ keys(gates)
        throw(ArgumentError("Operator must be a valid gate. "
            *"See PiccoloQuantumObjects.gates.GATES dict for available gates."))
    end
    return EmbeddedOperator(GATES[op], args...; kwargs...)
end

# TODO: TEST & DOCUMENT & FIX
function EmbeddedOperator(
    subspace_operator::AbstractMatrix{<:Number},
    composite_system::CompositeQuantumSystem,
    subsystem_indices::AbstractVector{Int};
    subspaces=fill(1:2, length(composite_system.subsystems)),
)
    @assert all(diff(subsystem_indices) .== 1) "subsystem_indices must be consecutive (for now)."

    if size(subspace_operator, 1) == prod(length.(subspaces[subsystem_indices]))
        Is = Matrix{ComplexF64}.(I.(length.(subspaces)))
        Is[subsystem_indices[1]] = subspace_operator
        deleteat!(Is, subsystem_indices[2:end])
        subspace_operator = kron(Is...)
    elseif size(subspace_operator, 1) != prod(length.(subspaces))
        throw(ArgumentError(
            "Operator size ($(size(subspace_operator, 1))) must match "
            *"the product of subsystem subspaces ($(prod(length.(subspaces))))."
        ))
    end

    return EmbeddedOperator(
        op,
        get_subspace_indices(subspaces, composite_system.subsystem_levels),
        composite_system.subsystem_levels
    )
end

EmbeddedOperator(
    subspace_operator::AbstractMatrix{<:Number},
    composite_system::CompositeQuantumSystem,
    subsystem_index::Int;
    kwargs...
) = EmbeddedOperator(subspace_operator, composite_system, [subsystem_index]; kwargs...)

function EmbeddedOperator(
    subspace_operators::AbstractVector,
    composite_system::CompositeQuantumSystem,
    subsystem_indices::AbstractVector{Int};
    kwargs...
)
    embedded_operators = [
        EmbeddedOperator(op, composite_system, i; kwargs...)
        for (op, i) ∈ zip(subspace_operators, subsystem_indices)
    ]
    return *(embedded_operators...)
end

# ----------------------------------------------------------------------------- #
#                           EmbeddedOperator operations                         #
# ----------------------------------------------------------------------------- #

"""
    embed(subspace_op::AbstractMatrix{<:Number}, embedded_op::EmbeddedOperator)

Embed the `subspace_op` in the subspace of a larger system defined by `op`.
"""
embed(
    subspace_op::AbstractMatrix{<:Number}, embedded_op::EmbeddedOperator
) = embed(subspace_op, embedded_op.subspace, prod(embedded_op.subsystem_levels))

"""
    unembed(embedded_op::EmbeddedOperator)::Matrix{ComplexF64}

Unembed an embedded operator, returning the original operator.
"""
unembed(
    op::EmbeddedOperator
)::Matrix{ComplexF64} = op.operator[op.subspace, op.subspace]

"""
    unembed(op::AbstractMatrix, embedded_op::EmbeddedOperator)

Unembed a sub-matrix from the `op` at the subspace defined by `embedded_op`.
"""
unembed(
    op::AbstractMatrix{<:Number}, embedded_op::EmbeddedOperator
) = op[embedded_op.subspace, embedded_op.subspace]

Base.size(op::EmbeddedOperator, args...) = size(op.operator, args...)

function Base.:*(op1::EmbeddedOperator, op2::EmbeddedOperator)
    @assert size(op1) == size(op2) "Operators must be of the same size."
    @assert op1.subspace == op2.subspace "Operators must have the same subspace."
    @assert op1.subsystem_levels == op2.subsystem_levels "Operators must have the same subsystem levels."
    return EmbeddedOperator(
        unembed(op1) * unembed(op2),
        op1.subspace,
        op1.subsystem_levels
    )
end

function Base.kron(op1::EmbeddedOperator, op2::EmbeddedOperator)
    levels = [size(op1, 1), size(op2, 2)]
    indices = get_subspace_indices(
        [op1.subspace, op2.subspace], levels
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
    get_subspace_indices(subspace::AbstractVector{Int}, levels::Int)    
    get_subspace_indices(subspaces::Vector{<:AbstractVector{Int}}, subsystem_levels::AbstractVector{Int})
    get_subspace_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2)
    get_subspace_indices(op::EmbeddedOperator)

Get the indices for the provided subspace of the quantum system.

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

function get_subspace_indices(subspace::AbstractVector{Int}, levels::Int)
    @assert all(1 ≤ s ≤ levels for s ∈ subspace)
    return subspace
end

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

get_subspace_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2) =
    get_subspace_indices(fill(subspace, length(subsystem_levels)), subsystem_levels)

get_subspace_indices(op::EmbeddedOperator) = op.subspace

"""
    get_subspace_enr_indices(excitation_restriction::Int, subsystem_levels::AbstractVector{Int})

Get the indices for the subspace of the quantum system with an excitation restriction.

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
    get_leakage_indices(subspace::AbstractVector{Int}, levels::Int)
    get_leakage_indices(subspaces::AbstractVector{<:AbstractVector{Int}}, subsystem_levels::AbstractVector{Int})
    get_leakage_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2)
    get_leakage_indices(op::EmbeddedOperator)

Get the indices for the states that are outside of the provided subspace of the quantum
system.

# Example

```jldoctest

julia> subspaces = [1:2, 1:2]
julia> subsystem_levels = [3, 3]
julia> get_leakage_indices(subspaces, subsystem_levels)
5-element Vector{Int64}:
    3
    6
    7
    8
    9

julia> subspace = 1:2
julia> levels = 3
julia> get_leakage_indices(subspace, levels)
1-element Vector{Int64}:
    3

```

"""
function get_leakage_indices end

get_leakage_indices(subspace::AbstractVector{Int}, levels::Int) =
   setdiff(1:levels, subspace)

function get_leakage_indices(
    subspaces::AbstractVector{<:AbstractVector{Int}},
    subsystem_levels::AbstractVector{<:Int}
)
    return get_leakage_indices(
        get_subspace_indices(subspaces, subsystem_levels),
        prod(subsystem_levels)
    )
end

get_leakage_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2) =
    get_leakage_indices(fill(subspace, length(subsystem_levels)), subsystem_levels)

get_leakage_indices(op::EmbeddedOperator) =
    get_leakage_indices(op.subspace, size(op, 1))

"""
    get_iso_vec_subspace_indices(subspace::AbstractVector{Int}, subsystem_levels::AbstractVector{Int})
    get_iso_vec_subspace_indices(op::EmbeddedOperator)

Get the indices for the subspace in the isomorphic vector space for operators.

```jldoctest

julia> subspace = 1:2
julia> levels = 3
julia> get_iso_vec_subspace_indices(subspace, levels)
8-element Vector{Int64}:
    1
    2
    4
    5
    7
    8
    10
    11

```

"""
function get_iso_vec_subspace_indices end

function get_iso_vec_subspace_indices(
    subspace::AbstractVector{Int},
    levels::Int
)
    iso_subspace_indices = Int[]
    for sⱼ ∈ subspace
        for sᵢ ∈ subspace
            push!(iso_subspace_indices, 2levels * (sⱼ - 1) + sᵢ)
        end
        for sᵢ ∈ subspace
            push!(iso_subspace_indices, 2levels * (sⱼ - 1) + sᵢ + levels)
        end
    end
    return iso_subspace_indices
end

function get_iso_vec_subspace_indices(
    subspaces::AbstractVector{<:AbstractVector{Int}},
    subsystem_levels::AbstractVector{Int}
)
    subspace = get_subspace_indices(subspaces, subsystem_levels)
    levels = prod(subsystem_levels)
    return get_iso_vec_subspace_indices(subspace, levels)
end

get_iso_vec_subspace_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2) =
    get_iso_vec_subspace_indices(fill(subspace, length(subsystem_levels)), subsystem_levels)

get_iso_vec_subspace_indices(op::EmbeddedOperator) =
    get_iso_vec_subspace_indices(op.subspace, size(op, 1))


"""
    get_iso_vec_leakage_indices(subspace::AbstractVector{Int}, levels::Int)
    get_iso_vec_leakage_indices(subspaces::AbstractVector{<:AbstractVector{Int}}, subsystem_levels::AbstractVector{Int})
    get_iso_vec_leakage_indices(subsystem_levels::AbstractVector{Int}; subspace=1:2)
    get_iso_vec_leakage_indices(op::EmbeddedOperator)

Get the indices for the leakage in the isomorphic vector space for operators.

# Example

```jldoctest
julia> get_iso_vec_leakage_indices(1:2, 3)
4-element Vector{Int64}:
    3
    6
    9
    12

```

"""
function get_iso_vec_leakage_indices end

function get_iso_vec_leakage_indices(
    subspace::AbstractVector{Int},
    levels::Int;
    ignore_pure_leakage::Bool=true
)
    leakage = get_leakage_indices(subspace, levels)
    iso_leakage_indices = Int[]
    
    rows = ignore_pure_leakage ? subspace : 1:levels

    for j ∈ 1:levels
        # Real
        if j ∈ leakage
            # (subspace, leakage) x leakage
            for i ∈ rows
                push!(iso_leakage_indices, 2levels * (j - 1) + i)
            end
        elseif j ∈ subspace
            # leakage x subspace
            for lᵢ ∈ leakage
                push!(iso_leakage_indices, 2levels * (j - 1) + lᵢ)
            end
        end

        # Imaginary
        if j ∈ leakage
            for i ∈ rows
                push!(iso_leakage_indices, 2levels * (j - 1) + i + levels)
            end
        elseif j ∈ subspace
            for lᵢ ∈ leakage
                push!(iso_leakage_indices, 2levels * (j - 1) + lᵢ + levels)
            end
        end
    end

    return iso_leakage_indices
end

function get_iso_vec_leakage_indices(
    subspaces::AbstractVector{<:AbstractVector{Int}},
    subsystem_levels::AbstractVector{Int};
    kwargs...
)
    subspace = get_subspace_indices(subspaces, subsystem_levels)
    levels = prod(subsystem_levels)
    return get_iso_vec_leakage_indices(subspace, levels; kwargs...)
end

function get_iso_vec_leakage_indices(
    subsystem_levels::AbstractVector{Int}; subspace=1:2, kwargs...
)
    return get_iso_vec_leakage_indices(
        fill(subspace, length(subsystem_levels)), subsystem_levels; kwargs...
    )
end

get_iso_vec_leakage_indices(op::EmbeddedOperator; kwargs...) =
    get_iso_vec_leakage_indices(op.subspace, size(op, 1); kwargs...)

# ****************************************************************************** #

@testitem "Embedded operator" begin
    using LinearAlgebra: I

    # Embed X
    op = Matrix{ComplexF64}([0 1; 1 0])
    embedded_op = Matrix{ComplexF64}([0 1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0])
    @test embed(op, 1:2, 4) == embedded_op
    embedded_op_struct = EmbeddedOperator(op, 1:2, 4)
    @test embedded_op_struct.operator == embedded_op
    @test embedded_op_struct.subspace == 1:2
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
    subspace = get_subspace_indices([1:2, 1:2], [3, 3])
    @test embed(op2, subspace, 9) == embedded_op2
    embedded_op2_struct = EmbeddedOperator(op2, subspace, [3, 3])
    @test embedded_op2_struct.operator == embedded_op2
    @test embedded_op2_struct.subspace == subspace
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
        subspace=get_subspace_indices([1:2, 1:2], [3, 3])
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
    leakage = get_leakage_indices([1:2, 1:2], [3, 3])
    @test leakage == [3, 6, 7, 8, 9]

    leakage = get_leakage_indices([1, 2], 3)
    @test leakage == [3]
end

@testitem "Single system iso vec indices" begin
    U_real = [1 4 7; 2 5 8; 3 6 9]
    U = U_real - im * U_real
    Ũ⃗ = operator_to_iso_vec(U)

    subspace = get_subspace_indices(1:2, 3)
    iso_vec_subspace = get_iso_vec_subspace_indices(1:2, 3)
    all_leakage = get_iso_vec_leakage_indices(1:2, 3, ignore_pure_leakage=false)
    ignore_pure_leakage = get_iso_vec_leakage_indices(1:2, 3, ignore_pure_leakage=true)
    ignore_pure_leakage_expect = [3, 6, 9, 12, 13, 14, 16, 17]
    all_leakage_expect = [3, 6, 9, 12, 13, 14, 15, 16, 17, 18]

    @test Ũ⃗[iso_vec_subspace] == operator_to_iso_vec(unembed(U, subspace))
    @test setdiff(setdiff(1:length(Ũ⃗), iso_vec_subspace), all_leakage) == []
    @test ignore_pure_leakage == ignore_pure_leakage_expect
    @test all_leakage == all_leakage_expect

    @test issorted(ignore_pure_leakage)
    @test issorted(all_leakage)
end

@testitem "Composite system iso vec indices" begin
    U_real = kron(rand(3, 3), rand(3, 3))
    U = U_real - im * U_real
    Ũ⃗ = operator_to_iso_vec(U)
    
    subspace = get_subspace_indices([3, 3])
    iso_vec_subspace = get_iso_vec_subspace_indices([3, 3])
    all_leakage = get_iso_vec_leakage_indices([3, 3], ignore_pure_leakage=false)
    ignore_pure_leakage = get_iso_vec_leakage_indices([3, 3], ignore_pure_leakage=true)

    @test Ũ⃗[iso_vec_subspace] == operator_to_iso_vec(unembed(U, subspace))
    @test setdiff(setdiff(1:length(Ũ⃗), iso_vec_subspace), all_leakage) == []

    @test issorted(all_leakage)
    @test issorted(ignore_pure_leakage)
end

@testitem "Embedded operator subspace" begin
    op = EmbeddedOperator([1 2; 3 4] + im * [5 6; 7 8], 1:2, 3)
    
    @test get_subspace_indices(op) == [1, 2]
    @test get_leakage_indices(op) == [3]

    @test get_iso_vec_subspace_indices(op) == [1, 2, 4, 5, 7, 8, 10, 11]
    
    # drop the pure leakage state at op.operator[3, 3]
    @test get_iso_vec_leakage_indices(op) == [3, 6, 9, 12, 13, 14, 16, 17]
end

end
