module DirectSums

export add_suffix
export get_suffix
export remove_suffix
export get_suffix_label
export direct_sum

using SparseArrays
using TestItems
using NamedTrajectories

using ..QuantumSystems


"""
    direct_sum(A::AbstractMatrix, B::AbstractMatrix)

Returns the direct sum of two matrices.
"""
function direct_sum(A::AbstractMatrix, B::AbstractMatrix)
    return [A spzeros((size(A, 1), size(B, 2))); spzeros((size(B, 1), size(A, 2))) B]
end

"""
    direct_sum(A::SparseMatrixCSC, B::SparseMatrixCSC)

Returns the direct sum of two sparse matrices.
"""
function direct_sum(A::SparseMatrixCSC, B::SparseMatrixCSC)
    return blockdiag(A, B)
end

"""
    direct_sum(Ã⃗::AbstractVector, B̃⃗::AbstractVector)

Returns the direct sum of two iso_vec operators.
"""
function direct_sum(Ã⃗::AbstractVector, B̃⃗::AbstractVector)
    return operator_to_iso_vec(
        direct_sum(
            iso_vec_to_operator(Ã⃗),
            iso_vec_to_operator(B̃⃗)
        )
    )
end

"""
    direct_sum(sys1::QuantumSystem, sys2::QuantumSystem)

Returns the direct sum of two `QuantumSystem` objects.
"""
function direct_sum(sys1::QuantumSystem, sys2::QuantumSystem)
    @assert sys1.n_drives == sys2.n_drives "System 1 drives ($(sys1.n_drives)) must equal System 2 drives ($(sys2.n_drives))"
    n_drives = sys1.n_drives
    H = a -> direct_sum(sys1.H(a), sys2.H(a))
    direct_sum_params = Dict{Symbol, Dict{Symbol, Any}}()
    if haskey(sys1.params, :system_1)
        n_systems = length(keys(sys1.params))
        direct_sum_params = sys1.params
        if haskey(sys2.params, :system_1)
            for i = 1:length(keys(sys2.params))
                direct_sum_params[Symbol("system_$(n_systems + i)")] =
                    sys2.params[Symbol("system_$(i)")]
            end
        else
            direct_sum_params[Symbol("system_$(n_systems + 1)")] = sys2.params
        end
    else
        direct_sum_params[:system_1] = sys1.params
        if haskey(sys2.params, :system_1)
            n_systems = length(keys(sys2.params))
            for i = 1:length(keys(sys2.params))
                direct_sum_params[Symbol("system_$(1 + i)")] =
                    sys2.params[Symbol("system_$(i)")]
            end
        else
            direct_sum_params[:system_2] = sys2.params
        end
    end
    return QuantumSystem(H, n_drives; params=direct_sum_params)
end

direct_sum(systems::AbstractVector{<:QuantumSystem}) = reduce(direct_sum, systems)

# Add suffix utilities
# -----------------------
Base.startswith(symb::Symbol, prefix::AbstractString) = startswith(String(symb), prefix)
Base.startswith(symb::Symbol, prefix::Symbol) = startswith(String(symb), String(prefix))

add_suffix(symb::Symbol, suffix::String) = Symbol(string(symb, suffix))
add_suffix(symbs::Tuple, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    Tuple(s ∈ exclude ? s : add_suffix(s, suffix) for s ∈ symbs)
add_suffix(symbs::AbstractVector, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    [s ∈ exclude ? s : add_suffix(s, suffix) for s ∈ symbs]
add_suffix(d::Dict{Symbol, Any}, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    typeof(d)(k ∈ exclude ? k : add_suffix(k, suffix) => v for (k, v) ∈ d)

function add_suffix(nt::NamedTuple, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[])
    symbs = Tuple(k ∈ exclude ? k : add_suffix(k, suffix) for k ∈ keys(nt))
    return NamedTuple{symbs}(values(nt))
end

function add_suffix(components::Union{Tuple, AbstractVector}, traj::NamedTrajectory, suffix::String)
    return add_suffix(get_components(components, traj), suffix)
end

function add_suffix(traj::NamedTrajectory, suffix::String)
    # TODO: Inplace
    # Timesteps are appended because of bounds and initial/final constraints.
    component_names = vcat(traj.state_names..., traj.control_names...)
    components = add_suffix(component_names, traj, suffix)
    controls = add_suffix(traj.control_names, suffix)
    return NamedTrajectory(
        components;
        controls=controls,
        timestep=traj.timestep isa Symbol ? add_suffix(traj.timestep, suffix) : traj.timestep,
        bounds=add_suffix(traj.bounds, suffix),
        initial=add_suffix(traj.initial, suffix),
        final=add_suffix(traj.final, suffix),
        goal=add_suffix(traj.goal, suffix)
    )
end

# function add_suffix(sys::QuantumSystem, suffix::String)
#     return QuantumSystem(
#         sys.H_drift,
#         sys.H_drives
#     )
# end

# get suffix label utilities
# --------------------

function get_suffix_label(s::String, pre::String)::String
    if startswith(s, pre)
        return chop(s, head=length(pre), tail=0)
    else
        error("Prefix '$pre' not found at the start of '$s'")
    end
end

get_suffix_label(symb::Symbol, pre::Symbol) = get_suffix_label(String(symb), String(pre))


# remove suffix utilities
# -----------------------

function remove_suffix(s::String, suffix::String)
    if endswith(s, suffix)
        return chop(s, tail=length(suffix))
    else
        error("Suffix '$suffix' not found at the end of '$s'")
    end
end

remove_suffix(symb::Symbol, suffix::String) = Symbol(remove_suffix(String(symb), suffix))
remove_suffix(symbs::Tuple, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    Tuple(s ∈ exclude ? s : remove_suffix(s, suffix) for s ∈ symbs)
remove_suffix(symbs::AbstractVector, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    [s ∈ exclude ? s : remove_suffix(s, suffix) for s ∈ symbs]
remove_suffix(d::Dict{Symbol, Any}, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[]) =
    typeof(d)(k ∈ exclude ? k : remove_suffix(k, suffix) => v for (k, v) ∈ d)

function remove_suffix(nt::NamedTuple, suffix::String; exclude::AbstractVector{<:Symbol}=Symbol[])
    symbs = Tuple(k ∈ exclude ? k : remove_suffix(k, suffix) for k ∈ keys(nt))
    return NamedTuple{symbs}(values(nt))
end

function get_suffix(nt::NamedTuple, suffix::String; remove::Bool=false)
    names = Tuple(remove ? remove_suffix(k, suffix) : k for (k, v) ∈ pairs(nt) if endswith(k, suffix))
    values = [v for (k, v) ∈ pairs(nt) if endswith(k, suffix)]
    return NamedTuple{names}(values)
end

function get_suffix(d::Dict{<:Symbol, <:Any}, suffix::String; remove::Bool=false)
    return Dict(remove ? remove_suffix(k, suffix) : k => v for (k, v) ∈ d if endswith(k, suffix))
end

function get_suffix(traj::NamedTrajectory, suffix::String; remove::Bool=false)
    state_names = Tuple(s for s ∈ traj.state_names if endswith(s, suffix))

    # control names
    if traj.timestep isa Symbol
        if endswith(traj.timestep, suffix)
            control_names = Tuple(s for s ∈ traj.control_names if endswith(s, suffix))
            timestep = remove ? remove_suffix(traj.timestep, suffix) : traj.timestep
            exclude = Symbol[]
        else
            # extract the shared timestep
            control_names = Tuple(s for s ∈ traj.control_names if endswith(s, suffix) || s == traj.timestep)
            timestep = traj.timestep
            exclude = [timestep]
        end
    else
        control_names = Tuple(s for s ∈ traj.control_names if endswith(s, suffix))
        timestep = traj.timestep
        exclude = Symbol[]
    end

    component_names = Tuple(vcat(state_names..., control_names...))
    components = get_components(component_names, traj)
    if remove
        components = remove_suffix(components, suffix; exclude=exclude)
    end

    return NamedTrajectory(
        components,
        controls=remove ? remove_suffix(control_names, suffix; exclude=exclude) : control_names,
        timestep=timestep,
        bounds=get_suffix(traj.bounds, suffix, remove=remove),
        initial=get_suffix(traj.initial, suffix, remove=remove),
        final=get_suffix(traj.final, suffix, remove=remove),
        goal=get_suffix(traj.goal, suffix, remove=remove)
    )
end


# =========================================================================== #

@testitem "Apply suffix to trajectories" begin
    include("../test/test_utils.jl")

    traj = named_trajectory_type_1(free_time=false)
    suffix = "_new"
    new_traj = add_suffix(traj, suffix)

    @test new_traj.state_names == add_suffix(traj.state_names, suffix)
    @test new_traj.control_names == add_suffix(traj.control_names, suffix)

    same_traj = add_suffix(traj, "")
    @test traj == same_traj
end

@testitem "Merge systems" begin
    H_drift = 0.01 * GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    T = 50
    sys_1 = QuantumSystem(H_drift, H_drives)
    sys_2 = deepcopy(sys_1)

    # direct sum of systems
    sys_sum = direct_sum(sys_1, sys_2)
    @info sys_sum.n_drives


    @test sys_sum.levels == sys_1.levels * 2
    @test isempty(symdiff(keys(sys_sum.params), [:system_1, :system_2]))

    sys_sum_2 = direct_sum(sys_sum, deepcopy(sys_1))

    @test sys_sum_2.levels == sys_1.levels * 3
    display(sys_sum_2.params)
    @test isempty(symdiff(keys(sys_sum_2.params), [:system_1, :system_2, :system_3]))

end

# TODO: fix broken test
@testitem "Get suffix" begin
    @test_broken false

    # using NamedTrajectories

    # sys = QuantumSystem(0.01 * GATES[:Z], [GATES[:X], GATES[:Y]])
    # T = 50
    # Δt = 0.2
    # ip_ops = IpoptOptions(print_level=1)
    # pi_ops = PiccoloOptions(verbose=false, free_time=false)
    # prob1 = UnitarySmoothPulseProblem(sys, GATES[:X], T, Δt, piccolo_options=pi_ops, ipopt_options=ip_ops)
    # prob2 = UnitarySmoothPulseProblem(sys, GATES[:Y], T, Δt, piccolo_options=pi_ops, ipopt_options=ip_ops)

    # # Direct sum problem with suffix extraction
    # # Note: Turn off control reset
    # direct_sum_prob = UnitaryDirectSumProblem([prob1, prob2], 0.99, drive_reset_ratio=0.0, ipopt_options=ip_ops)
    # # TODO: BROKEN HERE
    # prob1_got = get_suffix(direct_sum_prob, "1")
    # @test prob1_got.trajectory == add_suffix(prob1.trajectory, "1")

    # # Mutate the direct sum problem
    # update!(prob1_got.trajectory, :a1, ones(size(prob1_got.trajectory[:a1])))
    # @test prob1_got.trajectory != add_suffix(prob1.trajectory, "1")

    # # Remove suffix during extraction
    # prob1_got_without = get_suffix(direct_sum_prob, "1", remove=true)
    # @test prob1_got_without.trajectory == prob1.trajectory
end

# TODO: fix broken test
@testitem "Append to default integrators" begin
    @test_broken false
    # sys = QuantumSystem(0.01 * GATES[:Z], [GATES[:Y]])
    # T = 50
    # Δt = 0.2
    # ip_ops = IpoptOptions(print_level=1)
    # pi_false_ops = PiccoloOptions(verbose=false, free_time=false)
    # pi_true_ops = PiccoloOptions(verbose=false, free_time=true)
    # prob1 = UnitarySmoothPulseProblem(sys, GATES[:Y], T, Δt, piccolo_options=pi_false_ops, ipopt_options=ip_ops)
    # prob2 = UnitarySmoothPulseProblem(sys, GATES[:Y], T, Δt, piccolo_options=pi_true_ops, ipopt_options=ip_ops)

    # suffix = "_new"
    # # UnitaryPadeIntegrator
    # # TODO: BROKEN HERE
    # prob1_new = add_suffix(prob1.integrators, suffix)
    # @test prob1_new[1].unitary_symb == add_suffix(prob1.integrators[1].unitary_symb, suffix)
    # @test prob1_new[1].drive_symb == add_suffix(prob1.integrators[1].drive_symb, suffix)

    # # DerivativeIntegrator
    # @test prob1_new[2].variable == add_suffix(prob1.integrators[2].variable, suffix)

    # # UnitaryPadeIntegrator with free time
    # prob2_new = add_suffix(prob2.integrators, suffix)
    # @test prob2_new[1].unitary_symb == add_suffix(prob2.integrators[1].unitary_symb, suffix)
    # @test prob2_new[1].drive_symb == add_suffix(prob2.integrators[1].drive_symb, suffix)

    # # DerivativeIntegrator
    # @test prob2_new[2].variable == add_suffix(prob2.integrators[2].variable, suffix)
end

# TODO: fix broken test
@testitem "Free time get suffix" begin
    # using NamedTrajectories

    # sys = QuantumSystem(0.01 * GATES[:Z], [GATES[:Y]])
    # T = 50
    # Δt = 0.2
    # ops = IpoptOptions(print_level=1)
    # pi_false_ops = PiccoloOptions(verbose=false, free_time=false)
    # pi_true_ops = PiccoloOptions(verbose=false, free_time=true)
    # suffix = "_new"
    # timestep_name = :Δt

    # prob1 = UnitarySmoothPulseProblem(sys, GATES[:Y], T, Δt, piccolo_options=pi_false_ops, ipopt_options=ops)
    # traj1 = direct_sum(prob1.trajectory, add_suffix(prob1.trajectory, suffix), free_time=true)

    # # Direct sum (shared timestep name)
    # @test get_suffix(traj1, suffix).timestep == timestep_name
    # @test get_suffix(traj1, suffix, remove=true).timestep == timestep_name

    # prob2 = UnitarySmoothPulseProblem(sys, GATES[:Y], T, Δt, ipopt_options=ops, piccolo_options=pi_true_ops)
    # traj2 = add_suffix(prob2.trajectory, suffix)

    # # Trajectory (unique timestep name)
    # @test get_suffix(traj2, suffix).timestep == add_suffix(timestep_name, suffix)
    # @test get_suffix(traj2, suffix, remove=true).timestep == timestep_name
end

end # module
