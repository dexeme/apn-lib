struct ExperimentSpec
    id::Symbol
    description::String
    columns::Vector{Symbol}
    compute::Function
end

function ExperimentSpec(id::Symbol, compute::Function;
                        description::AbstractString = "",
                        columns = Symbol[])
    return ExperimentSpec(id, String(description), collect(Symbol, columns), compute)
end

struct ExperimentContext
    id::Symbol
    catalogue::Catalogue
    description::String
    source::String
    root_dir::String
    artifacts_dir::String
    fixtures::Dict{Symbol, Any}
    parameters::Dict{Symbol, Any}
    expected::Dict{Symbol, Any}
    metadata::Dict{Symbol, Any}
end

function ExperimentContext(id::Symbol, catalogue::Catalogue;
                           description::AbstractString = "",
                           source::AbstractString = "",
                           root_dir::AbstractString = pwd(),
                           artifacts_dir::AbstractString = joinpath(String(root_dir), "artifacts"),
                           fixtures = Dict{Symbol, Any}(),
                           parameters = Dict{Symbol, Any}(),
                           expected = Dict{Symbol, Any}(),
                           metadata = Dict{Symbol, Any}())
    return ExperimentContext(id,
                             catalogue,
                             String(description),
                             String(source),
                             String(root_dir),
                             String(artifacts_dir),
                             Dict{Symbol, Any}(fixtures),
                             Dict{Symbol, Any}(parameters),
                             Dict{Symbol, Any}(expected),
                             Dict{Symbol, Any}(metadata))
end

struct RunMetadata
    started_at::String
    julia_version::String
    threads::Int
    project_path::String
    hostname::String
    environment::Dict{String, String}
end

function RunMetadata(; started_at::AbstractString = string(Dates.now()),
                     project_path::AbstractString = Base.active_project() === nothing ? pwd() : Base.active_project(),
                     hostname::AbstractString = gethostname(),
                     env_keys = String[])
    environment = Dict{String, String}()
    for key in env_keys
        haskey(ENV, key) && (environment[String(key)] = ENV[key])
    end

    return RunMetadata(String(started_at),
                       string(VERSION),
                       Threads.nthreads(),
                       String(project_path),
                       String(hostname),
                       environment)
end

struct ExperimentResult
    spec::ExperimentSpec
    context_id::Symbol
    rows::Vector{Any}
    metadata::RunMetadata
end

catalogue_key(function_::APNFunction) = (function_.n, function_.id)

function expected_values(context::ExperimentContext, experiment_id::Symbol)
    return context.expected[experiment_id]
end

function expected_value(context::ExperimentContext,
                        experiment_id::Symbol,
                        function_::APNFunction;
                        key = catalogue_key(function_))
    return expected_values(context, experiment_id)[key]
end

function fixture(context::ExperimentContext, name::Symbol)
    return context.fixtures[name]
end

function parameter(context::ExperimentContext, name::Symbol)
    return context.parameters[name]
end

function run_experiment(context::ExperimentContext,
                        spec::ExperimentSpec;
                        metadata::RunMetadata = RunMetadata())
    rows = Any[]

    for function_ in context.catalogue.functions
        row = spec.compute(context, function_)
        if row isa AbstractVector
            append!(rows, row)
        else
            push!(rows, row)
        end
    end

    return ExperimentResult(spec, context.id, rows, metadata)
end
