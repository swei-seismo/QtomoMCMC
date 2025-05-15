using JLD, HDF5
using Distributed 
using Glob,YAML, Plots
@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

file_name = "TongaLau_Spherical.yml"
@time par = load_par_from_yml(file_name)
@time (dataStruct, RayTraces) = load_data_Tonga(par)
make_dir(par, file_name)
println("--------Data Loaded-------")

(m, n) = size(RayTraces["rayX"])
observed_traveltime = zeros(n)
predicted_traveltime = zeros(n)

Threads.@threads for i = 1:n
    rayl = RayTraces["rayL"][:,i]
    index = findfirst(isnan, rayl)
    if index == nothing
        predicted_traveltime[i] = sum(rayl .* RayTraces["rayU"][:,i])
    else
        predicted_traveltime[i] = sum(rayl[1:index-1] .* RayTraces["rayU"][1:index-1,i])
    end
    observed_traveltime[i] = 1000 * dataStruct["tS"][i] / dataStruct["allaveatten"][i]
end

# plot predicted traveltime vs observed traveltime with 1:1 line
plot(observed_traveltime, predicted_traveltime, seriestype = :scatter, label = "Predicted vs Observed", xlabel = "Observed Traveltime (s)", ylabel = "Predicted Traveltime (s)",legend = false)
# maxval = maximum([observed_traveltime; predicted_traveltime])
maxval = 130
plot!([0:maxval], [0:maxval], color = :black, linestyle = :dash)
savefig("traveltime_check.png")