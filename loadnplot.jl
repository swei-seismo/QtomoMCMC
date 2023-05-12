using HDF5, JLD
using Plots, Distributed

@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")
@everywhere include("MCsub.jl")
@everywhere include("load_data_Tonga.jl")

ENV["GKSwstype"] = "nul"

cmax = 30.0

FigDirLst = ["figures", "voronoi", "./figures/nCells",
        "./figures/xzMap", "./figures/xzUncertainty", "./figures/xzMask", "./figures/xzVoronoi", "./figures/xzMean",
        "./figures/xyMap", "./figures/xyUncertainty", "./figures/xyMask", "./figures/xyVoronoi", "./figures/xyMean",
        "./figures/phi", "./figures/nCells"
            ]

for iDir in FigDirLst
    if isdir(iDir) == false
        mkdir(iDir)
    end
end

TD_parameters = define_TDstructrure()
@time dataStruct = load_data_Tonga(TD_parameters)
println("--------data loaded-------")

num_models_per_chain = Int64((TD_parameters.n_iter - TD_parameters.burn_in) / TD_parameters.keep_each)
models = load("model.jld")["model"]
println("--------models loaded--------")

plot_model_hist(models, dataStruct, TD_parameters, cmax)

# plot nCells/phi vs time
plotnCells = false
if plotnCells == true
    for i in 1:length(models)
        nCells  = []
        phi     = []
        for j in 1:length(models[i])
            append!(nCells,models[i][j].nCells)
            append!(phi,models[i][j].phi)
        end
        p1 = plot(1:length(models[i]),nCells)
        title!(p1,"nCells of saved models in Chain" *string(i))
        xlabel!(p1,"iteration")
        ylabel!(p1,"nCells")
        p2 = plot(1:length(models[i]),phi)
        title!(p2,"phi of saved models in Chain" *string(i))
        xlabel!(p2,"iteration")
        ylabel!(p2,"phi")
        savefig(p1,"./figures/nCells/nCells_chain" * string(i))
        savefig(p2,"./figures/nCells/phi_chain" * string(i))
    end
end

# plot all the nCells/phi vs iterations
makeplot = true
if makeplot == true
    for i in 1:length(models)
        nCells  = []
        phi     = []
        for j in 1:length(models[i])
            append!(nCells,models[i][j].nCells)
            append!(phi,models[i][j].phi)
        end
        p1 = plot(1:length(models[i]),nCells)
        title!(p1,"nCells of saved models in Chain" *string(i))
        xlabel!(p1,"iteration")
        ylabel!(p1,"nCells")
        p2 = plot(1:length(models[i]),phi)
        title!(p2,"phi of saved models in Chain" *string(i))
        xlabel!(p2,"iteration")
        ylabel!(p2,"phi")
        savefig(p1,"./figures/nCells/nCells_chain" * string(i))
        savefig(p2,"./figures/nCells/phi_chain" * string(i))
    end
end
