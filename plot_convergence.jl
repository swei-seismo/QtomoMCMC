using JLD, HDF5
using Distributed 
using Glob
@everywhere include("DefStruct.jl")
@everywhere include("define_TDstructure.jl")

@everywhere include("MCsub.jl")
# @everywhere include("load_data_Tonga.jl")
# @everywhere include("TD_inversion_function.jl")


@time TD_parameters = define_TDstructrure()
make_dir(TD_parameters)

function plot_convergence(TD_parameters::Dict{String,Any},start_ind)
    # plot the convergence of nCells and phi value over iterations for each chain
    # saved in ./figures/nCells and ./figures/phi
        ENV["GKSwstype"] = "nul"

        phi_all = []
        nCells_all = []
        p_phi = plot()
        p_nCells = plot()
        for chain in 1:TD_parameters["n_chains"]
            model_checkpoint_lists = glob("./models/chain" * string(chain) * "_*")
            if length(model_checkpoint_lists) == 0
                println("ERROR: Couldn't Find Models For Chain" * string(chain))
            else
                # load the newest model
                split1      = split.(model_checkpoint_lists,"%")
                split2      = [split1[i][1] for i in 1:length(split1)]
                split3      = split.(split2,"_")
                split4      = [split3[i][end] for i in 1:length(split3)]
                model_ind   =  findmax(parse.(Float64,split4))[2]
                load_model  = load(model_checkpoint_lists[model_ind])
                for irm in 1:length(model_checkpoint_lists)
                    if irm == model_ind
                        continue
                    end
                    rm(sort(model_checkpoint_lists)[irm])
                end
                # plot nCells and phi over iterations
                cellnumber_list = load_model["nCells"][start_ind:end]
                phi_list        = load_model["phi"][start_ind:end]
                # push!(phi_all, phi_list)
    
                p1 = plot(1:length(cellnumber_list),cellnumber_list)
                title!(p1,"nCells from" *string(start_ind)* " in Chain" *string(chain))
                xlabel!(p1,"Iterations")
                ylabel!(p1,"Number of Cells")
                p2 = plot(1:length(phi_list),phi_list)
                title!(p2,"Phi from" *string(start_ind)* " in Chain" *string(chain))
                xlabel!(p2,"Iterations")
                ylabel!(p2,"Phi")

                plot!(p_phi,1:length(phi_list),phi_list)
                plot!(p_nCells,1:length(cellnumber_list),cellnumber_list)

                savefig(p1,"./figures/nCells/nCells_from"* string(start_ind) *"_chain" * string(chain))
                savefig(p2,"./figures/phi/phi_from"* string(start_ind) *"_chain" * string(chain))
     
            end
        end
        title!(p_phi,"Phi from" *string(start_ind)* "in all the chains")
        xlabel!(p_phi,"Iterations")
        ylabel!(p_phi,"Phi")
        savefig(p_phi,"./figures/phi/phi_all_from" *string(start_ind)* ".png")

        title!(p_nCells,"nCells from" *string(start_ind)* "in all the chains")
        xlabel!(p_nCells,"Iterations")
        ylabel!(p_nCells,"Number of Cells")
        savefig(p_nCells,"./figures/nCells/nCells_all_from" *string(start_ind)* ".png")

end

@time plot_convergence(TD_parameters,1)
@time plot_convergence(TD_parameters,10000)
@time plot_convergence(TD_parameters,100000)
# @time plot_llh(models)



# @time save("model.jld","model",models)

# model_checkpoint_lists = glob("chain*jld")
# [rm(model_checkpoint_lists[i]) for i in 1:length(model_checkpoint_lists)]
