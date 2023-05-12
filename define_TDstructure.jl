function define_TDstructrure() 
    TD_parameters = Dict(
       #====basic parameters====# 
        "debug_prior"       => 0        ::Int64,    # decide whether the data will be used or not
        "plot_voronoi"      => false     ::Bool,     # plot each voronoi diagram in the generated models or not
        "add_yVec"          => 1        ::Int64,    # 0: only work in 2D: x- and z- axis; 1: 3D tomography
        
        #====Voronoi diagram parameters====#
        "sig"               => 10       ::Int64,    # in percent of scales for all parameters.
        "zeta_scale"        => 50       ::Int64,    # std of the prior for Normal. +/- bounds for uniform. 
        "max_cells"         => 250      ::Int64,    # don't go over 100
        "min_cells"         => 5        ::Int64,    # minimum for a box shaped anomaly
        "max_sig"           => 0.1      ::Float64,  # on t*
        "interp_style"      => 1        ::Int64,    # nearest:1, Inverse Distance Weighting:2. Similar results for each
        "enforce_discon"    => 0        ::Int64,    # load a discontinuity, and use it in the inversion
        "prior"             => 1        ::Int64,    # Uniform:1, normal:2, or exponential:3. exponential not recommended
    # check whether they are used
        "event_statics"     => 1        ::Int64,    # much better if you use this if using relative measurements. 
        "demean"            => 1        ::Int64,    # for plotting only.  
    ### 

        #=====Monte Carlo parameters=====#
        # "n_chains"          => 20       ::Int64,    # Per script, not in total. See iter field in run_batch.pbs
        # "n_iter"            => 1e6      ::Float64,  # How many times you iterate each chain
        # "burn_in"           => 5e5      ::Float64,  # Don't save until you get to this iteration
        # "keep_each"         => 5e3      ::Float64,  # how often to save a model post burnin
        # "print_each"        => 1e5      ::Float64,  # how often to print to the screen if you are impatient and chekcing it (like I do)
        # "save_percent"      => 2        ::Int64,
        # "n_chains"          => 5        ::Int64,    # Per script, not in total. See iter field in run_batch.pbs
        # "n_iter"            => 1e3      ::Float64,  # How many times you iterate each chain
        # "burn_in"           => 5e2      ::Float64,  # Don't save until you get to this iteration
        # "keep_each"         => 1e1      ::Float64,  # how often to save a model post burnin
        # "print_each"        => 1e2      ::Float64,  # how often to print to the screen if you are impatient and chekcing it (like I do)
        # "save_percent"      => 5        ::Int64,
        
         "n_chains"          => 10       ::Int64,    # Per script, not in total. See iter field in run_batch.pbs
         "n_iter"            => 1e5      ::Float64,  # How many times you iterate each chain
         "burn_in"           => 5e4      ::Float64,  # Don't save until you get to this iteration
        "keep_each"         => 1e2      ::Float64,  # how often to save a model post burnin
        "print_each"        => 1e4     ::Float64,  # how often to print to the screen if you are impatient and chekcing it (like I do)        
        "save_percent"      => 5       ::Int64,

#        "n_chains"          => 5       ::Int64,    # Per script, not in total. See iter field in run_batch.pbs
 #       "n_iter"            => 5e5      ::Float64,  # How many times you iterate each chain
  #      "burn_in"           => 1e5      ::Float64,  # Don't save until you get to this iteration
   #     "keep_each"         => 1e2      ::Float64,  # how often to save a model post burnin
        # currently 7e5/1e2 = 7e3 models
    #    "print_each"        => 2e4      ::Float64,  # how often to print to the screen if you are impatient and chekcing it (like I do)        
     #   "save_percent"      => 5       ::Int64,

        #=====map parameters=====#
        "max_depth"         => 660.0    ::Float64,  # in km
        "min_depth"         => 0.0      ::Float64,  # in km, should just be zero
    #check whether it's used
        "rotation"          => 20       ::Int64,    # used to rotate the grid. Helpful when using linear arrays so everthing is on one axis
    ###
        "ZnodeSpacing"      => 20       ::Int64,    #10   # in km, how often in the z-direction do you grid the ray paths
        "buffer"            => 100      ::Int64,    # in km. If zero, edges of the modeling domain are [max_x min_x max_y min_y] station locations. Add buffer to broader the imaging domain. 
        "XYnodeSpacing"     => 20       ::Int64,    #5  # not actually part of the inversion, but used for synthetic models.
        
        #====Plot Cross section parameters====#
        "cmax"              => 20.0     ::Float64,  # max 1000/Q on colorbar
        "xyMap"             => true     ::Bool,     # work on x-y cross section with a given z value
        "z0"                => [50,100,300,500]     ::Array{Int64,1},     # the given z value (vertical axis)
        "xzMap"             => true     ::Bool,      # work on x-z cross section with a given y value
        "y0"                => [550,700,850]        ::Array{Int64,1}     # the given y value (horizontal axis, E-W afer rotation of β)
        )
    
    return TD_parameters

end
