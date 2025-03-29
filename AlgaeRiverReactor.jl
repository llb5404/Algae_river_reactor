
module AlgaeRiverReactor
    __precompile__()

    using DifferentialEquations, OrdinaryDiffEq, Sundials
    using CSV, BeepBeep, Logging, Interpolations, Printf, Tables, Debugger, LinearAlgebra
    using DataFrames
    #using JDL2
    using Plots; gr()
    using SparseArrays
    import Statistics


    break_on(:error)
   
   

    # include is used to import Julia files (effectively like copy/pasting the code in place)

    include("LoadParameters.jl")        # LoadDefaultParameters
    include("PDE_AlgaeBiomass.jl")      # PDE_AlgaeBiomass!
    include("PDE_Temperature.jl")       # HeatTransfer!
    include("MakePlots.jl")             # Plot_Biomass_Profile, Plot_Temperature_Profile, etc.
    include("PDE_CO2.jl")               #Carbon Flux\
    include("PDE_Flowrate.jl")          # Flowrate
   
    


    function Main_PDE!(dX, X, params, t)
        Nx = params.num_odes_x
        Ny = params.num_odes_y
        Nelements = (Nx+1)

        C_biomass = X[1:Nelements]
        Temperature = X[Nelements+1:2*Nelements]
        CO2 = X[2*Nelements+1:3*Nelements]
        DIC = X[3*Nelements+1:4*Nelements]
        A =   X[4*Nelements+1:5*Nelements]

        PDE_AlgaeBiomass!(dX, C_biomass, DIC, CO2, Temperature, A, params, t)      # changes dX
        HeatTransfer!(dX, Temperature, A, params, t)          # changes dX
        PDE_CO2!(dX, CO2, DIC, C_biomass, Temperature, A, params, t)            # changes dX
        PDE_Flowrate!(dX, A, Temperature, params, t)

        #dX is changed directly by the above functions
        
        nothing
    end
    
    
    function Run(l, q)
       
        include("LoadParameters.jl")        # LoadDefaultParameters

        ## Load Parameters
        params = LoadDefaultParameters(l, q)

        Time_end = params.time_end

        Nx = params.num_odes_x
        Ny = params.num_odes_y
        Nelements = (Nx+1)
        pos2idx(y) = (y.+1)
        idx2pos(pos) = [Integer(pos - 1 - (Nx+1) * floor( (pos-1) ./ (Nx.+1))), Integer(floor( (pos-1) ./ (Nx.+1)))]

        C_biomass_in = params.input_biomass_concentration #g/L
        Temperature_in = params.input_temperature #K
        CO2_in = params.co2_init #g/m3
        DIC_in = params.DIC_init #umol/L, remember this excludes CO2
        A_in = params.A_o

        Ny = params.num_odes_y
        Nx = params.num_odes_x
     
        #ICs: initial conditions at t = 0
        C_biomass_o = ones( (Nx+1))*C_biomass_in
        Temperature_o = ones( (Nx+1)) .* Temperature_in
        CO2_o = ones((Nx+1)) .* CO2_in
        CO2_o[pos2idx(0)] = CO2_in
        DIC_o = (ones((Nx+1)) .* DIC_in)
        A_o = ones(Nx+1)*A_in

        tspan = (0.0, Time_end)

        Xo = [C_biomass_o; Temperature_o; CO2_o;DIC_o;A_o]
        @show tspan
    
        prob = ODEProblem(Main_PDE!, Xo, tspan, params)
        
        @time sol = DifferentialEquations.solve(prob, CVODE_BDF(linear_solver=:GMRES), saveat=1.0, progress = true)
        
        T = sol.t
        TL = length(sol.t)

        Mout = zeros(TL, Nelements)
        for i in 1:TL
            Mout[i,:] .= sol.u[i][1:Nelements]
        end
        
        Tout = zeros(TL, Nelements)
        for i in 1:TL
            Tout[i,:] .= sol.u[i][Nelements+1:2*Nelements]
        end

        CO2_out = zeros(TL, Nelements)

        for i in 1:TL
            CO2_out[i,:] .= sol.u[i][2*Nelements+1:3*Nelements]
          
        end

        DIC_out = zeros(TL, Nelements)

        for i in 1:TL
            DIC_out[i,:] .= sol.u[i][3*Nelements+1:4*Nelements]
            
        end

        A_out = zeros(TL, Nelements)

        for i in 1:TL
            A_out[i,:] .= sol.u[i][4*Nelements+1:5*Nelements]
            
        end

        filesuffix = "v_v1"
        
        return Plot_Biomass_Profile(Mout, CO2_out,DIC_out,Tout, A_out,T, params, filesuffix)

    end
    export Run
end