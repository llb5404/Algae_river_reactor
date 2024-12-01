
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
    include("PDE_CO2.jl")               #Carbon Flux
   
    


    function Main_PDE!(dX, X, params, t)
        Ny = params.num_odes_y
        Nz = params.num_odes_z
        Nelements = (Ny+1) * (Nz+1)

        C_biomass = X[1:Nelements]
        Temperature = X[Nelements+1:2*Nelements]
        CO2 = X[2*Nelements+1:3*Nelements]
        DIC = X[3*Nelements+1:4*Nelements]

        PDE_AlgaeBiomass!(dX, C_biomass, DIC, CO2, Temperature, params, t)      # changes dX
        HeatTransfer!(dX, Temperature, params, t)          # changes dX
        PDE_CO2!(dX, CO2, DIC, C_biomass, Temperature, params, t)            # changes dX

        #dX is changed directly by the above functions
        
        nothing
    end
    
    
    function Run(l, q)
       
        include("LoadParameters.jl")        # LoadDefaultParameters

        ## Load Parameters
        params = LoadDefaultParameters(l, q)

        Time_end = params.time_end

        Ny = params.num_odes_y
        Nz = params.num_odes_z
        Nelements = (Ny+1) * (Nz+1)
        pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
        idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]

        C_biomass_in = params.input_biomass_concentration #g/L
        Temperature_in = params.input_temperature #K
        CO2_in = params.co2_init #g/m3
        DIC_in = params.DIC_init #umol/L, remember this excludes CO2

        Nz = params.num_odes_z
        Ny = params.num_odes_y
     
        #ICs: initial conditions at t = 0
        C_biomass_o = zeros( (Ny+1) * (Nz+1), 1)
        C_biomass_o[pos2idx(0,0:Nz)] .= C_biomass_in
        Temperature_o = ones( (Ny+1) * (Nz+1), 1) .* Temperature_in
        CO2_o = ones( (Ny+1) * (Nz+1), 1) .* 0.4401
        CO2_o[pos2idx(0,1:Nz)] .= CO2_in
        DIC_o = (ones( (Ny+1) * (Nz+1), 1) .* DIC_in)
        tspan = (0.0, Time_end)

        Xo = [C_biomass_o; Temperature_o; CO2_o;DIC_o]
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
        
        return Plot_Biomass_Profile(Mout, CO2_out,DIC_out,Tout, T, params)

    end
    export Run
end