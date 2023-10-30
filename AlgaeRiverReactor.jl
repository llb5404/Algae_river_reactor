
module AlgaeRiverReactor
    __precompile__()

    using DifferentialEquations, ModelingToolkit
    using CSV, BeepBeep, Logging, Interpolations, Printf, Tables
    using DataFrames
    #using JDL2
    using Plots; gr()
    import Statistics
    using Debugger

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
        
    
        # more ....

        PDE_AlgaeBiomass!(dX, C_biomass, CO2, Temperature, params, t)      # changes dX
        HeatTransfer!(dX, Temperature, params, t)          # changes dX
        PDE_CO2!(dX, CO2, C_biomass, Temperature, params, t)            # changes dX
   



        #dX is changed directly by the above functions
        
        nothing
    end

    function Run()
        #Units for params can be found in LoadParameters.jl
        filesuffix = "v1"

        ## Load Parameters
        params = LoadDefaultParameters(filesuffix)

        Time_end = params.time_end
        Time_interval = params.time_interval
        Ny = params.num_odes_y
        Nz = params.num_odes_z

    

        Nelements = (Ny+1) * (Nz+1)
        pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
        idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]

        C_biomass_in = params.input_biomass_concentration
        Temperature_in = params.input_temperature
        S_co2 = params.solubility_co2_water         #mole fraction of dissolved CO2 in water at equilibrium, depends on Temperature.
        density_water = params.density_water        # kg / m^3, depends on Temperature.
        molecular_weight_co2 = params.molecular_weight_co2
        molecular_weight_water = params.molecular_weight_water
        S_in = params.salinity_in
        CO2_in = params.initial_co2_g(Temperature_in,S_in)
        #S_co2(Temperature_in) .* density_water(Temperature_in) .* molecular_weight_co2 ./ molecular_weight_water
       

        
    

        #ICs: initial conditions at t = 0
        C_biomass_o = zeros( (Ny+1) * (Nz+1), 1)
        C_biomass_o[pos2idx(0,0:Nz)] .= C_biomass_in
        Temperature_o = ones( (Ny+1) * (Nz+1), 1) .* Temperature_in
        CO2_o = (ones( (Ny+1) * (Nz+1), 1) .* CO2_in)
        tspan = (0.0, Time_end)

    
        

        Xo = [C_biomass_o; Temperature_o; CO2_o]
        @show tspan

        prob = ODEProblem(Main_PDE!, Xo, tspan, params)
        sol = solve(prob, saveat=8.0)

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

        
        Plot_Biomass_Profile(Mout, CO2_out, Tout, T, params, filesuffix)
        Plot_Temperature_Profile(Tout, T, params, filesuffix)
        Plot_CO2_Profile(CO2_out, Tout, T, params, filesuffix)
        Plot_Height_Profile(Tout, T, params, filesuffix)
        Plot_Salinity_Profile(Tout, T, params, filesuffix)

    end

    export Run
end