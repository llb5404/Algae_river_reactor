cd("C:/Users/slant/OneDrive/Desktop/Julia_2")
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
    include("MakePlots.jl")             # Plot_Biomass_Profile, Plot_Temperature_Profile

    include("PDE_CO2.jl")
    include("PDE_Height.jl")

    #include("PDE_Water.jl")
    #include("ExportData.jl")

    function Main_PDE!(dX, X, params, t)
        Ny = params.num_odes_y
        Nz = params.num_odes_z
        Nelements = (Ny+1) * (Nz+1)

        M_biomass = X[1:Nelements]
        Temperature = X[Nelements+1:2*Nelements]
        CO2 = X[2*Nelements+1:3*Nelements]
        Height = X[3*Nelements+1:4*Nelements]
    
        # more ....

        PDE_AlgaeBiomass!(dX, M_biomass, Temperature, Height, params, t)      # changes dX
        HeatTransfer!(dX, Temperature, Height, params, t)          # changes dX
        PDE_CO2!(dX, CO2, Height, M_biomass, Temperature, params, t)            # changes dX
        HeightChange!(dX, Height, Temperature, params, t)
        #dX is changed directly by the above functions

        #@show dX
        #@show size(X)
        #@show size(C_biomass)
        #@show size(Temperature)
        #@show size(CO2)
        
        nothing
    end

    function Run()

        filesuffix = "v1"

        ## Load Parameters
        params = LoadDefaultParameters(filesuffix)

        Time_end = params.time_end
        Time_interval = params.time_interval
        Ny = params.num_odes_y
        Nz = params.num_odes_z
        L = params.reactor_length
        H = params.reactor_initial_liquid_level
        W = params.reactor_width                   # m
        dy = L/Ny
        dz = H/Nz

        Nelements = (Ny+1) * (Nz+1)
        pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
        idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]

        C_biomass_in = params.input_biomass_concentration
        Temperature_in = params.input_temperature
        S_co2 = params.solubility_co2_water         #mole fraction of dissolved CO2 in water at equilibrium, depends on Temperature.
        density_water = params.density_water        # kg / m^3, depends on Temperature.
        molecular_weight_co2 = params.molecular_weight_co2
        molecular_weight_water = params.molecular_weight_water
        CO2_in = S_co2(Temperature_in) .* density_water(Temperature_in) .* molecular_weight_co2 ./ molecular_weight_water
        H_init = params.reactor_initial_liquid_level #reactor initial liquid level, m
    

        #ICs: initial conditions at t = 0
        M_biomass_o = zeros( (Ny+1) * (Nz+1), 1)
        M_biomass_o[pos2idx(0,0:Nz)] .= C_biomass_in*(W*dy*dz)
        Temperature_o = ones( (Ny+1) * (Nz+1), 1) .* Temperature_in
        CO2_o = ones( (Ny+1) * (Nz+1), 1) .* CO2_in
        H_o = ones((Ny+1)*(Nz+1), 1).*H_init


        #H_in(x) = H_init(x)

        #for i = 0:Ny
            
            #H_o[pos2idx(i,0:Nz)] .= H_in(i)
            
        #end
        
        tspan = (0.0, Time_end)

        #@show size(C_biomass_o)
        #@show size(Temperature_o)
        @show size(H_o)
        

        Xo = [M_biomass_o; Temperature_o; CO2_o; H_o]
        #@show size(Xo)
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

        H_out = zeros(TL, Nelements)
        for i in 1:TL
            H_out[i,:] .= sol.u[i][3*Nelements+1:4*Nelements]
        end
        
        Plot_Biomass_Profile(Mout, Tout, T, params, filesuffix)
        Plot_Temperature_Profile(Tout, T, params, filesuffix)
        Plot_CO2_Profile(CO2_out, Tout, T, params, filesuffix)
        Plot_Height_Profile(Tout, T, params, filesuffix)

    end

    export Run
end