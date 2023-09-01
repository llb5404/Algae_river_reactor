

function Plot_Biomass_Profile(Mout, Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    Q(H,T) = params.volumetric_flow_rate(H,T)       # m^3/hour
    Re(H,T) = params.reynolds_number(H,T)
    V_prof = params.velocity_profile
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    dy = L/Ny
    
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
   
    
    Ieff(Cz, z, Hgt) = (1 - min(sum(Cz[1:z]) * Hgt / Nz / Cmax, 1.0))
    phiL(Cz, zlist, Hgt) = [Ieff(Cz, z, Hgt) * exp(1 - Ieff(Cz, z, Hgt)) for z in zlist]
    
    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end
    T_days = T ./ 24.0
    

    #V_profile(x,zlist,H) =([V_prof(x,z*(H/Int(Nz_max),H)) for z in zlist])/3600
    #velocity_profile(x,z,H) = max_flow_velocity(x)*(1-(z/H)) #z = dz(x)*z

    H_o = params.reactor_initial_liquid_level
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb)
    Pa_out = zeros(TL, 1)
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i])
    end

    Sal = params.salinity
    Sout = zeros(TL, Ny+1)

    for i in 1:TL
        for j in 0:Ny
            Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],Tambout[i],WNDSPDout[i],RHout[i],Pa_out[i],j)
        end
    end
    


    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o)
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz

    Cout = zeros(TL, Nelements)

    for i= 1:TL
        for j = 0:Ny
            for k = 0:Nz
                Cout[i,pos2idx(j,k)] = Mout[i,pos2idx(j,k)]/(dy*W*dz(j,Tout[i,pos2idx(j,k)],0,WNDSPDout[i],RHout[i],Pa_out[i])) #kg/m3
            end
        end
    end

    Prod_out = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                Prod_out[i,pos2idx(j,k)] = ((Cout[i,pos2idx(j,k)]- Cinit)* Q(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])* 24.0* 1000)/(W*L)
            end
        end
    end

    V_prof_out = zeros(TL,Nelements)
    for i = 1:TL
        for j = 0:Ny
            for k = 0:Nz
                V_prof_out[i,pos2idx(j,k)] = V_prof(k*dz(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tout[i,pos2idx(j,k)])
            end
        end
    end
    
    Re_out = zeros(TL, Ny+1)
    for i = 1:TL
        for j = 0:Ny
            Re_out[i,pos2idx(j,0)] = Re(Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i]),Tambout[i])
        end
    end

    (maxval, maxpos) = findmax(Cout[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)

    P(C) = Statistics.mean(C[:, pos2idx(Ny,0:Nz)], dims = 2)
    Q(C) = Statistics.mean(C[:, pos2idx(0:Ny,0)], dims = 2)

    #P(C) = (Statistics.mean(C[:, pos2idx(Ny,0:Nz)], dims=2) .- Cinit).* Q ./ L ./ W .* 24.0 .* 1000; # mass flow rate [g/day] / surface area [m^2]



    p1 = plot(Y,Cout[TL, pos2idx(0:Ny,0)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Surface Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,Cout[TL, pos2idx(0:Ny,Nz)] * 1000.0, xlabel = "Length [m]", ylabel = "Algae Cell Density (g/m^3)", title="Algae Floor Cell Density", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(Z,V_prof_out[TL, pos2idx(0,0:Nz)]/3600.0, xlabel = "Height [m]", ylabel = "Fluid Velocity (m/s)", title="Velocity Profile at Inlet", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T,Q(Re_out), xlabel = "Time [hours]", ylabel = "Reynold's #", title="Average Reynold's # Over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T, P(Prod_out), xlabel = "Time [hours]", ylabel = "Algae Productivity (g/m^2/day)", title = "Net Continuous Biomass Productivity", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p6 = plot(T, GHIout, xlabel = "Time [hours]", ylabel = "Global Horizontal Irradiance (GHI) [W/m^2]", title = "Solar Energy", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, p6, layout=(6,1), legend=false, size=(1200,1200))
    png("Biomass_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Biomass_AlgaeRiverReactor_$filesuffix.ps")
    Average_Continuous_Productivity = Statistics.mean(P(Mout)[max(1,TL-20):TL])
    # Mout2D(z, y), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Cout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Cout2D[j+1,i+1] = Cout[TL, pos2idx(i,j)]
        end
    end
    q = heatmap(Y, Z, log10.(Cout2D * 1000.0),
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Algae Cell Density (g/L)",
            size=(800,400)
            )
    savefig(q, "BiomassProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("BiomassProfile_AlgaeRiverReactor_$filesuffix")
end
function Plot_Temperature_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin
    TL = length(T)
    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end
    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(Tout[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)
    O(T) = Statistics.mean(T[:, pos2idx(Ny-1,0:Nz)], dims=2); #avg inlet
    I(T) = Statistics.mean(T[:, pos2idx(1,0:Nz)], dims=2); #avg outlet
    S(T) = Statistics.mean(T[:,pos2idx(0:Ny,0)], dims=2); #avg surface
    F(T) = Statistics.mean(T[:,pos2idx(0:Ny,Nz)], dims=2); #avg floor
    D(T) = O(T)-I(T)
    p1 = plot(T, S(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Surface Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T, F(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Floor Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T, O(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Outlet Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T, I(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Inlet Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T, D(Tout), xlabel = "Time [hours]", ylabel = "Temperature (K)", title = "Average Temperature Difference (Outlet - Inlet)", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p6 = plot(T, GHIout, xlabel = "Time [hours]", ylabel = "Global Horizontal Irradiance (GHI) [W/m^2]", title = "Solar Energy", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p7 = plot(T, WNDSPDout, xlabel = "Time [hours]", ylabel = "Wind Speed (WNDSPD) [m/s]", title = "Wind Speed", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p8 = plot(T, RHout, xlabel = "Time [hours]", ylabel = "Relative Humidity (RH) [%]", title = "Relative Humidity", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p9 = plot(T, Tambout, xlabel = "Time [hours]", ylabel = "Ambient Temperature (Tamb) [K]", title = "Ambient Temperature", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(9,1), legend=false, size=(1600,1600))
    png("Thermal_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Thermal_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Tout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Tout2D[j+1,i+1] = Tout[TL, pos2idx(i,j)]
        end
    end
    q = heatmap(Y, Z, Tout2D,
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Temperature (K)",
            size=(800,400)
            )
    savefig(q, "TemperatureProfile_AlgaeRiverReactor_$filesuffix.ps")
    png("TemperatureProfile_AlgaeRiverReactor_$filesuffix")
end
function Plot_CO2_Profile(CO2_out, Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)

    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end

    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
    dy = L/Ny

    Ieff(Cz, z, Hgt) = (1 - min(sum(Cz[1:z]) * Hgt / Nz / Cmax, 1.0))
    phiL(Cz, zlist,Hgt) = [Ieff(Cz, z, Hgt) * exp(1 - Ieff(Cz, z, Hgt)) for z in zlist]
    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(CO2_out[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)



    H_o = params.reactor_initial_liquid_level

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb)
    Pa_out = zeros(TL, 1)
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i])
    end

    Sal = params.salinity
    Sout = zeros(TL, Ny+1)

    for i in 1:TL
        for j in 0:Ny
            Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],Tambout[i],WNDSPDout[i],RHout[i],Pa_out[i],j)
        end
    end

    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o)
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz

    CO2C_out = zeros(TL, Nelements)

    for i= 1:TL
        for j = 0:Ny
            for k = 0:Nz
                CO2C_out[i,pos2idx(j,k)] = CO2_out[i,pos2idx(j,k)]/(dy*W*(Hght(j,Tout[i,pos2idx(j,k)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])/Nz))
            end
        end
    end
    p1 = plot(Y,CO2C_out[TL, pos2idx(0:Ny,0)] * 1000.0, xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Surface Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(Y,CO2C_out[TL, pos2idx(0:Ny,Nz)] * 1000.0, xlabel = "Length [m]", ylabel = "Dissolved CO2 (g/m^3)", title="CO2 Floor Concentration", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("CO2_AlgaeRiverReactor_$filesuffix")
    savefig(p, "CO2_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    Cout2D = zeros(Nz+1, Ny+1)
    for i in 0:Ny
        for j in 0:Nz
            Cout2D[j+1,i+1] = max(0, CO2C_out[TL, pos2idx(i,j)]) # don't allow negative values
        end
    end
    q = heatmap(Y, Z, log10.(Cout2D * 1000.0),
            yflip=true,
            c=cgrad([:blue, :white,:red, :yellow]),
            xlabel="Length (0 is entry) [m]", ylabel="Liquid Level (0 is surface) [m]",
            title="Dissolved CO2 (g/m^3)",
            size=(800,400)
            )
    savefig(q, "CO2Profile_AlgaeRiverReactor_$filesuffix.ps")
    png("CO2Profile_AlgaeRiverReactor_$filesuffix")
end

function Plot_Height_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)

    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end
    
    Evap_C = zeros(TL, Ny+1)
    X_air = zeros(TL, Ny+1)
    X_surface = zeros(TL, Ny+1)
    M_Evap = zeros(TL, Ny+1)
    T_new = zeros(TL,Ny+1)
    Cum_Evap_Loss = zeros(TL, Ny+1)

    P_atm = params.reference_pressure
    P_w = params.saturated_vapor_pressure_water_air

    H_o = params.reactor_initial_liquid_level

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb)

    Pa_out = zeros(TL, 1)
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i])
    end
    
    for i in 1:TL
        Evap_C[i,pos2idx(0:Ny,0)] .= 25 + 19*WNDSPDout[i]
        X_air[i,pos2idx(0:Ny,0)] .= 0.62198*(P_w(Tambout[i])*(RHout[i]/100))/(P_atm-(P_w(Tambout[i])*(RHout[i]/100)))
    end
    for i in 1:TL
        T_new[i,pos2idx(0:Ny,0)] .= T[i]
    end

    for i in 1:TL
        for j in 0:Ny
           
            X_surface[i,pos2idx(j,0)] = 0.62198*(P_w(Tout[i,pos2idx(j,0)])/(P_atm-P_w(Tout[i,pos2idx(j,0)]))) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j,0)] = Evap_C[i,pos2idx(j,0)]*(X_surface[i,pos2idx(j,0)]-X_air[i, pos2idx(j,0)])
            Cum_Evap_Loss[i,pos2idx(j,0)] = W*L*sum(M_Evap[1:i,pos2idx(j,0)]) #kg
        end
    end

    Sal = params.salinity
    Sout = zeros(TL, Ny+1)

    for i in 1:TL
        for j in 0:Ny
            Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],Tambout[i],WNDSPDout[i],RHout[i],Pa_out[i],j)
        end
    end

    K = params.dVavgdx
    rho_sol = params.density_solution
    Kterm_out = zeros(TL, Ny+1)
    for i in 1:TL
        for j in 0:Ny
            Kterm_out[i,pos2idx(j,0)] = K(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],RHout[i],WNDSPDout[i],P_w(Tambout[i]))*rho_sol(Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)])
        end
    end
    

    Hght(x,T,S,WNDSPD,RH,P) = params.height(x,T,S,WNDSPD,RH,P,H_o)
    dz(x,T,S,WNDSPD,RH,P) = Hght(x,T,S,WNDSPD,RH,P)/Nz


    Hght_out = zeros(TL,Ny+1)

    for i in 1:TL
        for j in 0:Ny

            Hght_out[i,pos2idx(j,0)] = Hght(j,Tout[i,pos2idx(j,0)],Sout[i,pos2idx(j,0)],WNDSPDout[i],RHout[i],Pa_out[i])
        end
    end
    
  
    T_days = T ./ 24.0
    (maxval, maxpos) = findmax(Hght_out[TL,:])
    (Ny_max, Nz_max) = idx2pos(maxpos)
    S(Hgt) = Statistics.mean(Hgt[:,pos2idx(0:Ny,0)], dims=2)
    Q(M) = Statistics.mean(M[:,pos2idx(0:Ny,0)], dims=2);
  
    p1 = plot(Y,Hght_out[TL, pos2idx(0:Ny,0)], xlabel = "Length [m]", ylabel = "Height [m]", title="Height vs Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T,S(Hght_out), xlabel = "Time [hours]", ylabel = "Average Height [m]", title="Average height over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p3 = plot(T,Q(M_Evap), xlabel = "Time [hours]", ylabel = "Average Evaporation Rate [kg/m2-hr]", title="Average evaporation rate over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p4 = plot(T,Q(Kterm_out), xlabel = "Time [hours]", ylabel = "Average Counter-Evaporation Rate [kg/m2-hr]", title="Average counter-evaporation rate over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p5 = plot(T,Q(Cum_Evap_Loss), xlabel = "Time [hours]", ylabel = "Average Cumulative Water Loss [kg]", title="Average cumulative water loss over time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p = plot(p1, p2, p3, p4, p5, layout=(5,1), legend=false, size=(1200,1200))
    png("Height_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Height_AlgaeRiverReactor_$filesuffix.ps")
    # Cout(y,z), z = 0 is the surface and z = H is the bottom. For plotting, we will invert the z-scale so that z = 0 is the bottom and z = H is the surface.
    
end

function Plot_Salinity_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    TL = length(T)
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    V_prof = params.velocity_profile
    L = params.reactor_length                   # m
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg/m^3
    Cmax = params.max_biomass_concentration     # kg/m^3
    dy = L/Ny
    
    Y = LinRange(0,L,Ny+1)
    Z = LinRange(0,H,Nz+1)
   

    
    GHIout = zeros(TL,1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        GHIout[i] = GHI_Data[data_begin + t_hour1] * (T[i]-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - T[i])
    end
    WNDSPDout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        WNDSPDout[i] = WNDSPD_Data[data_begin + t_hour1]*(T[i]-t_hour1)+WNDSPD_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    RHout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        RHout[i] = RH_Data[data_begin + t_hour1]*(T[i]-t_hour1)+RH_Data[data_begin + t_hour2]*(t_hour2-T[i])
    end
    Tambout = zeros(TL, 1)
    for i in 1:TL
        t_hour1 = floor(Int64, T[i])
        t_hour2 = floor(Int64, T[i])+1
        Tambout[i] = (Tamb_Data[data_begin + t_hour1]*(T[i]-t_hour1)+Tamb_Data[data_begin + t_hour2]*(t_hour2-T[i]))
    end
    T_days = T ./ 24.0
    

    #V_profile(x,zlist,H) =([V_prof(x,z*(H/Int(Nz_max),H)) for z in zlist])/3600
    #velocity_profile(x,z,H) = max_flow_velocity(x)*(1-(z/H)) #z = dz(x)*z

    H_o = params.reactor_initial_liquid_level
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb)
    Pa_out = zeros(TL, 1)
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i])
    end

    Sal = params.salinity
    Sout = zeros(TL, Ny+1)

    for i in 1:TL
        for j in 0:Ny
            Sout[i,pos2idx(j,0)] = Sal(Tout[i,pos2idx(j,0)],Tambout[i],WNDSPDout[i],RHout[i],Pa_out[i],j)
        end
    end
    

    P(S) = Statistics.mean(S[:,pos2idx(0:Ny,0)], dims = 2)


    p1 = plot(Y,Sout[TL, pos2idx(0:Ny,0)], xlabel = "Length [m]", ylabel = "Salinity (kg/m^3)", title="Reactor Salinity with Length", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)
    p2 = plot(T, P(Sout), xlabel = "Time [hours]", ylabel = "Average Salinity (kg/m^3)", title = "Average Reactor Salinity over Time", plot_titlefontsize=8, labelfontsize=7,tickfontsize=6, grid = false)

    p = plot(p1, p2, layout=(2,1), legend=false, size=(1200,1200))
    png("Salinity_AlgaeRiverReactor_$filesuffix")
    savefig(p, "Salinity_AlgaeRiverReactor_$filesuffix.ps")
end