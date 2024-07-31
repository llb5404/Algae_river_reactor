
using CSV, Tables

function Plot_Biomass_Profile(Mout, CO2_out, DIC_out,Tout, T, params, filesuffix)
    DIC_out = max.(DIC_out,1E-09)
    DIC_out = min.(DIC_out,70000)

    Ny = params.num_odes_y
    Nz = params.num_odes_z
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny + 1
    pos2idx1(z) = z.+1
    tpos2idx(t,y) = Nelements1*t + pos2idx1(y)
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    tpos2idx2(t,y,z) = Nelements*t + pos2idx(y,z)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    
    
    TL = length(T)
    NelementsT2 = (Nelements)*(TL+1)
  
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin
    bm = params.max_biomass_specific_growth_rate
 
    dy = L/Ny
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Cinit = params.input_biomass_concentration  # kg
    
    Y = LinRange(0,dy*Ny_new,Ny_new+1)
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
    
    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) #Pa
    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #pressure at ambient temp, Pa
    end
   
    Q(T,WNDSPD,RH,P,x) = params.volumetric_flow_rate(T,WNDSPD,RH,P,x)      # m^3/hour v

  
  

    #Initialize vectors
   
    Prod_vec = zeros(TL,Nz+1) #g/m2/day

    for i in 1:TL
            for k in 0:Nz
                    Prod_vec[i,pos2idx1(k)] = ((Mout[i,pos2idx(Ny_new,k)]*Q(Tout[i,pos2idx(Ny_new,0)],WNDSPDout[i],RHout[i],Pa_out[i],Ny_new) - Cinit*params.volumetric_flow_rate_o)* 24.0* 1000)/(W*L)
            end
    end


    P(C) = Statistics.mean(C[:, pos2idx1(1:Nz)], dims = 2) #Used for average productivity
    
    @show Average_Continuous_Productivity = Statistics.mean(P(Prod_vec)[max(1,TL-100):TL])
   

    return Average_Continuous_Productivity
end

function Plot_Height_Profile(Tout, T, params, filesuffix)
    Ny = params.num_odes_y
    L = params.reactor_length                   # m
    RL = params.real_length #e.g. real length is 1/10 of reactor length
    length_frac = RL/L
    Ny_new = trunc(Int,(length_frac)*Ny)
    Nz = params.num_odes_z
    Nelements = (Ny+1) * (Nz+1)
    Nelements1 = Ny+1
    pos2idx(y,z) = (y.+1) .+ z.*(Ny.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Ny+1) * floor( (pos-1) ./ (Ny.+1))), Integer(floor( (pos-1) ./ (Ny.+1)))]
    pos2idx1(z) = z.+1
    tpos2idx(t,y) = Nelements1*t + pos2idx1(y)
    
    TL = length(T)

    NelementsT = (Ny+1)*(TL+1)
    

    dy = L/Ny
    W = params.reactor_width                    # m
    H = params.reactor_initial_liquid_level     # m
    Y = LinRange(0,Ny_new*dy,Ny_new+1)
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
    
    Evap_C = zeros(TL, Ny+1) #evaporation constant, dimless
    X_air = zeros(TL, Ny+1) #ambient humidity ratio, dimless kg H2O/kg dry air
    X_surface = zeros(TL, Ny+1) #surface humidity ratio, dimless
    M_Evap = zeros(TL, Ny+1) #evaporation mass flux, kg H2O/m2/hr
    T_new = zeros(TL,Ny+1) #time with position index, hr
    Cum_Evap_Loss_vec = zeros(TL, Ny_new+1) #cumulative loss of water due to evaporation, kg
    Cum_Evap_Loss = zeros(TL+1)

    P_atm = params.reference_pressure #Pa
    P_w = params.saturated_vapor_pressure_water_air #Pa

    P_a(Tamb) = params.saturated_vapor_pressure_water_air(Tamb) # Pa

    Pa_out = zeros(TL, 1) #Pa
    for i in 1:TL
        Pa_out[i] = P_a(Tambout[i]) #Pa
    end
    
    for i in 1:TL
        Evap_C[i,pos2idx(0:Ny,0)] .= 25 + 19*WNDSPDout[i]
        X_air[i,pos2idx(0:Ny,0)] .= 0.62198*(P_w(Tambout[i])*(RHout[i]/100))/(P_atm-(P_w(Tambout[i])*(RHout[i]/100)))
    end
   

    for i in 1:TL
        for j in 0:Ny
           
            X_surface[i,pos2idx(j,0)] = 0.62198*(P_w(Tout[i,pos2idx(j,0)])/(P_atm-P_w(Tout[i,pos2idx(j,0)]))) #kg H2O/kg dry air
            M_Evap[i,pos2idx(j,0)] = Evap_C[i,pos2idx(j,0)]*(X_surface[i,pos2idx(j,0)]-X_air[i, pos2idx(j,0)]) #kg H2O/m2-hr
        end
    end

    for i in 1:TL
        for j in 0:Ny_new
            Cum_Evap_Loss_vec[i,pos2idx(j,0)] = W*dy*M_Evap[i,pos2idx(j,0)] #kg H20/hr
        end
    end
    T_new = zeros(TL+1)
    for i = 1:TL
        T_new[i+1] = T[i]
    end

    Cum_Evap_Loss_T = zeros(TL)
    for i = 1:TL
        for j = 0:Ny_new
            Cum_Evap_Loss_T[i] = Cum_Evap_Loss_T[i] + Cum_Evap_Loss_vec[i,pos2idx(j,0)] #kg H20/hr
        end
        Cum_Evap_Loss[i+1] = Cum_Evap_Loss[i] + Cum_Evap_Loss_T[i]*(T_new[i+1]-T_new[i])
    end

    
    popfirst!(Cum_Evap_Loss)
 
    return Cum_Evap_Loss[TL]
end
