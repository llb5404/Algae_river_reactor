
function PDE_Flowrate!(dX, A, C_biomass,Temperature, params, t)
    #Units for params can be found in LoadParameters.jl

    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data
    data_begin = params.data_begin

    #Timestep
    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    #Adjusts Environmental Data for each timestep
    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)


     #Geometric Properties
     L = params.reactor_length
     Nx = params.num_odes_x
     dx = L / Nx
 
     #Reindexing of 2D position in reactor to be contained in 1D vector
     pos2idx(y) = (y.+1)

     #Number of total elements in vector
     Nelements = (Nx+1) #number of elements in vectors containing info in y-direction
     A = max.(A, 0.0) #don't allow negative values
     dA = zeros(Nelements, 1) #dT/dt
     H = zeros(Nelements,1)
     W = params.reactor_width

    for i in 0:Nx
        H[pos2idx(i)] = A[pos2idx(i)]/(W)
    end


      #Calculation of light adjustment factor for biomass growth
    wavelength = zeros(301,1)
    absorbance = zeros(301,1)
    percent_light = zeros(301,1)
    
    csv_reader1 = CSV.File("WAVE_ABS.csv", header=["col1", "col2","col3"])
    for (i,row) in enumerate(csv_reader1)
        wavelength[i] = convert(Float64, row.col1) #wavelength of light, nm 400-700
        absorbance[i] = convert(Float64, row.col2) #absorbance of wavelength of light by P celeri, m2/mol
        percent_light[i] = convert(Float64,row.col3) #percentage of sunlight in each wavelength
    end

    watt_to_umolm2s = 0.425*4.6 #watt to umolm2s conversion

    #average light intensity for each wavelength
    I_avg_vec = zeros(301,Nelements)
    I_avg = zeros(Nelements)

    for i = 1:301
        for j = 0:Nx
            I_avg_vec[i,pos2idx(j)] = GHI*watt_to_umolm2s*percent_light[i]*exp(-absorbance[i]*(C_biomass[pos2idx(j)])*H[pos2idx(j)]/2) 
        end
    end

    for i = 0:Nx
            I_avg[pos2idx(i)] = sum(I_avg_vec[1:301,pos2idx(i)])
    end

 
     #Other parameters
     P_a = params.saturated_vapor_pressure_water_air(Tamb) 
     q(Temp,phiL) = params.lat_flow(Temp,WNDSPD,RH,P_a,phiL)

    alpha = params.alpha
    m = params.m

    for i=1:Nx
            dA[pos2idx(i)] = q(Temperature[pos2idx(i)],I_avg[pos2idx(i)]) -alpha*m*(((A[pos2idx(i)] + A[pos2idx(max(i-1,0))])/2)^(m-1))*(A[pos2idx(i)]-A[pos2idx(max(i-1,0))])/dx
    end

    @views dX[1+4*Nelements:5*Nelements] .= dA        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end
