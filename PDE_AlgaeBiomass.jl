

function PDE_AlgaeBiomass!(dX, C, DIC,CO2,Temperature,A, params, t)
    #Units for params can be found in LoadParameters.jl
    
    #Environmental Data
    GHI_Data = params.global_horizontal_irradiance_data
    Tamb_Data = params.ambient_temperature_data
    WNDSPD_Data = params.wind_speed_data
    RH_Data = params.relative_humidity_data

    #Timestep
    t_hour1 = floor(Int64, t)
    t_hour2 = floor(Int64, t)+1

    #Ensures that CO2 and DIC do not reach 0
    CO2 = max.(CO2, 1E-09)
    DIC = max.(DIC, 1E-09)
    
    data_begin = params.data_begin

    #Adjusts Environmental Data for each timestep
    GHI = GHI_Data[data_begin + t_hour1] * (t-t_hour1) + GHI_Data[data_begin + t_hour2] * (t_hour2 - t)
    Tamb = max.(((Tamb_Data[data_begin + t_hour1] * (t-t_hour1) + Tamb_Data[data_begin + t_hour2] * (t_hour2 - t))),0)
    RH = RH_Data[data_begin + t_hour1] * (t-t_hour1) + RH_Data[data_begin + t_hour2] * (t_hour2 - t)
    WNDSPD = max.((WNDSPD_Data[data_begin + t_hour1] * (t-t_hour1) + WNDSPD_Data[data_begin + t_hour2] * (t_hour2 - t)),0)

    #Geometric Properties
    L = params.reactor_length
    W = params.reactor_width
    Nx = params.num_odes_x 
    Ny = params.num_odes_y
    dx = L / Nx
 
    #Reindexing of 2D position in reactor to be contained in 1D vector
    pos2idx(y) = (y.+1)
    idx2pos(pos) = [Integer(pos - 1 - (Nx+1) * floor( (pos-1) ./ (Nx.+1))), Integer(floor( (pos-1) ./ (Nx.+1)))]

    #Number of total elements in vector
    Nelements = (Nx+1)
   
    C = max.(C, 1E-09) #biomass (g/L) don't allow negative values
    dC = zeros(Nelements, 1) #dC/dt
    dA = zeros(Nelements)
    dA .= dX[1+4*Nelements:5*Nelements]

    #Other parameters
    P_a = params.saturated_vapor_pressure_water_air(Tamb) 
    S(Area) = params.salinity(Area)
    bm = params.max_biomass_specific_growth_rate
    mu(T,S) = params.biomass_specific_growth_rate(T, S)
    alpha = params.alpha
    m = params.m

    #Initialize vectors for vertical increments, salinity, height
    H = zeros(Nelements)
    Sal = zeros(Nelements)
    Vavg = zeros(Nelements) #average velocity

    #Initialize Mass and Avg Velocity

    for i in 0:Nx

        Vavg[pos2idx(i)] = alpha*A[pos2idx(i)]^(m-1)
        H[pos2idx(i)] = A[pos2idx(i)]/(W)
        #salinity
        Sal[pos2idx(i)] = S(A[pos2idx(i)]) #kg/m3
        
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

    for i = 1:301
        for j = 0:Nx
            I_avg_vec[i,pos2idx(j)] = GHI*watt_to_umolm2s*percent_light[i]*exp(-absorbance[i]*(C[pos2idx(j)])*H[pos2idx(j)]/2) 
        end
    end

    #average light intensity summed across all wavelengths
    I_avg = zeros(Nelements,1)

    for i = 0:Nx   
            I_avg[pos2idx(i)] = sum(I_avg_vec[1:301,pos2idx(i)])
    end

    Y_xphm = 1.38/24 #maximal biomass yield on light, molx/mol
    a_x = 4.81 #absorbance cross section, m2/mol

    #Growth adjustment factor - light
    phiL = zeros(Nelements,1)

    
    for i = 0:Nx
            phiL[pos2idx(i)] = tanh((Y_xphm*a_x*I_avg[pos2idx(i)]*1E-06)/(bm*(1/3600)))
    end


    #Growth adjustment factor - CO2
    phiCO2 = zeros(Nelements, 1)

    #Piecewise linear function obtained from experiments on P celeri (see documentation)
    for i = 0:Nx
            if  CO2[pos2idx(i)] <= 0.1782/0.0315 && CO2[pos2idx(i)] >= 0.689655
                phiCO2[pos2idx(i)] = (0.0261*CO2[pos2idx(i)] - 0.018)/0.129383
            elseif CO2[pos2idx(i)] > 0.1782/0.0315 && CO2[pos2idx(i)] <= 29.67
                phiCO2[pos2idx(i)] = (-0.0054*CO2[pos2idx(i)] + 0.1602)/0.129383
            else
                phiCO2[pos2idx(i)] = 0
            end
    end


    #Biomass specific growth term
    mu_v = zeros(Nelements, 1)
    for i in 0:Nx
            mu_v[pos2idx(i)] = phiCO2[pos2idx(i)]*phiL[pos2idx(i)]*mu(Temperature[pos2idx(i)],Sal[pos2idx(i)]) -0.003621 
    end

    # ==========  ALGAE GROWTH ========

    #BC1: C(y,z) at y = 0 is constant [Cin]
    dC[pos2idx(0)] = 0.0

    q(Temp) = params.lat_flow(Temp,WNDSPD,RH,P_a)

    for i=0:Nx
        dA[pos2idx(i)] = q(Temperature[pos2idx(i)]) -alpha*m*(((A[pos2idx(i)] + A[pos2idx(max(i-1,0))])/2)^(m-1))*(A[pos2idx(i)]-A[pos2idx(max(i-1,0))])/dx
    end

    BM_growth = zeros(Nelements)
    BM_conv = zeros(Nelements)
    Term_adj = zeros(Nelements)


    for i=1:Nx
            BM_growth[pos2idx(i)] = mu_v[pos2idx(i)]* C[pos2idx(i)]
            BM_conv[pos2idx(i)] = -(A[pos2idx(i)]*Vavg[pos2idx(i)]*C[pos2idx(i)] - A[pos2idx(i-1)]*Vavg[pos2idx(i-1)]*C[pos2idx(i-1)])/(A[pos2idx(i)]*dx)
            Term_adj[pos2idx(i)] = -dA[pos2idx(i)]*C[pos2idx(i)]/A[pos2idx(i)]

            dC[pos2idx(i)] = (BM_conv[pos2idx(i)]
                                +BM_growth[pos2idx(i)]
                                +Term_adj[pos2idx(i)]
                               )
    end

    C = max.(C, 1E-09)
    
    @views dX[1:Nelements] .= dC        # order matters! The @views operator takes a slice out of an array without making a copy.

    nothing
end