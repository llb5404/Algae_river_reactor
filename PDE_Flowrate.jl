
function PDE_Flowrate!(dX, A, Temperature, params, t)
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
 
     #Other parameters
     P_a = params.saturated_vapor_pressure_water_air(Tamb) 
     q(Temp) = params.lat_flow(Temp,WNDSPD,RH,P_a)

    alpha = params.alpha
    m = params.m

    for i=1:Nx
            dA[pos2idx(i)] = q(Temperature[pos2idx(i)]) -alpha*m*(((A[pos2idx(i)] + A[pos2idx(max(i-1,0))])/2)^(m-1))*(A[pos2idx(i)]-A[pos2idx(max(i-1,0))])/dx
    end

    @views dX[1+4*Nelements:5*Nelements] .= dA        # order matters!  The @views operator takes a slice out of an array without making a copy.
    nothing
end
