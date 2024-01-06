##test height
using Roots

n = 0.05 #wood
Q = (85*(3.28)^3)/3600 #vol flow rate, ft3/s
S = tan(0.0174) #slope angle (1 deg)
W = 0.5*3.28 #ft

f(x) = W*x*((W*x)/(2*x + W))^(2/3) - (Q*n)/(1.49*S^0.5)
@show height = (find_zero(f, (0,10), Bisection()))/3.28

@show V = Q/(height*3.28*W)
density_water = 1000*(1/3.38)^3 #kg/ft3
viscosity = 0.001002*(1/3.28) #kg/ft/s

@show reynolds_number = V .* density_water * ((W*height*3.28)/(2*height*3.28 + W)) ./ (viscosity)


