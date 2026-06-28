module parameters

function waves()
    @show L = 1e-2
    @show T0 = 1.9
    @show rho0 = 145.4684
    @show s0 = 725.5
    @show xn0 = 61.03/rho0
    @show xs0 = 1.0 - xn0 
    @show u1 = 40.0
    @show u2 = 18.83
    @show t_end = 2*sqrt(2)*L/u1
    @show dt_frame = t_end/100
    @show C = T0*s0*s0*xs0/(xn0*u2*u2)
    @show x_prime = 1.17*T0/C
    #@show dr = L/N
    #@show m = dr*dr*rho0
    #@show h = 6.0*dr
    #@show dt = dt_hat*h/u1
end

function fountain()
    @show N = 20
    @show r_cap = 1.55e-3/2
    @show dr = r_cap/N       # average particle distance
    @show h = 2.8*dr        # size of kernel support
    @show rho0 = 145.2352   # fluid density
    @show m = rho0*dr^2     # particle mass
    @show c1 = 40.0         # numerical speed of first sound
    @show c2 = 20.37        # [m/s]
    #@show g = -9.8*VECY         #gravitational acceleration
    @show mu = 1.295e-4       #dynamic viscosity of helium

    @show Xn = 28.09/rho0
    @show Xs = 1.0 - Xn

    @show T = 1.65          # [K]
    @show S = 335.0        # [J/(kg*K)]
    @show X_prime = 0.66*(Xn/Xs)*(c2/S)^2
    @show beta = 1e-2*(r_cap*rho0*S*S/c2)

    ##temporal parameters
    @show dt = 0.2*h/c1
    @show t_end = 0.3
    @show dt_frame = max(dt, t_end/200)

    @show CELL_W = 25e-3
    @show CELL_Z1 =  120e-3
    @show CELL_Z0 = -10e-3
    @show HE_LEVEL = 23e-3
    @show wwall = 3.0*dr

    ##superleak parameters
    @show superleak_w = 15e-3
    @show superleak_y = 1e-3
    @show superleak_d = 2e-3
    @show superleak_k = 1e6 # [1/s]

    ##heater parameters
    @show heater_r = 4.0e-3
    @show heater_y = 8.0e-3
    @show heater_t = 0.02
    @show w_dot = 100.0 # [W/m]

    @show C = T*S*S*Xs/(Xn*c2*c2)
end

end
