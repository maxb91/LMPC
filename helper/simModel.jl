function simModel_s(z::Array{Float64},u::Array{Float64},dt::Float64,coeff::Array{Float64},modelParams::classes.ModelParams)

    local zNext::Array{Float64}
    L_a = modelParams.l_A
    L_b = modelParams.l_B
    c = ([z[1]^4 z[1]^3 z[1]^2 z[1] 1]*coeff)[1]        # Polynomial to approximate curvature kappa

    bta = atan(L_b/(L_a+L_b)*tan(u[2])) #?? warum to abs ^2 is nich in formel
    dsdt = z[4]*cos(z[3]+bta)/(1-z[2]*c)

    zNext = z
    zNext[1] = z[1] + dt*dsdt
    zNext[2] = z[2] + dt*z[4] * sin(z[3] + bta)
    zNext[3] = z[3] + dt*(z[4]/L_b*sin(bta)-dsdt*c)
    zNext[4] = z[4] + dt*(u[1] - 0.63*abs(z[4])*z[4])#sign)=

    return zNext
end


function simModel_x(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::classes.ModelParams)

   # kinematic bicycle model
   # u[1] = acceleration
   # u[2] = steering angle

    local zNext::Array{Float64}
    l_A = modelParams.l_A
    l_B = modelParams.l_B

    bta = atan(l_B/(l_A+l_B)*tan(u[2]))

    zNext = z
    zNext[1] = z[1] + dt*(z[4]*cos(z[3]+bta))       # x
    zNext[2] = z[2] + dt*(z[4]*sin(z[3] + bta))      # y
    zNext[3] = z[3] + dt*(z[4]/l_B*sin(bta))        # psi
    zNext[4] = z[4] + dt*(u[1])#- 0.23 * z[4]^2 * sign(z[4]))                #0.63     # v

    return zNext
end

function simModel_dyn_x(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::classes.ModelParams)

    local zNext::Array{Float64}
    l_A = modelParams.l_A
    l_B = modelParams.l_B
    zNext = z
    x = z[1]
    y = z[2]
    v_x = z[3]
    v_y = z[4]
    psi = z[5]
    psi_dot = z[6]
    # u[1] = acceleration
    # u[2] = steering angle

    m = 1.98 # kg
    mu  = 0.8
    g = 9.81 # m/s^2
    I_z = 0.03 # kg * m^2
    B = 1.0
    C = 1.25




    F_xr = m*u[1]    
    FxMax = mu*m*g / 2.0    
    if F_xr > FxMax      
        F_xr = FxMax    
    end

    # determine slip angles      
    if v_x < 0.1        
        alpha_f = 0.0        
        alpha_r = 0.0      
    else        
        alpha_f = atan( (v_y+l_A*psi_dot) / v_x ) - u[2]        
        alpha_r = atan( (v_y-l_B*psi_dot) / v_x)      
    end
    
    F_yf = -mu*m*g/2.0 * sin(C*atan(B*alpha_f))
    F_yr = -mu*m*g/2.0 * sin(C*atan(B*alpha_r))
    # F_yf = -(mu*m*g / 2.0) * 1.02*alpha_f      
    # F_yr = -(mu*m*g / 2.0) * 1.02*alpha_r
    if F_yr > sqrt(FxMax^2 - F_xr^2)        
        F_yr = sqrt(FxMax^2 - F_xr^2)      
    end


    zNext[1] = x + dt*(v_x*cos(psi) - v_y*sin(psi))
    zNext[2] = y + dt*(v_x*sin(psi) + v_y*cos(psi))
    zNext[3] = v_x + dt*(psi_dot*v_y+1/m*(F_xr-F_yf*sin(u[2])))
    zNext[4] = v_y + dt*(-psi_dot*v_x+1/m*(F_yf*cos(u[2])+F_yr))
    zNext[5] = psi + dt*psi_dot
    zNext[6] = psi_dot + dt*(1/I_z*(l_A*F_yf*cos(u[2])-l_B*F_yr))

    return zNext
end









