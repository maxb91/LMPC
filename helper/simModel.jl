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
    zNext[4] = z[4] + dt*(u[1])#- 0.23 * z[4]^2 * sign(z[4]))   # v

    return zNext
end

function simModel_dyn_x(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::classes.ModelParams,exact_sim_i::Int64)

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

    max_alpha =modelParams.max_alpha
    m = modelParams.mass
    mu  = modelParams.mu
    g = modelParams.g
    I_z = modelParams.I_z
    B = modelParams.B
    C = modelParams.C

    F_xr = m*u[1]    
    FMax = mu*m*g / 2.0    
    if F_xr > FMax      
        F_xr = FMax    
    elseif F_xr < -FMax  
        F_xr = -FMax 
    end

    # determine slip angles      
    if v_x < 0.1        
        alpha_f = 0.0        
        alpha_r = 0.0      
    else        
        alpha_f = atan( (v_y+l_A*psi_dot) / v_x ) - u[2]        
        alpha_r = atan( (v_y-l_B*psi_dot) / v_x)       
    end
    if exact_sim_i==1 && max(abs(alpha_f),abs(alpha_r))>max_alpha/180*pi
        warn("Large slip angles: alpha_f = $(alpha_f*180/pi)°, alpha_r = $(alpha_r*180/pi)° , x =$x, y = $y")
    end
    
    F_yf = -FMax * sin(C*atan(B*alpha_f))
    F_yr = -FMax * sin(C*atan(B*alpha_r))
    # F_yf = -(mu*m*g / 2.0) * 1.02*alpha_f      
    # F_yr = -(mu*m*g / 2.0) * 1.02*alpha_r
    if F_yr > sqrt(FMax^2 - F_xr^2)        
        F_yr = sqrt(FMax^2 - F_xr^2)    
    elseif  F_yr < -sqrt(FMax^2 - F_xr^2)  
        F_yr = -sqrt(FMax^2 - F_xr^2) 
    end


    zNext[1] = x + dt*(v_x*cos(psi) - v_y*sin(psi))
    zNext[2] = y + dt*(v_x*sin(psi) + v_y*cos(psi))
    zNext[3] = v_x + dt*(psi_dot*v_y+1/m*(F_xr-F_yf*sin(u[2])))
    zNext[4] = v_y + dt*(-psi_dot*v_x+1/m*(F_yf*cos(u[2])+F_yr))
    zNext[5] = psi + dt*psi_dot
    zNext[6] = psi_dot + dt*(1/I_z*(l_A*F_yf*cos(u[2])-l_B*F_yr))

    return zNext
end

function simModel_exact_dyn_x(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::classes.ModelParams)
    dtn = dt/10
    t = 0:dtn:dt
    z_final = copy(z)
    for i=1:length(t)-1
        z_final = simModel_dyn_x(z_final,u,dtn,modelParams,i)
    end
    return z_final
end







