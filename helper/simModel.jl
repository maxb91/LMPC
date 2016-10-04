# z[1] = xDot
# z[2] = yDot
# z[3] = psiDot
# z[4] = ePsi
# z[5] = eY
# z[6] = s

function simKinModel(z::Array{Float64},u::Array{Float64},dt::Float64,coeff::Array{Float64},modelParams::ModelParams)

    zNext::Array{Float64}
    L_a = modelParams.l_A
    L_b = modelParams.l_B
    c = ([z[1]^4 z[1]^3 z[1]^2 z[1] 1]*coeff)[1]        # Polynomial

    bta = atan(L_a/(L_a+L_b)*tan(u[2]+abs(u[2])*u[2]))
    dsdt = z[4]*cos(z[3]+bta)/(1-z[2]*c)

    zNext = z
    zNext[1] = z[1] + dt*dsdt
    zNext[2] = z[2] + dt*z[4] * sin(z[3] + bta)
    zNext[3] = z[3] + dt*(z[4]/L_a*sin(bta)-dsdt*c)
    zNext[4] = z[4] + dt*(u[1] - modelParams.c_f*abs(z[4])*z[4])

    return zNext
end

function simDynModel_exact(z::Array{Float64},u::Array{Float64},dt::Float64,coeff::Array{Float64},modelParams::ModelParams)
    # This function uses smaller steps to achieve higher fidelity than we would achieve using longer timesteps
    z_final = copy(z)
    u[1] = min(u[1],3)
    u[1] = max(u[1],-3)
    u[2] = min(u[2],pi/6)
    u[2] = max(u[2],-pi/6)
    dtn = dt/100
    println("0. u = $u")
    println("1. z = $(z_final)")
    for i=1:100
        z_final = simDynModel(z_final,u,dtn,coeff,modelParams)
    end
    println("2. z = $(z_final)")
    return z_final
end

function simDynModel(z::Array{Float64},u::Array{Float64},dt::Float64,coeff::Array{Float64},modelParams::ModelParams)

    zNext::Array{Float64}
    L_f = modelParams.l_A
    L_r = modelParams.l_B
    c0  = modelParams.c0
    m   = modelParams.m
    I_z = modelParams.I_z
    c_f = modelParams.c_f

    #z[7] = u[1]
    #z[8] = u[2]

    a_F = 0
    a_R = 0
    if abs(z[1]) >= 0.02
        a_F     = atan((z[2] + L_f*z[3])/z[1]) - z[8]
        a_R     = atan((z[2] - L_r*z[3])/z[1])
    end
    if max(abs(a_F),abs(a_R))>20/180*pi
        warn("Large tire angles: a_F = $a_F, a_R = $a_R, xDot = $(z[1]), d_F = $(z[8])")
    end
    
    FyF = -pacejka(a_F)
    FyR = -pacejka(a_R)
    
    c = ([z[1]^4 z[1]^3 z[1]^2 z[1] 1]*coeff)[1]                        # Polynomial for curvature
    
    dsdt = (z[1]*cos(z[4]) - z[2]*sin(z[4]))/(1-z[5]*c)

    zNext = copy(z)
    zNext[1] = z[1] + dt * (z[7] + z[2]*z[3] - c_f*z[1]^2*sign(z[1]))       # xDot
    zNext[2] = z[2] + dt * (2/m*(FyF*cos(z[8]) + FyR) - z[3]*z[1])          # yDot
    zNext[3] = z[3] + dt * (2/I_z*(L_f*FyF - L_r*FyR))                      # psiDot
    zNext[4] = z[4] + dt * (z[3]-dsdt*c)                                    # ePsi
    zNext[5] = z[5] + dt * (z[1]*sin(z[4]) + z[2]*cos(z[4]))                # eY
    zNext[6] = z[6] + dt * dsdt                                             # s
    zNext[7] = z[7] + dt * (u[1] - z[7]) * 0.1/dt                          # a
    zNext[8] = z[8] + dt * (u[2] - z[8]) * 0.1/dt                         # d_f
    #zNext[8] = z[8] + dt * sign(u[2] - z[8]) * 0.01/dt                     # d_f

    return zNext
end

function pacejka(a)
    B = 0.3#20
    C = 1.25
    mu = 0.234
    m = 1.98
    g = 9.81
    D = mu * m * g/2
    D = D*100.0

    C_alpha_f = D*sin(C*atan(B*a))
    return C_alpha_f
end
