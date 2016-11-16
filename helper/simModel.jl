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