function simModel(z::Array{Float64},u::Array{Float64},dt::Float64,coeff::Array{Float64},modelParams::ModelParams)

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
    zNext[4] = z[4] + dt*(u[1] - 0.63*abs(z[4])*z[4])

    return zNext
end