function simModel(z,u,dt,coeff,modelParams)

    L_a = modelParams.l_A
    L_b = modelParams.l_B
    c0  = modelParams.c0

    c = ([z[1]^4 z[1]^3 z[1]^2 z[1] 1]*coeff)[1]        # Polynomial

    bta = atan(L_a/(L_a+L_b)*tan(c0[3]*u[2]+c0[4]*abs(u[2])*u[2]))
    dsdt = z[4]*cos(z[3]+bta)/(1-z[2]*c)

    zNext = z
    zNext[1] = z[1] + dt*dsdt
    zNext[2] = z[2] + dt*z[4] * sin(z[3] + bta)
    zNext[3] = z[3] + dt*(z[4]/L_a*sin(bta)-dsdt*c)
    zNext[4] = z[4] + dt*(c0[1]*u[1] - c0[2]*abs(z[4])*z[4])

    zNext = zNext + 0*randn(1,4)*0.001
    return zNext
end