function InitModel()
    # read parameters
    dt              = p.dt
    L_a             = p.L_a
    L_b             = p.L_b
    coeffCurvature  = p.coeffCurvature
    N               = p.N
    v_ref           = p.v_ref
    c0              = p.c0

    # create model
    mdl     = Model(solver = IpoptSolver(print_level=3,max_cpu_time=0.1))

    @defVar( mdl, s[1:(N+1)] )
    @defVar( mdl, ey[1:(N+1)] )
    @defVar( mdl, epsi[1:(N+1)] )
    @defVar( mdl, v[1:(N+1)] )
    @defVar( mdl, 0.1 >= a[1:N] >= 0 )
    @defVar( mdl, -pi/6 <= d_f[1:N] <= pi/6 )

    # define objective function
    @setNLObjective(mdl, Min, sum{ey[i]^2+(v[i]-v_ref)^2,i=1:N}+5*(ey[N+1]^2+(v[N+1]-v_ref)^2))

    # define constraints
    @defNLParam(mdl, s0      == 0);     @addNLConstraint(mdl, s[1]      == s0);
    @defNLParam(mdl, ey0     == 0.1);   @addNLConstraint(mdl, ey[1]     == ey0);
    @defNLParam(mdl, epsi0   == 0);     @addNLConstraint(mdl, epsi[1]   == epsi0 );
    @defNLParam(mdl, v0      == 0);     @addNLConstraint(mdl, v[1]      == v0);

    @defNLParam(mdl, coeff[i=1:length(coeffCurvature)]==coeffCurvature[i]);

    @defNLExpr(mdl, c[i = 1:N],    coeff[1]*s[i]^4+coeff[2]*s[i]^3+coeff[3]*s[i]^2+coeff[4]*s[i]+coeff[5])
    @defNLExpr(mdl, bta[i = 1:N],  atan( L_a / (L_a + L_b) * tan( c0[3]*d_f[i] + c0[4]*abs(d_f[i])*d_f[i]) ) )
    @defNLExpr(mdl, dsdt[i = 1:N], v[i]*cos(epsi[i]+bta[i])/(1-ey[i]*c[i]))

    # define system dynamics
    for i in 1:N
        @addNLConstraint(mdl, s[i+1]     == s[i]       + dt*dsdt[i]  )
        @addNLConstraint(mdl, ey[i+1]    == ey[i]      + dt*v[i]*sin(epsi[i]+bta[i])  )
        @addNLConstraint(mdl, epsi[i+1]  == epsi[i]    + dt*(v[i]/L_a*sin(bta[i])-dsdt[i]*c[i])  )
        @addNLConstraint(mdl, v[i+1]     == v[i]       + dt*(c0[1]*a[i] - c0[2]*abs(v[i]) * v[i])  )
    end

    # return model
    return mdl
end

function changeObjective()
    v_ref = p.v_ref
    # define objective function
    @setNLObjective(mdl, Min, sum{ey[i]^2+(v[i]-v_ref)^2,i=1:N}+2*(ey[N+1]^2+(v[N+1]-v_ref)^2))
    return mdl
end


type mpcParams
    Hp
    nz
    OrderCostCons
end

type modelParams
    L_a
    L_b
    dt
    coeffCurvature
    N
    v_ref
    c0
end