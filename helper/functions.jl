# oldTraj.oldTraj contains all state information of one lap. It is structured like follows:
# The first <prebuf> values are the end of one lap before the next lap starts (necessary for inter-lap-system-ID)
# The last value of <prebuf> ends with the last value *before* the finish line
# The s-values within <prebuf> need to be below zero (-> before the finish line). Otherwise they can be mistaken as the end of the current (not previous) lap.
# After <prebuf> follow the recorded states of the trajectory (as many as there were during one lap)
# After the recorded trajectory, the rest of the vector (until <buffersize>) is filled up with constant values
# oldCost is the number of values between the start and finish line (0 <= s < s_target)

function saveOldTraj(oldTraj::OldTrajectory,zCurr::Array{Float64},uCurr::Array{Float64},lapStatus::LapStatus,posInfo::PosInfo,buffersize::Int64,dt::Float64)
                
                i               = lapStatus.currentIt-1         # i = number of points for 0 <= s < s_target (= cost of this lap)
                prebuf          = oldTraj.prebuf                # so many points of the end of the previous old traj will be attached to the beginning
                zCurr_export    = zeros(buffersize,6)
                uCurr_export    = zeros(buffersize,2)

                zCurr_export    = cat(1,oldTraj.oldTraj[oldTraj.oldCost[1]+1:oldTraj.oldCost[1]+prebuf,:,1],
                                        zCurr[1:i,:], NaN*ones(buffersize-i-prebuf,6))#[ones(buffersize-i-prebuf,1)*zCurr[i,1:5] zCurr[i,6]+collect(1:buffersize-i-prebuf)*dt*zCurr[i,1]])
                uCurr_export    = cat(1,oldTraj.oldInput[oldTraj.oldCost[1]+1:oldTraj.oldCost[1]+prebuf,:,1],
                                        uCurr[1:i,:], NaN*ones(buffersize-i-prebuf,2))#zeros(buffersize-i-prebuf,2))

                zCurr_export[1:prebuf,6] -= posInfo.s_target       # make the prebuf-values below zero
                costLap                   = i                      # the cost of the current lap is the time it took to reach the finish line
                # println("zCurr_export:")
                # println(zCurr_export)
                # println("uCurr_export:")
                # println(uCurr_export)

                # Save all data in oldTrajectory:
                if lapStatus.currentLap <= 2                        # if it's the first or second lap
                    oldTraj.oldTraj[:,:,1]  = zCurr_export          # ... just save everything
                    oldTraj.oldInput[:,:,1] = uCurr_export
                    oldTraj.oldTraj[:,:,2]  = zCurr_export
                    oldTraj.oldInput[:,:,2] = uCurr_export
                    oldTraj.oldCost = [costLap,costLap]
                else                                                # idea: always copy the new trajectory in the first array!
                    if oldTraj.oldCost[1] < oldTraj.oldCost[2]      # if the first old traj is better than the second
                        oldTraj.oldTraj[:,:,2]  = oldTraj.oldTraj[:,:,1]    # ... copy the first in the second
                        oldTraj.oldInput[:,:,2] = oldTraj.oldInput[:,:,1]   # ... same for the input
                        oldTraj.oldCost[2] = oldTraj.oldCost[1]
                    end
                    oldTraj.oldTraj[:,:,1]  = zCurr_export                 # ... and write the new traj in the first
                    oldTraj.oldInput[:,:,1] = uCurr_export
                    oldTraj.oldCost[1] = costLap
                end
end

function InitializeModel(m::MpcModel,mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff,z_Init::Array{Float64,1})

    dt   = modelParams.dt
    L_a  = modelParams.l_A
    L_b  = modelParams.l_B
    c0   = modelParams.c0
    u_lb = modelParams.u_lb
    u_ub = modelParams.u_ub
    z_lb = modelParams.z_lb
    z_ub = modelParams.z_ub

    N    = mpcParams.N

    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
    
    m.mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=0.5))#,check_derivatives_for_naninf="yes"))#,linear_solver="ma57",print_user_options="yes"))

    @variable( m.mdl, m.z_Ol[1:(N+1),1:6])      # z = s, ey, epsi, v
    @variable( m.mdl, m.u_Ol[1:N,1:2])
    @variable( m.mdl, 0 <= m.ParInt[1:1] <= 1)

    z_lb_6s = ones(mpcParams.N+1,1)*[0.1 -Inf -Inf -Inf -Inf -Inf]                    # lower bounds on states
    z_ub_6s = ones(mpcParams.N+1,1)*[3.0  Inf  Inf  Inf  Inf  Inf]                    # upper bounds
    u_lb_6s = ones(mpcParams.N,1) * [-1.0  -pi/6]                    # lower bounds on steering
    u_ub_6s = ones(mpcParams.N,1) * [3.0    pi/6]                  # upper bounds

    for i=1:2       # I don't know why but somehow the short method returns errors sometimes
        for j=1:N
            setlowerbound(m.u_Ol[j,i], u_lb_6s[j,i])
            setupperbound(m.u_Ol[j,i], u_ub_6s[j,i])
        end
    end
    for i=1:6       # I don't know why but somehow the short method returns errors sometimes
        for j=1:N
            setlowerbound(m.z_Ol[j,i], z_lb_6s[j,i])
            setupperbound(m.z_Ol[j,i], z_ub_6s[j,i])
        end
    end
    #for i=1:N+1
    #    setlowerbound(m.z_Ol[i,1],0.1)
    #end

    @NLparameter(m.mdl, m.z0[i=1:6] == z_Init[i])
    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i]);
    @NLparameter(m.mdl, m.c_Vx[i=1:3]  == 1)
    @NLparameter(m.mdl, m.c_Vy[i=1:4]  == 1)
    @NLparameter(m.mdl, m.c_Psi[i=1:3] == 1)
    #setvalue(m.c_Vx[3],1)
    #setvalue(m.c_Psi[3],1)

    @NLconstraint(m.mdl, [i=1:6], m.z_Ol[1,i] == m.z0[i])

    @NLexpression(m.mdl, m.c[i = 1:N], sum{m.coeff[j]*m.z_Ol[i,6]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1])
    #@NLexpression(m.mdl, m.c[i = 1:N], 0)

    @NLexpression(m.mdl, m.dsdt[i = 1:N], (m.z_Ol[i,1]*cos(m.z_Ol[i,4]) - m.z_Ol[i,2]*sin(m.z_Ol[i,4]))/(1-m.z_Ol[i,5]*m.c[i]))
    #@NLexpression(m.mdl, m.dsdt[i = 1:N], 1)
    
    println("Initializing model...")
    #println("z0   = $(getvalue(m.z0))")
    #println("z_Ol = $(getvalue(m.z_Ol))")
    #println("m.coeff = $(getvalue(m.coeff))")
    #println("m.c     = $(getvalue(m.c))")
    #println("n_poly_curv = $n_poly_curv")
    #println("dt = $dt")
    # System dynamics
    for i=1:N
        @NLconstraint(m.mdl, m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*(m.z_Ol[i,1] + m.c_Vx[1]*m.z_Ol[i,2] + m.c_Vx[2]*m.z_Ol[i,3] + m.c_Vx[3]*m.u_Ol[i,1]))  # xDot
        @NLconstraint(m.mdl, m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*(m.z_Ol[i,2] + m.c_Vy[1]*m.z_Ol[i,2]/m.z_Ol[i,1] + m.c_Vy[2]*m.z_Ol[i,1]*m.z_Ol[i,3] +
                                                               m.c_Vy[3]*m.z_Ol[i,3]/m.z_Ol[i,1] + m.c_Vy[4]*m.u_Ol[i,2]))
        @NLconstraint(m.mdl, m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,3] + m.c_Psi[1]*m.z_Ol[i,3]/m.z_Ol[i,1] + m.c_Psi[2]*m.z_Ol[i,2]/m.z_Ol[i,1] +
                                                               m.c_Psi[3]*m.u_Ol[i,2]))                       # psiDot
        @NLconstraint(m.mdl, m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.z_Ol[i,3]-m.dsdt[i]*m.c[i]))                                                   # ePsi
        @NLconstraint(m.mdl, m.z_Ol[i+1,5]  == m.z_Ol[i,5] + dt*(m.z_Ol[i,1]*sin(m.z_Ol[i,4])+m.z_Ol[i,2]*cos(m.z_Ol[i,4])))                    # eY
        @NLconstraint(m.mdl, m.z_Ol[i+1,6]  == m.z_Ol[i,6] + dt*m.dsdt[i]  )

        # @NLconstraint(m.mdl, m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )                                                    # s
        # @NLconstraint(m.mdl, m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                        # ey
        # @NLconstraint(m.mdl, m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_a*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )             # epsi
        # @NLconstraint(m.mdl, m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1] - 0.63*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))            # v
        # @NLconstraint(m.mdl, m.z_Ol[i+1,5]  == m.z_Ol[i,5] + dt*m.dsdt[i]  )                                                    # s
        # @NLconstraint(m.mdl, m.z_Ol[i+1,6]  == m.z_Ol[i,6] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                        # ey
    end

    @NLexpression(m.mdl, costZ, 0.5*sum{(m.z_Ol[j,5])^2+(m.z_Ol[j,1]-0.2)^2,j=2:N+1})
    # solve model once
    println("solving...")
    solve(m.mdl)
    #solve(m.mdl)
    #solve(m.mdl)
    #solve(m.mdl)
    #println("Solution: $(getvalue(m.z_Ol))")
    println("finished")

end


function InitializeModel_pathFollow(m::MpcModel_pF,mpcParams::MpcParams,modelParams::ModelParams,trackCoeff::TrackCoeff,z_Init::Array{Float64,1})
    dt   = modelParams.dt
    L_a  = modelParams.l_A
    L_b  = modelParams.l_B
    c0   = modelParams.c0
    u_lb = modelParams.u_lb
    u_ub = modelParams.u_ub
    z_lb = modelParams.z_lb
    z_ub = modelParams.z_ub

    N    = mpcParams.N

    n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation
    
    m.mdl = Model(solver = IpoptSolver(print_level=0,max_cpu_time=0.1))#,linear_solver="ma57",print_user_options="yes"))

    @variable( m.mdl, m.z_Ol[1:(N+1),1:4])      # z = s, ey, epsi, v
    @variable( m.mdl, m.u_Ol[1:N,1:2])

    for i=1:2       # I don't know why but somehow the short method returns errors sometimes
        for j=1:N
            setlowerbound(m.u_Ol[j,i], modelParams.u_lb[j,i])
            setupperbound(m.u_Ol[j,i], modelParams.u_ub[j,i])
        end
    end
    for i=1:4       # I don't know why but somehow the short method returns errors sometimes
        for j=1:N
            setlowerbound(m.z_Ol[j,i], modelParams.z_lb[j,i])
            setupperbound(m.z_Ol[j,i], modelParams.z_ub[j,i])
        end
    end
    #@variable( m.mdl, 1 >= m.ParInt >= 0 )

    @NLparameter(m.mdl, m.z0[i=1:4] == z_Init[i])
    @NLconstraint(m.mdl, [i=1:4], m.z_Ol[1,i] == m.z0[i])

    @NLparameter(m.mdl, m.coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i]);

    #@NLexpression(m.mdl, m.c[i = 1:N],    m.coeff[1]*m.z_Ol[1,i]^4+m.coeff[2]*m.z_Ol[1,i]^3+m.coeff[3]*m.z_Ol[1,i]^2+m.coeff[4]*m.z_Ol[1,i]+m.coeff[5])
    @NLexpression(m.mdl, m.c[i = 1:N],    sum{m.coeff[j]*m.z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + m.coeff[n_poly_curv+1])
    @NLexpression(m.mdl, m.bta[i = 1:N],  atan( L_a / (L_a + L_b) * ( m.u_Ol[i,2] ) ) )
    @NLexpression(m.mdl, m.dsdt[i = 1:N], m.z_Ol[i,4]*cos(m.z_Ol[i,3]+m.bta[i])/(1-m.z_Ol[i,2]*m.c[i]))

    # System dynamics
    for i=1:N
        @NLconstraint(m.mdl, m.z_Ol[i+1,1]  == m.z_Ol[i,1] + dt*m.dsdt[i]  )                                             # s
        @NLconstraint(m.mdl, m.z_Ol[i+1,2]  == m.z_Ol[i,2] + dt*m.z_Ol[i,4]*sin(m.z_Ol[i,3]+m.bta[i])  )                     # ey
        @NLconstraint(m.mdl, m.z_Ol[i+1,3]  == m.z_Ol[i,3] + dt*(m.z_Ol[i,4]/L_a*sin(m.bta[i])-m.dsdt[i]*m.c[i])  )            # epsi
        @NLconstraint(m.mdl, m.z_Ol[i+1,4]  == m.z_Ol[i,4] + dt*(m.u_Ol[i,1] - modelParams.c_f*abs(m.z_Ol[i,4]) * m.z_Ol[i,4]))  # v
    end
    solve(m.mdl)
end

function InitializeParameters(mpcParams::MpcParams,mpcParams_pF::MpcParams,trackCoeff::TrackCoeff,modelParams::ModelParams,
                                posInfo::PosInfo,oldTraj::OldTrajectory,mpcCoeff::MpcCoeff,lapStatus::LapStatus,buffersize::Int64)
    mpcParams.N                 = 5
    mpcParams.Q                 = [0.0,10.0,0.0,1.0]      # put weights on ey, epsi and v
    mpcParams.Q_term            = 1.0*[1.0,1.0,1.0,1.0,1.0]        # weights for terminal constraints (LMPC, for e_y, e_psi, and v)
    mpcParams.R                 = 0*[1.0,1.0]             # put weights on a and d_f
    mpcParams.QderivZ           = 0.0*[0,0,0.1,0,0,0]       # cost matrix for derivative cost of states
    mpcParams.QderivU           = 0.1*[1,1]               # cost matrix for derivative cost of inputs

    mpcParams_pF.N              = 5
    mpcParams_pF.Q              = [0.0,10.0,0.0,1.0]
    mpcParams_pF.R              = 0*[1.0,1.0]             # put weights on a and d_f
    mpcParams_pF.QderivZ        = 0.0*[0,0,0.1,0]       # cost matrix for derivative cost of states
    mpcParams_pF.QderivU        = 0.1*[1,1]               # cost matrix for derivative cost of inputs
    mpcParams_pF.vPathFollowing = 0.3                     # reference speed for first lap of path following

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.coeffCurvature   = zeros(trackCoeff.nPolyCurvature+1)         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 0.6                     # width of the track (0.5m)

    modelParams.u_lb            = ones(mpcParams.N,1) * [-0.2  -pi/6]                    # lower bounds on steering
    modelParams.u_ub            = ones(mpcParams.N,1) * [1.2   pi/6]                  # upper bounds
    modelParams.z_lb            = ones(mpcParams.N+1,1)*[-Inf -Inf -Inf -0.1]                    # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1)*[ Inf  Inf  Inf  2.0]                    # upper bounds
    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125
    modelParams.dt              = 0.1
    modelParams.m               = 1.98
    modelParams.I_z             = 0.24
    modelParams.c_f             = 0.63                 # friction coefficient: xDot = - c_f*xDot² (aerodynamic+tire)

    posInfo.s_start             = 0.0
    posInfo.s_target            = 5.0

    oldTraj.oldTraj             = NaN*ones(buffersize,6,2)
    #oldTraj.oldTraj[:,6,1]      = 1:buffersize
    #oldTraj.oldTraj[:,6,2]      = 1:buffersize
    oldTraj.oldInput            = NaN*ones(buffersize,2,2)
    oldTraj.oldCost             = ones(Int64,2)                   # dummies for initialization
    oldTraj.prebuf              = 30
    oldTraj.postbuf             = 30

    mpcCoeff.order              = 5
    mpcCoeff.coeffCost          = zeros(mpcCoeff.order+1,2)
    mpcCoeff.coeffConst         = zeros(mpcCoeff.order+1,2,5)
    mpcCoeff.pLength            = 4*mpcParams.N        # small values here may lead to numerical problems since the functions are only approximated in a short horizon
    mpcCoeff.c_Vx               = zeros(3)
    mpcCoeff.c_Vy               = zeros(4)
    mpcCoeff.c_Psi              = zeros(3)

    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap
end

function extendOldTraj(oldTraj::OldTrajectory,posInfo::PosInfo,zCurr::Array{Float64,2},uCurr::Array{Float64,2})
    #println(size(zCurr[1:21,:]))
    #println(size(oldTraj.oldTraj[oldTraj.oldCost[1]:oldTraj.oldCost[1]+20,1:4,1]))
    postbuf = oldTraj.postbuf
    for i=1:6
        for j=1:postbuf
            oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+j,i,1]  = zCurr[j,i]
        end
    end
    for i=1:2
        for j=1:postbuf
            oldTraj.oldInput[oldTraj.oldCost[1]+oldTraj.prebuf+j,i,1] = uCurr[j,i]
        end
    end
    oldTraj.oldTraj[oldTraj.oldCost[1]+oldTraj.prebuf+1:oldTraj.oldCost[1]+oldTraj.prebuf+postbuf,6,1] += posInfo.s_target
    println("oldTraj extended:")
    println(oldTraj.oldTraj[:,:,1])
end