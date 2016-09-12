#!/usr/bin/env julia

#=
 Licensing Information: You are free to use or extend these projects for 
 education or reserach purposes provided that (1) you retain this notice
 and (2) you provide clear attribution to UC Berkeley, including a link 
 to http://barc-project.com

 Attibution Information: The barc project ROS code-base was developed
 at UC Berkeley in the Model Predictive Control (MPC) lab by Jon Gonzales
 (jon.gonzales@berkeley.edu). The cloud services integation with ROS was developed
 by Kiet Lam  (kiet.lam@berkeley.edu). The web-server app Dator was 
 based on an open source project by Bruce Wootton
=# 

using RobotOS
@rosimport barc.msg: ECU, pos_info, Encoder, Ultrasound, Z_KinBkMdl, Logging
@rosimport data_service.msg: TimeData
@rosimport geometry_msgs.msg: Vector3
rostypegen()
using barc.msg
using data_service.msg
using geometry_msgs.msg
using JuMP
using Ipopt


# define model parameters
L_a     = 0.125         # distance from CoG to front axel
L_b     = 0.125         # distance from CoG to rear axel
dt      = 0.1           # time step of system
coeffCurvature   = [0,0,0,0,0,0,0]

# preview horizon
N       = 10

# define targets [generic values]
v_ref   = 0.2

c0 = [0.5431, 1.2767, 2.1516, -2.4169]

# c_1 = 0.5431, c_2 = 1.2767, c_3 = 2.1516, c_4 = -2.4169

# define decision defVars 
# states: position (x,y), yaw angle, and velocity
# inputs: acceleration, steering angle 
println("Creating kinematic bicycle model ....")
mdl     = Model(solver = IpoptSolver(print_level=3,max_cpu_time=0.1))

@defVar( mdl, s[1:(N+1)] )
@defVar( mdl, ey[1:(N+1)] )
@defVar( mdl, epsi[1:(N+1)] )
@defVar( mdl, v[1:(N+1)] )
@defVar( mdl, 0.1 >= a[1:N] >= 0 )
@defVar( mdl, -pi/6 <= d_f[1:N] <= pi/6 )

# define objective function
# @setNLObjective(mdl, Min, (x[N+1] - x_ref)^2 + (y[N+1] - y_ref)^2 )
@setNLObjective(mdl, Min, sum{ey[i]^2+(v[i]-v_ref)^2,i=1:N}+5*(ey[N+1]^2+(v[N+1]-v_ref)^2))

# define constraints
# define system dynamics
# Reference: R.Rajamani, Vehicle Dynamics and Control, set. Mechanical Engineering Series,
#               Spring, 2011, page 26
@defNLParam(mdl, s0      == 0); @addNLConstraint(mdl, s[1]      == s0);
@defNLParam(mdl, ey0     == 0.1); @addNLConstraint(mdl, ey[1]     == ey0);
@defNLParam(mdl, epsi0   == 0); @addNLConstraint(mdl, epsi[1]   == epsi0 );
@defNLParam(mdl, v0      == 0); @addNLConstraint(mdl, v[1]      == v0);
@defNLParam(mdl, coeff[i=1:length(coeffCurvature)]==coeffCurvature[i]);
@defNLExpr(mdl, c[i = 1:N],    coeff[1]*s[i]^6+coeff[2]*s[i]^5+coeff[3]*s[i]^4+coeff[4]*s[i]^3+coeff[5]*s[i]^2+coeff[6]*s[i]+coeff[7])
# @defNLExpr(mdl, c[i = 1:N],    sum{coeff[j]*s[i]^(7-j),j=1:7})
# @defNLExpr(mdl, bta[i = 1:N],    atan( L_a / (L_a + L_b) * tan( 1.3*d_f[i]) ) ) # workaround since there's no sign function
@defNLExpr(mdl, bta[i = 1:N],    atan( L_a / (L_a + L_b) * tan( c0[3]*d_f[i] + c0[4]*abs(d_f[i])*d_f[i]) ) )
# @defNLExpr(mdl, bta[i = 1:N],    atan( L_a / (L_a + L_b) * tan( 2.1516*d_f[i]) ) )

@defNLExpr(mdl, dsdt[i = 1:N], v[i]*cos(epsi[i]+bta[i])/(1-ey[i]*c[i]))
for i in 1:N
    @addNLConstraint(mdl, s[i+1]     == s[i]       + dt*dsdt[i]  )
    @addNLConstraint(mdl, ey[i+1]    == ey[i]      + dt*v[i]*sin(epsi[i]+bta[i])  )
    @addNLConstraint(mdl, epsi[i+1]  == epsi[i]    + dt*(v[i]/L_a*sin(bta[i])-dsdt[i]*c[i])  )
    @addNLConstraint(mdl, v[i+1]     == v[i]       + dt*(c0[1]*a[i] - c0[2]*abs(v[i]) * v[i])  )
    # @addNLConstraint(mdl, v[i+1]     == v[i]       + dt*(0.4*a[i] - 0.3 * v[i])  )
    # @addNLConstraint(mdl, v[i+1]     == v[i]       + dt*(0.5431*a[i])  )
end


# status update
println("initial solve ...")
solve(mdl)
println("finished initial solve!")

function SE_callback(msg::pos_info)
    # update mpc initial condition 
    setValue(s0,     msg.s)
    setValue(ey0,    msg.ey)
    setValue(epsi0,  msg.epsi)
    setValue(v0,     msg.v)
    setValue(coeff,  msg.coeffCurvature)
end

function main()
    # initiate node, set up publisher / subscriber topics
    init_node("mpc_traj")
    pub = Publisher("ecu", ECU, queue_size=10)
    pub2 = Publisher("logging", Logging, queue_size=10)
    s1  = Subscriber("pos_info", pos_info, SE_callback, queue_size=10)
    loop_rate = Rate(10)
    cmdcount = 0
    failcnt = 0
    while ! is_shutdown()
        # run mpc, publish command
        status = solve(mdl)
        if status == Symbol("Optimal")
            cmd = ECU(a_opt, d_f_opt)
        if status != Symbol("Optimal")
            failcnt = failcnt + 1
            if failcnt == 10
                cmd = ECU(0,0)
        else
            failcnt = 0
        end
        # get optimal solutions
        a_opt   = getValue(a[1])
        d_f_opt = getValue(d_f[1])
        # TO DO: transform to PWM signals
        cmd = ECU(a_opt, d_f_opt)
        loginfo = Logging(getObjectiveValue(mdl),getValue(s[N+1]),getValue(ey[N+1]),getValue(epsi[N+1]),getValue(v[N+1]),getValue(s),getValue(ey),getValue(epsi),getValue(v),getValue(a),getValue(d_f))
        # publish commands
        if cmdcount>10      # ignore first 10 commands since MPC often stagnates during the first seconds
            publish(pub, cmd)
        end
        publish(pub2, loginfo)
        cmdcount = cmdcount + 1
        rossleep(loop_rate)
    end
end

if ! isinteractive()
    main()
end
