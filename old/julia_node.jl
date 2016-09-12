#!/usr/bin/env julia

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

include("status.jl")
include("coeffConstraintCost.jl")


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
    pub     = Publisher("ecu", ECU, queue_size=10)
    pub2    = Publisher("logging", Logging, queue_size=10)
    s1      = Subscriber("pos_info", pos_info, SE_callback, queue_size=10)
    loop_rate = Rate(10)

    oldTraj         = OldTrajectory()
    lapStatus       = LapStatus(1,1)
    mpcCoeff        = MpcCoeff()
    posInfo         = PosInfo()
    mpcParams       = MpcParams()
    stateIn         = zeros(6,1)
    inputIn         = zeros(2,1)

    # Initialize Parameters


    while ! is_shutdown()

        # check current position and current lap
        # load previous trajectories
        # find coefficients
            # load trajectories
            # interpolate
            # find q-function
            # evaluate coefficients
        mpcCoeff = coeffConstraintCost(oldTraj, lapStatus, mpcCoeff, posInfo, mpcParams, stateIn, inputIn)
        # find optimal control input (MPC part)
        # send this input to car
        # if current lap is finished
            # save new trajectory
        
        rossleep(loop_rate)
    end
end

if ! isinteractive()
    main()
end