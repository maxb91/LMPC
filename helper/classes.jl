module classes# VARIOUS TYPES FOR CALCULATIONS
using JuMP
type LapStatus
    currentLap::Int64       # current lap number
    currentIt::Int64        # current iteration in current lap
end

# Structure of coeffConst:
# 1st dimension is the polynomial coefficient
# 2nd dimension dimension specifies the lap numbers between which are iterated
# 3rd dimension specifies the state (1 = eY, 2 = ePsi, 3 = v)

type MpcCoeff           # coefficients for trajectory approximation
    coeffCost::Array{Float64}
    coeffConst::Array{Float64}
    order::Int64
    pLength::Int64      # small values here may lead to numerical problems since the functions are only approximated in a short horizon
                        # "small" values are about 2*N, good values about 4*N
                        # numerical problems occur at the edges (s=0, when v is almost 0 and s does not change fast and at s=s_target)
    MpcCoeff(coeffCost=Float64[], coeffConst=Float64[], order=4, pLength=0) = new(coeffCost, coeffConst, order, pLength)
end

type OldTrajectory      # information about previous trajectories
    n_oldTraj::Int64
    oldTraj::Array{Float64}
    oldTrajXY::Array{Float64}
    distance2obst::Array{Float64}
    curvature::Array{Float64}
    oldInput::Array{Float64}
    oldNIter
    costs::Array{Float64}
    lambda_sol::Array{Float64}
    z_pred_sol::Array{Float64}
    u_pred_sol::Array{Float64}
    ssInfOn_sol::Array{Float64}
    eps::Array{Float64}
    cost2Target::Array{Float64}
    OldTrajectory(n_oldTraj = 0, oldTraj=Float64[],oldTrajXY=Float64[],distance2obst=Float64[],curvature=Float64[],oldInput=Float64[],oldNIter=Float64[],
        costs=Float64[],lambda_sol=Float64[],z_pred_sol=Float64[],u_pred_sol=Float64[],ssInfOn_sol=Float64[],eps=Float64[],cost2Target= Float64[]) =
                 new(n_oldTraj, oldTraj,oldTrajXY,distance2obst,curvature,oldInput,oldNIter,costs,lambda_sol,z_pred_sol,u_pred_sol,ssInfOn_sol,eps,cost2Target)
end

#old definition to load old results
# type OldTrajectory      # information about previous trajectories
#     n_oldTraj::Int64
#     oldTraj::Array{Float64}
#     oldTrajXY::Array{Float64}
#     oldInput::Array{Float64}
#     oldNIter
#     costs::Array{Float64}
#     lambda_sol::Array{Float64}
#     z_pred_sol::Array{Float64}
#     u_pred_sol::Array{Float64}
#     ssInfOn_sol::Array{Float64}
#     eps::Array{Float64}
#     cost2Target::Array{Float64}
#     OldTrajectory(n_oldTraj = 0, oldTraj=Float64[],oldTrajXY=Float64[],oldInput=Float64[],oldNIter=Float64[],
#         costs=Float64[],lambda_sol=Float64[],z_pred_sol=Float64[],u_pred_sol=Float64[],ssInfOn_sol=Float64[],eps=Float64[],cost2Target= Float64[]) =
#                  new(n_oldTraj, oldTraj,oldTrajXY,oldInput,oldNIter,costs,lambda_sol,z_pred_sol,u_pred_sol,ssInfOn_sol,eps,cost2Target)
# end

type MpcParams          # parameters for MPC solver
    N::Int64
    nz::Int64
    OrderCostCons::Int64
    Q::Array{Float64,1}
    Q_term::Array{Float64,1}
    Q_cost::Float64
    Q_obstacle::Float64
    Q_obstacleNumer::Float64
    Q_lane::Float64
    Q_velocity::Float64
    R::Array{Float64,1}
    vPathFollowing::Float64
    QderivZ::Array{Float64,1}
    QderivU::Array{Float64,1}
    MpcParams(N=0,nz=0,OrderCostCons=0,Q=Float64[],Q_term=Float64[],Q_cost=1.0,Q_obstacle = 1.0,Q_obstacleNumer = 1.0,Q_lane = 1.0,Q_velocity=1.0, R=Float64[],vPathFollowing=1.0,QderivZ=Float64[],QderivU=Float64[]) = new(N,nz,OrderCostCons,Q,Q_term,Q_cost,Q_obstacle,Q_obstacleNumer,Q_lane, Q_velocity, R,vPathFollowing)
end

type PosInfo            # current position information
    s_start::Float64
    s::Float64
    s_target::Float64
    PosInfo(s_start=0,s=0,s_target=0) = new(s_start,s,s_target)
end

type MpcSol             # MPC solution output
    a_x::Float64
    d_f::Float64
    solverStatus::Symbol
    u::Array{Float64}
    z::Array{Float64}
    lambda::Array{Float64,1}
    ssInfOn::Array{Int64,1}
    eps::Array{Float64}
    cost::Array{Float64}
    MpcSol(a_x=0.0, d_f=0.0, solverStatus=Symbol(), u=Float64[], z=Float64[], lambda= Float64[], ssInfOn= Int64[],eps=Float64[],cost=Float64[]) = new(a_x,d_f,solverStatus,u,z,lambda,ssInfOn,eps,cost)
end

type Obstacle
    s_obstacle::Array{Float64}
    sy_obstacle::Array{Float64}
    rs::Float64
    ry::Float64
    v::Array{Float64}
    index::Array{Int64,1} #!!currently not used , if that stay like this delete it
    xy_vector::Array{Float64}
    axis_y_up::Array{Float64}
    axis_y_down::Array{Float64}
    axis_s_up::Array{Float64}
    axis_s_down::Array{Float64}
    Obstacle(s_obstacle = Float64[], sy_obstacle = Float64[], rs = 0.0, ry = 0.0,v= Float64[], index=Int64[], xy_vector=Float64[], axis_y_up=Float64[], axis_y_down=Float64[], axis_s_up=Float64[], axis_s_down=Float64[])= new(s_obstacle,sy_obstacle,rs,ry,v,index, xy_vector,axis_y_up, axis_y_down, axis_s_up, axis_s_down)
end

##########old definition o load old results
# type Obstacle
#     s_obstacle::Array{Float64}
#     sy_obstacle::Array{Float64}
#     rs::Float64
#     ry::Float64
#     v::Float64
#     index::Array{Int64,1} #!!currently not used , if that stay like this delete it
#     xy_vector::Array{Float64}
#     axis_y_up::Array{Float64}
#     axis_y_down::Array{Float64}
#     axis_s_up::Array{Float64}
#     axis_s_down::Array{Float64}
#     Obstacle(s_obstacle = Float64[], sy_obstacle = Float64[], rs = 0.0, ry = 0.0,v= 0.4, index=Int64[], xy_vector=Float64[], axis_y_up=Float64[], axis_y_down=Float64[], axis_s_up=Float64[], axis_s_down=Float64[])= new(s_obstacle,sy_obstacle,rs,ry,v,index, xy_vector,axis_y_up, axis_y_down, axis_s_up, axis_s_down)
# end

type TrackCoeff         # coefficients of track
    coeffAngle::Array{Float64,1}
    coeffCurvature::Array{Float64,1}
    nPolyCurvature::Int64      # order of the interpolation polynom
    nPolyXY::Int64              # order of the interpolation polynom of the x y coordinates
    width::Float64               # lane width -> is used in cost function as soft constraints (to stay on track)
    ds::Rational{Int}
    TrackCoeff(coeffAngle=Float64[], coeffCurvature=Float64[], nPolyCurvature=4, nPolyXY = 6, width=1.0, ds=1//10) = new(coeffAngle,coeffCurvature,nPolyCurvature,nPolyXY,ds)
end

type ModelParams
    l_A::Float64
    l_B::Float64
    dt::Float64
    u_lb::Array{Float64}        # lower bounds for u
    u_ub::Array{Float64}        # upper bounds
    z_lb::Array{Float64}
    z_ub::Array{Float64}
    c0::Array{Float64}
    ModelParams(l_A=0.25,l_B=0.25,dt=0.1,u_lb=Float64[],u_ub=Float64[],z_lb=Float64[],z_ub=Float64[],c0=Float64[]) = new(l_A,l_B,dt,u_lb,u_ub,z_lb,z_ub,c0)
end


end
