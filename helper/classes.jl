module classes# VARIOUS TYPES FOR CALCULATIONS
using JuMP
type LapStatus
    currentLap::Int64       # current lap number
    currentIt::Int64        # current iteration in current lap
end

# Structure of coeffConst:
# 1st dimension specifies the state (1 = eY, 2 = ePsi, 3 = v)
# 2nd dimension is the polynomial coefficient
# 3rd dimension is not used
# 4th dimension specifies one of the two lap numbers between which are iterated

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
    oldInput::Array{Float64}
    oldCost::Array{Int64}
    OldTrajectory(n_oldTraj = 0, oldTraj=Float64[],oldTrajXY=Float64[],oldInput=Float64[],oldCost=Float64[]) = new(n_oldTraj, oldTraj,oldTrajXY,oldInput,oldCost)
end

type MpcParams          # parameters for MPC solver
    N::Int64
    nz::Int64
    OrderCostCons::Int64
    Q::Array{Float64,1}
    Q_term::Array{Float64,1}
    Q_cost::Float64
    R::Array{Float64,1}
    vPathFollowing::Float64
    QderivZ::Array{Float64,1}
    QderivU::Array{Float64,1}
    MpcParams(N=0,nz=0,OrderCostCons=0,Q=Float64[],Q_term=Float64[],Q_cost=1.0,R=Float64[],vPathFollowing=1.0,QderivZ=Float64[],QderivU=Float64[]) = new(N,nz,OrderCostCons,Q,Q_term,Q_cost,R,vPathFollowing)
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
    cost::Array{Float64}
    MpcSol(a_x=0.0, d_f=0.0, solverStatus=Symbol(), u=Float64[], z=Float64[], lambda= Float64[],cost=Float64[]) = new(a_x,d_f,solverStatus,u,z,lambda,cost)
end

type Obstacle
    s_obstacle::Array{Float64}
    sy_obstacle::Array{Float64}
    rs::Float64
    ry::Float64
    index::Array{Int64,1} #!!currently not used , if that stay like this delete it
    xy_vector::Array{Float64}
    axis_y_up::Array{Float64}
    axis_y_down::Array{Float64}
    axis_s_up::Array{Float64}
    axis_s_down::Array{Float64}
    Obstacle(s_obstacle = Float64[], sy_obstacle = Float64[], rs = 0.0, ry = 0.0, index=Int64[], xy_vector=Float64[], axis_y_up=Float64[], axis_y_down=Float64[], axis_s_up=Float64[], axis_s_down=Float64[])= new(s_obstacle,sy_obstacle,rs,ry,index, xy_vector,axis_y_up, axis_y_down, axis_s_up, axis_s_down)
end

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

type MpcModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    uCurr::Array{JuMP.NonlinearParameter,1}
    sCoord_obst::Array{JuMP.NonlinearParameter,1}
    coeffTermConst::Array{JuMP.NonlinearParameter,3}
    coeffTermCost::Array{JuMP.NonlinearParameter,2}
    #s_startC::JuMP.NonlinearParameter

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    lambda::Array{JuMP.Variable,1}
    #t::Array{JuMP.Variable,1}

    dsdt::Array{JuMP.NonlinearExpression,1}
    bta::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    costPath::JuMP.NonlinearExpression
    costZ::JuMP.NonlinearExpression
    costZTerm::JuMP.NonlinearExpression
    constZTerm::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression
    laneCost::JuMP.NonlinearExpression
    costObstacle::JuMP.NonlinearExpression

    MpcModel(mdl=JuMP.Model(),
                z0=@NLparameter(mdl,z0[i=1:4]==0),
                coeff=@NLparameter(mdl,coeff[i=1:5]==0),
                uCurr=@NLparameter(mdl,zCurr[i=1:4]==0),
                sCoord_obst=@NLparameter(mdl,sCoord_obst[i=1:2]==0),
                coeffTermConst= @NLparameter(mdl,coeffTermConst[i=1:7,k=1:5,j=1:3]==0),
                coeffTermCost= @NLparameter(mdl,coeffTermCost[i=1:7,k=1:5]==0),
                #s_startC=@NLparameter(mdl, s_startC==0),
                z_Ol=@variable(mdl,[1:11, 1:4]),
                u_Ol=@variable(mdl,[1:10, 1:2]),
                lambda=@variable(mdl,[1:5]),
                #t=@variable(mdl,[1:11]),
                dsdt=@NLexpression(mdl,dsdt[1:10],0), 
                bta=@NLexpression(mdl,bta[1:10],0),
                c=@NLexpression(mdl,c[1:10],0),
                costPath = @NLexpression(mdl,costPath,0),
                costZ =@NLexpression(mdl,costZ,0),
                costZTerm=@NLexpression(mdl,costZTerm,0),
                constZTerm=@NLexpression(mdl,constZTerm,0),
                derivCost=@NLexpression(mdl,derivCost,0),
                controlCost=@NLexpression(mdl,controlCost,0),
                laneCost=@NLexpression(mdl,laneCost,0),
                costObstacle=@NLexpression(mdl,costObstacle,0))= new(mdl,
                                                        z0,
                                                        coeff,
                                                        uCurr,
                                                        sCoord_obst,
                                                        coeffTermConst,
                                                        coeffTermCost,
                                                        #s_startC,
                                                        z_Ol,
                                                        u_Ol,
                                                        lambda,
                                                        #t,
                                                        dsdt,
                                                        bta,
                                                        c,   
                                                        costPath,
                                                        costZ,
                                                        costZTerm,
                                                        constZTerm,
                                                        derivCost,
                                                        controlCost,
                                                        laneCost,
                                                        costObstacle)
end
end