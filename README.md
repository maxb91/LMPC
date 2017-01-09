# LMPC

This repo provides a framework to run simple ideal-model LMPC simulations in julia.

How to run simulations:
-------------------------
1. open julia terminal in this folder
2. julia> include("testNode.jl")
3. julia> run_sim()

Further info:
-------------------------
- As of Sep. 2016 there are only constant curvatures available.
- The model is defined by global variables.
- There's a second file ("testNode_profiling.jl") which can be used to run the julia profiler and check what functions take long time to evaluate.
- There's a folder ("Matlab and C") which contains a file with LMPC implemented in a C-function. It was only used to derive cost functions for the LMPC and should stay there (to look up the functions)

Branches:
---------
simple-sim: Kinematic model for path following and Learning MPC, no system ID ("perfect mathematical equivalent models")
dyn-sim: Kinematic model for path following, dynamic model for LMPC (same system for simulation and for prediction)
master: Kinematic model for path following, dynamic model with system ID for LMPC