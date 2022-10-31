lib functions

Run:
basinSPsimComp.m:             File used to numerically evaluate the basin of attraction of the nonlinear, closed-loop simple pendulum system at the upper equilibrium.

Helper functions:
plotTree.m:           Used to plot tree policies by plotting all trajectories and, optionally, their funnels
plotFunnels.m:	      Used to plot a trajectory and, optionally, their funnels
get2DellipsePoints.m: Used to generate points on the surface of a hyper-ellipse, i.e. to generate points on the node funnel hyopthesis surface
Rot.m                 Used for visualization of the 2D Quadrotor, this mainly converted the coordinates from body frame to inertia frame

mex:
Matlab Mex-C-Code and compiled executables that are used for finding a control policy given a tree, and other functions that are used in algorithm.