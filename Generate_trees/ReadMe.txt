Generate Tree

Use:
generateSP.m	        Run to generate an LQR-Tree policy for the simple pendulum.
generateCP.m	        Run to generate an LQR-Tree policy for the cart-pole.
generateQuad2D.m	Run to generate an LQR-Tree policy for the Quadrotor.

Functions:
dltv_lqr.m            -generates the TVLQR policy for a given trajectory.
dir_col.m             -used for motion planning using direct collocation with direct transcription.
simLQRtree.m	      -the algorithm that generates the trajectories (nominal states,nominal inputs, K(tk),S(tk),funnel adjustments....ect).
LQR_base_proximity.m  -the function approximate the closeness of the states which performing the RRT 
collocate_trajectory.m -the procedure is the proposed steering function used in the RRT
RRT.m                  - the kino-dynanamic RRT algorithm


Notice:
Add OptimTraj to the path before running the program with 100% goal bias
