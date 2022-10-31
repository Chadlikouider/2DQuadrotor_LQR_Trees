# 2DQuadrotor_LQR_Trees
This project is a MATLAB implementation of the simulation based LQR-Tree algorithm for a planner quadrotor.

## Prerequisites

In order to run this program you will need to have the following installed on your machine:
* [Optimtraj](https://github.com/MatthewPeterKelly/OptimTraj)

and extract the folder to Generate_trees folder

## steps
* Generate the LQR-trees by running generateQuad2D file in Generate_trees folder
* run Quadrotor2D_withObstacles or Quadrotor2D to visualize and test the LQR-Trees generated

## Future Work

* Improve code documentation
* Remove extranneous code
* improve the running time

## Issues for generating LQR-Trees
* the running time for a 2d state space model is around 30min
* runing time for a 4d state space model is around 60 hours
* runing time for a 6d state space model is over 200 hours

## Aknowledgments

This software was developed as a part of my final graduation project.
