% Generate a Sim-LQR Tree policy for the planar Quadrotor system with 
% obstacle avoidance

function [] = generateQuad2D()
% add mex files and other functions to path
addpath(genpath('./../lib'))

% System parameters of experimental planar Quadrotor setup (the measuements are taken from some paper online)

% m: mass of the Quadrotor (Kg);
% g: gravity (m/s^2);
% r: distance from the center to the base of the propellor (m);
% I: moment of inertial (kg.m^2)

g = 9.81;
% parameters of the Quadroror
m = 0.486; % Kg
r = 0.25;  % meter
I = 0.00383; % Kg.m^2
% state and input dims
nx = 6;
nu = 2;

% goal state:
goal.x = [-3,2,0,0,0,0]';             % nominal goal state
goal.u = [m*g/2,m*g/2]';             % nominal goal input

% TILQR weights for balancing the Quadrotor at the goal Xg
goal.Qtilqr = diag([10,10,10,1,1,r/(2*pi)]);  % state penalty
goal.Rtilqr = [0.1 0.05;0.05 0.1];            % input penalty
goal.dt = 0.15;                               % controller sampling time
goal.bias = 1;
% Get A,B of continuous-time dynamics
[~,df] = dynamicsgr(0,goal.x,goal.u);

% Get discrete-time LTI dynamics using exact discretization (matrix exp)
sysd = c2d( ss(df{1}(:,1+(1:nx)),...
            df{1}(:,nx+1+(1:nu)),...
            ones(1,nx),zeros(1,nu)),...
            goal.dt);

% get TILQR policy (discretize cost matrices)
[goal.Ktilqr, goal.Stilqr] = dlqr(sysd.A, sysd.B,...
                                goal.dt*goal.Qtilqr, goal.dt*goal.Rtilqr);
% LQR weights for TVLQR policy, obtained empirically by tuning controller
% for stabilizing the trajectory from a relatively far position from the
% goal state
options.tvlqr.Q = diag([20,20,5000,2000,2000,100]);
options.tvlqr.R = goal.Rtilqr;
options.tvlqr.Qf = goal.Stilqr; % Final cost matrix


% Init Motion Planning Options
dirCopts.Tmax = 15;  	 % Maximal trajectory time
dirCopts.Tmin = 0.01;    % Minimal     "
limits = [-4.9 4.9;         % state constraints. Shrunk state constraints
          0.25 6;          % in order to allow for some control authority of
         -pi/4 pi/4;   % TVLQR policy along the trajectories
         -2 2;        % order: x,z,theta,dx,dz,dtheta
         -1 1;
         -pi/3 pi/3];
dirCopts.Umax = [10;10];      % same for input constraints. Choose smaller   
dirCopts.Umin = [-10;-10];     % constraints for leaving TVLQR control authority
dirCopts.Xmin = limits(:,1);
dirCopts.Xmax = limits(:,2);
dirCopts.maxNutape = 30;   % max. # of samp. allowed in motion plan. traj.
dirCopts.minDT = 0.0001;     % min. samp. time allowed of motion plan. traj.
dirCopts.maxDT = 0.9;      % max.              "
dirCopts.defaultDT = goal.dt; % default samp. time of motion plan. traj.
dirCopts.resampleToDefaultDT = true; % activate resamp. of traj. (set false when doing RRT bias<100%)
dirCopts.nDircolTrajectory = 25;   % default number of samples for traj.
dirCopts.nColTraj=2;           % just used in RRT (collocated point in the steering function)
options.dirCopts = dirCopts; % copy motion planning options to options
options.maxNodesRRT=1e4;     % the maximum RRT nodes

% Set parameters of termination heuristic and max. number of iterations,
options.termCond.alpha = 0.01;  % event probability parameter
options.termCond.pAlpha = 0.99; % success rate parameter
options.maxIter =  1e4;   % max. number of random samples / iterations

% STATE/INPUT LIMITS used in simulations for funnel falsifications
options.xLimit = [-Inf*ones(6,1), Inf*ones(6,1)];
%add state constraints in between the two
options.uLimit = [-25, 25;
                  -25, 25];             % control input is limited to 50 Newton

if(1) % enable to disable state constraints
    options.xLimit(1,:) = [-5.5 5.5];
    options.xLimit(2,:) = [0.25 6.5];
    options.xLimit(3,:) = [-pi/2 pi/2];
end
if(1) % enable to add Obstacles
    options.obstacles = [2.5, 3.8, 1;    % x coordinate, y coordinate, radius
                         0.5,  1, 0.8;
                         -1.5, 3, 0.7];
end

% design region D to be covered by LQR-Tree policy:
xRange = [-4 4;              % x
          0.3 5;              % z
          -pi/4 pi/4;        % theta
          -2 2;              % x dot
          -1 1;              % z dot
          -pi/3 pi/3];       % theta dot

% Init plotting options
options.plotOpt = initPlotOpts();

% basin of attraction approximation parameters
goal.basin.nEvals = 2e3;    % number of samples to evaluate basin
                            % by approximation

% Set termination heuristic parameters for basin approximation
goal.basin.termCond.alpha = 0.01;    % event probability param.
goal.basin.termCond.pAlpha = 0.99;   % success rate parameter
% goal.basin.epsH = 1.1215; for xg=[1 0 0 0 0 0]'                            
% enable to pause after goal basin approximation
options.promptUserAfterGoalApprox = false;

% if(1)   % resume the planning    (uncomment to run)
%     seed = load('name_of_file.mat');
%     options.tvlqr.Q = seed.options.tvlqr.Q;
%     options.tvlqr.R = seed.options.tvlqr.R;
%     goal.Rtilqr = seed.goal.Rtilqr;
%     goal.Qtilqr = seed.goal.Qtilqr;
%     goal.Stilqr = seed.goal.Stilqr;
%     options.seedTrajs(1).tf = seed.tf;
%     options.seedTrajs(1).xdir = seed.xdir;
%     options.seedTrajs(1).udir = seed.udir;
% end

if(1) % enable to save files to folder
    options.nplotPDF = 50;
    options.saveFolder = [datestr(now,30),'_2DQuadrotor'];
    if (isfield(options,'obstacles'))
        options.saveFolder = [options.saveFolder,'_WithObstacles'];
    end
    if(all(all(isinf(options.xLimit)))) % adapt folder name if no state constraints
        options.saveFolder = [options.saveFolder,'_NoStateConstraints'];
    end
end


% reset random number generator to get the results presented in the thesis:
 reset(RandStream.getGlobalStream, 0); % A fixed seed

% get policy
ts = tic;
[T,opts] = simLQRtree(@dynamicsgr,@cost,xRange,goal,options);
ttree = toc(ts) % the time needed to generate the tree 

% save tree generation time
if(isfield(opts,'saveFolder'))
    save(['./', opts.saveFolder,'/mat/ttree.mat'],'ttree');
end

% Dynamics of the planar quadrotor
% xdot: state time derivative given time, state, and ctrl. input
% df:   gradients of dyn. with resp. to time, state, and ctrl. input
    function [xdot, df] = dynamicsgr(t,x,u)
        
        % Equation of motion
        % written in this form to be usable with fmincon
        th = x(3,:);
        dx = x(4,:);
        dz = x(5,:);
        dth = x(6,:);
%         u1 = u(1,:);  % rotor 1 force magnitude
%         u2 = u(2,:);  % rotor 2 force magnitude

        xdot=[dx;dz;dth;-(u(1,:)+u(2,:)).*sin(th)/m;(u(1,:)+u(2,:)).*cos(th)/m-g;...
            (u(1,:)-u(2,:))*r/I];

        if (nargout>1)
          
            dfdx = [0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;0,0,-(u(1,:)+u(2,:)).*cos(th)/m,...
                    0,0,0;0,0,-(u(1,:)+u(2,:)).*sin(th)/m,0,0,0;0,0,0,0,0,0];
                
            dfdu=[0,0;0,0;0,0;-sin(th)/m,-sin(th)/m;cos(th)/m,cos(th)/m;r/I,-r/I];
            
            dfdt = [ 0; 0; 0; 0; 0; 0];
            df{1} = [dfdt, dfdx, dfdu];
        end
    end

% instantaneous cost for motion planning
% output:
% g: instantaneous cost
% dg: gradient of cost with respect to time, state, and ctrl. input
    function [g,dg] = cost(t,x,u)
        dx = x-goal.x;
        du = u-goal.u;
        g = dx'*options.tvlqr.Q*dx + du'*options.tvlqr.R*du;
        
        if (nargout>1)
            dgdt = 0;
            dgdx = 2*dx'*options.tvlqr.Q;
            dgdu = 2*du'*options.tvlqr.R;
            dg{1} = [dgdt,dgdx,dgdu];
        end
    end

    function plotOpts = initPlotOpts()
        plotOpts.figNum = 122;
        plotOpts.subplot = [3,1];                    % subplot dimensions
        plotOpts.titles{1} = 'Planar Quadrotor ';    % subplot titles
        plotOpts.titles{2} = ' ';
        plotOpts.titles{3} = ' ';
        plotOpts.xlabels{1} = 'x (m)';             % axis labels
        plotOpts.ylabels{1} = 'dx/dt (m/s)';
        plotOpts.xlabels{2} = 'z (m)';
        plotOpts.ylabels{2} = 'dz/dt (m/s)';
        plotOpts.xlabels{3} = '\theta (rad)';
        plotOpts.ylabels{3} = 'd\theta/dt (rad/s)';
        plotOpts.fsize = 12;                    % font size for axis desc.
        plotOpts.plotDim = [[1,0,0,1,0,0]',[0,1,0,0,1,0]',... % inds. of state dims.
                            [0,0,1,0,0,1]'];                 % couples to plot
                                                    
        plotOpts.basinColor = [0 0.5 1]; % color of basins [r g b]
        plotOpts.trajColor = [1 0 0];% color of trajectories [r g b]
        plotOpts.linewidth = 2;      % width of trajectory lines
        plotOpts.plotBasins = false; % plot the funnels on/off
        options.plotIncrement = true;% avoid replotting the tree (slow)
                                     % by just incrementally adding
                                     % trajectories to plot
                                     % region to be plotted in plot of
                                     % trajectories in tree
        plotOpts.plotRange = [-6 6;              % x
          0 7;                                  % z
          -pi/3 pi/3;                            % theta
          -5 5;                                  % x dot
          -4 4;                                  % z dot
          -pi/2 pi/2];                           %theta dot
    end
end