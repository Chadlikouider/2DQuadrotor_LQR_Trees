% Plot all trajectories in a Sim-LQR Tree policy. Optionally also plot the
% funnels.
%
% Inputs:
%   T:          Tree struct to plot
%   options:    Plot options with elements listed below
%   .plotRange  State range to be plotted [xmin,xmax] (nx*2) nx: state dim
%   .plotDim:   For higher dimensional systems, determine which state
%               combination should be shown together in a subplot. E.g. for
%               the double pendulum, plotDim = [1 0 1 0; 0 1 0 1]': x(1), the
%               2nd rod position should be plotted together with x(3), the
%               2nd pole velocity, in the first subplot. And x(2), the 1st pendulum
%               angle, should be plotted together with x(4), the angular
%               velocity.
%   .subplot    subplot dimensions
%   .figNum:    Figure number to plot the tree
%   .nSamp:     Number of samples to use for funnel ellipses.
%   .titles     Titles for phase-plots. Cell array with same dimension as
%               size(plotDim,2).
%   .xlabels    Labels for phase-plot x/y-axis. Cell array with same dim.
%   .ylabels    size(plotDim,2).
%   .fsize      Custom font size for labels/titles
%   .basinColor Color for funnel fills, array: [r g b]
%   .trajColor  Color for trajectories, array: [r g b]
%   .linewidth  Width of trajectory lines
%   .plotBasins Enable/disable plotting of basins: true/false
%  highlightFunnel: Enable highlighting of a trajectory's funnel, where
%                   highlightFunnel is the trajectory's index
%

function [] = plotTree(T,options,highlightFunnel)
if(nargin < 3)
    highlightFunnel = 0;
end

% Number of elements in tree:
N = length(T.dt);
% state dimension:
nx = size(T.x0,1);
% visited logical array:
notVisited = [false, true(1,N-1)];

% plot all trajectories:
trajectories = [];
trajInd = 1;
while(sum(notVisited))
    ind = find(notVisited,1,'first');
    %get nominal tape and controller matrices
    trajectories(trajInd).x = T.x0(:,ind:T.endnode(ind));
    trajectories(trajInd).epsH = T.epsH(:,ind:T.endnode(ind));
    trajectories(trajInd).S = T.S(:,((ind-1)*nx + 1):(T.endnode(ind)*nx));
    notVisited(ind:T.endnode(ind)) = false;
    %add goal for connecting plot
    trajectories(trajInd).x(:,end+1) = T.x0(:,1);
    trajectories(trajInd).S = [trajectories(trajInd).S, T.S(:,1:nx)];
    trajectories(trajInd).epsH(end+1) = Inf;   % avoid replotting goal set for each traj.
    trajInd = trajInd + 1;
end

% also plot goal basin
trajectories(trajInd).x = T.x0(:,1);
trajectories(trajInd).S = [trajectories(trajInd).S, T.S(:,1:nx)];
trajectories(trajInd).epsH = T.epsH(1);

for i = 1:trajInd
    traj = trajectories(i);
    opts = options;
    % highlight funnel if enabled
    if(highlightFunnel == i)
        opts.basinColor = opts.basinColor*0.5;
        opts.trajColor = opts.trajColor*0.5;
    end
    for kTraj = 1:length(traj)
        plotFunnels(  traj(kTraj).x,...
            traj(kTraj).S,...
            traj(kTraj).epsH,...
            opts);
        figure(options.figNum)
        for j = 1:size(options.plotDim,2)
            if(~isempty(options.subplot))
                subplot(options.subplot(1),options.subplot(2),j);
            end
            if(~isempty(options.titles))
                title(options.titles{j});
            end
            if(~isempty(options.ylabels))
                ylabel(options.ylabels{j});
            end
            if(~isempty(options.xlabels))
                xlabel(options.xlabels{j});
            end
            xr = options.plotRange(repmat(options.plotDim(:,j),1,size(options.plotRange,2))>0);
            xr = reshape(xr,sqrt(numel(xr)),[]);
            axis(reshape(xr',1,[]))
        end
    end
end
end
