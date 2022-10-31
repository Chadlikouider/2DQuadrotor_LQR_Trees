% Run a Sim-LQR Tree policy generated for the planar Quadrotor system in
% simulation with random initial conditions.
%
% You can also evaluate the quality measures using this function by
% setting: 
% nsims = 2000;       
% showsim = 0;        % disable animation
% showFails = 0;      % (optional) show failed simulations (1) or (0)
% signif = 0.01;      % confidence interval significance level (see binofit)
%
%

function [] = Quadrotor2D_withObstacles()
% add helper functions and mex files
addpath(genpath('./../lib'))
% Function parameters:
nsims = 5;   % number of simulations to run
showsim = 1;    % show animation (1) or not (0)
showTree = 0;   % plot all trajectories in the tree
signif = 0.01;  % confidence interval significance level (see binofit)
showFails = 1;  % show failed simulations (1) or (0)
tryClosestNode = 0; % switch to enable using the policy of the closest node
                    % (=1) in case the sample is not in any funnel hyp.

% String of tree policy folder name to load. Folder should be located
% in Generate_Tree subfolder:
treeFile = '20210826T173851_2DQuadrotor_WithObstacles';

% parameters of the Quadroror

% m: mass of the Quadrotor (Kg);
% g: gravity (m/s^2);
% r: distance from the center to the base of the propellor (m);
% I: moment of inertial (kg.m^2)

g = 9.81;
m = 0.486;
r = 0.25;
I = 0.00383;
% The tree was generated for m = 0.486, g = 9.81, r = 0.25, I=0.00383;

% reset random number generator to get the results presented in the thesis:
% reset(RandStream.getGlobalStream, 0);

% Simulation speed parameter:
disp_dt = 0.030; % increase for faster simulation on slower computers
tUp = 3;         % number of seconds to keep simulating time-invariant goal
                 % state controller

% load the mat file containing the tree policy:
s = load(['./../Generate_Tree/',treeFile,'/mat/',treeFile,'.mat']);

% set appropriate title for plotting
s.options.plotOpt.titles{1} = 'LQR-Tree Policy for Planar Quadrotor with obstacles';

% plot the tree, if enabled
if(showTree)
    try close(s.options.plotOpt.figNum); catch end
    plotTree(s.Ts,s.options.plotOpt);
end

% init some global function parameters:
nX = 6;     % state dimension
nU = 2;     % input dimension
nT = length(s.Ts.dt);   % number of nodes in tree

% init failed samples storage arrays
xfail = [];

if(showsim)
    disp('Click/Press Key to Continue')
end

% init statistic counters
nIn = 0;    % number of samples in any basin hypothesis in tree
nSuc = 0;   % number of successful simulations of samples that are in
            % any basin hyp. in the tree
nOutSuc = 0;% number of samples not in any hypothesis, but succ. stab. to
            % goal by "closest" node
nFailInGoal = 0; % number of samples that were stabilized to goal set by
                 % tree policy, but TILQR policy did not provide
                 % convergence
nFailSC = 0;    % number of samples that failed due to state constraints

% start simulations
for i = 1:nsims
    % get random sample from design region D:
    xsamp = getRandSample(s.xRange);
    
        [ind,flg,indMin] = findPolicyInt(xsamp);
    
    % flag == 1  -> found tree node ind with minimal cost-to-go,
    %               whose funnel contains sample
    % flag == -1 -> the random sample is outside all funnels in the tree
    if(flg < 0)
            disp('sample outside all funnels in tree!')
    else
        nIn = nIn + 1;
    end
    
    if(flg > 0 || tryClosestNode)
        % run simulation with found policy, and time-invariant LQR policy after
        % the time-varying policy ends:
        [xtape,tsim,x0,ln] = simTree();
        % check for success and show fail, if necessary
        fail = norm(xtape(:,end) - s.Ts.x0(:,1)) > 0.1;
        violatedConstraints = ~testConstraints(xtape);
        collision=~testCollision(xtape(1:2,:));
        nFailSC = nFailSC + violatedConstraints + collision;
        fail = fail || violatedConstraints || collision;
    else
        % sample outside all funnels and closest node policy should not be
        % tried, continue on to next sample
        continue
    end
    % check if the system state was in the goal region, when the goal state
    % TILQR policy was activated, excluding state-constraints fails
    dxG = xtape(:,ln+1) - s.Ts.x0(:,1);
    failGoal = dxG'*s.Ts.S(:,1:nX)*dxG < s.Ts.epsH(1) && testConstraints(xtape);

    % show simulation if showsim or if we want to show the simulation of a
    % failed sample:
    if(showsim ||(fail && showFails))
        % resample state trajectory to regular time interval given by
        % disp_dt:
        xrs = zeros(nX,floor(tsim(end)/disp_dt)+1);
        x0 = s.Ts.x0(:,ind:s.Ts.endnode(ind));
        tdisp = 0:disp_dt:(floor(tsim(end)/disp_dt)*disp_dt);
        for dim = 1:nX
            xrs(dim,:) = spline(tsim,xtape(dim,:),tdisp);
        end
        
        % show animation
        for j = 1:length(tdisp)
            draw(tdisp(j),xrs(:,j));
        end
        
        %plot the resulting trajectories
        figure(31)
        subplot(3,1,1), hold off
        plot(x0(1,:),x0(4,:),'Color','r','LineWidth',s.options.plotOpt.linewidth), hold on
        plot(xtape(1,:),xtape(4,:),'LineWidth',s.options.plotOpt.linewidth)
        plot(xtape(1,1),xtape(4,1),'Color','g','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none')
        plot(x0(1,end),x0(4,end),'Color','k','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none'), hold off
        legend('Nominal Trajectory','Simulated Trajectory', 'Random Sample', 'Goal State')
        title(['Planar Quadrotor Phaseplot, click/press key to continue'])
        xlabel(s.options.plotOpt.xlabels{1});
        ylabel(s.options.plotOpt.ylabels{1});
        subplot(3,1,2), hold off
        plot(x0(2,:),x0(5,:),'Color','r','LineWidth',s.options.plotOpt.linewidth), hold on
        plot(xtape(2,:),xtape(5,:),'LineWidth',s.options.plotOpt.linewidth)
        plot(xtape(2,1),xtape(5,1),'Color','g','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none')
        plot(x0(2,end),x0(5,end),'Color','k','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none'), hold off
        legend('Nominal Trajectory','Simulated Trajectory', 'Random Sample', 'Goal State')
        xlabel(s.options.plotOpt.xlabels{2});
        ylabel(s.options.plotOpt.ylabels{2});
        subplot(3,1,3), hold off
        plot(x0(3,:),x0(6,:),'Color','r','LineWidth',s.options.plotOpt.linewidth), hold on
        plot(xtape(3,:),xtape(6,:),'LineWidth',s.options.plotOpt.linewidth)
        plot(xtape(3,1),xtape(6,1),'Color','g','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none')
        plot(x0(3,end),x0(6,end),'Color','k','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none'), hold off
        legend('Nominal Trajectory','Simulated Trajectory', 'Random Sample', 'Goal State')
        xlabel(s.options.plotOpt.xlabels{3});
        ylabel(s.options.plotOpt.ylabels{3});
        drawnow
    end
    
    %check if sim is successful and show IC in appropriate color
    if(fail)
        pcolor = 'r';
        % store failed sample and flag:
        xfail(:,end+1) = [xsamp;flg];
        % check if state was in goal after tree policy ended:
        nFailInGoal = nFailInGoal + double(failGoal);
    else
        pcolor = 'g';
        nSuc = nSuc + 1;
        nOutSuc = nOutSuc + double(flg < 0);
    end
    
    if(showsim)
        % plot IC on tree plot
        figure(s.options.plotOpt.figNum)
        subplot(3,1,1), hold on
        plot(xsamp(1),xsamp(4),'Color',pcolor,'Marker','o','LineWidth',2)
        xlabel(s.options.plotOpt.xlabels{1});
        ylabel(s.options.plotOpt.ylabels{1});
        title('Initial Conditions. Green = Success, Red = Fail, Cross = Not in Funnel ')
        subplot(3,1,2), hold on
        plot(xsamp(2),xsamp(5),'Color',pcolor,'Marker','o','LineWidth',2)
        xlabel(s.options.plotOpt.xlabels{2});
        ylabel(s.options.plotOpt.ylabels{2});
        subplot(3,1,3), hold on
        plot(xsamp(3),xsamp(6),'Color',pcolor,'Marker','o','LineWidth',2)
        xlabel(s.options.plotOpt.xlabels{3});
        ylabel(s.options.plotOpt.ylabels{3});
        if(flg < 0)
            subplot(3,1,1)
            plot(xsamp(1),xsamp(4),'Color',pcolor,'Marker','x','LineWidth',2)
            subplot(3,1,2)
            plot(xsamp(2),xsamp(5),'Color',pcolor,'Marker','x','LineWidth',2)
            subplot(3,1,3)
            plot(xsamp(3),xsamp(6),'Color',pcolor,'Marker','x','LineWidth',2)
        end
        waitforbuttonpress();
    end
    if(nsims > 50 && mod(i,floor(nsims/10)) == 0)
        disp(['i = ',int2str(i),'/',int2str(nsims),...
            '    nsuc = ',int2str(nSuc),' of nIn = ',int2str(nIn)])
    end
end

% show failed initial conditions:
if(~isempty(xfail) && showFails)
    figure(111)
    subplot(3,1,1), hold off
    plot(xfail(1,xfail(7,:)>0),xfail(4,xfail(7,:)>0),'Color','r','Marker','x','LineWidth',2,'LineStyle','none','MarkerSize',7), hold on
    plot(xfail(1,xfail(7,:)<0),xfail(4,xfail(7,:)<0),'Color','b','Marker','o','LineWidth',2,'LineStyle','none','MarkerSize',7), hold off
    xlabel(s.options.plotOpt.xlabels{1});
    ylabel(s.options.plotOpt.ylabels{1});
    legend('Failed, IC was inside funnel', 'Failed, IC was outside funnel')
    title('Failed Samples')
    subplot(3,1,2), hold off
    plot(xfail(2,xfail(7,:)>0),xfail(5,xfail(7,:)>0),'Color','r','Marker','x','LineWidth',2,'LineStyle','none','MarkerSize',7), hold on
    plot(xfail(2,xfail(7,:)<0),xfail(5,xfail(7,:)<0),'Color','b','Marker','o','LineWidth',2,'LineStyle','none','MarkerSize',7), hold off
    xlabel(s.options.plotOpt.xlabels{2});
    ylabel(s.options.plotOpt.ylabels{2});
    subplot(3,1,3), hold off
    plot(xfail(3,xfail(7,:)>0),xfail(6,xfail(7,:)>0),'Color','r','Marker','x','LineWidth',2,'LineStyle','none','MarkerSize',7), hold on
    plot(xfail(3,xfail(7,:)<0),xfail(6,xfail(7,:)<0),'Color','b','Marker','o','LineWidth',2,'LineStyle','none','MarkerSize',7), hold off
    xlabel(s.options.plotOpt.xlabels{3});
    ylabel(s.options.plotOpt.ylabels{3});
end

% calculate statistics:
[covRatio, covRatioCI] = binofit(nIn, nsims, signif);
[sucRate, sucRateCI] = binofit(min(nSuc,nIn), nIn, signif);
disp(['Estimated Success Rate = ',num2str(sucRate),...
    ', with CI (',num2str(sucRateCI(1)),', ', num2str(sucRateCI(2)),...
    ') (alpha = ',num2str(signif),').'])
disp(['Estimated Coverage Ratio = ',num2str(covRatio),...
    ', with CI (',num2str(covRatioCI(1)),', ', num2str(covRatioCI(2)),...
    ') (alpha = ',num2str(signif),').'])
disp(['Term. Heuristic Params: pAlpha = ',...
    num2str(s.options.termCond.pAlpha),', alpha = ',...
    num2str(s.options.termCond.alpha)])
disp(['Tree File: ', treeFile])
disp(['Number of samples where system state ends up in goal, ',...
    'but convergence fails = ',int2str(nFailInGoal)])
disp(['Number of samples that failed due to state constraints = ',int2str(nFailSC)])

    function xSample = getRandSample(xRange)
        %uniformly sample from xRange:
            collision_free=1; % 1 sample in obstacle space
                              % 0 sample in the free space
            num_obs = size(s.options.obstacles,1); % number of obstacles
            radi = zeros(num_obs,1);
            while collision_free
                xSample = xRange(:,1) ...
                    + rand(size(xRange,1),1).*(xRange(:,2) - xRange(:,1));
                pos = xSample(1:2);  % just the x and z configuration for collision detection
                for p=1:num_obs
                    radi(p)=sqrt((pos(1)-s.options.obstacles(p,1))^2+(pos(2)-s.options.obstacles(p,2))^2);
                end
                if all(radi > s.options.obstacles(:,3)+r)
                    collision_free =0;
                end
            end
    end

    function [xtape,tsim,x0,lnom] = simTree()

        % simulate random sample with policy starting at ind
        % return resulting state trajectory (xtape) and nominal trajectory
        % to follow (x0).
        
        % check if sample is in goal basin:
        if(ind < 2)
            % we simulate for tUp:
            nsim = floor(tUp/s.Ts.dt(1));   % no. of time steps in sim. traj.
            lnom = 0;
            tsim = 0:s.Ts.dt(1):(nsim * s.Ts.dt(1)); % simulation time
            x0 = repmat(s.Ts.x0(:,1),1,nsim + 1); % nominal traj = goal state
        else
            % if there is a nonzero length trajectory to follow, prepare
            % it, and the LTV policy that stabilizes it:
            x0 = s.Ts.x0(:,ind:s.Ts.endnode(ind));
            u0 = s.Ts.u0(:,ind:s.Ts.endnode(ind));
            lnom = size(u0,2);
            % TVLQR compensator matrices
            Ksim = s.Ts.K(((ind-1)*nU+1):(s.Ts.endnode(ind)*nU),:);
            % Number of time steps in sim. traj. (we add tUp at the end to
            % check for actual convergence to the goal state).
            nsim = size(u0,2) + floor(tUp/s.Ts.dt(1));
            tsim = cumsum([0,s.Ts.dt(1,ind:s.Ts.endnode(ind)),ones(1,floor(tUp/s.Ts.dt(1)))*s.Ts.dt(1)]); % simulation time
            % append nominal state
            x0 = [x0, repmat(s.Ts.x0(:,1),1,nsim - size(x0,2) + 1)];
        end
        
        % input constraints:
        umax = Inf;
        umin = -umax;
        if(isfield(s.options,'uLimit'))
            umax = s.options.uLimit(:,2);
            umin = s.options.uLimit(:,1);
        end
        
            % init xtape, utape, where simulation state and input
            % trajectories are stored
            xtape = zeros(size(x0));
            utape = zeros(nU,nsim);
            % set first element to random state sample
            xtape(:,1) = xsamp;
            
            % simulate system for lnom steps
            for k = 1:nsim
                if(k > lnom)
                    % get TILQR input
                    u = s.Ts.u0(:,1) - s.Ts.K(1:nU,1:nX)*(xtape(:,k)-s.Ts.x0(:,1));
                else
                    % get TVLQR input:
                    u = u0(:,k)- s.Ts.K(1:nU,1:nX)*(xtape(:,k)-x0(:,k));% - Ksim(((k-1)*nU+1):(k*nU),:)*(xtape(:,k)-x0(:,k));
                end
                % enforce input constraints
                u = min(umax,max(u,umin));
                utape(:,k) = u;
                [~,xode] = ode45(@(t,x)fofx(t,x),tsim(k):tsim(k+1)-tsim(k):tsim(k+1),xtape(:,k));
                xtape(:,k+1) = xode(end,:)';
            end
        
        if(lnom == 0)
            lnom = floor(nsim/2);
            Ksim = repmat(s.Ts.K(1:nU,1:nX),lnom,1);
            u0 = repmat(s.Ts.u0(:,1),lnom,1);
        end
        % create init trajectory:
        init.xtape = xtape(:,1:lnom+1);
        init.tsim = tsim(1:lnom+1);
        init.utape = zeros(nU,lnom);
        %recover input tape (for init trajectories):
        for k = 1:lnom
            u = u0(:,k) - s.Ts.K(1:nU,1:nX)*(xtape(:,k)-x0(:,k));%- Ksim(((k-1)*nU + 1):k*nU,:)*(xtape(:,k) - x0(:,k));
            if(~isempty(s.options.uLimit))
                u = min([u,s.options.uLimit(:,2)],[],2);
                u = max([u,s.options.uLimit(:,1)],[],2);
            end
            init.utape(:,k) = u;
        end
        
        function [xdot] = fofx(t,x)
            xdot=dynamicsgr(t,x,u);
        end
        
    end

    function draw(t,x)
        persistent hFig len base wb pin prop obstacle_one obstacle_two obstacle_three
        if (isempty(hFig))
             hFig = figure(33);
             set(hFig,'DoubleBuffer', 'on');
             len = 0.25;     % moment arm (meters)
             hb = 0.025;     % hight of the base
             wb = 1.2*len;   % width of the base
             wp=0.005;       % the width of the pin
             hp=0.1;         % the hight of the pin
             base = [wb*[1 -1 -1 1 1];hb*[1 1 -1 -1 1]]'; % stack vertically
             pin = [wp*[1 1 -1 -1 1];hp*[1 0 0 1 1]]';
             a = linspace(0,2*pi,50);
             prop = [len/1.5*cos(a); hp+0.02*sin(2*a)]'; % make the propeler
             if(1)
             obstacle_one=[s.options.obstacles(1,1)+s.options.obstacles(1,3)*cos(a);s.options.obstacles(1,2)+s.options.obstacles(1,3)*sin(a)]';
             obstacle_two=[s.options.obstacles(2,1)+(s.options.obstacles(2,3)-0.4)*cos(a);s.options.obstacles(2,2)+(s.options.obstacles(2,3)-0.4)*sin(a)]';
             obstacle_three=[s.options.obstacles(3,1)+s.options.obstacles(3,3)*cos(a);s.options.obstacles(3,2)+s.options.obstacles(3,3)*sin(a)]';
             end
        end
        figure(hFig); cla; hold on;
        set(hFig,'Position', [100 100 640 480]);
        R=[cos(x(3)) -sin(x(3)); sin(x(3)) cos(x(3))];      % rotation matrix
        R_base = Rot(R,base);                               % the dot product
        patch(x(1)+R_base(:,1),x(2)+R_base(:,2),'k','FaceColor',[.6 .6 .6])
        R_lpin = Rot(R,[-len+pin(:,1),pin(:,2)]);                           % rotate left pin 
        patch(x(1)+R_lpin(:,1), x(2)+R_lpin(:,2),'k','FaceColor',[0 0 0])   % left pin
        R_rpin = Rot(R,[len+pin(:,1),pin(:,2)]);                            % rotate right pin 
        patch(x(1)+R_rpin(:,1), x(2)+R_rpin(:,2),'k','FaceColor',[0 0 0])   % right pin
        R_lprop = Rot(R,[-len+prop(:,1),prop(:,2)]);                        % rotate left prop
        patch(x(1)+R_lprop(:,1), x(2)+R_lprop(:,2),'k','FaceColor',[0 0 1]) % left prop
        R_rprop = Rot(R,[len+prop(:,1),prop(:,2)]);                         % rotate right prop
        patch(x(1)+R_rprop(:,1), x(2)+R_rprop(:,2),'k','FaceColor',[0 0 1]) % right prop
        if(1)
            patch(obstacle_one(:,1),obstacle_one(:,2),0*obstacle_one(:,1),'k');
            patch(obstacle_two(:,1),obstacle_two(:,2),0*obstacle_two(:,1),'k');
            patch(obstacle_three(:,1),obstacle_three(:,2),0*obstacle_three(:,1),'k');
        end
        plot(-3,2,'Color','g','MarkerSize',10,'LineWidth',2,'Marker','x','LineStyle','none')
        title(['t = ', num2str(t,'%.2f'),' sec'],'FontSize',22);%, x0str]);
        set(gca,'XTick',[],'YTick',[])
        axis image; axis([-5 6 0 6.5]);
        drawnow;
    end

    function [ind,flag,indMinAll] = findPolicyInt(x)
        % init flag, which indicates if the sample is within any funnel or
        % not (1: within a funnel, -1: outside all funnels)
        flag = -1;
        % evaluate TILQR cost-to-go and use it to init the overall minimal
        % cost variable (if the sample is outside all funnels, use policy
        % starting at node with minimal cost-to-go).
        jminall = x'*s.Ts.S(:,1:nX)*x;
        indMinAll = 1;
        
        % init minimal cost of nodes whose funnel contains the sample x
        jmin = Inf;
        
        % set starting node index to goal node
        ind = 1;
        
        % if the sample is within the goal basin, finish.
        if(jminall < s.Ts.epsH(1))
            flag = 1;
            return
        end
        
        % now check all other nodes in the tree
        for node = 2:nT
            % get state deviation to nominal node state
            dx = x-s.Ts.x0(:,node);
            % evaluate its cost to go
            J = dx'*s.Ts.S(:,((node-1)*nX+1):(node*nX))*dx;
            % if the sample is within the local funnel description of the
            % node, and the cost-to-go is the lowest so far of all nodes
            % whose funnel contains the sample, update best starting node:
            if(J < s.Ts.epsH(node))
                if(J < jmin)
                    ind = node;
                    jmin = J;
                    flag = 1;
                end
            end
            % keep track of overall minimal cost to go if no node has been
            % found so far whose funnel contains the sample
            if(J < jminall)
                indMinAll = node;
                jminall = J;
            end
        end
        if(flag < 0)
            ind = indMinAll;
        end
    end

% Function to check state constraints on sampled state trajectory
% return true if inside, and false if outside
    function inside = testConstraints(xtape)
        if(isempty(s.options.xLimit))
            inside = true;
            return;
        end
        upperLimit = repmat(s.options.xLimit(:,2),1,size(xtape,2));
        lowerLimit = repmat(s.options.xLimit(:,1),1,size(xtape,2));
        inside = all(all(xtape<upperLimit & xtape>lowerLimit));
    end
% Function to check collision on sampled state trajectory
% return true if no collision, and false if collision
     function No_collid = testCollision(config_tape)
        num_obs = size(s.options.obstacles,1); % number of obstacles
        ridius = zeros(num_obs,size(config_tape,2));
        for p=1:size(config_tape,2)
            pos = config_tape(:,p);
           for q=1:num_obs
                    ridius(q,p)=sqrt((pos(1)-s.options.obstacles(q,1))^2+(pos(2)-s.options.obstacles(q,2))^2);
           end 
        end
        obstacle_radius = repmat(s.options.obstacles(:,3)+r,1,size(config_tape,2));
        No_collid = all(all(ridius>obstacle_radius));
     end
% Dynamics of the 2D Quadrotor
% xdot: state time derivative given time, state, and ctrl. input
% df:   gradients of dyn. with resp. to time, state, and ctrl. input
    function [xdot, df] = dynamicsgr(t,x,u)
        
        xdot=[x(4),x(5),x(6),-(u(1)+u(2))*sin(x(3))/m,((u(1)+u(2))*cos(x(3))/m)-g,...
            (u(1)-u(2))*r/I]';
        
        
        if (nargout>1)
            dfdx=[0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;0,0,-(u(1)+u(2))*cos(x(3))/m,...
                0,0,0;0,0,-(u(1)+u(2))*sin(x(3)),0,0,0;0,0,0,0,0,0];
            
            dfdu=[0,0;0,0;0,0;-sin(x(3))/m,-sin(x(3))/m;cos(x(3))/m,cos(x(3))/m;r/I,-r/I];
            
            dfdt = [ 0; 0; 0; 0; 0; 0];
            df{1} = [dfdt, dfdx, dfdu];
        end
    end
end