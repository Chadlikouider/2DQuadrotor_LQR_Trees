function [T,options] = simLQRtree(PLANTFUN, COSTFUN, xRange, goal, options)

%  --SIMLQRTREE-----------------------------------------------------------
%   Generates a simulation-based LQR-Tree feedback policy for a nonlinear
%   dynamic system Based on RRT algorithm.
%
%   [T, options] = simLQRtree(PLANTFUN,COSTFUN,xRange,goal,options)
%   --OUTPUT---------------------------------------------------------------
%   T:      Struct with LQR-Tree feedback policy
%           Fields: (nx: dim state, nu: dim input)
%           .S: Discrete-time (DT) TVLQR Cost-To-Go Matrices
%               Size nx * (nx*|T|) (nx: state dim, |T|: no. tree nodes)
%           .K: DT TVLQR compensator matrices
%               Size (nu*|T|) * nx
%           .epsH: Local funnel description parameter (set epsilonHat in paper of Tedrake)
%               Size  1 * |T|
%           .x0: Nominal state of node
%               Size nx * |T|
%           .u0: Nominal input of node
%               Size nu * |T|
%           .dt: Sampling time of node
%               Size  1 * |T|
%           .endnode: Pointer to last element of nominal trajectory leading
%                     to the goal, starting at this node.
%               Size  1 * |T|
%           .tTrajID: Integer that is ID for what trajectory the node
%                     belongs to (used only for plotting and bookkeeping).
%               Size  1 * |T|
%
%   options:  Struct with algorithm options that may have been completed or
%             altered at initialization.
%
%   --INPUT----------------------------------------------------------------
%   PLANTFUN: Function handle to CONTINUOUS plant dynamics. Needs to return
%             two arguments:
%             [xdot,df] = PLANTFUN(t,x,u);
%             xdot must be dx/dt, Size nx * 1
%             df{1} must be [df/dt,df/dx,df/du] if xdot = f(t,x,u)
%             Size nx * (1 + nx + nu)
%   COSTFUN:  Function handle to user-defined instantaneous cost function
%             [g,dg] = COSTFUN(t,x,u,xf)
%             t: time, x: state, u: input, xf: goal state
%             g: cost
%             dg: gradient of cost function with respect to t, x, u
%   XRANGE:   Box-shaped design set to cover, described by [xmin,xmax]
%             Size: nx * 2
%   GOAL: Struct containing:
%         .x: Goal state, size nx*1
%         .u: Goal Input, size nu*1
%         .Stilqr: Goal state TILQR cost-to-go matrix, must be positive
%                  definite
%         .Ktilqr: Goal state TILQR compensator matrix
%         .dt: Sampling time of TILQR policy
%         .Qtilqr: Goal state TILQR design state penalty matrix
%         .Rtilqr: Goal state TILQR design input penalty matrix
%         .epsH: Parameter that determines goal region, a state q is inside
%               the goal region if (q-x)'*Stilqr*(q-x) < epsH
%               if no epsH is provided, the algorithm determines epsH using
%               the procedure outlined in the thesis
%         .basin: Struct containing:
%                 .termCond struct containing param. for term. heuristic
%                 for goal basin approximation, with elements:
%                 .alpha:  event probability parameter in (0,1)
%                 .pAlpha: Bernoulli process probability parameter in (0,1)
%   OPTIONS: Options, elaborated below
%
%   --OPTIONS--------------------------------------------------------------
%   options.ARG where ARG can be:
%   REQURIED FIELDS:
%   termCond: Struct that contains two fields for the paramteres of
%             the termintation heuristic:
%             .alpha:  event probability parameter (0,1)
%             .pAlpha: Bernoulli process probability parameter (0,1)
%   maxIter: Iteration limit of algorithm (max number of random samples)
%   tvlqr: TVLQR cost matrices
%         .Q: state penalty matrix (nx*nx)
%         .R: input penalty matrix (nu*nu)
%         .Qf: init cost matrix for Riccati iteration backwards
%              usually chosen equal goal.Stilqr
%   dirCopts: Options for direct-collocation motion planning(RRT based DC):
%           .resampleToDefaultDT: Enables resampling to default DT
%           .defaultDT: default dt for new trajectories
%           .Tmin, .Tmax: min/max final time
%           .Umin, .Umax: min/max input (nu,1)
%           .Xmin, .Xmax: min/max states (nx,1)
%           .nDircolTrajectory: default number of nodes in a traj.
%           .maxNutape: maximal length of trajectory obtained
%   simAllNodes: Enable testing of all nodes whose funnel hypotheses
%                contain the sample xSample. Slower, and results in more
%                trajectories, but the funnels are approximated better.
%
%   OPTIONAL FIELDS:
%   saveFolder: Folder for where to save snapshots of tree evolution and
%               mat file with policy. Leave empty if no saving desired.
%   uLimit: Input Constraints [umin, umax], size: nu * 2, elem = Inf
%   xLimit: State Constaints [xmin, xmax], size: nx*2
%   plotOpt: Struct with plot options for plotting tree, see plotTree
%            function for details
%   initTree: path to folder containing an existing tree that should be
%             used to initialize policy generation (use for "warm" start or
%             resuming)
%   nplotPDF: Limit number of JPG/PDFs saved for tree evolution snapshots
%   seedTrajs: Struct array with seed trajectories. Each element idx of the
%              struct array must contain:
%              .seedTrajs(idx).xdir: nominal trajectory (nx * no. el. traj)
%              .seedTrajs(idx).udir: nominal input (nu * no. elem. traj)
%              .seedTrajs(idx).tf: Time stamp nodes (1, no. elem. traj)
%              .xLimit: State constraints [xmin, xmax] (nx*2)
%              .nNodesInit: Number of nodes to init storage arrays with
%              .plotIncrement in {0,1}: incremental plotting off/on
%                                        helpful for larger trees
%   promptUserAfterGoalApprox: Enable user button press prompt after
%                              approximating goal set
%

% get state and input space size
nX = size(goal.x,1);
nU = size(goal.u,1);

% check provided options:
options = checkOptions(options);

% calculate success number M threshold for termination condition
nSuccTerm = ceil(log(options.termCond.alpha)...
                    /log(options.termCond.pAlpha));

%% Init Tree:
% init data arrays
tS = zeros(nX,nX*options.nNodesInit); % tvlqr cost to go matrices
tK = zeros(nU*options.nNodesInit,nX); % tvlqr comp. matrices
tEpsH = zeros(1,options.nNodesInit);  % funnel parameters
tX0 = zeros(nX,options.nNodesInit);   % nominal states
tU0 = zeros(nU,options.nNodesInit);   % nominal inputs
tDt = tEpsH;                          % sampling times
tTrajID = tEpsH;                      % trajectory ID array
tEndNode = tEpsH;                     % trajectory end-node indices
svtime = zeros(options.maxIter,5);    % array to save iteration stats
xNoSuc = [];                          % empty array to store unplann. samp.

% initialize arrays, iteration counter, and load initTree if necessary
initTree();

try % try to close the figure window where the tree will be plotted
    close(options.plotOpt.figNum)
catch
end

% initialize required file system folders, if saving is enabled
if(~isempty(options.saveFolder))
    [saveFile, fileno] = initFolders();
else
    saveFile = [];
end

%% build tree
% (re)set success counter
nSuccess = 0;

while(nSuccess < nSuccTerm && nIterations <= options.maxIter)
    nIterations = nIterations + 1; % iteration counter
    
    % draw random sample from uniform distribution on xRange
    xSample = getRandSample(xRange);
    tic % measure sample simulation time
    [init,flag] = testSample(); % check if current tree policy can get
    tSimTree = toc;             % the sample to the goal
    % possible flags:
    % -2: outside all funnels in current tree (=> policy unchanged)
    % -1: sample was in one or more hypotheses, but stabilization failed or
    %     collision detection
    %  1: successfully brought to goal
    %  2: successfully brought to goal in first attempt (=> policy unchanged)
    %  3: inside goal basin, Lyapunov test ok (=> policy unchanged)
 
    % if sample failed to reach goal, grow the tree using motion planning
    % algorithm RRT
    if(flag < 0)
        [xdir] = growTree(init);
        % if xdir is empty, motion planning failed
        if(isempty(xdir))
            % record unsuccessful xsample
            xNoSuc(:,end+1) = [xSample;0];
            % don't reset nSuccess if no trajectory could be found and no
            % funnel was adjusted. If a funnel was adjusted, reset nSucc.
            if(flag < -1)
                % nSuccess = nSuccess;
            else
                nSuccess = 0;
            end
        else
            % record successfully grown tree from sample
            xNoSuc(:,end+1) = [xSample;1];
            % reset success counter, since policy changed
            nSuccess = 0;
            % plot tree:
            if(options.plotIncrement)
                plotIncrementalTree(xdir);
            else
                plotFullTree();
            end
            % save tree
            if(~isempty(options.saveFolder))
                saveTree();
            end
        end
    elseif(flag == 1)
        % if the sample reached the goal after a funnel was adjusted,
        % reset success counter, since policy changed
        nSuccess = 0;
    else % implies 
        % If the sample was successful, increase success counter
        nSuccess = nSuccess + 1;
    end

    % keep stats on tree evolution:
    % #nodes timeToSimSample sampleTestResult #successStreak #trajInTree
    svtime(nIterations,:) = [tIdx tSimTree flag nSuccess tTrajID(tIdx)];
    
    % display progress
    disp([saveFile,': Iteration ', int2str(nIterations),...
        ', nSuccess = ',int2str(nSuccess) ,', numel(tree) = ',...
        int2str(tIdx),', time to sim tree: ',num2str(tSimTree),...
        's @ testFlag = ',int2str(flag),', nTraj = ',int2str(tTrajID(tIdx))])
end

% plot full tree, if it is not incrementally plotted
if(~options.plotIncrement)
    plotFullTree();
end

% save final tree
T = saveTree();
disp(['Termination condition M = ',...
    int2str(nSuccess),' reached.']);

%% NESTED FUNCTION DEFINITONS
    % Function for optimal-control-based motion planning.
    % The trajectory optimizer is initialized with a trajectory from the
    % failed simulations (initTraj), if available
    function [xdir] = growTree(initTraj)
         % Set options
        dirCopts = options.dirCopts;
        % set length of trajectory
        nUtape = options.dirCopts.nDircolTrajectory;
        dirCopts.xf = goal.x;   % set the goal state as state constraint
        if (isfield(options,'obstacles'))
            dirCopts.obstacles=options.obstacles;
        end
        % check if sample is not in state constraints:
        if(~testConstraints(xSample))
            % if the sample is not in the state constraints, cannot plan a
            % trajectory (this is not caught by the dir_col function)
            xdir = [];
            return
        end
        if(isempty(initTraj))
            % if there is no init. trajectory, use random inputs
            udir = randn(nU,nUtape);
            tf0 = nUtape*options.dirCopts.defaultDT;
        else
            % otherwise, use init trajectory
            tf0 = initTraj.tsim(end);
            
            % re-sample trajectory to default number of elements
            tnew = linspace(initTraj.tsim(1),initTraj.tsim(end),nUtape+1);
            
            % interpolate states:
            dirCopts.xtape0 = zeros(nX,nUtape+1);
            for ix = 1:nX
                dirCopts.xtape0(ix,:) = spline(initTraj.tsim,...
                                               initTraj.xtape(ix,:), tnew);
            end
            % interpolate inputs
            udir = zeros(nU,nUtape);
            for iu = 1:nU
                if nU == 1
                    udir(1:(length(tnew) - 1)) = spline(initTraj.tsim,...
                                     [initTraj.utape 0], tnew(1:end-1));
                else
                    udir(iu,1:(length(tnew) - 1)) = spline(initTraj.tsim,...
                                    [initTraj.utape(iu,:) 0], tnew(1:end-1));
                end

            end
        end
        
        % Dircol loop that makes sure that trajectory is within specs in
        % terms of maximal allowed sampling times and maximal number of
        % trajectory elements
        dirCol = 1;
        while(dirCol)
            % run dircol (see dir_col function file for definitions)
            [tf,xdir,udir,info] = dir_col(PLANTFUN, COSTFUN,...
                            @finalCostDummy,xSample,udir,tf0,dirCopts);
            % extract time steps and samping time achieved
            tf = linspace(0,tf,size(xdir,2));
            dtDircol = tf(end) - tf(end-1);
            % check if achieved samp. time and traj. length is acceptable
            if((dtDircol > options.dirCopts.maxDT...
                    || dtDircol < options.dirCopts.minDT)...
                    && info == 1 ... % check if success
                    && nUtape ~= options.dirCopts.maxNutape... 
                    && nUtape ~= 10)
                % resample trajectory:
                tf0 = tf(end);
                % enforce max length
                nUtape = min(options.dirCopts.maxNutape,...
                                ceil(tf(end)/options.dirCopts.defaultDT));
                nUtape = max(10,nUtape); % enforce minimum length
                % and resample the trajectory
                ts = timeseries([xdir;[udir goal.u]],tf);
                tnew = linspace(tf(1),tf(end),nUtape+1);
                tsnew = resample(ts,tnew,'linear');
                ndata = squeeze(tsnew.Data);
                dirCopts.xtape0 = ndata(1:nX,:);
                udir = ndata(nX+1:end,1:nUtape);
                disp(['Setting trajectory length to ',int2str(nUtape)])
            else
                dirCol = 0;
                % check if maximal dt is not violated
                if(dtDircol > options.dirCopts.maxDT ...
                        || tf(end) == options.dirCopts.Tmax)
                    info = -1;
                end
            end
        end
        % check for success of dircol or RRT planning
        if(info == 1)
            if(options.dirCopts.resampleToDefaultDT)
                % resample trajectory to default dt
                tf0 = tf(end); %get final time
                %get new element number
                nNew = floor(tf0/options.dirCopts.defaultDT);
                ts = timeseries([xdir;[udir tU0(:,1)]],tf);
                tnew = linspace(0,nNew*options.dirCopts.defaultDT,nNew+1);
                tsnew = resample(ts,tnew,'linear');
                ndata = squeeze(tsnew.Data);
                xdir = [ndata(1:nX,:), xdir(:,end)];
                udir = ndata(nX+1:end,:);
                % adapt last sampling time in trajectory to match total
                % trajectory time.
                if(tf0 == tnew(end))
                    dtDircol = [ones(1,nNew)*options.dirCopts.defaultDT, tDt(1)];
                else
                    dtDircol = [ones(1,nNew)*options.dirCopts.defaultDT, tf0-tnew(end), tDt(1)];
                end
            end
            
            %get discrete LTV regulator for new trajectory
            [Ki,Si] = dltv_lqr(PLANTFUN,options.tvlqr.Q,options.tvlqr.R,dtDircol,xdir,[udir,tU0(:,1)],options.tvlqr);
            %add all except final node
            addnodes(xdir(:,1:end-1),udir,Si(:,:,1:end-1),Ki(:,:,1:end-1),inf,dtDircol(1));
        else
            % motion planning failed
            xdir = [];
        end
    end

    % Generate random sample drawing from a uniform PDF over xRange
    function xSample = getRandSample(xRange)
        % uniformly sample from xRange (= the design set in the thesis)
        if(~isfield(options,'obstacles'))
        xSample = xRange(:,1) ...
                    + rand(size(xRange,1),1).*(xRange(:,2) - xRange(:,1));
        else
            collision_free=1; % 1 sample in obstacle space
                              % 0 sample in the free space
            num_obs = size(options.obstacles,1); % number of obstacles
            r = zeros(num_obs,1);
            while collision_free
                xSample = xRange(:,1) ...
                    + rand(size(xRange,1),1).*(xRange(:,2) - xRange(:,1));
                pos = xSample(1:2);  % just the x and z configuration for collision detection
                for i=1:num_obs
                    r(i)=sqrt((pos(1)-options.obstacles(i,1))^2+(pos(2)-options.obstacles(i,2))^2);
                end
                if all(r > options.obstacles(:,3)+0.25)
                    collision_free =0;
                end
            end
        end
    end

    % function to test a random sample xSample
    function [init,flag] = testSample()
        %check if sample is already in final basin:
        [inGoal] = testGoal(xSample);
        if(inGoal)
            init = [];   % return empty motion planning init trajectory
            flag = 3;    % if xSample in goal, return 3
            return
        end
        
        % initialize success state and test result flag
        success = -1;   % needed in case all nodes whose funnels contain
                        % the sample are tested (simAllNodes in options)
        flag = -1;
        
        % check if sample is within any funnel in the tree, and also get
        % closest node to xSample in tree if it's not inside
        [inside, closest] = inBasin(xSample, tX0, tS, tEpsH, tIdx);
        
        % if xSample is in any funnel, start simulation(s)
        if(inside)
            % get distances from sample to all tree node states:
            ji = getC2G(xSample, tX0, tS, tIdx);
            % ignore goal node
            ji = ji(2:end);
            % check distances to find what node funnels the sample is in
            drji = tEpsH(2:tIdx) - ji; % drji > 0 if inside
            % init node index array
            inds = 2:tIdx;
            % get indices of nodes where the sample is within funnel disc.
            lind = drji > 0;
            % get relevant arrays of distances and indices
            inds = inds(lind);
            ji = ji(lind);
            % sort by distances, then start simulating in ascending
            % order (lowest distance first).
            [~,ind] = sort(ji,'ascend');
            jbest = inf; % init variable to keep track best init traj.
            i = 1;  % init index to go through all possible node policies
            while(i <= length(ind)) % go through all nodes
                % if node is outside funnel due to earlier sim, continue
                if((ji(ind(i)) - tEpsH(inds(ind(i)))) >= 0)
                    i = i + 1;
                    continue;
                end
                
                %simulate using the policy of current Node policy
                [xtape,utape,tsim] = simTree(inds(ind(i)));
                
                %check if simulation is successful
                stateConstraintsOK = true;
                collisionFreeOK = true;
                [inGoal] = testGoal(xtape(:,end));
                %if the sample ended up in the goal,
                % check if state constraints were violated:
                if(inGoal)
                    stateConstraintsOK = testConstraints(xtape);
                    collisionFreeOK = testCollision(xtape(1:2,:)); % comment it if the the model isn't a quadrotor
                end
                
                if(inGoal && stateConstraintsOK && collisionFreeOK)
                    % only set flag the first time a successful sim is
                    % encountered (needed in case simAllNodes is set)
                    if(success < 0)
                        %we're good, set success
                        success = 1 + (i == 1);
                        % return 2 if the first node was successful
                        % and 1 if it was a lter node
                        if(~options.simAllNodes)
                            % return if not all nodes whose funnel contains
                            % the sample should be tested
                            flag = success;
                            init = [];
                            return
                        end
                    end
                else
                    % Failure, adjust epsH
                    adjustFunnels(xtape,inds(ind(i)));  % adjust funnels
                    
                    % store the trajectory for initializing motion-planning
                    cost = getCost(xtape,utape,tsim);   % get cost of traj.
                    % if the failed trajectory has the lowest cost
                    % store it
                    if(cost < jbest)
                        init.xtape = xtape;
                        init.utape = utape;
                        init.tsim = tsim;
                        jbest = cost;
                    end
                end
                i = i+1;
            end
        else
            % if the sample is not in any funnel, obtain an initialization
            % trajectory using the closest node's feedback policy
            flag = -2;
            % simulate from closest node to get initial condition for
            % dircol, only if tree has more than one node (=goal), and
            % closest node is not goal node
            if(tIdx > 1 && closest > 1)
                [xsim,usim,tsim] = simTree(closest);
                init.xtape = xsim;
                init.utape = usim;
                init.tsim = tsim;
            else
                init = [];   % return empty motion planning init trajectory
            end
        end
        % in case that simAllNodes is set, make sure flag is set before
        % function returns
        if(success > 0)
            flag = success;
            init = [];
        end
    end

    % Evaluate cost for determining what failed simulation trajectory
    % should be used for initializing motion planning.
    function cost = getCost(xtape,utape,tsim)
        dt = tsim(2) - tsim(1);
        xdes = repmat(goal.x,1,size(xtape,2)-1);
        udes = repmat(goal.u,1,size(utape,2));
        xerr = xtape(:,1:(end-1))-xdes;
        uerr = utape(:,1:end)-udes;
        % get instantaneous cost
        if size(utape,1)==1
         cost = sum(sum((options.tvlqr.Q*xerr*dt).*xerr,1)...
           + (options.tvlqr.R*uerr*dt).*uerr)... % and final cost
           + (xtape(:, end) - goal.x)'*tS(:, 1:nX)*(xtape(:,end) - goal.x);
        else
         cost = sum(sum((options.tvlqr.Q*xerr*dt).*xerr,1)...
           + sum((options.tvlqr.R*uerr*dt).*uerr,1))... % and final cost
           + (xtape(:, end) - goal.x)'*tS(:, 1:nX)*(xtape(:,end) - goal.x);
        end
    end

    % Function that checks box-shaped state contraints defined by xLimit
    function inside = testConstraints(xtape)
        if(isempty(options.xLimit))
            inside = true;
            return;
        end
        upperLimit = repmat(options.xLimit(:,2),1,size(xtape,2));
        lowerLimit = repmat(options.xLimit(:,1),1,size(xtape,2));
        inside = all(all(xtape<upperLimit & xtape>lowerLimit));
    end
    % Function that checks for collision (comment the function if not a quadrotor model)
    function No_collid = testCollision(config_tape)
        if(~isfield(options,'obstacles'))
           No_collid = true;
           return;
        end
        num_obs = size(options.obstacles,1); % number of obstacles
        r = zeros(num_obs,size(config_tape,2));
        for i=1:size(config_tape,2)
            pos = config_tape(:,i);
           for j=1:num_obs
                    r(j,i)=sqrt((pos(1)-options.obstacles(j,1))^2+(pos(2)-options.obstacles(j,2))^2);
           end 
        end
        obstacle_radius = repmat(options.obstacles(:,3)+0.25,1,size(config_tape,2));
        No_collid = all(all(r>obstacle_radius));
    end
    % Simulates the random sample with the policy of node ind
    function [xtape,utape,tsim,x0] = simTree(ind)
        % get open-loop trajectory and TV LQR compensator matrices
        x0 = tX0(:,ind:tEndNode(ind));
        u0 = tU0(:,ind:tEndNode(ind));
        Ksim = tK(((ind-1)*nU+1):(tEndNode(ind)*nU),:);
        tsim = tDt(1,ind:tEndNode(ind));
        % init simulation time:
        tsim = [0,cumsum(tsim)];
        % Simulate
        if(options.ode45MaxStep)
            % enforce maximal time step (may improve sim. accuracy)
            opts = odeset('MaxStep',options.dirCopts.defaultDT);
            % closed-loop dynamics are in fofx
            [~,xtape] = ode45(@fofx, tsim, xSample, opts);
        else
            [~,xtape] = ode45(@fofx, tsim, xSample);
        end
        
       if(numel(tsim) == 2) % catch degenerate cases of 2 step sims
            xtape = [xtape(1,:); xtape(end,:)];
        end
        
        xtape = xtape';
        %recover input tape (for init trajectories):
        utape = zeros(nU,size(xtape,2)-1);
        for i = 1:size(utape,2)
             u = u0(:,i) - Ksim(((i-1)*nU + 1):i*nU,:)*(xtape(:,i) - x0(:,i));
%               u = u0(:,i) - tK(1:nU,:)*(xtape(:,i) - x0(:,i));
            if(~isempty(options.uLimit))
                u = min([u,options.uLimit(:,2)],[],2);
                u = max([u,options.uLimit(:,1)],[],2);
            end
            utape(:,i) = u;
        end
        
        % Nested function that evaluates closed-loop dynamics
        % again this is not the clean way to simulate the CT system, since
        % the control input may change within a time-step.
        
        function [xdot] = fofx(t,x)
            tind = find(tsim <= t,1,'last');
            if(tind > length(u0(1,:)))
                % apply tilqr policy if tvlqr policy ended
                uDT = tU0(:,1) - tK(1:nU, :)*(x - tX0(:,1));
            else
                dx = (x - x0(:,tind));
                uDT = u0(:,tind) - Ksim(((tind-1)*nU + 1):tind*nU,:)*dx;
%                   uDT = u0(:,tind) - tK(1:nU,:)*dx;
            end
            
            if(~isempty(options.uLimit))
                uDT = min([uDT,options.uLimit(:,2)],[],2);
                uDT = max([uDT,options.uLimit(:,1)],[],2);
            end
            % get derivative from system dynamics function:s
            xdot = PLANTFUN(t,x,uDT);
        end
    end

    % Function that adjusts the funnels after a failed simulation
    function [] = adjustFunnels(xtape,ind)
        nel = tEndNode(ind) - ind + 1;
        % get new values for the epsH values from the mex function:
        newEps = getNewRho(xtape(:,1:nel), tX0(:,ind:tEndNode(ind)),...
                           tS(:,((ind-1)*nX+1):tEndNode(ind)*nX),...
                           tEpsH(1,ind:tEndNode(ind)));
        % read out current values of the epsH:
        tempRho = tEpsH(1,ind:tEndNode(ind));
        % the negative entries of newEps indicate funnel parameters that
        % should be left unchanged:
        newEps(newEps < 0) = tempRho(newEps < 0);
        tEpsH(1,ind:tEndNode(ind)) = newEps;
    end

    % Function to add new nodes to the tree after a successful motion plan
    function [] = addnodes(xtape, utape, Si, Ki, epsHinit, dt)
        nel = size(xtape,2);
        % make sure the trajectory is not just a single point:
        assert(nel >= 2,...
           'tried to add trajectory with less than 2 elements in addnodes')

        % add trajectory to tree
        tS(:,(tIdx*nX + 1):(tIdx+nel)*nX) = reshape(Si,nX,[]);
        tK((tIdx*nU + 1):(tIdx+nel)*nU,:) = reshape(Ki,nX,[])';
        tEpsH(tIdx + 1 : tIdx + nel) = epsHinit;
        tX0(:,(tIdx + 1 : tIdx + nel)) = xtape;
        tU0(:,(tIdx + 1 : tIdx + nel)) = utape;
        tDt(tIdx + 1 : tIdx + nel) = dt;
        tTrajID(tIdx + 1 : tIdx + nel) =  tTrajID(tIdx) + 1;
        tEndNode(tIdx + 1 : tIdx + nel) = tIdx + nel;
        tIdx = tIdx + nel;
    end

    % function that checks if a state x is inside goal set G
    function [inside] = testGoal(x)
        dx = x - goal.x;
        % cost of x
        ji = dx'*tS(:,1:nX)*dx;
        % if outside, return false
        inside = ji < tEpsH(1);
    end

    % function that determines goal set G by approximating the invariant
    % set of the closed-loop goal-state dynamics.
    function epsH = findGoalBasinSim()
        disp('Approximating goal basin of attraction...')
        tgoal = tic;
        epsH = Inf;
        % determine success limit
        termCond = ceil(log(goal.basin.termCond.alpha)...
            /log(goal.basin.termCond.pAlpha));
        % determine eigenvalues and eigenvectors of goal.Stilqr
        [L,ev] = eig(goal.Stilqr);
        ev = diag(ev);
        % init state that controls process of first estimating epsH and
        % then estimating the success rate, given epsH
        estimate = 2;
        
        while(estimate > 0)
            % if estimate = 2, get initial estimate epsH
            if(estimate > 1)
                maxIter = options.maxIter;
            else
                % if estimate = 1, estimate rate of samples for which
                % Lyapunov condition holds. Use nEvals samples.
                maxIter = goal.basin.nEvals;
                % make sure term cond does not kick in
                termCond = goal.basin.nEvals; 
                disp(['Found epsH = ',num2str(epsH),' in ',...
                    num2str(toc(tgoal)),'s.']);
                disp(['Now evaluating success rate with ',...
                    int2str(maxIter),' samples']);
            end
            nGoalSuccess = 0;
            iter = 0;
            while(nGoalSuccess < termCond && iter < maxIter)
                iter = iter + 1;
                % get samples until a sample is inside the ellipse
                inside = 0;
                while(~inside)
                    if(isinf(epsH))
                        x = getRandSample(xRange);
                        dx = (x-goal.x);
                        ji = dx'*goal.Stilqr*dx;
                        inside = ji < epsH;
                    else
                        lim = sqrt(epsH./ev);
                        X = -lim + rand(nX,1).*2.*lim;
                        inside = sum(X.^2.*ev) < epsH;
                        if(inside)
                            % calculate x
                            x = L*X + goal.x;
                            % calculate dx and check if x is inside of ellipse
                            % with radius epsH
                            dx = x - goal.x;
                            ji = dx'*goal.Stilqr*dx;
                            inside = ji < epsH;
                        end
                    end
                end

                if(testConstraints(x))
                    % get control input for a single step
                    u = goal.u - tK(1:nU,:)*dx;
                    % apply input constraints
                    if(~isempty(options.uLimit))
                        u = min([u,options.uLimit(:,2)],[],2);
                        u = max([u,options.uLimit(:,1)],[],2);
                    end
                    % simulate system
                    [~,xtape] = ode45(@(t,x)PLANTFUN(t,x,u),[0,tDt(1)],x);
                    x = xtape(end,:)';
                    % evaluate Lyapunov function    
                    if((x-goal.x)'*goal.Stilqr*(x-goal.x) - ji < 0)
                        % decreasing, increase success count
                        nGoalSuccess = nGoalSuccess + 1;
                    else
                        % increasing, shrink epsH; do this also while
                        % estimating success rate, since it's good to have
                        % a conservative estimate
                        epsH = ji;
                        disp(['Lyapunov increasing by ',num2str((x-goal.x)'*goal.Stilqr*(x-goal.x) - ji),', shrunk goal basin epsH to ',num2str(epsH)])
                        % only reset success counter during initial estim.
                        nGoalSuccess = (estimate < 2) * nGoalSuccess;
                    end
                else
                    % if state constraints are violated:
                    epsH = ji;
                    disp(['State constraints violated, shrunk goal basin epsH to ',num2str(epsH)])
                    % only reset success counter during initial estim.
                    nGoalSuccess = (estimate < 2) * nGoalSuccess;
                end
            end
            % go from initial estimate to success rate estimate to done
            estimate = estimate - 1;
        end
        % make sure the goal is not almost just a point:
        assert( epsH > eps,...
            'Goal ellipse shrunk to zero, double check parameters');
        disp(['Approximated a goal set with epsH = ',...
            num2str(epsH)])
        % calculate success rate to display:
        [p,pci] = binofit(nGoalSuccess, maxIter);
        disp(['With ', int2str(maxIter),' samples from Goal basin, ',...
            'the success rate is at least ',num2str(p),...
            ' with confidence int. ', num2str(pci)]);
        if(options.promptUserAfterGoalApprox)
            input('Press the Return key to start generating the tree policy');        
        end
    end

    % The motion-planning function we use requires a final cost function.
    % Since the final state is introduced as a constraint, this function is
    % not used and we pass this dummy function.
    function [h,dh] = finalCostDummy(t,x)
        h = 0;
        % dh = [dcost/dt, dcost/dx] => the gradient with respect to time
        %                              and state, respectively 
        dh{1} = zeros(1, 1 + nX);
    end

    % Initialize the tree data structures depending on whether a tree
    % generation should be resumed or a fresh start is given:
    function [] = initTree()
        % check if we should resume with building a tree:
        if(~isempty(options.initTree))
            % if yes, extract the data for resuming:
            s = load([options.initTree,'/mat/',options.initTree,'.mat']);
            tIdx = length(s.Ts.epsH);
            xNoSuc = s.xNoSuc;
            nIterations = size(s.svtimeSave,1);
            svtime(1:nIterations,:) = s.svtimeSave;
            tS(:,1:(nX*tIdx)) = s.Ts.S;
            tK(1:(nU*tIdx),:) = s.Ts.K;
            tEpsH(:,1:tIdx) = s.Ts.epsH;
            tX0(:,1:tIdx) = s.Ts.x0;
            tU0(:,1:tIdx) = s.Ts.u0;
            tDt(:,1:tIdx) = s.Ts.dt;
            tTrajID(:,1:tIdx) = s.Ts.trajInd;
            tEndNode(:,1:tIdx) = s.Ts.endnode;
        else
            % if no tree is provided, start fresh:
            tS(:,1:nX) = goal.Stilqr;
            tK(1:nU,:) = goal.Ktilqr;
            tEpsH(:,1) = Inf; % init goal basin hypothesis to cover D
            tX0(:,1) = goal.x;
            tU0(:,1) = goal.u;
            tDt(:,1) = goal.dt;
            tEndNode(:,1) = 1;
            tIdx = 1;
            
            % init total number of samples counter
            nIterations = 0;
            
            if(isfield(goal.basin,'epsH'))
                tEpsH(1) = goal.basin.epsH;
            else
                % determine epsH of goal
                tEpsH(1) = findGoalBasinSim();
                goal.basin.epsH = tEpsH(1);
            end
            
            % extract and add seed trajectories:
            if(~isempty(options.seedTrajs))
                for idx = 1:length(options.seedTrajs)
                    xdir = options.seedTrajs(idx).xdir;
                    udir = options.seedTrajs(idx).udir;
                    tf = options.seedTrajs(idx).tf;
                    dtDircol = tf(2) - tf(1);
                    
                    if(options.dirCopts.resampleToDefaultDT)
                        % resample trajectory to default dt
                        tf0 = options.seedTrajs(idx).tf(end); %get duration
                        %get new element number
                        nNew = floor(tf0/options.dirCopts.defaultDT);
                        ts = timeseries([xdir;[udir tU0(:,1)]], tf);
                        tnew = linspace(0,...
                            nNew*options.dirCopts.defaultDT, nNew+1);
                        tsnew = resample(ts,tnew,'linear');
                        ndata = squeeze(tsnew.Data);
                        xdir = [ndata(1:nX,:),xdir(:,end)];
                        udir = [ndata(nX+1:end,:)];
                        if(tf0 == tnew(end))
                            dtDircol = ...
                         [ones(1,nNew)*options.dirCopts.defaultDT, tDt(1)];
                        else
                            dtDircol = ...
                               [ones(1,nNew)*options.dirCopts.defaultDT,...
                                tf0-tnew(end), tDt(1)];
                        end
                    end
                    
                    % generate tvlqr policy for trajectory
                    [Ki,Si] = dltv_lqr( PLANTFUN,...
                        options.tvlqr.Q,...
                        options.tvlqr.R,...
                        dtDircol,...
                        xdir,...
                        [udir,tU0(:,1)],...
                        options.tvlqr);
                    
                    %add all except final node (equal to goal state)
                    if(     options.dirCopts.resampleToDefaultDT ...
                            && dtDircol(1) ~= options.dirCopts.defaultDT)
                        error('stop, resampling failed');
                    end
                    % add nodes to tree:
                    addnodes( xdir(:,1:end-1),...
                        udir,Si(:,:,1:end-1),...
                        Ki(:,:,1:end-1),...
                        Inf,... % set ball radii to Inf
                        dtDircol(1));
                    % plot stuff, if necessary
                    if(options.plotIncrement)
                        % do not save picture: 2nd arg = 0
                        plotIncrementalTree(xdir, 0);
                    end
                end
            end
        end
    end

    % Create folders for saving the policy
    function [saveFile, fileno] = initFolders()
        % create saving directories, if necessary
        if(~isempty(options.saveFolder))
            % if an tree is resumed, create a new folder for it, and append
            % folder name accordingly.
            if(~isempty(options.initTree))
                saveFile = [options.initTree,...
                    '_resumedAt_',datestr(now,30)];
                options.saveFolder = saveFile;
            else
                saveFile = options.saveFolder;
            end
            %setup folders and subfolders
            if(~exist(options.saveFolder,'file'))
                mkdir(options.saveFolder)
                fileno = 1;
            else
                dirlist = dir(options.saveFolder);
                a = [dirlist.isdir];
                fileno = length(a) - sum(a) + 1;
            end
            if(~exist(['./',options.saveFolder,'/pdf'],'file'))
                mkdir(['./',options.saveFolder,'/pdf'])
            end
            if(~exist(['./',options.saveFolder,'/mat'],'file'))
                mkdir(['./',options.saveFolder,'/mat'])
            end
        else
            saveFile = 'notSavingData';
        end
    end
    
    % Incremental plotting of the tree, trajectory by trajectory
    function [] = plotIncrementalTree(xtraj, sv)
        if(nargin < 2)
            % by default save if enabled
            sv = 1;
        end
        figure(options.plotOpt.figNum);
        for dim = 1:size(options.plotOpt.plotDim, 2)
            if(~isempty(options.plotOpt.subplot))
                subplot(options.plotOpt.subplot(1),...
                        options.plotOpt.subplot(2), dim);
            end
            if(~isempty(options.plotOpt.titles))
                title(options.plotOpt.titles{dim},...
                    'FontSize',options.plotOpt.fsize);
            end
            if(~isempty(options.plotOpt.ylabels))
                ylabel(options.plotOpt.ylabels{dim},...
                    'FontSize', options.plotOpt.fsize);
            end
            if(~isempty(options.plotOpt.xlabels))
                xlabel(options.plotOpt.xlabels{dim},...
                    'FontSize', options.plotOpt.fsize);
            end 
            xi = xtraj(options.plotOpt.plotDim(:,dim)>0,:);
            plot(xi(1,:),xi(2,:), 'Color', options.plotOpt.trajColor,...
                'LineWidth', options.plotOpt.linewidth), hold on
            xm = options.plotOpt.plotRange(...
                                options.plotOpt.plotDim(:,dim) > 0, 1);
            xp = options.plotOpt.plotRange(...
                                options.plotOpt.plotDim(:,dim) > 0, 2);
            axis([xm(1) xp(1) xm(2) xp(2)])
        end
        if(~isempty(options.saveFolder) && sv && fileno < options.nplotPDF)
            fileOrig= sprintf(['./',options.saveFolder,'/pdf/%08i.pdf'],...
                                fileno);
            print('-dpdf', '-painters',fileOrig); %creates 1 frame
            fileno = fileno+1;
        end
    end

    % Plots the full tree, including funnels
    function plotFullTree()
        try
            close(options.plotOpt.figNum)
        catch
        end
        Ts.x0 = tX0(:,1:tIdx);
        Ts.u0 = tU0(:,1:tIdx);
        Ts.dt = tDt(:,1:tIdx);
        Ts.S = tS(:,1:(tIdx*nX));
        Ts.K = tK(1:(tIdx*nU),:);
        Ts.epsH = tEpsH(:,1:tIdx);
        Ts.endnode = tEndNode(:,1:tIdx);
        plotTree(Ts,options.plotOpt);
        if(~isempty(options.saveFolder) && fileno < options.nplotPDF)
            fileOrig = sprintf(...
                ['./',options.saveFolder,'/pdf/%08i.pdf'],fileno);
            print('-dpdf', '-painters',fileOrig); %creates 1 frame
            fileno = fileno+1;
        end
    end
    
    % function that checks the options struct upon startup
    function opts = checkOptions(opts)
        if (~isfield(opts,'uLimit'))
            opts.uLimit = [];
        end
        if (~isfield(opts,'xLimit'))
            opts.xLimit = [];
        end
        if (~isfield(opts,'plotOpt'))
            opts.plotOpt = [];
        end
        if (~isfield(opts,'saveFolder'))
            opts.saveFolder = [];
        end
        if (~isfield(opts,'maxIter'))
            error('You need to provide a maximal iterations limit');
        end
        if (~isfield(opts,'initTree'))
            opts.initTree = [];
        end
        if (~isfield(opts,'nplotPDF'))
            opts.nplotPDF = [];
        end
        if (~isfield(opts,'seedTrajs'))
            opts.seedTrajs = [];
        end
        if (~isfield(opts,'nNodesInit'))
            opts.nNodesInit = 2e4;
        end
        if (~isfield(opts,'plotIncrement'))
            opts.plotIncrement = false;
        end
        if(~isfield(opts,'simAllNodes'))
           opts.simAllNodes = 0;
        end
        if (~isfield(opts,'tvlqr'))
            error('You must provide tvlqr parameters');
        end
        if (~isfield(opts,'dirCopts'))
            error('You need to provide direct collocation options');
        end
        if (~isfield(opts,'ode45MaxStep'))
            opts.ode45MaxStep = 0;
        end
        if (~isfield(opts,'promptUserAfterGoalApprox'))
            opts.promptUserAfterGoalApprox = 1;
        end
        if (~isfield(opts,'termCond'))
            error('You must provide parameters for the term. heuristic')
        end
        if(~isfield(goal,'basin'))
            error('You must provide goal basin parameters')
        end
        assert( opts.dirCopts.minDT < opts.dirCopts.defaultDT...
            && opts.dirCopts.maxDT > opts.dirCopts.defaultDT,...
            'Default DT for Dircol must be in min/max time interval')
    end

    % Function that saves tree to mat file in save folder
    function Ts = saveTree()
        Ts.x0 = tX0(:,1:tIdx);
        Ts.u0 = tU0(:,1:tIdx);
        Ts.dt = tDt(:,1:tIdx);
        Ts.S = tS(:,1:(tIdx*nX));
        Ts.K = tK(1:(tIdx*nU),:);
        Ts.trajInd = tTrajID(:,1:tIdx);
        Ts.epsH = tEpsH(:,1:tIdx);
        Ts.endnode = tEndNode(:,1:tIdx);
        if(~isempty(options.saveFolder))
            svtimeSave = svtime(1:nIterations,:);
            save(['./',options.saveFolder,'/mat/',saveFile,'.mat'],...
                'Ts','options','xNoSuc','xRange','svtimeSave','goal');
        end
    end
end