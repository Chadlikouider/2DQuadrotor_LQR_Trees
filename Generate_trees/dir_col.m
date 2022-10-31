function [tf,xtape,utape,info] = dir_col(PLANTFUN,COSTFUN,FINALCOSTFUN,x0,utape0,tf0,options)
% Open-loop trajectory optimization using the direct collocation method via
% Optimtraj.  which is based on nonlinear programming that is constructed 
% from fmincon function
%
% Inputs are:
%    [x,df] = PLANTFUN(t,x,u)    - the plant dynamics, including gradients
%    [g,dg] = COSTFUN(t,x,u)     - the cost function, including gradients
%    [h,dh] = FINALCOSTFUN(t,x)  - the final time cost, including gradients
%    x0                          - the initial conditions (a constraint)
%    utape0                      - an initial guess for the open-loop control
%    
%
%    Tmin,Tmax,Xmax,Xmin,Umax,Umin - impose hard constraints
%     
%    xf                          - impose final value constraint
%    xtape0                      - initial guess for the xtape.  Default is
%                                  the result of simulating utape0 from x0,
%                                  or (if xf is set) a straight line from
%                                  x0 to xf.
%    maxDT                       - a limit on the step-size, currently this
%                                  limits the horizon time, but in future 
%                                  implementations, dircol will increase the 
%                                  length of the tape as necessary.
% Outputs are:
%     tf                         - a scalar corresponding to the final time 
%    xtape                       - a matrix that holds the history of
%                                  the states in the trajectory.
%    utape                       - a matrix holding the history of inputs
%    info
% Kouider CHADLI (chadli_kouider@hotmail.com)
if (nargin<7)
    options = struct(); 
end
[nU,nT] = size(utape0);           % dimension of input and # of time steps 
nX = length(x0);                  % dimension of the state space
dt = tf0/(nT+1);
tvlqr_R=[500 0;0 500];
%% implement options
if (~isfield(options,'Tmax')) 
    options.Tmax = inf;
end
if (~isfield(options,'Tmin')) 
    options.Tmin = 0;  
end
if (~isfield(options,'Xmax'))
    options.Xmax = inf(nX,1); 
end
if (~isfield(options,'Xmin')) 
    options.Xmin = repmat(-inf,nX,1); 
end
if (~isfield(options,'Umax'))
    options.Umax = inf(nU,1); 
end
if (~isfield(options,'Umin')) 
    options.Umin = repmat(-inf,nU,1);
end
if (~isfield(options,'xf')) 

    options.xf = []; 
end 
if (~isfield(options,'xtape0'))
    
 for i=1:nX, options.xtape0(i,:) = linspace(x0(i),options.xf(i),nT+1); end

end
if (~isfield(options,'maxDT'))
    options.maxDT = 0.1; 
end
if (~isfield(options,'transcription'))
             options.transcription = 'trapezoid'; %  <-- this is robust, but less accurate
%            options.transcription = 'rungeKutta';
%         options.transcription = 'trapGrad';   % Trapezoid with gradient
% options.transcription = 'hermiteSimpson'; approximated using a medium-order spline
%  options.transcription = 'hermiteSimpsonGrad';% hermite Simpson with gradient
%               options.transcription ='rungeKuttaGrad';
end
if (exist('optimtraj')==0) 
  error('dir_col needs the matlab folder for optimtraj to be in the matlab path  if you don''t have it, you can download it in https://www.mathworks.com/matlabcentral/fileexchange/54386-optimtraj-trajectory-optimization-library');
end
%% set up optimtraj inputs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         SET UP FUNCTION HANDLES                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.func.dynamics = @(t,x,u) DYNFUN(t,x,u,PLANTFUN);
problem.func.pathObj = @(t,x,u) OBJFUN(t,x,u,tvlqr_R,COSTFUN);
if (isfield(options,'obstacles'))
problem.func.pathCst = @(t,x,u)pathCst(x,options.obstacles);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              SET UP BOUNDS ON TIME, STATE AND CONTROLS                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low =options.Tmin;
problem.bounds.initialTime.upp =options.Tmin;
problem.bounds.finalTime.low =options.Tmin;
problem.bounds.finalTime.upp =options.Tmax;

problem.bounds.state.low =options.Xmin;
problem.bounds.state.upp =options.Xmax;

problem.bounds.initialState.low = x0;
problem.bounds.initialState.upp = x0;
problem.bounds.finalState.low = options.xf;
problem.bounds.finalState.upp = options.xf;

problem.bounds.control.low = options.Umin;
problem.bounds.control.upp = options.Umax;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

                       %%%% Method-independent options:
problem.options.nlpOpt = optimset(...
            'Display','off',...   % {'iter','final','off'}
            'TolFun',1e-5,...
            'MaxFunEvals',2e5,...
            'MaxIter',2e3);   %options for fmincon
%  problem.options(2).nlpOpt = optimset(...
%             'Display','off',...   % {'iter','final','off'}
%             'TolFun',1e-6,...
%             'MaxFunEvals',1e5);   %options for fmincon
switch options.transcription
    case 'trapezoid'
        % First iteration: get a more reasonable guess
        problem.options.verbose = 3; % How much to print out? 3=>debug
        problem.options.method = 'trapezoid'; % Select the transcription method
        problem.options.trapezoid.nGrid = nT;  %method-specific options
        
        
%         % Second iteration: refine guess to get precise soln
%         problem.options(2).verbose = 3; % How much to print out? 3=>debug
%         problem.options(2).method = 'trapezoid'; % Select the transcription method
%         problem.options(2).trapezoid.nGrid = nT;  %method-specific options
    case 'trapGrad'  %trapezoid with analytic gradients
        
        problem.options.method = 'trapezoid'; % Select the transcription method
        problem.options.trapezoid.nGrid = nT;  %method-specific options
        problem.options.nlpOpt.GradConstr = 'on';
        problem.options.nlpOpt.GradObj = 'on';
        problem.options.nlpOpt.DerivativeCheck = 'off';
        
%         problem.options(2).method = 'trapezoid'; % Select the transcription method
%         problem.options(2).trapezoid.nGrid = nT+1;  %method-specific options
%         problem.options(2).nlpOpt.GradConstr = 'on';
%         problem.options(2).nlpOpt.GradObj = 'on';
     case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = nT;  %method-specific options
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = nT;  %method-specific options
        
    case 'hermiteSimpsonGrad'  %hermite simpson with analytic gradients
        
        problem.options.method = 'hermiteSimpson'; % Select the transcription method
        problem.options.hermiteSimpson.nSegment = nT;  %method-specific options
        problem.options.nlpOpt.GradConstr = 'on';
        problem.options.nlpOpt.GradObj = 'on';
        problem.options.nlpOpt.DerivativeCheck = 'off';
        
%         problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
%         problem.options(2).hermiteSimpson.nSegment = nT;  %method-specific options
%         problem.options(2).nlpOpt.GradConstr = 'on';
%         problem.options(2).nlpOpt.GradObj = 'on';
    case 'rungeKutta'
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'medium';
        problem.options(1).rungeKutta.nSegment=5;
        problem.options(1).rungeKutta.nSubStep=2;
    case 'rungeKuttaGrad'
      
        problem.options.method = 'rungeKutta'; % Select the transcription method
        problem.options.defaultAccuracy = 'medium';
        problem.options.nlpOpt.GradConstr = 'on';
        problem.options.nlpOpt.GradObj = 'on';
        problem.options.rungeKutta.nSegment=nT;
        problem.options.nlpOpt.DerivativeCheck = 'off';
    otherwise
        error('Invalid method!');
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%           Create an initial guess for the trajectory                    %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.guess.time = 0:dt:(dt*nT);
problem.guess.state = options.xtape0;
problem.guess.control = [utape0,zeros(nU,1)];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       SOLVE with direct collocation                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = optimTraj(problem);

tf = soln.grid.time(end);        % value of the final time
% t = linspace(soln.grid.time(1), soln.grid.time(end), 2000);
% % x_t = soln.interp.state(t);
% x_s = zeros(nX,2000);
% for ix = 1:nX
%      x_s(ix,:) = spline(soln.grid.time,...
%                            soln.grid.state(ix,:), t);
% end
% u_t = soln.interp.control(t);
xtape = [x0, soln.grid.state];   % size(xtape)=[nX,1+Nt]
utape = soln.grid.control;       % size(xtape)=[nU,Nt]
% xtape = [x0, x_s];   % size(xtape)=[nX,1+Nt]
% utape = u_t;       % size(xtape)=[nU,Nt]
info = soln.info.exitFlag;            % 1 ,2 correspond to successful outcomes
                                      % -1,-2 not feasible solution of exceed
                                      % the number of iterations.
info
if(info == 1 || info == 2)
    info = 1;
end
end
function [j,dj] = OBJFUN(t,x,u,R,COSTFUN)
[nX,nT]=size(x);
m = 0.486;
g=9.81;
ug=[m*g/2;m*g/2];
j = zeros(1,nT);
for i=1:nT
      j(i)=1+(u(:,i)-ug)'*R*(u(:,i)-ug);%COSTFUN(t,x(:,i),u(:,i));%1+u(:,i)'*R*u(:,i);
end
% specific for cart pole model
if nargout == 2  % Analytic gradients
    nU=length(u(:,1));
    grad_dim=1+nX+nU;
    dj = zeros(grad_dim,nT); %1+nX+nU=[time,elements of the state,elements of the inputs]
    for k=1:nT
    [~, dg]=COSTFUN(t,x(:,k),u(:,k));
    for p=1:grad_dim
        dj(p,k)=dg{1}(p);
    end
    end
end
end
function [dx, grad] = DYNFUN(t,x,u,PLANTFUN)
[nX,nT]=size(x);
dx = PLANTFUN(t,x,u);

if nargout == 2  % Analytic gradients
    nU=length(u(:,1));
    grad_dim=1+nU+nX;     %1+nU+nX=[time,elements of the state,elements of the inputs]
    dxgrad = zeros(nX/2,grad_dim,nT); 
    ddxgrad = zeros(nX/2,grad_dim,nT); 
    for i=1:nT
    [~,df] = PLANTFUN(t,x(:,i),u(:,i));
     for j=1:nX/2
         for k=1:grad_dim
             dxgrad(j,k,i) = df{1}(j,k);
             ddxgrad(j,k,i)=df{1}(nX/2+j,k);
         end
     end
    end
    grad = cat(1,dxgrad,ddxgrad);
end
end
function [C, Ceq] = pathCst(x,obstacles)
q1=x(1,:);
q2=x(2,:);
C=[obstacles(1,3)-((q1-obstacles(1,1)).^2+(q2-obstacles(1,2)).^2).^(0.5);...
    obstacles(2,3)-((q1-obstacles(2,1)).^2+(q2-obstacles(2,2)).^2).^(0.5);...
    obstacles(3,3)-((q1-obstacles(3,1)).^2+(q2-obstacles(3,2)).^2).^(0.5)];
Ceq=[];
end
