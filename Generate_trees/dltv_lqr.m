%   DLTV LQR (Discrete-linear-time-varying LQR)
%   Returns the discrete-time time-varying TVLQR policy for the continuous-
%   time system linearized and discretized around a discrete-time
%   trajectory x(k),u(k).
%
%   [Ki,Si] = ltv_lqr(PLANTFUN,Q,R,dt,x0,u0,options)
%
%   Inputs:
%   PLANTFUN:   Function handle to continuous-time (CT) dynamics of system
%               to be stabilized. The function needs to have the following
%               I/O: [dx,dfdx] = PLANTFUN(t,x,u), where dx is the state
%               derivative given time t, state x, and input u; and dfdx is
%               a 1x1 cell array with the cell corresponding to the
%               gradients of the dynamics with respect to time, state, and
%               input: dfdx{1} = [df/dt, df/dx, df/du], where size(df/dt)
%               is [nx,1], size(df/dx) is [nx,nx], and size(df/du) is
%               [nx,nu], where nx is the state dim. and nu is the input dim
%   Q and R:    The CT penalty matrices for the TVLQR design.
%               They are discretized in this function using the sampling
%               time provided with the nominal trajectory
%   dt:         The (possibly non-homogenenous) array of sampling times
%               associated with the nominal trajectory. dt(i) is sampling
%               time for i-th element of state and input trajectory. May
%               also be a scalar, if sampling times are homogeneous.
%   x0,u0:      Nominal state and input trajectory of system described by
%               PLANTFUN. size(x0) = [nx,N], size(u0) = [nu,N], where N is
%               the length of the DT nominal trajectory.
%   Options:    Optional struct with the following elements:
%               -bInfHorizon:  if true, seek infinite horizon soulution
%                              where K(:,:,end) = infinite horizon
%                              compensator from dlqr
%               - Qf:          final cost, if unset when finite horizon
%                              solution is desired, Qf is set to Q

%   Outputs: Ki, Si, Ai, Bi are three dimensional arrays with
%   Ki(:,:,i): corresponding to the TVLQR compensator matrix for discrete
%              time k
%   Si(:,:,i): discrete time riccati matrix for cost-to-go matrix S obtained from the TVLQR design
%   Ai(:,:,i), Bi(:,:,i): discrete-time dynamics of system about trajectory
%
%


function [Ki,Si,Ai,Bi] = dltv_lqr(PLANTFUN,Q,R,dt,x0,u0,options)

if (nargin < 7)
    options = struct();
end

% get length of nominal trajectory:
N = size(x0,2);
% get state and input dims:
nX = size(x0,1);
nU = size(u0,1);

if(nargout > 2)
    Ai = zeros(nX,nX,N);
    Bi = zeros(nX,nU,N);
end


%check length of u0

if (length(u0)~=(N))
    error('length(u0) should equal length(x0)');
end
if (length(dt)~=(N) && length(dt) ~= 1)
    error('length(u0) should equal length(x0)'); 
end
if(length(dt) == 1)
    dt = dt*ones(1,N);
end

if (~isfield(options,'Qf'))
    options.Qf = []; 
end
if (~isfield(options,'Kf')) 
    options.Kf = zeros(nU,nX); 
end
if (~isfield(options,'dMethod'))
    options.dMethod = 'zoh'; 
end

Ki = zeros(nU,nX,N);
Si = zeros(nX,nX,N);

if (isempty(options.Qf))
    %get linearized plant dynamics:
    [xdot,df] = PLANTFUN(0,x0(:,end),u0(:,end));
    
    % make sure the tape ends at an equilibrium (stabilizability is checked
    % by DLQR).
    if (max(abs(xdot))>0.01)
        error('can''t do infinite horizon stabilization except at a equilibrium');
    end
    %discretize dynamics
    sysd = c2d(ss(df{1}(:,1+(1:nX)),df{1}(:,nX+1+(1:nU)),ones(1,nX),zeros(1,nU)),dt(end),options.dMethod);
    [Ki(:,:,end),Si(:,:,end)] = dlqr(sysd.A,sysd.B,dt(end)*Q,dt(end)*R);
    if(nargout > 2)
        Ai(:,:,end) = sysd.A;
        Bi(:,:,end) = sysd.B;
    end
else
    Si(:,:,end) = options.Qf;
    Ki(:,:,end) = options.Kf;
end

%go backwards in time and update Ki and Si
for n = N-1:-1:1
    %get continuous dynamics
    [~,df] = PLANTFUN(0,x0(:,n),u0(:,n));
    sysd = c2d(ss(df{1}(:,1+(1:nX)),df{1}(:,nX+1+(1:nU)),ones(1,nX),zeros(1,nU)),dt(n),options.dMethod);
    Ad = sysd.A; % discrete A
    Bd = sysd.B; % discrete B
    if(nargout > 2)
        Ai(:,:,n) = sysd.A;
        Bi(:,:,n) = sysd.B;
    end
    Qp = Si(:,:,n+1); % cost to go i+1 (p for plus)
    Ki(:,:,n) = (dt(n)*R + Bd'*Qp*Bd)\Bd'*Qp*Ad;
    Si(:,:,n) = dt(n)*Q + Ad'*Qp*Ad - Ad'*Qp*Bd*Ki(:,:,n);
end
end