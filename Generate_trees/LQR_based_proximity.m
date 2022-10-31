function [id,x_near] = LQR_based_proximity(x_rand,tree,A,B,c,R,t_horizon,count)
% This procedure is a measure of ”closeness” between points in state space.
% it follow the heuristic proposed in "LQR-Based Heuristics for Rapidly
% Exploring State Space" proposed by Elena Glassman and Russ Tedrake
%
% The inputs:
%           x_rand          - The random sampled state
%           tree            - the RRT
%
%                    xdot=A(x-x_rand)+Bu+c
%
%           A               - The jacobian of the nonlinear system with
%                             respect to the states at x=x_rand and Uf=0.
%           B               - the The jacobian of the nonlinear system with
%                             respect to the inputs at x=x_rand and Uf=0.
%           c               - the result of the nonlinear dynamics at  
%                             x=x_rand and Uf=0.
%
%           R               - The input penalty matrix
%
%           t_horizon         - The final time horizon 
%
%
% The outputs:
%           id              - the id of the nearest state on the tree
%           x_near          - the coordinates of the state
% Author: Kouider CHADLI
% Email: chadli_kouider@hotmail.com

P_inv = get_cost(A,B,R,t_horizon);       % find P^(-1)(t) for each t 
                                         % from t=0 to t=t_horizon
    
    j = inf;
    for q =1:count
        d = get_d(A,c,tree(q).x-x_rand,t_horizon); 
        cost = t_horizon+.5*d'*P_inv(t_horizon)*d;
        if cost < j
            j = cost;           % The LQR-based proximity heuristic’s value for the distance
                                % from x_tree to x_rand
            id = q;             % the id on the tree
            x_near = tree(id).x; % the closest state
        end
        
    end
end

function [P_inv] = get_cost(A,B,R,t_f)
    P_0 = A*0;
    P0 = reshape(P_0, [numel(P_0) 1]); % convert the P_0 matrix to column
                                          % vector to be usable in ode45 
    
    % Solve the following matrix differential equation:
    % Pdot(t) = AP(t)+P(t)A'+B*inv(R)*B'   with P(t_start)=P_0=0
    
    [t, P] = ode45(@(t,P) Pdot(P, A,B,R), [0 t_f], P0);
    

    P = spline(t, P'); 

    P_inv = @(t) reshape_inv(ppval(P,t)); % decouple the results 
                                          % with respective time
                                          % and make it piecewise
     
end

function [P_dot] = Pdot(P,A,B,R)
    P = reshape(P,[sqrt(length(P)) sqrt(length(P))]); % column vector to matrix
    
    P_dot = A*P + P*A'+B*(R\B');
    
    P_dot = reshape(P_dot,[numel(P_dot) 1]); % matrix to column vector
    
end

function P_inverse = reshape_inv(P)
    % function that reshape and inverse the input argument
    Psquare = reshape(P, [sqrt(length(P)) sqrt(length(P))]);
    P_inverse = inv(Psquare);
end

function [d] = get_d(A,c,x_err,t_f)
    d = expm(A*t_f)*x_err;
    %+integral(@(t) exp(A*(t_max - t))*c,t_0,t_f);
   % d_int = d*0;
    N = 10; % Number of sub-intervals
    d_int=zeros(2,N+1);
    delta_t = t_f/N;
    for i = 1:N + 1
      t = (i-1)*delta_t;
      d_int(:,i) = expm(A*(t_f - t))*c;
    end
    int = delta_t*cumtrapz(d_int')';
    d = d+int(:,end);
        
end

% function [u] = open_loop(A,B,R,Pf_inv,d,tf)
% % this function generates an optimal open-loop control 
% % that depends only on initial and final states
%    u = @(t) -R\B'*expm(A'*(tf-t))*Pf_inv\d;
% end