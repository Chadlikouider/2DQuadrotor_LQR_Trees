% Plot a trajectory in a Sim-LQR Tree policy.
% Optionally also plot its funnel.
%
% Inputs:    (nx: state dim, N: trajectory length)
%   x0:      nominal state trajectory to plot, (nx*N)
%   Si:      TVLQR cost-to-go matrices for each state in x0 (nx*(nx*N))
%   epsH:    funnel parameters for each element in x0, (1*N)
%   options: plot options, see plotTree for explanation
%

function plotFunnels(x0,Si,epsH,options)
N = size(epsH,2);

if(options.figNum > 0)
    figure(options.figNum);
else
    figure;
end

%get dimension of problem:
nx = size(x0,1);
ndp = size(options.plotDim,2);

ztp = linspace(0,2*pi,options.nSamp);
angles = zeros(nx-1,options.nSamp^(nx-1));
%smaple angles on n-sphere circle:
for i = 1:nx-1
    angles(i,:) = repmat(reshape(repmat(ztp,options.nSamp^(i-1),1),1,[]),1,options.nSamp^(nx-1-i));
end
%sample hypersphere surface
x = ones(nx,length(angles));
for j = 1:nx-1
    x(j,:) = x(j,:).*cos(angles(j,:));
    x(j+1:nx,:) = x(j+1 : end,:).*sin(repmat(angles(j,:),nx-j,1));
end

%set up axes:
minD = options.plotRange(:,1);
maxD = options.plotRange(:,2);

options.plotDim = options.plotDim > 0;
if(options.plotBasins)
    for i = N:-1:1
        if(epsH(i)>0 && epsH(i) ~= Inf)
            %project hypersphere onto ellipsis of epsH = x'Sx
            z = sqrt(epsH(i))*(Si(:,((i-1)*nx+1):i*nx)^(-1/2))*x + repmat(x0(:,i),1,length(x));
            % plot it
            for j = 1:ndp
                if(ndp>1) subplot(options.subplot(1),options.subplot(2),j); end
                hold on
                zi = z(options.plotDim(:,j),:);
                cvInd = convhull(zi(1,:),zi(2,:));
                fill3(zi(1,cvInd),zi(2,cvInd),-0.0001*ones(1,length(cvInd)),options.basinColor,'EdgeColor','none')
            end
        end
    end
end
for j = 1:ndp
    if(ndp>1)
        subplot(options.subplot(1),options.subplot(2),j);
        hold on;
    end
    %set axis
    if(size(x0,2) > 1)
        mi = minD(options.plotDim(:,j));
        ma = maxD(options.plotDim(:,j));
        axis([mi(1),ma(1),mi(2),ma(2)]);
    end
    if(options.trajColor(1) >= 0 && N > 1)
        xi = x0(options.plotDim(:,j),:);
        plot(xi(1,:),xi(2,:),'Color',options.trajColor,'LineWidth',options.linewidth)
    end
end
end