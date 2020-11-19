function meanErr = RE_ProblemSet5Problem5_Solution()
% meanErr = RE_ProblemSet5Problem5_Solution()
%
% Example Solution to Problem 5 of the Particle Filtering Problem Set.
%
% Outputs a performance measure of the Particle Filter: The mean of the
% average Euclidean distance error of the 50% closest particles to the
% person. You may use this to fine-tune your implementation.
%
% Class:
% Recursive Estimation
% Spring 2018
% Problem Set 5, Particle Filtering
%
% --
% Revision history
% [19.05.11, ST]    First version by Sebastian Trimpe for Recitation
% [28.05.13, PR]    Adapted for Particle Filter Problem Set

%% PF Setup
% Select Number of Particles
N = 300;
% Select Number of Simulation Steps:
Nsim = 90;

%% Animation Setup
% Delay in between time steps
pauseSec = 0.01;
% Activate step-by-step animation, press any key on the keyboard to proceed
% by one time step. Abort by pressing any mouse button.
activateButtonPress = 0;

%% System Parameters
% Process model
theta = 1/12;       % Horizontal velocity of person, in km (per time step)
Q = (theta/4)^2*eye(2); % Process noise variance in kilometers^2
s0 = [7.5; 7.5];        % Mean of initial state distribution in kilometers
P0 = 1/25*eye(2); % Variance of initial state distribution in kilometers^2

% Measurement model
R = (0.05)^2;   % Variance of measurement noise in kilometers^2

%% Initialize Storage Arrays for Person State, Measurements, and Particles

% Person State:
personStates = zeros(2,Nsim + 1); % [s(0), s(1), s(2),...,s(Nsim)]
% personStates(1,n): person location x at time k = n - 1
% personStates(2,n): person location y at time k = n - 1

% Measurements:
measurements = zeros(1,Nsim + 1);
% measurements(n): altitude measurement z at time k = n - 1

% Posterior particles:
posteriorParticles = zeros(2,N,Nsim + 1);
% posteriorParticles(1,:,n): posterior particles of location x at k = n - 1
% posteriorParticles(2,:,n): posterior particles of location y at k = n - 1

% If error measure requested, initialize error array:
if(nargout > 0)
    errs = zeros(1,Nsim+1);
end

%% CHANGE CODE HERE: Simulation of person, Part a)
% Run this function after you implemented Part a) to check the results.

% Draw random initial state s(0) here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
personStates(:,1) = s0 + chol(P0)*randn(2,1);

for k = 2:(Nsim+1) % for k = 1:Nsim (2:Nsim+1 because Matlab is 1-based indexing)
    % Implement the update of the person's state here: %%%%%%%%%%%%%%%%%%%%
    % Get current altitude and terrain gradient
    [~,grad] = altitudeAndGradient(personStates(1,k-1),personStates(2,k-1));
    % Update person location:
    processNoiseSample = chol(Q)*randn(2,1); % v(k-1), in km
    % s(k) = s(k-1) - theta * grad / |grad| + v(k-1), in km
    personStates(:,k) = personStates(:,k-1) - theta*grad/norm(grad) + processNoiseSample;
    % Generate a measurement here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    measurementNoiseSample = sqrt(R)*randn; % meas. noise sample w(k), in km
    % Get current altitude of person:
    altitude = altitudeAndGradient(personStates(1,k), personStates(2,k));
    measurements(k) =  altitude + measurementNoiseSample; % = z(k), in km
end

%% CHANGE CODE HERE: Particle Filter, Parts b) - d)
% You can run and test your code after implementing each part.

% If you want to check your solution to parts b) and c), indicate this
% here, for the solution with roughening, set part = 'd'.
part = 'd';

% Implement particle initialization here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
posteriorParticles(:,:,1) = repmat(s0,1,N) + chol(P0)*randn(2,N);

for k = 2:(Nsim + 1)
    % b) Implement process update here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get gradients at particle positions:
    [~,grads] = altitudeAndGradient(posteriorParticles(1,:,k-1),posteriorParticles(2,:,k-1));
    % get gradient norms:
    normGrads = sqrt(sum(grads.^2));
    % draw noise samples:
    processNoiseSamples = chol(Q)*randn(2,N);
    % calculate prior particles:
    priorParticles = posteriorParticles(:,:,k-1) - theta * grads./[normGrads;normGrads] + processNoiseSamples;
    
    if(part > 'b')
        % c) Implement measurement update here: %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get heights of particles:
        [hparticles] = altitudeAndGradient(priorParticles(1,:),priorParticles(2,:));
        % get measurement likelihood of particles, beta_i = f_w(z - h_i) with
        % h_i the height of particle i, and z the given measurement
        betas = exp(-(measurements(k) - hparticles).^2/(2*R));
        % the normalization factor 1/sqrt(2*pi*R) is not neccessary, since the
        % betas will be normalized later anyway. Therefore, the above yields
        % the same result, after normalization, as using the full definition of
        % the Gaussian distribution:
        % betas = 1/sqrt(2*pi*R)*exp(-(measurements(k) - hparticles).^2/(2*R))
        if(sum(betas) > 0)
            % normalize weights
            betas = betas/sum(betas);
        else
            % In theory, this should never happen since the measurement noise
            % is Gaussian and nonzero for all z - h(x,y). In practice however,
            % if the particles are really far away from the actual person, the
            % particle weights can be set to zero due to numerical round-off.
            % In this case, we just assign equal probabilities to all
            % particles, which is not a good solution to this problem (i.e. use
            % something smarter in the programming exercise).
            betas = ones(1,N)*1/N;
            warning('Particle weights were all zero')
        end
        % build cumulative sum of particles (similar to a CDF)
        betaCumSum = cumsum(betas);
        % Resample according to lecture notes:
        for i = 1:N
            randNumber = rand;
            ind = find(betaCumSum >= randNumber,1,'first');
            posteriorParticles(:,i,k) = priorParticles(:,ind);
        end
        
        if(part > 'c')
            % d) Implement roughening here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K = 0.5; % Roughening parameter, tuned for the specific number of particles used
            D = 2;    % Dimension of the state space
            % Find maximal inter-sample variability
            Ei = [max(posteriorParticles(1,:,k)) - min(posteriorParticles(1,:,k));
                max(posteriorParticles(2,:,k)) - min(posteriorParticles(2,:,k))];
            % Build diagonal matrix of standard deviations for drawing roughening samples:
            StdRough = K*diag(Ei)*N^(-1/D);
            % Get roughening samples from Gaussian distribution with StdRough stand. dev.
            deltaX = StdRough*randn(2,N);
            % and add them to posterior particles:
            posteriorParticles(:,:,k) = posteriorParticles(:,:,k) + deltaX;
        end 
    else % if part b), set posterior particles to prior ones
        posteriorParticles(:,:,k) = priorParticles;
    end
    % Leave this part - used for performance assessment and optional output
    if(nargout > 0)
        % evaluate mean error of 50% closest particles:
        dists = sqrt(sum([posteriorParticles(1,:,k) - personStates(1,k);
            posteriorParticles(2,:,k) - personStates(2,k)].^2));
        dsorted = sort(dists,'ascend'); % sort distances
        errs(k) = mean(dsorted(1:floor(length(dsorted)/2)));
    end
end

% If output requested, do not show animation
if(nargout > 0)
    meanErr = mean(errs);
    return
end

%% PLOT RESULTS (You should not need to change anything here).
% Setup
gridding = 0:0.2:10;
% Axis limits
xLim = [0 gridding(end)];
yLim = [0 gridding(end)];
zLim = [-1 5];

% 3D plot
try close(1); catch end % try to close figure 1 first
figure(1)
subplot(1,2,1)
[x_mesh,y_mesh] = meshgrid(gridding, gridding);
[h_mesh,~] = altitudeAndGradient(x_mesh,y_mesh);
surf(x_mesh,y_mesh,h_mesh);
caxis([-1 2.2])
colormap winter
% set view
view([-30,60]);
set(gca,'xlim',xLim);
set(gca,'ylim',yLim);
set(gca,'zlim',zLim);
xlabel('Location x (km)')
ylabel('Location y (km)')
zlabel('Altitude a (km)')
hold on

% 2D plot
subplot(1,2,2)
[x_mesh,y_mesh] = meshgrid(gridding, gridding);
[h_mesh,~] = altitudeAndGradient(x_mesh,y_mesh);
contour(x_mesh,y_mesh,h_mesh);
set(gca,'xlim',xLim);
set(gca,'ylim',yLim);
caxis([-1 2.2])
title('Level-Curves Plot')
xlabel('Location x (km)')
ylabel('Location y (km)')
hold on;

% Set position and width of figure
pos = get(gcf,'Position');
pos(1) = 100;
pos(3) = 1000;
set(gcf,'Position',pos);

hperson = [];
hparticles = [];

% plot particles, if any of the first x location particles are nonzero (for cases b-d)
plotParticles = any(posteriorParticles(1,:,1));

% Show Animation
for k = 1:(Nsim + 1)
    if(plotParticles)
        % plot particles as well
        % 3D:
        subplot(1,2,1)
        try delete(hparticles(1)); catch end
        [hval,~] = altitudeAndGradient(posteriorParticles(1,:,k),posteriorParticles(2,:,k));    % h value
        hparticles(1) = plot3(posteriorParticles(1,:,k),posteriorParticles(2,:,k),hval,'c.','MarkerSize',15);
        % 2D:
        subplot(1,2,2)
        try delete(hparticles(2)); catch end
        hparticles(2) = plot(posteriorParticles(1,:,k),posteriorParticles(2,:,k),'c.','MarkerSize',15);
    end
    % plot actual location of person
    % delete previous location:
    try delete(hperson(1)); catch end
    subplot(1,2,1)
    % Plot person's true state
    [hval,~] = altitudeAndGradient(personStates(1,k),personStates(2,k));    % h value
    % plot person 100m above actual location, so that the person marker is
    % still visible among all PF particles.
    hperson(1) = plot3(personStates(1,k),personStates(2,k),hval+0.1,'r.','Markersize',20);
    title(['Step: k=',int2str(k-1),'   Particles: N=',int2str(N),'']);
    try delete(hperson(2)); catch end
    subplot(1,2,2)
    hperson(2) = plot(personStates(1,k),personStates(2,k),'r.','Markersize',20);
    legend(hperson(2),'Person','Location','SouthWest')
    if(activateButtonPress)
        if(~waitforbuttonpress())
            break;
        end
    else
        pause(pauseSec)
    end
end

%% Altitude and gradient function: (you should not need to change this)

% [a,grad] = altitudeAndGradient(x,y)
% Function that returns the altitude and its gradient given a location
% INPUT:
% x,y:  Coordinates of person, in km
% OUTPUT:
% a:    altitude at location given by x,y, in km
% grad: gradient of terrain at location x,y: [dh/dx; dh/dy](x,y), w.o. unit
    function [a,grad] = altitudeAndGradient(x, y)
        % values and gradient of l (level/elevation).
        coeff = [   -1.394079911146759e+001
            4.318449838008606e+001
            -5.873035185058867e-002
            -1.400961327221415e+001
            4.393262455274894e+001
            -6.301283314674033e-002
            1.936241014498105e+002
            -5.881117127815230e+002
            -5.803913633249534e+002
            1.558126588975868e+003
            1.575320828668508e+002
            1.468258796264717e+002
            2.000000000000000e+003];
        
        
        h = @(x,y) 1*x.^3.*y.^3 + coeff(1)*x.^2.*y.^3 + coeff(2)*x.*y.^3 + coeff(3).*y.^3 ...
            + coeff(4)*x.^3.*y.^2 + coeff(5)*x.^3.*y + coeff(6)*x.^3 ...
            + coeff(7)*x.^2.*y.^2 ...
            + coeff(8)*x.^2.*y + coeff(9)*x.*y.^2 + coeff(10)*x.*y ...
            + coeff(11)*x + coeff(12)*y + coeff(13);
        
        
        dhdx = @(x,y) 3*x.^2.*y.^3 + 2*coeff(1)*x.*y.^3 + coeff(2)*y.^3 + 0*coeff(3) ...
            + 3*coeff(4)*x.^2.*y.^2 + 3*coeff(5)*x.^2.*y + 3*coeff(6)*x.^2 ...
            + 2*coeff(7)*x.*y.^2 ...
            + 2*coeff(8)*x.*y + coeff(9)*y.^2 + coeff(10)*y ...
            + coeff(11) + 0*coeff(12) + 0*coeff(13);
        
        dhdy = @(x,y) 3*x.^3.*y.^2 + 3*coeff(1)*x.^2.*y.^2 + 3*coeff(2)*x.*y.^2 + 3*coeff(3).*y.^2 ...
            + 2*coeff(4)*x.^3.*y + coeff(5)*x.^3 + 0*coeff(6) ...
            + 2*coeff(7)*x.^2.*y ...
            + coeff(8)*x.^2 + 2*coeff(9)*x.*y + coeff(10)*x ...
            + 0*coeff(11) + coeff(12) + 0*coeff(13);
        
        % convert to KM here
        a = 0.001 * h(x,y);
        grad = 0.001 * [dhdx(x,y); dhdy(x,y)];
    end
end