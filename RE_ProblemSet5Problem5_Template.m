function meanErr = RE_ProblemSet5Problem5_Template()
% meanErr = RE_ProblemSet5Problem5_Template()
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
N = 10;
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
theta = 1/12; % Horizontal velocity of person, in km (per time step)
Q = (theta/4)^2*eye(2); % Process noise variance in kilometers^2
s0 = [7.5; 7.5]; % Mean of initial state distribution in kilometers
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
% Run this function after you implemented this to check the results.

% Draw random initial state s(0) of the person here: %%%%%%%%%%%%%%%%%%%%%%
personStates(:,1) = s0;

for k = 2:(Nsim+1) % for k = 1:Nsim (2:Nsim+1 because Matlab is 1-based indexing)
    % Implement the update of the person's state here: %%%%%%%%%%%%%%%%%%%%
    % Hint, use [a,grad] = altitudeAndGradient(x,y) to get the altitude and
    % gradient at a location x,y (see definition of this nested function at
    % the end of this file.
    personStates(:,k) = personStates(:,k-1);
    % Generate a measurement here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    measurements(k) = 0;
end

%% CHANGE CODE HERE: Particle Filter, Parts b) - d)
% You can run and test your code after implementing each part.

% Implement particle initialization here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
posteriorParticles(:,:,1) = zeros(2,N);

for k = 2:(Nsim + 1)
    % b) Implement process update here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    priorParticles = posteriorParticles(:,:,k-1); % replace this
    
    % c) Implement measurement update here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    posteriorParticles(:,:,k) = priorParticles; % replace this
        
    % d) Implement roughening here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw random samples to add to posterior particles
    deltaX = zeros(2,N); % replace this
    % Add them to posterior particles
    posteriorParticles(:,:,k) = posteriorParticles(:,:,k) + deltaX;
 
    % Do not change this - used for performance assessment
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