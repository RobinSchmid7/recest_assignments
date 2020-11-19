% Observer-based (LQG) controller design for the Balancing Cube.
%
% The observer-based control system consists of 
%   - a steady-state Kalman Filter, and
%   - an LQR state-feedback controller.
% 
% Implementation/tasks:
% - Design steady-state KF
% - Design an LQR controller
% - Analysis: estimator poles, controller poles, control performance (RMS 
%   of states), control effort (RMS of inputs), input peak values
% - Simulation: observer-based controller vs. state-feedback controller
%
% Course: Recursive Estimation, Spring 2013
% Problem Set: Observer-Based Control and the Separation Principle
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% S. Trimpe
% strimpe@ethz.ch
% 2013
%
% --
% Revision history
% [28.05.13, ST]    First version
%

close all;
clear;


%% Configuration
% Simulation time:
T = 1000;    % []=samples (i.e. simulation time T/100 seconds)

% Select what states to plot (select 4 indices).
% See below (next cell) for meanings of state indices.
sel = [5 6 13 14];


%% Load Balancing Cube model.
% Contained in sys_d.
% Discrete-time model with sampling rate of 100 Hz.
%
% States:
%   x_1: angle arm 1
%   x_2: angular velocity arm 1
%   x_3: angle arm 2
%   ...
%   x_11: angle arm 6
%   x_12: angular velocity arm 6
%   x_13: cube roll angle
%   x_14: cube roll angular velocity
%   x_15: cube pitch angle
%   x_16: cube pitch angular velocity
%
% Measurements:
%   z_1: angle arm 1
%   z_2: cube roll angular velocity
%   z_3: cube pitch angular velocity
%   z_4: angle arm 2
%   ...
%   z_16: angle arm 6
%   z_17: cube roll angular velocity
%   z_18: cube pitch angular velocity
%
% Inputs:
%   u_1: torque applied at arm 1
%   u_2: torque applied at arm 2
%   ...
%   u_6: torque applied at arm 6
load CubeModel;

A = sys_d.a;
B = sys_d.b;
H = sys_d.c;

n = length(A);      % Number of states
n_meas = size(H,1); % Number of measurements
n_arms = 6;         % Number of arms (control agents) = number of inputs


%% Initial condition and noise variance

% Initial condition
x0 = zeros(n,1);
x0(13) = 1/180*pi;  % assume start at 0, except pitch angle initially at 1 deg
P0 = 0*eye(n);      % variance = 0 -> perfect knowledge of initial state

% Noise and initial state statistics
R = 1e-6*diag(repmat([1 0 0],1,n_arms)) ...  % absolute encoder noise (small number)
    + 2e-5*diag(repmat([0 1 1],1,n_arms));      % rate gyro noise
Q = 1e-7*eye(16);


%% Design a steady-state Kalman Filter.
% Compute the steady-state Kalman Filter gain from the solution of the
% Discrete Algebraic Riccati Equations.

P = dare(A',H',Q,R);
K_inf = P*H'*inv(H*P*H'+R);

% Steady-state KF A-matrix
A_sskf = (eye(n)-K_inf*H)*A;

% Steady-stte KF B-matrix
B_sskf = (eye(n)-K_inf*H)*B;


%% Design an LQR controller

% Weights
Q_bar = 10*eye(n);
R_bar = 0.1*eye(n_arms);    
% Change to R_bar = 1*eye(n_arms), for example, to roughly meet the 
% actuation limitations.

% Discrete-time LQR design:
F = -dlqr(A,B,Q_bar,R_bar);


%% Simulation and Kalman Filter Implementation
% Variable to store state and state estimate
x1 = zeros(n,T+1);              % Observer-based control
x1(:,1) = x0+sqrt(P0)*randn(n,1); 
x2 = zeros(n,T+1);              % State-feedback control
x2(:,1) = x1(:,1); 

xHat = zeros(n,T+1);   % Steady-state KF
xHat(:,1) = x0;

% Control input
u1 = zeros(n_arms,T+1);    % for observer-based control implementation
u2 = zeros(n_arms,T+1);    % for state-feedback control implementation

Q_sq = sqrt(Q);
R_sq = sqrt(R);

for k=1:T
    % I: Observer-based control
    %
    % (i) Simulate system:
    x1(:,k+1) = A*x1(:,k) + B*u1(:,k) + Q_sq*randn(16,1);
    z = H*x1(:,k+1) + R_sq*randn(n_meas,1);
    %
    % (ii) Steady-state Kalman Filter:
    xHat(:,k+1) = A_sskf*xHat(:,k) + B_sskf*u1(:,k) + K_inf*z;
    %
    % (iii) Observer-based control:
    u1(:,k+1) = F*xHat(:,k+1);
    
    % II: State-feedback control
    % Notice that this is a hypothetical implementation since, in reality,
    % we do not have access to perfect state measurements.
    %
    % (i) Simulate system:
    x2(:,k+1) = A*x2(:,k) + B*u2(:,k) + Q_sq*randn(16,1);
    %
    % (ii) State-feedback control:
    u2(:,k+1) = F*x2(:,k+1);
   
end;


%% Analysis
% Compute the poles of the error dynamics for the steady-state KF.
poles_estimator = eig(A_sskf);
disp(['Magnitude of estimator eigenvalues: ']);
disp(sort(abs(poles_estimator)));

% Compute poles of LQR state-feedback controller
poles_controller = eig(A+B*F);
disp(['Magnitude of controller eigenvalues: ']);
disp(sort(abs(poles_controller)));

% Compute control performance.
disp(['Control performance (RMS of states).']);
disp(['Observer-based control: ',num2str(sqrt(trace(x1*x1')/T))]);
disp(['State-feedback control: ',num2str(sqrt(trace(x2*x2')/T))]);
disp(' ');

% Compute control effort.
disp(['Control effort (RMS of inputs).']);
disp(['Observer-based control: ',num2str(sqrt(trace(u1*u1')/T))]);
disp(['State-feedback control: ',num2str(sqrt(trace(u2*u2')/T))]);
disp(' ');

% Maximum control input.
disp(['Peak values of input:']);
disp(['Observer-based control: ',num2str(max(abs(u1),[],2)')]);
disp(['State-feedback control: ',num2str(max(abs(u2),[],2)')]);


%% Plots
% - some states (indices in SEL)
% - all inputs

% Figure 1: states
figure;
for i=1:4
    subplot(4,1,i);
    plot(0:T,x1(sel(i),:)/pi*180,0:T,x2(sel(i),:)/pi*180);
    grid;
    if i==1
        title('States and state estimates (in deg or deg/s)');
    end;
    ylabel(['x_{',int2str(sel(i)),'}']);
end;
xlabel('Discrete-time step k');
legend('Observer-based control','State-feedback control');

% Figure 2: control inputs
figure;
for i=1:6
    subplot(6,1,i);
    plot(0:T,u1(i,:),0:T,u2(i,:));
    grid;
    if i==1
        title('Control input: applied torque (in N*m)');
    end;
    ylabel(['u_{',int2str(i),'}']);
end;
xlabel('Discrete-time step k');
legend('Observer-based control','State-feedback control');
