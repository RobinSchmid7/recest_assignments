% Counterexample to the separation principle for LTI system and DISTRIBUTED
% control implementation without perfect knowledge of the control input.
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
% [30.05.13, ST]    First version
%

clear;


%% Problem parameters

A = [1 0.3; 1.2 1.1];
B = [-1 1; 1 0];
H = [1 0];
K = [-0.2; 1.4];
F = [-1.2 -0.2; -1.3 -0.5];


%% Stability tests

% State-feedback controller
disp('Magnitude of eigenvalues of A+BF:');
disp(abs(eig( A+B*F )));

% State observer
disp('Magnitude of eigenvalues of (I-KH)A:');
disp(abs(eig( (eye(length(A))-K*H)*A )));

% "Extra" eigenvalues of distributed implementation
disp('Magnitude of eigenvalues of (I-KH)(A+BF):');
disp(abs(eig( (eye(length(A))-K*H)*(A+B*F) )));