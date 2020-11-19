% Problem Set: Observers and the Separation Principle - Problem 2a
%
% Observer design
%
% Recursive Estimation
% Spring 2011
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Sebastian Trimpe
% strimpe@ethz.ch
%
% --
% Revision history
% [30.05.11, ST]    first version
% [31.05.11, ST]    changed to new simulation task

clc;
clear;
close all;

%% System matrices
% State space model:
%   x(k) = A*x(k-1) + B*u(k)
%   z(k) = C*x(k)
A = [0 -6; 1 5];
B = [1; 1];
H = [0 1];



%% Check observability
% You can check the rank of the matrix obsv(A,C): if it is full column
% rank, the system is observable.  For practical problems however, its
% better to look at the singular values of obsv(A,C) and check if some of
% them are close to zero (relative to the other singular values).
Ob = obsv(A,H);

disp('Singular values of observability matrix: (system is observable if all are different from 0)');
disp(svd(Ob));


%% Design observer

% Desired poles of observer dynamics: (I-K*H)*A = A-K*H*A
% Notice that place places the poles of A-K*H; hence we have to call it
% with H*A as the second argument.  Furthermore, place is written for pole
% placement for static controller feedback (the dual problem of the static
% observer design), therefore input and output matrices are transposed 
% (also see place documentation).
p = [0.8 0.7];              % Desired pole locations.
K = place(A',(H*A)',p)';

disp('Observer gain K:');
disp(K);

% Display the eigenvalues
disp('Eigenvalues of observer dynamics:');
disp(eig((eye(2)-K*H)*A));


%% Simulation
% The closed loop system (see Problem Set for derivation).
% As outputs, we define the actual states and the observer error.
Acl = [A, zeros(2,2); K*H*A, (eye(2)-K*H)*A];
Bcl = [B; B];
Ccl = [eye(2), zeros(2,2); eye(2), -eye(2)];
sys_cl = ss(Acl,Bcl,Ccl,[],1);

% State transformation to (x(k), e(k)).  Doesn't change the input-output
% behaviour of the system; if we simulate with the original system however, 
% numerical issues occur.
T = [eye(2), zeros(2,2); eye(2), -eye(2)];      % transformation matrix
sys_cl = ss2ss(sys_cl,T);

% Response to x(0)=(1,1) and \hat{x}(0)=0, hence e(0)=(1,1).
initial(sys_cl,[1 1 1 1]');