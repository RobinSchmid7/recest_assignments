% Problem Set: Observers and the Separation Principle - Problem 2b
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
%

clc;
clear;
close all;

%% System matrices
% State space model:
%   x(k) = A*x(k-1) + B*u(k)
%   z(k) = C*x(k)
A = [6 0.5 1 6.5; 5.25 -0.5 1.25 6.5; -5 0.5 0 -5.5; -1 -0.5 0 -0.5];
B = [0.5; 0.5; 0.5; -0.5];
H = [1 0 1 2; 0 0 1 1];


%% Check observability
% You can check the rank of the matrix obsv(A,C): if it is full column
% rank, the system is observable.  For practical problems however, its
% better to look at the singular values of obsv(A,C) and check if some of
% them are close to zero (relative to the other singular values).
Ob = obsv(A,H);

disp('Singular values of observability matrix: (system is observable if all are different from 0)');
disp(svd(Ob));


%% Design observer

% Part 1: decompose the system into observable and unobservable subsystems.
% Can use the command 'obsvf' to bring the system to observability
% staircase form.
[Abar,Bbar,Hbar,T,k] = obsvf(A,B,H);

% The observable subsystem is represented by the lower right 2x2 block of 
% Abar and the right 2x2 block of Hbar.  We will design an observer for
% this part.  The unobservable part we can't do anything about, but
% fortunately, it is stable.

p = [0.8 0.7];              % Desired pole locations.

Ao = Abar(3:4,3:4);
Ho = Hbar(1:2,3:4);
K = place(Ao',(Ho*Ao)',p)';

K = [zeros(2,2); K];        % Add zeros for states not used in observer.
K = T'*K;                   % Transform back to original state order.

disp('Observer gain K:');
disp(K);

% Display the eigenvalues
disp('Eigenvalues of observer dynamics:');
disp(eig((eye(4)-K*H)*A));


%% Simulation
% The closed loop system (see Problem Set for derivation).
% As outputs, we define the actual states and the observer error.
Acl = [A, zeros(4,4); K*H*A, (eye(4)-K*H)*A];
Bcl = [B; B];
Ccl = [eye(4), zeros(4,4); eye(4), -eye(4)];
sys_cl = ss(Acl,Bcl,Ccl,[],1);

% State transformation to (x(k), e(k)).  Doesn't change the input-output
% behaviour of the system; if we simulate with the original system however, 
% numerical issues occur.
T = [eye(4), zeros(4,4); eye(4), -eye(4)];      % transformation matrix
sys_cl = ss2ss(sys_cl,T);

% Response to x(0)=(1,1,1,1) and \hat{x}(0)=0, hence e(0)=(1,1,1,1).
initial(sys_cl,[1 1 1 1 1 1 1 1]');