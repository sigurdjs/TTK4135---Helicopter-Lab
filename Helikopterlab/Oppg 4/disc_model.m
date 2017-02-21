Ac = [0 1 0 0 0 0;
      0 0 -K_2 0 0 0;
      0 0 0 1 0 0;
      0 0 -K_1*K_pp -K_1*K_pd 0 0;
      0 0 0 0 0 1;
      0 0 0 0 -K_3*K_ep -K_3*K_ed];
  
Bc = [0 0; 
      0 0;
      0 0;
      K_1*K_pp 0;
      0 0;
      0 K_3*K_ep];

A1 = eye(6) + delta_t.*Ac;
B1 = delta_t.*Bc;
  

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
% x = [lambda, r = lambda_dot, p, p_dot, e, e_dot]
x0 = [pi 0 0 0 0 0]'; 

% Time horizon and initialization
N  = 80;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization
z0(1,1) = pi;

% Bounds
pitch_lim_pos = 30*pi/180;
pitch_lim_neg = -30*pi/180; 
ul 	    = [pitch_lim_neg; -Inf];                  % Lower bound on control -- u1
uu 	    = [pitch_lim_pos; Inf];                    % Upper bound on control -- u1

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = pitch_lim_neg;                           % Lower bound on state x3
xu(3)   = pitch_lim_pos;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu); % hint: genbegr2
vlb(N*mx+M*mu)  = 0;                         % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                         % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = diag([1 0 0 0 0 0]);
P1 = diag([0.8 0.03]);                                 % Weight on input
Q = genq2(Q1,P1,N,M,mu);              % Generate Q
%c = zeros(N*mx+M*mu,1);                 % Generate c