% TTK4135 - Helicopter lab
% Updated spring 2017, Andreas L. Fl?ten

%% Initialization and model definition
init;
run = true; % set to true to auto-compile and run model
delta_t	= 0.25; % sampling time
h = delta_t;
% Discrete time system model. x = [lambda r p p_dot]'
A1 = [1 h 0 0; 0 1 -(K_2*h) 0; 0 0 1 h; ...
    0 0 -(h*K_1*K_pp) 1-(h*K_1*K_pd)];
B1 = [0; 0; 0; (h*K_1*K_pp)];

% Number of states and inputs
mx = size(A1,2); % Number of states
mu = size(B1,2); % Number of inputs

% Initial values
x1_0 = pi;                       % Lambda
x2_0 = 0;                        % r
x3_0 = 0;                        % p
x4_0 = 0;                        % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0]';     % Initial values

% Time horizon and initialization
N  = 100;                        % Time horizon for states
M  = N;                          % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);         % Initialize z
z0 = z;                          % Initial value
z0(1,1) = pi;

% Bounds
ul 	    = -30*pi/180;   % Lower bound on control -- u1
uu 	    = 30*pi/180;    % Upper bound on control -- u1

xl      = -Inf*ones(mx,1);  % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);   % Upper bound on states (no bound)
xl(3)   = ul;               % Lower bound on state x3
xu(3)   = uu;               % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = genbegr2(N,M,xl,xu,ul,uu);
vlb(N*mx+M*mu)  = 0;       % We want the last input to be zero
vub(N*mx+M*mu)  = 0;       % We want the last input to be zero

% Generate the matrix Q and the vector c
% (objecitve function weights in the QP problem)
Q1 = zeros(mx,mx);
Q1(1,1) = 1;               % Weight on state x1
Q1(2,2) = 0;               % Weight on state x2
Q1(3,3) = 0;               % Weight on state x3
Q1(4,4) = 0;               % Weight on state x4
P1 = 0.5;                  % Weight on input
Q = 2*genq2(Q1,P1,N,M,mu); % Generate Q
c = zeros(N*mx+M*mu,1);    % Generate c

%% Generate system matrixes for linear model
Aeq = gena2(A1,B1,N,mx,mu);      % Generate A
beq = zeros(N*mx,1);        	 % Generate b
beq(1:mx) = A1*x0;               % Initial value

%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub);
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end
%% Find K matrix with LQR
pitch_lim = uu;
Q_lq = diag([10,16/(pi^2),16/(pi^2),4/(pi^2)]);
R_lq = 1/(pitch_lim^2);
K = dlqr(A1,B1,Q_lq,R_lq);

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)];

x1 = [x0(1);z(1:mx:N*mx)];  % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];  % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];  % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];  % State x4 from solution


num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

t = 0:delta_t:delta_t*(length(u)-1);

%% Simulation
if run == true
    qc_build_model;
    qc_start_model;
    pause(30);
    qc_stop_model;
    pause(5);
end

%% Plotting
plotting_3;   % Small script for plotting simulations
