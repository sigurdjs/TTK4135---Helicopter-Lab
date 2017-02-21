% % TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2017, Andreas L. Fl?ten
% 
%% Initialization and model definition
init; % NB: Change this to the init file corresponding to your helicopter
global N mx mu beta alpha lambda_t 
delta_t	= 0.25; % sampling time
run = true; % set to true to auto-compile and run model

%% Definition of discrete time system model. x = [lambda r p p_dot]'
disc_model;

%% Generate system matrixes for linear model
Aeq = gena2(A1,B1,N,mx,mu);           % Generate A, hint: gena2
beq = zeros(N*mx,1);        	  % Generate b
beq(1:mx) = A1*x0; % Initial value

%% Generate nonlinear inequality constraints
alpha = 0.5;
beta = 20;
lambda_t = (2*pi)/3;

%% Solve QP problem with linear model
obj_fun = @(X) X'*Q*X;
options = optimoptions('fmincon','Algorithm','sqp','Display','iter','UseParallel',true);
tic
z = fmincon(obj_fun,z0,[],[],Aeq,beq,vlb,vub,@nonlinear_constraint,options);
% [z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub); % hint: quadprog
t1=toc;

%% Find K matrix with LQR
Q_lq = diag([10 1/(16*pi^2) 16/(pi^2) 4/(pi^2) 4/(pi^2) 4/(pi^2)]);
R_lq = diag([(2*pitch_lim_pos)^2 pitch_lim_pos^2]);
K = dlqr(A1,B1,Q_lq,R_lq);

%% Extract control inputs and states
u1  = z(N*mx+1:mu:N*mx+M*mu-1); % Control input from solution
u2  = z(N*mx+2:mu:N*mx+M*mu);

x1 = [x0(1);z(7:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(8:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(9:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(10:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(11:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(12:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

t = 0:delta_t:delta_t*(length(u1)-1);

%% Simulation
if run == true
    qc_build_model;
    qc_start_model;
    pause(30);
    qc_stop_model;
    pause(5);
end

%% Plotting
global trav_offset pitch_offset 
trav_offset = 0;
pitch_offset  = 0;

plotting;   % Small script for plotting simulations


