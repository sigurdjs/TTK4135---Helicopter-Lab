%% Some globals to make life easy
global trav_offset pitch_offset

%% Plots of calculated values
figure(1)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','mo'),grid
ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5,'mo'),grid
ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6,'mo'),grid
xlabel('tid (s)'), ylabel('e_dot')

%% Plots of measured values
data = load('data.mat');
figure(2);
plot(t,x1,data.data(1,:),data.data(2,:) + trav_offset);
xlab = xlabel('Time[s]');
ylab = ylabel('Travel[rad]');
h_legend = legend('\lambda_{calculated}','\lambda_{measured}');
set(h_legend,'Fontsize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);

figure(3);
plot(t,x3,t,u1,data.data(1,:),data.data(4,:) + pitch_offset);
xlab = xlabel('Time[s]');
ylab = ylabel('Pitch[rad]');
h_legend = legend('p_{calculated}','u = p_{c}','p_{measured}');
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14)

figure(4);
plot(t,x5,data.data(1,:),data.data(6,:));
xlab = xlabel('Time[s]');
ylab = ylabel('Elevation[rad]');
h_legend = legend('e_{calculated}','e_{measured}');
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14)