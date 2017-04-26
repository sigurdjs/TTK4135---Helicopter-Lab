% Script for plotting results in task 10.2
load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/pitch_des_q10.mat');
load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/pitch_meas_q10.mat');
load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/u_q10.mat');

load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/travel_des_q10.mat');
load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/travel_meas_q10.mat');
load('/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/Helikopterlab_ny/Oppg 2/t.mat');

clf;
fig1 = figure(1);
plot(t(1:80),(x1(1:80)),trav_meas(1,:) ...
    ,(trav_meas(2,:).*(pi/180)) + pi);
xlab = xlabel('Time [s]');
ylab = ylabel('Travel [deg]');
h_legend = legend('\lambda_{calculated}' ...
    ,'\lambda_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);

fig2 = figure(2);
plot(t(1:80),(x3(1:80)),t(1:80), ...
    (u(1:80)),pitch_meas(1,:),(pitch_meas(2,:).*(pi/180)));
xlab = xlabel('Time [s]');
ylab = ylabel('Pitch [deg]');
h_legend = legend('p_{calculated}','u','p_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);

cleanfigure('handle',fig1);
cleanfigure('handle',fig2);

matlab2tikz('filename','/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/report/fig/2_plots/travel_q10.tikz','figurehandle',fig1,'width', '5.3cm');
matlab2tikz('filename','/home/sigurdjs/Documents/TTK4135---Helicopter-Lab/report/fig/2_plots/pitch_q10.tikz','figurehandle',fig2,'width', '5.3cm');

close all;