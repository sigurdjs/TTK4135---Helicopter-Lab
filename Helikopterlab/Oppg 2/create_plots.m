load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/u_q1.mat');
load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/travel_meas_q1.mat');
load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/travel_des_q1.mat');
load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/t.mat');
load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/pitch_meas_q1.mat');
load('/Users/Sigurd/Google Drive/1. Skole/3.2 TTK4135/Helicopter Lab/Helikopterlab/Oppg 2/pitch_des_q1.mat');


clf;
figure(1);
plot(t(1:80),(x1(1:80).*(180/pi)),trav_meas(1,:),trav_meas(2,:)+165);
xlab = xlabel('Time [s]');
ylab = ylabel('Travel [deg]');
h_legend = legend('\lambda_{calculated}','\lambda_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);

figure(2);
plot(t(1:80),(x3(1:80).*(180/pi)),t(1:80),(u(1:80).*(180/pi)),pitch_meas(1,:),pitch_meas(2,:));
xlab = xlabel('Time [s]');
ylab = ylabel('Pitch [deg]');
h_legend = legend('p_{calculated}','u','p_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);