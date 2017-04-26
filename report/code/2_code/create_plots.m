% Script for plotting results in task 10.2
load('path to where .mat data files are located');

clf;
figure(1);
plot(t(1:80),x1(1:80),trav_meas(1,:) ...
    ,trav_meas(2,:));
xlab = xlabel('Time [s]');
ylab = ylabel('Travel [rad]');
h_legend = legend('\lambda_{calculated}' ...
    ,'\lambda_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);

figure(2);
plot(t(1:80),x3(1:80),t(1:80), ...
    u(1:80),pitch_meas(1,:),pitch_meas(2,:));
xlab = xlabel('Time [s]');
ylab = ylabel('Pitch [rad]');
h_legend = legend('p_{calculated}','u','p_{measured}');
set(h_legend,'FontSize',14);
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14);
