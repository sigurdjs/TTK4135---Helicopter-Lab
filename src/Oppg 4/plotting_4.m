figure(4);
plot(t,x5,data.data(1,:),data.data(6,:));
xlab = xlabel('Time[s]');
ylab = ylabel('Elevation[rad]');
h_legend = legend('e_{calculated}','e_{measured}');
set(xlab,'Fontsize',14);
set(ylab,'Fontsize',14)