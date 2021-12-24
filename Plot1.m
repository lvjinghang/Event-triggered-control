
Design_K_P;
dynamic_ETM;
dynamic_relative_ETM;

figure(1)
plot(t,x1(:,1),'b-','linewidth',1.5)
hold on
plot(t,x2(:,1),'r--','linewidth',1.5)
hold off
xlabel('Time(s)');
mlstr1 = {'$x_1(t)$'};
ylabel(mlstr1,'interpreter','latex');
axis([0 tm -5 3])
% mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')

figure(2)
plot(t,x1(:,2),'b-','linewidth',1.5)
hold on
plot(t,x2(:,2),'r--','linewidth',1.5)
hold off
xlabel('Time(s)');
mlstr1 = {'$x_2(t)$'};
ylabel(mlstr1,'interpreter','latex');
axis([0 tm -2 3])
mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')


figure(3)
plot(t,x1(:,3),'b-','linewidth',1.5)
hold on
plot(t,x2(:,3),'r--','linewidth',1.5)
hold off
xlabel('Time(s)');
mlstr1 = {'$x_3(t)$'};
ylabel(mlstr1,'interpreter','latex');
axis([0 tm -1 5])
mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')

figure(4)
plot(t,u1,'b-','linewidth',1.5)
hold on
plot(t,u2,'r--','linewidth',1.5)
hold off
xlabel('Time(s)');
mlstr1 = {'Control input u(t)'};
ylabel(mlstr1,'interpreter','latex');
axis([0 tm -10 10])
mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')

figure(5)
subplot(3,1,1)
plot(t,phi1,'b-','linewidth',1.0)
axis([0 tm 0.05 0.25])
xlabel('Time(s)');
mlstr = {'$\varphi_{1}(t)$'};
legend(mlstr,'interpreter','latex')
subplot(3,1,2)
plot(t,phi2,'r-','linewidth',1.0)
axis([0 tm 0.03 0.14])
xlabel('Time(s)');
mlstr = {'$\varphi_{2}(t)$'};
legend(mlstr,'interpreter','latex')
subplot(3,1,3)
plot(t,delta_r,'k-','linewidth',1.0)
axis([0 tm 0.07 0.13])
xlabel('Time(s)');
mlstr = {'$\vartheta(t)$'};
legend(mlstr,'interpreter','latex')

figure(6)
stem(t,intervals1,'b-','linewidth',0.5)
hold on
stem(t,intervals2,'r.','linewidth',0.5)
axis([0 tm 0 1]) 
xlabel('Time(s)');
ylabel('Release intervals');
mlstr = {'Dynamic-threshold ETM';'Dynamic-threshold relative ETM'};
legend(mlstr,'interpreter','latex')

figure(7)
plot(t,d1,'b-','linewidth',1.5)
hold on
plot(t,d1_es,'m-.','linewidth',1.5)
hold off
xlabel('Time(s)');
ylabel('Disturbance');
axis([0 tm -6 6])
mlstr = {'$\omega_{11}$';'The estimation of $\omega_{11}$'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')

figure(8)
plot(t,d12,'b-','linewidth',1.5)
hold on
plot(t,d1_es2,'m-.','linewidth',1.5)
hold off
xlabel('Time(s)');
ylabel('Disturbance');
axis([0 tm -6 6])
mlstr = {'$\omega_{11}$';'The estimation of $\omega_{11}$'};
%mlstr = {'Fixed threshold';'Relative threshold';'Switching threshold';'The method in [20]'};
legend(mlstr,'interpreter','latex')