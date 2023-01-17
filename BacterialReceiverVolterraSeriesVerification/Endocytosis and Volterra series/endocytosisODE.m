% Parameters from Wattis-paper (for LDL) % mol = molecules = particles
ke = 0; % [mol/ml/sec]
Vmedium = 10; %[ml]
Vcell = 0.0066; %[ml]
alpha = Vmedium/Vcell;
a = 6.64e-17; % [mL/mol/sec]
pm = 200;
N = 180;
ki = 0.0027; % [1/sec]
kd = 0.0002; % [1/sec]
% Time domain ODE
f = @(t,x) [ke - a/alpha*x(1)*(pm*N - x(2));
    a*x(1)*(pm*N - x(2)) - ki*x(2);
    ki*x(2) - kd*x(3)
    ];
% % ODE with initial condition for the biding
% tspan = [0,1e4];
% xinit = [1.17e13 1.6e13 0]; % [mol/ml]
% [t, y] = ode45(f, tspan, xinit);
% figure, ode45(f, tspan, xinit);
% leg1 = legend('$c_e(t)$', '$c_b(t)$', '$c_r(t)$');
% set(leg1, 'Interpreter','latex','FontSize',15);
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
% ylabel('Exosome concentration [molecules/ml]', 'Interpreter', 'latex', ...
%     'FontSize', 15);
% title('With $C_{e_0} = 1.17\times10^{13}$ and $C_{b_0} = 1.6\times10^{13}$ [mol/ml]', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% ODE without initial condition for the biding
tspan = [0,6e4];
InConst = 1.17e13;
xinit = [InConst 0 0]; % [mol/ml]
[t, y] = ode45(f, tspan, xinit);
figure, ode45(f, tspan, xinit);
leg2 = legend('$c_e(t)$', '$c_b(t)$', '$c_r(t)$');
set(leg2, 'Interpreter','latex','FontSize',15);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Exosome concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('With $C_{e_0} = 1.17\times10^{13}$ and $C_{b_0} = 0$ [mol/ml]', ...
    'Interpreter', 'latex', 'FontSize', 15);

% % Linear system c_e --> c_b (First-order kernels) (from our paper)
% num1 = [a*pm*N];
% den1 = [1 ki];
% sys1 = tf(num1, den1);
% % Linear system c_b --> c_r (Linear is anyway) (from our paper)
% num2 = [ki];
% den2 = [1 kd];
% sys2 = tf(num2, den2);
% % Linear endocytosis (First-order kernels)
% sysLin = sys1*sys2;
% [num,den] = tfdata(sys2);
% 
% input = y(:,1)';
% time = 0:1:t(end);
% inputInterp = interp1(t, input, time);
% figure, lsim(sysLin, inputInterp, time);

% Test - impulse responses - first three orders 
% (valildated for the 1st-order impulse response)
% Test for the Heaviside input (approximating c_e from ODEs)

InConst;

% 1-order Volterra
nilt1 = NILT1_Volterra(firstOrderFuncEndocytosis(a, pm, N, ki, kd, InConst), tspan(2));
time1 = linspace(0,tspan(2),length(nilt1));
figure, subplot(3,1,1), plot(time1, nilt1, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('1D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);
% 2-order Volterra
nilt2 = NILT2_Volterra(secondOrderFuncEndocytosis(a, pm, N, ki, kd, InConst), tspan(2));
time2 = linspace(0,tspan(2),length(nilt2));
subplot(3,1,2), plot(time2, nilt2, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('2D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);
% 3-order Volterra
nilt3 = NILT3_Volterra(thirdOrderFuncEndocytosis(a, pm, N, ki, kd, InConst), tspan(2));
time3 = linspace(0,tspan(2),length(nilt3));
subplot(3,1,3), plot(time3, nilt3, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('3D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);

%% Comparison of ODE outputs and Volterra series
outputAux = y(:,3)';
time = 0:1:t(end);
outputReal = interp1(t, outputAux, time);
figure, plot(time, outputReal, 'b', 'LineWidth', 2);
outputVolterra1 = interp1(time1, (nilt1), time);
outputVolterra2 = interp1(time2, (nilt1 + nilt2), time);
outputVolterra3 = interp1(time3, (nilt1 + nilt2 + nilt3), time);
hold on, plot(time, outputVolterra1, 'k-.', 'LineWidth', 2);
hold on, plot(time, outputVolterra2, 'k--', 'LineWidth', 2);
hold on, plot(time, outputVolterra3, 'k-', 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('Comparison of ODE output and output of Volterra series', ...
    'Interpreter', 'latex', 'FontSize', 15);
leg3 = legend('ODE', 'Volterra series ($n=1$)', ...
    'Volterra series ($n=2$)', 'Volterra series ($n=3$)');
set(leg3, 'Interpreter','latex','FontSize',15);
