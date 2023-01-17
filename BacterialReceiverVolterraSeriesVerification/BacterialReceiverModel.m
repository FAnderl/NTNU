% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %                                               
% Bacterial Receiver Model                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% ResultsComputation
% System of Ordinary, non-linear differential equations
% describing bacterial receiver 
% adopted from http://dx.doi.org/10.1063/1.4884519
% Bacterial population is modeled

% Gram - negative bacteria


% Model needs to include:
% 1) Production of surface/ membrane receptors
% 2) Production of internal, unphosphorylated AGrA; 
% 3) External AIP binding to receptors
% 4) Phosphorylation of interal Agr A
% 5) 
%


% ------------- Substances ---------- % 
% AgrA ...  A
%
% Phosphorylated AgrA ... A_Pi
%
% Membrance Receptors ... R
%
% External AIP Concentration ... P
%
% Receptor+P Complex ... C_CP
%
% Transcription Product ... (RNA) RNA
%
% (Bio-Luminescent Protein BLp)
% 

% ------------- Parameters ---------- % 


% k1 ... receptor-AIP complex formation/association rate
%
% k2 ... receptor-AIP complex DISassociation rate
%
% deltaCcp ... receptor-AIP complex degradation rate
%
% kappaB ... basal transctiption rate
%
% alphaC ... effective factor of AgrC synthesis
%
% alphaA ... effective factor of AgrA synthesis
%
% kP ... Phosphorylation rate og AgrA

% k_Dp ... DEphosphorylation rate og AgrA
%
% deltaA ... degradation of A
%
% deltaApi ... degradation of APi
%
% deltaC ... degradation of C
%
% kApi ... Saturation konstant of Api_induced RNA synthesis
%
% kappaApi ... maximum rate of APi-induced RN promotion & synthesis rate 
 
% Parameters taken from original Publication
%% Parameter Defintion 
kappaApi = 10;    % amount time^-1
kApi = 1;       % concentration % TODO: Verify!
deltaC = 2 ;      % time^-1
deltaA = 2 ;      % time^-1
kDp = 1 ;         % time^-1
kP = 10 ;         % concentration^-1 time^-1
alphaA = 1 ; 
alphaC = 1 ; 
kappaB = 1 ;   % amount time^-1
deltaCcp = 2 ; 
k1 =  1 ;       % concentration^-1 time^-1
k2 = 0.1 ;      % concentration^-1 
deltaRx = 2 ;   % time^-1
deltaApi = 2;   % time^-1

n = 1000;        
v=1;            % Volume of Bacterial Cell ~1 um^3 = 10^-15 l
 
deltaE = 1;     % time^-1
betaE = 0.1;      % time^-1


hill_exp = 1;

% ---->
a_tilde = (n*alphaA*kappaB)/deltaA;

c_tilde = (n*alphaC*kappaB)/deltaC;

%% Scaled System Parameter Defition - TODO: Remove (Deprecated)

Gamma  = ((kDp+deltaApi)/(c_tilde*k1));

Chi = (k2+deltaCcp)/(c_tilde*k1);

epsilon = 1/(c_tilde*k1);


%% Taylor Factor Definition 
% ------------------------------------

% Definition of Taylor Series - Volterra Factors % Hard Coded for a= k/2;
% TODO -> Deprecated; Remove
sigma1 = 4/(9*kApi);
sigma2 = -8/(27*kApi^2);
sigma3 = 16/(81*kApi^3);
sigma4 = -32/(243*kApi^4);



sigma1 = (1/(kApi/2 + kApi*n*v) - kApi/(2*(kApi/2 + kApi*n*v)^2)) / factorial(1);
sigma2 = (kApi/(kApi/2 + kApi*n*v)^3 - 2/(kApi/2 + kApi*n*v)^2)   / factorial(2);
sigma3 = (6/(kApi/2 + kApi*n*v)^3 - (3*kApi)/(kApi/2 + kApi*n*v)^4) / factorial(3);
sigma4 = ((12*kApi)/(kApi/2 + kApi*n*v)^5 - 24/(kApi/2 + kApi*n*v)^4) / factorial(4);

% Expansion Point for Taylor Series is defined and used here 
z = 0; % default is 0 

q1 = sigma1 - 2*z*sigma2 + 3*z^2*sigma3 - 4*z^3*sigma4;
q2 = sigma2 - 3*z*sigma3 + 6*z^2*sigma4;
q3 = sigma3 - 4*sigma4*z;


% ------------- ODE system ----------- % 


% RNA transcription follows saturation kinetics (Michaelis-Menten)



% (1) - Concentration of membrane bound receptors AgrC




% (2) - Concentration of UN-phosphorylated AgrA




% (3) - Concentration of AgrC-AIP complex




% (4) - Concentration of phosphorylated AgrA




% (5) - Model of Promoter Complex Binding




% (6) - "Do I need another formula to model the process from promoter to transcription" 






%% Regular System
% The Bacterial Receiver System including all 6 involved state variables 

dstat_vars_dt= zeros(6,1);

% Time scale for simulation
tspan = [0  20];
options = odeset('NonNegative',[1,2,3,4,5,6]);

% Time Scale for external AIP concentration
time = 0:0.01:10; 

dirac_const = 60;  % 100nM 

% Define AIP function - Dirac
AIP =@(t) dirac_const*(dirac(t)==inf) ;


% Define AIP function - Heaviside
AIP =@(t) dirac_const*heaviside(t);



% Define AIP function - Heaviside
AIP =@(t) dirac_const*heaviside(t) - dirac_const*heaviside(t-2);

% Define AIP function - Sinuoid
%sAIP= @(t)sin(2*pi*2*t)*(heaviside(t)-heaviside(ti-2));
%AIP= @(t)awgn(sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2)), 3)  ;
figure(99), fplot(AIP)

% Constant supply
%AIP =@(t) ((10.*t)./(1.+t)) ;
%figure(1), plot(time, feval(AIP,time))

% Initial Conditons

% Units mol/l
% stat_vars = [
%     C_Cp
%     APi
%     R(x)
%     C
%     A
%     ELux
% ]; 

init_cond = [0.0;0.0;0.0; (n*alphaC*kappaB)/deltaC; (n*alphaA*kappaB)/deltaA; 0];
%init_cond = [0.0;0.0;0.0; 0; 0; 0];




[t,stat_vars] = ode45(@(t,stat_vars) bacterial_receiver_sys(t, stat_vars,AIP,  k1, k2,...
        kappaApi, kApi,deltaC,deltaA,deltaRx ,deltaApi, ...
    kDp,kP,alphaA,alphaC,kappaB,deltaCcp,deltaE, betaE ,n,v, hill_exp), tspan, init_cond, options);


varNameStrings = ["Bound Surface Receptor", "Phosphorylized Promoter", "mRNA","Unbound Receptor", "De-Phosphorylized Promoter", "Target Protein"];
lgd_strings = ["C_{CP}","APi","R(x)","C", "A", "E_{Lux}"];
color_strings = ["r","b","g","c", "m", "r"];

%close all

size_stat_vec = size(stat_vars);

for i = 1:1:size_stat_vec(2)
figure(1), subplot(2,3,i)
plot(t,stat_vars(:,i), strcat("-x",color_strings(i), ""));
%ylim([0  0.01])
%set(gca, 'YScale', 'log')
xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title("QS Signal: " + num2str(dirac_const) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(lgd_strings(i));
end




%% Reduced System

% Simplified and Reduced Model for Bacterial Receiver 
dstat_vars_dt= zeros(4,1);

% Time scale for simulation
tspan = [0  5];
options = odeset('NonNegative',[1,2,3,4]);

% Time Scale for external AIP concentration
time = 0:0.01:5; 

dArr = 6 ;  % Array of Input Concentrations
InArr = cell(size(dArr,2));  
OutArr = cell(size(dArr,2),1);
tArr = cell(size(dArr,2),1);

for k = 1:size(dArr,2)

dirac_const = dArr(k); 

% Define AIP function - Dirac
AIP =@(t) dirac_const*(dirac(t)==inf);


% Define AIP function - Heaviside
AIP =@(t) dirac_const*heaviside(t);

%AIP =@(t) dirac_const*exp(-(t-2).^2);
AIP =@(t) (dirac_const/10).*t;% - heaviside(t-5).*((dirac_const/5)*(t-5));
AIP =@(t) dirac_const*triangularPulse(0,1,2,t);
AIP =@(t) dirac_const*heaviside(t) - dirac_const*heaviside(t-2);

InArr{k} = AIP;

% Define AIP function - Heaviside
%AIP =@(t) dirac_const*heaviside(t) - dirac_const*heaviside(t-2);

% Define AIP function - Sinusoid
% AIP= @(t)sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2));
% AIP= @(t)awgn(sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2)), 3)  ;




% Initial Conditons

% Units mol/l
% stat_vars = [
%     C_Cp
%     APi
%     R(x)
%     ELux
% ]; 

init_cond = [0.0;0.0;0.0;0];

% Steady State Solution for C & A
a_tilde = (n*alphaA*kappaB)/deltaA;

c_tilde = (n*alphaC*kappaB)/deltaC;


 [t,stat_vars] = ode45(@(t,stat_vars) ReducedBacterialReceiverSys(t, stat_vars, ...
         AIP, k1, k2,...
        kappaApi, kApi,deltaRx,...
        deltaApi, kDp, kP, deltaCcp,...
        deltaE, betaE, ...
        n,v, hill_exp, c_tilde, a_tilde),...
        tspan, init_cond, options);


OutArr{k} = stat_vars;
tArr{k} = t;


varNameStrings = ["Bound Surface Receptor", "Phosphorylized Promoter", "mRNA", "Target Protein"];
lgd_strings = ["C_{CP}","APi","R(x)", "E_{Lux}"];
color_strings = ["r","b","g", "r"];

%close all

size_stat_vec = size(stat_vars);

for i = 1:1:size_stat_vec(2)
figure(k), subplot(2,3,i)
plot(t,stat_vars(:,i), strcat("-x",color_strings(i), ""));
%ylim([0  0.01])
%set(gca, 'YScale', 'log')
xlabel('Time', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title("QS Signal: " + num2str(dirac_const) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(lgd_strings(i));
end

end

RedSysRes = {t,stat_vars};

%% Max - Reduced System 
% NB! Max Reduced system can be used for assumption 
% Deploys quasi steady-state assumption for C_Cp

% Simplified and Reduced Model for Bacterial Receiver 
dstat_vars_dt= zeros(2,1);

% Time scale for simulation
tspan = [0  5];
options = odeset('NonNegative',[1,2]);

% Time Scale for external AIP concentration
time = 0:0.01:5; 


dirac_const = 6; 

% Define AIP function - Dirac
% Ugly re-definition; FIXME 
AIP =@(t) dirac_const*(dirac(t)==inf);
AIP =@(t) dirac_const*heaviside(t);

%AIP =@(t) (dirac_const/6000).*t;% - heaviside(t-5).*((dirac_const/5)*(t-5));
AIP =@(t) dirac_const*triangularPulse(0,1,2,t);
AIP =@(t) dirac_const*heaviside(t) - dirac_const*heaviside(t-2);

% MICHAELIS-MENTEN-KINETICS -> Steady-State Solutions
k_AIP = (deltaCcp+k2)/(k1);


k_CCP = ((kDp+deltaApi)/(kP*(n*v)^-1));



% Quasi-Steady-State Solution for Surface - Recptors
CCp_qss =@(t) (c_tilde.*AIP(t))./((k_AIP) + AIP(t));
A_Pi_qss = @(t) (a_tilde*c_tilde.*CCp_qss(t))./...
                                         (c_tilde.*CCp_qss(t) + k_CCP.*(k_AIP+CCp_qss(t)));

init_cond = [0.0;0.0];

[t,stat_vars] = ode45(@(t,stat_vars) MaxReducedBacterialReceiverSys(t, stat_vars, ...
         CCp_qss, k1, k2,...
        kappaApi, kApi,deltaRx,...
        deltaApi, kDp, kP, deltaCcp,...
        deltaE, betaE, ...
        n,v, hill_exp, c_tilde, a_tilde),...
        tspan, init_cond, options);




varNameStrings = [ "mRNA", "Target Protein"];
RefvarNameStrings = ["R_{X} - Reference", "E_{Lux} - Reference"];
lgd_strings = ["R_{(x)}", "E_{Lux}"];
color_strings = ["g", "r"];

%close all

size_stat_vec = size(stat_vars);

for i = 1:1:size_stat_vec(2)
figure(10), subplot(1,4,i)
plot(t,stat_vars(:,i), strcat("-x",color_strings(i), ""), "DisplayName",lgd_strings(i));
%ylim([0  0.01])
%set(gca, 'YScale', 'log')
xlabel('Time', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title("QS Signal: " + num2str(dirac_const) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
%legend(lgd_strings(i));

hold on
subplot(1,4,i)
t_ref = RedSysRes{1};
sv_ref =  RedSysRes{2};
plot(t_ref,sv_ref(:,i+2), strcat("-x","b", ""),"DisplayName",RefvarNameStrings(i));
%legend(RefvarNameStrings(i));
legend show
hold off

end

subplot(1,4,3)
plot(time, CCp_qss(time), strcat("-x","b", ""));
legend("C_{Cp} - QSS");

subplot(1,4,4)
plot(time, A_Pi_qss(time), strcat("-x",color_strings(i), ""));
legend("A_{Pi} - QSS");

%% DIMENSIONLESS & SCALED SYSTEM


dstat_vars_dt= zeros(4,1);

time_uLim = 500;

% Time scale for simulation
tspan = [0  time_uLim];
options = odeset('NonNegative',[1,2,3,4]);

% Time Scale for external AIP concentration
time = 0:0.001:time_uLim; 

dirac_const = 10000/(c_tilde*n*v); 

% Define AIP function - Dirac
P =@(t) dirac_const*(dirac(t)==inf) ;


% Define AIP function - Heaviside
P =@(t) dirac_const*heaviside(t);


% Define AIP function - Heaviside
%AIP =@(t) dirac_const*heaviside(t) - dirac_const*heaviside(t-2);

% Define AIP function - Sinuoid
%AIP= @(t)sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2));
%AIP= @(t)awgn(sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2)), 3)  ;


%figure(99), fplot(AIP)

% Constant supply
%AIP =@(t) ((10.*t)./(1.+t)) ;
%figure(1), plot(time, feval(AIP,time))

% Initial Conditons

% Units mol/l
% stat_vars = [
%     C_Cp
%     APi
%     R(x)
%     ELux
% ]; 


init_cond = [0.0;0.0;0.0;0];



 [t,stat_vars] = ode45(@(t,stat_vars) ScaledBacterialReceiverSystem(t, stat_vars, P, k1,...
        kappaApi, kApi,deltaRx , ...
    kP,deltaE, betaE , n,v, hill_exp, c_tilde, a_tilde, Gamma, Chi, epsilon),...
    tspan, init_cond, options);


varNameStrings = ["Bound Surface Receptor", "Phosphorylized Promoter", "mRNA", "Target Protein"];
lgd_strings = ["C_{CP}","APi","R(x)", "E_{Lux}"];
color_strings = ["r","b","g", "r"];

%close all

size_stat_vec = size(stat_vars);

for i = 1:1:size_stat_vec(2)
figure(2), subplot(2,3,i)
plot(t,stat_vars(:,i), strcat("-x",color_strings(i), ""));
%ylim([0  0.01])
%set(gca, 'YScale', 'log')
xlabel('Time Unit ', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Normalized Molecule Number ', 'Interpreter', 'latex', ...
    'FontSize', 15);
title("QS Signal: " + num2str(dirac_const) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(lgd_strings(i));
end






%% Volterra Verification - Saturation ONLY   %%%
% Verification of Volterra Kernel for saturation term ONLY (NB!)   

close all

tspan = [0  20];
time = 0:0.01:20; 
test_const = dirac_const;
        

% Isolated Differential Equation with non-linear saturation term 
xinit = [0  0];


disp("k-Parameter: " +num2str(dirac_const/(kApi*n*v)));
disp("Point of Development: " +num2str(z));

if (test_const/(n*v)) > kApi
    warning("Conditions for Series Extension NOT fulfilled")
end

iFunc = @(t)(test_const/5)*t;

f =  @(t,x) [(test_const/5);
    n*kappaApi*((x(1))/ ((kApi*n*v)+(x(1)))) ...
- deltaRx*x(2)];  % [Total Amount [molecules]]
[t, y] = ode45(f, tspan, xinit);
figure,  plot(t, y(:,2), 'r-o');
title('Reduced System Output')

% Input is Heaviside
%input_fct = @(s) (dirac_const./s);

syms tau
test_input = laplace(iFunc(tau));
input_fct = matlabFunction(test_input);
%input_fct = @(s) test_const./s;



% 1-order Volterra
nilt1 = NILT1_Volterra(ReducedBacRecvFirstOrderVolterra(n,v,q1, kappaApi, kApi,deltaRx, input_fct), tspan(2));
time1 = linspace(0,tspan(2),length(nilt1));
figure, subplot(3,1,1), plot(time1, nilt1, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('1D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);

% 2-order Volterra
nilt2 = NILT2_Volterra(ReducedBacRecvSecondOrderVolterra(n,v,q2, kappaApi, kApi,deltaRx, input_fct), tspan(2));
time2 = linspace(0,tspan(2),length(nilt2));
subplot(3,1,2), plot(time2, nilt2, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('2D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);


% 3-order Volterra
nilt3 = NILT3_Volterra(ReducedBacRecvThirdOrderVolterra(n,v, q3,kappaApi, kApi,deltaRx, input_fct), tspan(2));
time3 = linspace(0,tspan(2),length(nilt3));
subplot(3,1,3), plot(time3, nilt3, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('3D Numerical Inverse Laplace Transform', ...
    'Interpreter', 'latex', 'FontSize', 15);




% Comparison of ODE outputs and Volterra series
outputAux = y(:,end)';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Volterra Verification - Full (1)

% Steady-  Solutions for A and C for equation system under assumption
% of constant a_tilde & c_tilde

a_tilde = (n*alphaA*kappaB)/deltaA;
 

c_tilde = (n*alphaC*kappaB)/deltaC;


%figure(101), plot(stat_vars(:,2)./(kApi*n*v))
for k=1:size(dArr,2)

% Input is Heaviside
%input_fct = @(s) (dirac_const./s);
syms tau s
helperFctHdl = InArr{k};
test_input = laplace(helperFctHdl(tau));

%test_input = laplace(AIP(tau));
input_fct = matlabFunction(test_input);
printIFunc = latex(test_input);


% Binary Value that determines to which step of the processing chain in the
% bacertial receiver the Volterra kernels are computed 
%
% 0001: C_Cp
% 0011: A_Pi
% 0111: R_(X)
% 1111: E_Lux

order = 0b1111;

% 1-order Volterra
nilt1 = NILT1_Volterra(BacRecvFirstOrderVolterra(order,n,v, q1, kappaApi, kApi,deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP, kDp, deltaApi, betaE, deltaE,input_fct), tspan(2));
time1 = linspace(0,tspan(2),length(nilt1));
figure(size(dArr,2)+k), subplot(3,1,1), plot(time1, nilt1, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('1D Numerical Inverse Laplace Transform  $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);
% 2-order Volterra
nilt2 = NILT2_Volterra(BacRecvSecondOrderVolterra(order,n,v, q2,kappaApi, kApi,deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP, kDp, deltaApi, betaE, deltaE,input_fct), tspan(2));
time2 = linspace(0,tspan(2),length(nilt2));
subplot(3,1,2), plot(time2, nilt2, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('2D Numerical Inverse Laplace Transform - $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);

%% Volterra Verification - V3
% 3-order Volterra
nilt3 = NILT3_Volterra(BacRecvThirdOrderVolterra(order,n,v, q3,kappaApi, kApi, deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP, kDp, deltaApi, betaE, deltaE,input_fct), tspan(2));
time3 = linspace(0,tspan(2),length(nilt3));
subplot(3,1,3), plot(time3, nilt3, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('3D Numerical Inverse Laplace Transform -  $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);



%% Volterra Verification - Full (2)
% Comparison of ODE outputs and Volterra series
outputAux = OutArr{k}'; %TODO changed from end -> 1 -> 2 -> make automatic
outputAux = outputAux(sum(bitget(order,1:4)),:);
outputReal = interp1(tArr{k}, outputAux, time);
figure(2*size(dArr,2)+(k)), plot(time, outputReal, 'b', 'LineWidth', 2);
outputVolterra1 = interp1(time1, (nilt1), time);
outputVolterra2 = interp1(time2, (nilt1 + nilt2), time);
outputVolterra3 = interp1(time3, (nilt1 + nilt2 + nilt3), time);
hold on, plot(time, outputVolterra1, 'r-.', 'LineWidth', 2);
hold on, plot(time, outputVolterra2, 'g--', 'LineWidth', 2);
hold on, plot(time, outputVolterra3, 'c-', 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
%ylim([0,max(outputVolterra1)]);
title('Comparison of ODE output and output of Volterra series- $' + convertCharsToStrings(printIFunc)+"$",...
    'Interpreter', 'latex', 'FontSize', 15);
leg3 = legend('ODE', 'Volterra series ($n=1$)', ...
    'Volterra series ($n=2$)', 'Volterra series ($n=3$)'); 
set(leg3, 'Interpreter','latex','FontSize',15);
hold off

% figure(100)
% diff1 = outputVolterra1-outputReal;
% diff2 = (outputVolterra1 + outputVolterra2)-outputReal;
% diff3 = outputVolterra1 + outputVolterra2 + outputVolterra3-outputReal;
% hold on, plot(time, diff1, 'r-.', 'LineWidth', 2);
% hold on, plot(time, diff2, 'g--', 'LineWidth', 2);
% hold on, plot(time, diff3, 'c-', 'LineWidth', 2);
% leg4 = legend('diff1', 'diff2', ...
%     'diff3');
% set(leg4, 'Interpreter','latex','FontSize',15);
end


%% Save Figures
FolderName = "./figures/01";   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
FigHandle = FigList(iFig);
FigAxes   = findall(FigHandle, 'type', 'axes');
FigName = get(FigAxes(1), 'title').String;
FigName = strrep(FigName, '/', "\");
savefig(FigHandle, fullfile(FolderName, FigName+ '.fig'));
end

%% ODE System Definition
function dstat_vars_dt = bacterial_receiver_sys(t, stat_vars, AIP, k1, k2,...
        kappaApi, kApi,deltaC,deltaA,deltaRx ,deltaApi, ...
    kDp,kP,alphaA,alphaC,kappaB,deltaCcp,deltaE, betaE , n,v, hill_exp)

% (3) - Concentration of AgrC-AIP complex


% d_CCp_dt
dstat_vars_dt(1,1) = k1* stat_vars(4) * AIP(t) - k2*stat_vars(1) - deltaCcp* ...
    stat_vars(1); 


% (4) - Concentration of phosphorylated AgrA

% d_APi_dt
dstat_vars_dt(2,1) =(1/(n*v)) *kP* stat_vars(1) * stat_vars(5) - kDp* stat_vars(2) - ...
        deltaApi* stat_vars(2); 


% (5) - Model of Promoter Complex Binding

% d_R(x)_dt
dstat_vars_dt(3,1) = n*kappaApi*(((stat_vars(2)^(hill_exp))/(n*v))/ (kApi^(hill_exp)+((stat_vars(2)^(hill_exp))/(n*v)))) ...
- deltaRx*stat_vars(3); 

% ----------- Uncoupled -------------



% (1) - Concentration of membrane bound receptors AgrC

% d_C_dt - Non - Linear
dstat_vars_dt(4,1)  = n*alphaC*kappaB - k1*AIP(t)*stat_vars(4) + k2*stat_vars(1)...
    - deltaC*stat_vars(4); 



% (2) - Concentration of UN-phosphorylated AgrA

% d_A_dt - Non Linear
dstat_vars_dt(5,1)  = n*alphaA*kappaB - (kP/(n*v))*stat_vars(1)*stat_vars(5) + kDp*stat_vars(2)...
    - deltaA*stat_vars(5)  ; 

% d_ELux_dt- Linear
dstat_vars_dt(6,1)  = betaE * stat_vars(3) - deltaE*stat_vars(6);



end


function dstat_vars_dt = ReducedBacterialReceiverSys(t, stat_vars, ...
         AIP, k1, k2,...
        kappaApi, kApi,deltaRx,...
        deltaApi, kDp, kP, deltaCcp,...
        deltaE, betaE, ...
        n,v, hill_exp, c_tilde, a_tilde)




        % d_CCp_dt
        dstat_vars_dt(1,1) = k1 * (c_tilde-stat_vars(1))* AIP(t) - k2*stat_vars(1) - deltaCcp* ...
            stat_vars(1); 
        
        
        %Concentration of phosphorylated AgrA
        % d_APi_dt
        dstat_vars_dt(2,1) =(1/(n*v)) * kP * stat_vars(1) * (a_tilde-stat_vars(2)) - kDp* stat_vars(2) - ...
                deltaApi* stat_vars(2); 
        
        
       
        
        % d_R(x)_dt
        dstat_vars_dt(3,1) = n*kappaApi*(((stat_vars(2)^(hill_exp))/(n*v))/ (kApi^(hill_exp)+((stat_vars(2)^(hill_exp))/(n*v)))) ...
        - deltaRx*stat_vars(3) - betaE * stat_vars(3);         % TODO: Verify that, i.e., Translation degrades mRNA
        
     
        
         
        
        % d_ELux_dt- Linear
        dstat_vars_dt(4,1)  = betaE * stat_vars(3) - deltaE*stat_vars(4);


        

    

end


function dstat_vars_dt = MaxReducedBacterialReceiverSys(t, stat_vars, ...
         A_Pi_qss, k1, k2,...
        kappaApi, kApi,deltaRx,...
        deltaApi, kDp, kP, deltaCcp,...
        deltaE, betaE, ...
        n,v, hill_exp, c_tilde, a_tilde)




        
        %Concentration of phosphorylated AgrA
        % d_APi_dt
%         dstat_vars_dt(1,1) =(1/(n*v)) * kP * CCp_qss(t) * (a_tilde-stat_vars(1)) - kDp* stat_vars(1) - ...
%                 deltaApi* stat_vars(1); 
        
        
       
        
        % d_R(x)_dt
        dstat_vars_dt(1,1) = n*kappaApi*(((A_Pi_qss(t)^(hill_exp))/(n*v))/ (kApi^(hill_exp)+((A_Pi_qss(t)^(hill_exp))/(n*v)))) ...
        - deltaRx*stat_vars(1)-betaE * stat_vars(1); 
        
     
        
         
        
        % d_ELux_dt- Linear
        dstat_vars_dt(2,1)  = betaE * stat_vars(1) - deltaE*stat_vars(2);


        

end


% TODO: Moved in Vanilla Form to another .m-file
% Implementation of DIMENSIONLESS & Scaled Bacterial Receiver System 

function dstat_vars_dt = ScaledBacterialReceiverSystem(t, stat_vars, P, k1,...
        kappaApi, kApi,deltaRx , ...
    kP,deltaE, betaE , n,v, hill_exp, c_tilde, a_tilde, Gamma, Chi, epsilon)

% (3) - Concentration of AgrC-AIP complex


% d_CCp_dt - SCALED
dstat_vars_dt(1,1) = P(t) - stat_vars(1)*P(t) - Chi*stat_vars(1); 

% (4) - Concentration of phosphorylated AgrA

% d_APi_dt
dstat_vars_dt(2,1) =(1/(n*v)) *(kP/k1)*stat_vars(1) - (1/(n*v)) *(kP/k1)*stat_vars(1)* stat_vars(2) - Gamma*stat_vars(2);

% (5) - Model of Promoter Complex Binding

% d_R(x)_dt
dstat_vars_dt(3,1) = n*(kappaApi/(k1*c_tilde^2))*(((stat_vars(2)^(hill_exp))/(n*v))/ ((kApi/a_tilde)^(hill_exp)+((stat_vars(2)^(hill_exp))/(n*v)))) ...
- deltaRx*epsilon*stat_vars(3); 



 % d_ELux_dt- Linear
dstat_vars_dt(4,1)  = betaE * epsilon * stat_vars(3) - deltaE*epsilon*stat_vars(4);



end

