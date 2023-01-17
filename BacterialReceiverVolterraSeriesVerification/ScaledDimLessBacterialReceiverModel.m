%% Parameter Defintion 
kappaApi = 10;    % amount time^-1
kApi = 1;       % concentration
deltaC = 2 ;      % time^-1
deltaA = 2 ;      % time^-1
kDp = 1 ;         % time^-1
kP = 10 ;         % concentration^-1 time^-1
alphaA = 1 ; 
alphaC = 1 ; 
kappaB = 0.10 ;   % amount time^-1
deltaCcp = 1 ; 
k1 =  1 ;       % concentration^-1
k2 = 0.1 ;      % concentration^-1
deltaRx = 2 ;   % time^-1
deltaApi = 2;   % time^-1

n = 1000;        
v=1;            % Volume of Bacterial Cell ~1 um^3 = 10^-15 l
 
deltaE = 1;     % time^-1
betaE = 1;      % time^-1


hill_exp = 1;

a_tilde = (n*alphaA*kappaB)/deltaA;

c_tilde = (n*alphaC*kappaB)/deltaC;

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


% Dependent Variable Scale Constant
constCcp = c_tilde;
constApi = a_tilde;
constRx  = c_tilde;
constElux = c_tilde;

% Input Magnitude / Scale Factor
IScale = 100;

% Time Scale Constant
tc = 1/(IScale*c_tilde*k1);
tc = (n*v)/(kP*c_tilde);
constAIP  = c_tilde*IScale;



dstat_vars_dt= zeros(4,1);

time_uLim = 3;

% Time scale for simulation
tspan = [0  time_uLim];
options = odeset('NonNegative',[1,2,3,4]);

% Time Scale for external AIP concentration
time = 0:0.001:time_uLim; 

InputSignal = 0.001; 

% Define AIP function - Dirac
AIP =@(t) InputSignal * (dirac(t)==inf) ;


% Define AIP function - Heaviside
AIP =@(t) InputSignal * heaviside(t);
AIP =@(t) InputSignal*heaviside(t) - InputSignal*heaviside(t-2); % TODO Verifiy Time Units for Input Time Shift
%AIP =@(t) InputSignal*exp(-(t-2).^2)
AIP =@(t) (InputSignal/10).*t;% - heaviside(t-5).*((dirac_const/5)*(t-5));
AIP =@(t) InputSignal*triangularPulse(0,1,2,t);


% Define AIP function - Heaviside
%   

% Define AIP function - Sinuoid
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



[t,stat_vars] = ode45(@(t,stat_vars) VanillaScaledBacterialReceiverSystem(t, stat_vars, AIP, k1, k2,...
        kappaApi, kApi,deltaRx , deltaCcp, deltaApi, ...
        kP,kDp,deltaE, betaE , n,v, ...
        hill_exp, c_tilde, a_tilde, ...
        constCcp , constApi , constRx , ...
        constElux, constAIP, tc),...
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
title("QS Signal: " + num2str(InputSignal) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
legend(lgd_strings(i));
end

scaleArr = [constCcp, constApi, constRx, constElux];

 
for i = 1:1:size_stat_vec(2)
figure(99),subplot(2,3,i);
plot(t*tc,scaleArr(i)*stat_vars(:,i), strcat("-x",color_strings(i), ""));
%ylim([0  0.01])
%set(gca, 'YScale', 'log')
sgtitle("ODE Solution - Re-Scaled to Real Units" + num2str(InputSignal) +" [1/Bacteria Volume]", ...
    'Interpreter', 'latex', 'FontSize', 15);
xlabel('Time Unit ', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Normalized Molecule Number ', 'Interpreter', 'latex', ...
    'FontSize', 15);
legend(lgd_strings(i));
end




%% Scaled Volterra Verification - Full (1)



%figure(101), plot(stat_vars(:,2)./(kApi*n*v))

% Input is Heaviside
%input_fct = @(s) (dirac_const./s);
syms tau s
test_input = laplace(AIP(tau));

% Replace any occurance of erfc(.) with 1-erfz() to allow for complex input
% TODO: Fix
% if contains(string(test_input),"erfc")
%     
% syms symArg    
% symERFZ = sym(sin(symArg));
% 
% [extrArgIdxStart,extrArgIdxEnd] = regexp(string(test_input), "erfc\(.*[^\)]\)");
%      
%     extrArg = extractBetween(string(test_input),extrArgIdxStart+5,extrArgIdxEnd-1);
% %     test_input= sym(...
% %         replaceBetween(...
% %         string(test_input),extrArgIdxStart, extrArgIdxEnd, "(1-erfz("+extrArg+"))"...
% %         ));
%         test_input= subs(...
%         test_input,str2sym("erfc("+(extrArg)+")"), str2sym("(1-erfz("+extrArg+"))")...
%         );
% 
% end

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
% 

order = 0b0001;

% 1-order Volterra
nilt1 = NILT1_Volterra(ScaledBacRecvFirstOrderVolterra(order,n,v, q1, kappaApi, ...
    kApi,deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP ...
    , kDp, deltaApi, betaE, deltaE,input_fct, ...
    tc, constCcp, constApi, constAIP, constRx, constElux), ...
    tspan(2));
time1 = linspace(0,tspan(2),length(nilt1));
figure(2), subplot(3,1,1), plot(time1, nilt1, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('1D Numerical Inverse Laplace Transform  $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);
% 2-order Volterra
nilt2 = NILT2_Volterra(ScaledBacRecvSecondOrderVolterra(order,n,v, q2,kappaApi, ...
    kApi,deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP, ...
    kDp, deltaApi, betaE, deltaE,input_fct, ...
    tc, constCcp, constApi, constAIP, constRx, constElux), ...
    tspan(2));
time2 = linspace(0,tspan(2),length(nilt2));
subplot(3,1,2), plot(time2, nilt2, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('2D Numerical Inverse Laplace Transform - $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);

%% Volterra Verification - V3
% 3-order Volterra
nilt3 = NILT3_Volterra(ScaledBacRecvThirdOrderVolterra(order,n,v, q3,kappaApi, kApi, ...
    deltaRx, c_tilde, k1,k2, deltaCcp, a_tilde,kP, kDp, deltaApi, ...
    betaE, deltaE,input_fct, ...
    tc, constCcp, constApi, constAIP, constRx, constElux), ...
    tspan(2));
time3 = linspace(0,tspan(2),length(nilt3));
subplot(3,1,3), plot(time3, nilt3, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Concentration [molecules/ml]', 'Interpreter', 'latex', ...
    'FontSize', 15);
title('3D Numerical Inverse Laplace Transform -  $' + convertCharsToStrings(printIFunc)+"$", ...
    'Interpreter', 'latex', 'FontSize', 15);



%% Scaled Volterra Verification - Full (2)
% Comparison of ODE outputs and Volterra series
outputAux = stat_vars(:,sum(bitget(order,1:4)));
outputReal = interp1(t, outputAux, time);
figure(3), plot(time, outputReal, 'b', 'LineWidth', 2);
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



%% Scaled ODE Definition

function dstat_vars_dt = VanillaScaledBacterialReceiverSystem(t, stat_vars, AIP, k1, k2,...
        kappaApi, kApi,deltaRx , deltaCcp, deltaApi,...
    kP, kDp,deltaE, betaE , n,v, hill_exp, c_tilde, a_tilde, constCcp , constApi , constRx , constElux, constAIP, tc)

% (3) - Concentration of AgrC-AIP complex


% d_CCp_dt - SCALED
dstat_vars_dt(1,1) = tc * (c_tilde)/(constCcp) * k1 * (AIP(t) * constAIP) -... 
                    k1*tc*stat_vars(1) * (AIP(t) * constAIP) - ...
                    (k2 + deltaCcp)*stat_vars(1)*tc;

% (4) - Concentration of phosphorylated AgrA

% d_APi_dt
dstat_vars_dt(2,1) =(1/(n*v)) * tc * kP * a_tilde * (constCcp)/(constApi) * stat_vars(1) - ...
                    (1/(n*v)) * (kP) * tc * constCcp * stat_vars(1)* stat_vars(2) -...
                    (kDp+deltaApi) * tc * stat_vars(2);

% (5) - Model of Promoter Complex Binding

% d_R(x)_dt
dstat_vars_dt(3,1) = n*kappaApi *...
    (((stat_vars(2))/(n*v))^(hill_exp)/ ((kApi/constApi)^(hill_exp)+((stat_vars(2))/(n*v)^(hill_exp)))) * tc  / (constRx)  ...
    - deltaRx * tc * stat_vars(3) - betaE * stat_vars(3) * tc ; 



 % d_ELux_dt- Linear
dstat_vars_dt(4,1)  = betaE * constRx * stat_vars(3) * tc * 1/(constElux) -...
                    deltaE * tc * stat_vars(4);



end

