% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %                                               
% Bacterial Receiver Model                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% Model of XIP Signaling Pathway in Strepptococcus mutans


% Parameter Typical Order of Magnitudes
%
%
% mRNA Degradation: 0.02 1/min 
% Protein Degradation: 0.001 1/min
% Transcription Rate: 
%  
%
%
%


% Parameter Declaration & Definition


Df = 0.54288*1e-2;                % Facilitated active diffusion constant / transfer rate / conductivity
qD = Df;                         % Renaming Routine
delta_XIP = 0.00035;
delta_XIP_int = delta_XIP;       % Internal XIP Degradation Rate Constant
delta_XIP_ext = delta_XIP;       % External XIP Degradation Rate Constant
k_TF_m = 1e-3;                     % Transcription Factor Synthesis Rate Constant
delta_tf = 0.01;                 % Free Transcription Factor (ComR) Degradation Rate Constant
k_TF_f = 0.00000001;                  % Transcription Factor Dematuration Rate Constant; has to be really, really smal
delta_mtf = 1e3;                 % Mature Transcription Factor Degradation Rate

alpha_comR = 0.5;                % Translation Efficieny of comR mRNA to ComR
kappaB = 1e-2;                   % Basal RNA synthesis rate

kappa_comR =  200;               % Maximum ComR-XIP_int-complex-induced sigX transcription rate
k_comR = 0.001;                      % Michaelis-Menten Constant

alpha_sigX = 1;                % SigX Translation Efficiency

kappaX = 1000;                      % Maximum SigX induced comR Transcription Rate

k_SigX = 0.1;                   % 0.5 (TRANSIENT EXPLANATION); Michaelis-Menten Constant

deltaSigX = 10;                 % 10 (TRANSIENT EXPLANATION)

deltaR = 0.02;                   % mRNA Degradation Rate

alpha_L = 1;

kL = 0.05;                       % kEl: Luciferase Translation Rate
deltaEL = 0.00385;               % Luciferase Decay Rate Constant 

v = 1;                          % Bacteria Volume ~1e-15l or 1um^3

mu0 = 0.0090222;                % Bacterial Growth Constant           
Kn = 7.5602;                    %1e8 (SCALE)  % Carrying Capacity (for bacterial gorowth function)

od600_coeff = 1.64e-9;

Ve = 2.2e11;                    % UNITS [um^3] Extracellular Volume; arbitrary 10ml; TODO  -> How to handle this?

tEnd = 900;                   % Time Scale of Experiment


n_init = 0.016  * (1/od600_coeff);      % Initial Bacterial Cell Count DEFAULT:0.016 
%n_init = 1;          % Initial Bacterial Cell Count

n = n_init;                     % Population Size; remnant albeit necessary 



% Steady-State comR-concentration (due to basal production)
comR_tilde = (n*alpha_comR*kappaB)/delta_tf;




% Extracellular XIP - Driving CONCENTRATION
XIP_Sig_Num = 60;  % 100nM, i.e., ~60 molecules per bacterial cell volume
%% Model for Gene Expression





% Time scale for simulation
tspan = [0  tEnd];


% Time Scale for external AIP concentration
time = 0:0.1:tEnd; 


% Comment/Uncomment XIP Signal Input 

% Define XIP function - Heaviside
XIP_ext =@(t) XIP_Sig_Num * heaviside(t);


% Define XIP function - Rect Pulse
XIP_ext =@(t) XIP_Sig_Num*heaviside(t) - XIP_Sig_Num*heaviside(t-2);

% Define XIP function - Sinusoid
% AIP= @(t)sin(2*pi*2*t)*(heaviside(t)-heaviside(t-2));


% Define XIP realistically
%XIP_ext=@(t) XIP_Sig_Num*exp(-(t).^2);


% Initial Conditons

% Units mol/l
% stat_vars = [
%     XIP_int
%     fcomR
%     mcomR
%     RNA
% ]; 

% Model with(0) or without(1) mRNA
SIMPLIFY_FLAG = 1; 

% Dynamic Bacterial Growth Model  
DPS_FLAG = 1;

if SIMPLIFY_FLAG
    MODEL_SZ = 6; 
    options = odeset('NonNegative',[1,2,3,4,5,6]);
    init_cond = [ 0.0;
              comR_tilde;
              0.0;
              0.0;
              0.0;
              XIP_Sig_Num;
              ];

elseif ~SIMPLIFY_FLAG
    MODEL_SZ = 7;
    options = odeset('NonNegative',[1,2,3,4,5,6,7]);
    init_cond = [ 0.0;
              comR_tilde;
              0.0;
              0.0;
              0.0;
              XIP_Sig_Num;
              0];
end




dn_dt = zeros(1,1);

% Additional Scaling Parameters for Bacterial Growth
a =  1;
b =  1;
tau = 80.5504;
n_init_DDE  = n_init;

% Solve for Bacterial Growth
%[t,n_t] = ode45(@(t,n_t) BacterialGrowth(t,n_t,mu,Kn, od600_coeff,), tspan, n_init);


% Solve Delay Differential Equation for Heuristic Fit 
sol = dde23(@(t,n_t,Z) BacGrowthDelayDE(t, n_t, Z, mu0,Kn, a,b),...
        tau,@(t) BacHist(t,mu0,n_init_DDE),[0, 900]);


figure(1);
plot(sol.x,sol.y,'-o', 'LineWidth',3.0,'Color','m')
title('Bacterial Growth - Modelled as Delay Differential Equation','FontSize',15);
xlabel("Time",'FontSize',15);
ylabel("Population size [Number of Cells]",'FontSize',15)

% Use Dynamic Bacterial Growth Model
if DPS_FLAG
  XIP_model = @SM_XIP_DPS;
  n = sol; % n_sol

elseif ~DPS_FLAG
  XIP_model = @SM_XIP;

end


% ---------- Solve for mRNA -----------

% Initialize State Variable Vector
dstat_vars_dt= zeros(MODEL_SZ,1);




[t,stat_vars] = ode15s(@(t,stat_vars) XIP_model(t, stat_vars,XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX,n ,v, Ve, kL, deltaEL, alpha_L)... 
    , tspan, init_cond, options);




varNameStrings = ["Internal XIP","Free comR" ,"Mature TF comR", "Luciferase", "SigX", "External XIP", "(mRNA)"];
lgd_strings = ["XIP_{int}","comR","comR-XIP", "Luciferase", "SigX", "XIP_{ext}", "(mRNA)"];
color_strings = ["r","b","g","k", "c", "m", "r"];


size_stat_vec = size(stat_vars);

figure(2);
for i = 1:1:size_stat_vec(2)
    subplot(2,4,i)
    plot(t,stat_vars(:,i), strcat("-x",color_strings(i), ""));
    xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
        'FontSize', 15);
    title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);
    legend(lgd_strings(i));
end





%% Scaled Model for Bio Luminescence 


n_c = Kn;
v_c = 1;

 % Scaling Parameters
 XIP_int_c =   XIP_Sig_Num*(Kn*1e8);
 XIP_ext_c =   XIP_Sig_Num;                 % Concentration
 ComRf_c   =   XIP_Sig_Num*(Kn*1e8);        % This is for sure wrong
 ComRm_c   =   XIP_Sig_Num*(Kn*1e8); 
 sigX_c    =   XIP_Sig_Num*(Kn*1e8); 
 SigX_c    =   XIP_Sig_Num*(Kn*1e8);
 E_L_c     =   XIP_Sig_Num*(Kn*1e8) ;   
 tc        =   sigX_c/kappa_comR;

tEnd_SCALED = 0.001;

comR_tilde_SCALED = 1;     % UNITY

XIP_ext_init_SCALED = 1;   % UNITY


% Time scale for simulation
SCALED_tspan = [0  tEnd_SCALED];


% Time Scale for external AIP concentration
SCALED_time = 0:0.1:tEnd_SCALED; 


SCALED_options = odeset('NonNegative',[1,2,3,4,5,6,7]);


SCALED_init_cond = [ 0.0;
          0;
          0.0;
          0.0;
          0.0;
          XIP_ext_init_SCALED;
          0];



[t_SCALED,stat_vars_SCALED] = ode15s(@(t_SCALED,stat_vars_SCALED) SCALED_DIMLESS_SM_XIP(t_SCALED, stat_vars_SCALED ...
     ,qD, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR,  kappaX, k_SigX, deltaSigX,sol,v, Ve,kL, deltaEL, ...
     XIP_int_c, ...
     XIP_ext_c, ...
     ComRf_c,   ...
     ComRm_c,   ...
     sigX_c,    ...
     SigX_c,    ...
     E_L_c,     ...
     n_c,       ...
     v_c,       ...
     tc),  SCALED_tspan, SCALED_init_cond, SCALED_options); 



varNameStrings = ["Internal XIP","Free comR" ,"Mature TF comR", "mRNA", "SigX", "External XIP", "Luciferase"];
lgd_strings = ["XIP_{int}","comR","comR-XIP", "mRNA", "SigX", "XIP_{ext}", "E_L"];
color_strings = ["r","b","g","k", "c", "m", "r"];



SCALED_size_stat_vec = size(stat_vars_SCALED);
figure('Name','SCALED MODEL','NumberTitle','off')
for i = 1:1:SCALED_size_stat_vec(2)
    subplot(2,4,i)
    plot(t_SCALED,stat_vars_SCALED(:,i), strcat("-x",color_strings(i), ""));
    xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
        'FontSize', 15);
    title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);
    legend(lgd_strings(i));
end






%% Model for Bioluminscence resulting from mRNA 
% 
% Model captures the photon emission rate (implicitely by deriving the 
% photon count from the last state variable in the model below)


betaE = 0.1;    % kEl
deltaEp = 0.063;     
ks = 80;      
kp = 10; 
s_zero = 10;    



eta = 0.5;             % Quantum Efficient
%reacVol = n * 1e-15;  % Reaction Chamber Volume
N_Avo = 6.022e23;     % Avogadro Constant
RLUconst = 1e5;       % RLU constant


mRNA =stat_vars(:,end-2);



% Time scale for simulation
tspan1 = [0  tEnd];

options1 = odeset('NonNegative',[1,2,3],'RelTol',1e-5);

init_cond1 = [ 0.0;
               0.0;
               0.0  ];


[t1,stat_vars1] = ode15s(@(t1,stat_vars1) SMLuxLuminscence(t1, stat_vars1, t,...
    betaE, deltaEp, ks, kp, s_zero,n,v, mRNA)... 
    , tspan1, init_cond1, options1);



% Computing
% (A) Number of generated photons
% (B) Photon emission rate


% Formula as proposed in original paper
% NB! I am directly modeling the molecule count 
% N_ph = eta * N_Avo * reacVol * stat_vars1(:,end);

N_ph = eta * stat_vars1(:,end); 

% NOTE: This apparently does not work; needs to be fixed!!!
% Derivative 
dstat_vars_dt = zeros(size(t1,1), size(stat_vars1,2) );
for k = 1:size(t1,1)
dstat_vars_dt(k,:) = SMLuxLuminscence(t1(k), stat_vars1(k,:), t,...
    betaE, deltaEp, ks, kp, s_zero,n,v, mRNA);
end

dp_dt = dstat_vars_dt(:,end);

I_ph = eta * dp_dt;

RLU = RLUconst * I_ph;


% Plot Bioluminescence Measures

varNameStrings1 = ["Luciferase","Reductase" ,"Fatty Acid", "mRNA", "Photon Count"];
lgd_strings1 = ["E_L","E_P","P", "N_{Ph}"];
color_strings1 = ["k", "c", "r", "m"];


size_stat_vec1 = size(stat_vars1);

figure;
for i = 1:1:size_stat_vec1(2)
    subplot(2,3,i)
    plot(t1,stat_vars1(:,i), strcat("-x",color_strings1(i), ""));
    xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Molecule Number [molecules]', 'Interpreter', 'latex', ...
        'FontSize', 15);
    title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);
    legend(lgd_strings1(i));
end

subplot(2,3,4)
plot(t1,N_ph, strcat("-x","g ", ""), DisplayName='N_ph');
xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Photon Count', 'Interpreter', 'latex', ...
        'FontSize', 15);
title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);

subplot(2,3,5)
plot(t1,I_ph, strcat("-x","b", ""), DisplayName='I_ph');
xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$I_{Ph}$', 'Interpreter', 'latex', ...
        'FontSize', 15);
title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);



subplot(2,3,6)
plot(t1,RLU, strcat("-x","b", ""), DisplayName='RLU');
xlabel('Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('RLU', 'Interpreter', 'latex', ...
        'FontSize', 15);
title("QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
        'Interpreter', 'latex', 'FontSize', 15);




%% Parameter Fitting using WetLab Experimental Data


% Load WetLab Experimental Data
opts=detectImportOptions('../ExperimentalData/SM091 Luminescence 1st experiment.csv');
lum_data_01 = readtable('../ExperimentalData/SM091 Luminescence 1st experiment.csv'...
    ,opts, 'ReadVariableNames', true);

% Manually coded from Input CSV Table; dim = nM (TODO: Automatize)
dataXIP_In = [ 
    2000;
    2000;
    2000;
    1000;
    1000;
    1000;
    500;
    500;
    500;
    250;
    250;
    250;
    125;
    125;
    125;
    62.5;
    62.5;
    62.5;
    31.25;
    31.25;
    31.25;
    15.625;
    15.625;
    15.625;
    7.8;
    7.8;
    7.8;
    3.9;
    3.9;
    3.9;
    1.9;
    1.9;
    1.9        ];



% dataXIP_in scaled to AMOUNT[1]/VOLUME[um^3]
dataXIP_In_scaled =  (1e-9*N_Avo*1e-15)*dataXIP_In;


% Load Time Data from dataset
time_data = lum_data_01.Var1;

% Date String
time_conv_str = datestr(time_data, 'HH:MM');

% Number of days relative to January 1st XXXX (TODO: look up)
time_conv = datenum(time_conv_str,'HH:mm');

time_conv2 = datetime(time_conv_str,'InputFormat','HH:mm');

data_duration = diff([time_conv2(1) ,time_conv2(end)]);

data_duration_min = minutes(data_duration);


% Load Luminescence Data 
[dataset_size_row, dataset_size_cols] = size(lum_data_01);


% The data starts from column 3
data_set = zeros(dataset_size_row, dataset_size_cols-2);


for i=1:size(data_set,2)
    disp(i);
    dummy = lum_data_01(:, i + 2);
    data_set(:,i) = table2array(dummy);
end


% Prepare Optimization variable
r = optimvar('r',24, "LowerBound", 0.01, "UpperBound",1000);


testSolPts = RtoODE(r,...
     XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX, kx ,n ,v, ...
    tspan, init_cond,options, ...
    betaE, deltaEp, ks, kp, s_zero, ...
    tspan1, init_cond1, options1);


disp("Testing Optimization");


% Convert ODE Function to Optimization Expression
%myfcn = fcn2optimexpr(@RtoODE,r,tspan,init_cond);


% Create Objective Function
%obj = sum(sum((myfcn - yvalstrue).^2));



%optim_prob = optimproblem("Objective",obj);





% Parameter Vector Initial Guess


r0.r = [
            Df;
            delta_XIP_int;
            delta_XIP_int;
            k_TF_m;
            delta_tf ;
            k_TF_f;
            delta_mtf;
            alpha_comR;
            kappaB;
            kappa_comR;
            k_comR;
            alpha_SigX;
            kappaX;
            k_SigX;
            deltaSigX;
            deltaR;
            betaE;
            deltaEp;
            ks;
            kp;
            s_zero
       ].';



function solpts = RtoODE(r,...
     XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX, kx ,n ,v,...
    tspan, init_cond,options, ...
    betaE, deltaEp, ks, kp, s_zero, ...
    tspan1, init_cond1, options1)


        
        %sol = ode45(@(t,y)diffun(t,y,r),tspan,init_cond);
        
        %solpts = deval(sol,tspan);


    [t,stat_vars] = ode45(@(t,stat_vars) SM_XIP(t, stat_vars,XIP_ext  ...
        ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
        ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,kappaB ,kappa_comR ,k_comR ... 
        ,deltaR, kappaX, k_SigX,deltaSigX, kx ,n ,v,Ve, mu, Kn)... 
        , tspan, init_cond, options);
    
    
    mRNA =stat_vars(:,end-2);
    
    
    sol1 = ode45(@(t1,stat_vars1) SMLuxLuminscence(t1, stat_vars1, t,...
        betaE, deltaEp, ks, kp, s_zero, mRNA)... 
        , tspan1, init_cond1, options1);
    
    
    solpts = deval(sol1,tspan1);

end




% Model Definition - XIP Pathway Strepptococcus Mutans

% function dstat_vars_dt = SM_XIP(t, stat_vars,XIP_ext  ...
%     ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
%     ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
%     ,deltaR,  kappaX, k_SigX, deltaSigX ,n ,v, Ve, mu, Kn)
% 
% 
%         % XIP_int
%         dstat_vars_dt(1,1) = - Df *  n * (stat_vars(1)/(n*v) - stat_vars(6)) - delta_XIP_int*stat_vars(1)...
%                 - 1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) + k_TF_f* stat_vars(3); 
%         
%         
%         % Free comR
%         dstat_vars_dt(2,1) =  n * alpha_comR * kappaB  + k_TF_f* stat_vars(3)  - ...
%                               1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
%                               delta_tf * stat_vars(2)   + ...
%                               n * alpha_comR * kappaX * (((stat_vars(5))/(n*v))/ (k_SigX+((stat_vars(5))/(n*v))));   % SigX induced gene expression
%         
%         
%         
%         % Mature comR (XIP-comR Complex)
%         dstat_vars_dt(3,1) =  1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
%                               k_TF_f* stat_vars(3) - delta_mtf* stat_vars(3); 
%         
%         
%         % mRNA - sigX
%         dstat_vars_dt(4,1) =  n * kappa_comR *  (((stat_vars(3))/(n*v))/ (k_comR+((stat_vars(3))/(n*v)))) ...
%                 - deltaR * stat_vars(4); 
%         
%         
%         % SigX - Alternative Sigma Factor (Translated)
%         dstat_vars_dt(5,1) = n * kappa_comR * alpha_sigX *  (((stat_vars(3))/(n*v))/ (k_comR+((stat_vars(3))/(n*v)))) ...
%                             - deltaSigX * stat_vars(5);
%         
%         
%         
%         % XIP_ext
%         % NB! Note the addition of the factor 1/Ve that implies that a
%         % CONCENTRATION is modeled here
%         % TDODO: Check whether modeling this as concentration actually
%         % makes sense, for optimization purposes this should probably be
%         % modeled as AMOUNT (VERIFY!!!)
%         dstat_vars_dt(6,1) = Df *  n * 1/Ve * (stat_vars(1)/(n*v) -  stat_vars(6)) - delta_XIP_ext*stat_vars(6);
% 
% 
%         % Bacterial Growth
%         dstat_vars_dt(7,1) = mu * stat_vars(7)*(1-stat_vars(7)/Kn);
% 
% 
% end
%  
% 
% 
% 
% 
% function dstat_vars_dt = SMLuxLuminscence(t1, stat_vars1, t,...
%     betaE, deltaEp, ks, kp, s_zero, mRNA )
% 
%     % EL  
%     dstat_vars_dt(1,1) = betaE * interp1(t,mRNA,t1); 
%     
%     
%     % EP 
%     dstat_vars_dt(2,1) = betaE * interp1(t,mRNA,t1)- deltaEp * stat_vars1(2); 
%     
%     
%     % P : Product (Fatty Acid)
%     dstat_vars_dt(3,1)  = ks * s_zero * stat_vars1(1) - (ks * stat_vars1(1) + ...
%         kp * stat_vars1(2)) * stat_vars1(3);
% 
% end