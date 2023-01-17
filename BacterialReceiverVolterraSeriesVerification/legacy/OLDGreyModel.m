function [dx,y] = NLMod4GreyBox(t,x,u,  ...
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    ,deltaR, kappaX, k_SigX,deltaSigX, Ve,v,...
    betaE, deltaEp, ks, kp, s_zero, eta, RLUconst,deltaEL, varargin)

% NLMod4GreyBox Contains the Non-linear GreyBox Model Specification for the
% XIP Signaling Pathway in Strepptococcus mutans (Implementation following 
% https://se.mathworks.com/help/ident/ug/estimating-nonlinear-grey-box-models.html)
% 
% 
% x -> x0
%
% u -> Input NB! Time Series data which implies that input u = []; i.e., NO
% control input
% 
% NB! System does not have a classical input in the traditional sense
% 


params = [
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    ,deltaR, kappaX, k_SigX,deltaSigX, Ve,v,...
    betaE, deltaEL,deltaEp, ks, kp, s_zero, eta, RLUconst];

% DBG
global y_glob glob_i glob_stor glob_params glob_switch ...
t_glob stor_glob_t dirty_switch DBG





% Use Globally stored bacterial 
global n_sol_glob
n_sol = n_sol_glob;

try

x0 = x(1:7);
x1 = x(8:10);

% (1) Solving for mRNA

 % XIP_int
dx0(1,1) = - Df *  deval(n_sol,t) * (x0(1)/(deval(n_sol,t)*v) - x0(6)) - delta_XIP_int*x0(1)...
        - 1/(deval(n_sol,t)*v)  * k_TF_m * x0(1) * x0(2) + k_TF_f* x0(3); 


% Free comR
dx0(2,1) =  deval(n_sol,t) * alpha_comR * kappaB  + k_TF_f* x0(3)  - ...
                      1/(deval(n_sol,t)*v)  * k_TF_m * x0(1) * x0(2) - ...
                      delta_tf * x0(2)   + ...
                      deval(n_sol,t) * alpha_comR * kappaX * (((x0(5))/(deval(n_sol,t)*v))/ (k_SigX+((x0(5))/(deval(n_sol,t)*v))));   % SigX induced gene expression



% Mature comR (XIP-comR Complex) aka Activated Response Regulator
dx0(3,1) =  1/(deval(n_sol,t)*v)  * k_TF_m * x0(1) * x0(2) - ...
                      k_TF_f* x0(3) - delta_mtf* x0(3); 


% mRNA - sigX
dx0(4,1) =  deval(n_sol,t) * kappa_comR *  (((x0(3))/(deval(n_sol,t)*v))/ (k_comR+((x0(3))/(deval(n_sol,t)*v)))) ...
        - deltaR * x0(4); 


% SigX - Alternative Sigma Factor (Translated)
dx0(5,1) = deval(n_sol,t) * kappa_comR * alpha_sigX *  (((x0(3))/(deval(n_sol,t)*v))/ (k_comR+((x0(3))/(deval(n_sol,t)*v)))) ...
                    - deltaSigX * x0(5);



% XIP_ext - Concentration measures in units [1/um^3]
dx0(6,1) = Df *  deval(n_sol,t) * 1/Ve * (x0(1)/(deval(n_sol,t)*v) -  x0(6)) - delta_XIP_ext*x0(6);


% NEW - Added Bio-Luminescence Directly (Enzyme, i.e., Luciferase is assumed constant)
dx0(7,1) = betaE * x0(4) - x0(7) * deltaEL;


% EL  
dx1(1,1) = betaE * x0(4); 


% EP 
dx1(2,1) = betaE * x0(4)- deltaEp * x1(2); 


% P : Product (Fatty Acid)
dx1(3,1)  = ks * s_zero * x1(1) - (ks * x1(1)/(deval(n_sol,t)*v) + ...
kp * x1(2)/(deval(n_sol,t)*v)) * x1(3);



dx = [dx0;dx1];
%dx = dx0;

%dp_dt = dx(end);

I_ph = eta * x0(7);

RLU = (RLUconst*1e6) * I_ph;

y = RLU;

if DBG
    y_glob = [y_glob; y];
    t_glob = [t_glob; t];
    
    
    if((t == 900) && glob_switch)
           y_glob = [];
           t_glob = [];
    end
    
    
    if (t == 15 && glob_switch)
        glob_switch = false;
    end
    
    
    if (t == 900 && ~glob_switch)
        glob_switch = true;
    
        glob_stor = [glob_stor  [-1;y_glob]];
        stor_glob_t = [stor_glob_t [-1;t_glob]];
    
        % At FIRST (ONLY) iteration remove -1 
        if(~dirty_switch)
            dirty_switch = true;
            glob_stor = glob_stor(2:end);
            stor_glob_t = stor_glob_t(2:end);
        end
    
    
        glob_params = [glob_params; params];
        glob_i = glob_i +1;
        y_glob = [];
        t_glob = [];
    
    end
end

catch MExc

warning("OOps..Something went wrong");

end


end\



 %%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This files serves as solver for the SM XIP GreyBox Model
%
%
% F.Anderl - 16.11.2022



% Parameter Declaration & Definition - TODO: Implement "Load Parameters from CSV-file"


Df = 2.36e4;                    % Facilitated active diffusion transfer rate;
delta_XIP_int = 0.02;           % Internal XIP Degradation Rate
delta_XIP_ext = 0.02;           % External XIP Degradation Rate
k_TF_m = 0.01;                  % Transcription Factor Synthesis Rate
delta_tf = 0.02;                 % Transcription Factor Degradation Rate
k_TF_f = 0.1;                 % Transcription Factor Dematuration Rate
delta_mtf = 0.2;                % Mature Transcription Factor Degradation Rate
alpha_comR = 1;               % Translation Efficieny of sigX mRNA to ComR
kappaB = 0.1;                 % Basal RNA synthesis rate
kappa_comR =  0.01;
k_comR = 1;
alpha_sigX = 1;               % SigX Translation Efficiency
kappaX = 0.05;                 % Maximum SigX induced Transcription Rate
k_SigX = 1;
deltaSigX = 0.02;
deltaR = 0.02;                   % mRNA Degradation Rate
v = 1;                        % Bacteria Volume ~1e-15l or 1um^3
mu0 = 0.0090222;              % Bacterial Growth Constant
tau = 80.5504;                % Delay for Delay-Differential Equation
Kn = 7.5602;                  %1e8 (SCALE)  % Carrying Capacity (for bacterial gorowth function)
od600_coeff = 1.64e-9;
Ve = 2.2e11;                    % UNITS [um^3] Extracellular Volume; arbitrary 10ml; TODO  -> How to handle this?
tEnd = 900;                   % Time Scale of Experiment
n_init = 0.016 * (1/od600_coeff);      % Initial Bacterial Cell Count
%n_init = 1;                  % Initial Bacterial Cell Count
n = n_init;                   % Population Size; remnant albeit necessary
betaE = 0.05;                  % kEl
deltaEL = 0.00385;
deltaEp = 0.1;
ks = 0.01;
kp = 0.01;
s_zero = 200;


% Scaling Parameters for Bacterial Population
a = 1;
b = 1;

% Further Luminescence-related parameters
eta = 0.5;                     % Quantum Efficiency
RLUconst = 9.9394e-04;         %* 1e6 RLU constant

N_Avo = 6.022e23;

n_init_DDE  = n;

global n_sol_glob
% Solve for Bacterial Growth and store globally
n_sol_glob = dde23(@(tb,n_t,Z) BacGrowthDelayDE(tb, n_t, Z, mu0,Kn, a,b),...
    tau,@(tb) BacHist(tb,mu0,n_init_DDE),[0, 900]);



% Definition of added XIP supply as driving force of model
XIP_ext = 60;

% Steady-State comR-concentration (due to basal production)
comR_tilde = (n*alpha_comR*kappaB)/delta_tf;


% Load Experimental WetLab Data for Model Specifications
% Load WetLab Experimental Data
opts=detectImportOptions('../ExperimentalData/SM091 Luminescence 1st experiment.csv');
lum_data_01 = readtable('../ExperimentalData/SM091 Luminescence 1st experiment.csv'...
    ,opts, 'ReadVariableNames', true);

% Manually coded from Input CSV Table; dim = nM (TODO: Automatize)
dataXIP_In = [
    0;
    0;
    0;
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


% Read data into array
for i=1:size(data_set,2)
    dummy = lum_data_01(:, i + 2);
    data_set(:,i) = table2array(dummy);
end






% Define data which is to be investigated; `0 ->
inv_data = data_set(:,10);



sampling_times = linspace(0,data_duration_min,(size(inv_data,1)));


Ts = data_duration_min/(size(inv_data,1)-1);




% Setting up GreyBox Model


model_order = [1 0 10];

params = [
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    ,deltaR, kappaX, k_SigX,deltaSigX, Ve,v,...
    betaE ,  deltaEp, ks, kp, s_zero, eta, RLUconst, deltaEL];



init_states = [ 0;  % XIP_Int (A)
    comR_tilde;     % Free comR (A)
    0;              % comR-XIP_int Dimer (A)
    0;              % sigX_mRNA (A)
    0;              % SigX Protein/Product (A)
    XIP_ext;        % Extracellular XIP (C)
    0;              % Luciferase (New); TODO: Remove second lucerifase entry
    0;              % Luciferase (A)
    0;              % Reductase (A)
    0];             % Lux Reaction Product (A)




%model = idnlgrey('NLMod4GreyBox',model_order, params.', init_states, 0,'TimeUnit','minute');


% Interpolation from sparse data to get densly spaced virtual measurement
% data; this cannot be the end of the story though
t_pop = linspace(0, data_duration_min, 1000);
y_data_pop = interp1(sampling_times,inv_data,t_pop);


model_updated = idnlgrey('FullyIntegratedModel4GreyBox',model_order, params.', init_states, 0,'TimeUnit','minute');

model_updated = setinit(model_updated, 'Name', {'Internal XIP' 'Free comR' 'XIP-comR-Dimer' 'sigX mRNA' 'SigX' 'External XIP' 'Luciferase (NEW)' 'Luciferase' 'Reductase' 'Fatty Acid'});


% Define ODE Solver for GreyBox Model
model_updated.SimulationOptions.Solver = 'ode15s';



model_updated.InitialStates(1).Fixed =true;
model_updated.InitialStates(2).Fixed =true;
model_updated.InitialStates(3).Fixed =true;
model_updated.InitialStates(4).Fixed =true;
model_updated.InitialStates(5).Fixed =true;
model_updated.InitialStates(6).Fixed =true;
model_updated.InitialStates(7).Fixed =true;
model_updated.InitialStates(8).Fixed =true;
model_updated.InitialStates(9).Fixed =true;
model_updated.InitialStates(10).Fixed =true;

% Specify Parameers

model_updated.Parameters(1).Name = 'Transfer Rate';
model_updated.Parameters(2).Name = 'Int. XIP Degradation Constant';
model_updated.Parameters(3).Name = 'Ext. XIP Degradation Constant';
model_updated.Parameters(4).Name = 'Dimer Synthesis Rate Constant';
model_updated.Parameters(5).Name = 'Dimer Unbinding Rate Constant';
model_updated.Parameters(6).Name = 'comR Degradation Constant';
model_updated.Parameters(7).Name = 'comR-XIP-dimer Degradation Constant';
model_updated.Parameters(8).Name = 'ComR Translation Efficiency';
model_updated.Parameters(9).Name = 'SigX Translation Efficiency';
model_updated.Parameters(10).Name = 'Basal comR Transcription Rate';
model_updated.Parameters(11).Name = 'comR-activated sigX Transcription Rate';
model_updated.Parameters(12).Name = 'ComR MM Constant (ComR as TF)';
model_updated.Parameters(13).Name = 'sigX Degradation Constant';
model_updated.Parameters(14).Name = 'SigX-activated comR Transcription Rate';
model_updated.Parameters(15).Name = 'SigX MM Constant (SigX as TF)';
model_updated.Parameters(16).Name = 'SigX Degradation Constant';
model_updated.Parameters(17).Name = 'Test Volume MINUS Total Bacterial Population Volume';
model_updated.Parameters(18).Name = 'Bacterial Volume [um^3]';
model_updated.Parameters(19).Name = 'Luciferase/Reductase Synthesis Rate Constant';
model_updated.Parameters(20).Name = 'Reductase Degradation Constant';
model_updated.Parameters(21).Name = 'Product Formation Rate Constant (/w Substrate)';
model_updated.Parameters(22).Name = 'Product Turn-over (Reduction) Rate Constant';
model_updated.Parameters(23).Name = 'Initial Substrate (Luciferine) Supply';
model_updated.Parameters(24).Name = 'Reaction Quantum Yield';
model_updated.Parameters(25).Name = 'RLU Constant (LAB Equipment)';


model_updated.Parameters(1).Fixed = false; %"Transfer Rate"
model_updated.Parameters(2).Fixed = false; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Fixed = false; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Fixed = false; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Fixed = false; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Fixed = false; %"comR Degradation Constant"
model_updated.Parameters(7).Fixed = false; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Fixed = false; %"ComR Translation Efficiency"
model_updated.Parameters(9).Fixed = false; %"SigX Translation Efficiency"
model_updated.Parameters(10).Fixed = false; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Fixed = false; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Fixed = false; %"ComR MM Constant (ComR as TF)"
model_updated.Parameters(13).Fixed = false; %"sigX Degradation Constant"
model_updated.Parameters(14).Fixed = false; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(15).Fixed = false; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(16).Fixed = false; %"SigX Degradation Constant"
model_updated.Parameters(17).Fixed = true; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(18).Fixed = true; %" Bacterial Volume"
model_updated.Parameters(19).Fixed = false; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(20).Fixed = false; %"Reductase Degradation Constant"
model_updated.Parameters(21).Fixed = false; %"Product Formation Rate Constant (/w Substrate)"
model_updated.Parameters(22).Fixed = false; %"Product Turn-over (Reduction) Rate Constant"
model_updated.Parameters(23).Fixed = false; %"Initial Substrate (Luciferine) Supply"
model_updated.Parameters(24).Fixed = false; %"Reaction Quantum Yield"
model_updated.Parameters(25).Fixed = false; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(26).Fixed = false; %"RLU Constant (LAB Equipment)"


% Parameter Boundraries - Min
model_updated.Parameters(1).Minimum = 1e-22; %"Transfer Rate"
model_updated.Parameters(2).Minimum = 1e-22; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Minimum = 1e-22; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Minimum = 1e-22; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Minimum = 1e-22; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Minimum = 1e-22; %"comR Degradation Constant"
model_updated.Parameters(7).Minimum = 1e-22; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Minimum = 1e-22; %"ComR Translation Efficiency"
model_updated.Parameters(9).Minimum = 1e-22; %"SigX Translation Efficiency"
model_updated.Parameters(10).Minimum = 1e-22; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Minimum = 1e-22; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Minimum = 1e-22; %"ComR MM Constant (ComR as TF)"
model_updated.Parameters(13).Minimum = 1e-22; %"sigX Degradation Constant"
model_updated.Parameters(14).Minimum = 1e-22; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(15).Minimum = 1e-22; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(16).Minimum = 1e-22; %"SigX Degradation Constant"
model_updated.Parameters(17).Minimum = 1e-22; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(18).Minimum = 1; %"Bacterial Volume"
model_updated.Parameters(19).Minimum = 1e-22; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(20).Minimum = 1e-22; %"Reductase Degradation Constant"
model_updated.Parameters(21).Minimum = 1e-22; %"Product Formation Rate Constant (/w Substrate)"
model_updated.Parameters(22).Minimum = 1e-22; %"Product Turn-over (Reduction) Rate Constant"
model_updated.Parameters(23).Minimum = 1e-22; %"Initial Substrate (Luciferine) Supply"
model_updated.Parameters(24).Minimum = 0; %"Reaction Quantum Yield"
model_updated.Parameters(25).Minimum =1e-10; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(26).Minimum =1e-10; %"Luciferase Degradation Rate"


% Paramter Boundraries - Max


model_updated.Parameters(1).Maximum = 1e8; %"Transfer Rate"
model_updated.Parameters(2).Maximum = 10; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Maximum = 10; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Maximum = 10; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Maximum = 10; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Maximum = 10; %"comR Degradation Constant"
model_updated.Parameters(7).Maximum = 10; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Maximum = 10; %"ComR Translation Efficiency"
model_updated.Parameters(9).Maximum = 10; %"SigX Translation Efficiency"
model_updated.Parameters(10).Maximum = 10; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Maximum = 10; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Maximum = 10; %"ComR MM Constant (ComR as TF)"
model_updated.Parameters(13).Maximum = 10; %"sigX Degradation Constant"
model_updated.Parameters(14).Maximum = 10; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(15).Maximum = 10; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(16).Maximum = 10; %"SigX Degradation Constant"
model_updated.Parameters(17).Maximum = 1e20; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(18).Maximum = 1; %"Bacterial Volume"
model_updated.Parameters(19).Maximum = 10; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(20).Maximum = 10; %"Reductase Degradation Constant"
model_updated.Parameters(21).Maximum = 10; %"Product Formation Rate Constant (/w Substrate)"
model_updated.Parameters(22).Maximum = 10; %"Product Turn-over (Reduction) Rate Constant"
model_updated.Parameters(23).Maximum = 2000; %"Initial Substrate (Luciferine) Supply (Concentration)"
model_updated.Parameters(24).Maximum = 1; %"Reaction Quantum Yield"
model_updated.Parameters(25).Maximum = 1e2; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(26).Maximum = 1; %"Luciferase Degradation Rate"


% -------------------------------------------------------------------------------------------------

% NB!Excluding first element (t=0) so that ODE solvers in model do not throw
% an error

y_data = iddata(inv_data(2:end),[],[], 'SamplingInstances',(sampling_times(2:end)).', 'TimeUnit', 'minutes');

% Trying if densily sampled 'virtual' data works better than
% y_data = iddata(y_data_pop.',[],[], 'SamplingInstances',t_pop.', 'TimeUnit', 'minutes');


% Output Model Information
get(model_updated)

% Updated, fuly integrated model and plot
tic
figure;sim(model_updated,y_data);
toc

% Define Estimator Options
opt = nlgreyestOptions;
opt.Display = 'on';
opt.SearchOptions.Advanced.UseParallel = true;
opt.SearchMethod = 'fmincon';
opt.SearchOptions.Algorithm = 'trust-region-reflective';


clear global y_glob glob_i glob_stor glob_params
global y_glob glob_i glob_stor glob_params glob_switch ...
    t_glob stor_glob_t dirty_switch DBG

glob_i = 1;
dirty_switch = false;
glob_switch = false;
DBG = true;

try

% Perform NL Parameter Estimation
estm_model = nlgreyest(y_data,model_updated,opt);

catch MExc

warning("OOps..Something went wrong");

end




% GreyBox-Model based Dual Estimation (based on Extended Kalman Filter)

