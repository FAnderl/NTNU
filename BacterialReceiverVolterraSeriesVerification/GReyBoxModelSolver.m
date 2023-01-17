
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file implements the wrapper for the solver of the SM XIP GreyBox Model
%
%
% F.Anderl - 16.11.2022


RED_SYS_FLAG = 1;



% Parameter Declaration & Definition - TODO: Implement "Load Parameters from CSV-file"


Df = 1e-6;                    % *1e6 ; Facilitated active diffusion transfer rate; ALTERNATIVELY: 0.54288 VS. 0.0236
delta_XIP_int = 0.07;         % Internal XIP Degradation Rate
delta_XIP_ext = 0.07;         % External XIP Degradation Rate
k_TF_m = 0.0091;              % Transcription Factor Synthesis Rate
delta_tf = 0.01;              % Transcription Factor Degradation Rate
k_TF_f = 1e-3;                % Transcription Factor Dematuration Rate    -> Likely a very small numerical value
delta_mtf = 0.01;             % Mature Transcription Factor Degradation Rate
alpha_comR = 1;               % Translation Efficieny of sigX mRNA to ComR
kappaB = 0.01;                % Basal RNA synthesis rate
kappa_comR =  9.575;
k_comR = 1.807;
alpha_sigX = 1;               % SigX Translation Efficiency
kappaX = 1.807;               % Maximum SigX induced Transcription Rate
k_SigX = 27.09;
deltaSigX = 0.2;
deltaR = 0.00152;             % mRNA Degradation Rate
v = 1;                        % Bacteria Volume ~1e-15l or 1um^3
mu0 = 0.0090222;              % Bacterial Growth Constant
tau = 80.5504;                % Delay for Delay-Differential Equation
Kn = 7.5602;                  % 1e8 (SCALE)  % Carrying Capacity (for bacterial gorowth function)
od600_coeff = 1.64e-9;
Ve = 0.22;                    % *1e12 -> UNITS [um^3] Extracellular Volume; arbitrary 10ml; TODO  -> How to handle this?
tEnd = 900;                   % Time Scale of Experiment
n_init = 0.016 * (1/od600_coeff);      % Initial Bacterial Cell Count
%n_init = 1;                  % Initial Bacterial Cell Count
n = n_init;                   % Population Size; remnant albeit necessary
betaE = 0.05;                 % kEl
deltaEL = 0.015;
alpha_L = 1;


% Scaling Parameters for Bacterial Population
a = 1;
b = 1;

% Further Luminescence-related parameters
eta = 0.98;                     % Quantum Efficiency
RLUconst = 1e-5;                % RLU constant

N_Avo = 6.022e23;

n_init_DDE  = n;

global n_sol_glob bac_gro_mode n_vec t_vec
% Solve for Bacterial Growth and store globally
n_sol_glob = dde23(@(tb,n_t,Z) BacGrowthDelayDE(tb, n_t, Z, mu0,Kn, a,b),...
    tau,@(tb) BacHist(tb,mu0,n_init_DDE),[0, 900]);



% NEW: 29.12.2022 - Use Experimental Bacterial Growth Curves for Estimation
load("bac_growth_data.mat");
n_vec = (1/od600_coeff) .* data_set;
t_vec = time_vec;



% Setting Bacterial Growth Mode (1: Use empirical data; 0: Use Bacterial Growth Model)
bac_gro_mode = 1;


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

data_cell_arr = cell(1,size(data_set,2));

for i=1:size(data_set,2)
    dummy = lum_data_01(:, i + 2);
    data_set(:,i) = table2array(dummy);
    data_cell_arr{i} = data_set(:,i); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation & Processing %%%%%%%%%%%%%%%%%%%%%%

% Select which data to use for estimation (Single Experiment Set-Up)
data_idx = 19;

% Cut corresponding slice from bacterial_growth experimental data
n_vec = n_vec(:,data_idx);


% Define data which is to be investigated; `0 ->
inv_data = data_set(:,data_idx);



% Create data_object that contains ALL experimental data



sampling_times = linspace(0,data_duration_min,(size(inv_data,1)));



Ts = data_duration_min/(size(inv_data,1)-1);




% Setting up GreyBox Model

if ~RED_SYS_FLAG


model_order = [1 0 7];


params = [
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    ,deltaR, kappaX, k_SigX,deltaSigX, Ve,v,...
    betaE , eta, RLUconst, deltaEL];

init_states = [ 0;  % XIP_Int (A)
    comR_tilde;     % Free comR (A)
    0;              % comR-XIP_int Dimer (A)
    0;              % Luciferase (A)
    0;              % SigX Protein/Product (A)
    dataXIP_In_scaled(data_idx);           % Extracellular XIP (C)
    0               % (mRNA)
     ];

elseif RED_SYS_FLAG

model_order = [1 0 6]; 

params = [
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    , kappaX, k_SigX,deltaSigX, Ve,v,...
     eta, RLUconst, deltaEL, alpha_L];


init_states = [ 0;  % XIP_Int (A)
    comR_tilde;     % Free comR (A)
    0;              % comR-XIP_int Dimer (A)
    0;              % Luciferase (A)
    0;              % SigX Protein/Product (A)
    dataXIP_In_scaled(data_idx);           % Extracellular XIP (C)
     ];


end

% Creating initital states matrix for all_experiments set-up

data_offset = 18;   % should always be data_idx-1 to maintain consistency

if ~RED_SYS_FLAG
    init_states_arr = zeros(7,size(data_set,2)-data_offset);
elseif RED_SYS_FLAG
    init_states_arr = zeros(6,size(data_set,2)-data_offset);
end

for i = 1 : (size(data_set,2) - data_offset) 
    idx = data_offset + i;


    if ~RED_SYS_FLAG 
    init_states_arr(:,i) = [ 0;    % XIP_Int (A)
    comR_tilde;                    % Free comR (A)
    0;                             % comR-XIP_int Dimer (A)
    0;                             % Luciferase (New)
    0;                             % SigX Protein/Product (A)
    dataXIP_In_scaled(idx);        % Extracellular XIP (C)
    0                              % sigX_mRNA (A)
     ];     

    elseif RED_SYS_FLAG

     init_states_arr(:,i) = [ 0;    % XIP_Int (A)
    comR_tilde;                    % Free comR (A)
    0;                             % comR-XIP_int Dimer (A)
    0;                             % Luciferase (New)
    0;                             % SigX Protein/Product (A)
    dataXIP_In_scaled(idx);        % Extracellular XIP (C)
     ];     


    end

end


% Populate 'sampling_times_cell_arr'
sampling_times_cell_arr = cell(1,size(data_set,2)-data_offset);
for i=1:(size(data_set,2)-data_offset)
    sampling_times_cell_arr{i} = sampling_times;
end






% Interpolation from sparse data to get densly spaced virtual measurement
% data; this cannot be the end of the story though

t_pop = linspace(0, data_duration_min, 1000);
y_data_pop = interp1(sampling_times,inv_data,t_pop);


% Model Specification

% NB! Init_states must be adapted to all_experiment vs single experiment
% set_up
if ~RED_SYS_FLAG
    model_updated = idnlgrey('FullyIntegratedModel4GreyBox',model_order, params.', init_states, 0,'TimeUnit','minute');
elseif RED_SYS_FLAG
    model_updated = idnlgrey('RED_FullyIntegratedModel4GreyBox',model_order, params.', init_states, 0,'TimeUnit','minute');
end 



if ~RED_SYS_FLAG
cell__name_arr = {'Internal XIP' 'Free comR' 'XIP-comR-Dimer' 'Luciferase (NEW)' 'SigX' 'External XIP' 'sigX mRNA' };
elseif RED_SYS_FLAG
cell__name_arr = {'Internal XIP' 'Free comR' 'XIP-comR-Dimer' 'Luciferase (NEW)' 'SigX' 'External XIP'  };
end

model_updated = setinit(model_updated, 'Name', cell__name_arr);

% Define ODE Solver for GreyBox Model
model_updated.SimulationOptions.Solver = 'ode15s';


% Non-Negative Options for ODE Solver




if ~RED_SYS_FLAG    % Including mRNA 

% Implies: DO NOT ESTIMATE INITIAL STATES
model_updated.InitialStates(1).Fixed =true;
model_updated.InitialStates(2).Fixed =true;
model_updated.InitialStates(3).Fixed =true;
model_updated.InitialStates(4).Fixed =true;
model_updated.InitialStates(5).Fixed =true;
model_updated.InitialStates(6).Fixed =true;
model_updated.InitialStates(7).Fixed =true;


% Specify Parameters
model_updated.Parameters(1).Name =  'qD - Transfer Rate';
model_updated.Parameters(2).Name =  'delta_XIP(int) - Int. XIP Degradation Constant';
model_updated.Parameters(3).Name =  'delta_XIP(ext) - Ext. XIP Degradation Constant';
model_updated.Parameters(4).Name =  'k_tf_m - Dimer Synthesis Rate Constant';
model_updated.Parameters(5).Name =  'k_tf_f - Dimer Unbinding Rate Constant';
model_updated.Parameters(6).Name =  'delta_tf - comR Degradation Constant';
model_updated.Parameters(7).Name =  'delta_mtf - comR-XIP-dimer Degradation Constant';
model_updated.Parameters(8).Name =  'alpha_comR - ComR Translation Efficiency';
model_updated.Parameters(9).Name =  'alpha_sigX - SigX Translation Efficiency';
model_updated.Parameters(10).Name = 'kappaB - Basal comR Transcription Rate';
model_updated.Parameters(11).Name = 'kappa_comR - comR-activated sigX Transcription Rate';
model_updated.Parameters(12).Name = 'k_comR - ComR MM Constant (ComR as TF)';
model_updated.Parameters(13).Name = 'delta_sigX - sigX Degradation Constant';
model_updated.Parameters(14).Name = 'kappa_SigX - SigX-activated comR Transcription Rate';
model_updated.Parameters(15).Name = 'k_SigX - SigX MM Constant (SigX as TF)';
model_updated.Parameters(16).Name = 'delta_SigX - SigX Degradation Constant';
model_updated.Parameters(17).Name = 'VE - Test Volume MINUS Total Bacterial Population Volume';
model_updated.Parameters(18).Name = 'V - Bacterial Volume [um^3]';
model_updated.Parameters(19).Name = 'k_L - Luciferase/Reductase Synthesis Rate Constant';
model_updated.Parameters(20).Name = 'eta - Reaction Quantum Yield';
model_updated.Parameters(21).Name = 'RLUconst - RLU Constant (LAB Equipment)';
model_updated.Parameters(22).Name = 'deltaEL - Luciferase Degradation Constant';


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
model_updated.Parameters(20).Fixed = true; %"Reaction Quantum Yield"
model_updated.Parameters(21).Fixed = false; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(22).Fixed = false; %"Luciferase Degradation Constant"


% Parameter Boundraries - Min
model_updated.Parameters(1).Minimum = 1e-10; %"Transfer Rate"
model_updated.Parameters(2).Minimum = 1e-5; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Minimum = 1e-5; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Minimum = 1e-5; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Minimum = 1e-9; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Minimum = 1e-5; %"comR Degradation Constant"
model_updated.Parameters(7).Minimum = 1e-5; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Minimum = 1e-5; %"ComR Translation Efficiency"
model_updated.Parameters(9).Minimum = 1e-5; %"SigX Translation Efficiency"
model_updated.Parameters(10).Minimum = 1e-5; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Minimum = 1e-10; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Minimum = 1e-5; %"ComR MM Constant (ComR as TF)"
model_updated.Parameters(13).Minimum = 1e-5; %"sigX Degradation Constant"
model_updated.Parameters(14).Minimum = 1e-5; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(15).Minimum = 1e-5; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(16).Minimum = 1e-5; %"SigX Degradation Constant"
model_updated.Parameters(17).Minimum = 1e-5; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(18).Minimum = 1; %"Bacterial Volume"
model_updated.Parameters(19).Minimum = 1e-5; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(20).Minimum = 0.001; %"Reaction Quantum Yield"
model_updated.Parameters(21).Minimum =1e-10; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(22).Minimum =1e-5; %"Luciferase Degradation Rate"


% Paramter Boundraries - Max


model_updated.Parameters(1).Maximum = 1; %"Transfer Rate"
model_updated.Parameters(2).Maximum = 1; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Maximum = 1; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Maximum = 1; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Maximum = 1; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Maximum = 1; %"comR Degradation Constant"
model_updated.Parameters(7).Maximum = 1; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Maximum = 1; %"ComR Translation Efficiency"
model_updated.Parameters(9).Maximum = 1; %"SigX Translation Efficiency"
model_updated.Parameters(10).Maximum = 1; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Maximum = 100; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Maximum = 100; %"ComR MM Constant (ComR as TF)"
model_updated.Parameters(13).Maximum = 1; %"sigX Degradation Constant"
model_updated.Parameters(14).Maximum = 100; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(15).Maximum = 100; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(16).Maximum = 1; %"SigX Degradation Constant"
model_updated.Parameters(17).Maximum = 100; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(18).Maximum = 1; %"Bacterial Volume"
model_updated.Parameters(19).Maximum = 1; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(20).Maximum = 1; %"Reaction Quantum Yield"
model_updated.Parameters(21).Maximum = 1; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(22).Maximum = 1; %"Luciferase Degradation Rate"



elseif RED_SYS_FLAG

% Implies: DO NOT ESTIMATE INITIAL STATES
model_updated.InitialStates(1).Fixed =true;
model_updated.InitialStates(2).Fixed =true;
model_updated.InitialStates(3).Fixed =true;
model_updated.InitialStates(4).Fixed =true;
model_updated.InitialStates(5).Fixed =true;
model_updated.InitialStates(6).Fixed =true;




% Specify Parameters
model_updated.Parameters(1).Name =  'qD - Transfer Rate';
model_updated.Parameters(2).Name =  'delta_XIP(int) - Int. XIP Degradation Constant';
model_updated.Parameters(3).Name =  'delta_XIP(ext) - Ext. XIP Degradation Constant';
model_updated.Parameters(4).Name =  'k_tf_m - Dimer Synthesis Rate Constant';
model_updated.Parameters(5).Name =  'k_tf_f - Dimer Unbinding Rate Constant';
model_updated.Parameters(6).Name =  'delta_tf - comR Degradation Constant';
model_updated.Parameters(7).Name =  'delta_mtf - comR-XIP-dimer Degradation Constant';
model_updated.Parameters(8).Name =  'alpha_comR - ComR Translation Efficiency';
model_updated.Parameters(9).Name =  'alpha_sigX - SigX Translation Efficiency';
model_updated.Parameters(10).Name = 'kappaB - Basal comR Transcription Rate';
model_updated.Parameters(11).Name = 'kappa_comR - comR-activated sigX Transcription Rate';
model_updated.Parameters(12).Name = 'k_comR - ComR MM Constant (ComR as TF)';
%model_updated.Parameters(13).Name = 'delta_sigX - sigX Degradation Constant';
model_updated.Parameters(13).Name = 'kappa_SigX - SigX-activated comR Transcription Rate';
model_updated.Parameters(14).Name = 'k_SigX - SigX MM Constant (SigX as TF)';
model_updated.Parameters(15).Name = 'delta_SigX - SigX Degradation Constant';
model_updated.Parameters(16).Name = 'VE - Test Volume MINUS Total Bacterial Population Volume';
model_updated.Parameters(17).Name = 'V - Bacterial Volume [um^3]';
%model_updated.Parameters(19).Name = 'k_L - Luciferase/Reductase Synthesis Rate Constant';
model_updated.Parameters(18).Name = 'eta - Reaction Quantum Yield';
model_updated.Parameters(19).Name = 'RLUconst - RLU Constant (LAB Equipment)';
model_updated.Parameters(20).Name = 'deltaEL - Luciferase Degradation Constant';
model_updated.Parameters(21).Name = 'alphaL - luxAB Translation Efficieny';

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
%model_updated.Parameters(13).Fixed = false; %"sigX Degradation Constant"
model_updated.Parameters(13).Fixed = false; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(14).Fixed = false; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(15).Fixed = false; %"SigX Degradation Constant"
model_updated.Parameters(16).Fixed = true; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(17).Fixed = true; %" Bacterial Volume"
%model_updated.Parameters(19).Fixed = false; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(18).Fixed = true; %"Reaction Quantum Yield"
model_updated.Parameters(19).Fixed = false; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(20).Fixed = false; %"Luciferase Degradation Constant"
model_updated.Parameters(21).Fixed =  false; %'alphaL - luxAB Translation Efficieny';


% Parameter Boundraries - Min
model_updated.Parameters(1).Minimum = 1e-10; %"Transfer Rate"
model_updated.Parameters(2).Minimum = 1e-5; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Minimum = 1e-5; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Minimum = 1e-5; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Minimum = 1e-9; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Minimum = 1e-5; %"comR Degradation Constant"
model_updated.Parameters(7).Minimum = 1e-5; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Minimum = 1e-5; %"ComR Translation Efficiency"
model_updated.Parameters(9).Minimum = 1e-5; %"SigX Translation Efficiency"
model_updated.Parameters(10).Minimum = 1e-5; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Minimum = 1e-5; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Minimum = 1e-5; %"ComR MM Constant (ComR as TF)"
%model_updated.Parameters(13).Minimum = 1e-5; %"sigX Degradation Constant"
model_updated.Parameters(13).Minimum = 1e-5; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(14).Minimum = 1e-5; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(15).Minimum = 1e-5; %"SigX Degradation Constant"
model_updated.Parameters(16).Minimum = 1e-5; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(17).Minimum = 1; %"Bacterial Volume"
%model_updated.Parameters(19).Minimum = 1e-5; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(18).Minimum = 0.001; %"Reaction Quantum Yield"
model_updated.Parameters(19).Minimum =1e-10; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(20).Minimum =1e-5; %"Luciferase Degradation Rate"
model_updated.Parameters(21).Minimum = 1e-5; %'alphaL - luxAB Translation Efficieny';



% Paramter Boundraries - Max
model_updated.Parameters(1).Maximum = 1; %"Transfer Rate"
model_updated.Parameters(2).Maximum = 1; %"Int. XIP Degradation Constant"
model_updated.Parameters(3).Maximum = 1; %"Ext. XIP Degradation Constant"
model_updated.Parameters(4).Maximum = 1; %"Dimer Synthesis Rate Constant"
model_updated.Parameters(5).Maximum = 1; %"Dimer Unbinding Rate Constant"
model_updated.Parameters(6).Maximum = 1; %"comR Degradation Constant"
model_updated.Parameters(7).Maximum = 1; %"comR-XIP-dimer Degradation Constant"
model_updated.Parameters(8).Maximum = 1; %"ComR Translation Efficiency"
model_updated.Parameters(9).Maximum = 1; %"SigX Translation Efficiency"
model_updated.Parameters(10).Maximum = 1; %"Basal comR Transcription Rate"
model_updated.Parameters(11).Maximum = 100; %"comR-activated sigX Transcription Rate"
model_updated.Parameters(12).Maximum = 1e3; %"ComR MM Constant (ComR as TF)"
%model_updated.Parameters(13).Maximum = 1; %"sigX Degradation Constant"
model_updated.Parameters(13).Maximum = 100; %"SigX-activated comR Transcription Rate"
model_updated.Parameters(14).Maximum = 1e3; %"SigX MM Constant (SigX as TF)"
model_updated.Parameters(15).Maximum = 1; %"SigX Degradation Constant"
model_updated.Parameters(16).Maximum = 100; %"Test Volume MINUS Total Bacterial Population Volume"
model_updated.Parameters(17).Maximum = 1; %"Bacterial Volume"
%model_updated.Parameters(19).Maximum = 1; %"Luciferase/Reductase Synthesis Rate Constant"
model_updated.Parameters(18).Maximum = 1; %"Reaction Quantum Yield"
model_updated.Parameters(19).Maximum = 1; %"RLU Constant (LAB Equipment)"
model_updated.Parameters(20).Maximum = 1; %"Luciferase Degradation Rate"
model_updated.Parameters(21).Maximum = 1; %'alphaL - luxAB Translation Efficieny';


end



% -------------------------------------------------------------------------

% NB!Excluding first element (t=0) so that ODE solvers in model do not throw
% an error

y_data = iddata(inv_data(2:end),[],[], 'SamplingInstances',(sampling_times(2:end)).', 'TimeUnit', 'minutes');

% Data object when optimizing over all available experiments
data_cell_arr = data_cell_arr((data_offset+1):end);
%y_data = iddata(data_cell_arr,[],[], 'SamplingInstances',sampling_times_cell_arr, 'TimeUnit', 'minutes');



%DBG
% Set initial parameter values to Max/Min  Valuesto test feasability of
% bounds
% for i = 1:size(params,2)
%     model_updated.Parameters(i).Value = model_updated.Parameters(i).Minimum;
% end

% Output Model Information
get(model_updated)



%% Pre Estimation Simulation
% Simulate Model
preEstm_model_pred = sim(model_updated,y_data);
preEstm_model_pred_data = preEstm_model_pred.OutputData;

pre_estm_fig = figure('Name','Pre Estimation Fit','NumberTitle','off');
tg = uitabgroup(pre_estm_fig); % tabgroup


if size(preEstm_model_pred_data,2) > 1

for i = 1 : 1:(size(data_set,2) - data_offset)
    idx = data_offset + i;
    thistab = uitab(tg, "Title",preEstm_model_pred.ExperimentName{i});
    axes('Parent',thistab);
    plot(sampling_times_cell_arr{i}, data_cell_arr{i}, "o");
    hold on
    plot(sampling_times_cell_arr{i}, preEstm_model_pred_data{i});

end


elseif size(preEstm_model_pred_data,2) == 1
    figure;
    plot(sampling_times, inv_data, "o");
    hold on
    plot(sampling_times(2:end), preEstm_model_pred_data);


end

% TODO: Implement Pre-Estimation Fit Plot for Single Experiment Set-Up









%% Estimation
try

% Define Estimator Options
opt = nlgreyestOptions;
opt.Display = 'on';
opt.SearchOptions.Advanced.UseParallel = true;
opt.SearchMethod = 'fmincon';
opt.SearchOptions.Algorithm = 'active-set';   % sqp, trust-region-reflective,interior-point                                    
%opt.SearchOptions.Advanced.FinDiffRelStep = 1e-1;    

% Perform NL Parameter Estimation
estm_model = nlgreyest(y_data,model_updated,opt);

catch MExc

warning("OOps..Something went wrong -> " + MExc.message);
disp(getReport(MExc,'extended'));
end

%% Simulate Estimated Model
tic
figure;sim(estm_model,y_data);
hold on
toc



%% Plot Estimated Model against Experimental Data

model_pred = sim(estm_model,y_data);
model_pred_data = model_pred.OutputData;

fig = figure('Name','Post Estimation Fit','NumberTitle','off');
tg = uitabgroup(fig); % tabgroup

if size(model_pred_data,2) > 1

for i = 1 : 1:(size(data_set,2) - data_offset)
    idx = data_offset + i;
    thistab = uitab(tg, "Title",model_pred.ExperimentName{i});
    axes('Parent',thistab);
    plot(sampling_times_cell_arr{i}, data_cell_arr{i}, "o");
    hold on
    plot(sampling_times_cell_arr{i}, model_pred_data{i});

end

elseif size(model_pred_data,2) == 1
    figure;
    plot(sampling_times, inv_data, "o");
    hold on
    plot(sampling_times(2:end), model_pred_data);

end

disp("dbg");



%% Further Analysis
% Load Specific Estiamted Model for closer Analysis
if ~exist("estm_model","var")
    load("estm_model_FULL_EXPERIMENT_ALG_active_set.mat");
end
