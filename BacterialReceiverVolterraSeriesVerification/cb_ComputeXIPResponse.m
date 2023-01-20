function cb_ComputeXIPResponse(app,...
qD,...
delta_XIP ,...
k_TF_m ,...
delta_tf ,...
k_TF_f ,...
delta_mtf ,...
alpha_comR ,...
kappaB ,...
kappa_comR ,...
k_comR ,...
alpha_sigX ,...
kappaX ,...
k_SigX ,...
deltaSigX ,...
deltaR ,...
alpha_L ,...
deltaEL, ...
comR_tilde,...
XIP_input)

XTENDED_MODEL_FLAG = 1;   % USE MODEL CONSIDERING STOCHIOMETRY OF DIMERIZATION (1)


%cb_ComputeXIPResponse Computes the response of the XIP system with
%parameters set in parent app
%   Detailed explanation goes here

% Df = 0.54288*1e-2;               % Facilitated active diffusion constant / transfer rate / conductivity
% qD = Df;                         % Renaming Routine
% delta_XIP = 0.00035;
 delta_XIP_int = delta_XIP;        % Internal XIP Degradation Rate Constant
 delta_XIP_ext = delta_XIP;        % External XIP Degradation Rate Constant
% k_TF_m = 1e-3;                   % Transcription Factor Synthesis Rate Constant
% delta_tf = 0.01;                 % Free Transcription Factor (ComR) Degradation Rate Constant
% k_TF_f = 0.00000001;             % Transcription Factor Dematuration Rate Constant; has to be really, really smal
% delta_mtf = 1e3;                 % Mature Transcription Factor Degradation Rate 
% alpha_comR = 0.5;                % Translation Efficieny of comR mRNA to ComR
% kappaB = 1e-2;                   % Basal RNA synthesis rate
% kappa_comR =  200;               % Maximum ComR-XIP_int-complex-induced sigX transcription rate
% k_comR = 0.001;                  % Michaelis-Menten Constant
% alpha_sigX = 1;                  % SigX Translation Efficiency 
% kappaX = 1000;                   % Maximum SigX induced comR Transcription Rate 
% k_SigX = 0.1;                    % 0.5 (TRANSIENT EXPLANATION); Michaelis-Menten Constant
% deltaSigX = 10;                  % 10 (TRANSIENT EXPLANATION)
% deltaR = 0.02;                   % mRNA Degradation Rate 
% kL = 0.05;                       % kEl: Luciferase Translation Rate
% deltaEL = 0.00385;               % Luciferase Decay Rate Constant 

v = 1;                             % Bacteria Volume ~1e-15l or 1um^3

% Fitted Parameters for Bacterial Growth @ 62.6nM
% mu: 0.0082041
% Kn: 7.6964
% a: 1
% b: 1
% tau: 83.3785

mu0 = 0.0082041;%0.0090222;                % Bacterial Growth Constant           
Kn = .6964;%7.5602;                    % 1e8 (SCALE)  % Carrying Capacity (for bacterial gorowth function)

od600_coeff = 1.64e-9;

Ve = 2.2e11;                    % UNITS [um^3] Extracellular Volume; arbitrary 10ml; TODO  -> How to handle this?

tEnd = 900;                     % Time Scale of Experiment


n_init = 0.016  * (1/od600_coeff);      % Initial Bacterial Cell Count DEFAULT:0.016 
%n_init = 1;                    % Initial Bacterial Cell Count

n = n_init;                     % Population Size; remnant albeit necessary 


% Steady-State comR-concentration (due to basal production)
%comR_tilde = (n*alpha_comR*kappaB)/delta_tf;

% Imposed ZERO initital condition on comR concentration 
%comR_tilde = 0;


% comR_tilde overwritten from App GUI

N_Avo = 6.022e23;


dim_stoch_exp = 2.5;


% Extracellular XIP - Driving CONCENTRATION
XIP_Sig_Num = 60;  % 100nM, i.e., ~60 molecules per bacterial cell volume
XIP_Sig_Num = XIP_input*(1e-9*N_Avo*1e-15);  % Overwritten by external parameter

%% Model for Gene Expression


% Incorporating Bacterial Growth




global n_sol_glob bac_gro_mode n_vec t_vec


%  Use Experimental Bacterial Growth Curves for Estimation (Note: Copied from )
load("bac_growth_data.mat");
n_vec = (1/od600_coeff) .* data_set;
t_vec = time_vec;
% Select which data to use for estimation (Single Experiment Set-Up)
data_idx = 19;
% Cut corresponding slice from bacterial_growth experimental data
n_vec = n_vec(:,data_idx);



% Setting Bacterial Growth Mode (1: Use empirical data; 0: Use Bacterial Growth Model)
bac_gro_mode = 1;



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
tau = 83.3785;% 80.5504;
n_init_DDE  = n_init;

% Solve for Bacterial Growth
%[t,n_t] = ode45(@(t,n_t) BacterialGrowth(t,n_t,mu,Kn, od600_coeff,), tspan, n_init);


% Solve Delay Differential Equation for Heuristic Fit 
sol = dde23(@(t,n_t,Z) BacGrowthDelayDE(t, n_t, Z, mu0,Kn, a,b),...
        tau,@(t) BacHist(t,mu0,n_init_DDE),[0, 900]);
n_sol_glob = sol;

% figure(1);
% plot(app.UIAxes,sol.x,sol.y,'-o', 'LineWidth',3.0,'Color','m')
% title('Bacterial Growth - Modelled as Delay Differential Equation','FontSize',15);
% xlabel("Time",'FontSize',15);
% ylabel("Population size [Number of Cells]",'FontSize',15)

% Use Dynamic Bacterial Growth Model
if DPS_FLAG
  XIP_model = @SM_XIP_DPS;

if XTENDED_MODEL_FLAG 
    XIP_model = @SM_XIP_DPS_XTENDED;
end

  n = sol; % n_sol

elseif ~DPS_FLAG
  XIP_model = @SM_XIP;

end


% ---------- Solve for mRNA -----------

% Initialize State Variable Vector
dstat_vars_dt= zeros(MODEL_SZ,1);



[t,stat_vars] = ode15s(@(t,stat_vars) XIP_model(t, stat_vars,XIP_ext  ...
    ,qD, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX,n ,v, Ve, deltaEL, alpha_L,dim_stoch_exp)... 
    , tspan, init_cond, options);




varNameStrings = ["Internal XIP","Free comR" ,"Mature TF comR", "Luciferase", "SigX", "External XIP", "(mRNA)"];
lgd_strings = ["XIP_{int}","comR","comR-XIP", "Luciferase", "SigX", "XIP_{ext}", "(mRNA)"];
color_strings = ["r","b","g","k", "c", "m", "r"];


size_stat_vec = size(stat_vars);

%figure(2);
axesLst = [
           app.XIP_intUIAxes,
           app.comRUIAxes,
           app.comR_XIPUIAxes,
           app.LuciferaseUIAxes,
           app.SigXUIAxes,
           app.XIP_extUIAxes,
           app.mRNAUIAxes];
for i = 1:1:size_stat_vec(2)
    %subplot(2,4,i)
    plot(axesLst(i),t,stat_vars(:,i), strcat("-x",color_strings(i), ""));
%     xlabel(axesLst(i),'Time [Unit Time]', 'Interpreter', 'latex', 'FontSize', 15);
%     ylabel(axesLst(i),'Molecule Number [molecules]', 'Interpreter', 'latex', ...
%         'FontSize', 15);
%     title(axesLst(i),"QS Signal: " + num2str(XIP_Sig_Num) +" [1/Bacteria Volume]", ...
%         'Interpreter', 'latex', 'FontSize', 15);
    legend(axesLst(i),lgd_strings(i));
end



end