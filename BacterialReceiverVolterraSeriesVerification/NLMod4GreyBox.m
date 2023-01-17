function [dx,y] = NLMod4GreyBox(t,x,u,  ...
    XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX,n ,v, Ve,...
    betaE, deltaEp, ks, kp, s_zero, mu0, Kn, a, b, tau, eta, RLUconst, varargin)

% NLMod4GreyBox Contains the Non-linear GreyBox Model Specification for the
% XIP Signaling Pathway in Strepptococcus mutans (Implementation following 
% https://se.mathworks.com/help/ident/ug/estimating-nonlinear-grey-box-models.html)
% 
% 
% x -> stat_vars
%
% u -> Input NB! Time Series data which implies that input u = []; i.e., NO
% control input
% 
% NB! System does not have a classical input in the traditional sense
% 

%disp(t)

% Global variables
global t_minus_one prev_t


if abs(t-prev_t) > 1e-5 
    t_minus_one = prev_t;
    disp(t_minus_one);
end


% Computation of Model Outputs 



% (0) - Bacterial Growth - Delay-Differential Equation Model
dn_dt = zeros(1,1);

n_init_DDE  = n;


% Solve Delay Differential Equation for Heuristic Fit up to 't'
n_sol = dde23(@(tb,n_t,Z) BacGrowthDelayDE(tb, n_t, Z, mu0,Kn, a,b),...
        tau,@(tb) BacHist(tb,mu0,n_init_DDE),[0, 900]);
 
n = n_sol;


% (1) Solving for mRNA





options = odeset('NonNegative',[1,2,3,4,5,6]);
% init_cond = [ 0.0;
%           comR_tilde;
%           0.0;
%           0.0;
%           0.0;
%           XIP_Sig_Num];

% NOTE: Needs to be only solved to supply the dx1_dt equation system with
% the required input (mRNA )
[t0,x0] = ode45(@(t0,x0) SM_XIP_DPS(t0, x0,XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX,n ,v, Ve)... 
    , [t_minus_one t], x(1:6), options);



% (2) Solving for Luminescence 

% betaE = 0.5;      % kEl
% deltaEp = 1;     
% ks = 15;        
% kp = 1;         
% s_zero = 1;    





% Part of Solution of Subsystem 1 is used as input in Subsystem 2
mRNA =x0(:,end-2);


options1 = odeset('NonNegative',[1,2,3],'RelTol',1e-5);



% Solving this ODE system is NOT strictly necessary
[t1,x1] = ode45(@(t1,x1) SMLuxLuminscence(t1, x1, t0,...
    betaE, deltaEp, ks, kp, s_zero,n,v, mRNA)... 
    , [t_minus_one t], x(7:end), options1);


% x - State Variables 

%x = [x0;x1];



% dx - Computation of State-Variable Derivatives


dx0 = SM_XIP_DPS(t, x(1:6),XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR, kappaX, k_SigX,deltaSigX,n ,v, Ve);



dx1 = SMLuxLuminscence(t, x(7:end), t0,...
    betaE, deltaEp, ks, kp, s_zero,n,v, mRNA);



dx = [dx0;dx1];


dp_dt = dx(end);

I_ph = eta * dp_dt;

RLU = RLUconst * I_ph;

y = RLU;


%disp(dx)

% Update Previous time_step
prev_t = t; 



end