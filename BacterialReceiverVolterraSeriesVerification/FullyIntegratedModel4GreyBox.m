function [dx,y] = FullyIntegratedModel4GreyBox(t,x,u,  ...
    Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR ,alpha_sigX,kappaB ,kappa_comR ,k_comR ...
    ,deltaR, kappaX, k_SigX,deltaSigX, Ve,v,...
     betaE, eta, RLUconst,deltaEL, varargin)

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

try

dim_stoch_exp = varargin{1}{1};    % Parameter defined internally without handle to make adjustable/fittable (maybe change that in the long term) 

% Use Globally stored bacterial model
global n_sol_glob bac_gro_mode t_vec n_vec

if bac_gro_mode == 0
    n_sol = n_sol_glob;
elseif bac_gro_mode==1 
    n_sol = {t_vec, n_vec};
end

try EvalBacGrowth(t,bac_gro_mode,n_sol);

catch MExc
 
warning("Bacteria Solution Evaluation Error");

end

% Compute Bacterial Population size at t
bac_gro_at_t = EvalBacGrowth(t,bac_gro_mode,n_sol);


x0 = x(1:7);



% XIP_int
dx0(1,1) = - (Df*1e6) *  bac_gro_at_t * (x0(1)/(bac_gro_at_t*v) - x0(6)) - delta_XIP_int*x0(1)...
        - k_TF_m * (1/(bac_gro_at_t*v) *x0(1) * x0(2))^dim_stoch_exp + k_TF_f* x0(3); 


% Free comR
dx0(2,1) =  bac_gro_at_t * alpha_comR * kappaB  + k_TF_f* x0(3)  - ...
                      k_TF_m * (1/(bac_gro_at_t*v) *x0(1) * x0(2))^dim_stoch_exp - ...
                      delta_tf * x0(2)   + ...
                      bac_gro_at_t * alpha_comR * kappaX * (((x0(5))/(bac_gro_at_t*v))/ (k_SigX+((x0(5))/(bac_gro_at_t*v))));   % SigX induced gene expression



% Mature comR (XIP-comR Complex) aka Activated Response Regulator
dx0(3,1) =  1/2 * k_TF_m * (1/(bac_gro_at_t*v) *x0(1) * x0(2))^dim_stoch_exp - ...
                      k_TF_f* x0(3) - delta_mtf* x0(3); 


% mRNA - sigX
dx0(4,1) =  bac_gro_at_t * kappa_comR *  (((x0(3))/(bac_gro_at_t*v))/ (k_comR+((x0(3))/(bac_gro_at_t*v)))) ...
        - deltaR * x0(4); 
 


% SigX - Alternative Sigma Factor (Translated)
dx0(5,1) = bac_gro_at_t * kappa_comR * alpha_sigX *  (((x0(3))/(bac_gro_at_t*v))/ (k_comR+((x0(3))/(bac_gro_at_t*v)))) ...
                    - deltaSigX * x0(5);


% XIP_ext - Concentration measures in units [1/um^3]
dx0(6,1) = (Df*1e6) *  bac_gro_at_t * 1/(Ve*1e12) * (x0(1)/(bac_gro_at_t*v) -  x0(6)) - delta_XIP_ext*x0(6);


% NEW - Added Bio-Luminescence Directly (Enzyme, i.e., Luciferase is assumed constant)
dx0(7,1) = betaE * x0(4) - x0(7) * deltaEL;




dx = dx0;

I_ph = eta * x0(7);

RLU = (RLUconst) * I_ph;


y = RLU;

if (isnan(y) || isinf(y)) || ((sum(isnan(dx))>0) || (sum(isinf(dx))>0))

    warning("Y is Nan or Inf")
end


catch MExc

    warning("OOps..Something went wrong");

end



end