% ODE Model Incorporating Bacterial Growth


function dstat_vars_dt = SM_XIP_DPS(t, stat_vars,XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR,  kappaX, k_SigX, deltaSigX,n_sol,v, Ve, deltaEL, alpha_L)



     

        % XIP_int
        dstat_vars_dt(1,1) = - Df *  deval(n_sol,t) * (stat_vars(1)/(deval(n_sol,t)*v) - stat_vars(6)) - delta_XIP_int*stat_vars(1)...
                - 1/(deval(n_sol,t)*v)  * k_TF_m * stat_vars(1) * stat_vars(2) + k_TF_f* stat_vars(3); 
   

        % Free comR
        dstat_vars_dt(2,1) =  deval(n_sol,t) * alpha_comR * kappaB  + k_TF_f* stat_vars(3)  - ...
                              1/(deval(n_sol,t)*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
                              delta_tf * stat_vars(2)   + ...
                              deval(n_sol,t) * alpha_comR * kappaX * (((stat_vars(5))/(deval(n_sol,t)*v))/ (k_SigX+((stat_vars(5))/(deval(n_sol,t)*v))));   % SigX induced gene expression
        
        
        
        % Mature comR (XIP-comR Complex) aka Activated Response Regulator
        dstat_vars_dt(3,1) =  1/(deval(n_sol,t)*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
                              k_TF_f* stat_vars(3) - delta_mtf* stat_vars(3); 
        
        


        % NEW: Model Luciferase directly /wo taking mRNA into account
        dstat_vars_dt(4,1) =  deval(n_sol,t) * alpha_L * kappa_comR *  (((stat_vars(3))/(deval(n_sol,t)*v))/ (k_comR+((stat_vars(3))/(deval(n_sol,t)*v)))) ...
                - deltaEL * stat_vars(4); 
        
        
        % SigX - Alternative Sigma Factor (Translated)
        dstat_vars_dt(5,1) = deval(n_sol,t) * kappa_comR * alpha_sigX *  (((stat_vars(3))/(deval(n_sol,t)*v))/ (k_comR+((stat_vars(3))/(deval(n_sol,t)*v)))) ...
                            - deltaSigX * stat_vars(5);
        
        
        
        % XIP_ext - Concentration measures in units [1/um^3]
        dstat_vars_dt(6,1) = Df *  deval(n_sol,t) * 1/Ve * (stat_vars(1)/(deval(n_sol,t)*v) -  stat_vars(6)) - delta_XIP_ext*stat_vars(6);

%         % Moved here: mRNA - sigX
%         dstat_vars_dt(7,1) =  deval(n_sol,t) * kappa_comR *  (((stat_vars(3))/(deval(n_sol,t)*v))/ (k_comR+((stat_vars(3))/(deval(n_sol,t)*v)))) ...
%                 - deltaR * stat_vars(7); 


%         % NEW/OLD - Added Bio-Luminescence Directly (Enzyme, i.e., Luciferase is assumed constant)
%         dstat_vars_dt(7,1) = betaE * stat_vars(4) - stat_vars(7) * deltaEL;


end
 