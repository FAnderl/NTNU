% ODE Model Incorporating Bacterial Growth


function dstat_vars_dt = SM_XIP_DPS_XTENDED(t, stat_vars,XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR,  kappaX, k_SigX, deltaSigX,n_sol,v, Ve, deltaEL, alpha_L,...
    dim_stoch_exp)


        % Use Globally stored bacterial model
        global n_sol_glob bac_gro_mode t_vec n_vec
        
        if bac_gro_mode == 0
            n_sol = n_sol;   % NO CHANGE
        elseif bac_gro_mode==1 
            n_sol = {t_vec, n_vec};
        end
        
        try EvalBacGrowth(t,bac_gro_mode,n_sol);
        
        catch MExc
         
        warning("Bacteria Solution Evaluation Error");
        
        end

        % Compute Bacterial Population size at t
        bac_gro_at_t = EvalBacGrowth(t,bac_gro_mode,n_sol);   % replaces deval(n_sol,t)

     
        if(any(stat_vars<0) )
            warning("Solution is Negative");
        end

        % Force Non-Negative Solution
        %stat_vars(stat_vars<0) = 0;

        disp(t)

        % XIP_int
        dstat_vars_dt(1,1) = - Df *  bac_gro_at_t * (stat_vars(1)/(bac_gro_at_t*v) - stat_vars(6)) - delta_XIP_int*stat_vars(1)...
                -   k_TF_m * (1/(bac_gro_at_t*v) *stat_vars(1) * stat_vars(2))^dim_stoch_exp + k_TF_f* stat_vars(3); 
   

        % Free comR
        dstat_vars_dt(2,1) =  bac_gro_at_t * alpha_comR * kappaB  + k_TF_f* stat_vars(3)  - ...
                               k_TF_m * (1/(bac_gro_at_t*v) *stat_vars(1) * stat_vars(2))^dim_stoch_exp - ...
                              delta_tf * stat_vars(2)   + ...
                              bac_gro_at_t * alpha_comR * kappaX * (((stat_vars(5))/(bac_gro_at_t*v))/ (k_SigX+((stat_vars(5))/(bac_gro_at_t*v))));   % SigX induced gene expression
        
        
        
        % Mature comR (XIP-comR Complex) aka Activated Response Regulator
        dstat_vars_dt(3,1) =   1/2 * k_TF_m * (1/(bac_gro_at_t*v) *stat_vars(1) * stat_vars(2))^dim_stoch_exp - ...
                              k_TF_f* stat_vars(3) - delta_mtf* stat_vars(3); 
        
        
        % NEW: Model Luciferase directly /wo taking mRNA into account
        dstat_vars_dt(4,1) =  bac_gro_at_t * alpha_L * kappa_comR *  (((stat_vars(3))/(bac_gro_at_t*v))/ (k_comR+((stat_vars(3))/(bac_gro_at_t*v)))) ...
                - deltaEL * stat_vars(4); 
        
        
        % SigX - Alternative Sigma Factor (Translated)
        dstat_vars_dt(5,1) = bac_gro_at_t * kappa_comR * alpha_sigX *  (((stat_vars(3))/(bac_gro_at_t*v))/ (k_comR+((stat_vars(3))/(bac_gro_at_t*v)))) ...
                            - deltaSigX * stat_vars(5);
        
        
        
        % XIP_ext - Concentration measures in units [1/um^3]
        dstat_vars_dt(6,1) = Df *  bac_gro_at_t * 1/Ve * (stat_vars(1)/(bac_gro_at_t*v) -  stat_vars(6)) - delta_XIP_ext*stat_vars(6);

%         % Moved here: mRNA - sigX
%         dstat_vars_dt(7,1) =  bac_gro_at_t * kappa_comR *  (((stat_vars(3))/(bac_gro_at_t*v))/ (k_comR+((stat_vars(3))/(bac_gro_at_t*v)))) ...
%                 - deltaR * stat_vars(7); 


%         % NEW/OLD - Added Bio-Luminescence Directly (Enzyme, i.e., Luciferase is assumed constant)
%         dstat_vars_dt(7,1) = betaE * stat_vars(4) - stat_vars(7) * deltaEL;




        


end
 