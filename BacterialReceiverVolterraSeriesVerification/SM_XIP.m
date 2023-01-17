function dstat_vars_dt = SM_XIP(t, stat_vars,XIP_ext  ...
    ,Df, delta_XIP_int, delta_XIP_ext ,k_TF_m  ...
    ,k_TF_f ,delta_tf ,delta_mtf ,alpha_comR,alpha_sigX,kappaB ,kappa_comR ,k_comR ... 
    ,deltaR,  kappaX, k_SigX, deltaSigX ,n ,v, Ve)


        % XIP_int
        dstat_vars_dt(1,1) = - Df *  n * (stat_vars(1)/(n*v) - stat_vars(6)) - delta_XIP_int*stat_vars(1)...
                - 1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) + k_TF_f* stat_vars(3); 
        
       

        % Free comR
        dstat_vars_dt(2,1) =  n * alpha_comR * kappaB  + k_TF_f* stat_vars(3)  - ...
                              1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
                              delta_tf * stat_vars(2)   + ...
                              n * alpha_comR * kappaX * (((stat_vars(5))/(n*v))/ (k_SigX+((stat_vars(5))/(n*v))));   % SigX induced gene expression
        
        
        
        % Mature comR (XIP-comR Complex)
        dstat_vars_dt(3,1) =  1/(n*v)  * k_TF_m * stat_vars(1) * stat_vars(2) - ...
                              k_TF_f* stat_vars(3) - delta_mtf* stat_vars(3); 
        
        
        % mRNA - sigX
        dstat_vars_dt(4,1) =  n * kappa_comR *  (((stat_vars(3))/(n*v))/ (k_comR+((stat_vars(3))/(n*v)))) ...
                - deltaR * stat_vars(4); 
        
        
        % SigX - Alternative Sigma Factor (Translated)
        dstat_vars_dt(5,1) = n * kappa_comR * alpha_sigX *  (((stat_vars(3))/(n*v))/ (k_comR+((stat_vars(3))/(n*v)))) ...
                            - deltaSigX * stat_vars(5);
        
        
        
        % XIP_ext
        % NB! Note the addition of the factor 1/Ve that implies that a
        % CONCENTRATION is modelled here
        % TDODO: Check whether modeling this as concentration actually
        % makes sense, for optimization purposes this should probably be
        % modeled as AMOUNT (VERIFY!!!)
        dstat_vars_dt(6,1) = Df *  n * 1/Ve * (stat_vars(1)/(n*v) -  stat_vars(6)) - delta_XIP_ext*stat_vars(6);


end
 