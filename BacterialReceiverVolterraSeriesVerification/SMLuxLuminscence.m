
function dstat_vars_dt = simple_SMLuxLuminscence(t1, stat_vars1, t,...
    betaE, deltaEp, ks, kp, s_zero, n,v, mRNA )

    interp_mRNA = interp1(t,mRNA,t1);

    
    % EL  
    dstat_vars_dt(1,1) = betaE * interp_mRNA; 
    
    
    % Luminescence   
   dstat_vars_dt(2,1) = betaE * interp_mRNA;
end