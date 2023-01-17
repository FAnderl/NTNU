
function dstat_vars_dt = ScaledSMLuxLuminscence(t1, stat_vars1, t,...
    betaE, deltaEp, ks, kp, s_zero, mRNA )

    % EL  
    dstat_vars_dt(1,:) = betaE * interp1(t,mRNA,t1); 
    
    
    % EP 
    dstat_vars_dt(2,:) = betaE * interp1(t,mRNA,t1)- deltaEp * stat_vars1(2); 
    
    
    % P : Product (Fatty Acid)
    dstat_vars_dt(3,:)  = ks * s_zero * stat_vars1(1) - (ks * stat_vars1(1) + ...
        kp * stat_vars1(2)) * stat_vars1(3);
end