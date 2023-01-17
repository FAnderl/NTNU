function bac_gro_at_t = EvalBacGrowth(t,mode,n_sol, varargin)
%UNTITLED Evaluates Bacterial Grwoth Curve at argument t
%   mode 0: Growth Model
%   mode 1: Use Experimental Data for Growth
if mode == 0
    
    bac_gro_at_t = deval(n_sol,t);


elseif mode == 1 % Use Experimental Data

    % n_sol = {t_vec, n_vec}
    t_vec = n_sol{1};
    n_vec = n_sol{2};

    
    bac_gro_at_t = interp1(t_vec,n_vec,t);


end







end