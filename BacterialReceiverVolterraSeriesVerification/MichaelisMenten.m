function out = MichaelisMenten(in,k, kappa, hill_coeff)
% MichaelisMenten computes the conventional MM saturation term for  
%   in: Enzyme Input CONCENTRATION
%   k : Michaelis-Menten-constan
%   kappa : Maximum Activity

out = kappa * (...
    (in)^(hill_coeff)...
    /...
    (k+in^(hill_coeff))     );
end 