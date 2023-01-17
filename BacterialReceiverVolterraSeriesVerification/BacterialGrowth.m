function dn_dt = BacterialGrowth(t,n_t, mu,Kn, a,b )
%   BacterialGrowth ODE Model for Bacterial Growth
%   Bacterial Growth over time measured as OD600 Light Absorption
% %DBG
% disp("Iteration:")
% r0_str = ["mu0", "Kn0", "a0", "b0"];
% rInit = [mu, Kn, a, b];
% for i = 1:size(r0_str,2)
%     disp(strcat(r0_str(i), ": ", num2str( rInit(i))));
% end

dn_dt(1,1) =  mu * n_t(1).^a * (1-n_t(1)/(Kn*1e8)).^b;

dn_dt(2,1) = 1.6408e-9.*dn_dt(1,1);

end