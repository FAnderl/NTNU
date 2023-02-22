function dn_dt = BacterialGrowth_wXIP(t,n_t, mu,Kn, a,b, c_XIP, theta_XIP)
%   BacterialGrowth ODE Model for Bacterial Growth
%   Bacterial Growth over time measured as OD600 Light Absorption
% %DBG
% disp("Iteration:")
% r0_str = ["mu0", "Kn0", "a0", "b0"];
% rInit = [mu, Kn, a, b];
% for i = 1:size(r0_str,2)
%     disp(strcat(r0_str(i), ": ", num2str( rInit(i))));
% end


mu_XIP = (mu-(theta_XIP*log(c_XIP)));

dn_dt(1,1) =  mu_XIP * n_t(1).^a * (1-n_t(1)/(Kn*1e8)).^b;

dn_dt(2,1) = 1.6408e-9.*dn_dt(1,1);

end