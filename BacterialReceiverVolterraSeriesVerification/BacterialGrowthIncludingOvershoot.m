function dn_dt = BacterialGrowthIncludingOvershoot(t,n_t, mu,Kn, a,b)
%   BacterialGrowth ODE Model for Bacterial Growth
%   Bacterial Growth over time measured as OD600 Light Absorption


dn_dt(1,1) =  mu * n_t(1).^a * (1-n_t(1)/(Kn*1e8)).^b;

dn_dt(2,1) = 1.6408e-9.*dn_dt(1,1);


end