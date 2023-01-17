
function dn_dt = BacGrowthDelayDE(t, n_t,Z, mu,Kn, a,b) 

nlag = Z(:,1);

dn_dt =  (mu * n_t(1).^a * (1-nlag(1)/(Kn*1e8)).^b);

%dn_dt(2,1) = 1.6408e-9.*dn_dt(1,1);

end