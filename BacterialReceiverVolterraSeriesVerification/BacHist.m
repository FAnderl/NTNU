function s = BacHist(t, mu0,n_init_DelayDE)
% History for Bacterial Population

% Calculate approximate end using growth ODE soltuion WITHOUT quadratic
% damping
tspan_end = log(n_init_DelayDE/1)*(1/mu0);


time = 0:1:tspan_end;
% Solving equation in interval [0 tspan_end]
s = (1 .* exp(mu0.*t)).';

% TODO Constant History of Bacterial Population Size
s = n_init_DelayDE;

end