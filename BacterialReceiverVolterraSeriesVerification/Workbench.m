% Fit XIP-affected bacterial growth constants to logarithmic decay model
close all

XIP = [1.900000
    3.9000000
7.8000000
15.6250000
31.2500000
62.5000000
125.0000000
250.0000000
500.0000000
1000.0000000
2000.0000000];

mu = [0.012306118
0.011822689
0.011805228
0.011165575
0.010984436
0.009485443
0.0064108
0.004883943
0.003239925
0.002381929
0.002011366];




tbl = table(log(XIP), mu, 'VariableNames',{'XIP', 'mu'});

figure(1),scatter(XIP, mu);
hold on 
figure(2),scatter(log(XIP), mu);
hold on 

mdl = fit(log(XIP),mu, 'poly1')
figure(2);plot(mdl);

coeffs = coeffvalues(mdl);


figure(1), plot(XIP, coeffs(1)*log(XIP)+coeffs(2));


%% Bacterial Growth Model - Validation


  






function solpts = OptimizationODESolver(r,optim_tspan, init_vec)
    [mu] = r(1);
    [Kn] = r(2);
    [a] = r(3);
    [b] = r(4);
    %[od600_coeff] = r(3);
 
    sol= ode45(@(t,n_t) BacterialGrowth(t,n_t,mu,Kn, a,b), optim_tspan, init_vec);
    dummy_sol=deval(sol,optim_tspan);
    
    % Force output to be stricly real; TODO: check if that is really legit
    solpts = real(dummy_sol(2,:));
    % Check output for complex values
    %disp(dummy_sol(1,:));
end



function solpts = OptimizationODESolver_wOvershoot(r,optim_tspan, init_vec)
    [mu] = r(1);
    [Kn] = r(2);
    [a] = r(3);
    [b] = r(4);
    [tau] = r(5);
 
    sol= dde23(@(t,n_t,Z) BacGrowthDelayDE(t, n_t, Z, mu,Kn, a,b),...
        [tau],@(t) BacHist(t,mu,init_vec),[0, 900]);

    dummy_sol = deval(sol,optim_tspan);
    solpts = real(dummy_sol);
end




function solpts = OptimizationODESolver_wXIP(r,c_XIP,optim_tspan, init_vec)
    [mu] = r(1);
    [Kn] = r(2);
    [a] = r(3);
    [b] = r(4);
    [theta_XIP] = r(5);
    %[od600_coeff] = r(3);
 
    sol= ode45(@(t,n_t) BacterialGrowth_wXIP(t,n_t,mu,Kn, a,b, c_XIP, theta_XIP), optim_tspan, init_vec);
    dummy_sol=deval(sol,optim_tspan);
    
    % Force output to be stricly real; TODO: check if that is really legit
    solpts = real(dummy_sol(2,:));
    % Check output for complex values
    %disp(dummy_sol(1,:));
end


