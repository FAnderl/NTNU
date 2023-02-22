close all

% Investigation of Influence of XIP on Bacterial Growth
dataXIP_In = [ 
    0;
    0;
    0;
    2000;
    2000;
    2000;
    1000;
    1000;
    1000;
    500;
    500;
    500;
    250;
    250;
    250;
    125;
    125;
    125;
    62.5;
    62.5;
    62.5;
    31.25;
    31.25;
    31.25;
    15.625;
    15.625;
    15.625;
    7.8;
    7.8;
    7.8;
    3.9;
    3.9;
    3.9;
    1.9;
    1.9;
    1.9        ];


tspan = [0 1000];

mu0 = 0.0090222; 

Kn = 7.921;               % (*1e8)


theta_XIP = 1;

od600_coeff = 1.64e-9;


n_init = 0.016 * (1/od600_coeff);   % Initial Value before experiments start (from Gabriela's Lab Journal) 0.089 vs 0.016

a = 1;               % Scaling Parameters

b = 1;                  % Scaling Parameters

g = 500;                % Overshoot Parameter; TODO Remove

time = 0 : 1: 1000; 

dn_dt = zeros(2,1);

% Solve for Bacterial Growth
[t,n_t] = ode45(@(t,n_t) BacterialGrowth(t,n_t,mu0,Kn, a,b), tspan, [n_init,0.016]);




figure, plot(t/60, n_t(:,1), '-o','Color','m', 'LineWidth',3.0);
xlabel("Time",'FontSize',15);
ylabel("Population size [Number of Cells]",'FontSize',15);
%set(gca, 'YScale', 'log')
title("Bacterial Growth: $\mu_0=$" + num2str(mu0)+" OD600coeff="+num2str(od600_coeff), 'Interpreter','latex', 'FontSize',15);
figure, plot(t/60, n_t(:,2), '-o', 'Color','r','LineWidth',3.0);
xlabel("Time",'FontSize',15);
ylabel("Population size [OD600]",'FontSize',15);
title("Bacterial Growth: $\mu_0=$" + num2str(mu0)+" OD600coeff="+num2str(od600_coeff), 'Interpreter','latex', 'FontSize',15);



%% (SKIP) Heuristic for XIP Influence on Bacterial Growth
% Hypothesis: XIP leads to slower cell metabolism

d = 0.30;  % Coefficient capturing the effect of XIP on cell metabolism 

XIP = 0.250;   % muM

% Model for how XIP affects bacterial growth
mu = mu0 - d*XIP^2;

% Solve for Bacterial Growth
[t1,n_t1] = ode45(@(t1,n_t1) BacterialGrowth(t1,n_t1,mu,Kn, a,b), tspan, [n_init,0]);



figure, plot(t1/60, n_t1(:,1), '-o', 'LineWidth',3.0);
xlabel("Time",'FontSize',15);
ylabel("Population size [Number of Cells]",'FontSize',15);
title("Bacterial Growth under XIP : $XIP=$" + num2str(XIP)+"nM", 'Interpreter','latex', 'FontSize',15);





%% Data Preparation 

% Parameter Fitting

cfu_opts=detectImportOptions('../ExperimentalData/SM091 Absorbance OD600 first experiment.csv');
cfu_data = readtable('../ExperimentalData/Cell number CFU per mililiter.csv'...
    ,cfu_opts, 'ReadVariableNames', true);
% Load Time Data from dataset
cfu_time_data = cfu_data.Var1;

% Date String
cfu_time_conv_str = datestr(cfu_time_data, 'HH:MM');

% Number of days relative to January 1st XXXX (TODO: look up)
cfu_time_conv = datenum(cfu_time_conv_str,'HH:mm');

cfu_time_conv2 = datetime(cfu_time_conv_str,'InputFormat','HH:mm');

cfu_data_duration = diff([cfu_time_conv2(1) ,cfu_time_conv2(end)]);

cfu_data_duration_min = minutes(cfu_data_duration);

% Load Luminescence Data 
[dataset_size_row, dataset_size_cols] = size(cfu_data);


% The data starts from column 3
cfu_data_set = zeros(dataset_size_row, dataset_size_cols-2);


for i=1:size(cfu_data_set,2)
    dummy = cfu_data(:, i + 2);
    cfu_data_set(:,i) = table2array(dummy);
end

cfu_time_vec = linspace(0, cfu_data_duration_min, size(cfu_time_conv_str,1));


opts=detectImportOptions('../ExperimentalData/SM091 Absorbance OD600 first experiment.csv');
od600_data_01 = readtable('../ExperimentalData/SM091 Absorbance OD600 first experiment.csv'...
    ,opts, 'ReadVariableNames', true);



% Load Time Data from dataset
time_data = od600_data_01.Var1;

% Date String
time_conv_str = datestr(time_data, 'HH:MM');

% Number of days relative to January 1st XXXX (TODO: look up)
time_conv = datenum(time_conv_str,'HH:mm');

time_conv2 = datetime(time_conv_str,'InputFormat','HH:mm');

data_duration = diff([time_conv2(1) ,time_conv2(end)]);

data_duration_min = minutes(data_duration);


time_vec = linspace(0, data_duration_min, size(time_conv_str,1));

% Load Luminescence Data 
[dataset_size_row, dataset_size_cols] = size(od600_data_01);


% The data starts from column 3
data_set = zeros(dataset_size_row, dataset_size_cols-2);


for i=1:size(data_set,2)
    dummy = od600_data_01(:, i + 2);
    data_set(:,i) = table2array(dummy);
end



% Save Bacterial Growth Curves for later use
save("bac_growth_data.mat", "time_vec", "data_set");




%% OD600 - CFU Standard Curve Construction
sol_arr = [];
for i = 4:(length(dataXIP_In))

% Choose which bacterial growth curve to pick 
XIP_idx = i;

% % Plot OD & CFU in same plot
% figure, yyaxis left
% plot(time_vec/60, data_set(:,XIP_idx), '-.', 'LineWidth',3.0);
% xlabel("Time [h]",'FontSize',15);
% ylabel("OD600 ",'FontSize',15);
% title("Bacterial Growth (Control)", 'Interpreter','latex', 'FontSize',15);
% yyaxis right
% plot(cfu_time_vec/60, cfu_data_set(:,XIP_idx), '^', 'LineWidth',3.0, 'Color','r');
% ylabel("CFU ",'FontSize',15)

% Construct Standard Curve
standard_curve_cfu_data = cfu_data_set(:,XIP_idx);

standard_curve_od600_data = data_set(:,XIP_idx);
idx_dummy = ismember(time_vec, cfu_time_vec);
match_idxs = find(idx_dummy);

% Linear Fit to Data
p = polyfit(standard_curve_cfu_data,standard_curve_od600_data(match_idxs),1);
fit = polyval(p,linspace(0,standard_curve_cfu_data(end)));

% figure
% plot(standard_curve_cfu_data, standard_curve_od600_data(match_idxs), 'o', 'LineWidth',3.0);
% xlabel("CFU/ml",'FontSize',15);
% ylabel("OD600 ",'FontSize',15);
% hold on
% plot(linspace(0,standard_curve_cfu_data(end)), fit, 'LineWidth',3.0);
% hold off
% txt = {"Slope Linear Fit: " + num2str(p(1))};
% text(6e7,max(fit)/2, txt,'BackgroundColor',[0.9290 0.6940 0.1250])
% title("ODE600 CFU Standard Curve")





% Overwrite od600 constant
od600_coeff = p(1);


% Overwrite n_init
n_init = 0.016 * (1/od600_coeff);



%%% Bacterial Growth Paramter Estimation (/wo XIP Influence)

% Create Optimization Variable
r = optimvar('r',4, Type='continuous', LowerBound=[0.002,Kn,a,0.5],UpperBound=[0.05,Kn,a,5]);


init_vec = [n_init 0.016];


myfcn = fcn2optimexpr(@OptimizationODESolver,r,time_vec,init_vec);


obj = sum(sum((myfcn - data_set(:,XIP_idx).').^2));


prob = optimproblem("Objective",obj);


r0.r= [mu0 Kn a b];

options = optimoptions("lsqnonlin",...
    'MaxFunctionEvaluations',1500,'Display','iter');
[rsol,sumsq] = solve(prob,r0, 'Options',options);

sol_arr = [sol_arr ; (rsol.r).'];

figure, plot(time_vec/60, data_set(:,XIP_idx), '-o', 'LineWidth',3.0, 'Color','b', 'DisplayName','Experimental Data');
ylabel("Population size [Number of Cells]",'FontSize',15);
hold on 
fittedSol = OptimizationODESolver(rsol.r,time_vec, init_vec);
plot(time_vec/60, fittedSol, '--', 'LineWidth',3.0, 'Color','r', 'DisplayName','Fitted Model');
xlabel("Time [h]",'FontSize',15);
ylabel("Population size [OD600]",'FontSize',15);
title("Experimental Data vs. Fitted Model : $XIP=$" + num2str(dataXIP_In(XIP_idx))+"nM", 'Interpreter','latex', 'FontSize',15);
% legend;
% txt = {'Fitted Parameters:', '$\mu  K_n$ ', num2str(rsol.r(1)) + " " + ...
%     num2str(rsol.r(2))};
% text(2,1.2,txt, 'Interpreter','latex')
hold off


r0_str = ["mu0", "Kn0", "a0", "b0"];
rInit = [mu0, Kn, a, b];
for i = 1:size(r0_str,2)
    disp(strcat(r0_str(i), ": ", num2str( rInit(i))));
end
disp('-------------------------');
r_str = ["mu", "Kn", "a", "b"];
for i = 1:size(r_str,2)
    disp(strcat(r_str(i), ": ", num2str( rsol.r(i))));
end

end
%% Bacterial Growth Paramter Estimation (/w XIP Influence)
sol_arr = [];
for i = 1

% Choose which bacterial growth curve to pick 
XIP_idx = i;


mu0 = 0.013773; %%%0.014879802841075; 

Kn = 7.921;               % (*1e8)

theta_XIP = 0.001543;  %0.001724958390316;

od600_coeff = 1.64e-9;

n_init = 0.016 * (1/od600_coeff);   % Initial Value before experiments start (from Gabriela's Lab Journal) 0.089 vs 0.016

a = 1;                  % Scaling Parameters

b = 1;                  % Scaling Parameters



% Create Optimization Variable
r = optimvar('r',5, Type='continuous', LowerBound=[mu0-mu0/2,Kn-Kn/2,1,0.5, theta_XIP-theta_XIP/2] ...
    ,UpperBound=[mu0+mu0/2,Kn+Kn/2,1,4,theta_XIP+theta_XIP/2]);


init_vec = [n_init 0.016];


c_XIP =  dataXIP_In(4:18); 


% Test
%mu_XIP = (mu0-(theta_XIP*log(c_XIP)));


myfcn = fcn2optimexpr(@OptimizationODESolver_wXIP,r, c_XIP ,time_vec,init_vec);


obj = sum(sum((myfcn - data_set(:,4:18).').^2));


prob = optimproblem("Objective",obj);


r0.r= [mu0 Kn a b theta_XIP];

options = optimoptions("lsqnonlin",...
    'MaxFunctionEvaluations',10000,'Display','iter');
[rsol,sumsq] = solve(prob,r0, 'Options',options);



sol_arr = [sol_arr ; (rsol.r).'];

for j = 1: length(c_XIP)
figure, plot(time_vec/60, data_set(:,3+j), '-o', 'LineWidth',3.0, 'Color','b', 'DisplayName','Experimental Data');
ylabel("Population size [Number of Cells]",'FontSize',15);
hold on 
fittedSol = OptimizationODESolver_wXIP(sol_arr,c_XIP(j), time_vec, init_vec);
plot(time_vec/60, fittedSol, '--', 'LineWidth',3.0, 'Color','r', 'DisplayName','Fitted Model');
xlabel("Time [h]",'FontSize',15);
ylabel("Population size [OD600]",'FontSize',15);
title("Experimental Data vs. Fitted Model : $XIP=$" + num2str(dataXIP_In(XIP_idx))+"nM", 'Interpreter','latex', 'FontSize',15);
% legend;
% txt = {'Fitted Parameters:', '$\mu  K_n$ ', num2str(rsol.r(1)) + " " + ...
%     num2str(rsol.r(2))};
% text(2,1.2,txt, 'Interpreter','latex')
hold off
end

r0_str = ["mu0", "Kn0", "a0", "b0", "theta0"];
rInit = [mu0, Kn, a, b, theta_XIP];
for i = 1:size(r0_str,2)
    disp(strcat(r0_str(i), ": ", num2str( rInit(i))));
end
disp('-------------------------');
r_str = ["mu", "Kn", "a", "b", "theta"];
for i = 1:size(r_str,2)
    disp(strcat(r_str(i), ": ", num2str( rsol.r(i))));
end

end
%% Parameter Inference Including Overshoot (Delay DE Approach)

tau = 80.555;  % Initial Guess for the DELAY-DE delay parameter

n_init_DelayDE = n_init;

dn_dt = zeros(1,1);


% Solve Delay Differential Equation for Heuristic Fit 
sol = dde23(@(t,n_t,Z) BacGrowthDelayDE(t, n_t, Z, mu0,Kn, a,b),...
        [tau],@(t) BacHist(t,mu0,n_init_DelayDE),[0, 900]);

%Test
solpts = deval(sol,time_vec);

figure;
plot(sol.x,sol.y,'-o', 'LineWidth',3.0,'Color','m')
title('Bacterial Growth - Modelled as Delay Differential Equation','FontSize',15);
xlabel('time t','FontSize',15);
ylabel('solution y','FontSize',15);


% Non-Linear Least-Squares Fit to Data
% Create Optimization Variable
% r = [
% mu;
% kn; * 1e8
% a 
% b
% tau
% ]


r = optimvar('r',5, Type='continuous', LowerBound=[0.006,7,1,1,50],UpperBound=[0.01,8,1,1,120]);


init_vec = [n_init];


myfcn = fcn2optimexpr(@OptimizationODESolver_wOvershoot,r,time_vec,init_vec);


obj = sum(sum((myfcn*od600_coeff - data_set(:,XIP_idx).').^2));


prob = optimproblem("Objective",obj);


r0.r= [mu0 Kn a b tau];

options = optimoptions("lsqnonlin",...
    'MaxFunctionEvaluations',10000,'Display','iter');
[rsol,sumsq] = solve(prob,r0, 'Options',options);


figure, plot(time_vec/60, data_set(:,XIP_idx), '-o', 'LineWidth',3.0, 'Color','b', 'DisplayName','Experimental Data');
ylabel("Population size [Number of Cells]",'FontSize',15);
hold on 
fittedSol = od600_coeff.* OptimizationODESolver_wOvershoot(rsol.r,time_vec, init_vec);
plot(time_vec/60, fittedSol, '--', 'LineWidth',3.0, 'Color','r', 'DisplayName','Fitted Model');
xlabel("Time [h]",'FontSize',15);
ylabel("Population size [OD600]",'FontSize',15);
title("Experimental Data vs. Fitted Model : $XIP=$" + num2str(dataXIP_In(XIP_idx))+"nM", 'Interpreter','latex', 'FontSize',15);
% legend;
% txt = {'Fitted Parameters:', '$\mu  K_n$ ', num2str(rsol.r(1)) + " " + ...
%     num2str(rsol.r(2))};
% text(2,1.2,txt, 'Interpreter','latex')
hold off


r0_str = ["mu0", "Kn0", "a0", "b0", "tau0"];
rInit = [mu0, Kn, a, b, tau];
for i = 1:size(r0_str,2)
    disp(strcat(r0_str(i), ": ", num2str( rInit(i))));
end
disp('-------------------------');
r_str = ["mu", "Kn", "a", "b", "tau"];
for i = 1:size(r_str,2)
    disp(strcat(r_str(i), ": ", num2str( rsol.r(i))));
end



%% Local Functions


function dn_dt = BacGrowthDelayDE(t, n_t,Z, mu,Kn, a,b) 

nlag = Z(:,1);

dn_dt =  (mu * n_t(1).^a * (1-nlag(1)/(Kn*1e8)).^b);

%dn_dt(2,1) = 1.6408e-9.*dn_dt(1,1);

end

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
 
    solpts = [];
    for i = 1 : length(c_XIP)
        
        c_XIP_i= c_XIP(i);
        sol= ode45(@(t,n_t) BacterialGrowth_wXIP(t,n_t,mu,Kn, a,b, c_XIP_i, theta_XIP), optim_tspan, init_vec);
        dummy_sol=deval(sol,optim_tspan);
    
    % Force output to be stricly real; TODO: check if that is really legit
    solpts = [ solpts; real(dummy_sol(2,:))];

    end


end

