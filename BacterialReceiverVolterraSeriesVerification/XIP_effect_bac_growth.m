% In this script the influence of the XIP concentration 

Kn = 1e5;

mu = 0.0090222; 

n0 = 1e3;


t = linspace(0,900,1000);

% Bacterial Growth Formula
n =  Kn ./ (1+ ((Kn-n0)/n0 ) * exp(-mu*t));



figure, plot(t,n, "r");

