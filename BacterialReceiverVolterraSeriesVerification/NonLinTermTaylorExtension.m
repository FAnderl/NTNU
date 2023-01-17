% Derivation of Polynpmial Series Representation of the Saturation Term 
% syms x
% g = exp(x*sin(x));
% t = taylor(g, 'ExpansionPoint', 2, 'Order', 12);
% xd = 1:0.05:3;
% yd = subs(g,x,xd);
% fplot(t, [1, 3])
% hold on
% plot(xd, yd, 'r-.')
% title('Taylor approximation vs. actual function')
% legend('Taylor','Function')



U = @(x) x -1;
L = @(x) U(x)+1; 

L(10)




close all

t = 0:0.001:0.99;
beta = 0.2:0.2:1;
k = t./(t+3); 

figure(1),plot(t,k, 'r.-')
title("Saturation Function")

syms x k n v a_tilde diffTerm


term = x/(x+(k*n*v)); 
diffTerm = term;

kEval = 1.0;
nEval = 1000.0;
vEval = 1.0;
a_tildeEval = 50;

term =subs(term,k,kEval);
exp_point = (kEval*nEval*vEval);
repr = taylor(term, 'Order',4, 'ExpansionPoint', exp_point); 

figure(2),fplot(subs(repr,{v,n}, {vEval,nEval}) ,[kEval*nEval*vEval - (kEval*nEval*vEval)/2,...
    kEval*nEval*vEval + (kEval*nEval*vEval)/2], 'c-o', 'LineWidth',2.0)

%+ "Saturation Constant:" + num2str(beta(i))
hold on
xd = kEval*nEval*vEval - (kEval*nEval*vEval)/2 : 1: kEval*nEval*vEval + (kEval*nEval*vEval)/2;
yd = double(subs(term,{x,v,n},{xd,vEval,nEval})); 
plot(xd, yd, 'r-.');


title("Taylor Series Representation; Expansion Point: " + num2str(exp_point))
legend('Taylor Series','Original Function')



order = 4;
Dfs = sym(zeros(order,1));



symevalDfs = sym(zeros(order,1));
evalDfs = zeros(order,1);

syms symExpPt

symExpPt = (k*n*v)/2; % NB! Needs to be adjusted to line 35 when changes are made 

for i=1 : size(Dfs,1)
    Dfs(i) = diff(diffTerm,x,i);
    symevalDfs(i,1) = subs(Dfs(i),x,symExpPt);
    evalDfs(i) = subs(symevalDfs(i), {k,n,v}, {kEval,nEval,vEval});
end

