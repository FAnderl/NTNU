function H = firstOrderFuncEndocytosis(a, pm, N, ki, kd, InConst)
% Input is Heaviside with C0 = 1.17e13;
% sysLin
H = @(s) (a*pm*N./(ki + s)).*(ki./(kd + s))...
    *InConst./(s);
end