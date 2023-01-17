
syms s q q0 kh kb kr ki t


q = q0/((s+kb+kh)-(kr*kb)/(s+kr+ki));
pretty(q)
qPC=partfrac(q,s, "FactorMode","full")
pretty(qPC)
qt = ilaplace(qPC);
pretty(qt)


