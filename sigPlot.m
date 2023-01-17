t = linspace(0,100,10000);
sigma=1;
v=0;
sig = (t/sigma.^2) .* exp(-(t.^2+v.^2)/2*sigma.^2);
close all
figure(1); plot(t,sig );

rx_sig = awgn(sig,20);
figure(2); plot(t, rx_sig);


