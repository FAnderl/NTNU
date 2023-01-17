function H = ReducedBacRecvSecondOrderVolterra(n,v, q2, kappaApi, kApi,deltaRx, input_fct)
%UNTITLED Second Order Volterra Transfer FUnction/ Frequency Response
%   see Bacterial Receiver
H =@(s1, s2)  (n*kappaApi*q2)./(((s1+s2 + deltaRx))).*...
    input_fct(s1) .* input_fct(s2);  % Input is scaled Heaviside
    
end



