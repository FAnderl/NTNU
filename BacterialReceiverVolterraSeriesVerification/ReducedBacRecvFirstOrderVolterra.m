function H = ReducedBacRecvFirstOrderVolterra(n,v, q1,kappaApi, kApi,deltaRx,input_fct)
%UNTITLED First Order Volterra Transfer Function/ Frequency Response
%   see Bacterial Receiver
H =@(s) ((n*kappaApi*q1)./((s + deltaRx))).*...   % Saturation Terms, i.e., A_Pi -> R_X
    input_fct(s);     % Input is scaled Heaviside
 %H =@(s) ((n*kappaApi)./(kApi*n*v.*(s + deltaRx))).*...   % Saturation Terms, i.e., A_Pi -> R_X 

end

