function H = ReducedBacRecvThirdOrderVolterra(v,n, q3, kappaApi, kApi,deltaRx,input_fct)
%UNTITLED First Order Volterra Transfer FUnction/ Frequency Response
%   see Bacterial Receiver
H =@(s1, s2, s3)  ((n*kappaApi*q3)/((s1+s2+s3 + deltaRx))).*...        % Saturation Term
    input_fct(s1).*input_fct(s2).*input_fct(s3); % Input is scaled Heaviside


end

