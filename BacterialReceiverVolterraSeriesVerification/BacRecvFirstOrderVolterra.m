function H = BacRecvFirstOrderVolterra(order,n,v, q1,kappaApi, kApi,deltaRx, c_tilde, k1,k2, delta_ccp, a_tilde,kP, k_Dp, delta_Api, betaE, deltaE,input_fct)
%UNTITLED First Order Volterra Transfer Function/ Frequency Response
%   see Bacterial Receiver



H =@(s) ((k1*c_tilde) ./ (s + k2 + delta_ccp)) .*...      % C_CP = C + AIP / P
    input_fct(s);     % Input is scaled Heaviside




if(bitget(order,2))

    H =@(s) H(s) .* ...                                          % A_PI = C_CP + A 
                ((kP*a_tilde) ./ (n*v*(s + k_Dp + delta_Api)));             
end


if(bitget(order,3))                                                         % A_Pi -> R_(X)

    H =@(s) H(s) .* ...                                         
                ((n*kappaApi*q1)./((s + deltaRx)));
end




if(bitget(order,4))                                                         % R_(x) -> E_Lux

    H =@(s) H(s) .*...   
                  (betaE ./ (s+deltaE));
end


  

%H =@(s) ((n*kappaApi)./((kApi*n*v).*(s + deltaRx))).*...   % Saturation Terms, i.e., A_Pi -> R_X


end



% H =@(s) ((n*kappaApi*q1)./((s + deltaRx))).*...   % Saturation Terms, i.e., A_Pi -> R_X
%     order((k1*c_tilde) ./ (s + k2 + delta_ccp)) .*...      % C_CP = C + AIP / P
%     ((kP*a_tilde) ./ (n*v*(s + k_Dp + delta_Api))) .*...    % A_PI = C_CP + A 
%     1.*...%(betaE ./ (s+deltaE)) .*...     % R_X -> Bioluminescent Protein
%     input_fct(s);     % Input is scaled Heaviside

