function H = BacRecvThirdOrderVolterra(order, n,v, q3, kappaApi, kApi,deltaRx, c_tilde, ...
    k1,k2, delta_ccp, a_tilde,kP, k_Dp, delta_Api, betaE, deltaE,input_fct, ...
    tc, constCcp, constApi, constAIP, constRx, constElux)
%UNTITLED First Order Volterra Transfer FUnction/ Frequency Response
%   see Bacterial Receiver
H =@(s1, s2, s3) -(k1/3)*tc*constAIP*( ( ...
                ( -(1/2)*k1.*tc*constAIP* ((((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s1 + tc*(k2 + delta_ccp)))) + ((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s2 + tc*(k2 + delta_ccp))))) ./...      % C_CP = C + AIP / P
                  ((s1 +s2+tc*(k2 + delta_ccp)))) ...
                ) +...                
                (-(1/2)*tc*constAIP*k1.*((((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s1 + tc*(k2 + delta_ccp)))) + ((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s3 + tc*(k2 + delta_ccp))))) ./...     
                 ((s1 +s3 + tc*(k2 + delta_ccp)))) ...
                 ) +...
                (-(1/2)*tc*constAIP*k1.*((((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s2 + tc*(k2 + delta_ccp)))) + ((k1*c_tilde*constAIP*(tc./constCcp)) ./ ((s3 + tc*(k2 + delta_ccp))))) ./...     
                 (s2 +s3 + tc*(k2 + delta_ccp))))...                        
               )./...
                ((s1 + s2 + s3 + tc*(k2 + delta_ccp))) ...
            ).*...
    input_fct(s1) .* input_fct(s2) .*input_fct(s3); 



if(bitget(order,2))

    H =@(s1, s2, s3) H(s1,s2, s3) .* ...  
   -(1/3) * kP * constCcp*tc.* ( ( ( -(1/2)*kP*constCcp*tc.*( (((kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s1 + tc*(k_Dp + delta_Api)))) +  ((kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s2 + tc*(k_Dp + delta_Api)))))./...  % A_PI = C_CP + A
                    (n*v*(s1 + s2 + tc*(k_Dp + delta_Api))))) ...
                    +...              
                (-(1/2)*kP*constCcp*tc.*( (((kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s1 + tc*(k_Dp + delta_Api)))) +  (kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s3 + tc*(k_Dp + delta_Api))))./...   
                    (n * v * (s1 + s3 + tc*(k_Dp + delta_Api)))))  ...
                    +...
                (-(1/2)*kP*constCcp*tc.*( (((kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s2 + tc*(k_Dp + delta_Api)))) +  (kP*a_tilde*tc*(constCcp/constApi)) ./ (n*v*(s3 + tc*(k_Dp + delta_Api))))./...   
                 (n*v*(s2 + s3 + tc*(k_Dp + delta_Api)))))... 
                ) ./...
             (n*v*(s1 + s2 + s3 + tc*(k_Dp + delta_Api)))); 
end


if(bitget(order,3))                                                         % A_Pi -> R_(X)

    H =@(s1, s2, s3) H(s1,s2, s3) .* ...                                         
  (q3*n*kappaApi*(tc/constRx))./((s1+s2+s3 + tc*deltaRx));
end




if(bitget(order,4))                                                         % R_(x) -> E_Lux

     H =@(s1, s2, s3) H(s1,s2, s3) .* ...    
                 ((betaE*constRx*(tc/constElux)) ./ ((s1+s2+s3+tc*deltaE)));
end







%H =@(s1, s2, s3)  ((2)*(n*kappaApi)/((kApi*n*v)^3*(s1+s2+s3 + deltaRx))).*...        % Saturation Term


end




% H =@(s1, s2, s3) ... (q3*(n*kappaApi)/((s1+s2+s3 + deltaRx))).*...        % Saturation Term
%     -(k1/3)*( ( ( -(1/2)*k1.* ((((k1*c_tilde) ./ (s1 + k2 + delta_ccp)) + ((k1*c_tilde) ./ (s2 + k2 + delta_ccp))) ./...      % C_CP = C + AIP / P
%                   (s1 +s2+k2 + delta_ccp)) ...
%                 ) +...                
%                 (-(1/2)*k1.*((((k1*c_tilde) ./ (s1 + k2 + delta_ccp)) + ((k1*c_tilde) ./ (s3 + k2 + delta_ccp))) ./...     
%                  (s1 +s3 + k2 + delta_ccp)) ...
%                  ) +...
%                 (-(1/2)*k1.*((((k1*c_tilde) ./ (s2 + k2 + delta_ccp)) + ((k1*c_tilde) ./ (s3 + k2 + delta_ccp))) ./...     
%                  (s2 +s3 + k2 + delta_ccp)))...                        
%                )./...
%                 (s1 + s2 + s3 + k2 + delta_ccp) ...
%             ).*...     
%   -(1/3) * kP .* ( ( ( -(1/2)*kP.*( (((kP*a_tilde) ./ (n*v*(s1 + k_Dp + delta_Api))) +  ((kP*a_tilde) ./ (n*v*(s2 + k_Dp + delta_Api))))./...  % A_PI = C_CP + A
%                     (n*v*(s1 + s2 + k_Dp + delta_Api)))) +...              
%                 (-(1/2)*kP.*( (((kP*a_tilde) ./ (n*v*(s1 + k_Dp + delta_Api))) +  (kP*a_tilde) ./ (n*v*(s3 + k_Dp + delta_Api)))./...   
%                     (n * v * (s1 + s3 + k_Dp + delta_Api))))  +...
%                 (-(1/2)*kP.*( (((kP*a_tilde) ./ (n*v*(s2 + k_Dp + delta_Api))) +  (kP*a_tilde) ./ (n*v*(s3 + k_Dp + delta_Api)))./...   
%                  (n*v*(s2 + s3 + k_Dp + delta_Api))))... 
%                 ) ./...
%              (n*v*(s1 + s2 + s3 + k_Dp + delta_Api))) .*...
%     1.*...(betaE ./ (s1+s2+s3+deltaE)).*...            % R_X -> Bioluminescent Protein
%     (input_fct(s1) .* input_fct(s2) .* input_fct(s3)); % Input is scaled Heaviside

