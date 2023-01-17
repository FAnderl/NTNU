% Prettify System Transform 

syms n kappaApi kApi delta  c_tilde k1 k2 delta_ccp a_tilde kP  k_Dp  delta_Api betaE deltaE input_fct

H =@(s) ((n*kappaApi)./(kApi.*(s + delta))).*...   % Saturation Terms, i.e., A_Pi -> R_X
    ((k1*c_tilde) ./ (s + k2 + delta_ccp)) .*...      % C_CP = C + AIP / P
    ((kP*a_tilde) ./ (s + k_Dp + delta_Api)) .*...    % A_PI = C_CP + A 
    (betaE ./ (s+deltaE)) %.*...     % R_X -> Bioluminescent Protein
    %(1./s);     % Input is scaled Heaviside

% syms s
pretty(H(s))
ltx = latex(H(s))

imp_response = ilaplace(H(s));

pretty(imp_response)
impr_resp_ltx = latex(imp_response)