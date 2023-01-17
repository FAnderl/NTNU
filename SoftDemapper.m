function [out,dummy] = SoftDemapper(symbol_in,symMap)
% Performs LLR soft demapping -> DECODING??
%   symbol_in: symbols to demap
%   symMap: Corresponding symbol map ]


% Perform LLR Soft Demapping
partitionPoints = (-1.5:0.5:1.5);


LLRData = pskdemod(symbol_in,4,OutputType="llr");



out = LLRData;
dummy = 0;
end