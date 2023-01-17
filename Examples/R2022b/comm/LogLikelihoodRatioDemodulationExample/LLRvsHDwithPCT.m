function LLRvsHDwithPCT(EbNo,numLLRRuns)
% LLR vs. Hard Decision Demodulation Example with PCT 
%
%    This file runs parallel simulations at each of the Eb/No values listed
%    in the hard decision demodulation regime and at each of the first
%    numLLRRuns Eb/No values in the LLR demodulation regime. This file is a
%    modified version of the simLLRvsHD.m file, which is split into
%    LLRvsHDwithPCT.m and simLLRvsHDPCT.m functions to achieve
%    parallelization . The data required for the simulation is initialized
%    here, while the System objects are created in the simLLRvsHDPCT.m,
%    which is set up to simulate the communication system for the BER
%    calculations of a single Eb/No value.

% Copyright 2013-2021 The MathWorks, Inc.

%%
% Opening the parallel pool
if isempty(gcp('nocreate'))
    parpool;
end
pool = gcp;
numWorkers = pool.NumWorkers;

lenEbNo = length(EbNo);      % Length of the EbNo range

% Check for validity of numLLRRuns
if (gt(numLLRRuns,lenEbNo))
    error(message('comm:simLLRvsHD:invalidnumLLRRuns1'));
end

% Modulation properties
M = 4;                       % Modulation Order
k = log2(M);
seed = 599*(1:numWorkers);   % Random number generator seed

% Code properties
codeRate = 1/2;              % Code rate 
constlen = 7;                % Constraint length 
codegenpoly  = [171 133];
trellis  = poly2trellis(constlen, codegenpoly);
dSpect   = distspec(trellis,14);

% Simulation parameters
bitsPerIter = 1.2e4;         % Number of bits per iteration
maxNumIters = 1e3;           % Maximum number of iterations
maxNumErrs  = 300;           % Maximum number of bit errors to collect

% Adjust signal-to-noise ratio for coded bits and multi-bit symbols
adjSNR = convertSNR(EbNo,"ebno","BitsPerSymbol",k,"CodingRate",codeRate);

[errs_HD, bits_HD] = deal(zeros(numWorkers,lenEbNo));
[errs_UNQUANTIZED, bits_UNQUANTIZED] = deal(zeros(numWorkers,lenEbNo));
[errs_SD, bits_SD] = deal(zeros(numWorkers,lenEbNo));

% Set up a figure for visualizing BER results
fig = figure;
grid on; hold on;
set(fig.CurrentAxes,'yscale','log','xlim', ...
    [EbNo(1)-1, EbNo(end)+1],'ylim',[1e-6 1]);
xlabel('Eb/No (dB)'); ylabel('BER');
set(fig, 'renderer', 'zbuffer');
title('LLR vs. Hard Decision Demodulation');
grid on; hold on;

% Use 'bercoding' to calculate theoretical upper bounds on BER
theoryBER_HD = bercoding( ...
    adjSNR,'conv','hard',codeRate,dSpect);
theoryBER_LLR = bercoding( ...
    adjSNR(1:numLLRRuns),'conv','soft',codeRate,dSpect);

semilogy(EbNo,theoryBER_HD,'mo-',EbNo(1:numLLRRuns),theoryBER_LLR,'go-');
legend('Hard Decision: Theoretical Upper Bound',...
    'LLR: Theoretical Upper Bound', ...
    'Location', 'SouthWest');

% Running the simulations in parallel
parfor n = 1: numWorkers
    for idx = 1:lenEbNo
        [ber_hd, ber_unquantized, ber_sd] = simLLRvsHDPCT( ...
            EbNo,numLLRRuns,codeRate,trellis,seed(n),bitsPerIter, ...
            maxNumIters/numWorkers,maxNumErrs/numWorkers,idx,k);
        errs_HD(n,idx) = ber_hd(idx,2);
        bits_HD(n,idx) = ber_hd(idx,3);
        errs_UNQUANTIZED(n,idx) = ber_unquantized(idx,2);
        bits_UNQUANTIZED(n,idx) = ber_unquantized(idx,3);
        errs_SD(n,idx) = ber_sd(idx,2);
        bits_SD(n,idx) = ber_sd(idx,3);
    end
end
berHD = sum(errs_HD,1)./sum(bits_HD,1);
berUNQUANTIZED = sum(errs_UNQUANTIZED,1)./sum(bits_UNQUANTIZED,1);
berSD = sum(errs_SD,1)./sum(bits_SD,1);

% Plot results
semilogy(EbNo, berHD, 'r*-', ...
    EbNo(1:numLLRRuns), berUNQUANTIZED(1:numLLRRuns), 'b*-',...
    EbNo(1:numLLRRuns), berSD(1:numLLRRuns), 'k*-');
legend('Hard Decision: Theoretical Upper Bound', ...
    'LLR with unquantized decoding: Theoretical Upper Bound', ...
    'Hard Decision: Simulation' , ...
    'LLR with unquantized decoding: Simulation',...
    'LLR with Soft Decision: Simulation',...
    'Location', 'SouthWest');
