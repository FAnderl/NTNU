function [ber_HD, ber_UNQUANTIZED, ber_SD] = simLLRvsHD(EbNo, numLLRRuns)
% LLR vs. Hard Decision Demodulation Example
%
%    [BER_HD, BER_LLR, BER_SD] = SIMLLRVSHD(EbNo, numLLRRuns) illustrates
%    how to improve BER performance of a coded communication system by
%    using log-likelihood ratios (LLR) instead of hard decision
%    demodulation in conjunction with a Viterbi decoder. For the LLR case,
%    it also shows how to use the Viterbi decoder in the unquantized mode
%    as well as the more practically viable soft-decision mode. The channel
%    is assumed to be an AWGN channel. EbNo is a vector that contains the
%    signal-to-noise ratios per information bit.This file runs a simulation
%    at each of the EbNo values listed in the hard-decision demodulation
%    regime and at each of the first numLLRRuns EbNo values in the LLR
%    demodulation regime. The BER results are plotted as they are
%    generated.
% 
%   See also comm.ConvolutionalEncoder, comm.QPSKModulator,
%   comm.QPSKDemodulator, comm.ViterbiDecoder, and comm.ErrorRate, AWGN,
%   BERCODING, BERFIT.

% Copyright 2006-2021 The MathWorks, Inc.

% Check for validity of numLLRRuns
if (gt(numLLRRuns,length(EbNo))) 
    error(message('comm:simLLRvsHD:invalidnumLLRRuns1'));
end

% BERFIT requires at least four values for curve-fitting
if (gt(4,numLLRRuns)) 
    error(message('comm:simLLRvsHD:invalidnumLLRRuns2'));
end
% Modulation properties
M = 4;

% Create a rate 1/2, constraint length 7 comm.ConvolutionalEncoder System
% object. This encoder takes one-bit symbols as inputs and generates 2-bit
% symbols as outputs.
codeRate = 1/2;
constlen = 7;
codegenpoly = [171 133];    
trellis = poly2trellis(constlen, codegenpoly);
enc = comm.ConvolutionalEncoder(trellis);
dSpect = distspec(trellis,14);

% Adjust signal-to-noise ratio for coded bits and multi-bit symbols.
adjSNR = convertSNR( ...
    EbNo,"ebno", ...
    "BitsPerSymbol",log2(M), ...
    "CodingRate",codeRate);

% Create comm.ViterbiDecoder System objects to act as the hard-decision,
% unquantized and soft-decision decoders.
vitDecHD = comm.ViterbiDecoder(trellis, ...
    'InputFormat','Hard',...
    'TracebackDepth',32); 
vitDecUNQUANT = comm.ViterbiDecoder(trellis, ...
    'InputFormat','Unquantized', ...
    'TracebackDepth', 32); 
vitDecSD = comm.ViterbiDecoder(trellis, ...
    'InputFormat','Soft', ...
    'SoftInputWordLength',3, ...
    'TracebackDepth',32); 

% Create comm.ErrorRate System objects to compare the decoded bits to the
% original transmitted bits. The Viterbi decoder creates a delay in the
% output decoded bit stream equal to the traceback length. To account for
% this delay set the 'ReceiveDelay' property of the comm.ErrorRate objects
% to 32.
errorCalcHD  = comm.ErrorRate('ReceiveDelay', 32);
errorCalcUNQUANT = comm.ErrorRate('ReceiveDelay', 32);
errorCalcSD  = comm.ErrorRate('ReceiveDelay', 32);

% Before using a comm.ViterbiDecoder object in the 'soft decision' mode,
% the output of the comm.QPSKDemodulator object needs to be quantized. This
% example uses a comm.ViterbiDecoder object with a 'SoftInputWordLength'
% value of 3.
 
% Since the AWGN Channel as well as the RANDI function use the default
% random stream, the following commands are executed so that the results
% will be repeatable, i.e. same results will be obtained for every run of
% the example. The default stream will be restored at the end of the
% example.
s0 = RandStream.getGlobalStream;
s = RandStream.create('mt19937ar', 'seed',12345);
RandStream.setGlobalStream(s);

% Number of bits per iteration
bitsPerIter = 1.2e4;
% Maximum number of iterations
maxNumIters = 100;
% Maximum number of bit errors to collect
maxNumErrs  = 300;

% Initialize variables for storing BER results
ber_HD = zeros(3,length(EbNo));
ber_UNQUANTIZED = zeros(3,numLLRRuns);
ber_SD = zeros(3,numLLRRuns); 

% Set up a figure for visualizing BER results
fig = figure;
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');
ax.YScale = 'log';
xlim(ax, [EbNo(1)-1, EbNo(end)+1]); ylim(ax, [1e-6 1]);
xlabel(ax,'Eb/No (dB)'); ylabel(ax, 'BER');
title(ax,'LLR vs. Hard Decision Demodulation');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
set(fig,'DefaultLegendAutoUpdate','off');

% Use 'bercoding' to calculate theoretical upper bounds on BER
theoryBER_HD = bercoding(adjSNR,'conv','hard',codeRate,dSpect);
theoryBER_LLR = bercoding(adjSNR(1:numLLRRuns), ...
    'conv','soft',codeRate,dSpect);

semilogy(EbNo,theoryBER_HD,'mo-',EbNo(1:numLLRRuns),theoryBER_LLR,'go-');
legend('Hard Decision: Theoretical Upper Bound', ...
    'LLR: Theoretical Upper Bound','Location','SouthWest');

% This is the common loop for hard decision and LLR with unquantized and
% soft decision decoding
for idx=1:numLLRRuns    

    % Reset objects with state at the start of an EbNo value
    reset(errorCalcHD)
    reset(errorCalcUNQUANT)
    reset(errorCalcSD)
    reset(enc)
    reset(vitDecHD)
    reset(vitDecUNQUANT)
    reset(vitDecSD)
    iter=1;
    
    % Get noise variance given SNR and set the NoiseVariance, required for
    % LLR demodulation
    [~,nVar] = awgn(0,adjSNR(idx));

    % Exit loop when either the number of bit errors exceeds 'maxNumErrs'
    % or the maximum number of iterations have completed
    while (ber_UNQUANTIZED(2,idx) < maxNumErrs) && (iter <= maxNumIters)
        data = randi([0 1], bitsPerIter, 1); % Generate message bits        
        encData = enc(data);                 % Convolutionally encode data            
        modData = pskmod(encData,M,pi/4, ...
            InputType="bit");                % Modulate encoded data
        [chOut,nVar] = awgn(modData,adjSNR(idx)); % Pass modulated signal
                                                  % through AWGN channel     
        demodDataHD = pskdemod(chOut,M,pi/4, ...
            OutputType="bit");                    % Hard decision demod
        demodDataLLR = pskdemod(chOut,M,pi/4, ...
            OutputType="llr",NoiseVariance=nVar); % 'LLR' demod
        % Hard-decision decoding: Pass the demodulated data through the
        % Viterbi decoder; and compute and accumulate errors     
        decDataHD = vitDecHD(demodDataHD);
        ber_HD(:,idx) = errorCalcHD(data, decDataHD);
        
        % Unquantized decoding: Pass the demodulated data through the
        % Viterbi decoder; and compute and accumulate errors
        decDataUNQUANT = vitDecUNQUANT(demodDataLLR);
        ber_UNQUANTIZED(:,idx) = errorCalcUNQUANT(data, decDataUNQUANT);
        
        % Soft-decision decoding: The demodulated data must pass through a
        % quantizer before being fed to the Viterbi decoder. However the
        % output from the demodulator must be sign-reversed before being
        % fed to the quantizer. This is because in the soft-decision
        % decoding mode, the comm.ViterbiDecoder object assumes that
        % positive numbers correspond to 1s and negative numbers to 0s.
        % Thus the decoding operation consists of feeding in the sign
        % reversed data from the comm.PSKDemodulator object to the
        % |quantiz| function and feeding in the output to the
        % comm.ViterbiDecoder. The quantization partitions are fine tuned
        % according to the noise SNR. Compute and accumulate errors.
        quantizedValue = quantiz(-demodDataLLR, ...
            (-2.1:.7:2.1)/nVar);
        decDataSD = vitDecSD(double(quantizedValue));
        ber_SD(:,idx) = errorCalcSD(data, decDataSD);
        
        iter = iter+1;
    end
    
    % Plot results
    semilogy(ax, EbNo(1:idx), ber_HD(1,1:idx), 'r*', ...
             EbNo(1:numLLRRuns), ber_UNQUANTIZED(1,1:numLLRRuns), 'b*',...
             EbNo(1:numLLRRuns), ber_SD(1,1:numLLRRuns), 'k*');
    legend('Hard Decision: Theoretical Upper Bound', ...
           'LLR with unquantized decoding: Theoretical Upper Bound', ...
           'Hard Decision: Simulation' , ...
           'LLR with unquantized decoding: Simulation',...
           'LLR with Soft Decision: Simulation',...
           'Location', 'SouthWest');
    drawnow;
end

% If LLR decoding is done for fewer values than hard decision decoding,
% this loop performs hard decision decoding for the remaining EbNo values.
if (numLLRRuns < length(EbNo))
    for idx=numLLRRuns+1:length(EbNo)
             
        % Reset
        reset(errorCalcHD)
        reset(enc)
        reset(vitDecHD)
        iter=1;
        
        while (ber_HD(2,idx) < maxNumErrs) && (iter <= maxNumIters)            
            data = randi([0 1],bitsPerIter,1);
            encData = enc(data);
            modData = pskmod(encData,M,pi/4,InputType="bit");
            chOut = awgn(modData,adjSNR(idx));
            demodDataHD = pskdemod(chOut,M,pi/4,OutputType="bit");
            decDataHD = vitDecHD(demodDataHD);
            ber_HD(:,idx) = errorCalcHD(data,decDataHD);
            iter = iter+1;
        end
        
        semilogy(ax, EbNo(1:idx),ber_HD(1,1:idx),'r*', ...
            EbNo(1:numLLRRuns),ber_UNQUANTIZED(1,1:numLLRRuns),'b*', ...
            EbNo(1:numLLRRuns),ber_SD(1,1:numLLRRuns),'k*');
        drawnow;       
    end
end

% Perform curve fitting and plot the results
fitBER_HD  = berfit(EbNo,ber_HD(1,:));
fitBER_LLR = berfit(EbNo(1:numLLRRuns),ber_UNQUANTIZED(1,:));
fitBER_SD = berfit(EbNo(1:numLLRRuns),ber_SD(1,:));
semilogy(ax,EbNo,fitBER_HD,'r*-',EbNo(1:numLLRRuns),fitBER_LLR,'b*-', ...
        EbNo(1:numLLRRuns),fitBER_SD,'k*-');
hold(ax,'off');

% Restore default stream
RandStream.setGlobalStream(s0);

% [EOF]