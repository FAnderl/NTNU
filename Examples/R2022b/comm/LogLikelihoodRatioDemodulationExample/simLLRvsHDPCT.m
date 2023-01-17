function [ber_HD, ber_UNQUANTIZED, ber_SD] = simLLRvsHDPCT(EbNo, numLLRRuns, codeRate, trellis, seed, bitsPerIter, maxNumIters, maxNumErrs, idx, k)
% SIMLLRVSHDPCT Simulate the communication system for the BER calculations
% of a single Eb/No value
%
%    [ber_HD, ber_UNQUANTIZED, ber_SD] = SIMLLRVSHDPCT(EbNo, numLLRRuns,
%    codeRate, trellis, seed, bitsPerIter, maxNumIters, maxNumErrs, idx, k)
%    simulates a coded communication system that uses log-likelihood ratios
%    (LLR) and hard decision demodulation in conjunction with a Viterbi
%    decoder. For the LLR case, it also shows how to use the Viterbi
%    decoder in the unquantized mode as well as the more practically viable
%    soft-decision mode. The channel is assumed to be an AWGN channel. EbNo
%    is a vector that contains the signal-to-noise ratios per information
%    bit. numLLRRuns is the number of LLR runs, codeRate is the code rate,
%    seed is the random number generator seed, bitsPerIter is the number of
%    bits per iteration, maxNumIters is the maximum number of iterations,
%    maxNumErrs is the maximum number of errors, idx is the index of the
%    current Eb/No value and k is the number of bits in the modulated
%    signal.
%
%  See also comm.ConvolutionalEncoder, comm.QPSKModulator,
%  comm.QPSKDemodulator, comm.ViterbiDecoder, comm.ErrorRate, AWGN,
%  QUANTIZ.

% Copyright 2013-2021 The MathWorks, Inc.

%%
% Create a comm.QPSKModulator System object and two comm.QPSKDemodulator
% System objects one each for hard decision and LLR demodulation.
qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemodHD = comm.QPSKDemodulator('BitOutput',true,...
    'DecisionMethod', 'Hard decision');
qpskDemodLLR = comm.QPSKDemodulator('BitOutput',true,...
    'DecisionMethod', 'Log-likelihood ratio');

% Create a rate 1/2, constraint length 7 comm.ConvolutionalEncoder System
% object. This encoder takes one-bit symbols as inputs and generates 2-bit
% symbols as outputs.
enc = comm.ConvolutionalEncoder(trellis);


% Adjust signal-to-noise ratio for coded bits and multi-bit symbols.
adjSNR = convertSNR(EbNo,"ebno","BitsPerSymbol",k,"CodingRate",codeRate);

% Create comm.ViterbiDecoder System objects to act as the hard-decision,
% unquantized and soft-decision decoders.
vitDecHD  = comm.ViterbiDecoder(trellis, 'InputFormat', 'Hard',...
    'TracebackDepth', 32); 
vitDecUNQUANT = comm.ViterbiDecoder(trellis, 'InputFormat', 'Unquantized', ...
    'TracebackDepth', 32); 
vitDecSD  = comm.ViterbiDecoder(trellis, 'InputFormat', 'Soft',...
  'SoftInputWordLength',3, 'TracebackDepth', 32); 

% Create comm.ErrorRate System objects to compare the decoded bits to the
% original transmitted bits. The Viterbi decoder creates a delay in the
% output decoded bit stream equal to the traceback length. To account for
% this delay set the 'ReceiveDelay' property of the comm.ErrorRate objects
% to 32.
errorCalcHD  = comm.ErrorRate('ReceiveDelay', 32);
errorCalcUNQUANT = comm.ErrorRate('ReceiveDelay', 32);
errorCalcSD  = comm.ErrorRate('ReceiveDelay', 32);

% Since the AWGN Channel as well as the RANDI function use the default
% random stream, the following commands are executed so that the results
% will be repeatable, i.e. same results will be obtained for every run of
% the example. The default stream will be restored at the end of the
% example. Use the combined multiple recursive generator since it
% supports substreams.
s0 = RandStream.getGlobalStream;
s = RandStream.create('mrg32k3a', 'seed',seed);
RandStream.setGlobalStream(s);


% Initialize variables for storing BER results
ber_HD = zeros(length(EbNo),3);
ber_UNQUANTIZED = zeros(length(EbNo),3);
ber_SD = zeros(length(EbNo),3);


% This is the common loop for hard decision and LLR with unquantized and
% soft decision decoding
if idx<=numLLRRuns    

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
    while (ber_UNQUANTIZED(idx,2) < maxNumErrs) && (iter <= maxNumIters)
        
        data = randi([0 1], bitsPerIter, 1);      % Generate message bits        
        encData = enc(data);                  % Convolutionally encode data        
        modData = qpskMod(encData);               % Modulate encoded data        
        [chOut,nVar] = awgn(modData,adjSNR(idx)); % Pass modulated signal 
                                                  % through an AWGN channel        
        demodDataHD = pskdemod(chOut,M,pi/4, ...
            OutputType="bit");                    % Hard decision demod
        demodDataLLR = pskdemod(chOut,M,pi/4, ...
            OutputType="llr",NoiseVariance=nVar); % 'LLR' demod
        % Hard-decision decoding: Pass the demodulated data through the
        % Viterbi decoder; and compute and accumulate errors     
        decDataHD = vitDecHD(demodDataHD);
        ber_HD(idx,:) = errorCalcHD(data, decDataHD);
        
        % Unquantized decoding: Pass the demodulated data through the
        % Viterbi decoder; and compute and accumulate errors
        decDataUNQUANT = vitDecUNQUANT(demodDataLLR);
        ber_UNQUANTIZED(idx,:) = errorCalcUNQUANT(data, decDataUNQUANT);
        
        % Soft-decision decoding: The demodulated data must pass through a
        % quantizer before being fed to the Viterbi decoder. However the
        % output from the demodulator must be sign-reversed before being
        % fed to the quantizer. This is because in the soft-decision
        % decoding mode, the comm.ViterbiDecoder object assumes that
        % positive numbers correspond to 1s and negative numbers to 0s.
        % Thus the decoding operation consists of feeding in the sign
        % reversed data from the comm.PSKDemodulator object to the quantiz
        % function and feeding in the output from this quantiz function to
        % the comm.ViterbiDecoder. Compute and accumulate errors. Use a
        % fine tuned quantizer range according to the noise variance.
        quantizedValue = quantiz(-demodDataLLR, (-2.1:.7:2.1)/nVar);
        decDataSD = vitDecSD(double(quantizedValue));
        ber_SD(idx,:) = errorCalcSD(data, decDataSD);
        iter=iter+1;
    end
end

% If LLR decoding is done for fewer values than hard decision decoding,
% this loop performs hard decision decoding for the remaining EbNo values.
if idx>numLLRRuns
        % Reset
        reset(errorCalcHD)
        reset(enc)
        reset(vitDecHD)
        iter=1;
        
        while (ber_HD(idx,2) < maxNumErrs) && (iter <= maxNumIters)
            
            data = randi([0 1], bitsPerIter, 1);
            encData = enc(data);
            modData = qpskMod(encData);
            chOut = awgn(modData,adjSNR(idx));
            demodDataHD = pskdemod(chOut,M,pi/4,OutputType="bit");
            decDataHD = vitDecHD(demodDataHD);
            ber_HD(idx,:) = errorCalcHD(data, decDataHD);
            iter=iter+1;
            
        end
end

% Restore default stream
RandStream.setGlobalStream(s0);

