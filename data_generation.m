close all;
clc;
clear;


mod_order=2; %for bpsk
bitsPerSymbol= log2(mod_order); 

numBits = 1e5;  %total number of bits sent for each SNR

snrDb = -10:5:20;  % SNR range in dB
numTx = 4;  % Number of transmit antennas
numRx = 4;  % Number of receive antennas
symbolsPerAntenna=100; % Number of symbols being transmitted per antenna

numOfIterations=numBits/(symbolsPerAntenna*bitsPerSymbol); %we are sending..
% symbolsPerAntenna*bitsPerSymbol bits in each iteration

bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator;

% Initialize arrays to store results
berValues = zeros(length(snrDb), 1);

tx_data=zeros(length(snrDb)*numOfIterations,symbolsPerAntenna*bitsPerSymbol*numTx, 1);
rx_data=zeros(length(snrDb)*numOfIterations,symbolsPerAntenna*bitsPerSymbol*numTx, 1);

snr_values=zeros(length(snrDb)*numOfIterations,1);
cdModulated = comm.ConstellationDiagram('Name', 'Modulated Data');
cdReceived = comm.ConstellationDiagram('Name', 'Received Signal');
cdEqualized = comm.ConstellationDiagram('Name', 'Equalized Signal');
cdDemodulated = comm.ConstellationDiagram('Name', 'Demodulated Signal');



for i = 1:length(snrDb)
    
    snr = 10^(snrDb(i)/10);
    
    totalErrors = 0;
    totalBits = 0;

    iteration_count=0;
    % In each iteration we will be sending symbolsPerAntenna*bitsPerSymbol bits


    while totalBits < numBits  %this is the loop for multiple iterations..
                               % of data generation and BER calculation
        % Generate random binary data
        data = randi([0 1], symbolsPerAntenna*bitsPerSymbol*numTx, 1);  % Column vector
        
        
        % BPSK modulation
        modulatedData = bpskMod(data);
      
        % Reshape modulated data for MIMO transmission
        modulatedDataReshaped = reshape(modulatedData, numTx, []);

        % Generate Rayleigh fading channel
        H = (randn(numRx, numTx) + 1j*randn(numRx, numTx))/sqrt(2);
        
        % Pass through channel
        rxSignal = H * modulatedDataReshaped;
        
        % Add white Gaussian noise
        rxSignalNoisy=awgn(rxSignal,snrDb(i),'measured');
       
        % Perform ZF equalization
        equalized = pinv(H) * rxSignalNoisy;
        
        % Reshape equalized signal for demodulation
        equalizedReshaped = reshape(equalized, [], 1);
        
        % Demodulate
        receivedBits = bpskDemod(equalizedReshaped);
        
        % Calculate bit errors
        numErrors = biterr(data, receivedBits);
        totalErrors = totalErrors + numErrors;
        totalBits = totalBits + (symbolsPerAntenna*bitsPerSymbol);
        
        iteration_count= iteration_count +1;

        % saving transmitted and recieved bits to compare later
        tx_data((i-1)*numOfIterations+iteration_count,:,:)=data;
        rx_data((i-1)*numOfIterations+iteration_count,:,:)=receivedBits;
        
        % saving SNR values
        snr_values((i-1)*numOfIterations+iteration_count,:)=snrDb(i);
    end
    
    % Calculate BER for this SNR
    berValues(i) = totalErrors / totalBits;

    disp(['total errors at ',   num2str(snrDb(i)),  'db =', num2str(totalErrors)]);

    disp(['total bits =', num2str(totalBits)]);
end

% Plot results
semilogy(snrDb, berValues, 'bo-');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for 4x4 MIMO with BPSK and Rayleigh Fading');

%%
% %  Display constellation diagram for modulated data
% cdModulated(modulatedDataReshaped(:));
% 
% 
% % Display constellation diagram for received signal
% cdReceived(rxSignalNoisy(:));
% 
% % Display constellation diagram for equalized signal
% cdEqualized(equalized(:));
% 
% 
% % Display constellation diagram for demodulated signal
% cdDemodulated(receivedBits);
%%
% Save results
save('mimo_data_surya_final.mat', 'tx_data','rx_data','snr_values','-v7.3');
