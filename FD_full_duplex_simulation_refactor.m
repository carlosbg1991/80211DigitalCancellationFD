% ============================ IMPORTANT NOTE ============================
% This code evaluates two Digital Cancellation methods for a Full Duplex
% system, in the Frequency and Time domain respectivelly. To that end, two
% WiFi frames will be generated, passed through a TGn channel Model A,
% and AWGN. Impairments will be incorporated at the receiver such as
% Frequency Offset or IQ Imbalanced, in order to emulate a real radio
% behavior.
% 
% Since no WLAN System Toolbox was available at the moment of implementing 
% this code, we took the work in [1] as a reference, and adapted their C++ 
% code into Matlab code accordingly.
% 
% Pleace specify there the following parameters:
% (0) SNR            - SNR for the desired signal
% (1) Prate          - Power Ratio between the desired signal and the 
%                      interference in linear scale.
% (2) offset         - Time Offset between the Desired and the 
%                      Interfering signal.
% (3) can_meth       - Cancellation method: 'freq' or 'time'
% (4) FO             - Impairment at RX: Frequency offset
% (5) IQImb_gain     - Impairment at RX: IQ Imbal - Gain missmatch
% (6) IQImb_phase    - Impairment at RX: IQ Imbal - Phase Missmatch

clear; close all;

%% --------------------------- VARIABLES ------------------------------- %%
SNR             = 15;               % Signal-to-Noise Ratio
Prate           = 0.5;              % Relationship betweeen received sgn
                                    % and Interference in linear scale
offset          = 320;              % Offset Desired signal vs Interference
can_meth        = 'freq';           % Cancellation method (time or freq)
FO              = 20e3;             % Frequency Offset in Hz
IQImb_gain      = 0.002;            % Receiver Gain Missmatch between I-Q Branches
IQImb_phase     = 0.2;              % Receiver Phase Missmatch between I-Q Branches
% ----------------------------------------------------------------------- %

Fs              = 10e6;          % Sampling Freq used in USRP
nsamplesXblock  = 80;            % Samples / Block 802.11a
thresholdUSRP   = 4e-3;          % Threshold TX USRP
thresholdUSRP2  = 4e-3;          % Threshold RX USRP
threshold2round = 4e-3;          % Threshold Desired Packet
noise_ref       = -50;           % Reference Noise Power in dBm measured in USRP

load('shortp.mat','short_preamble');                % Known Short Preamble in 802.11a
load('longp.mat','long_preamble');                  % Known Long Preamble in 802.11a
tx_frame       = load('frame.mat','frame');         % WiFi presotred Frame using 802.11a
tx_frame       = tx_frame.frame;
input_samples  = [zeros(10000,1);tx_frame;zeros(10000,1)];

%% --------------------------- FIXED VARIABLES ------------------------- %%
% This parameters are fixed since the WiFi frame is a prestored frame.
% Please, leave them as is.
message         = 'Hello World';    % Message expected to be received
rate            = 1/2;           % 1/2,       3/4,       2/3
M               = 2;             % 2  (BPSK),  4 (QPSK),  16 (16QAM)
n_bpsc          = 1;             % 1  (BPSK),  2 (QPSK),  4 (16QAM)
n_cbps          = 48;            % 48 (BPSK), 96 (QPSK), 192 (16QAM)
n_dbps          = n_cbps*rate;   % Bits per OFDM Symbol
psdu_size       = 39;
n_sym           = ceil((16 + 8*psdu_size + 6)/n_dbps);
n_data          = n_sym * n_dbps;
n_pad           = n_data - (16 + 8*psdu_size + 6);
n_encoded_bits  = n_sym*n_cbps;
nblocks         = n_sym + 1;     %OFDM Symbols + Signal Field

%========================================================================%
%***************************FRAME CREATION*******************************%
%========================================================================%

interference   = load('frame.mat','frame');
interference   = interference.frame;
desired        = interference;

% RAYLEIGH FADING CHANNEL
Ts = 1/Fs;                                             % Sample time
fd = 0.1;                                              % Maximum Doppler shift
% MODEL A
tauA = [0 10 20 30 40 50 60 70 80]*1e-9;               % Vector of path delays
pdbA = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9];    % Vector of average path gains
h2A1 = rayleighchan(Ts,fd,tauA,pdbA);                  % Rayleigh Channel
h2A2 = rayleighchan(Ts,fd,tauA,pdbA);                  % Rayleigh Channel
h2A1.StoreHistory = 1;
h2A2.StoreHistory = 1;
int_after_rayleighA = filter(h2A1, interference);    % Signal after Rician Channel
des_after_rayleighA = filter(h2A2, desired);         % Signal after Rician Channel
Pi = fFD_Power(int_after_rayleighA);
Pd = fFD_Power(des_after_rayleighA);

Pawgn_desired = 10^((noise_ref-30)/10);
awgn_noise = randn(length(des_after_rayleighA) + offset,1) + 1i*randn(length(des_after_rayleighA) + offset,1);
Pawgn_actual = fFD_Power(awgn_noise);
awgn_noise = sqrt(Pawgn_desired/Pawgn_actual).*awgn_noise;
Eawgn = Pawgn_desired;
awgn = awgn_noise;

ErxD = fFD_Power(des_after_rayleighA);
ErxI = fFD_Power(int_after_rayleighA);

SNR_lin = 10^(SNR/10);
des_after_rayleighA1 = (sqrt((SNR_lin*Eawgn)/ErxD)).*des_after_rayleighA;
int_after_rayleighA1 = (sqrt((SNR_lin*Eawgn)/(ErxI*Prate))).*int_after_rayleighA;

des_after_rayleighA1 = [zeros(offset,1);des_after_rayleighA1];
int_after_rayleighA1 = [int_after_rayleighA1;zeros(offset,1)];
out = des_after_rayleighA1 + int_after_rayleighA1;

out_awgn = out + awgn;            % With AWGN noise
aw_append = [];
for i=1:5
    aw_append = [awgn;aw_append];       %#ok<AGROW>
end
rx_frame = [aw_append;out_awgn;aw_append];

% FREQUENCY AND PHASE OFFSET IN THE RX
output_samples1  = fFD_frequency_offset(Fs,FO,rx_frame);

% IQ IMBALANCE
real = (output_samples1+conj(output_samples1))./2;
imag = (output_samples1-conj(output_samples1))./2i;
imag = (1+IQImb_gain)*(sin(IQImb_phase).*real + cos(IQImb_phase).*imag);
output_samples = real + imag*1i;

%========================================================================%
%***************************CANCELLATION ROUND***************************%
%========================================================================%

[tx_pak]   = fFD_preamble_detect_final(Fs, input_samples, short_preamble, thresholdUSRP);
[rx_pak,rx_folp] = fFD_preamble_detect_final(Fs, output_samples, short_preamble, thresholdUSRP2);
[tx_packet,tx_lptot,rx_fosp,tx_num_pkt] = fFD_symbolalign2(Fs, tx_pak, long_preamble);
[rx_packet,rx_lptot,~      ,rx_num_pkt] = fFD_symbolalign2(Fs, rx_pak, long_preamble);
[tx_bits,~,~,tx_message,~] = fFD_decode_packet( tx_packet , nsamplesXblock, tx_num_pkt, nblocks, n_bpsc, n_cbps, n_dbps, psdu_size);

%     if(rx_num_pkt ~= tx_num_pkt)
%         fprintf('Num packet Missmatch. Double check the threshold - Missmatch \n');
%     end

cancelled = zeros((nblocks+10)*nsamplesXblock,min(rx_num_pkt,tx_num_pkt));
cancelled_t = zeros(8000,min(rx_num_pkt,tx_num_pkt));
for packet = 1:min(rx_num_pkt,tx_num_pkt)
    tx_preamble1 = tx_lptot(end-63:end,packet);
    tx_preamble2 = tx_lptot(end-64-63:end-64,packet);
    tx_preamble  = (tx_preamble1 + tx_preamble2)./2;
    fft_tx_preamble = fft(tx_preamble,64);
    rx_preamble1 = rx_lptot(end-63:end,packet);
    rx_preamble2 = rx_lptot(end-64-63:end-64,packet);
    rx_preamble  = (rx_preamble1 + rx_preamble2)./2;
    fft_rx_preamble = fft(rx_preamble,64);

    %NOISE POWER ESTIMATION
    Pn = 10*log10(fFD_Power(rx_preamble-rx_preamble2)) + 30;
    Pnr = 10*log10(fFD_Power(rx_packet(4000:5000,packet))) + 30;
    Prx = 10*log10(fFD_Power(rx_preamble2) - fFD_Power(rx_packet(4000:5000,packet))) + 30;
    Prxr = 10*log10(fFD_Power(rx_preamble2) - fFD_Power(rx_packet(4000:5000,packet))) + 30;
    SNRest = Prx - Pn;
    SNRestr = Prxr - Pnr;

    %CHANNEL ESTIMATION
    H_EST = fFD_channel_est(fft_tx_preamble,fft_rx_preamble);
    H_ESTnew = [H_EST(39:end);H_EST(2:27)];

    %NULL SUBCARRIERS INTERPOLLATION - LINEAR
    for k = 28:38
        alpha = (k-27)/(39-27);
        H_EST(k) = H_EST(27)+( (H_EST(39)-H_EST(27) )*alpha );
    end

    %RECEIVED BLOCKS
    tx_blocks = fFD_blocks_extraction2(tx_packet(:,packet),nsamplesXblock,nblocks+10);
    rx_blocks = fFD_blocks_extraction2(rx_packet(:,packet),nsamplesXblock,nblocks+10);

    if(size(strfind(can_meth,'freq'))>0)
       % CANCELLATION FREQ DOMAIN
        cancellation = fFD_cancellation2(tx_blocks,rx_blocks,H_EST,nblocks+10);
        cancelled(:,packet) = reshape(cancellation, 1, nsamplesXblock*(nblocks+10));
	elseif(size(strfind(can_meth,'time'))>0)
        % CANCELLATION TIME DOMAIN IIR & PA (ASYNCHRONOUS)
        Me = 10;
        Nn = 128;
        xprn = tx_lptot(end-64-63:end,packet).';
        frame = [tx_lptot(end-64-63:end,packet);tx_packet(1:end-128,packet)].';
        K = size(frame,2);
        y4_awgn = [rx_lptot(end-64-63:end,packet);rx_packet(1:end-128,packet)];

        xp4_10 = [zeros(1,Me-1).';frame.'].';
        xp4_01 = [zeros(1,Me-1).';conj(frame).'].';
        xp4_03 = [zeros(1,Me-1).';xp4_01.^3.'].';
        xp4_12 = [zeros(1,Me-1).';((abs(xp4_10).^2).*xp4_01).'].';
        xp4_21 = [zeros(1,Me-1).';((abs(xp4_01).^2).*xp4_10).'].';
        xp4_30 = [zeros(1,Me-1).';xp4_10.^3.'].';

        xpre4_10 = [zeros(1,Me-1).';frame(1:Nn).'].';
        xpre4_01 = [zeros(1,Me-1).';conj(frame(1:Nn)).'].';
        xpre4_03 = [zeros(1,Me-1).';xpre4_01.^3.'].';
        xpre4_12 = [zeros(1,Me-1).';((abs(xpre4_10).^2).*xpre4_01).'].';
        xpre4_21 = [zeros(1,Me-1).';((abs(xpre4_01).^2).*xpre4_10).'].';
        xpre4_30 = [zeros(1,Me-1).';xpre4_10.^3.'].';

        Ap4_10 = zeros(Nn,Me);
        Ap4_01 = zeros(Nn,Me);
        Ap4_03 = zeros(Nn,Me);
        Ap4_12 = zeros(Nn,Me);
        Ap4_21 = zeros(Nn,Me);
        Ap4_30 = zeros(Nn,Me);
        for i = 1:Nn
            Ap4_10(i,:) = fliplr(xpre4_10(i:Me-1+i));
            Ap4_01(i,:) = fliplr(xpre4_01(i:Me-1+i));
            Ap4_03(i,:) = fliplr(xpre4_03(i:Me-1+i));
            Ap4_12(i,:) = fliplr(xpre4_12(i:Me-1+i));
            Ap4_21(i,:) = fliplr(xpre4_21(i:Me-1+i));
            Ap4_30(i,:) = fliplr(xpre4_30(i:Me-1+i));
        end

        Z0mag = pinv(Ap4_10'*Ap4_10)*Ap4_10';
        h0test = Z0mag*y4_awgn(1:Nn);

        X1 = [];
        X1 = cat(2, X1, Ap4_01);
        X1 = cat(2, X1, Ap4_10);

        Z1mag = pinv(X1'*X1)*X1';
        h1est = Z1mag*y4_awgn(1:Nn);
        h1est01 = h1est(0*Me+1:1*Me);
        h1est10 = h1est(1*Me+1:2*Me);

        X2 = [];
        X2 = cat(2, X2, Ap4_01);
        X2 = cat(2, X2, Ap4_10);
        X2 = cat(2, X2, Ap4_03);
        X2 = cat(2, X2, Ap4_12);
        X2 = cat(2, X2, Ap4_21);
        X2 = cat(2, X2, Ap4_30);

        Zmag = pinv(X2'*X2)*X2';
        h2test = Zmag*y4_awgn(1:Nn);
        h2est01 = h2test(0*Me+1:1*Me);
        h2est10 = h2test(1*Me+1:2*Me);
        h2est03 = h2test(2*Me+1:3*Me);
        h2est12 = h2test(3*Me+1:4*Me);
        h2est21 = h2test(4*Me+1:5*Me);
        h2est30 = h2test(5*Me+1:6*Me);

        A4_10 = zeros(K,Me);
        A4_01 = zeros(K,Me);
        A4_03 = zeros(K,Me);
        A4_12 = zeros(K,Me);
        A4_21 = zeros(K,Me);
        A4_30 = zeros(K,Me);
        for i = 1:K
            A4_10(i,:) = fliplr(xp4_10(i:Me-1+i));
            A4_01(i,:) = fliplr(xp4_01(i:Me-1+i));
            A4_03(i,:) = fliplr(xp4_03(i:Me-1+i));
            A4_12(i,:) = fliplr(xp4_12(i:Me-1+i));
            A4_21(i,:) = fliplr(xp4_21(i:Me-1+i));
            A4_30(i,:) = fliplr(xp4_30(i:Me-1+i));
        end

        cancelled_t(:,packet) = A4_10*h0test - y4_awgn;
    end
end

%========================================================================%
%***************************DECODING ROUND*******************************%
%========================================================================%

if(size(strfind(can_meth,'freq'))>0)
    disp('CANCELLATION METHOD - FREQ');
    cancelled_sel = cancelled;
elseif(size(strfind(can_meth,'time'))>0)
    disp('CANCELLATION METHOD - TIME');
    cancelled_sel = cancelled_t;
else
    disp('ERROR - Wrong Cancellation method selected');
end

awgn_noise1 = randn(1e5,1) + 1i*randn(1e5,1);
Pawgn = fFD_Power(awgn_noise1);
Prx = fFD_Power(cancelled_sel(offset:offset+200,packet));
SNR_lin = 10^(SNR(1,packet)/10);
awgn1 = sqrt(Prx/(SNR_lin*Pawgn)).*awgn_noise1;
cancelled_new = zeros(2*length(awgn_noise1)+length(cancelled_sel(:,packet)),1);
for packet = 1:size(cancelled_sel,2)
    cancelled_new(:,packet) = [awgn1;cancelled_sel(:,packet);awgn1];
end
cancelled_new_serie = reshape(cancelled_new, 1, size(cancelled_new,1)*size(cancelled_new,2)).';
[can_pak, fo_sp] = fFD_preamble_detect_final(Fs, cancelled_new_serie, short_preamble, threshold2round);
[can_packet , can_lptot , fo_lp, can_num_pkt] = fFD_symbolalign2(Fs, can_pak, long_preamble);

if(fFD_Power(can_packet)>0)
    disp('RESULT - Packet detected');
    detected = 1;
else
    disp('RESULT - Packet Not detected');
    detected = 0;
end

%NOISE POWER ESTIMATION
SNR_can = zeros(1,min(can_num_pkt,tx_num_pkt));
for packet = 1:can_num_pkt
    can_preamble1 = can_lptot(end-63:end,packet);
    can_preamble2 = can_lptot(end-64-63:end-64,packet);
    can_preamble  = (can_preamble1 + can_preamble2)./2;

    Pn = 10*log10(fFD_Power(rx_packet(4000:5000,packet))) + 30;
    Prx = 10*log10(fFD_Power(can_preamble2) - fFD_Power(can_packet(4000:5000,packet))) + 30;
    SNRest_can = Prx - Pn;
end

% DECODING
[rx_bits,~,~,can_message,~] = ...
    fFD_decode_packet(can_packet , nsamplesXblock, can_num_pkt,...
                       nblocks, n_bpsc, n_cbps, n_dbps, psdu_size);

% BER (BIT ERROR RATE)
BER = [];
for packet = 1:can_num_pkt
    error = find(rx_bits(:,packet)~=tx_bits(:,packet));
    BER = [BER 100*length(error)/n_data];                   %#ok<AGROW>

    fprintf('Pkt %d - Tx Message: %s \n',packet,tx_message{packet});
    fprintf('Pkt %d - Rx Message: %s \n',packet,can_message{packet});
    fprintf('Pkt %d - BER = %d \n',packet,BER(packet));

    if(size(strfind(can_message{packet},message))>0); decoded = 1;
    else; decoded = 0;
    end
end

% % TABLE RESULTS
% T = table((1:min(rx_num_pkt,tx_num_pkt)).',SNRestr.',SNRest_can',...
%     BER.', 'VariableNames',...
%     {'Packet' 'SNR' 'SNR_can' 'BER'})