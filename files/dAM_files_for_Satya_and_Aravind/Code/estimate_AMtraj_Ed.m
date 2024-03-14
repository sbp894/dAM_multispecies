clear;
clc;

AMfreqEst= load('..\..\AMfreqvec_est.mat');
AMfreqEst.time= (1:length(AMfreqEst.AMfreqvec_est))/AMfreqEst.fs;

fs=48828.125; %TDT sampling rate
T_half= 0.5;
tvec=1/fs:1/fs:1; %time in s
tvec2=1/fs:1/fs:.5; %duration of each piece
AMmodfq1(1:length(tvec2))=linspace(4,6.5,length(tvec2));
AMmodfq1(length(tvec2)+1:length(tvec))=linspace(6.5,10.25,length(tvec2));
AMfreqvec=2.^AMmodfq1;


%%

stim_data= load('DOD_sweeps_Nov2021.mat');
fs_stim= numel(stim_data.AMmod1);
tStim_ms= (1:length(stim_data.amfmtonevec))/fs_stim*1e3;

figure(2)
clf;
subplot(121);
hold on;
plot_spectrogram(detrend(stim_data.AMmod1), fs_stim)
line(AMfreqEst.time*1e3, AMfreqEst.AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
title('AMmod1');

subplot(122);
hold on;
plot_spectrogram(detrend(stim_data.AMFMtoneup_piece), fs_stim)
line(AMfreqEst.time*1e3, 8+AMfreqEst.AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
line(AMfreqEst.time*1e3, 8-AMfreqEst.AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
title('AMFMtoneup_piece');
% line(AMfreqEst.time*1e3, 8+AMfreqEst.AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
% line(AMfreqEst.time*1e3, 8-AMfreqEst.AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
% title('amfmtonevec');

%%
AMfreqvec_est_1= AMfreqvec(1:length(tvec2))/(2*pi) .* ( 1 + log(2)*(6.5-4)/T_half*tvec2)*3;
AMfreqvec_est_2= AMfreqvec(length(tvec2)+1:length(tvec))/(2*pi) .* ( log(2)*2.4/T_half*tvec(length(tvec2)+1:length(tvec))).^2.22;
AMfreqvec_est= [AMfreqvec_est_1, AMfreqvec_est_2];


figure(3);
clf;
hold on;
plot_spectrogram(detrend(stim_data.AMFMtoneup_piece), fs_stim, 40e-3)
line(tvec*1e3, 8+AMfreqvec_est/1e3, 'color', 'r', 'LineStyle', '--');
ylim([6 10])

save('approx_am_est.mat', 'AMfreqvec_est', 'fs')
