% function [power_SUM_FMFR, power_DIFF_FMFR, fm_depth]= get_fmfr_power(file2use, FM_carrier, FMrate)
% Assumes
function [power_SUM_FMFR, power_DIFF_FMFR, fm_depth]= get_fmfr_power(fmfr_filename, FM_carrier, FMrate, doSaveMAT)

%% Important parameters
mat_filename= ['DataOut' filesep fmfr_filename '_PowerValues.mat'];
flag_usePMTM1_timedomain0= 1;
filtWidth_Hz= 5; % if left empty, defaults to 2/stim_duration
FMFR_filter_BW= [300 2.6e3];

if ~exist('doSaveMAT', 'var')
    doSaveMAT= 1;
end

%% Load data
fmfr_data= load(['data' filesep fmfr_filename '.mat']);
fs= fmfr_data.fs;
fm_depth= fmfr_data.fmdepthstemp(:)';

hpFilter= get_filter_designfilt('bp', FMFR_filter_BW, fs); % filter for data

power_SUM_FMFR= nan(size(fm_depth));
power_DIFF_FMFR= nan(size(fm_depth));

plotPSD= 0;

for fmVar= 1:length(fm_depth)
    cur_fm= fm_depth(fmVar);
    
    cur_sum_data= fmfr_data.tdatsumall(fmVar, :);
    cur_diff_data= fmfr_data.tdatdiffall(fmVar, :);
    
    % filter
    cur_sum_data= filtfilt(hpFilter, cur_sum_data);
    cur_diff_data= filtfilt(hpFilter, cur_diff_data);
    
    tVec= (1:length(cur_diff_data))/fs;
    cur_fm_traj= FM_carrier + cur_fm*sin(2*pi*FMrate*tVec + pi/2);
    
    if flag_usePMTM1_timedomain0==1
        [sum_power_2Fc, ~]= helper.get_freq_trajectory_power(cur_sum_data, fs, 2*cur_fm_traj, plotPSD, filtWidth_Hz);
        
        [diff_power_1Fc, ~]= helper.get_freq_trajectory_power(cur_diff_data, fs, cur_fm_traj, plotPSD, filtWidth_Hz);
        [diff_power_3Fc, ~]= helper.get_freq_trajectory_power(cur_diff_data, fs, 3*cur_fm_traj, plotPSD, filtWidth_Hz);
        
        power_SUM_FMFR(fmVar)= sum_power_2Fc;
        power_DIFF_FMFR(fmVar)= diff_power_1Fc + diff_power_3Fc;
        
    elseif flag_usePMTM1_timedomain0==0
        
        [~, filtSignal_sum_2Fc]= helper.get_trajectory_signal(cur_sum_data, fs, 2*cur_fm_traj, filtWidth_Hz);
        
        [~, filtSignal_diff_1Fc]= helper.get_trajectory_signal(cur_diff_data, fs, cur_fm_traj, filtWidth_Hz);
        [~, filtSignal_diff_3Fc]= helper.get_trajectory_signal(cur_diff_data, fs, 3*cur_fm_traj, filtWidth_Hz);
        
        power_SUM_FMFR(fmVar)= var(filtSignal_sum_2Fc);
        power_DIFF_FMFR(fmVar)= var(filtSignal_diff_1Fc) + var(filtSignal_diff_3Fc);
    end
end

%%
fmfr_data= struct('FMFR_filename', fmfr_filename, 'SUM_power', power_SUM_FMFR, 'DIFF_power', power_DIFF_FMFR, 'FM_depth', fm_depth, 'FM_carrier', FM_carrier, 'FMrate', FMrate);
if doSaveMAT
    save(mat_filename, 'fmfr_data');
end