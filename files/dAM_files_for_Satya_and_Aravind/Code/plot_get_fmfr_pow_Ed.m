clear;
clc;

%% Update these based on real trajectory & optimize freq window
AMfreqEst= load('approx_am_est.mat');
AMfreqEst.time= (1:length(AMfreqEst.AMfreqvec_est))/AMfreqEst.fs;
freq_half_window_co_Hz= 15;

doPlotFig= 1;
doSaveFig= 1;
doSaveData= 1;
doDetrend= 1;

tStart_ms= 0;
tEnd_ms= inf;

freq_low_Hz= 80;
freq_high_Hz= inf;

%%
CodeDir= '..\..\';
addpath(CodeDir);

inDir= '..\inData\';
all_in_files= dir([inDir '*.mat']);
all_in_files= all_in_files(~contains({all_in_files.name}', 'peaks'));

OutFigDir= sprintf('../outData/Figures_Time%.0f_%.0fms__Filter_%.0f_%.0fHz/', tStart_ms, tEnd_ms, freq_low_Hz, freq_high_Hz);
OutDataDir= '../outData/';

if ~isfolder(OutFigDir)
    mkdir(OutFigDir);
end

fs= 24414.0625;
ch1ind.pos= 1;
ch1ind.neg= 3;
ch2ind.pos= 2;
ch2ind.neg= 4;
cbar_range= 20;

all_pow_data= repmat( ...
    struct('filename', '', 'ch1sum_pow', nan, 'ch1diff_pow', nan, 'ch2sum_pow', nan, 'ch2diff_pow', nan, ...
    'ch1sum_pow_frac', nan, 'ch1diff_pow_frac', nan, 'ch2sum_pow_frac', nan, 'ch2diff_pow_frac', nan), ...
    length(all_in_files), 1);

for fileVar=1:length(all_in_files)
    cur_fStruct= all_in_files(fileVar);
    cur_fName_in= [cur_fStruct.folder filesep cur_fStruct.name];
    cur_data= load(cur_fName_in);
    cur_data= cur_data.tdat;
    
    %% read data
    t_data_ms= (1:size(cur_data,2))/fs*1e3;
    am_est_kHz= interp1(AMfreqEst.time(:), AMfreqEst.AMfreqvec_est(:), t_data_ms/1e3)/1e3;
    
    ch1_sum= (cur_data(ch1ind.pos,:)+cur_data(ch1ind.neg,:))/2;
    ch1_diff= (cur_data(ch1ind.pos,:)-cur_data(ch1ind.neg,:))/2;
    
    ch2_sum= (cur_data(ch2ind.pos,:)+cur_data(ch2ind.neg,:))/2;
    ch2_diff= (cur_data(ch2ind.pos,:)-cur_data(ch2ind.neg,:))/2;
    
    valid_time_inds= (t_data_ms>tStart_ms) & (t_data_ms<tEnd_ms);
    
    %% read data in appropriate time window
    t_data_ms= t_data_ms(valid_time_inds);
    am_est_kHz= am_est_kHz(valid_time_inds);
    ch1_sum= ch1_sum(valid_time_inds);
    ch1_diff= ch1_diff(valid_time_inds);
    ch2_sum= ch2_sum(valid_time_inds);
    ch2_diff= ch2_diff(valid_time_inds);
    t_line_ms= (1:length(am_est_kHz))/fs*1e3;
    
    %% detrend
    if doDetrend
        ch1_sum = detrend(ch1_sum);
        ch1_diff = detrend(ch1_diff);
        ch2_sum = detrend(ch2_sum);
        ch2_diff = detrend(ch2_diff);
    end
    
    %% filter data if needed
    applyFilter= 1;
    if (freq_low_Hz>0) && (freq_high_Hz<(fs/2)) % BP
        cur_filt= helper.get_filter_designfilt('bp', [freq_low_Hz, freq_high_Hz], fs);
    elseif (freq_low_Hz>0) && (freq_high_Hz>=(fs/2)) % HP
        cur_filt= helper.get_filter_designfilt('hp', freq_low_Hz, fs);
    elseif (freq_low_Hz<=0) && (freq_high_Hz<(fs/2)) % LP
        cur_filt= helper.get_filter_designfilt('lp', freq_high_Hz, fs);
    else % no filt
        applyFilter= 0;
    end
    
    if applyFilter
        ch1_sum= filtfilt(cur_filt, ch1_sum);
        ch1_diff= filtfilt(cur_filt, ch1_diff);
        ch2_sum= filtfilt(cur_filt, ch2_sum);
        ch2_diff= filtfilt(cur_filt, ch2_diff);
    end
    
    
    %%
    [ch1sum_pow, totPower_ch1, ~, ~]= helper.get_freq_trajectory_power(ch1_sum, fs, am_est_kHz*1e3, 0, freq_half_window_co_Hz);
    [ch1diff_pow, ~, ~, ~]= helper.get_freq_trajectory_power(ch1_diff, fs, am_est_kHz*1e3, 0, freq_half_window_co_Hz);
    [ch2sum_pow, totPower_ch2, ~, ~]= helper.get_freq_trajectory_power(ch2_sum, fs, am_est_kHz*1e3, 0, freq_half_window_co_Hz);
    [ch2diff_pow, ~, ~, ~]= helper.get_freq_trajectory_power(ch2_diff, fs, am_est_kHz*1e3, 0, freq_half_window_co_Hz);
    
    all_pow_data(fileVar)= struct('filename', cur_fStruct.name, ...
        'ch1sum_pow', ch1sum_pow, 'ch1diff_pow', ch1diff_pow, 'ch2sum_pow', ch2sum_pow, 'ch2diff_pow', ch2diff_pow, ...
        'ch1sum_pow_frac', ch1sum_pow/totPower_ch1, 'ch1diff_pow_frac', ch1diff_pow/totPower_ch1, ...
        'ch2sum_pow_frac', ch2sum_pow/totPower_ch2, 'ch2diff_pow_frac', ch2diff_pow/totPower_ch2);
    %%
    
    if doPlotFig
        figSize_cm= [55 5 25 15];
        figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
        figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
        figure(1);
        set(gcf,figure_prop_name,figure_prop_val);
        
        figure(1);
        clf;
        
        sp_ax(1)= subplot(321);
        hold on;
        plot(t_data_ms, ch1_sum);
        plot(t_data_ms, ch1_diff);
        xlabel('Time (ms)')
        ylabel('Amp (au)');
        title(cur_fStruct.name(1:end-4), 'Interpreter', 'none')
        
        sp_ax(2)= subplot(322);
        hold on;
        Pxx= helper.plot_dpss_psd(ch1_sum, fs, 'xunit', 'khz');
        helper.plot_dpss_psd(ch1_diff, fs, 'xunit', 'khz');
        xlim([max(2, .8*freq_low_Hz) 3e3]/1e3)
        ylim([-40 2]+max(Pxx))
        set(gca, 'XTick', [50 100 500 1e3 2e3]/1e3)
        legend('Sum', 'Diff');
        
        sp_ax(3)= subplot(323);
        hold on;
        helper.plot_spectrogram(ch1_sum, fs, [], [], [], cbar_range);
        line(t_line_ms, am_est_kHz, 'color', 'r', 'LineStyle', '--');
        title('SUM | Ch1')
        xlabel(sprintf('Time re. %.0f (ms)', max(0, tStart_ms)))
        
        sp_ax(4)= subplot(324);
        hold on;
        helper.plot_spectrogram(ch1_diff, fs, [], [], [], cbar_range);
        line(t_line_ms, am_est_kHz, 'color', 'r', 'LineStyle', '--');
        title('DIFF | Ch1')
        xlabel(sprintf('Time re. %.0f (ms)', max(0, tStart_ms)))
        
        sp_ax(5)= subplot(325);
        hold on;
        helper.plot_spectrogram(ch2_sum, fs, [], [], [], cbar_range);
        line(t_line_ms, am_est_kHz, 'color', 'r', 'LineStyle', '--');
        title('SUM | Ch2')
        xlabel(sprintf('Time re. %.0f (ms)', max(0, tStart_ms)))
        
        sp_ax(6)= subplot(326);
        hold on;
        helper.plot_spectrogram(ch2_diff, fs, [], [], [], cbar_range);
        line(t_line_ms, am_est_kHz, 'color', 'r', 'LineStyle', '--');
        title('DIFF | Ch2')
        xlabel(sprintf('Time re. %.0f (ms)', max(0, tStart_ms)))
        
        linkaxes(sp_ax(3:6), 'xy')
        ylim(sp_ax(3), [0 1.2])
        
        %%
        xc=.05;
        xw= .4;
        xs= .08;
        yc= .07;
        ys= .1;
        yw= .23;
        
        set(sp_ax(5), 'Units', 'normalized', 'Position', [xc, yc, xw, yw]);
        set(sp_ax(6), 'Units', 'normalized', 'Position', [xc+xs+xw, yc, xw, yw]);
        set(sp_ax(3), 'Units', 'normalized', 'Position', [xc, yc+ys+yw, xw, yw]);
        set(sp_ax(4), 'Units', 'normalized', 'Position', [xc+xs+xw, yc+ys+yw, xw, yw]);
        set(sp_ax(1), 'Units', 'normalized', 'Position', [xc, yc+2*ys+2*yw, xw, yw]);
        set(sp_ax(2), 'Units', 'normalized', 'Position', [xc+xs+xw, yc+2*ys+2*yw, xw, yw]);
        %%
        fig_name= sprintf('%sFig_%s', OutFigDir, cur_fStruct.name(1:end-4));
        if doSaveFig
            print(fig_name, '-dpng', '-r300')
        end
    end
end

%%
out_fName= sprintf('%sPowerValuesTime%.0f_%.0fms__Filter_%.0f_%.0fHz.mat', OutDataDir, tStart_ms, tEnd_ms, freq_low_Hz, freq_high_Hz);
if doSaveData
    save(out_fName, 'all_pow_data');
end