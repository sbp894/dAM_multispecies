<<<<<<< HEAD
% This script quantifies the effect of filter BW on freq-demod Gerbil data 
clear;
clc;

%% 
% 'FMFR_GER36_unfilt' is a really good/clean example to look at the effect of filter BW
doSaveFig= 0;
addpath('..\');
all_files= dir(['GrandAvg' filesep '*.mat']);
all_filt_halfBW= [2:20 30 40 50];
print_freq= 12;
print_ind= dsearchn(all_filt_halfBW(:), print_freq);

for fileVar= 1:length(all_files)
    file2use= all_files(fileVar);
    
    fmfr_data= load(['GrandAvg' filesep file2use.name]);
    
    f_traj= fmfr_data.dam_traj_Hz;
    valid_inds= ~isnan(f_traj);
    f_traj= f_traj(valid_inds);
    y_t= fmfr_data.mean_data(valid_inds);
    fs= fmfr_data.fs_ffr;
    t= (1:length(y_t))/fs;
    dur= max(t);
    
    figSize_cm= [5 5 25 10];
    figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
    figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
    figure(fileVar);
    set(gcf,figure_prop_name,figure_prop_val);
    clf;
    
    sp_ax(1)= subplot(2,3,1:2);
    hold on;
    helper.plot_spectrogram(y_t, fs);
    colorbar off;
    ylim([.4 2.5]);
    % ylim([.4 .6]);
    line(t*1e3, f_traj/1e3, 'color', 'r', 'linestyle', '--')
    line(t*1e3, 3*f_traj/1e3, 'color', 'r', 'linestyle', '--')
    
    all_power_sigFun_prct= nan(size(all_filt_halfBW));
    all_power_powFun_prct= nan(size(all_filt_halfBW));
    
    sp_ax(2)= subplot(2,3,4:5);
    
    for filtVar=1:length(all_filt_halfBW)
        cur_filtHalfWidth= all_filt_halfBW(filtVar);
        [fracSignal_1, filtSignal_1]= helper.get_trajectory_signal(y_t(:), fs, f_traj, cur_filtHalfWidth);
        [fracSignal_3, filtSignal_3]= helper.get_trajectory_signal(y_t(:), fs, 3*f_traj, cur_filtHalfWidth);
        
        %
        plotPSD= 0;
        
        [outPower_1, ~]= helper.get_freq_trajectory_power(y_t(:), fs, f_traj, plotPSD, cur_filtHalfWidth);
        [outPower_3, totPower]= helper.get_freq_trajectory_power(y_t(:), fs, 3*f_traj, plotPSD, cur_filtHalfWidth);
        %     all_power_sigFun(filtVar)= rms(filtSignal).^2/rms(y_tr).^2;
        all_power_sigFun_prct(filtVar)= 100*(rms(fracSignal_1).^2+rms(fracSignal_3).^2);
        all_power_powFun_prct(filtVar)= 100*(outPower_1/totPower + outPower_3/totPower);
        
        subplot(2,3,4:5);
        hold on;
        plot(t*1e3, fracSignal_1)
        
    end
    axes(sp_ax(2));
    xlabel('Time (ms)');
    ylabel('Fractional signal');
    ylim([0 1.2]);
    
    
    figure(fileVar);
    sp_ax(3)= subplot(133);
    hold on;
    plot(1:length(all_filt_halfBW), all_power_sigFun_prct, 'd', 'MarkerSize', 8, 'linew', 1.5);
    plot(1:length(all_filt_halfBW), all_power_powFun_prct, 'x', 'MarkerSize', 8, 'linew', 1.5);
    plot([1 1]*2/dur, ylim(), 'm--');
    set(gca, 'xtick', 1:length(all_filt_halfBW), 'XTickLabel', num2str(all_filt_halfBW(:)));
    xlabel('Filter half-BW');
    ylabel('Power captured (%)');
    title([file2use.name '|% @12Hz=' num2str(round(all_power_powFun_prct(print_ind)))], 'Interpreter', 'none');
    
    legend('Using time-domain filter', 'using multi-taper PSD', 'Min theoretical BW', 'Location', 'best');
    linkaxes(sp_ax(1:2), 'x')
    xlim(sp_ax(1), [-1 dur*1e3+5])
    set(findall(gcf,'-property','FontSize'),'FontSize', 11);
    
    %%
    spa.xc= .055;
    spa.xw= .42;
    spa.xs= .08;
    spa.yc= .125;
    spa.yw= .34;
    spa.ys= .13;
    
    set(sp_ax(2), 'Units', 'normalized', 'Position', [spa.xc, spa.yc, spa.xw, spa.yw]);
    set(sp_ax(1), 'Units', 'normalized', 'Position', [spa.xc, spa.yc+spa.ys+spa.yw, spa.xw, spa.yw]);
    set(sp_ax(3), 'Units', 'normalized', 'Position', [spa.xc+spa.xs+spa.xw, spa.yc, spa.xw, 2*spa.yw+1*spa.ys]);
    
    
    if doSaveFig
        figName2Save= sprintf('FigOut%sEffectofBW%s%s_FM%.0fHz', filesep, filesep, file2use, f_depth);
        print(figName2Save, '-dpng', '-r600')
    end
end
=======
% This script quantifies the effect of filter BW on freq-demod Gerbil data 
clear;
clc;

%% 
% 'FMFR_GER36_unfilt' is a really good/clean example to look at the effect of filter BW
doSaveFig= 0;
addpath('..\');
all_files= dir(['GrandAvg' filesep '*.mat']);
all_filt_halfBW= [2:20 30 40 50];
print_freq= 12;
print_ind= dsearchn(all_filt_halfBW(:), print_freq);

for fileVar= 1:length(all_files)
    file2use= all_files(fileVar);
    
    fmfr_data= load(['GrandAvg' filesep file2use.name]);
    
    f_traj= fmfr_data.dam_traj_Hz;
    valid_inds= ~isnan(f_traj);
    f_traj= f_traj(valid_inds);
    y_t= fmfr_data.mean_data(valid_inds);
    fs= fmfr_data.fs_ffr;
    t= (1:length(y_t))/fs;
    dur= max(t);
    
    figSize_cm= [5 5 25 10];
    figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
    figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
    figure(fileVar);
    set(gcf,figure_prop_name,figure_prop_val);
    clf;
    
    sp_ax(1)= subplot(2,3,1:2);
    hold on;
    helper.plot_spectrogram(y_t, fs);
    colorbar off;
    ylim([.4 2.5]);
    % ylim([.4 .6]);
    line(t*1e3, f_traj/1e3, 'color', 'r', 'linestyle', '--')
    line(t*1e3, 3*f_traj/1e3, 'color', 'r', 'linestyle', '--')
    
    all_power_sigFun_prct= nan(size(all_filt_halfBW));
    all_power_powFun_prct= nan(size(all_filt_halfBW));
    
    sp_ax(2)= subplot(2,3,4:5);
    
    for filtVar=1:length(all_filt_halfBW)
        cur_filtHalfWidth= all_filt_halfBW(filtVar);
        [fracSignal_1, filtSignal_1]= helper.get_trajectory_signal(y_t(:), fs, f_traj, cur_filtHalfWidth);
        [fracSignal_3, filtSignal_3]= helper.get_trajectory_signal(y_t(:), fs, 3*f_traj, cur_filtHalfWidth);
        
        %
        plotPSD= 0;
        
        [outPower_1, ~]= helper.get_freq_trajectory_power(y_t(:), fs, f_traj, plotPSD, cur_filtHalfWidth);
        [outPower_3, totPower]= helper.get_freq_trajectory_power(y_t(:), fs, 3*f_traj, plotPSD, cur_filtHalfWidth);
        %     all_power_sigFun(filtVar)= rms(filtSignal).^2/rms(y_tr).^2;
        all_power_sigFun_prct(filtVar)= 100*(rms(fracSignal_1).^2+rms(fracSignal_3).^2);
        all_power_powFun_prct(filtVar)= 100*(outPower_1/totPower + outPower_3/totPower);
        
        subplot(2,3,4:5);
        hold on;
        plot(t*1e3, fracSignal_1)
        
    end
    axes(sp_ax(2));
    xlabel('Time (ms)');
    ylabel('Fractional signal');
    ylim([0 1.2]);
    
    
    figure(fileVar);
    sp_ax(3)= subplot(133);
    hold on;
    plot(1:length(all_filt_halfBW), all_power_sigFun_prct, 'd', 'MarkerSize', 8, 'linew', 1.5);
    plot(1:length(all_filt_halfBW), all_power_powFun_prct, 'x', 'MarkerSize', 8, 'linew', 1.5);
    plot([1 1]*2/dur, ylim(), 'm--');
    set(gca, 'xtick', 1:length(all_filt_halfBW), 'XTickLabel', num2str(all_filt_halfBW(:)));
    xlabel('Filter half-BW');
    ylabel('Power captured (%)');
    title([file2use.name '|% @12Hz=' num2str(round(all_power_powFun_prct(print_ind)))], 'Interpreter', 'none');
    
    legend('Using time-domain filter', 'using multi-taper PSD', 'Min theoretical BW', 'Location', 'best');
    linkaxes(sp_ax(1:2), 'x')
    xlim(sp_ax(1), [-1 dur*1e3+5])
    set(findall(gcf,'-property','FontSize'),'FontSize', 11);
    
    %%
    spa.xc= .055;
    spa.xw= .42;
    spa.xs= .08;
    spa.yc= .125;
    spa.yw= .34;
    spa.ys= .13;
    
    set(sp_ax(2), 'Units', 'normalized', 'Position', [spa.xc, spa.yc, spa.xw, spa.yw]);
    set(sp_ax(1), 'Units', 'normalized', 'Position', [spa.xc, spa.yc+spa.ys+spa.yw, spa.xw, spa.yw]);
    set(sp_ax(3), 'Units', 'normalized', 'Position', [spa.xc+spa.xs+spa.xw, spa.yc, spa.xw, 2*spa.yw+1*spa.ys]);
    
    
    if doSaveFig
        figName2Save= sprintf('FigOut%sEffectofBW%s%s_FM%.0fHz', filesep, filesep, file2use, f_depth);
        print(figName2Save, '-dpng', '-r600')
    end
end
>>>>>>> 033e1e8208ba0bc9e95b7a551303ab1999f126dc
rmpath('..\');