clear;
clc;

do_save_fig= 0;

code_dir= fileparts(pwd);
addpath(code_dir)

%% with actual trajectory
fs= 20e3;
stimdur= 1;
stimlen= round(fs*stimdur);
dam_f_start= 2^4;
dam_f_end= 2^10.25;

fc= 3e3;
mod_depth= 1;
tStim= (1:stimlen)/stimlen*stimdur;

dam_traj_Hz= 2.^linspace(log2(dam_f_start), log2(dam_f_end), stimlen);

phi_instantaneous= -cumtrapz(dam_traj_Hz)/fs; % Formula for phase in terms of frequency
stimMod= cos(2*pi*phi_instantaneous);
stimCar= sin(2*pi*fc*tStim);
stim= (1-mod_depth*stimMod(:))/2.*stimCar(:);
stim_env= abs(hilbert(stim));
stim_env= stim_env-mean(stim_env);

figSize_cm= [7 3 17 8];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(4);
set(gcf,figure_prop_name,figure_prop_val);

clf;

nSProws= 2;
nSPcols= 3;

sp_ax= nan(nSProws*nSPcols, 1);
for sp_var=1:nSProws*nSPcols
    sp_ax(sp_var)= subplot(nSProws, nSPcols, sp_var);
end

%%
all_time_windows_ms= [10, 20, 30, 50, 75, 100, 200];
example_time_windows_ms= [30, 100];
window_dur_tick= [10 30 50 100 200];
STres_tick= [1 3 10 30 100 300];

total_ST_resolution= nan(length(all_time_windows_ms), 1);
total_power_captured= nan(length(all_time_windows_ms), 1);
txt_han= nan(length(all_time_windows_ms), 1);

for win_var=1:length(all_time_windows_ms)
    tWindow_ms= all_time_windows_ms(win_var);

    if ismember(tWindow_ms, example_time_windows_ms)
        axes(sp_ax(tWindow_ms==example_time_windows_ms))
        do_plot_boxes= 1;

        hold on
        helper.plot_spectrogram(stim_env, fs, 40e-3, .95)
        xlim([-50 max(tStim)*1e3+50]);
        colorbar off;
        ylim([0 dam_f_end*1.2]/1e3);

    else
        do_plot_boxes= 0;
    end

    [total_ST_resolution(win_var), total_power_captured(win_var), txt_han(win_var)]= helper.get_STresolution(stim_env, fs, dam_traj_Hz, tWindow_ms/1e3, do_plot_boxes);

end

axes(sp_ax(3))
mrk_alpha= .75;
scatter(all_time_windows_ms, total_ST_resolution, 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, 'MarkerFaceColor', helper.get_color('gray'), 'MarkerEdgeColor', helper.get_color('gray'));
xlim([.9*min(all_time_windows_ms) 1.1*max(all_time_windows_ms)]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', window_dur_tick, 'YTick', STres_tick);
ylim([1.5 1.2*max(total_ST_resolution)]);
xlabel('Window duration (ms)');
ylabel('ST resolution');

[best_st_res, best_st_ind]= min(total_ST_resolution);
best_st_window= all_time_windows_ms(best_st_ind);

%%
all_filt_halfBW= [.5 1:0.25:2 2:4];
all_power_powFun_prct= nan(size(all_filt_halfBW));
plotPSD= 0;

for filtVar=1:length(all_filt_halfBW)
    cur_filtHalfWidth= all_filt_halfBW(filtVar);

    %
    [outPower, totPower, Pxx_dB, Freq_xx]= helper.get_freq_trajectory_power(stim_env(:), fs, dam_traj_Hz, plotPSD, cur_filtHalfWidth, fs/4);
    all_power_powFun_prct(filtVar)= 100*outPower/totPower;

end
optimal_ST_res= 2*interp1(all_power_powFun_prct(:)+.00001*randn(size(all_power_powFun_prct(:))), all_filt_halfBW(:), 95);

axes(sp_ax(4));
hold on
helper.plot_spectrogram(stim_env, fs, 40e-3, .95)
xlim([-50 max(tStim)*1e3+50]);
line(tStim*1e3, dam_traj_Hz/1e3, 'color', helper.get_color('r'), 'linew', 1.5, 'linestyle', '--')
colorbar off;
ylim([0 dam_f_end*1.2]/1e3);

axes(sp_ax(5));
hold on;
plot((Freq_xx-fs/4), Pxx_dB, 'Color', helper.get_color('b'), 'LineWidth', 1.5);
plot(optimal_ST_res*[-0.5 0.5], max(Pxx_dB)+[2 2], 'Color', helper.get_color('r'), 'LineWidth', 1.5)
xlim([-20 20])
ylim([-50 +5]+max(Pxx_dB));
xlabel('Frequency (Hz)');
ylabel('DFT amp (dB)');

axes(sp_ax(6));
plot(2*all_filt_halfBW, all_power_powFun_prct, 'x', 'Color', helper.get_color('r'), 'LineWidth', 1.5);
ylim([0 105]);
xlabel('LP filter bandwidth (Hz)');
ylabel('% power')

axes(sp_ax(3));
hold on
plot([min(all_time_windows_ms), max(all_time_windows_ms)], optimal_ST_res*ones(1,2), 'LineStyle', '-', 'Color', helper.get_color('r'), 'LineWidth', 2)

fprintf('Traditional best = %.1f (for %.0f ms) | Spectrally specific = %.1f \n', best_st_res, best_st_window, optimal_ST_res)

set(findall(gcf,'-property','FontSize'),'FontSize', 8);
set(txt_han(~isnan(txt_han)),'FontSize', 11);
set(findall(gcf,'-property','box'),'box', 'off');
set(findall(gcf,'-property','TickDir'),'TickDir', 'both');

txtHan= helper.add_subplot_letter(2, 3, 11, -.17, 1.12,  'abcdef');

%%
xc= .055;
xw= .27;
xs= .06;
yc= .11;
yw= .35;
ys= .13;

set(sp_ax(4), 'Units', 'normalized', 'Position', [xc, yc, xw, yw])
set(sp_ax(5), 'Units', 'normalized', 'Position', [xc+xw+xs, yc, xw, yw])
set(sp_ax(6), 'Units', 'normalized', 'Position', [xc+2*xw+2*xs, yc, xw, yw])

set(sp_ax(1), 'Units', 'normalized', 'Position', [xc, yc+yw+ys, xw, yw])
set(sp_ax(2), 'Units', 'normalized', 'Position', [xc+xw+xs, yc+yw+ys, xw, yw])
set(sp_ax(3), 'Units', 'normalized', 'Position', [xc+2*xw+2*xs, yc+yw+ys, xw, yw])

%%
fig_dir= ['figures' filesep];
fig_name= [fig_dir 'Fig2_STres_SpctSpcfc'];
if do_save_fig
    print(fig_name, '-dpng', '-r600')
end