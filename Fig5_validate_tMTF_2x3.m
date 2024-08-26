clear;
clc;

%%
do_save_fig= 0;

nSProws= 2;
nSPcols= 3;
sp_ax= get_axes(5, nSProws, nSPcols);

% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\dAM_power_data\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\0ms_12Hz\';
Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\5ms_12Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\10ms_12Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\40ms_12Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\40ms_20Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\25ms_20Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\25ms_8Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\25ms_5Hz\';


%% 

[human_dAM_data, human_SAM_data, lin_mdl_human] = plot_human_sAM_dAM_data(sp_ax([1, 6]), Root_dAM_power_dir);

[gerbil_dAM_data, gerbil_SAM_data, lin_mdl_gerbil] = plot_gerbil_sAM_dAM_data(sp_ax([2, 6]), Root_dAM_power_dir);

[mice_dAM_data, mice_SAM_data, lin_mdl_mice] = plot_mice_sAM_dAM_data(sp_ax([3, 6]), Root_dAM_power_dir);

[rat_vert_dAM_data, rat_vert_SAM_data, lin_mdl_rat_vert] = plot_rat_sAM_dAM_data(sp_ax([4, 6]), 'vertical', Root_dAM_power_dir);

[rat_horz_dAM_data, rat_horz_SAM_data, lin_mdl_rat_horz] = plot_rat_sAM_dAM_data(sp_ax([5, 6]), 'horizontal', Root_dAM_power_dir);

set(findall(gcf,'-property','FontSize'),'FontSize', 8);
set(findall(gcf,'-property','TickDir'),'TickDir', 'both');
set(findall(gcf,'-property','box'),'box', 'off');

%%
SP_letters= 'ABCDEF';
panel_han= nan(nSProws*nSPcols, 1);
for sp_var=1:nSProws*nSPcols
    axes(sp_ax(sp_var));
    panel_han(sp_var)= text(-.05, 1.1, SP_letters(sp_var), 'Units', 'normalized', 'FontWeight','bold', 'FontSize', 11);
end

%% 
% overall Rsq 
% all_dam_data= [normalize(human_dAM_data.mu_young), normalize(gerbil_dAM_data.mu_19wk), normalize(rat_vert_dAM_data.mu_control), ...
%     normalize(rat_horz_dAM_data.mu_control), normalize(mice_dAM_data.mu_control)];
% all_SAM_data= [normalize(human_SAM_data.mu_young), normalize(gerbil_SAM_data.mu_19wk), normalize(rat_vert_SAM_data.mu_control), ...
%     normalize(rat_horz_SAM_data.mu_control), normalize(mice_SAM_data.mu_control)];
all_dam_data= [normalize(human_dAM_data.mu_young), normalize(gerbil_dAM_data.mu_19wk), normalize(rat_vert_dAM_data.mu_control), ...
    normalize(rat_horz_dAM_data.mu_control), normalize(mice_dAM_data.mu_control)];
all_SAM_data= [normalize(human_SAM_data.mu_young), normalize(gerbil_SAM_data.mu_19wk), normalize(rat_vert_SAM_data.mu_control), ...
    normalize(rat_horz_SAM_data.mu_control), normalize(mice_SAM_data.mu_control)];

lin_mdl_all= fitlm(all_SAM_data, all_dam_data);
lin_fit_x= [min(all_SAM_data); max(all_SAM_data)];
lin_fit_y= predict(lin_mdl_all, lin_fit_x);

mrk_alpha= .75;

axes(sp_ax(6));
hold on
lHan(1)= scatter(normalize(human_SAM_data.mu_young), normalize(human_dAM_data.mu_young), 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, ...
    'MarkerFaceColor', get_color('m'), 'DisplayName', 'Human');
lHan(2)= scatter(normalize(gerbil_SAM_data.mu_19wk), normalize(gerbil_dAM_data.mu_19wk), 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, ...
    'MarkerFaceColor', get_color('b'), 'DisplayName', 'Gerbil');
lHan(3)= scatter(normalize(rat_vert_SAM_data.mu_control), normalize(rat_vert_dAM_data.mu_control), 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, ...
    'MarkerFaceColor', get_color('r'), 'DisplayName', 'Rat (vertical)');
lHan(4)= scatter(normalize(rat_horz_SAM_data.mu_control), normalize(rat_horz_dAM_data.mu_control), 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, ...
    'MarkerFaceColor', get_color('prp'), 'DisplayName', 'Rat (horizontal)');
lHan(5)= scatter(normalize(mice_SAM_data.mu_control), normalize(mice_dAM_data.mu_control), 'filled', 'MarkerFaceAlpha', mrk_alpha, 'MarkerEdgeAlpha', mrk_alpha, ...
    'MarkerFaceColor', get_color('g'), 'DisplayName', 'Mice');
plot(lin_fit_x, lin_fit_y, 'k-', 'LineWidth', 2)
text(.05, .9, sprintf('$p<10^{-16}$'), 'Units', 'normalized', 'Interpreter','latex')
text(.05, .75, sprintf('$R^2=%.2f$', lin_mdl_all.Rsquared.Adjusted), 'Units', 'normalized', 'Interpreter','latex')

% legend(lHan, {'Human', 'Gerbil', 'Rat (vertical)', 'Rat (horizontal)', 'Mice'}, 'Location', 'southeast', 'box', 'off');
xlabel('Transformed SAM power (au)');
ylabel('Transformed dAM power (au)');

%% axes tuning 
linkaxes(sp_ax(1:5), 'x');
xlim(sp_ax(1), [12 1.25e3])


axes(sp_ax(1));
legend('dAM', 'SAM', 'box', 'off', 'Location', 'west')
yyaxis left;
ylabel('dAM power (dB)')

axes(sp_ax(3));
yyaxis right;
ylabel('SAM power (dB)')

axes(sp_ax(4));
yyaxis left;
ylabel('dAM power (dB)')

axes(sp_ax(5));
yyaxis right;
ylabel('SAM power (dB)')


%%
if strcmp(getenv('COMPUTERNAME'), 'NB-VATS-LAB1')
    fig_dir= 'C:\Users\sap245\Google Drive\PostDoc\dAM_MultiSpecies\Figures\';
else
    fig_dir= 'G:\My Drive\PostDoc\dAM_MultiSpecies\Figures\';
end

fig_name= [fig_dir 'Fig5_validate'];
if do_save_fig
    print(fig_name, '-dpng', '-r600')
end


%% print all Rsq 
figure(98)
plot([lin_mdl_human.Rsquared.Adjusted, lin_mdl_gerbil.Rsquared.Adjusted, lin_mdl_mice.Rsquared.Adjusted, lin_mdl_rat_vert.Rsquared.Adjusted, lin_mdl_rat_horz.Rsquared.Adjusted, lin_mdl_all.Rsquared.Adjusted,], 'd-');
ylim([0, 1]);
set(gca, 'XTick', 1:6, 'XTickLabel', {'human', 'gerbil', 'mice', 'rat-V', 'rat-H', 'all'})
%% 
function [gerbil_dAM_data, gerbil_SAM_data, lin_mdl] = plot_gerbil_sAM_dAM_data(sp_ax, Root_dAM_power_dir)
freq_tick= [16 50 200 800];
lw1= 1;
lw2= 1.25;

% dAM data 
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\dAM_power_data\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\';
gerbil_fName= [Root_dAM_power_dir 'gerbil_three_group_power.mat'];
gerbil_dAM_data= load(gerbil_fName);
AMfreqs_Hz= gerbil_dAM_data.AM_Hz_2use;
if ~gerbil_dAM_data.use_dB
    error('Need dB, but not dB');
end
gerbil_dAM_data.mu_19wk= nanmean(gerbil_dAM_data.pow_val_young_frac);
gerbil_dAM_data.sem_19wk= nanstd(gerbil_dAM_data.pow_val_young_frac)/sqrt(size(gerbil_dAM_data.pow_val_young_frac,2));
gerbil_dAM_data.mu_44wk= nanmean(gerbil_dAM_data.pow_val_midaged_frac);
gerbil_dAM_data.sem_44wk= nanstd(gerbil_dAM_data.pow_val_midaged_frac)/sqrt(size(gerbil_dAM_data.pow_val_midaged_frac,2));
gerbil_dAM_data.mu_75wk= nanmean(gerbil_dAM_data.pow_val_old_frac);
gerbil_dAM_data.sem_75wk= nanstd(gerbil_dAM_data.pow_val_old_frac)/sqrt(size(gerbil_dAM_data.pow_val_old_frac,2));

% Discrete data 
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\Discrete_FFTdata\';

discrete_freq_data_19wk= readtable([Discrete_FFT_Dir 'Gerbilefr3k_19wk.xlsx']);
discrete_freqs_19wk= discrete_freq_data_19wk.Level;
discrete_freq_data_19wk= table2array(discrete_freq_data_19wk);
discrete_freq_data_19wk= discrete_freq_data_19wk(ismember(discrete_freqs_19wk, AMfreqs_Hz), 2:end);
discrete_freq_data_19wk= discrete_freq_data_19wk'; % make columns ~ freq
discrete_freq_data_19wk= mag2db(discrete_freq_data_19wk);

discrete_freq_data_44wk= readtable([Discrete_FFT_Dir 'Gerbilefr3k_44wk.xlsx']);
discrete_freqs_44wk= discrete_freq_data_44wk.Level;
discrete_freq_data_44wk= table2array(discrete_freq_data_44wk);
discrete_freq_data_44wk= discrete_freq_data_44wk(ismember(discrete_freqs_44wk, AMfreqs_Hz), 2:end);
discrete_freq_data_44wk= discrete_freq_data_44wk'; % make columns ~ freq
discrete_freq_data_44wk= mag2db(discrete_freq_data_44wk);

discrete_freq_data_75wk= readtable([Discrete_FFT_Dir 'Gerbilefr3k_75wk.xlsx']);
discrete_freqs_75wk= discrete_freq_data_75wk.Level;
discrete_freq_data_75wk= table2array(discrete_freq_data_75wk);
discrete_freq_data_75wk= discrete_freq_data_75wk(ismember(discrete_freqs_44wk, AMfreqs_Hz), 2:end);
discrete_freq_data_75wk= discrete_freq_data_75wk'; % make columns ~ freq
discrete_freq_data_75wk= mag2db(discrete_freq_data_75wk);

if ~isequal(discrete_freqs_19wk, discrete_freqs_44wk) & ~isequal(discrete_freqs_19wk, discrete_freqs_75wk) 
    error('AM freqs differ');
end

gerbil_SAM_data= struct(...
    'discrete_freq_data_19wk', discrete_freq_data_19wk, 'discrete_freq_data_44wk', discrete_freq_data_44wk, 'discrete_freq_data_75wk', discrete_freq_data_75wk, 'AM_Hz_2use', AMfreqs_Hz, ...
    'mu_19wk', nanmean(discrete_freq_data_19wk), 'mu_44wk', nanmean(discrete_freq_data_44wk), 'mu_75wk', nanmean(discrete_freq_data_75wk), ...
    'sem_19wk', nanstd(discrete_freq_data_19wk)/sqrt(size(discrete_freq_data_19wk,2)), 'sem_44wk', nanstd(discrete_freq_data_44wk)/sqrt(size(discrete_freq_data_44wk,2)), 'sem_75wk', nanstd(discrete_freq_data_75wk)/sqrt(size(discrete_freq_data_75wk,2)));

xlim_val= [.9*min(gerbil_dAM_data.AM_Hz_2use), 1.1*max(gerbil_dAM_data.AM_Hz_2use)];

axes(sp_ax(1));
yyaxis left;
hold on;
lHan(1)= errorbar(gerbil_dAM_data.AM_Hz_2use, gerbil_dAM_data.mu_19wk, gerbil_dAM_data.sem_19wk, 'LineWidth', lw2, 'Color', get_color('b'));
% lHan(2)= errorbar(gerbil_dAM_data.AM_Hz_2use, gerbil_dAM_data.mu_44wk, gerbil_dAM_data.sem_44wk, 'LineWidth', lw2, 'Color', get_color('b'));
% lHan(3)= errorbar(gerbil_dAM_data.AM_Hz_2use, gerbil_dAM_data.mu_75wk, gerbil_dAM_data.sem_75wk, 'LineWidth', lw2, 'Color', get_color('lb'));
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
text(.05, .1, 'Gerbil', 'Color', get_color('b'), 'Units', 'normalized');
xlim(xlim_val)
% legend(lHan, {'19 week', '45 week', '75 week'}, 'Location','southwest');

yyaxis right;
hold on;
errorbar(gerbil_SAM_data.AM_Hz_2use, gerbil_SAM_data.mu_19wk, gerbil_SAM_data.sem_19wk, 'LineWidth', lw2, 'Color', get_color('b'), 'LineStyle','--')
% errorbar(gerbil_SAM_data.AM_Hz_2use, gerbil_SAM_data.mu_44wk, gerbil_SAM_data.sem_44wk, 'LineWidth', lw2, 'Color', get_color('b'))
% errorbar(gerbil_SAM_data.AM_Hz_2use, gerbil_SAM_data.mu_75wk, gerbil_SAM_data.sem_75wk, 'LineWidth', lw2, 'Color', get_color('lb'))
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
xlabel('AM frequency (Hz)');
xlim(xlim_val)

% axes(sp_ax(3));
% hold on;
% errorbar(gerbil_SAM_data.mu_19wk, gerbil_dAM_data.mu_19wk, gerbil_dAM_data.sem_19wk/2, gerbil_dAM_data.sem_19wk/2, gerbil_SAM_data.sem_19wk/2, gerbil_SAM_data.sem_19wk/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', 'b')
% % errorbar(gerbil_SAM_data.mu_44wk, gerbil_dAM_data.mu_44wk, gerbil_dAM_data.sem_44wk/2, gerbil_dAM_data.sem_44wk/2, gerbil_SAM_data.sem_44wk/2, gerbil_SAM_data.sem_44wk/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('b'))
% % errorbar(gerbil_SAM_data.mu_75wk, gerbil_dAM_data.mu_75wk, gerbil_dAM_data.sem_75wk/2, gerbil_dAM_data.sem_75wk/2, gerbil_SAM_data.sem_75wk/2, gerbil_SAM_data.sem_75wk/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('lb'))
% xlabel('SAM power, dB');
% % ylabel('dAM power, dB');
% 
% lin_mdl= fitlm([gerbil_SAM_data.mu_19wk(:); gerbil_SAM_data.mu_44wk(:); gerbil_SAM_data.mu_75wk(:)], [gerbil_dAM_data.mu_19wk(:); gerbil_dAM_data.mu_44wk(:); gerbil_dAM_data.mu_75wk(:)]);
fprintf("Gerbil:");
lin_mdl= fitlm(gerbil_SAM_data.mu_19wk(:), gerbil_dAM_data.mu_19wk(:))
% text(.02, .9, sprintf('R^2=%.2f', lin_mdl.Rsquared.Adjusted), 'Units','normalized');
% xest= xlim();
% yest= predict(lin_mdl, xest(:));
% plot(xest, yest, 'k-', 'LineWidth', 1.5)

end

function [mice_dAM_data, mice_SAM_data, lin_mdl] = plot_mice_sAM_dAM_data(sp_ax, Root_dAM_power_dir)
freq_tick= [16 50 200 800];
lw1= 1;
lw2= 1.25;

% dAM data 
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\dAM_power_data\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\';
gerbil_fName= [Root_dAM_power_dir 'mice_two_group_power.mat'];
mice_dAM_data= load(gerbil_fName);
AMfreqs_Hz= mice_dAM_data.AM_Hz_2use;
if ~mice_dAM_data.use_dB
    error('Need dB, but not dB');
end
mice_dAM_data.mu_control= nanmean(mice_dAM_data.pow_val_control_frac);
mice_dAM_data.sem_control= nanstd(mice_dAM_data.pow_val_control_frac)/sqrt(size(mice_dAM_data.pow_val_control_frac,2));
mice_dAM_data.mu_exposed= nanmean(mice_dAM_data.pow_val_exposed_frac);
mice_dAM_data.sem_exposed= nanstd(mice_dAM_data.pow_val_exposed_frac)/sqrt(size(mice_dAM_data.pow_val_exposed_frac,2));

% Discrete data 
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Mice\';

discrete_freq_data_sem_control= readtable([Discrete_FFT_Dir 'MousediscreteEFRs_Control.xlsx']);
discrete_freqs_control= discrete_freq_data_sem_control.Rate;
discrete_freq_data_sem_control= table2array(discrete_freq_data_sem_control);
discrete_freq_data_sem_control= discrete_freq_data_sem_control(ismember(discrete_freqs_control, AMfreqs_Hz), 2:end);
discrete_freq_data_sem_control= discrete_freq_data_sem_control'; % make columns ~ freq
discrete_freq_data_sem_control= mag2db(discrete_freq_data_sem_control);

discrete_freq_data_exposed= readtable([Discrete_FFT_Dir 'MousediscreteEFRs_Exposed.xlsx']);
discrete_freqs_exposed= discrete_freq_data_exposed.Rate;
discrete_freq_data_exposed= table2array(discrete_freq_data_exposed);
discrete_freq_data_exposed= discrete_freq_data_exposed(ismember(discrete_freqs_exposed, AMfreqs_Hz), 2:end);
discrete_freq_data_exposed= discrete_freq_data_exposed'; % make columns ~ freq
discrete_freq_data_exposed= mag2db(discrete_freq_data_exposed);

if ~isequal(discrete_freqs_control, discrete_freqs_exposed)
    error('AM freqs differ');
end

mice_SAM_data= struct(...
    'discrete_freq_data_control', discrete_freq_data_sem_control, 'discrete_freq_data_exposed', discrete_freq_data_exposed, 'AM_Hz_2use', AMfreqs_Hz, ...
    'mu_control', nanmean(discrete_freq_data_sem_control), 'mu_exposed', nanmean(discrete_freq_data_exposed), ...
    'sem_control', nanstd(discrete_freq_data_sem_control)/sqrt(size(discrete_freq_data_sem_control,2)), 'sem_exposed', nanstd(discrete_freq_data_exposed)/sqrt(size(discrete_freq_data_exposed,2)));


axes(sp_ax(1));
hold on;
lHan(1)= errorbar(mice_dAM_data.AM_Hz_2use, mice_dAM_data.mu_control, mice_dAM_data.sem_control, 'LineWidth', lw2, 'Color', get_color('dg'));
% lHan(2)= errorbar(mice_dAM_data.AM_Hz_2use, mice_dAM_data.mu_exposed, mice_dAM_data.sem_exposed, 'LineWidth', lw2, 'Color', get_color('lg'));
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
text(.05, .1, 'Mice', 'Color', get_color('dg'), 'Units', 'normalized')
% legend(lHan, {'Control', 'Noise-exposed'}, 'Location','southwest')

% axes(sp_ax(2));
yyaxis right;
hold on;
errorbar(mice_SAM_data.AM_Hz_2use, mice_SAM_data.mu_control, mice_SAM_data.sem_control, 'LineWidth', lw2, 'Color', get_color('dg'), 'LineStyle', '--')
% errorbar(mice_SAM_data.AM_Hz_2use, mice_SAM_data.mu_exposed, mice_SAM_data.sem_exposed, 'LineWidth', lw2, 'Color', get_color('lg'))
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
xlabel('AM frequency (Hz)');

% axes(sp_ax(3));
% hold on;
% errorbar(mice_SAM_data.mu_control, mice_dAM_data.mu_control, mice_dAM_data.sem_control/2, mice_dAM_data.sem_control/2, mice_SAM_data.sem_control/2, mice_SAM_data.sem_control/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('dg'))
% % errorbar(mice_SAM_data.mu_exposed, mice_dAM_data.mu_exposed, mice_dAM_data.sem_exposed/2, mice_dAM_data.sem_exposed/2, mice_SAM_data.sem_exposed/2, mice_SAM_data.sem_exposed/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('lg'))
% xlabel('SAM power, dB');
% % ylabel('dAM power, dB');
% 
% lin_mdl= fitlm([mice_SAM_data.mu_control(:); mice_SAM_data.mu_exposed(:)], [mice_dAM_data.mu_control(:); mice_dAM_data.mu_exposed(:)]);
% lin_mdl= fitlm(mice_SAM_data.mu_control(:), mice_dAM_data.mu_control(:));
% text(.02, .9, sprintf('R^2=%.2f', lin_mdl.Rsquared.Adjusted), 'Units','normalized');
% xest= xlim();
% yest= predict(lin_mdl, xest(:));
% plot(xest, yest, 'k-', 'LineWidth', 1.5)

fprintf("Mice: ")
lin_mdl= fitlm(mice_SAM_data.mu_control(:), mice_dAM_data.mu_control(:))


end

function [rat_dAM_data, rat_SAM_data, lin_mdl] = plot_rat_sAM_dAM_data(sp_ax, configuration, Root_dAM_power_dir)
freq_tick= [16 50 200 800];
lw1= 1;
lw2= 1.25;

% dAM data 
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\dAM_power_data\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\';
rat_fName= [Root_dAM_power_dir 'rat_three_group_power_noise_' configuration '.mat'];
rat_dAM_data= load(rat_fName);
AMfreqs_Hz= rat_dAM_data.AM_Hz_2use;
if ~rat_dAM_data.use_dB
    error('Need dB, but not dB');
end
rat_dAM_data.mu_control= nanmean(rat_dAM_data.pow_val_ctrl_frac);
rat_dAM_data.sem_control= nanstd(rat_dAM_data.pow_val_ctrl_frac)/sqrt(size(rat_dAM_data.pow_val_ctrl_frac,2));

% Discrete data 
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Rat_galactose\';

if strcmp(configuration, 'vertical')
    discrete_freq_data_sem_control= readtable([Discrete_FFT_Dir 'NAMfmod_historical_ctrl_ch1_vert.xlsx']);
    discrete_freqs_control= discrete_freq_data_sem_control.fmod_channel_1;
    ttl_str= 'Rat (vertical)';
    deep_color= get_color('r');
elseif strcmp(configuration, 'horizontal')
    discrete_freq_data_sem_control= readtable([Discrete_FFT_Dir 'NAMfmod_historical_ctrl_ch2_horz.xlsx']);
    discrete_freqs_control= discrete_freq_data_sem_control.fmodChannel2;
    ttl_str= 'Rat (horizontal)';
    deep_color= get_color('prp');
end

discrete_freq_data_sem_control= table2array(discrete_freq_data_sem_control);
discrete_freq_data_sem_control= discrete_freq_data_sem_control(ismember(discrete_freqs_control, AMfreqs_Hz), 2:end);
discrete_freq_data_sem_control= discrete_freq_data_sem_control'; % make columns ~ freq
discrete_freq_data_sem_control= mag2db(discrete_freq_data_sem_control);

rat_SAM_data= struct(...
    'discrete_freq_data_control', discrete_freq_data_sem_control, 'AM_Hz_2use', AMfreqs_Hz, ...
    'mu_control', nanmean(discrete_freq_data_sem_control), 'sem_control', nanstd(discrete_freq_data_sem_control)/sqrt(size(discrete_freq_data_sem_control,2)));

xlim_val= [.9*min(rat_dAM_data.AM_Hz_2use), 1.1*max(rat_dAM_data.AM_Hz_2use)];

axes(sp_ax(1));
hold on;
lHan= errorbar(rat_dAM_data.AM_Hz_2use, rat_dAM_data.mu_control, rat_dAM_data.sem_control, 'LineWidth', lw2, 'Color', deep_color);
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
text(.05, .1, ttl_str, 'Color', deep_color, 'Units', 'normalized')
xlim(xlim_val)
% legend(lHan, 'Control', 'Location','southwest');


yyaxis right;
% axes(sp_ax(2));
hold on;
errorbar(rat_SAM_data.AM_Hz_2use, rat_SAM_data.mu_control, rat_SAM_data.sem_control, 'LineWidth', lw2, 'Color', deep_color, 'LineStyle', '--')
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
xlabel('AM frequency (Hz)');
xlim(xlim_val)

% axes(sp_ax(3));
% hold on;
% errorbar(rat_SAM_data.mu_control, rat_dAM_data.mu_control, rat_dAM_data.sem_control/2, rat_dAM_data.sem_control/2, rat_SAM_data.sem_control/2, rat_SAM_data.sem_control/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', deep_color)
% xlabel('SAM power, dB');
% % ylabel('dAM power, dB');
% lin_mdl= fitlm(rat_SAM_data.mu_control(:), rat_dAM_data.mu_control(:));
% text(.02, .9, sprintf('R^2=%.2f', lin_mdl.Rsquared.Adjusted), 'Units','normalized');
% xest= xlim();
% yest= predict(lin_mdl, xest(:));
% plot(xest, yest, 'k-', 'LineWidth', 1.5)

fprintf("Rats")
lin_mdl= fitlm(rat_SAM_data.mu_control(:), rat_dAM_data.mu_control(:))

end

%%
function [human_dAM_data, human_SAM_data, lin_mdl] = plot_human_sAM_dAM_data(sp_ax, Root_dAM_power_dir)
freq_tick= [16 50 200 800];
lw1= 1;
lw2= 1.25;

% dAM data 
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\dAM_power_data\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\';
human_fName= [Root_dAM_power_dir 'human_two_group_power.mat'];
human_dAM_data= load(human_fName);
AMfreqs_Hz= human_dAM_data.AM_Hz_2use;
if ~human_dAM_data.use_dB
    error('Need dB, but not dB');
end
human_dAM_data.mu_young= nanmean(human_dAM_data.pow_val_young_frac);
human_dAM_data.sem_young= nanstd(human_dAM_data.pow_val_young_frac)/sqrt(size(human_dAM_data.pow_val_young_frac,2));
human_dAM_data.mu_ma= nanmean(human_dAM_data.pow_val_midaged_frac);
human_dAM_data.sem_ma= nanstd(human_dAM_data.pow_val_midaged_frac)/sqrt(size(human_dAM_data.pow_val_midaged_frac,2));

% Discrete data 
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Humans\EFR_Discrete\';

discrete_freq_data_sem_young= readtable([Discrete_FFT_Dir 'EFR_Human_Young.xlsx']);
discrete_freqs_young= discrete_freq_data_sem_young.Level;
discrete_freq_data_sem_young= table2array(discrete_freq_data_sem_young);
discrete_freq_data_sem_young= discrete_freq_data_sem_young(ismember(discrete_freqs_young, AMfreqs_Hz), 2:end);
discrete_freq_data_sem_young= discrete_freq_data_sem_young'; % make columns ~ freq
discrete_freq_data_sem_young= mag2db(discrete_freq_data_sem_young);

discrete_freq_data_ma= readtable([Discrete_FFT_Dir 'EFR_Human_MA.xlsx']);
discrete_freqs_ma= discrete_freq_data_ma.Level;
discrete_freq_data_ma= table2array(discrete_freq_data_ma);
discrete_freq_data_ma= discrete_freq_data_ma(ismember(discrete_freqs_ma, AMfreqs_Hz), 2:end);
discrete_freq_data_ma= discrete_freq_data_ma'; % make columns ~ freq
discrete_freq_data_ma= mag2db(discrete_freq_data_ma);


if ~isequal(discrete_freqs_young, discrete_freqs_ma)
    error('AM freqs differ');
end

human_SAM_data= struct(...
    'discrete_freq_data_young', discrete_freq_data_sem_young, 'discrete_freq_data_ma', discrete_freq_data_ma, 'AM_Hz_2use', AMfreqs_Hz, ...
    'mu_young', nanmean(discrete_freq_data_sem_young), 'mu_ma', nanmean(discrete_freq_data_ma), ...
    'sem_young', nanstd(discrete_freq_data_sem_young)/sqrt(size(discrete_freq_data_sem_young,2)), 'sem_ma', nanstd(discrete_freq_data_ma)/sqrt(size(discrete_freq_data_ma,2)));

xlim_val= [.9*min(human_dAM_data.AM_Hz_2use), 1.1*max(human_dAM_data.AM_Hz_2use)];


axes(sp_ax(1));
yyaxis left;

lHan(1)= errorbar(human_dAM_data.AM_Hz_2use, human_dAM_data.mu_young, human_dAM_data.sem_young, 'LineWidth', lw2, 'Color', get_color('m'));
% lHan(2)= errorbar(human_dAM_data.AM_Hz_2use, human_dAM_data.mu_ma, human_dAM_data.sem_ma, 'LineWidth', lw2, 'Color', get_color('lavender'));
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
text(.05, .1, 'Human', 'Color', get_color('m'), 'Units', 'normalized');
xlim(xlim_val);
% ylabel('dAM power, dB');
% legend(lHan, {'Young', 'Middle-aged'}, 'Location','southwest');

yyaxis right;
hold on;
errorbar(human_SAM_data.AM_Hz_2use, human_SAM_data.mu_young, human_SAM_data.sem_young, 'LineStyle', '--', 'LineWidth', lw2, 'Color', get_color('m'))
% errorbar(human_SAM_data.AM_Hz_2use, human_SAM_data.mu_ma, human_SAM_data.sem_ma, 'LineWidth', lw2, 'Color', get_color('lavender'))
set(gca, 'XScale', 'log', 'XTick', freq_tick, 'YColor', get_color('k'))
xlabel('AM frequency (Hz)');
xlim(xlim_val);
% ylabel('SAM power, dB');

% axes(sp_ax(end));
% hold on;
% errorbar(human_SAM_data.mu_young, human_dAM_data.mu_young, human_dAM_data.sem_young/2, human_dAM_data.sem_young/2, human_SAM_data.sem_young/2, human_SAM_data.sem_young/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('m'))
% % errorbar(human_SAM_data.mu_ma, human_dAM_data.mu_ma, human_dAM_data.sem_ma/2, human_dAM_data.sem_ma/2, human_SAM_data.sem_ma/2, human_SAM_data.sem_ma/2, 'LineWidth', lw1, 'LineStyle', 'none', 'Color', get_color('lavender'))
% xlabel('SAM power, dB');
% ylabel('dAM power, dB');
% 
% % lin_mdl= fitlm([human_SAM_data.mu_young(:); human_SAM_data.mu_ma(:)], [human_dAM_data.mu_young(:); human_dAM_data.mu_ma(:)]);
% lin_mdl= fitlm(human_SAM_data.mu_young(:), human_dAM_data.mu_young(:));
% text(.02, .9, sprintf('R^2=%.2f', lin_mdl.Rsquared.Adjusted), 'Units','normalized');
% xest= xlim();
% yest= predict(lin_mdl, xest(:));
% plot(xest, yest, 'k-', 'LineWidth', 1.5)

fprintf("Humans: ")
lin_mdl= fitlm(human_SAM_data.mu_young(:), human_dAM_data.mu_young(:))

end


%%
function sp_ax =get_axes(fig_TF, nSProws, nSPcols)

figure(fig_TF);
figSize_cm= [3, 3, 17, 9];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
set(gcf, figure_prop_name, figure_prop_val);
clf;

xc= .065;
yc= .01;
xw= .245;
xs= .07;
yw= .36;
ys= .12;

sp_ax= nan(nSProws*nSPcols, 1);
for spNum=1:(nSProws*nSPcols)
    n_yskips= nSProws - ceil(spNum/nSPcols);
    yc_current= yc + n_yskips*(yw+ys);
    if spNum<=2*nSPcols
        yc_current= yc_current + .75*ys;
    end

    n_xskips= rem(spNum-1, nSPcols);
    xc_current= xc + n_xskips*(xw+xs);

    sp_ax(spNum) = subplot(nSProws, nSPcols, spNum);
    set(gca, 'Units', 'normalized', 'Position', [xc_current, yc_current, xw, yw])

end

pos_sp_ax_end= get(sp_ax(end), 'Position');
pos_sp_ax_end(1)= pos_sp_ax_end(1)+.03;
pos_sp_ax_end(3)= pos_sp_ax_end(3)+.02;
set(sp_ax(end), 'Position', pos_sp_ax_end);
end
