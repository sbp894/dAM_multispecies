clear;
clc;

num_cols= 4;
figure(1);
clf;
set(gcf, 'Units', 'centimeters', 'Position', [5, 3, 17, 5])
h_map= tiledlayout(1, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
sp_ax= zeros(1, num_cols);
for ax_var=1:num_cols
    sp_ax(ax_var)= nexttile;
end

col= struct('human',helper.get_color('m'), 'gerbil', helper.get_color('b'), 'mice', helper.get_color('g'));

% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\0ms_12Hz\';
Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\5ms_12Hz\';
% Root_dAM_power_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\dAM_power_data\10ms_12Hz\';


hold on;
[human_dAM_data, human_SAM_data, lin_mdl_human] = plot_human_sAM_dAM_data(sp_ax(1), Root_dAM_power_dir);
[gerbil_dAM_data, gerbil_SAM_data, lin_mdl_gerbil] = plot_gerbil_sAM_dAM_data(sp_ax(2), Root_dAM_power_dir);
[mice_dAM_data, mice_SAM_data, lin_mdl_mice] = plot_mice_sAM_dAM_data(sp_ax(3), Root_dAM_power_dir);

fprintf("humans=%d, gerbils=%d,mice=%d\n", length(human_dAM_data.dam_human_names), length(gerbil_dAM_data.dam_gerbil_names), length(mice_dAM_data.dam_mice_names))

all_dam_data= [normalize(human_dAM_data.dam_freq_data_ctrl(:)); normalize(gerbil_dAM_data.dam_freq_data_19wk(:)); normalize(mice_dAM_data.dam_freq_data_ctrl(:))];
all_sam_data= [normalize(human_SAM_data.sam_freq_data_ctrl(:)); normalize(gerbil_SAM_data.sam_freq_data_19wk(:)); normalize(mice_SAM_data.sam_freq_data_ctrl(:))];

axes(sp_ax(1));
plot(human_SAM_data.sam_freq_data_ctrl(:), human_dAM_data.dam_freq_data_ctrl(:), '.', 'DisplayName', 'human', 'Color', col.human)
assert(lin_mdl_human.Coefficients.pValue(2)<1e-4)
text(.95, .17, sprintf('$p<10^{-4}$'), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
text(.95, .05, sprintf('$R^2=%.2f$', lin_mdl_human.Rsquared.Adjusted), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
title('Human', 'Color', col.human)
ylabel('Frac. dAM power (dB)');
text(-.2, 1.14, 'A', 'FontSize', 10, 'FontWeight','bold', 'Units','normalized')

axes(sp_ax(2));
plot(gerbil_SAM_data.sam_freq_data_19wk(:), gerbil_dAM_data.dam_freq_data_19wk(:), '.', 'DisplayName', 'gerbil', 'Color', col.gerbil)
assert(lin_mdl_gerbil.Coefficients.pValue(2)<1e-4)
text(.95, .17, sprintf('$p<10^{-4}$'), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
text(.95, .05, sprintf('$R^2=%.2f$', lin_mdl_gerbil.Rsquared.Adjusted), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
title('Gerbil', 'Color', col.gerbil)
xlabel('SAM power (dB)');

axes(sp_ax(3));
plot(mice_SAM_data.sam_freq_data_ctrl(:), mice_dAM_data.dam_freq_data_ctrl(:), '.', 'DisplayName', 'mice', 'Color', col.mice)
assert(lin_mdl_mice.Coefficients.pValue(2)<1e-4)
text(.95, .17, sprintf('$p<10^{-4}$'), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
text(.95, .05, sprintf('$R^2=%.2f$', lin_mdl_mice.Rsquared.Adjusted), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
title('Mice', 'Color', col.mice)

axes(sp_ax(4))
hold on;

plot(normalize(human_SAM_data.sam_freq_data_ctrl(:)), normalize(human_dAM_data.dam_freq_data_ctrl(:)), '.', 'DisplayName', 'human', 'Color', col.human)
plot(normalize(gerbil_SAM_data.sam_freq_data_19wk(:)), normalize(gerbil_dAM_data.dam_freq_data_19wk(:)), '.', 'DisplayName', 'gerbil', 'Color', col.gerbil)
plot(normalize(mice_dAM_data.dam_freq_data_ctrl(:)), normalize(mice_SAM_data.sam_freq_data_ctrl(:)), '.', 'DisplayName', 'mice', 'Color', col.mice)
full_mdl= fitlm(all_dam_data, all_sam_data);
assert(full_mdl.Coefficients.pValue(2)<1e-4)
lin_fit_x= [min(all_sam_data); max(all_sam_data)];
lin_fit_y= predict(full_mdl, lin_fit_x);
% plot(lin_fit_x, lin_fit_y, 'k-', 'LineWidth', 2)

text(.95, .17, sprintf('$p<10^{-4}$'), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
text(.95, .05, sprintf('$R^2=%.2f$', full_mdl.Rsquared.Adjusted), 'Units', 'normalized', 'Interpreter','latex', 'HorizontalAlignment', 'right')
text(-.1, 1.14, 'B', 'FontSize', 10, 'FontWeight','bold', 'Units','normalized')
title('Summary')
xlabel('Transformed SAM power (au)');
ylabel('Transformed dAM power (au)');

set(findall(gcf,'-property','box'),'box', 'off');


%%
if strcmp(getenv('COMPUTERNAME'), 'NB-VATS-LAB1')
    fig_dir= 'C:\Users\sap245\Google Drive\PostDoc\dAM_MultiSpecies\Figures\';
else 
    fig_dir= 'G:\My Drive\PostDoc\dAM_MultiSpecies\Figures\';
end

do_save_fig= 0;
fig_name= [fig_dir 'Fig6_sam_vs_dam_indv'];
if do_save_fig
    print(fig_name, '-dpng', '-r600')
end

%%
function [gerbil_dAM_data, gerbil_SAM_data, lin_mdl] = plot_gerbil_sAM_dAM_data(sp_ax, Root_dAM_power_dir)
gerbil_fName= [Root_dAM_power_dir 'gerbil_three_group_power.mat'];
all_gerbil_dAM_data= load(gerbil_fName);
AMfreqs_Hz_dAM= all_gerbil_dAM_data.AM_Hz_2use;
if ~all_gerbil_dAM_data.use_dB
    error('Need dB, but not dB');
end
dam_gerbil_names= all_gerbil_dAM_data.all_names(all_gerbil_dAM_data.cur_young0_ma1_old2_flag==0);
dam_gerbil_names= cell2mat(cellfun(@(x) str2double(x(4:min(length(x), 5))), dam_gerbil_names, 'UniformOutput', false));

% Discrete data
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\ARO2023_files\Discrete_FFTdata\';

sam_freq_data_19wk= readtable([Discrete_FFT_Dir 'Gerbilefr3k_19wk.xlsx']);
AMfreqs_Hz_SAM= sam_freq_data_19wk.Level;
sam_gerbil_names= sam_freq_data_19wk.Properties.VariableNames(2:end);
sam_gerbil_names= cell2mat(cellfun(@(x) str2double(x(4:min(length(x), 5))), sam_gerbil_names, 'UniformOutput', false));

sam_freq_data_19wk= table2array(sam_freq_data_19wk);
sam_freq_data_19wk= sam_freq_data_19wk(ismember(AMfreqs_Hz_SAM, AMfreqs_Hz_dAM), 2:end);
sam_freq_data_19wk= sam_freq_data_19wk'; % make columns ~ freq
sam_freq_data_19wk= mag2db(sam_freq_data_19wk);


common_animals= intersect(dam_gerbil_names, sam_gerbil_names);
valid_dAM_animal_inds= cellfun(@(x) find(x==dam_gerbil_names), num2cell(common_animals));
valid_sam_animal_inds= cellfun(@(x) find(x==sam_gerbil_names), num2cell(common_animals));

common_freqs= intersect(AMfreqs_Hz_dAM, AMfreqs_Hz_SAM);
valid_dAM_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_dAM), num2cell(common_freqs));
valid_sam_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_SAM), num2cell(common_freqs));



gerbil_SAM_data= struct('sam_freq_data_19wk', sam_freq_data_19wk(valid_sam_animal_inds, valid_sam_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'sam_gerbil_names', common_animals);
gerbil_dAM_data= struct('dam_freq_data_19wk', all_gerbil_dAM_data.pow_val_young_frac(valid_dAM_animal_inds,valid_dAM_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'dam_gerbil_names', common_animals);

fprintf("Gerbils --------------- \n")
lin_mdl= fitlm(gerbil_SAM_data.sam_freq_data_19wk(:), gerbil_dAM_data.dam_freq_data_19wk(:))

end

%%
function [human_dAM_data, human_SAM_data, lin_mdl] = plot_human_sAM_dAM_data(sp_ax, Root_dAM_power_dir)
human_fName= [Root_dAM_power_dir 'human_two_group_power.mat'];
all_human_dAM_data= load(human_fName);
AMfreqs_Hz_dAM= all_human_dAM_data.AM_Hz_2use;
if ~all_human_dAM_data.use_dB
    error('Need dB, but not dB');
end
dam_human_names= all_human_dAM_data.allfiles(all_human_dAM_data.cur_young0_ma1_flag==0);
dam_human_names= {dam_human_names.name}';
dam_human_names= cell2mat(cellfun(@(x) str2double(x(3:min(length(x), 4))), dam_human_names, 'UniformOutput', false));

% Discrete data
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Humans\EFR_Discrete\';

sam_freq_data_young= readtable([Discrete_FFT_Dir 'EFR_Human_Young.xlsx']);
AMfreqs_Hz_SAM= sam_freq_data_young.Level;
sam_human_names= sam_freq_data_young.Properties.VariableNames(2:end);
sam_human_names= cell2mat(cellfun(@(x) str2double(x(3:min(length(x), 5))), sam_human_names, 'UniformOutput', false));

sam_freq_data_young= table2array(sam_freq_data_young);
sam_freq_data_young= sam_freq_data_young(ismember(AMfreqs_Hz_SAM, AMfreqs_Hz_dAM), 2:end);
sam_freq_data_young= sam_freq_data_young'; % make columns ~ freq
sam_freq_data_young= mag2db(sam_freq_data_young);


common_humans= intersect(dam_human_names, sam_human_names);
valid_dAM_human_inds= cellfun(@(x) find(x==dam_human_names), num2cell(common_humans));
valid_sam_human_inds= cellfun(@(x) find(x==sam_human_names), num2cell(common_humans));

common_freqs= intersect(AMfreqs_Hz_dAM, AMfreqs_Hz_SAM);
valid_dAM_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_dAM), num2cell(common_freqs));
valid_sam_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_SAM), num2cell(common_freqs));



human_SAM_data= struct('sam_freq_data_ctrl', sam_freq_data_young(valid_sam_human_inds, valid_sam_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'sam_human_names', common_humans);
human_dAM_data= struct('dam_freq_data_ctrl', all_human_dAM_data.pow_val_young_frac(valid_dAM_human_inds,valid_dAM_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'dam_human_names', common_humans);

fprintf("humans --------------- \n")
lin_mdl= fitlm(human_SAM_data.sam_freq_data_ctrl(:), human_dAM_data.dam_freq_data_ctrl(:))
end

%% mice
function [mice_dAM_data, mice_SAM_data, lin_mdl] = plot_mice_sAM_dAM_data(sp_ax, Root_dAM_power_dir)

mice_fName= [Root_dAM_power_dir 'mice_two_group_power.mat'];
all_mice_dAM_data= load(mice_fName);
AMfreqs_Hz_dAM= all_mice_dAM_data.AM_Hz_2use;
if ~all_mice_dAM_data.use_dB
    error('Need dB, but not dB');
end
dam_mice_names= all_mice_dAM_data.allfiles(all_mice_dAM_data.cur_ctrl0_exposed1_flag==0);
dam_mice_names= {dam_mice_names.name}';
dam_mice_names= cell2mat(cellfun(@(x) str2double(x(2:4)), dam_mice_names, 'UniformOutput', false));

% Discrete data
Discrete_FFT_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Mice\';
sam_freq_data_ctrl= readtable([Discrete_FFT_Dir 'MousediscreteEFRs_Control.xlsx']);
AMfreqs_Hz_SAM= sam_freq_data_ctrl.Rate;
sam_mice_names= sam_freq_data_ctrl.Properties.VariableNames(2:end);
sam_mice_names= cell2mat(cellfun(@(x) str2double(x(2:end)), sam_mice_names, 'UniformOutput', false));

sam_freq_data_ctrl= table2array(sam_freq_data_ctrl);
sam_freq_data_ctrl= sam_freq_data_ctrl(ismember(AMfreqs_Hz_SAM, AMfreqs_Hz_dAM), 2:end);
sam_freq_data_ctrl= sam_freq_data_ctrl'; % make columns ~ freq
sam_freq_data_ctrl= mag2db(sam_freq_data_ctrl);

bad_mice= 109;
common_animals= intersect(dam_mice_names, sam_mice_names);
common_animals= setdiff(common_animals, bad_mice);
valid_dAM_animal_inds= cellfun(@(x) find(x==dam_mice_names), num2cell(common_animals));
valid_sam_animal_inds= cellfun(@(x) find(x==sam_mice_names), num2cell(common_animals));

common_freqs= intersect(AMfreqs_Hz_dAM, AMfreqs_Hz_SAM);
valid_dAM_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_dAM), num2cell(common_freqs));
valid_sam_freq_inds= cellfun(@(x) find(x==AMfreqs_Hz_SAM), num2cell(common_freqs));



mice_SAM_data= struct('sam_freq_data_ctrl', sam_freq_data_ctrl(valid_sam_animal_inds, valid_sam_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'sam_mice_names', common_animals);
mice_dAM_data= struct('dam_freq_data_ctrl', all_mice_dAM_data.pow_val_control_frac(valid_dAM_animal_inds,valid_dAM_freq_inds), 'AM_Hz_2use', common_freqs, ...
    'dam_mice_names', common_animals);

fprintf("mice --------------- \n")
lin_mdl= fitlm(mice_SAM_data.sam_freq_data_ctrl(:), mice_dAM_data.dam_freq_data_ctrl(:))
end

