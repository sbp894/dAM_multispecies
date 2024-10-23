clear;
clc;

doSaveFig= 0;
doSaveGrandAvg= 0;

CodeDir= '..\';
addpath(CodeDir);


%% load estimate AM trajectory

fig_TF= 4;
nSProws= 5;
nSPcols= 5;
sp_ax_TF = get_axes(fig_TF, nSProws, nSPcols);

%% human plots: first column
verboseAxis= 1;
[dam_traj_Hz_human_NFvalid, all_pow_am_human, all_pow_am_human_NF, mu_pow_am_human, mu_pow_am_human_NF]= load_species_data('human', sp_ax_TF( 1 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg);

%% gerbil plots: second column
verboseAxis= 0;
[dam_traj_Hz_gerbil_NFvalid, all_pow_am_gerbil, all_pow_am_gerbil_NF, mu_pow_am_gerbil, mu_pow_am_gerbil_NF]= load_species_data('gerbil', sp_ax_TF( 2 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg);

%% rat sagittal plots: third column
[dam_traj_Hz_rat_NFvalid_vert, all_pow_am_rat_vert, all_pow_am_rat_NF_vert, mu_pow_am_rat_vert, mu_pow_am_rat_NF_vert]= load_species_data('rat_vert', sp_ax_TF( 3 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg); % rat_vert= rat + vertical (sagittal) 

%% rat plots: fourth column
[dam_traj_Hz_rat_NFvalid_horz, all_pow_am_rat_horz, all_pow_am_rat_NF_horz, mu_pow_am_rat_horz, mu_pow_am_rat_NF_horz]= load_species_data('rat_horz', sp_ax_TF( 4 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg);  % rat_vert= rat + vertical (interaural)

%% mice plots: fourth column
[dam_traj_Hz_mice_NFvalid, all_pow_am_mice, all_pow_am_mice_NF, mu_pow_am_mice, mu_pow_am_mice_NF]= load_species_data('mice', sp_ax_TF( 5 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg);

axes(sp_ax_TF(nSPcols*(nSProws-1)))
text(.02, .3, 'dAM', 'Color', helper.get_color('dg'), 'Units', 'normalized')
text(.02, .13, 'noise floor', 'Color', helper.get_color('k'), 'Units', 'normalized')

ylim(sp_ax_TF(18), [-39, 0]+ max(get(sp_ax_TF(18), 'YLim')))

%% panel letters 
% delete(panel_han);

axes(sp_ax_TF(1)); 
panel_han(1)= text(-.25, 1.075, "A", 'Units','normalized', 'FontWeight','bold');

axes(sp_ax_TF(6)); 
panel_han(2)= text(-.25, 1.075, "B", 'Units','normalized', 'FontWeight','bold');

axes(sp_ax_TF(11)); 
panel_han(3)= text(-.25, 1.075, "C", 'Units','normalized', 'FontWeight','bold');

axes(sp_ax_TF(16)); 
panel_han(4)= text(-.25, 1.075, "D", 'Units','normalized', 'FontWeight','bold');

axes(sp_ax_TF(21)); 
panel_han(5)= text(-.25, 1.075, "E", 'Units','normalized', 'FontWeight','bold');

%% fig properties 
set(findall(gcf,'-property','FontSize'),'FontSize', 7);
set(findall(gcf,'-property','TickDir'),'TickDir', 'both');
set(findall(gcf,'-property','box'),'box', 'off');

do_save_fig= 0;
fig_dir= ['figures' filesep];
fig_name= [fig_dir 'Fig4_robust_dAM_all_species'];
if do_save_fig
    print(fig_name, '-dpng', '-r600')
end
%%
rmpath(CodeDir);

%% all functions
function [dam_traj_Hz_NFvalid, all_pow_am, all_pow_am_NF, mu_pow_am, mu_pow_am_NF, lHan]= load_species_data(Species, sp_ax, verboseAxis, doSaveGrandAvg)

fs_env= 5e3;
bp_window_Hz= [5, 1.5e3];
Filter_HalfWindow_Hz= 12;

%% load AM trajectory
switch lower(Species)
    case 'human'

        % load dAM trajectory
        dam_traj_Hz_struct= load('AMfreqvec_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        Xorg= (1:length(dam_traj_Hz_struct.AMfreqvec_est))'/dam_traj_Hz_struct.fs;
        Yorg= dam_traj_Hz_struct.AMfreqvec_est(:);

        % load sitmuli
        fs_sig= 10e3;
        [stim_wav, fs_stim_org]= audioread('files\wavfiles\DoDAMPatternsNEW_3k_gerbils_ER3c2cc.wav');
        stim_wav= helper.gen_resample(stim_wav, fs_stim_org, fs_sig);

%         Species_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Humans\DoDAM3k_Human\YA\';
        Species_Dir= 'files\Exports_Humans\DoDAM3k_Human\YA\';
        all_files= dir([Species_Dir '*Fz-Rtip']);
        fs_ffr= 1/(61.035156e-6);
        freq_lim_kHz= [1, 5];
        mu_pow_lim_dB= 30;

        tMin_s= 0;
        tMax_s= 1;

        col_deep= helper.get_color('m');
        col_light= helper.get_color('lavender');
        FFR_gain= 1;

    case 'gerbil'
        % load dAM trajectory
        dam_traj_Hz_struct= load('AMfreqvec_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        Xorg= (1:length(dam_traj_Hz_struct.AMfreqvec_est))'/dam_traj_Hz_struct.fs;
        Yorg= dam_traj_Hz_struct.AMfreqvec_est(:);

        % load sitmuli
        fs_sig= 10e3;
%         [stim_wav, fs_stim_org]= audioread('D:\Dropbox\Pitts_files\Aravind_ExpAM\wavfiles\DoDAMPatternsNEW_3k_gerbils_ER3c2cc.wav');
        [stim_wav, fs_stim_org]= audioread('files\wavfiles\DoDAMPatternsNEW_3k_gerbils_ER3c2cc.wav');
        stim_wav= helper.gen_resample(stim_wav, fs_stim_org, fs_sig);

%         Species_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Gerbils\DoDAM3k\19wk\';
        Species_Dir= 'files\Exports_Gerbils\DoDAM3k\19wk\';
        all_files= dir([Species_Dir 'GER*']);
        fs_ffr= 1/(40.96e-6);
        freq_lim_kHz= [1, 5];
        mu_pow_lim_dB= 30;

        tMin_s= 0;
        tMax_s= 1;

        col_deep= helper.get_color('b');
        col_light= helper.get_color('lb');
        FFR_gain= 1;

    case {'rat', 'rat_vert', 'rat_horz'}
        
        %         tone_or_noise= 'tone';
        tone_or_noise= 'noise';
        % load dAM trajectory
%         dam_traj_Hz_struct= load('D:\Dropbox\Pitts_files\Aravind_ExpAM\dAM_files_for_Satya_and_Aravind\Code\approx_am_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        dam_traj_Hz_struct= load('files\dAM_files_for_Satya_and_Aravind\Code\approx_am_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        Xorg= (1:length(dam_traj_Hz_struct.AMfreqvec_est))'/dam_traj_Hz_struct.fs;
        Yorg= dam_traj_Hz_struct.AMfreqvec_est(:);

        % load sitmuli
        fs_sig= 24e3;
%         temp_stim_data= load('D:\Dropbox\Pitts_files\Aravind_ExpAM\dAM_files_for_Satya_and_Aravind\Code\DOD_sweeps_Nov2021.mat');
        temp_stim_data= load('files\dAM_files_for_Satya_and_Aravind\Code\DOD_sweeps_Nov2021.mat');
        if strcmp(tone_or_noise, 'tone')
            stim_wav= temp_stim_data.AMFMtoneup_piece(:);
        elseif strcmp(tone_or_noise, 'noise')
            stim_wav= temp_stim_data.AMFMnoiseup_piece(:);
        end
        fs_stim_org= temp_stim_data.amfmtonesweep.SampleRate;
        stim_wav= helper.gen_resample(stim_wav, fs_stim_org, fs_sig);

%         Species_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\dAM_files_for_Satya_and_Aravind\inData\';
        Species_Dir= 'files\dAM_files_for_Satya_and_Aravind\inData\';
        if strcmp(tone_or_noise, 'tone')
            all_files= dir([Species_Dir '*tone*.mat']);
            freq_lim_kHz= [6, 10];
        elseif strcmp(tone_or_noise, 'noise')
            all_files= dir([Species_Dir '*noise*.mat']);
            freq_lim_kHz= [0, fs_sig/2]/1e3;
        end
        all_files= all_files(~contains({all_files.name}', 'peaks'));

        fs_ffr= 24414.0625;
        mu_pow_lim_dB= 50;

        tMin_s= 0;
        tMax_s= 1;

        if strcmp(Species, 'rat_vert')
            col_deep= helper.get_color('r');
            col_light= helper.get_color('lr');
        elseif strcmp(Species, 'rat_horz')
            col_deep= helper.get_color('prp');
            col_light= helper.get_color('pink');
        end
%         FFR_gain= 1e6/5;
        FFR_gain= 1e6;

    case 'mice'
        % load dAM trajectory
        dam_traj_Hz_struct= load('AMfreqvec_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        Xorg= (1:length(dam_traj_Hz_struct.AMfreqvec_est))'/dam_traj_Hz_struct.fs;
        Yorg= dam_traj_Hz_struct.AMfreqvec_est(:);

        % load sitmuli
        fs_sig= 30e3;
%         [stim_wav, fs_stim_org]= audioread('D:\Dropbox\Pitts_files\Aravind_ExpAM\wavfiles\DoDAMPatternsNEW_12k_mouse_MF2FF2inch.wav');
        [stim_wav, fs_stim_org]= audioread('files\wavfiles\DoDAMPatternsNEW_12k_mouse_MF2FF2inch.wav');
        stim_wav= helper.gen_resample(stim_wav, fs_stim_org, fs_sig);

%         Species_Dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports_Mice\50wk\12k\';
        Species_Dir= 'files\Exports_Mice\50wk\12k\';
        all_files= dir([Species_Dir 'M*']);
        fs_ffr= 1/(40.96e-6);
        freq_lim_kHz= [10, 14];
        mu_pow_lim_dB= 50;

        tMin_s= 0;
        tMax_s= 1;

        col_deep= helper.get_color('g');
        col_light= helper.get_color('lg');
        FFR_gain= 1;
end

stim_wav= stim_wav/max(abs(stim_wav));

% hp_env_filter= helper.get_filter_designfilt('hp', 50, fs_env);

% Loop through all files to load all FFR data
all_data= cell(1,length(all_files));

all_data_frac_tracked= cell(1,length(all_files));
all_data_abs_tracked= cell(1,length(all_files));
all_data_frac_tracked_NF= cell(1,length(all_files));
all_data_abs_tracked_NF= cell(1,length(all_files));


nMin_ind_ffr= max(1, floor(fs_ffr*tMin_s));
nMax_ind_ffr= floor(fs_ffr*tMax_s);

t_ffr= (nMin_ind_ffr:nMax_ind_ffr)/fs_ffr;

% dAM trajectory for gerbils
dam_traj_Hz= interp1(Xorg(:), Yorg(:), t_ffr(:));
dam_traj_Hz(isnan(dam_traj_Hz))= dam_traj_Hz(find(~isnan(dam_traj_Hz), 1)); % make all preceding nan = min non-nan value
dam_traj_Hz(end)= nan;
% dam_traj_Hz_lower= dam_traj_Hz/dAM_NF_traj_factor;
dam_traj_Hz_NFupper= dam_traj_Hz + 3*Filter_HalfWindow_Hz;


bp_filter_ffr= helper.get_filter_designfilt('bp', bp_window_Hz, fs_ffr);

for fileVar=1:length(all_files)
    fStruct= all_files(fileVar);
    fName= [fStruct.folder filesep fStruct.name];

    switch Species
        case {'human', 'gerbil', 'mice'}
            temp_data= read_txt_data(fName);
            temp_data= temp_data(1:nMax_ind_ffr);

        case {'rat', 'rat_vert', 'rat_horz'}
            %             chan2use= 'sagittal';
            if strcmp(Species, 'rat_vert')
                ch1ind.pos= 1;
                ch1ind.neg= 3;
            elseif strcmp(Species, 'rat_horz')
                ch1ind.pos= 2;
                ch1ind.neg= 4;
            end

            cur_data= load(fName);
            cur_data= cur_data.tdat;
            temp_data= (cur_data(ch1ind.pos,:)+cur_data(ch1ind.neg,:))/2; % ch1 tracks high AM better
            temp_data= temp_data(nMin_ind_ffr:nMax_ind_ffr);
    end

    temp_data= FFR_gain * detrend(temp_data(:)); % should be column vector
%     temp_data=  detrend(temp_data(:)); % should be column vector
    all_data{fileVar}= temp_data;

    [all_data_frac_tracked{fileVar}, all_data_abs_tracked{fileVar}]= get_trajectory_hilbert_signal(temp_data, fs_ffr, dam_traj_Hz, Filter_HalfWindow_Hz);
    [all_data_frac_tracked_NF{fileVar}, all_data_abs_tracked_NF{fileVar}]= get_trajectory_hilbert_signal(temp_data, fs_ffr, dam_traj_Hz_NFupper, Filter_HalfWindow_Hz);

end
mean_data= nanmean(cell2mat(all_data), 2);
mean_data= filtfilt(bp_filter_ffr, mean_data);
mean_env= detrend(helper.gen_resample(mean_data, fs_ffr, fs_env));
mean_env_hp= mean_env;
mean_env= detrend(mean_env);
% mean_env= .98*mean_env/max(abs(mean_env));
t_ffr_env= (1:length(mean_env))/fs_env;

fprintf("%s: ptp=%f | rms=%f \n", Species, range(temp_data), std(mean_data))
fprintf("\t after gain: ptp=%f | rms=%f \n", FFR_gain*range(temp_data), FFR_gain*std(mean_data))

if doSaveGrandAvg
    outdir_GrandAvg= 'GrandAvg\';
    out_fName= sprintf('%sGA_%s.mat', outdir_GrandAvg, Species);
    text_note= 'mean_env_hp@fs_env. mean_data&dam_traj_Hz@fs_ffr';
    save(out_fName, 'mean_env_hp', 'fs_env', 'mean_data', 'dam_traj_Hz', 'fs_ffr', 'text_note');
end

min_AM_Hz= 15;
max_AM_Hz= 1202;
tMin_AM_ind= dsearchn(dam_traj_Hz(:), min_AM_Hz)-1;
tMin_AM_s= t_ffr(tMin_AM_ind);
tMax_AM_ind= dsearchn(dam_traj_Hz(:), max_AM_Hz)+1;
tMax_AM_s= t_ffr(tMax_AM_ind);

% get valid inds for gerbils
valid_NF_inds= (t_ffr>tMin_AM_s) & (t_ffr<tMax_AM_s);
dam_traj_Hz_NFvalid= dam_traj_Hz(valid_NF_inds);

if ismember(Species, {'rat', 'rat_vert', 'rat_horz'})
    %% have to do some book-keeping to avoid the weird jump for rats

    rat_seg1_end= find(diff(dam_traj_Hz_NFvalid)<-5); % do this to remove the weird jump down
    rat_seg2_start= find(dam_traj_Hz_NFvalid>dam_traj_Hz_NFvalid(rat_seg1_end), 1 );
    valid_rat_inds= [1:rat_seg1_end, rat_seg2_start:length(dam_traj_Hz_NFvalid)];


end

all_pow_am= cell2mat(all_data_abs_tracked);
all_pow_am_NF= cell2mat(all_data_abs_tracked_NF);

all_pow_am= db(all_pow_am(valid_NF_inds, :));
mu_pow_am= nanmean(all_pow_am, 2);

all_pow_am_NF= db(all_pow_am_NF(valid_NF_inds, :));
mu_pow_am_NF= nanmean(all_pow_am_NF, 2);

%%
SG_range_dB= 40;
doPlot= 1;
time_ticks_vals_ms= 0:500:1000;
time_ticks_labs= num2str(time_ticks_vals_ms(:)/1e3);
lw0= .75;
lw1= 2;
freq_tick_Hz_val= [16, 64, 256, 1200];
freq_tick_Hz_lab= {'16', '64', '256', '1.2k'};

if doPlot
    axes(sp_ax(1));
    helper.plot_spectrogram(stim_wav, fs_sig);
    title(sent_case(Species), 'Color', col_deep, 'Interpreter','none')
    xlabel('');
    if verboseAxis
        ylabel('Freq, kHz');
    else
        ylabel('');
    end
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', '');
    ylim(freq_lim_kHz);

    axes(sp_ax(2));
    plot(t_ffr_env*1e3, mean_env)
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', '');
%     ylim([-1 1])
    xlabel('');
    if verboseAxis
        ylabel('Amplitude (uV)');
    else
        ylabel('');
    end

    axes(sp_ax(3));
    hold on;
    helper.plot_spectrogram(mean_env_hp, fs_env);
    if range(get(colorbar, 'limit'))>SG_range_dB
        caxis([-SG_range_dB 0]+max(get(colorbar, 'limit')))
    end
    line(t_ffr(:)*1e3, (dam_traj_Hz(:))/1e3, 'linestyle', ':', 'linew', 1, 'color', 'r')
    %     line(t_ffr(:)*1e3, (dam_traj_Hz(:)*.8)/1e3, 'linestyle', ':', 'linew', 1, 'color', 'r')
    colorbar off;
    if verboseAxis
        ylabel('Freq, kHz');
    else
        ylabel('');
    end
    xlabel('Time (s)');
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', time_ticks_labs);
    ylim([0 1])

    %% 
    axes(sp_ax(4));
    hold on;
    if ismember(Species, {'human', 'gerbil', 'mice'})
        plot(dam_traj_Hz_NFvalid, all_pow_am, 'Color', col_light, 'LineWidth', lw0)
        % plot(dam_traj_Hz_gerbil_NFvalid, all_pow_am_gerbil_NF, 'Color', helper.get_color('lr'), 'LineWidth', lw0)
        lHan(1)= plot(dam_traj_Hz_NFvalid, mu_pow_am, 'Color', col_deep, 'linew', lw1);
        lHan(2)= plot(dam_traj_Hz_NFvalid, mu_pow_am_NF, 'Color', helper.get_color('lgray'), 'linew', lw1);
    elseif ismember(Species, {'rat', 'rat_vert', 'rat_horz'})
        plot(dam_traj_Hz_NFvalid(valid_rat_inds), all_pow_am(valid_rat_inds, :), 'Color', col_light, 'LineWidth', lw0)
        % plot(dam_traj_Hz_gerbil_NFvalid, all_pow_am_gerbil_NF, 'Color', helper.get_color('lr'), 'LineWidth', lw0)
        lHan(1)= plot(dam_traj_Hz_NFvalid(valid_rat_inds), mu_pow_am(valid_rat_inds), 'Color', col_deep, 'linew', lw1);
        lHan(2)= plot(dam_traj_Hz_NFvalid(valid_rat_inds), mu_pow_am_NF(valid_rat_inds), 'Color', helper.get_color('lgray'), 'linew', lw1);

    end
    set(gca, 'XScale', 'log', 'XTick', freq_tick_Hz_val, 'XTickLabel', '')
    ylim([-mu_pow_lim_dB 0]+max(mu_pow_am(:))+5)
    if verboseAxis
        ylabel('Frac. power, dB');
    else
        ylabel('');
    end

    %% 
    axes(sp_ax(5));
    hold on;
    tri_filt_dur= 25e-3;
    tri_filt_width= round(fs_env*tri_filt_dur);
    if ismember(Species, {'human', 'gerbil', 'mice'})
        plot(dam_traj_Hz_NFvalid, helper.trifilt(mu_pow_am-mu_pow_am_NF, tri_filt_width), 'Color', col_deep, 'LineWidth', lw1)
        plot(dam_traj_Hz_NFvalid, 0*mu_pow_am_NF, 'Color', helper.get_color('lgray'), 'linew', lw1, 'LineStyle', '--');
    elseif ismember(Species, {'rat', 'rat_vert', 'rat_horz'})
        plot(dam_traj_Hz_NFvalid(valid_rat_inds), helper.trifilt(mu_pow_am(valid_rat_inds)-mu_pow_am_NF(valid_rat_inds), tri_filt_width), 'Color', col_deep, 'linew', lw1);
        plot(dam_traj_Hz_NFvalid(valid_rat_inds), 0*mu_pow_am_NF(valid_rat_inds), 'Color', helper.get_color('lgray'), 'linew', lw1, 'LineStyle', '--');

    end
    
    if verboseAxis
        ylabel('SNR, dB');
    else
        ylabel('');
    end

    ylim([-3 30])

    set(gca, 'XScale', 'log', 'XTick', freq_tick_Hz_val, 'XTickLabel', freq_tick_Hz_lab)
    xlabel('dAM freq (Hz)');

end
end

%%
function word_out = sent_case(word_in)
word_out= [upper(word_in(1)), word_in(2:end)];
end

%%
function sp_ax_TF=get_axes(fig_TF, nSProws, nSPcols)

figure(fig_TF);
figSize_cm= [2, 4, 17, 13];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
set(gcf, figure_prop_name, figure_prop_val);
clf;

xc= .05;
yc= .06;
xw= .16;
xs= .03;
yw= .145;
ys= .035;

sp_ax_TF= nan(nSProws*nSPcols, 1);
for spNum=1:(nSProws*nSPcols)
    n_yskips= nSProws - ceil(spNum/nSPcols);
    yc_current= yc + n_yskips*(yw+ys);
    if spNum<=3*nSPcols
        yc_current= yc_current + 1.33*ys;
    end

    n_xskips= rem(spNum-1, nSPcols);
    xc_current= xc + n_xskips*(xw+xs);

    sp_ax_TF(spNum) = subplot(nSProws, nSPcols, spNum);
    set(gca, 'Units', 'normalized', 'Position', [xc_current, yc_current, xw, yw])

end
end