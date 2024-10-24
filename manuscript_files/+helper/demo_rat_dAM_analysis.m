clear;
clc;

doSaveFig= 0;
doSaveGrandAvg= 0;

if strcmp(getenv('COMPUTERNAME'), 'NB-VATS-LAB1')
    fig_dir= 'C:\Users\sap245\Google Drive\PostDoc\dAM_MultiSpecies\Figures\';

else
    fig_dir= 'G:\My Drive\PostDoc\dAM_MultiSpecies\Figures\';
end

if ~isfolder(fig_dir)
    fig_dir = './saved_figures/';
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
end

% CodeDir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\';
CodeDir= '..\';
addpath(CodeDir);

SPpathcode= '.\SP_path_code\';
addpath(SPpathcode);


%% load estimate AM trajectory

fig_TF= 4;
nSProws= 5;
nSPcols= 2;
sp_ax_TF= get_axes(fig_TF, nSProws, nSPcols);


%% human plots: first column
%% rat sagittal plots: first column
verboseAxis= 1;
[dam_traj_Hz_rat_NFvalid_vert, all_pow_am_rat_vert, all_pow_am_rat_NF_vert, mu_pow_am_rat_vert, mu_pow_am_rat_NF_vert]= load_species_data('rat_vert', sp_ax_TF( 1 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg); % rat_vert= rat + vertical (sagittal)

%% rat plots: second column
verboseAxis= 0;
[dam_traj_Hz_rat_NFvalid_horz, all_pow_am_rat_horz, all_pow_am_rat_NF_horz, mu_pow_am_rat_horz, mu_pow_am_rat_NF_horz]= load_species_data('rat_horz', sp_ax_TF( 2 + (0:(nSProws-1))*nSPcols), verboseAxis, doSaveGrandAvg);  % rat_vert= rat + vertical (interaural)

%% fig properties
set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(findall(gcf,'-property','TickDir'),'TickDir', 'both');
set(findall(gcf,'-property','box'),'box', 'off');

%%
rmpath(CodeDir);

%% all functions
function [dam_traj_Hz_NFvalid, all_pow_am, all_pow_am_NF, mu_pow_am, mu_pow_am_NF, lHan]= load_species_data(Species, sp_ax, verboseAxis, doSaveGrandAvg)

fs_env= 2e3;
bp_window_Hz= [14, 1.5e3];
Filter_HalfWindow_Hz= 12;

%% load AM trajectory
switch lower(Species)

    case {'rat', 'rat_vert', 'rat_horz'}

        %% load dAM trajectory
        %         tone_or_noise= 'tone';
        tone_or_noise= 'noise';
        dam_traj_Hz_struct= load('..\dAM_files_for_Satya_and_Aravind\Code\approx_am_est.mat'); % estimated using AMFMforDODpiecewise_APeditAug21
        Xorg= (1:length(dam_traj_Hz_struct.AMfreqvec_est))'/dam_traj_Hz_struct.fs;
        Yorg= dam_traj_Hz_struct.AMfreqvec_est(:);

        %% load sitmuli
        fs_sig= 24e3;
        temp_stim_data= load('..\dAM_files_for_Satya_and_Aravind\Code\DOD_sweeps_Nov2021.mat');
        if strcmp(tone_or_noise, 'tone')
            stim_wav= temp_stim_data.AMFMtoneup_piece(:);
        elseif strcmp(tone_or_noise, 'noise')
            stim_wav= temp_stim_data.AMFMnoiseup_piece(:);
        end
        fs_stim_org= temp_stim_data.amfmtonesweep.SampleRate;
        stim_wav= gen_resample(stim_wav, fs_stim_org, fs_sig);

        %% get data filenames to go through and specific data parameters like gain, sampling freq, color names to plot 
        Species_Dir= '..\dAM_files_for_Satya_and_Aravind\inData\';
        if strcmp(tone_or_noise, 'tone')
            all_files= dir([Species_Dir '*tone*.mat']);
            freq_lim_kHz= [6, 10];
        elseif strcmp(tone_or_noise, 'noise')
            all_files= dir([Species_Dir '*noise*.mat']);
            freq_lim_kHz= [0, fs_sig/2]/1e3;
        end
        all_files= all_files(~contains({all_files.name}', 'peaks'));

        fs_ffr= 24414.0625;
        FFR_gain= 1e6/5;
        mu_pow_lim_dB= 50;

        tMin_s= 0;
        tMax_s= 1;

        if strcmp(Species, 'rat_vert')
            col_deep= get_color('r');
            col_light= get_color('lr');
        elseif strcmp(Species, 'rat_horz')
            col_deep= get_color('prp');
            col_light= get_color('pink');
        end


end

% normalize stim 
stim_wav= stim_wav/max(abs(stim_wav));


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


bp_filter_ffr= get_filter_designfilt('bp', bp_window_Hz, fs_ffr);

for fileVar=1:length(all_files)
    fStruct= all_files(fileVar);
    fName= [fStruct.folder filesep fStruct.name];

    switch Species

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
    all_data{fileVar}= temp_data;

    [all_data_frac_tracked{fileVar}, all_data_abs_tracked{fileVar}]= get_trajectory_hilbert_signal(temp_data, fs_ffr, dam_traj_Hz, Filter_HalfWindow_Hz);
    [all_data_frac_tracked_NF{fileVar}, all_data_abs_tracked_NF{fileVar}]= get_trajectory_hilbert_signal(temp_data, fs_ffr, dam_traj_Hz_NFupper, Filter_HalfWindow_Hz);

end
mean_data= nanmean(cell2mat(all_data), 2);
mean_data= filtfilt(bp_filter_ffr, mean_data);
mean_env= detrend(gen_resample(mean_data, fs_ffr, fs_env));
mean_env_hp= mean_env;
mean_env= detrend(mean_env);
mean_env= .98*mean_env/max(abs(mean_env));
t_ffr_env= (1:length(mean_env))/fs_env;

min_AM_Hz= 15;
max_AM_Hz= 1202;
tMin_AM_ind= dsearchn(dam_traj_Hz(:), min_AM_Hz)-1;
tMin_AM_s= t_ffr(tMin_AM_ind);
tMax_AM_ind= dsearchn(dam_traj_Hz(:), max_AM_Hz)+1;
tMax_AM_s= t_ffr(tMax_AM_ind);

valid_NF_inds= (t_ffr>tMin_AM_s) & (t_ffr<tMax_AM_s);
dam_traj_Hz_NFvalid= dam_traj_Hz(valid_NF_inds);

if ismember(Species, {'rat', 'rat_vert', 'rat_horz'})
    %% have to do some book-keeping to avoid the jump for rats

    rat_seg1_end= find(diff(dam_traj_Hz_NFvalid)<-5); % do this to remove the jump down
    rat_seg2_start= find(dam_traj_Hz_NFvalid>dam_traj_Hz_NFvalid(rat_seg1_end), 1 );
    valid_rat_inds= [1:rat_seg1_end, rat_seg2_start:length(dam_traj_Hz_NFvalid)];

end


all_pow_am= cell2mat(all_data_abs_tracked);
all_pow_am_NF= cell2mat(all_data_abs_tracked_NF);

all_pow_am= db(all_pow_am(valid_NF_inds, :));
mu_pow_am= nanmean(all_pow_am, 2);

all_pow_am_NF= db(all_pow_am_NF(valid_NF_inds, :));
mu_pow_am_NF= nanmean(all_pow_am_NF, 2);


%% plot stim spectrogram, EFR waveform, EFR spectraogram, frac power, and SNR 
SG_range_dB= 40;
doPlot= 1;
time_ticks_vals_ms= 0:500:1000;
time_ticks_labs= num2str(time_ticks_vals_ms(:)/1e3);
lw0= .75;
lw1= 2;
freq_tick_Hz_val= [16, 64, 256, 1200];
freq_tick_Hz_lab= {'16', '64', '256', '1.2k'};

if doPlot
    % 1) stim spectrogram
    axes(sp_ax(1));
    plot_spectrogram(stim_wav, fs_sig);
    title(sent_case(Species), 'Color', col_deep, 'Interpreter','none')
    xlabel('');
    if verboseAxis
        ylabel('Freq, kHz');
    else
        ylabel('');
    end
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', '');
    ylim(freq_lim_kHz);

    % 2) EFR waveform
    axes(sp_ax(2));
    plot(t_ffr_env*1e3, mean_env)
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', '');
    ylim([-1 1])
    xlabel('');
    if verboseAxis
        ylabel('Amplitude (au)');
    else
        ylabel('');
    end

    % 3) EFR spectraogram
    axes(sp_ax(3));
    hold on;
    plot_spectrogram(mean_env_hp, fs_env);
    if range(get(colorbar, 'limit'))>SG_range_dB
        caxis([-SG_range_dB 0]+max(get(colorbar, 'limit')))
    end
    line(t_ffr(:)*1e3, (dam_traj_Hz(:))/1e3, 'linestyle', ':', 'linew', 1, 'color', 'r')
    colorbar off;
    if verboseAxis
        ylabel('Freq, kHz');
    else
        ylabel('');
    end
    xlabel('Time (s)');
    set(gca, 'XTick', time_ticks_vals_ms, 'XTickLabel', time_ticks_labs);
    ylim([0 1])

    % 4) frac power
    axes(sp_ax(4));
    hold on;

    plot(dam_traj_Hz_NFvalid(valid_rat_inds), all_pow_am(valid_rat_inds, :), 'Color', col_light, 'LineWidth', lw0)
    lHan(1)= plot(dam_traj_Hz_NFvalid(valid_rat_inds), mu_pow_am(valid_rat_inds), 'Color', col_deep, 'linew', lw1);
    lHan(2)= plot(dam_traj_Hz_NFvalid(valid_rat_inds), mu_pow_am_NF(valid_rat_inds), 'Color', get_color('lgray'), 'linew', lw1);

    set(gca, 'XScale', 'log', 'XTick', freq_tick_Hz_val, 'XTickLabel', '')
    ylim([-mu_pow_lim_dB 0]+max(mu_pow_am(:))+5)
    if verboseAxis
        ylabel('Frac. power, dB');
    else
        ylabel('');
    end

    % 5) EFR SNR along dAM 
    axes(sp_ax(5));
    hold on;
    tri_filt_dur= 25e-3;
    tri_filt_width= round(fs_env*tri_filt_dur);
    if ismember(Species, {'human', 'gerbil', 'mice'})
        plot(dam_traj_Hz_NFvalid, trifilt(mu_pow_am-mu_pow_am_NF, tri_filt_width), 'Color', col_deep, 'LineWidth', lw1)
        plot(dam_traj_Hz_NFvalid, 0*mu_pow_am_NF, 'Color', get_color('lgray'), 'linew', lw1, 'LineStyle', '--');
    elseif ismember(Species, {'rat', 'rat_vert', 'rat_horz'})
        plot(dam_traj_Hz_NFvalid(valid_rat_inds), trifilt(mu_pow_am(valid_rat_inds)-mu_pow_am_NF(valid_rat_inds), tri_filt_width), 'Color', col_deep, 'linew', lw1);
        plot(dam_traj_Hz_NFvalid(valid_rat_inds), 0*mu_pow_am_NF(valid_rat_inds), 'Color', get_color('lgray'), 'linew', lw1, 'LineStyle', '--');

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
figSize_cm= [5, 4, 12, 13];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
set(gcf, figure_prop_name, figure_prop_val);
clf;

xc= .1;
yc= .08;
xw= .39;
xs= .08;
yw= .14;
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