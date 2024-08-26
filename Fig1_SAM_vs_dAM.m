clear;
clc;

do_save_fig= 0;

%% Common parameters 
fs= 20e3;
stimdur_s= 1;
stimlen= round(fs*stimdur_s);
fc_Hz= 3e3;
mod_depth= 1;
t_stim_ms= (1:stimlen)/stimlen*stimdur_s*1e3;

%% freq params 
fm_Hz= .5e3;
dam_f_start= 2^4;
dam_f_end= 2^10.25;
dam_traj_Hz= 2.^linspace(log2(dam_f_start), log2(dam_f_end), stimlen);

t_dam_x_am_ms= t_stim_ms(dsearchn(dam_traj_Hz(:), fm_Hz));
t_zoom_edge_ms= t_dam_x_am_ms + [-3 3]/fm_Hz*1e3;
t_zoom_inds= t_stim_ms>min(t_zoom_edge_ms) & t_stim_ms<max(t_zoom_edge_ms);

%% create SAM 
sam_stim= (1+mod_depth*sin(2*pi*fm_Hz*t_stim_ms/1e3)).*sin(2*pi*fc_Hz*t_stim_ms/1e3);
sam_stim= sam_stim/max(abs(sam_stim));
sam_stim_env= abs(hilbert(sam_stim));
%% create dAM 

phi_instantaneous= -cumtrapz(dam_traj_Hz)/fs; % Formula for phase in terms of frequency
dam_Mod= cos(2*pi*phi_instantaneous);
dam_Car= sin(2*pi*fc_Hz*t_stim_ms/1e3);
dam_stim= (1-mod_depth*dam_Mod(:))/2.*dam_Car(:);
dam_stim_env= abs(hilbert(dam_stim));

%% define axes 
nSProws= 3;
nSPcols= 4;
carrier_freq_ticks_kHz= 0:5;
mod_freq_tick_kHz= 0:0.5:2;


lw1= 1;
lw2= 1.5;

figSize_cm= [27 3 15 8.5];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
set(gcf,figure_prop_name,figure_prop_val);

clf;

sp_ax(1)= subplot(nSProws, nSPcols, 1);
hold on
plot(t_stim_ms, sam_stim, 'LineWidth', lw1);
plot(t_stim_ms, sam_stim_env, 'LineWidth', lw2);
xlabel('Time (ms)');
ylabel('Amp. (au)');
ttl_han(1)= text(1.05, 1.2, 'SAM tone', 'Units', 'normalized', 'HorizontalAlignment', 'center');
panel_han(1)= text(-0.24, 1.22, 'A', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');


sp_ax(2)= subplot(nSProws, nSPcols, 2);
hold on
plot(t_stim_ms(t_zoom_inds), sam_stim(t_zoom_inds), 'LineWidth', lw1);
plot(t_stim_ms(t_zoom_inds), sam_stim_env(t_zoom_inds), 'LineWidth', lw2);
xlim(t_zoom_edge_ms)
xlabel('Time (ms)');

sp_ax(5)= subplot(nSProws, nSPcols, 5);
helper.plot_dAM_spectrogram(sam_stim, fs);
set(gca, 'YTick', carrier_freq_ticks_kHz, 'XTickLabel', '')
xlabel('');
panel_han(2)= text(-0.24, 1.1, 'B', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

sp_ax(6)= subplot(nSProws, nSPcols, 6);
hold on
[dft_sam_dB, Freq_sam_Hz]= helper.plot_dft(sam_stim, fs, 'plot', false);
plot(dft_sam_dB, Freq_sam_Hz/1e3);
xlim([-50 5]+max(dft_sam_dB));
set(gca, 'YTick', carrier_freq_ticks_kHz, 'XTickLabel', '')

sp_ax(9)= subplot(nSProws, nSPcols, 9);
helper.plot_dAM_spectrogram(sam_stim_env, fs);
set(gca, 'YTick', mod_freq_tick_kHz)
panel_han(3)= text(-0.24, 1.1, 'C', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

sp_ax(10)= subplot(nSProws, nSPcols, 10);
[dft_sam_env_dB, Freq_sam_env_Hz]= helper.plot_dft(sam_stim_env, fs, 'plot', false);
plot(dft_sam_env_dB, Freq_sam_env_Hz/1e3);
xlim([-50 5]+max(dft_sam_env_dB));
set(gca, 'YTick', mod_freq_tick_kHz)
xlabel('DFT Amp (dB)');

%% 
sp_ax(3)= subplot(nSProws, nSPcols, 3);
hold on
plot(t_stim_ms, dam_stim, 'LineWidth', lw1);
plot(t_stim_ms, dam_stim_env, 'LineWidth', lw2);
xlabel('Time (ms)');
ylabel('Amp. (au)');
ttl_han(2)= text(1.05, 1.2, 'dAM tone', 'Units', 'normalized', 'HorizontalAlignment', 'center');
panel_han(4)= text(-0.24, 1.22, 'D', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

sp_ax(4)= subplot(nSProws, nSPcols, 4);
hold on
plot(t_stim_ms(t_zoom_inds), dam_stim(t_zoom_inds), 'LineWidth', lw1);
plot(t_stim_ms(t_zoom_inds), dam_stim_env(t_zoom_inds), 'LineWidth', lw2);
xlim(t_zoom_edge_ms)
xlabel('Time (ms)');

sp_ax(7)= subplot(nSProws, nSPcols, 7);
helper.plot_dAM_spectrogram(dam_stim, fs);
set(gca, 'YTick', carrier_freq_ticks_kHz, 'XTickLabel', '')
xlabel('');
panel_han(5)= text(-0.24, 1.1, 'E', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

sp_ax(8)= subplot(nSProws, nSPcols, 8);
hold on
[dft_dam_dB, Freq_dam_Hz]= helper.plot_dft(dam_stim, fs, 'plot', false);
plot(dft_dam_dB, Freq_dam_Hz/1e3);
xlim([-50 5]+max(dft_dam_dB));
set(gca, 'YTick', carrier_freq_ticks_kHz, 'XTickLabel', '')

sp_ax(11)= subplot(nSProws, nSPcols, 11);
helper.plot_dAM_spectrogram(dam_stim_env, fs);
set(gca, 'YTick', mod_freq_tick_kHz)
panel_han(6)= text(-0.24, 1.1, 'F', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

sp_ax(12)= subplot(nSProws, nSPcols, 12);
[dft_dam_env_dB, Freq_dam_env_Hz]= helper.plot_dft(dam_stim_env, fs, 'plot', false);
plot(dft_dam_env_dB, Freq_dam_env_Hz/1e3);
xlim([-50 5]+max(dft_dam_env_dB));
set(gca, 'YTick', mod_freq_tick_kHz)
xlabel('DFT Amp (dB)');

%% link axis
linkaxes(sp_ax(1:2:11), 'x')
xlim(sp_ax(1), [0 stimdur_s*1e3])

linkaxes(sp_ax(1:4), 'y')
ylim(sp_ax(1), [-1.2 1.2]);

linkaxes(sp_ax(5:8), 'y')
ylim(sp_ax(5), [0 5])

linkaxes(sp_ax(9:12), 'y')
ylim(sp_ax(9), [0 2])

set(findall(gcf,'-property','FontSize'),'FontSize', 8);
set(findall(gcf,'-property','TickDir'),'TickDir', 'both');
set(findall(gcf,'-property','box'),'box', 'off');
set(ttl_han, 'FontSize', 11);
set(panel_han, 'FontSize', 11);
%% 
xc1= .057;
xs= .035;
xs2= .06;
xw= .201;
xc2= xc1+2*xw+xs+xs2;

yc= .11;
yw2= .25;
yw1= .15;
ys1= .05;
ys2= .14;

set(sp_ax(9), 'Units', 'normalized', 'Position', [xc1, yc, xw, yw2])
set(sp_ax(10), 'Units', 'normalized', 'Position', [xc1+xs+xw, yc, xw, yw2])
set(sp_ax(11), 'Units', 'normalized', 'Position', [xc2, yc, xw, yw2])
set(sp_ax(12), 'Units', 'normalized', 'Position', [xc2+xw+xs, yc, xw, yw2])

set(sp_ax(5), 'Units', 'normalized', 'Position', [xc1, yc+yw2+ys1, xw, yw2])
set(sp_ax(6), 'Units', 'normalized', 'Position', [xc1+xs+xw, yc+yw2+ys1, xw, yw2])
set(sp_ax(7), 'Units', 'normalized', 'Position', [xc2, yc+yw2+ys1, xw, yw2])
set(sp_ax(8), 'Units', 'normalized', 'Position', [xc2+xw+xs, yc+yw2+ys1, xw, yw2])

set(sp_ax(1), 'Units', 'normalized', 'Position', [xc1, yc+2*yw2+ys2+ys1, xw, yw1])
set(sp_ax(2), 'Units', 'normalized', 'Position', [xc1+xs+xw, yc+2*yw2+ys2+ys1, xw, yw1])
set(sp_ax(3), 'Units', 'normalized', 'Position', [xc2, yc+2*yw2+ys2+ys1, xw, yw1])
set(sp_ax(4), 'Units', 'normalized', 'Position', [xc2+xw+xs, yc+2*yw2+ys2+ys1, xw, yw1])

%% 
if strcmp(getenv('COMPUTERNAME'), 'NB-VATS-LAB1')
    fig_dir= 'C:\Users\sap245\Google Drive\PostDoc\dAM_MultiSpecies\Figures\';
elseif strcmp(getenv('COMPUTERNAME'), 'DESKTOP-SR5PJKK')
    fig_dir= 'D:\Dropbox\Pitts_files\Aravind_ExpAM\___manuscript\saved_figures\'
else 
    fig_dir= 'G:\My Drive\PostDoc\dAM_MultiSpecies\Figures\';
end

fig_name= [fig_dir 'Fig1_sam_vs_dam'];
if do_save_fig
    print(fig_name, '-dpng', '-r600')
    print([fig_name '_ink'], '-dsvg')
end