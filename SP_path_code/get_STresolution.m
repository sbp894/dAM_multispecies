function [total_ST_resolution, total_power_captured, txt_han]= get_STresolution(stim_env, fs, dam_traj_Hz, tWindow, do_plot_boxes)

stimdur= length(stim_env)/fs;

freq_res_default= 1/tWindow;
nSegs= floor(stimdur/tWindow);
min_power_ratio= 0.95;

current_pow_box= nan(nSegs, 1);
current_pow_total= nan(nSegs, 1);
all_freq_res= nan(nSegs, 1);

for segVar=1:nSegs
    freq_res_factor= 1;
    power_ratio= 0;
    cur_start_time_s= (segVar-1)*tWindow;
    cur_end_time_s= segVar*tWindow;

    cur_start_ind= max(1, round(cur_start_time_s*fs));
    cur_end_ind= min(length(stim_env), round(cur_end_time_s*fs));

    cur_freq_start_Hz= dam_traj_Hz(cur_start_ind);
    cur_freq_end_Hz= dam_traj_Hz(cur_end_ind);

    cur_sig= stim_env(cur_start_ind:cur_end_ind);
    [cur_psd_dB, psd_freq_Hz]= plot_dpss_psd(cur_sig, fs, 'nw', 1.5, 'plot', false);
    cur_psd_pow= db2pow(cur_psd_dB);

    while power_ratio<min_power_ratio
        freq_res= freq_res_default*freq_res_factor;

        freq_box_lower_Hz= max(0, cur_freq_start_Hz-freq_res);
        freq_box_upper_Hz= min(fs/2, cur_freq_end_Hz+freq_res);

        box_psd_inds= psd_freq_Hz>freq_box_lower_Hz & psd_freq_Hz<freq_box_upper_Hz;

        current_pow_box(segVar)= sum(cur_psd_pow(box_psd_inds));
        current_pow_total(segVar)= sum(cur_psd_pow);
        all_freq_res(segVar)= freq_res_factor*range(psd_freq_Hz(box_psd_inds));

        power_ratio= sum(cur_psd_pow(box_psd_inds))/sum(cur_psd_pow);
        freq_res_factor= freq_res_factor+0.5;
    end

    do_debug_plot= 0;
    if do_debug_plot
        figure(41);
        clf;
        hold on;
        plot(psd_freq_Hz, cur_psd_dB);
        plot([max(psd_freq_Hz(1), freq_box_lower_Hz), freq_box_upper_Hz], [1, 1]+max(cur_psd_dB), 'r');
        set(gca, 'XScale', 'log')
    end

    box_x= [cur_start_time_s, cur_end_time_s, cur_end_time_s, cur_start_time_s, cur_start_time_s]*1e3;
    box_y= [freq_box_lower_Hz, freq_box_lower_Hz, freq_box_upper_Hz, freq_box_upper_Hz, freq_box_lower_Hz]/1e3;

    if do_plot_boxes
        line(box_x, box_y, 'color', 'k', 'linew', .5, 'linestyle', '-')
    end
end

total_power_captured= sum(current_pow_box)./sum(current_pow_total);
total_ST_resolution= tWindow*sum(all_freq_res);

if do_plot_boxes
    %     title(sprintf('tWindow=%.0f ms, tot_pow=%.0f%%, STres= %s', ...
    %         tWindow*1e3, sum(total_power_captured)*100, num2str(total_ST_resolution)), 'Interpreter', 'none');
    txt_han= text(.075, .8, sprintf('Window=%.0f ms', tWindow*1e3), 'Units', 'normalized', 'Color', 'w');
else
    txt_han= nan;
end