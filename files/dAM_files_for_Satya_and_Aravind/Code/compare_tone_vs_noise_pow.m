clear;
clc;

doSaveFig = 0;
all_pow_data= load('../outData/PowerValuesTime500_Infms__Filter_80_InfHz.mat');
all_pow_data= all_pow_data.all_pow_data;

animal_ids= extractBefore({all_pow_data.filename}, '_amfm')';
uniq_animal_ids= unique(animal_ids);

tone_file_inds= contains({all_pow_data.filename}, 'tone')';
noise_file_inds= contains({all_pow_data.filename}, 'noise')';


%% estimate paired power
ch1data= repmat(struct('animal_id', '', 'tone_ch1pow_sum_frac', nan, 'noise_ch1pow_sum_frac', nan, 'tone2noise_sum_ch1pow', nan, ...
    'tone_ch1pow_diff_frac', nan, 'noise_ch1pow_diff_frac', nan, 'tone2noise_diff_ch1pow', nan), length(uniq_animal_ids), 1);
ch2data= repmat(struct('animal_id', '', 'tone_ch2pow_sum_frac', nan, 'noise_ch2pow_sum_frac', nan, 'tone2noise_sum_ch2pow', nan, ...
    'tone_ch2pow_diff_frac', nan, 'noise_ch2pow_diff_frac', nan, 'tone2noise_diff_ch2pow', nan), length(uniq_animal_ids), 1);
for uniqIDvar=1:length(uniq_animal_ids)
    cur_animal= uniq_animal_ids{uniqIDvar};
    cur_tone_ind= strcmp(animal_ids, cur_animal) & tone_file_inds;
    cur_noise_ind= strcmp(animal_ids, cur_animal) & noise_file_inds;

    cur_tone_data= all_pow_data(cur_tone_ind);
    cur_noise_data= all_pow_data(cur_noise_ind);

    ch1data(uniqIDvar)= struct('animal_id', cur_animal, ...
        'tone_ch1pow_sum_frac', cur_tone_data.ch1sum_pow_frac, 'noise_ch1pow_sum_frac', cur_noise_data.ch1sum_pow_frac, 'tone2noise_sum_ch1pow', cur_tone_data.ch1sum_pow_frac/cur_noise_data.ch1sum_pow_frac, ...
        'tone_ch1pow_diff_frac', cur_tone_data.ch1diff_pow_frac, 'noise_ch1pow_diff_frac', cur_noise_data.ch1diff_pow_frac, 'tone2noise_diff_ch1pow', cur_tone_data.ch1diff_pow_frac/cur_noise_data.ch1diff_pow_frac);

    ch2data(uniqIDvar)= struct('animal_id', cur_animal, ...
        'tone_ch2pow_sum_frac', cur_tone_data.ch2sum_pow_frac, 'noise_ch2pow_sum_frac', cur_noise_data.ch2sum_pow_frac, 'tone2noise_sum_ch2pow', cur_tone_data.ch2sum_pow_frac/cur_noise_data.ch2sum_pow_frac, ...
        'tone_ch2pow_diff_frac', cur_tone_data.ch2diff_pow_frac, 'noise_ch2pow_diff_frac', cur_noise_data.ch2diff_pow_frac, 'tone2noise_diff_ch2pow', cur_tone_data.ch2diff_pow_frac/cur_noise_data.ch2diff_pow_frac);
end

%%
figSize_cm= [55 5 25 15];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
set(gcf,figure_prop_name,figure_prop_val);

figure(1);
clf;

sp_ax(1)= subplot(231);
boxplot([[ch1data.tone_ch1pow_sum_frac]', [ch1data.noise_ch1pow_sum_frac]'], {'Tone', 'Noise'})
ylabel('Ch1 Frac Sum power');

sp_ax(2)= subplot(232);
boxplot([[ch2data.tone_ch2pow_sum_frac]', [ch2data.noise_ch2pow_sum_frac]'], {'Tone', 'Noise'})
ylabel('Ch2 Frac Sum power');

sp_ax(3)= subplot(233);
boxplot([[ch1data.tone2noise_sum_ch1pow]', [ch2data.tone2noise_sum_ch2pow]'], {'Ch1Sum', 'Ch2Sum'})
ylabel('SUM: Tone to noise power')

sp_ax(4)= subplot(234);
boxplot([[ch1data.tone_ch1pow_diff_frac]', [ch1data.noise_ch1pow_diff_frac]'], {'Tone', 'Noise'})
ylabel('Ch1 Frac Diff power');

sp_ax(5)= subplot(235);
boxplot([[ch2data.tone_ch2pow_diff_frac]', [ch2data.noise_ch2pow_diff_frac]'], {'Tone', 'Noise'})
ylabel('Ch2 Frac Diff power');

sp_ax(6)= subplot(236);
boxplot([[ch1data.tone2noise_sum_ch1pow]', [ch2data.tone2noise_sum_ch2pow]'], {'Ch1Sum', 'Ch2Sum'})
ylabel('DIFF: Tone to noise power');


%%
xc=.06;
xw= .25;
xs= .08;
yc= .07;
ys= .1;
yw= .4;

set(sp_ax(4), 'Units', 'normalized', 'Position', [xc, yc, xw, yw]);
set(sp_ax(5), 'Units', 'normalized', 'Position', [xc+xs+xw, yc, xw, yw]);
set(sp_ax(6), 'Units', 'normalized', 'Position', [xc+2*xs+2*xw, yc, xw, yw]);
set(sp_ax(1), 'Units', 'normalized', 'Position', [xc, yc+ys+yw, xw, yw]);
set(sp_ax(2), 'Units', 'normalized', 'Position', [xc+xs+xw, yc+ys+yw, xw, yw]);
set(sp_ax(3), 'Units', 'normalized', 'Position', [xc+2*xs+2*xw, yc+ys+yw, xw, yw]);

%%
OutFigDir= '../outData/';
fig_name= sprintf('%sFig_ToneVNoise_Pow', OutFigDir);
if doSaveFig
    print(fig_name, '-dpng', '-r600')
end