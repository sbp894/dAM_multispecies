clear;
clc;

doSaveFig = 0;
inDir= '../inData/';
all_pow_data= load('../outData/PowerValuesTime500_Infms__Filter_80_InfHz.mat');
all_pow_data= all_pow_data.all_pow_data;

all_files= {all_pow_data.filename}';
animal_ids= extractBefore({all_pow_data.filename}, '_')';
uniq_animal_ids= unique(animal_ids);

% For each animal, conditions can be 
% - either tone or noise 
% - either day#0 or day#28

% For each condition, we have the following power metrics 
% - ch1 or ch2 power  
% - sum or diff power  
% - abs or frac power 


tone_file_inds= contains({all_pow_data.filename}, 'tone')';
noise_file_inds= contains({all_pow_data.filename}, 'noise')';


%%

ch1pow_data= repmat(struct('animal_id', '', 'tone_day0_sum_frac', nan, 'noise_day0_sum_frac', nan, ...
    'tone_day28_sum_frac', nan, 'noise_day28_sum_frac', nan, ...
    'tone_day0_diff_frac', nan, 'noise_day0_diff_frac', nan, ...
    'tone_day28_diff_frac', nan, 'noise_day28_diff_frac', nan), ...
    length(uniq_animal_ids), 1);

ch2pow_data= repmat(struct('animal_id', '', 'tone_day0_sum_frac', nan, 'noise_day0_sum_frac', nan, ...
    'tone_day28_sum_frac', nan, 'noise_day28_sum_frac', nan, ...
    'tone_day0_diff_frac', nan, 'noise_day0_diff_frac', nan, ...
    'tone_day28_diff_frac', nan, 'noise_day28_diff_frac', nan), ...
    length(uniq_animal_ids), 1);

[ch1pow_data.animal_id]= uniq_animal_ids{:};
[ch2pow_data.animal_id]= uniq_animal_ids{:};

for fileVar=1:length(all_files)
    cur_fName= all_files{fileVar};
    cur_powdata= all_pow_data(fileVar);
    
    uniq_animal_ind= strcmp(uniq_animal_ids, animal_ids{fileVar});
    
    % tone or noise 
    if contains(cur_fName, 'tone')
        TN_fix= 'tone';
    elseif contains(cur_fName, 'noise')
        TN_fix= 'noise';
    else
        error('Tone or noise? ');
    end
    
    % day 0 or day 1
    if contains(cur_fName, 'day0')
        Day_fix= 'day0';
    elseif contains(cur_fName, 'day28')
        Day_fix= 'day28';
    else
        error('Day 0 or 28? ');
    end
    
    ch1pow_data(uniq_animal_ind).(sprintf('%s_%s_sum_frac', TN_fix, Day_fix))= cur_powdata.ch1sum_pow_frac;
    ch1pow_data(uniq_animal_ind).(sprintf('%s_%s_diff_frac', TN_fix, Day_fix))= cur_powdata.ch1diff_pow_frac; 
    
    ch2pow_data(uniq_animal_ind).(sprintf('%s_%s_sum_frac', TN_fix, Day_fix))= cur_powdata.ch2sum_pow_frac;
    ch2pow_data(uniq_animal_ind).(sprintf('%s_%s_diff_frac', TN_fix, Day_fix))= cur_powdata.ch2diff_pow_frac; 
end

%%
figure(1);
clf;

subplot(221);
boxplot([[ch1pow_data.tone_day0_sum_frac]', [ch1pow_data.tone_day28_sum_frac]'], {'day0', 'day28'})
title('Tone_sum', 'Interpreter', 'none')
text(1.0, 1.15, 'Channel #1', 'Units', 'normalized', 'FontSize', 14);

subplot(222);
boxplot([[ch1pow_data.tone_day0_diff_frac]', [ch1pow_data.tone_day28_diff_frac]'], {'day0', 'day28'})
title('Tone_diff', 'Interpreter', 'none')

subplot(223);
boxplot([[ch1pow_data.noise_day0_sum_frac]', [ch1pow_data.noise_day28_sum_frac]'], {'day0', 'day28'})
ylabel('Power');
title('Noise_sum', 'Interpreter', 'none')

subplot(224);
boxplot([[ch1pow_data.noise_day0_diff_frac]', [ch1pow_data.noise_day28_diff_frac]'], {'day0', 'day28'})
title('Noise_diff', 'Interpreter', 'none')

%%
figure(2);
clf;

subplot(221);
boxplot([[ch2pow_data.tone_day0_sum_frac]', [ch2pow_data.tone_day28_sum_frac]'], {'day0', 'day28'})
title('Tone_sum', 'Interpreter', 'none')
text(1.0, 1.15, 'Channel #2', 'Units', 'normalized', 'FontSize', 14);

subplot(222);
boxplot([[ch2pow_data.tone_day0_diff_frac]', [ch2pow_data.tone_day28_diff_frac]'], {'day0', 'day28'})
title('Tone_diff', 'Interpreter', 'none')

subplot(223);
boxplot([[ch2pow_data.noise_day0_sum_frac]', [ch2pow_data.noise_day28_sum_frac]'], {'day0', 'day28'})
ylabel('Power');
title('Noise_sum', 'Interpreter', 'none')

subplot(224);
boxplot([[ch2pow_data.noise_day0_diff_frac]', [ch2pow_data.noise_day28_diff_frac]'], {'day0', 'day28'})
title('Noise_diff', 'Interpreter', 'none')