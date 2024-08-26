% All units are in Hz/s unless stated otherwise
function [stim, dam_traj_Hz, tStim]= create_dAM_stim(fs, stimdur, fc, dam_f_start, dam_f_end)

stimlen= round(fs*stimdur); % stimulus length in samples 

mod_depth= 1; % modulation depth 
tStim= (1:stimlen)/stimlen; % time 

dam_traj_Hz= 2.^linspace(log2(dam_f_start), log2(dam_f_end), stimlen); % Desired dAM trajectory 

phi_instantaneous= -cumtrapz(dam_traj_Hz)/fs; % Formula for phase (in terms of frequency)
stimMod= cos(2*pi*phi_instantaneous);
stimCar= sin(2*pi*fc*tStim);
stim= (1-mod_depth*stimMod(:))/2.*stimCar(:);

end