<<<<<<< HEAD
% Filter object can be used for filtfilt
% function Hd= get_filter_designfilt(filtType, freqWindow, fs, filtOrder, plotYes)
function Hd= get_filter_designfilt(filtType, freqWindow, fs, filtOrder, plotYes)

if ~exist('filtOrder', 'var')
    filtOrder= 2;
end
if ~exist('plotYes', 'var')
    plotYes =0;
end

switch lower(filtType)
    case {'band', 'bandpass', 'bp'}
        Hd = designfilt('bandpassiir','FilterOrder',filtOrder, ...
            'HalfPowerFrequency1',freqWindow(1),'HalfPowerFrequency2',freqWindow(2), ...
            'SampleRate',fs);
    case {'low', 'lowpass', 'lp'}
        Hd = designfilt('lowpassiir','FilterOrder',filtOrder, ...
            'PassbandFrequency',freqWindow, 'SampleRate',fs);
    case {'high', 'highpass', 'hp'}
        Hd = designfilt('highpassiir','FilterOrder',filtOrder, ...
            'PassbandFrequency',freqWindow, 'SampleRate',fs);
end

if plotYes
    freqz(Hd, 2^12, fs);
=======
% Filter object can be used for filtfilt
% function Hd= get_filter_designfilt(filtType, freqWindow, fs, filtOrder, plotYes)
function Hd= get_filter_designfilt(filtType, freqWindow, fs, filtOrder, plotYes)

if ~exist('filtOrder', 'var')
    filtOrder= 2;
end
if ~exist('plotYes', 'var')
    plotYes =0;
end

switch lower(filtType)
    case {'band', 'bandpass', 'bp'}
        Hd = designfilt('bandpassiir','FilterOrder',filtOrder, ...
            'HalfPowerFrequency1',freqWindow(1),'HalfPowerFrequency2',freqWindow(2), ...
            'SampleRate',fs);
    case {'low', 'lowpass', 'lp'}
        Hd = designfilt('lowpassiir','FilterOrder',filtOrder, ...
            'PassbandFrequency',freqWindow, 'SampleRate',fs);
    case {'high', 'highpass', 'hp'}
        Hd = designfilt('highpassiir','FilterOrder',filtOrder, ...
            'PassbandFrequency',freqWindow, 'SampleRate',fs);
end

if plotYes
    freqz(Hd, 2^12, fs);
>>>>>>> 033e1e8208ba0bc9e95b7a551303ab1999f126dc
end