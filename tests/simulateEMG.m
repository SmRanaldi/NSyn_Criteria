function [emg]=simulateEMG(Fh,Fl,Fs,wTrue,s)

% [emg]=generaEMG(Fh,Fl,Fs,wTrue,s)
%
% Generates a modulated EMG trace using the Stulen & De Luca filter.
%
% INPUTS:
% 	Fh: Low pass cut-off frequency.
% 	Fl: High pass cut-off frequancy.
% 	Fs: Sampling frequency.
% 	wTrue: Modulating waveform.
% 	s: White noise series.
%
% OUTPUTS:
% 	emg: The simulated EMG signal.

%% Initialization.

wTrue=[zeros(1,5000),wTrue,zeros(1,5000)]; % Needed for the filtering.
s=[randn(1,5000),s,randn(1,5000)];

%% Filter parameters.

a = fdeluca(Fh,Fl,Fs); b=1;

%% EMG generation.

s=filtfilt(a,b,s);

emg=(s.*wTrue)/std(s);

emg=emg(5001:end-5000);
