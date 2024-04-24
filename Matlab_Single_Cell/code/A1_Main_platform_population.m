clear all
close all
clc

%% Parameters

freq              = 1;        % Pacing freq, Hz
duration          = 10e3;     % [ms]
ISO_param         = 0;        % ISO 0.1 or 0.02   % [uM] - SET LIGAND CONCENTRATION HERE 0-0.1
beat_analysis     = 1;        % Output APD90,CaT Amp, RMP...etc

%% Run Single Cell

frequencies = freq;
[result] = A0_Obtain_intial_conditions(duration,frequencies,ISO_param,beat_analysis);


%% Save Biomarkers

% save population_600_biomarkers result

