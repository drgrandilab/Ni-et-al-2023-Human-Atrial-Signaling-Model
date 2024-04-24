clear all; close all; clc

%% Parameters

freq              = 1;        % Pacing freq, Hz
duration          = 10e3;    % Duration of simulation (ms)
plot_currents     = 1;        % Plot all currents
plot_basic        = 1;        % Plot Vm, CaT, [Na]i
ISO_param         = 0;        % ISO 0.1 or 0.02   % [uM] - SET LIGAND CONCENTRATION HERE 0-0.1
Int_cond          = 0;        % Saves each simulations final conditions;
T_Y_value         = 1;        % Save Y values and Time
output_currents   = 1;        % Save all Currents
outout_biomarkers = 1;        % Output APD90, CaT Amp, RMP...etc

%% Run Single Cell

for i = 1:numel(freq)
    frequencies = freq(i);
    [result, t, y, currents, txt] = NH_single_cell(duration,frequencies, ISO_param, plot_currents,...
        plot_basic, Int_cond, T_Y_value, output_currents, outout_biomarkers);

    results_array_full(i,:) = result;
    results_array = array2table(results_array_full, 'VariableNames', txt);

end

%% Save Biomarkers

y_and_t = [y t];
% save 0p5-4Hz_0p0_ISO_biomarkers result_array
save 1Hz_0p0_ISO_y_values_extra y_and_t
save 1Hz_0p0_ISO_currents_extra currents

