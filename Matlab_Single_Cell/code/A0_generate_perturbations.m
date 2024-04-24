close all
clear
clc

%% Parameters

parameter_names = {'GNa' 'GNaL' 'GNaB' 'VNaK' 'Gto'...
    'GKur' 'GK2P' 'GKr' 'GKs' 'GK1' 'GKp' 'GKach' 'GKCa'...
    'GCaL' 'VNCX' 'ClCa' 'ClCa' 'Clb' 'Jrel' 'JSERCA' 'JLeak'...
    'Ca Buffer-Cleft' 'Ca Buffer-SR' 'GCab' ' GCap' 'Ca Buffer-Cyto'};
n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variation
variations = 10000; % number of trials

sigmaG = 0.1*ones(1,n_parameters); % standard deviation for parameters

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    scaling_2 = sigmaG(ii)*randn(1,variations) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end


aaa = mean(scaling_2);
aa2 = mean(scaling);
% all_parameters % size(all_parameters)
% columns: N parameters
% rows: N trials

%%
save parameter_matrix_600_trial all_parameters
