function outputs = function_beat_analysis(time,Vm,Ca,Na,dVm,period,AP_index)
% Beat analysis first or last AP (index = 1 or 2)

%% Check
if AP_index == 1,
    t1 = 0; t2 = t1+period-1.5; t3 = t1+2*period-1.5;
else
    t1 = time(end)-2*period; t2 = t1+period; t3 = t1+2*period;
end

% Alternans
% t1_roi = find(time>t1); t1_index = t1_roi(1); Vm1 = Vm(t1_index);
% t2_roi = find(time>t2); t2_index = t2_roi(1); Vm2 = Vm(t2_index);
% t3_roi = find(time>t3); t3_index = t3_roi(1); Vm3 = Vm(t3_index);
% Vmin1 = min(Vm(t1_index:t2_index));
% Vmin2 = min(Vm(t2_index:t3_index));
%[abs(Vm1-Vm2) abs(Vm2-Vm3) abs(Vm1-Vm3) abs(Vm2-Vmin1) abs(Vm2-Vmin2)]
% delta = [abs(Vm1-Vm2) abs(Vm2-Vm3) abs(Vm1-Vm3) abs(Vm2-Vmin1) abs(Vm2-Vmin2)];

% Triggered APs
% t4_roi = find(time>t2+15); t4_index = t4_roi(1);
% taps = find(dVm(t4_index:t3_index)>25);
%figure, plot(time(t4_index:t3_index),dVm(t4_index:t3_index))

% if max(delta) > 5 %| length(taps) > 0,
%     outputs = zeros(1,15);
% else

    %% ROI
    if AP_index == 1,
        t_in = 0; t_fin = t_in+period-5;
    else
        t_in = time(end)-period; t_fin = t_in+period-5;
    end

    t_in_roi = find(time>t_in); t_in_index = t_in_roi(1);%-1;
    t_fin_roi = find(time>t_fin); t_fin_index = t_fin_roi(1);

    %% Em
    [dVm_max, index] = max(dVm(t_in_index:t_fin_index)); index = index-1; % max slope

    [Vm_max, index_max] = max(Vm(t_in_index:t_fin_index)); index_max = index_max-1; % peak
    Vm_min = min(Vm(t_in_index:t_fin_index)); % Em resting
    AP_amp = Vm_max-Vm_min;

    APfraction = 0.90;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD90 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);
    
    P20_APD90 = 0.20*APD90; P20_APD90 = P20_APD90 + time(t_in_index+index);  Index_P20 = find(time>P20_APD90); Index_P20_first = Index_P20(1);
    Second_value = P20_APD90 + 5; Index_20_second = find(time>Second_value); Second_index = Index_20_second(1);
    Vm_in_window = Vm(Index_P20_first:Second_index); Vm_in_window_mean = mean(Vm_in_window); Mean_closest = find(Vm_in_window>Vm_in_window_mean); %ms %https://www.sciencedirect.com/science/article/pii/S0022282823000068?via%3Dihub
    VPLT = Vm_in_window(Mean_closest(end)); Index_mean_value = find(Vm == VPLT); Time_VPLT = time(Index_mean_value);


    start_point = time(t_in_index+index); index_val_1 = find(time == start_point); Vm_val_1 = Vm(index_val_1);
    last_point = time(t_in_index+index_max+APD_roi(1)); index_val_2 = find(time == last_point); Vm_val_2 = Vm(index_val_2);


%     figure(1)
%     subplot(3,1,1)
%     hold on, plot (start_point, Vm_val_1, 'bo'); plot (last_point, Vm_val_2, 'bo');
%     plot(Time_VPLT, VPLT, 'ko');


    APfraction = 0.7;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD70 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.5;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD50 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    start_point_50 = time(t_in_index+index); index_val_1_50 = find(time == start_point_50); Vm_val_1_50 = Vm(index_val_1_50);
    last_point_50 = time(t_in_index+index_max+APD_roi(1)); index_val_2_50 = find(time == last_point_50); Vm_val_2_50 = Vm(index_val_2_50);


%     figure(1)
%     subplot(3,1,1)
%     hold on, plot (start_point_50, Vm_val_1_50, 'ro'); plot (last_point_50, Vm_val_2_50, 'ro'); 

    APfraction = 0.3;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD30 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    APfraction = 0.2;
    APD_roi = find(Vm(t_in_index+index_max:t_fin_index)<Vm_max-APfraction*AP_amp)-1;
    APD20 = time(t_in_index+index_max+APD_roi(1))-time(t_in_index+index);

    % Ca
    [Ca_max, index_ca] = max(Ca(t_in_index:t_fin_index)); store_index = index_ca; index_ca = index_ca-1; % peak CaT
    ca_max_index = find(Ca == Ca_max); time_max_ca = time(ca_max_index);

    Ca_min = min(Ca(t_in_index:t_fin_index)); % diast Ca
    CaT_amp = Ca_max-Ca_min;
    
%     start_point_ca = min(Ca(t_in_index:t_fin_index)); index_val_1_ca = find(Ca == start_point_ca); time_val_1_ca = time(index_val_1_ca);
%     last_point_ca = max(Ca(t_in_index:t_fin_index)); index_val_2_ca = find(Ca == last_point_ca); time_val_2_ca = time(index_val_2_ca);
% 
%     figure(1)
%     subplot(3,1,2),
%     hold on, plot (time_val_2_ca, Ca_max ,'bo'); plot (time_val_1_ca, Ca_min ,'bo')

    CaT_rise = time(t_in_index+index_ca)-time(t_in_index);

    Cafraction = 0.5;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_50 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);

    Cafraction = 0.632;
    Ca_roi = find(Ca(t_in_index+index_ca:t_fin_index)<Ca_max-Cafraction*CaT_amp)-1;
    CaT_decay_63 = time(t_in_index+index_ca+Ca_roi(1))-time(t_in_index+index_ca);

    %% Na
    [Na_min, index_na] = min(Na(t_in_index:t_fin_index)); index_na = index_na-1; % diast Na


    %% Collect all outputs
    outputs = [dVm_max Vm_max -Vm_min AP_amp APD90 APD70 APD50 APD30...
        Ca_max Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min];

end