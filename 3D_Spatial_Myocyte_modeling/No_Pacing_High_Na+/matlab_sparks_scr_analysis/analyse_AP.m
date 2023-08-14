clc
clear
close all  

    
folder_name = {'AF.BCL.250.CKII_db.ISO.0', 'AF.BCL.250.CKII_db.ISO.1', 'AF.CL.250.CKII_inh.ISO.0', 'AF.CL.250.CKII_inh.ISO.1'}

% for id_folder = 1:length(folder_name)

%     f_n = folder_name{id_folder};

    % checker(f_n, t_range, x_select, y_select, NX, NY, z); 

    CKII_AP = load('../AF.BCL.250.CKII_db.ISO.0/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_AP(:,2), 1, 50);
% end

out_CKII.amplitude = amplitude;
out_CKII.duration = duration;
out_CKII.num_pks = num_pks;



CKII_ISO_AP = load('../AF.BCL.250.CKII_db.ISO.1/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_ISO_AP(:,2), 1, 50);
% end

out_CKII_ISO.amplitude = amplitude;
out_CKII_ISO.duration = duration;
out_CKII_ISO.num_pks = num_pks;




    CKII_inh_AP = load('../AF.CL.250.CKII_inh.ISO.0/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_inh_AP(:,2), 1, 50);
% end

out_CKII_inh.amplitude = amplitude;
out_CKII_inh.duration = duration;
out_CKII_inh.num_pks = num_pks;




CKII_inh_ISO_AP = load('../AF.CL.250.CKII_inh.ISO.1/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_inh_ISO_AP(:,2), 1, 50);
% end

out_CKII_inh_ISO.amplitude = amplitude;
out_CKII_inh_ISO.duration = duration;
out_CKII_inh_ISO.num_pks = num_pks;



figure(7); hold on
subplot(1,2,1);hold on
swarmchart(ones(1,out_CKII.num_pks), [out_CKII.amplitude], 36,[0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,out_CKII_ISO.num_pks)+1, [out_CKII_ISO.amplitude], 36,[0 0.4470 0.7410],'LineWidth',2)
title('num_peaks')
xlim([0.5,2.5])
% ylim([0,30])
xticks([1,2])
% yticks(0:10:30)
ylabel('DAD (mV)')
set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)
saveas(gcf,'F7','epsc')






    % CKII_AP = load('../AF.BCL.250.CKII_db.ISO.0/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_AP(:,3), 0.05, 50);
% end

Ca_out_CKII.amplitude = amplitude;
Ca_out_CKII.duration = duration;
Ca_out_CKII.num_pks = num_pks;



% CKII_ISO_AP = load('../AF.BCL.250.CKII_db.ISO.1/global_result/global_cai.txt');


[peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(CKII_ISO_AP(:,3), 0.05, 50);
% end

Ca_out_CKII_ISO.amplitude = amplitude;
Ca_out_CKII_ISO.duration = duration;
Ca_out_CKII_ISO.num_pks = num_pks;


Vm_ca_coupling_ISO = out_CKII_ISO.amplitude(4:end) ./ Ca_out_CKII_ISO.amplitude(4:end)
Vm_ca_coupling = out_CKII.amplitude(2:4) ./ Ca_out_CKII.amplitude(2:4)



   16.4301
   62.0100
   20.4084
   23.6047
   18.5109
   20.2029
   17.5921
   13.2663
   15.8587
   13.3917

     12.6580
   12.3492
    9.9077


%  figure, plot(reshape(mean(ci(1,:,:),2), 7501,[])*100-100); hold on;
% >> plot (CKII_ISO_AP(1:2:end,2))



figure(8); hold on
subplot(1,2,1);hold on

plot(CKII_inh_AP(:,2), 'color', [0 0.4470 0.7410],'LineWidth',2); 
plot(CKII_AP(:,2), 'color', [0.8500 0.3250 0.0980],'LineWidth',2); 

xlim([-100,5000])
ylim([-85,20])

ylabel('V_m (mV)')
set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)

subplot(1,2,2);hold on

plot(CKII_inh_ISO_AP(:,2), 'color', [0 0.4470 0.7410],'LineWidth',2); 
plot(CKII_ISO_AP(:,2), 'color', [0.8500 0.3250 0.0980],'LineWidth',2); 
xlim([-100,5000])
ylim([-85,20])

set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)
set(gcf,'Position',[108,250,1100,250])

saveas(gcf,'F8','epsc')
