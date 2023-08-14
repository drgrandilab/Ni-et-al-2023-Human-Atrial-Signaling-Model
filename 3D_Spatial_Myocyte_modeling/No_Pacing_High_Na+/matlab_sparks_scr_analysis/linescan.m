clc
clear
close all  

% tub_trails          = {'T1' 'T5' 'T11' 'T15' 'T20'};
% rate_trails         = {'2hz' '3hz' '4hz' '5hz'};
% % folder_name       = {'1018_paper1_AP4' '1013_newAP_python_AP4' '1106_mn_AP4'};

% len_tub             = length(tub_trails);
% len_rate            = length(rate_trails);
% len_folder          = length(folder_name); 
    
t_range     = 0:2:15000;
x_select    = 10:3:46;
y_select    = 1:17; 
NX  = 55; 
NY  = 17; 
z   = 6; 
% for id_folder = 1:len_folder
%     f_n = folder_name{id_folder};
%             for id_tub = 1:len_tub
%                 tub = tub_trails{id_tub};
%                 for id_rate = 1:len_rate
%                     file_rate   = rate_trails{id_rate};
%                     checker(f_n, tub, file_rate, t_range, x_select, y_select, NX, NY, z); 
%                 end
%             end
% end
    

    
folder_name = {'AF.BCL.250.CKII_db.ISO.0', 'AF.BCL.250.CKII_db.ISO.1', 'AF.CL.250.CKII_inh.ISO.0', 'AF.CL.250.CKII_inh.ISO.1'}

for id_folder = 1:length(folder_name)

    f_n = folder_name{id_folder};

    checker(f_n, t_range, x_select, y_select, NX, NY, z); 

end




out_CKII = analyze_linescanes(folder_name{1});
out_CKII_ISO = analyze_linescanes(folder_name{2});
out_CKII_inh = analyze_linescanes(folder_name{3});
out_CKII_inh_ISO = analyze_linescanes(folder_name{4});


figure(3); hold on
subplot(1,2,1);hold on
swarmchart(ones(1,13), [out_CKII.num_peaks], 36,[0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh.num_peaks], 36,[0 0.4470 0.7410],'LineWidth',2)
title('num_peaks')
xlim([0.5,2.5])
ylim([0,30])
xticks([1,2])
yticks(0:10:30)
ylabel('SCR incidence #')
set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)

subplot(1,2,2);hold on
swarmchart(ones(1,13), [out_CKII_ISO.num_peaks],36, [0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh_ISO.num_peaks],36, [0 0.4470 0.7410],'LineWidth',2)
title('num_peaks')
xlim([0.5,2.5])
ylim([0,30])
xticks([1,2])
yticks(0:10:30)
set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)

set(gcf,'Position',[108,236,500,326])

saveas(gcf,'F3','epsc')


figure(4)
subplot(1,2,1);hold on
swarmchart(ones(1,13), [out_CKII.mean_amplitude], 36,[0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh.mean_amplitude], 36,[0 0.4470 0.7410],'LineWidth',2)
title('mean_amplitude')
xlim([0.5,2.5])
ylim([0,0.4])
xticks([1,2])
yticks(0:0.1:0.4)
ylabel('Mean \Delta [Ca^{2+}]_i (\muM)')

set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)


subplot(1,2,2);hold on
swarmchart(ones(1,13), [out_CKII_ISO.mean_amplitude],36, [0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh_ISO.mean_amplitude],36, [0 0.4470 0.7410],'LineWidth',2)
title('mean_amplitude')
xlim([0.5,2.5])
ylim([0,0.4])
xticks([1,2])
yticks(0:0.1:0.4)


set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)
set(gcf,'Position',[108,236,500,326])
saveas(gcf,'F4','epsc')




figure(5)
subplot(1,2,1);hold on
swarmchart(ones(1,13), [out_CKII.mean_duration]*(t_range(2)-t_range(1)), 36,[0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh.mean_duration]*(t_range(2)-t_range(1)), 36,[0 0.4470 0.7410],'LineWidth',2)
title('meanduration')
xlim([0.5,2.5])
% ylim([0,0.4])
ylim([0,90])

xticks([1,2])
% yticks(0:0.1:0.4)
yticks(0:30:90)

ylabel('Mean duration (ms)')

set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)


subplot(1,2,2);hold on
swarmchart(ones(1,13), [out_CKII_ISO.mean_duration]*(t_range(2)-t_range(1)),36, [0.8500 0.3250 0.0980],'LineWidth',2)
swarmchart(ones(1,13)+1, [out_CKII_inh_ISO.mean_duration]*(t_range(2)-t_range(1)),36, [0 0.4470 0.7410],'LineWidth',2)
title('mean_duration')
xlim([0.5,2.5])
ylim([0,90])
xticks([1,2])
yticks(0:30:90)


set(gca,"FontSize",18,'fontname','arial')
set(gca, 'LineWidth',1.5)
set(gcf,'Position',[108,236,500,326])
saveas(gcf,'F5','epsc')



[h1,p1] = ttest2(out_CKII.num_peaks, out_CKII_inh.num_peaks)

[h2,p2] = ttest2(out_CKII_ISO.num_peaks, out_CKII_inh_ISO.num_peaks)

[h3,p3] = ttest2(out_CKII_ISO.mean_amplitude, out_CKII_inh_ISO.mean_amplitude)

[h4,p4] = ttest2(out_CKII_ISO.mean_duration, out_CKII_inh_ISO.mean_duration)

% figure(4)
% subplot(1,2,2);hold on
% swarmchart(ones(1,13), [out_CKII.num_peaks])

% swarmchart(ones(1,13)+2.5, [out_CKII_ISO.num_peaks])
% swarmchart(ones(1,13)+3.5, [out_CKII_inh_ISO.num_peaks])
% xlim([0.5,5.5])
% xticks([1,2,3.5,4.5])
% title('num peaks')


% figure(5); hold on
% subplot(1,2,1);hold on

% % duration should be scaled by t_range accordingly

% swarmchart(ones(1,13), [out_CKII.mean_duration]*(t_range(2)-t_range(1)))
% swarmchart(ones(1,13)+1, [out_CKII_inh.mean_duration]*(t_range(2)-t_range(1)))
% swarmchart(ones(1,13)+2.5, [out_CKII_ISO.mean_duration]*(t_range(2)-t_range(1)))
% swarmchart(ones(1,13)+3.5, [out_CKII_inh_ISO.mean_duration]*(t_range(2)-t_range(1)))
% xlim([0.5,5.5])
% xticks([1,2,3.5,4.5])
% title('mean duration')


function  output = analyze_linescanes(folder_name)

load([folder_name '_.mat'])




figure(1)
imagesc(reshape(ci(7,:,:),17,[]))
colormap(hot)

res = mean(ci,2); % get mean of the 17 CRUs for each x linescan
res = reshape(res, 13,[]);  % reshape to 2D matrix;
figure(2)
plot(res');


output.mean_amplitude = [];
output.num_peaks = [];
output.mean_duration = [];

for i = 1:13
    [peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(res(i,:), 0.05,4);

    output.mean_amplitude(end+1) = mean(amplitude);
    output.num_peaks(end+1) = num_pks;
    output.mean_duration(end+1) = mean(duration);
end


end
