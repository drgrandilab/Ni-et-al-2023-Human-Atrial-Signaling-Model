
clear; close all;

res_db_05Hz = load_res("Summarized_data/BCL.2000.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_inh_05Hz = load_res("Summarized_data/BCL.2000.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_inh_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_2Hz = load_res("Summarized_data/BCL.500.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_inh_2Hz = load_res("Summarized_data/BCL.500.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_inh_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_4Hz = load_res("Summarized_data/BCL.250.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_inh_4Hz = load_res("Summarized_data/BCL.250.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");





% Xticklabels = {'CaMKII\n inhb.', 'Normal', 'CaMKII\nX2','CaMKII\n inhb.', 'Normal', 'CaMKII\nX2'};
Xticklabels = {[]};


% figure(10)
% res_DAD_all = [res_inh.tAP_num; res.tAP_num; res_db.tAP_num; res_ISO_inh.tAP_num; res_ISO.tAP_num; res_ISO_db.tAP_num];
% [p,tbl,stats] = kruskalwallis(res_DAD_all',[],'off');
% c = multcompare(stats);
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])



% figure(11)
% res_DAD_all = [res_inh.DAD; res.DAD; res_db.DAD; res_ISO_inh.DAD; res_ISO.DAD; res_ISO_db.DAD];
% [p,tbl,stats] = kruskalwallis(res_DAD_all',[],'off');
% c = multcompare(stats);
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])



cell_num = length([res_db_05Hz.tAP_num]);
figure(1), hold on
title('tAPs')

% swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP],36, [31, 119, 180]/250)
% swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP],36, [214, 39, 40]/255)

swarmchart(2+ones(1,cell_num), [res_inh_1Hz.tAP],36, [31, 119, 180]/250)
swarmchart(2+ones(1,cell_num)+1, [res_db_1Hz.tAP],36, [214, 39, 40]/255)

swarmchart(4+ones(1,cell_num), [res_inh_2Hz.tAP],36, [31, 119, 180]/250)
swarmchart(4+ones(1,cell_num)+1, [res_db_2Hz.tAP],36, [214, 39, 40]/255)

swarmchart(6+ones(1,cell_num), [res_inh_3Hz.tAP],36, [31, 119, 180]/250)
swarmchart(6+ones(1,cell_num)+1, [res_db_3Hz.tAP],36, [214, 39, 40]/255)
swarmchart(8+ones(1,cell_num), [res_inh_4Hz.tAP],36, [31, 119, 180]/250)
swarmchart(8+ones(1,cell_num)+1, [res_db_4Hz.tAP],36, [214, 39, 40]/255)



% swarmchart(ones(1,cell_num), [res.tAP])

ylabel('\Delta Vm (mV)')

set(gca,'xtick',1:10);
% set(gca,'xticklabel',Xticklabels,'fontsize',12)
xtickangle(60)



figure(2), hold on
% title('\DeltaV_m/\DeltaCa^{2+} (for DADs)')


% title('tAPs')
title('tAP number in each cell')

% swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP_num],36, [31, 119, 180]/250)
% swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP_num],36, [214, 39, 40]/255)

swarmchart(2+ones(1,cell_num), [res_inh_1Hz.tAP_num],36, [31, 119, 180]/250)
swarmchart(2+ones(1,cell_num)+1, [res_db_1Hz.tAP_num],36, [214, 39, 40]/255)

swarmchart(4+ones(1,cell_num), [res_inh_2Hz.tAP_num],36, [31, 119, 180]/250)
swarmchart(4+ones(1,cell_num)+1, [res_db_2Hz.tAP_num],36, [214, 39, 40]/255)

swarmchart(6+ones(1,cell_num), [res_inh_3Hz.tAP_num],36, [31, 119, 180]/250)
swarmchart(6+ones(1,cell_num)+1, [res_db_3Hz.tAP_num],36, [214, 39, 40]/255)

swarmchart(8+ones(1,cell_num), [res_inh_4Hz.tAP_num],36, [31, 119, 180]/250)
swarmchart(8+ones(1,cell_num)+1, [res_db_4Hz.tAP_num],36, [214, 39, 40]/255)



% swarmchart(ones(1,cell_num), [res.tAP])
set(gca,'xtick',1:10);
% set(gca,'xticklabel',Xticklabels,'fontsize',12)
xtickangle(60)
% ylabel('\DeltaVm/\DeltaCa (mV/\muM)')
ylabel('tAP incidence (#)')




for i = 1:2

	h = figure(i);

	box off
	xlim([0.5+2,10.5])
set(gca,'xticklabel',Xticklabels,'fontsize',18,'fontname','arial', 'innerposition', [0.1300 0.1100 0.7750 0.8150])

% set(gca,'xticklabel',Xticklabels,'fontsize',18)
% set(gcf,'position',[100,100,width,height])

set(gcf,'position',[100,100,480*10/6,350])
set(gca,'linewidth',1.5)
	figname = sprintf('F_rate_dep_tAP_%d.eps', i);
saveas(gcf,figname,'epsc')
end

% 
% figure(9), hold on
% title('dCa_dt_at_v_threshold (for tAPs)')
% % title('tAPs')
% 
% swarmchart(ones(1,cell_num), [res_inh.dCa_dt_at_v_threshold])
% swarmchart(ones(1,cell_num)+1, [res.dCa_dt_at_v_threshold])
% swarmchart(ones(1,cell_num)+2, [res_db.dCa_dt_at_v_threshold])
% swarmchart(ones(1,cell_num)+3, [res_ISO_inh.dCa_dt_at_v_threshold])
% swarmchart(ones(1,cell_num)+4, [res_ISO.dCa_dt_at_v_threshold])
% swarmchart(ones(1,cell_num)+5, [res_ISO_db.dCa_dt_at_v_threshold])
% % swarmchart(ones(1,cell_num), [res.tAP])
% set(gca,'xtick',1:6);
% set(gca,'xticklabel',Xticklabels,'fontsize',12)
% xtickangle(60)
% ylabel('dCa_dt_at_v_threshold (uM/ms)')

% figure(3), hold on
% title('\DeltaVm/\DeltaCa2+')
% scatter([res_inh.Ca_DAD], [res_inh.DAD])
% scatter([res.Ca_DAD], [res.DAD])
% scatter([res_db.Ca_DAD], [res_db.DAD])

% scatter([res_ISO_inh.Ca_DAD], [res_ISO_inh.DAD])
% scatter([res_ISO.Ca_DAD], [res_ISO.DAD])
% scatter([res_ISO_db.Ca_DAD], [res_ISO_db.DAD])
% xlabel('Ca2+ (uM)')
% ylabel('DAD amplitude (mV)')

%% get_res: function description

function [outputs] = load_res(filename)

	load_res = analyze_single_file(filename);
	outputs = load_res.res;
end