
clear;



clear; close all;


res_db_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_norm_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
res_inh_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_norm_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
res_inh_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_norm_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
res_inh_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_norm_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
res_inh_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");





cell_num = length([res_db_1Hz.tAP]);
figure(1), hold on
title('tAP')

% swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP],18, [31, 119, 180]/250)
% swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP],18, [214, 39, 40]/255)

swarmchart(0+ones(1,cell_num)+0, [res_inh_1Hz.tAP],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(0+ones(1,cell_num)+1, [res_norm_1Hz.tAP],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(0+ones(1,cell_num)+2, [res_db_1Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(3.5+ones(1,cell_num)+0, [res_inh_2Hz.tAP],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(3.5+ones(1,cell_num)+1, [res_norm_2Hz.tAP],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(3.5+ones(1,cell_num)+2, [res_db_2Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(7+ones(1,cell_num)+0, [res_inh_3Hz.tAP],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(7+ones(1,cell_num)+1, [res_norm_3Hz.tAP],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(7+ones(1,cell_num)+2, [res_db_3Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(10.5+ones(1,cell_num)+0, [res_inh_4Hz.tAP],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(10.5+ones(1,cell_num)+1, [res_norm_4Hz.tAP],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(10.5+ones(1,cell_num)+2, [res_db_4Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\DeltaV_m (mV)')


cell_num = length([res_db_1Hz.tAP_num]);

figure(2), hold on

title('tAP number in each cell')

% swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP_num],36, [31, 119, 180]/250)
% swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP_num],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)

swarmchart(0+ones(1,cell_num)+0, [res_inh_1Hz.tAP_num],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(0+ones(1,cell_num)+1, [res_norm_1Hz.tAP_num],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(0+ones(1,cell_num)+2, [res_db_1Hz.tAP_num],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(3.5+ones(1,cell_num)+0, [res_inh_2Hz.tAP_num],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(3.5+ones(1,cell_num)+1, [res_norm_2Hz.tAP_num],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(3.5+ones(1,cell_num)+2, [res_db_2Hz.tAP_num],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(7+ones(1,cell_num)+0, [res_inh_3Hz.tAP_num],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(7+ones(1,cell_num)+1, [res_norm_3Hz.tAP_num],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(7+ones(1,cell_num)+2, [res_db_3Hz.tAP_num],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)



swarmchart(10.5+ones(1,cell_num)+0, [res_inh_4Hz.tAP_num],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(10.5+ones(1,cell_num)+1, [res_norm_4Hz.tAP_num],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
swarmchart(10.5+ones(1,cell_num)+2, [res_db_4Hz.tAP_num],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)


ylabel('tAP incidence (#)')





Xticklabels = {'1 Hz','2 Hz', '3 Hz', '4 Hz'}

for i = 1:2

	h = figure(i);

	box off
	xlim([0.5,14])

set(gca,'xtick',[2, 2+3.5, 2+3.5+3.5, 2+3.5+3.5+3.5],'tickdir', 'out');
set(gca,'xticklabel',Xticklabels,'fontsize',22,'fontname','arial', 'innerposition', [0.1300 0.1100 0.7750 0.8150])

% set(gca,'xticklabel',Xticklabels,'fontsize',18)
% set(gcf,'position',[100,100,width,height])
set(gcf,'position',[100,100,480*2,350])
set(gca,'linewidth',1.5)
	figname = sprintf('tAPs_F_%d.NoISO.svg', i);
saveas(gcf,figname,'svg')
end


function [outputs] = load_res(filename)

	load_res = analyze_single_file(filename);
	outputs = load_res.res;
end