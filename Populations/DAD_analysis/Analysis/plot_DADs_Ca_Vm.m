
clear;



clear; close all;


% res_db_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
% res_norm_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
% res_inh_1Hz = load_res("Summarized_data/BCL.1000.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

% res_db_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
% res_norm_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
% res_inh_2Hz = load_res("Summarized_data/BCL.500.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

res_db_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
res_norm_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
res_inh_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

% res_db_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
% res_norm_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
% res_inh_4Hz = load_res("Summarized_data/BCL.250.ISO.0.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");







figure(20), hold on;




scatter([res_db_3Hz.Ca_DAD], [res_db_3Hz.DAD],'ro')
alpha(.5)

scatter([res_db_3Hz.Ca_tAP], [res_db_3Hz.tAP],'ro')
alpha(.5)



scatter([res_norm_3Hz.Ca_DAD], [res_norm_3Hz.DAD],'ko')
alpha(.5)

scatter([res_norm_3Hz.Ca_tAP], [res_norm_3Hz.tAP],'ko')
alpha(.5)


scatter([res_inh_3Hz.Ca_DAD], [res_inh_3Hz.DAD],'bo')
alpha(.5)

scatter([res_inh_3Hz.Ca_tAP], [res_inh_3Hz.tAP],'bo')
alpha(.5)



ISO_res_db_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
ISO_res_norm_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.1.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
ISO_res_inh_3Hz = load_res("Summarized_data/BCL.333.333.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");

% ISO_res_db_4Hz = load_res("Summarized_data/BCL.250.ISO.0.1.CaMKII_inhb.0.CaMKII_db.1_DAD_only_60s.mat");
% ISO_res_norm_4Hz = load_res("Summarized_data/BCL.250.ISO.0.1.CaMKII_inhb.0.CaMKII_db.0_DAD_only_60s.mat");
% ISO_res_inh_4Hz = load_res("Summarized_data/BCL.250.ISO.0.1.CaMKII_inhb.1.CaMKII_db.0_DAD_only_60s.mat");
xlabel('\DeltaCa (\muM)')
ylabel('\DeltaVm (mV)')

figure(21), hold on;


p1=scatter([ISO_res_db_3Hz.Ca_DAD], [ISO_res_db_3Hz.DAD],'ro')
alpha(.5)
p2=scatter([ISO_res_db_3Hz.Ca_tAP], [ISO_res_db_3Hz.tAP],'ro')
alpha(.5)
p3=scatter([ISO_res_norm_3Hz.Ca_DAD], [ISO_res_norm_3Hz.DAD],'ko')%,'MarkerEdgeColor',[0.5 .5 .5])
alpha(.5)
p4=scatter([ISO_res_norm_3Hz.Ca_tAP], [ISO_res_norm_3Hz.tAP],'ko')%,'MarkerEdgeColor',[0.5 .5 .5])
alpha(.5)


scatter([ISO_res_inh_3Hz.Ca_DAD], [ISO_res_inh_3Hz.DAD],'bo')
alpha(.5)

scatter([ISO_res_inh_3Hz.Ca_tAP], [ISO_res_inh_3Hz.tAP],'bo')
alpha(.5)


alpha(.5)
xlabel('\DeltaCa (\muM)')
ylabel('\DeltaVm (mV)')
 % p1.Color(4) = 0.25;
 % p2.Color(4) = 0.75;
 % p3.Color(4) = 0.75;
 % p4.Color(4) = 0.75;
% % ylabel('\DeltaCa (\muM)')
% ylabel('\DeltaVm/\DeltaCa (mV/\muM)')



% Xticklabels = {'1 Hz','2 Hz', '3 Hz', '4 Hz'}

for i = 20:21

	h = figure(i);

	box off
	% xlim([0.5,14])

% set(gca,'xtick',[2, 2+3.5, 2+3.5+3.5, 2+3.5+3.5+3.5]);
set(gca, 'fontsize',22,'TickDir','out', 'fontname','arial', 'innerposition', [0.1300 0.1100 0.7750 0.8150])

set(gca,'fontsize',22)
xlim([-0.05,1.8])
ylim([-5,130])
% set(gcf,'position',[100,100,width,height])
set(gcf,'position',[100,100,480*2/1.2,480*1.7/1.2])
set(gca,'linewidth',1.5)
	figname = sprintf('R2_Vm_Ca_%d.NoISO.svg', i);
saveas(gcf,figname,'svg')
end



for i = 20:21

	h = figure(i);

	box off
	% xlim([0.5,14])

% set(gca,'xtick',[2, 2+3.5, 2+3.5+3.5, 2+3.5+3.5+3.5]);
set(gca, 'fontsize',22,'TickDir','out', 'fontname','arial', 'innerposition', [0.1300 0.1100 0.7750 0.8150])

set(gca,'fontsize',22)
xlim([-0.02,0.15])
ylim([-2,20])
% set(gcf,'position',[100,100,width,height])
set(gcf,'position',[100,100,480*2/1.2/1.9,480*1.7/1.2/1.9])
set(gca,'linewidth',1.5)
	figname = sprintf('R2_Vm_Ca_%d.NoISO.zoom.svg', i);
saveas(gcf,figname,'svg')
end



% cell_num = length([res_db_1Hz.DAD]);
% figure(1), hold on
% title('DAD')

% % swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP],18, [31, 119, 180]/250)
% % swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP],18, [214, 39, 40]/255)

% swarmchart(0+ones(1,cell_num)+0, [res_inh_1Hz.DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+1, [res_norm_1Hz.DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+2, [res_db_1Hz.DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(3.5+ones(1,cell_num)+0, [res_inh_2Hz.DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+1, [res_norm_2Hz.DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+2, [res_db_2Hz.DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(7+ones(1,cell_num)+0, [res_inh_3Hz.DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+1, [res_norm_3Hz.DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+2, [res_db_3Hz.DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(10.5+ones(1,cell_num)+0, [res_inh_4Hz.DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+1, [res_norm_4Hz.DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+2, [res_db_4Hz.DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% ylabel('\DeltaV_m (mV)')


% cell_num = length([res_db_1Hz.Ca_DAD]);

% figure(2), hold on
% title('Ca_DAD')

% % swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP],36, [31, 119, 180]/250)
% % swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)

% swarmchart(0+ones(1,cell_num)+0, [res_inh_1Hz.Ca_DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+1, [res_norm_1Hz.Ca_DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+2, [res_db_1Hz.Ca_DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(3.5+ones(1,cell_num)+0, [res_inh_2Hz.Ca_DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+1, [res_norm_2Hz.Ca_DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+2, [res_db_2Hz.Ca_DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(7+ones(1,cell_num)+0, [res_inh_3Hz.Ca_DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+1, [res_norm_3Hz.Ca_DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+2, [res_db_3Hz.Ca_DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(10.5+ones(1,cell_num)+0, [res_inh_4Hz.Ca_DAD],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+1, [res_norm_4Hz.Ca_DAD],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+2, [res_db_4Hz.Ca_DAD],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)

% ylabel('\DeltaCa (\muM)')



% figure(3), hold on
% title('\DeltaV_m/\DeltaCa^{2+} (for DADs)')

% % swarmchart(ones(1,cell_num), [res_inh_05Hz.tAP],36, [31, 119, 180]/250)
% % swarmchart(ones(1,cell_num)+1, [res_db_05Hz.tAP],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)

% swarmchart(0+ones(1,cell_num)+0, [res_inh_1Hz.DAD_Ca_ratio],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+1, [res_norm_1Hz.DAD_Ca_ratio],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(0+ones(1,cell_num)+2, [res_db_1Hz.DAD_Ca_ratio],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(3.5+ones(1,cell_num)+0, [res_inh_2Hz.DAD_Ca_ratio],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+1, [res_norm_2Hz.DAD_Ca_ratio],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(3.5+ones(1,cell_num)+2, [res_db_2Hz.DAD_Ca_ratio],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(7+ones(1,cell_num)+0, [res_inh_3Hz.DAD_Ca_ratio],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+1, [res_norm_3Hz.DAD_Ca_ratio],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(7+ones(1,cell_num)+2, [res_db_3Hz.DAD_Ca_ratio],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)



% swarmchart(10.5+ones(1,cell_num)+0, [res_inh_4Hz.DAD_Ca_ratio],36, [31, 119, 180]/250,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+1, [res_norm_4Hz.DAD_Ca_ratio],36, [127,127,127]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)
% swarmchart(10.5+ones(1,cell_num)+2, [res_db_4Hz.DAD_Ca_ratio],36, [214, 39, 40]/255,'XJitterWidth',0.7,'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9)

% % ylabel('\DeltaCa (\muM)')
% ylabel('\DeltaVm/\DeltaCa (mV/\muM)')



% Xticklabels = {'1 Hz','2 Hz', '3 Hz', '4 Hz'}

% for i = 1:3

% 	h = figure(i);

% 	box off
% 	xlim([0.5,14])

% set(gca,'xtick',[2, 2+3.5, 2+3.5+3.5, 2+3.5+3.5+3.5]);
% set(gca,'xticklabel',Xticklabels,'fontsize',22,'fontname','arial', 'innerposition', [0.1300 0.1100 0.7750 0.8150])

% % set(gca,'xticklabel',Xticklabels,'fontsize',18)
% % set(gcf,'position',[100,100,width,height])
% set(gcf,'position',[100,100,480*2,350])
% set(gca,'linewidth',1.5)
% 	figname = sprintf('DADs_F_%d.NoISO.svg', i);
% saveas(gcf,figname,'svg')
% end


function [outputs] = load_res(filename)

	load_res = analyze_single_file(filename);
	outputs = load_res.res;
end