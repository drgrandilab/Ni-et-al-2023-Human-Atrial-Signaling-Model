function [output, t, y, total_store_values, txt] = NH_single_cell(duration, frequency,...
    ISO_param, Currents_record_val, plot_par, save_int_cond,  output_y_values, output_currents, beat_analysis)

%% Parameters for external modules

load y_final.mat
y0n = y_final';

AF_index = 0; % nSR = 0, AF = 1 (AF-Model is not finished)
freq = frequency;
cycleLength = 1e3/freq;
ISO = ISO_param;
Currents_record = Currents_record_val;
PLOT = plot_par;

%Prot_index
%protocol = 'pace_cc';
%protocol = 'pace_cc_erp';
%protocol = 'v_step';
%protocol = 'DAD';

prot_index = 1;
gender_param = zeros(1,9);

%% Establish and define globals

if Currents_record  == 1 || output_currents == 1

    global tStep tArray ICa_store Ito_store INa_store IK1_store
    global Jserca_store IKs_store Jleak_store Incx_store
    global INaK_store Ikr_store INabk_store
    global Lmyo_store Fmyo_store Vmax_store
    global IKp_store Iclca_store Iclbk_store Ipca_store Icabk_store Cai_store
    global I_app_store I_NaLstore I_kurstore I_k2pstore I_kistore I_kachstore
    global I_skstore I_ClCFTRstore I_pcastore I_J_SRCarelstore

    tStep = 1; tArray = zeros(1,1e6); ICa_store = zeros(1,1e6); Ito_store = zeros(1,1e6);
    INa_store = zeros(1,1e6); IK1_store = zeros(1,1e6);
    Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
    Jleak_store = zeros(1e6,2); Incx_store = zeros(1,1e6);
    INaK_store = zeros(1,1e6); Ikr_store = zeros(1,1e6); INabk_store = zeros(1,1e6);
    Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6); Vmax_store = zeros(1,1e6);
    IKp_store = zeros(1,1e6); Iclca_store = zeros(1,1e6); Iclbk_store = zeros(1,1e6);
    Ipca_store = zeros(1,1e6); Icabk_store = zeros(1,1e6); Cai_store = zeros(1,1e6); I_J_SRCarelstore = zeros(1,1e6);
    I_app_store = zeros(1,1e6); I_NaLstore = zeros(1,1e6); I_kurstore = zeros(1,1e6);
    I_k2pstore = zeros(1,1e6); I_kistore = zeros(1,1e6); I_kachstore = zeros(1,1e6);
    I_skstore = zeros(1,1e6); I_ClCFTRstore = zeros(1,1e6); I_pcastore = zeros(1,1e6);

end

%% Run single simulation
prot_rate = freq;
period = 1000/prot_rate;
num_beats = floor(duration/period);   % can add +1 here for extra beat
duration = (num_beats*period);

tic

par_SA = ones(1,25);
Voltage_step = 0;
p = [cycleLength, AF_index, prot_index, ISO, Currents_record, output_currents, par_SA];
tspan = [0 duration]; % [ms]
options = odeset('RelTol',1e-6,'MaxStep',1);
[t,y] = ode15s(@NH_ODE,tspan,y0n,options,p);
yfinal = y(end,:);

toc
%% Rename outputs
if Currents_record == 1 || output_currents == 1

    tA = tArray(1:tStep); dVm = Vmax_store(1:tStep);
    ICa = ICa_store(1:tStep); Ito = Ito_store(1:tStep);
    INa = INa_store(1:tStep); IK1 = IK1_store(1:tStep);
    IKs = IKs_store(1:tStep); IKr = Ikr_store(1:tStep);
    IKp = IKp_store(1:tStep);
    IClCa = Iclca_store(1:tStep); IClbk = Iclbk_store(1:tStep);
    INCX = Incx_store(1:tStep); INaK = INaK_store(1:tStep);
    INabk = INabk_store(1:tStep); IPMCA = Ipca_store(1:tStep);
    ICabk = Icabk_store(1:tStep); Cai_tA = Cai_store(1:tStep);
    Jserca = Jserca_store(1:tStep); Jleak = Jleak_store(1:tStep,:);
    I_app = I_app_store(1:tStep); I_NaL = I_NaLstore(1:tStep); I_Kur = I_kurstore(1:tStep);
    I_K2P = I_k2pstore(1:tStep); I_Ki = I_kistore(1:tStep); I_Kach = I_kachstore(1:tStep);
    I_SK = I_skstore(1:tStep); I_ClFTR = I_ClCFTRstore(1:tStep); I_pca = I_pcastore(1:tStep);
    I_J_SRCarel = I_J_SRCarelstore(1:tStep); Jleak = Jleak(:,1);
end

if output_currents == 1
    total_store_values = [dVm; ICa; Ito; INa; IKs; IKr; IKp; IClCa; IClbk; INCX; INaK;...
        INabk; IPMCA; ICabk; Cai_tA; Jserca; I_app; I_NaL; I_Kur; I_K2P;...
        I_Ki; I_Kach; I_SK; I_ClFTR; I_pca; IK1; I_J_SRCarel; tA]';

    total_store_values(:,29) = Jleak(:,1);
%     total_store_values(:,30) = Jleak(:,2);
else
    total_store_values = 0;
end
%% Current Plot
if Currents_record == 1
    figure(2)
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,2,1); plot(tA,INa); title('INa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,3); plot(tA,ICa); title('ICa');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,5); plot(tA,INCX); title('INCX');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,7); plot(tA,INaK); title('INaK');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,2); plot(tA,IKr); title('IKr');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,4); plot(tA,IKs); title('IKs');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,6); plot(tA,Ito); title('Ito');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlim([duration-period*2 duration])
    subplot(4,2,8); plot(tA,IK1); title('IK1');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(3)
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1); plot(tA,IKp); title('IKp');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2); plot(tA,IClCa); title('IClCa');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3); plot(tA,IClbk); title('IClbk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4); plot(tA,INabk); title('INabk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5); plot(tA,IPMCA); title('IPMCA');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6); plot(tA,ICabk); title('ICabk');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7); plot(tA,Cai_tA); title('Cai_tA');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8); plot(tA,Jserca); title('Jserca');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9); plot(tA,Jleak); title('Jleak');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

    figure(4)
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,1); plot(tA,I_app); title('Iapp');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,2); plot(tA,I_NaL); title('INaL');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,3); plot(tA,I_Kur); title('IKur');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,4); plot(tA,I_K2P); title('IK2P');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,5); plot(tA,I_Ki); title('IKi');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,6); plot(tA,I_Kach); title('IKach');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,7); plot(tA,I_SK); title('ISK');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,8); plot(tA,I_ClFTR); title('IClFTR');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(5,2,9); plot(tA,I_pca); title('Ipca');
    xlim([duration-period*2 duration])
    set(gca,'box','off','tickdir','out','fontsize',12)

end
%% Plot Vm, CaT, Na

if PLOT == 1
    figure(1)
    subplot(3,1,1)

    hold on,
    plot(t,y(:,39),'LineWidth', 1.5)
    ylabel('Em (mV)')
    xlim([(duration-period*1)-10 duration])
    fontsize(gca,20,"pixels")
    ax = gca;
    ax.TickDir = 'out';
    set(gca,'linewidth',2, 'FontWeight','bold')

    subplot(3,1,2)
    hold on,
    plot(t,y(:,38)*1000000,'LineWidth', 1.5)
    ylabel('CaT (nM)')
    xlim([(duration-period*1)-10 duration])
    fontsize(gca,20,"pixels")
    ax = gca;
    ax.TickDir = 'out';
    set(gca,'linewidth',2, 'FontWeight','bold')

    subplot(3,1,3)
    hold on,
    plot(t,y(:,34),'LineWidth', 1.5)
    ylabel('[Na]_i (mM)')
    xlabel('Time (ms)')
    xlim([(duration-period*1)-10 duration])
    fontsize(gca,20,"pixels")
    ax = gca;
    ax.TickDir = 'out';
    set(gca,'linewidth',2, 'FontWeight','bold')
end

if beat_analysis == 1
    time = t; % (ms)
    Vm = y(:,39); % (mV)
    Ca = y(:,38); % (mM)
    Na = y(:,34); % (mM)
    dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
    dVm = [dVm_array(1); dVm_array];
    AP_index = 2;
    period = 1000/freq;

    outputs = function_beat_analysis(time,Vm,Ca,Na,dVm,period,AP_index);
    dVm_max = outputs(1);
    Vm_max = outputs(2);
    RMP = outputs(3);
    AP_amp = outputs(4);
    APD90 = outputs(5);
    APD70 = outputs(6);
    APD50 = outputs(7);
    APD30 = outputs(8);
    Ca_max = outputs(9)*1000000;
    Ca_min = outputs(10)*1000000;
    CaT_amp = outputs(11)*1000000;
    CaT_rise = outputs(12);
    CaT_decay_50 = outputs(13);
    CaT_decay_63 = outputs(14);
    Na_min = outputs(15);

    output = [dVm_max Vm_max -RMP AP_amp APD90 APD70 APD50 APD30 Ca_max...
        Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min];
    txt = {'dVm_max', 'Vm_max', 'RMP', 'AP_amp', 'APD90', 'APD70', 'APD50', 'APD30', 'Ca_max',...
        'Ca_min', 'CaT_amp', 'CaT_rise', 'CaT_decay_50', 'CaT_decay_63', 'Na_min'};
else
    output = [];
    txt = [];
end


if output_y_values ~= 1
    t = 0;
    y = 0;
    
end

if save_int_cond == 1
save yfin_nSR_1Hz_0p0_ISO yfinal
end
