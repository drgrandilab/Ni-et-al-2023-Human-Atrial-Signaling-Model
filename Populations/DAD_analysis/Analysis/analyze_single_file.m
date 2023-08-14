% % analyze_single_file: function description
function [outputs] = analyze_single_file(filename)
% outputs = ;
data = load(filename);

% threshold = 20; % Vogit et al., 2012
threshold = 50;

tAP = zeros(1, 600);
DAD = zeros(1, 600);
Ca_DAD = zeros(1, 600);
Ca_tAP = zeros(1, 600);
tAP_num = zeros(1, 600);
DAD_num = zeros(1, 600);
tAP_Ca_ratio = zeros(1, 600);
DAD_Ca_ratio = zeros(1, 600);


outputs = {};
for i = 1 : 600

    d = data.data2{i} .DAD_ind;
    res(i) = analyze_single_cell(d);
    out = analyze_single_cell(d);
    tAP(i) = out.tAP;
    Ca_tAP(i) = out.Ca_tAP;
    tAP_num(i) = out.tAP_num;
    tAP_Ca_ratio(i) = out.tAP_Ca_ratio;
    DAD(i) = out.DAD;
    Ca_DAD(i) = out.Ca_DAD;
    DAD_num(i) = out.DAD_num;
    DAD_Ca_ratio(i) = out.DAD_Ca_ratio;
end



outputs.res = res;
outputs. tAP = tAP;
outputs. Ca_tAP = Ca_tAP;
outputs. tAP_num = tAP_num;
outputs. tAP_Ca_ratio = tAP_Ca_ratio;
outputs. DAD = DAD;
outputs. Ca_DAD = Ca_DAD;
outputs. DAD_num = DAD_num;
outputs. DAD_Ca_ratio = DAD_Ca_ratio;

end



function [output] = analyze_single_cell(d)
threshold = 50; 
threshold_DAD = 1;
num  = length(d);

output.tAP = nan;
output.Ca_tAP = nan;
output.tAP_num = nan;
output.DAD = nan;
output.Ca_DAD = nan;
output.DAD_num = nan;
output.tAP_Ca_ratio = nan;
output.DAD_Ca_ratio = nan;
output.Arrhy = 0;
output.latency_1st = nan;
output.latency_DAD_mean = nan;
output.mean_Ca_at_v_threshold = nan;
output.dCa_dt_at_v_threshold = nan;

output.last_tAP_time = nan;

if length(d) >= 1

    tmp_tAP = [];
    tmp_peakCa_tAP = [];
    tmp_Thres_Ca_tAP = [];
    tmp_dCa_dt_at_v_threshold = [];
    tmp_DAD = [];
    tmp_peakCa_DAD = [];
    time = [];
    for jj = 1 : num

        if (d{jj} .AMP) > threshold
            tmp_tAP(end + 1) = d{jj} .AMP;
            tmp_peakCa_tAP(end + 1) = d{jj}.Ca_AMP_peak;
            tmp_Thres_Ca_tAP(end + 1) = d{jj}.Ca_at_v_threshold;
            tmp_dCa_dt_at_v_threshold(end + 1) = d{jj}.dCa_dt_at_v_threshold;
            time(end+1) = d{jj}.time;
            output.last_tAP_time = d{jj}.time;

        elseif (d{jj} .AMP) >= threshold_DAD

            tmp_DAD(end + 1) = d{jj} .AMP;
            % disp('aaa')
            tmp_peakCa_DAD(end + 1) = d{jj} .Ca_AMP_peak;
            time(end+1) = d{jj}.time;
        end
    end

    if(length(time)>=1)
    output.latency_1st = time(1);
    end 
    if length(time) >= 2
        output.latency_DAD_mean = (time(end) - time(1)) / (num -1);
    end
    output.tAP = mean(rmmissing(tmp_tAP));
    output.Ca_tAP = mean(rmmissing(tmp_peakCa_tAP));
    output.tAP_Ca_ratio = mean(rmmissing(tmp_tAP./ tmp_peakCa_tAP));
    output.tAP_num = length(tmp_tAP);

    output.mean_Ca_at_v_threshold = mean(rmmissing(tmp_Thres_Ca_tAP));
    output.dCa_dt_at_v_threshold = mean(rmmissing(tmp_dCa_dt_at_v_threshold));

    output.DAD = mean(rmmissing(tmp_DAD));
    output.Ca_DAD = mean(rmmissing(tmp_peakCa_DAD));
    output.DAD_Ca_ratio = mean(rmmissing(tmp_DAD./ tmp_peakCa_DAD));
    output.DAD_num = length(tmp_DAD);

    if(output.DAD_num >0 || output.tAP_num>0)
        output.Arrhy = 1;
    end
    if(output.DAD_num == 0)
        output.DAD_num = nan;
    end
    if(output.tAP_num == 0)
        output.tAP_num = nan;
    end
end

end
