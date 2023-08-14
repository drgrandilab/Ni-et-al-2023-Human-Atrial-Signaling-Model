function [outputs] = process_folder(folder)
% clear;clc;
data1 = {};
data2 = {};

% folder='../ISO';
stimulation_start_time = 0;

% stimulation_train = [10000:500:290000,301450];
% 2 Hz
% stimulation_train = [0:1000:600000,620990]; % full sim
% stimulation_train = [580000:1000:600000,620990]; % SS sim
% stimulation_train = [500000:1000:600000,620990]; % SS new


stimulation_train = [600000:500:600500,620990]; % 2 Hz

% stimulation_train = [280000:1000:289500,301450];
% stimulation_train = [270000:1000:289000,302000];

plot_figure = 0;

parfor id = 0:599
    
    file = sprintf('%s/AP.BCL.1000.ID.%d', folder, id);
    disp(file)
    data = load(file);
    data(:,1) = data(:,1)- stimulation_start_time;
    
    time1 = data(:,1);
    
    start_1 = find(time1 >= 0e3);   % starting from 0 sec after pacing.
    start_1 = start_1(1);
    % end_1 = find (time1 > 128e3-5);
    end_1 = length(time1);%end_1(1);
    
    time = data(start_1:end_1,1);
%     time = time - time(1);  # no need anymore with the pacing train
%     infor
     Vm = data(start_1:end_1,2);
    
    CaSR = data(start_1:end_1,5);
    
    Cai = data(start_1:end_1, 8);
    
    Nai = data(start_1:end_1,10);  % 33 in the full dataset, 10 in reduced dataset 
    % dV = data(:,3);
    
    dV=zeros(length(Vm),1);
    num = length(Vm);
    dV(1:num-1) = (Vm(2:num) - Vm(1:num-1)) ./((time(2:num) - time(1:num-1))) ;
    
    dV(end) = dV(end-1);
    
%     data1{id+1} =   function_beat_analysis_EAD(time,Vm,Cai,CaSR,Nai,dV,1000,1,stimulation_train);
    
    % data2{id+1} = function_beat_analysis_both_EAD_DAD(time,Vm,Cai,CaSR,Nai,dV,1000,1,stimulation_train);
    data2{id+1} = function_analyze_DAD_only(time,Vm,Cai,CaSR,Nai,dV,stimulation_train,plot_figure);

%     data2{id+1} =  function_beat_analysis_DAD(time,Vm,Cai,CaSR,Nai,dV,1000,1,stimulation_train);
    
end

% out1 = sprintf('%s_EAD_count.csv', folder);
out2 = sprintf('%s_DAD_only.mat', folder);

% dlmwrite(out1, data1)
% dlmwrite(out2, data2)
save (out2, 'data2');
% outputs = [data1, data2];

end
