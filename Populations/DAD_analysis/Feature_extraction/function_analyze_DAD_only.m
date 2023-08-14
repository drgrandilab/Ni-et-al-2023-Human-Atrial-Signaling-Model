function [outputs] = function_analyze_DAD_only(time,Vm,Ca,CaSR,Na,dVm,stimulation_train,plot_figure)

run_delta_Na = Na(end)-Na(1);

% Check for impaired repolarization
flag_abnormal_rep = 0;

run_Vm_min = min(Vm);
if run_Vm_min > -50
    flag_abnormal_rep = 1;
else
    for i = 1:(length(stimulation_train)-1)
        %         i
        
        ti_roi = find(time>stimulation_train(i),1,'first'); ti_index = ti_roi(1);
        
        
        te_roi = find(time<stimulation_train(i+1),1, 'last'); te_index = te_roi(end);
        
        time_roi = time(ti_index:te_index);
        %         time_roi(end)
        Vm_roi = Vm(ti_index:te_index);
        %Ca_roi = Ca(ti_index:te_index);
        %CaSR_roi = CaSR(ti_index:te_index);
        %Na_roi = Na(ti_index:te_index);
        %dVm_roi = dVm(ti_index:te_index);
        
        %         Vm_end_array(i) = Vm_roi(end);
        Vm_end_array(i) = min(Vm_roi);
    end
    
    
    %     test = 1;
    if max(Vm_end_array) > -50
        flag_abnormal_rep = 1;
    end
end



% plot_figure
% if(plot_figure)
%
% % figure(2), plot(time, Vm); hold on;
% end
% Check for EADs
flag_ead = 0;
DAD_ind = {};
total_DAD = 0;
for i = 1:(length(stimulation_train)-1)
    
    
    ti_roi = find(time>stimulation_train(i)-1,1,'first'); ti_index = ti_roi(1);
    
    
    te_roi = find(time<stimulation_train(i+1),1, 'last'); te_index = te_roi(end);
    
    new_index = find(Vm(ti_index:te_index) > -75, 1, 'last');
    
    % if(length(new_index) == 1)
        
    %     te_index = min(te_index, ti_index+new_index);
    % end
    
    
    time_roi = time(ti_index:te_index);
    Vm_roi = Vm(ti_index:te_index);
    Ca_roi = Ca(ti_index:te_index);
    %Ca_roi = Ca(ti_index:te_index);
    %CaSR_roi = CaSR(ti_index:te_index);
    %Na_roi = Na(ti_index:te_index);
    dVm_roi = dVm(ti_index:te_index);
    
    [pks,locs] = findpeaks(Vm_roi,'MinPeakDistance',10);
    %     [~, locs_dv]  = findpeaks(dVm_roi, 'MinPeakHeight',12.5,'MinPeakDistance',10);
    [~, locs_dv]  = findpeaks(dVm_roi, 'MinPeakHeight',8,'MinPeakDistance',10);
    
    if(plot_figure)
        
        figure(3), plot(time_roi, dVm_roi); hold on;
    end
    
    
    MDP = min(Vm_roi);

    Ca_diast = min(Ca_roi);
    if(MDP>-65)
        MDP=-65;
        disp('MDP higher than -65 mV');
    end
    
    
    DAD_num = 0;
    EAD_num = 0;
    dad_array = [];
    DAD_cutoff_amp = 2.0;  % DAD cutoff amplitude for analysis
    if ~isempty(pks)
        for j = 1:length(pks)
            
            tm = time_roi(locs(j));
            % // 14:01:36, Sun, 07-June-2020, By Haibo adding threshold for each cell instead of -65mV;
            if(tm - time_roi(1) > 10) && (pks(j) > MDP+DAD_cutoff_amp)  && (min(Vm_roi(locs(max(1, j-1)):locs(j))) < -65)  % DAD should be lower than -60
                DAD_num = DAD_num +1;
                total_DAD= total_DAD+1;
                incid.time = tm;
                incid.AMP = pks(j) - MDP;
%                 incid.Ca_AMP = max(Ca_roi(locs(max(1, j-1)):locs(j))) - Ca_diast;
                incid.Ca_AMP = Ca_roi(locs(j)) - Ca_diast;



                DAD_ind{total_DAD} = incid;

                if(plot_figure)
                    figure(2), plot(time_roi, Vm_roi), hold on;
                    plot(tm, Vm_roi(locs(j)), 'rs');
                    
                    figure(4), plot(time_roi, Ca_roi), hold on;
                    plot(tm, Ca_roi(locs(j)), 'rs');
                end
            end
        end
        %             i
        %             tm
        
        
        
        if(~isempty(locs_dv))  % there is at least one AP here.
            if( abs(time_roi(1)-time_roi(locs_dv(1)))>10)
                locs_dv = [1; locs_dv];
            end
            
            for j = 1:length(pks)
                
                tm = time_roi(locs(j));
                
                
                index_dv = find(time_roi(locs_dv) < tm,1,'last');  % see which AP dvdtmax is close to the current peak
                
                if(isempty(index_dv))
                    disp('index_dv = []');
                    index_dv = 1;
                    
                end
                if (    (tm - time_roi(locs_dv(index_dv)) > 150) && (pks(j) > -60)  ... 
                    && (min(Vm_roi(locs(max(1, j-1)):locs(j))) > -65 ) ...
                    && (min( Vm_roi(locs_dv(index_dv):locs(j)) ) > -65)  ) ...
                        || (   (tm - time_roi(locs_dv(index_dv)) > 50) && (pks(j) > -60) ...
                    && (min( Vm_roi(locs_dv(index_dv):locs(j)) ) > -65)  ...     %                         && (min(Vm_roi(locs(max(1, j-1)):locs(j))) > -65) && ...
                       && (pks(j) - min(Vm_roi(locs(max(1, j-1)):locs(j))) > 10) ...
                        && (max( dVm_roi(locs(max(1, j-1)) :locs(j)) ) > 1)   )   %% maxium dvdt > 1 for EAD 
                    EAD_num = EAD_num +1;
                    if(plot_figure)
                        % tm
                        figure(2), plot(time_roi, Vm_roi), hold on;
                        plot(tm, pks(j), 'bs');
                    end
                end
            end
        else
            EAD_num = 0;
        end
    end
    
    if DAD_num > 0
        dad_array(i) = DAD_num;  % num of DADs
    else
        dad_array(i) = 0;
    end
    if EAD_num > 0
        ead_array(i) = 1;
    else
        ead_array(i) = 0;
    end
    
end
if sum(dad_array)>0
    flag_ead = 1;
    %         sum(ead_array)
end

flag_normal_AP = 1;
% Beat analysis
if flag_abnormal_rep == 1 || flag_ead == 1
    flag_normal_AP = 0;
end

output = [flag_abnormal_rep flag_ead flag_normal_AP run_Vm_min ...
    run_delta_Na, sum(dad_array), sum(ead_array), length(ead_array), ead_array...
    ];


outputs.EAD_info = output;
outputs.DAD_ind = DAD_ind;
end