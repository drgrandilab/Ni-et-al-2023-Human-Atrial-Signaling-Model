function checker(folder_name, t_range, x_select, y_select, NX, NY, z)
    ci  = zeros(length(x_select), length(y_select), length(t_range));
    % cp  = zeros(length(x_select), length(y_select), length(t_range));
    for idt = 1:length(t_range)
        t   = t_range(idt);
        FileLoc1     = ['../' folder_name '/BinaryFiles/map_ci_' num2str(t) '.ID.0.bin'];
        fileID1      = fopen(FileLoc1);
        rawfile_ci  = fread(fileID1, 'float32');
        fclose(fileID1);
        
        % if length(strfind(folder_name,'large_cell')) == 1
        % FileLoc2     = ['../' folder_name '/BinaryFiles/map_cp_' num2str(t) '.ID.0.bin'];
        % fileID2      = fopen(FileLoc2);
        % rawfile_cp  = fread(fileID2, 'float32');
        % fclose(fileID2);
        % end

                
        disp(t);
        
        for idx = 1:length(x_select)
            x = x_select(idx);
            for idy = 1:length(y_select)
                y = y_select(idy);
                n = 1 + (x-1) + (y-1)*NX + (z-1)*NX*NY;
                ci(idx, idy, idt)       = rawfile_ci(n);
                % if length(strfind(folder_name,'large_cell')) == 1
                % cp(idx, idy, idt)       = rawfile_cp(n);
                % end
            end
        end
    end

    model   = [folder_name '_'];
    save([model '.mat'],'ci', '-v7.3');
end