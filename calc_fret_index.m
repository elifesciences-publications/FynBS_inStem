%% Calcuates FRET index from bleed through values and Acceptor Donor and FRET channel images
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% The code displays the progress by showing the iteration number.
% The generates the .mat of all the channels after calculating the FRET values.
 
clc
clear all

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\'; % Cell data
drt_nm = input('Directory name(Result)', 's');
sv_dirt = fullfile(moth_dirt, drt_nm);
mkdir(sv_dirt);

lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag); %List all Files in directory with Cell as starting

var_nm_arr = {'acceptor','donor','FRET_ch','FRET_Index'}; % Acceptor Donor FRETch FRET idx
a_bt = 0.028; % Acceptor Bleed through
d_bt = 0.55; % Donor Bleed through
md_filt = 2; % Median Filtering result

for i = 1 : length(lst) % Iterating through all the files
    disp(i) % Displaying Progress
    
    nm = lst(i).name;
    
    fld_nm = strcat(moth_dirt,'\',nm);
    save_fld_nm = fullfile(sv_dirt, nm);
    mkdir (save_fld_nm)
    dirt = dir(fld_nm);
    
    %% Setting up video files and parameters
    op_fret_comp = VideoWriter(fullfile(save_fld_nm,'FRET Computed.avi'),'Motion JPEG AVI');
    op_norm_fret_comp = VideoWriter(fullfile(save_fld_nm,'FRET Norm computed.avi'),'Motion JPEG AVI');
    op_norm_fret_idx = VideoWriter(fullfile(save_fld_nm,'FRET Norm Index.avi'),'Motion JPEG AVI');
    
    op_fret_comp.FrameRate = 5;
    op_norm_fret_comp.FrameRate = 5;
    op_norm_fret_idx.FrameRate = 5;
    
    open(op_fret_comp)
    open(op_norm_fret_comp)
    open(op_norm_fret_idx)
    %% Looping to calcualte the file names
    
    vid_files = cell(1);
    for j=1:length(var_nm_arr)
        for k=3:length(dirt)
            
            nm = dirt(k).name;
            if contains(nm,strcat(var_nm_arr{j}), 'Ignorecase', 1)==1
                
                info=imfinfo(strcat(fld_nm,'\',nm));
                num = numel(info);
                vid =cell(1);
                for g=1:num
                    vid{g,1} = imread(strcat(fld_nm,'\',nm),g);
                end
                vid_files{j,1} = vid;
                break
            end 
        end
    end
    
    
    if sum(cellfun(@isempty,vid_files))~=0 % Error checking step
        disp('One of the cell folder could not be read, no match found in files');
        pause
    end
    
    
    %% Calculating FRET Index
    
    num_im = length(vid_files{1});
    fret_comp = cell(1);
    norm_fret_comp = cell(1);
    norm_fret_idx = cell(1);
    
    
    for j = 1:num_im
        fret_comp{j,1} = vid_files{3}{j} - a_bt*vid_files{1}{j} - d_bt*vid_files{2}{j};
        fret_comp{j,1} = medfilt2(fret_comp{j,1},[md_filt, md_filt]);
        writeVideo(op_fret_comp,fret_comp{j,1})
        norm_fret_comp{j,1} = medfilt2(fret_comp{j,1}./vid_files{2}{j},[md_filt, md_filt]);
        writeVideo(op_norm_fret_comp,norm_fret_comp{j})
        norm_fret_idx{j,1} = medfilt2(vid_files{4}{j}./vid_files{2}{j},[md_filt, md_filt]);
        writeVideo(op_norm_fret_idx,norm_fret_idx{j})
    end
    
    vid_files{end+1} = fret_comp;
    vid_files{end+1} = norm_fret_idx;
    vid_files{end+1} = norm_fret_comp;
    
    save(strcat(save_fld_nm,'\','video_files.mat'),'vid_files'); % Acceptor Donor FRETch FRETidx FRETcomp NormFRETidx NormFRETcomp
    
    close(op_fret_comp)
    close(op_norm_fret_comp)
    close(op_norm_fret_idx)
end
