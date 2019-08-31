%% Analysis of Libe Scan Data Line Scan for membrane Fyn Localisation.
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% Performs the Scan Analysis and saves the data

clear all;
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';

var_nm_arr = {'Acceptor','Donor','FRET_ch','FRET index','FRET Computed','Norm FRET idx','Norm FRET computed'}; % This is very important, this is based on how the files has been saved in vid files.

lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);

lin_bin = 5;
band_extn = 20;
spat_dist = 100;
depth_memb = 10;

for drt_idx = 1:length(lst) % Looping through the directories
     disp(strcat('Directory is ',num2str(drt_idx)));
     
     dirt = strcat(moth_dirt,lst(drt_idx).name,'\');
     sv_dirt = strcat(dirt, 'Line_Scan\');
     l1 = load(strcat(sv_dirt,'Line_int_data.mat'));
     perm_pt_cell = l1.perm_pt_cell;
     intn_cell = l1.intn_cell;
     profl_cord_cell = l1.profl_cord_cell;
     
     l2 = load(strcat(dirt,'Cell_Intensity\Thresh_Video.mat'));
     bndr_img = l2.bndr_img;
     cent_arr = l2.cent_arr;
     
     %% Finding the most poalrised part
     empt_idx = cellfun(@isempty, intn_cell);
     intn_cell(empt_idx) = {nan(100,1)}; % Replacing empty places with nan elements, 
     near_memb_int = cellfun(@(x) x(1:depth_memb,1),intn_cell,'UniformOutput',false);
     near_memb_int = cellfun(@nanmean, near_memb_int);  % Calculating near membrane intensities
     [~,mx_idx] = max(near_memb_int'); % Maximum intensity line is roughly which line at each time
     
     mx_int_ln = histcounts(mx_idx,[0:lin_bin:size(intn_cell,2)]);
     [~,idx] = max(mx_int_ln); % finding the line "region(bin length 5)" which happens to the most polarised line most time
     mx_int_ln = idx*lin_bin; % finding the line "region(bin length 5)" which happens to the most polarised line most time
     band = [mx_int_ln-band_extn, mx_int_ln+band_extn]; % Defining a region based on most polarised band
     
     %% Computing_Intensities
     band_intn = intn_cell(:,band(1):band(2)); % Values at that band
     
     A = cellfun(@(x) [x;nan(spat_dist,1)], band_intn,'UniformOutput',false);  % Concating false values so as to each element is greater than a certain threshold
     A = cellfun(@(x) x(1:spat_dist,1), A, 'UniformOutput',false); % cutting till a point from the membrane
     mean_arr = []; % Calculating mean array
     for i = 1:size(A,1)
          M = cat(3,A{i,:});
          mean_arr(i,:) = nanmean(M,3);
     end
     
     %% Plotting Data
     
     % 3D plot
     colormap hot
     figure(1)
     clf
     surf(mean_arr)
     xlabel('Dist From Edge');
     ylabel('Time');
     zlabel('FRET');
     colorbar
     saveas(gcf, strcat(sv_dirt,'3Dgraph.fig'));
     
     % color coded plot
     figure(1)
     clf
     imagesc((mean_arr))
     xlabel('Distance From Membrane(Pixels)')
     ylabel('Time')     
     colorbar
     set(gca,'FontSize',16)
     saveas(gcf, strcat(sv_dirt,'Graph.tif'));
     
    % Signifying part of periphery picture
    figure(1)
    clf
    imshow(bndr_img{1})
    perm_pt = perm_pt_cell{1};
    hold on
    line([perm_pt(band(1),1),cent_arr(1,1)],[perm_pt(band(1),2), cent_arr(1,2)],'Color','b', 'LineWidth',3)
    hold on
    line([perm_pt(band(2),1),cent_arr(1,1)],[perm_pt(band(2),2), cent_arr(1,2)],'Color','b', 'LineWidth',3)
    saveas(gcf, strcat(sv_dirt,'Area_part.tif'));
    
    
    save(strcat(sv_dirt,'line_scan_analysis.mat'), 'band', 'band_intn', 'mx_int_ln', 'mean_arr', 'near_memb_int');
    xlswrite(strcat(sv_dirt,'Line_Scan.xlsx'), mean_arr);
    
end