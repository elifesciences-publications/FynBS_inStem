%% Quantifies Quadrant Wise Data for channels. Perfoms masking.
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% Computes and saves the time lapse of cellular masks, the quadrant wise
% trace of intensity for each channel.

clear all;
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';
var_nm_arr = {'Acceptor','Donor','FRET_ch','FRET_Index','FRET Computed','Norm FRET idx','Norm FRET computed'};
lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);
ch_for_bndr = 1; % Which Channel to extract the boundary
% Parameters for extraction
open_var = 200; 
thr_val = 14;
se = strel('disk',5);
set(0, 'DefaultLineLineWidth', 2);


for drt_idx = 1:length(lst)
     
     disp(drt_idx) % Displays Progress
     dirt = strcat(moth_dirt,lst(drt_idx).name,'\');
     load(strcat(dirt,'video_files.mat'))
     sv_dirt_one = strcat(dirt,'Cell_Intensity');
     mkdir(sv_dirt_one);
     
     %% Boundary Extraction for Masking
     
     op_thr_vid = VideoWriter(fullfile(sv_dirt_one,'\Thresh_Video.avi'));
     op_thr_vid.FrameRate = 5;
     open(op_thr_vid)
     
     img = vid_files {ch_for_bndr};
     cent_arr = [];
     bndr_img = cell(1);
     for i = 1:length(img)
        im = pre_procs(img{i},open_var,thr_val/255); % Over writing pre processed version
        im_procc =  imerode(im,se); % Erosion (might do imopen also) to remove small irregularities on the surface
        im_procc = bwareaopen(im_procc,open_var);
        im_procc = bwareafilt(im_procc,1);
        bndr_img{i} = im_procc;
        cent = regionprops(bndr_img{i},'Centroid');
        cent = cent.Centroid;
        cent_arr(i,:) = round(cent);
        writeVideo(op_thr_vid,uint8(bndr_img{i}.*255));
     end
     close(op_thr_vid);
     save(strcat(sv_dirt_one,'\Thresh_Video.mat'),'bndr_img','cent_arr');
   
     
     imshow(bndr_img{1},[])
     
     
    
     %% Channel wise quantification
     if length(vid_files) ~= length(var_nm_arr)
         disp('Error!') % Throw error if Value names does nto match with variable
     end
     
     for i = 1:size(vid_files,1) %Looping through different sets
         sv_dirt = strcat(sv_dirt_one,'\',var_nm_arr{i});
         mkdir(sv_dirt);
         curr_vid  = vid_files{i};
         quad_fret_data=[];
         for k = 1:length(curr_vid) % Looping through time frames
             quad_crp_dim = [[cent_arr(k,1),0,size(im_procc,2)-cent_arr(k,1),cent_arr(k,2)];...
                 [0,0,cent_arr(k,1),cent_arr(k,2)];[0,cent_arr(k,2),cent_arr(k,1),size(im_procc,1)-cent_arr(k,2)];...
                 [cent_arr(k,1),cent_arr(k,2),size(im_procc,2)-cent_arr(k,1),size(im_procc,1)-cent_arr(k,2)]];
             img = uint8(bndr_img{k}).*curr_vid{k};
             for quad = 1:4 % Looping through quads
                 crp_im = imcrop(img,quad_crp_dim(quad,:));
                 crp_procc = imcrop(bndr_img{k},quad_crp_dim(quad,:));
                 [~, mean_var,rept]=quad_procc(crp_im, crp_procc);
                 quad_fret_data(k,quad) = mean_var;
             end
         end
         quad_fret_data(:,end+1) = mean(quad_fret_data,2); %Averaging Whole Cell
         
         save(strcat(sv_dirt,'\Quad_data.mat'),'quad_fret_data');
         xlswrite(strcat(sv_dirt,'\Quad_data_',lst(drt_idx).name,'.xlsx'),quad_fret_data);
         f = figure('visible','off');
         clf
         plot(quad_fret_data);
         legend('Q1','Q2','Q3','Q4','Average')
         saveas(f,strcat(sv_dirt,'\Quad_Plots_',lst(drt_idx).name,'.tif'));
     end
end