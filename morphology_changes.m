%% Calcuates Morphology Changes From Masked Images
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India
% Creates and saves the morphology changes parameters along with a
% representastive video

clear all;
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';

set(0, 'DefaultLineLineWidth', 2);
lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);


for drt_idx = 1:length(lst) % Iterating through all the folders
    dirt = strcat(moth_dirt,lst(drt_idx).name,'\');
    l1 = load(strcat(dirt,'Cell_Intensity\Thresh_Video.mat')); % Loading the video
    bndr_img = l1.bndr_img;
    cent_arr = l1.cent_arr;
    sv_dirt = strcat(dirt,'Morph_Quant');
    mkdir (sv_dirt) % Making folder
    prot = [];
    retr = [];
    tot_ar = [];
    cell_ar = [];
    morph_vid = VideoWriter(fullfile(sv_dirt,'\Morph_Video.avi'));
    morph_vid.FrameRate = 5;
    open(morph_vid)
    sub_img = cell(1);
    for frm_idx = 1:length(bndr_img)-1 % Looping through all the images
        ar_img = bndr_img{frm_idx}; % Actual Cell Area image
        sub = bndr_img{frm_idx+1} - bndr_img{frm_idx}; % Subtracted Frame
        sub_img{frm_idx} = sub; % Storing subtracted frame in a cell 
        sub(sub==0) = 128; % Making the file visually better for writing
        sub(sub==-1) = 0;
        sub(sub==1) = 255;
        writeVideo(morph_vid,uint8(sub)); % Written Video
        quad_crp_dim = [[cent_arr(frm_idx,1),0,size(bndr_img{1},2)-cent_arr(frm_idx,1),cent_arr(frm_idx,2)];...
                 [0,0,cent_arr(frm_idx,1),cent_arr(frm_idx,2)];[0,cent_arr(frm_idx,2),cent_arr(frm_idx,1),size(bndr_img{1},1)-cent_arr(frm_idx,2)];...
                 [cent_arr(frm_idx,1),cent_arr(frm_idx,2),size(bndr_img{1},2)-cent_arr(frm_idx,1),size(bndr_img{1},1)-cent_arr(frm_idx,2)]]; % Cropping Dimensions
        for quad = 1:4 % Looping through Quadrants
            crp_sub = imcrop(sub_img{frm_idx}, quad_crp_dim(quad,:)); % Cropping Subtracted image
            crp_ar = imcrop(ar_img, quad_crp_dim(quad,:)); % Cropping Area image
            prot(frm_idx, quad) = length(find(crp_sub==1)); % Protrusions
            retr(frm_idx, quad) = length(find(crp_sub==-1)); % Retractions
            tot_ar(frm_idx, quad) = nnz(crp_sub); % Total Area of Movement (Protrusion and Retractions)
            cell_ar(frm_idx, quad) = length(find(crp_ar==1)); % Cell Area Quad wise
        end     
    end
    
    prot(:,end+1) = sum(prot,2); % Summing it along Rows accounting for whole cell
    retr(:,end+1) = sum(retr,2); % Summing it along Rows accounting for whole cell
    tot_ar(:,end+1) = sum(tot_ar,2); % Summing it along Rows accounting for whole cell
    cell_ar(:,end+1) = sum(cell_ar,2); % Summing it along Rows accounting for whole cell
    
    norm_prot = prot ./ cell_ar; % Normalising by dividing each value by respective Quadrant Area or Cell Area
    norm_retr = retr ./ cell_ar; % Normalising by dividing each value by respective Quadrant Area or Cell Area
    norm_tot_ar = tot_ar ./ cell_ar; % Normalising by dividing each value by respective Quadrant Area or Cell Area
    
    save(strcat(sv_dirt,'\Quad_morph_data.mat'),'prot','retr','tot_ar',...
        'cell_ar','norm_prot','norm_retr','norm_tot_ar');
    save(strcat(sv_dirt,'\Morph_Vid.mat'),'sub_img');
    
    tab_var = array2table(norm_prot,...
        'VariableNames',{'Q1','Q2','Q3','Q4','Average'});
    writetable(tab_var,strcat(sv_dirt,'\Quad_Morph_Data.xlsx'), 'Sheet', 1);
    
    tab_var = array2table(norm_retr,...
        'VariableNames',{'Q1','Q2','Q3','Q4','Average'});
    writetable(tab_var,strcat(sv_dirt,'\Quad_Morph_Data.xlsx'), 'Sheet', 2);
    
    tab_var = array2table(norm_tot_ar,...
        'VariableNames',{'Q1','Q2','Q3','Q4','Average'});
    writetable(tab_var,strcat(sv_dirt,'\Quad_Morph_Data.xlsx'), 'Sheet', 3);   
    
    figure(1)
    clf
    plot(norm_tot_ar);
    legend('Q1','Q2','Q3','Q4','Average')
    saveas(gcf,strcat(sv_dirt,'\Quad_Morph_Total_Area_Plots.tif'));
end

set(0, 'DefaultLineLineWidth', 1);