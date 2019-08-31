%% Performs Line Scan for membrane Fyn Localisation.
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India. 
% Creates a subfolder and saves line scan data. Analysis file will read the
% .mat file and save it

clear all;
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';

lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);
ch_for_quant = 5;


for drt_idx = 1:length(lst)
      dirt = strcat(moth_dirt,lst(drt_idx).name,'\');
      l1 = load(strcat(dirt,'video_files.mat'));
      vid_files = l1.vid_files;
      sv_dirt = strcat(dirt,'Line_Scan\');
      mkdir(sv_dirt)
      
      l2 = load(strcat(strcat(dirt,'\Cell_Intensity\Thresh_Video.mat')));
      cent_arr = l2.cent_arr;
      bndr_img = l2.bndr_img;
      perm_pt_cell = cell(1);
      intn_cell = cell(1);
      profl_cord_cell = cell(1);
      
      for num_frm = 1:size(cent_arr,1)
          disp(strcat('Number of Frame is ', num2str(num_frm), 'Directory ', num2str(drt_idx)));
          cent = cent_arr(num_frm,:);
          x_cent = cent(1);
          y_cent = cent(2);
          bw_img = bndr_img{num_frm};
          data_img = uint8(bw_img) .* vid_files{ch_for_quant}{num_frm};
          %[perm_pt, ~] = contour(bw_img);
          perm_pt = bwboundaries(bw_img);
          perm_pt = fliplr(perm_pt{1}); % flip is required as the the way improfile and bwboundaries works
          %perm_pt = 512 - perm_pt;
          perm_pt_cell{num_frm} = perm_pt;
          figure(1)
          imshow(bw_img);
          arr = zeros(size(perm_pt,1),4);
          x_perm_pt = perm_pt(:,1);
          y_perm_pt = perm_pt(:,2);
          
         for pt_num = 1:size(perm_pt,1)
              x_bnd = x_perm_pt(pt_num);
              y_bnd = y_perm_pt(pt_num);
              xi = [x_bnd x_cent];
              yi = [y_bnd y_cent];
              arr(pt_num,:) = [xi yi];
              prof = improfile(data_img, xi, yi);
              intn_cell{num_frm, pt_num} = prof;
          end
          profl_cord_cell{num_frm} = arr;
      end
      save(strcat(sv_dirt, 'Line_int_data.mat'), 'intn_cell','perm_pt_cell','profl_cord_cell');    
end