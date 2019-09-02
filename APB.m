%% Computes Acceptor PhotoBleaching
% Code developed by Shilpa Dilip Kumar(Email id: shilpadilip@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% Matlab code written for analyzing data obtained in acceptor-photobleaching experiment
close all;
clear all;
rep_nmb = 2;
moth_dirt = 'Z:\Imag_Data\Test_temp\';
img_sv = strcat(moth_dirt,'Cell_Outline\');
mkdir(img_sv);
lst = dir(strcat(moth_dirt,'*.tif'));
tot_fil=length(lst)/rep_nmb;
load(strcat(moth_dirt,'thr_res.mat'))
delete(strcat(moth_dirt,'APB.xls'))

for fil_idx=1:tot_fil
    clearvars -except lst tot_fil fil_idx moth_dirt thr_res_arr rep_nmb img_sv
    for i=1:rep_nmb
         nm=lst((fil_idx-1)*rep_nmb+i).name;
          if length(strfind(nm,'pre'))==1
              nm_pre=nm;
          elseif length(strfind(nm,'post'))==1
              nm_post=nm;
          end
    end
   pre_bleach=imread(strcat(moth_dirt,nm_pre));
   if size(pre_bleach,3)==1
       pre_bleach=pre_bleach(:,:,1);
   end
   post_bleach=imread(strcat(moth_dirt,nm_post));
   figure;imshow(pre_bleach,[]); title('Pre Bleach')
   
   if size(post_bleach,3)==1
       post_bleach=post_bleach(:,:,1);
   end
   figure;imshow(post_bleach,[]); title('Post Bleach')
    
   
   diff_img = post_bleach - pre_bleach;
   figure;imshow(diff_img,[]);
   T=thr_res_arr(fil_idx,1)/2^16;
   Im_bw = im2bw(post_bleach,T);

   figure;imshow(Im_bw,[])
   
   
   fh = imfill(Im_bw,'holes');

Im_L = bwlabeln(fh, 8);

S = regionprops(Im_L, 'Area');

Im_regProp = ismember(Im_L, find([S.Area] >= 2500));
Im_regProp=bwareafilt(Im_regProp,1);
figure;imshow(Im_regProp,[])

imwrite(Im_regProp.*255,strcat(img_sv,'Cell_',num2str(fil_idx),'Boundary.tif'));

 BW1 = edge(Im_regProp,'canny');

figure;imshow(BW1,[]);

S1 = regionprops(Im_regProp, 'Centroid','MajorAxisLength','PixelIdxList','PixelList');

centroid = S1.Centroid;

% pixels =S1.PixelList;
% figure;
[C,h]=contour(Im_regProp,1);
hold on;plot(centroid(1), centroid(2),'g+');
tic
figure;imshow(post_bleach,[]);
for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(post_bleach,xi,yi,200);

    int(:,i) = c;
    if mod(i,50)==0
        hold on;
        plot(xi,yi,'r')
    end
end

for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(pre_bleach,xi,yi,200);

    int1(:,i) = c;

end

for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(diff_img,xi,yi,200);

    int2(:,i) = c;

end

toc

mean_int_postbleach = mean(int,2);
mean_int_prebleach = mean(int1,2);
mean_int_diffimg = mean_int_postbleach - mean_int_prebleach;
% mean_int_diffimg = mean(int2,2);


figure; plot(mean_int_prebleach);
hold on; plot(mean_int_postbleach,'r');
hold on; plot(mean_int_diffimg,'g');

legend('pre bleach','post bleach','difference');
saveas(gcf,strcat(img_sv,'Cell_',num2str(fil_idx),'Graph.tif'));


xlswrite(strcat(moth_dirt,'APB.xls'),[mean_int_prebleach,mean_int_postbleach,mean_int_diffimg],fil_idx)

 close all
    
end

