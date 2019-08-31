%% Quantifies the Maximum FRET point - and other channels as well
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% Generates parameter of the time trace of highest FRET quadrant.


clear all;
clc
moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';
var_nm_arr = {'Acceptor','Donor','FRET Computed'};
lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);
load(fullfile(moth_dirt,'bef_idx.mat'));

inten_ch = 3;

fold = 'Cell';
val_arr = [];

for drt_idx = 1:length(lst)
       disp(drt_idx)
       dirt = fullfile(moth_dirt, lst(drt_idx).name, strcat(fold,'_Intensity'), var_nm_arr{inten_ch});
       l = load(fullfile(dirt,'Quad_data.mat'));
       quad_fret_data = l.quad_fret_data;
       bef_frame = before_idx(drt_idx);
       bef = quad_fret_data(1:bef_frame,:);
       aft = quad_fret_data(bef_frame+1:end,:);
       
       [frm_bf,q_bf_mx] = find(bef==max(bef(:)));
       bef_mx_qd = bef(:, q_bf_mx);
       [~, frm_bf_qd_mx] = max(bef_mx_qd);
       [~, frm_bf_qd_mn] = min(bef_mx_qd);
       
       [frm_af,q_af_mx] = find(aft==max(aft(:)));
       aft_mx_qd = aft(:, q_af_mx);
       [~, frm_af_qd_mx] = max(aft_mx_qd);
       [~, frm_af_qd_mn] = min(aft_mx_qd);
       
       % Making the index
       idx_bef_mx = [frm_bf_qd_mx, q_bf_mx];
       idx_bef_mn = [frm_bf_qd_mn, q_bf_mx];
       idx_aft_mx = [frm_af_qd_mx+bef_frame, q_af_mx];
       idx_aft_mn = [frm_af_qd_mn+bef_frame, q_af_mx];
       idx_arr(drt_idx,:) = [idx_bef_mx, idx_bef_mn, idx_aft_mx, idx_aft_mn];
       
       for ch_idx = 1:length(var_nm_arr)
           dirt = fullfile(moth_dirt, lst(drt_idx).name, strcat(fold,'_Intensity'), var_nm_arr{ch_idx});
           l = load(fullfile(dirt,'Quad_data.mat'));
           quad_fret_data = l.quad_fret_data;
        
           for i = 1:size(idx_arr(drt_idx,:),2)/2
               val_arr(drt_idx, i, ch_idx) = quad_fret_data(idx_arr(drt_idx,2*i-1),idx_arr(drt_idx,2*i));
           end
       end
end

for ch_idx = 1:length(var_nm_arr)
    tab = array2table([idx_arr, val_arr(:,:,ch_idx)], 'VariableNames', ...
        {'Bef_mx_tm', 'Bef_mx_qd', 'Bef_mn_tm', 'Bef_mn_qd','Aft_mx_tm', 'Aft_mx_qd', 'Aft_mn_tm', 'Aft_mn_qd', 'Bef_Mx', 'Bef_Min', 'Aft_Mx', 'Aft_Min'},...
        'RowNames', {lst(:).name});
    writetable(tab, fullfile(moth_dirt, 'HQ_FRET_parm.xlsx'), 'Sheet', var_nm_arr{ch_idx}, 'WriteRowNames', 1);
end

save(fullfile(moth_dirt, 'HQ_FRET_parm.mat'), 'val_arr', 'idx_arr');
