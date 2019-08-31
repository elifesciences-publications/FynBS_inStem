%% Quantifies Initial Data for channels.
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India.
% Initail Intensity of the whole cell.

clear all
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';
var_nm_arr = {'Acceptor','Donor','FRET_ch','FRET_index','FRET Computed','Norm FRET idx','Norm FRET computed'}; % This is very important, this is based on how the files has been saved in vid files.

lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);
set(0, 'DefaultLineLineWidth', 2);
init_data = [];
for drt_idx = 1:length(lst)
    drt_idx
    dirt = strcat(moth_dirt,lst(drt_idx).name,'\Cell_Intensity\');
    
     for i = 1:size(var_nm_arr,2) % Looping through different sets
         curr_dirt = strcat(dirt,var_nm_arr{i});
         l1 = load(strcat(curr_dirt,'\Quad_data.mat'));
         init_data(drt_idx,i) = l1.quad_fret_data(1,end); % end one is the average 1 is first time point
     end
end

kurt = kurtosis(init_data);
std = std(init_data/mean(init_data));

writing_var_nm = cell(1);
for i = 1:length(var_nm_arr)
    writing_var_nm{i} = strrep(var_nm_arr{i},' ','_');
end

tab = array2table(init_data,...
    'VariableNames',writing_var_nm,'RowNames',{lst(:).name});
 writetable(tab,strcat(moth_dirt,'Initial_Data.xlsx'),'WriteRowNames',true,'Sheet',1);
 
 tab = array2table(kurt,...
    'VariableNames',writing_var_nm,'RowNames',{'Kurtosis'});
 writetable(tab, strcat(moth_dirt,'Initial_Data.xlsx'),'WriteRowNames',true,'Sheet',2);
 
 save(strcat(moth_dirt,'Initial_Data.mat'),'init_data','kurt');