%% Corelation of FRET and Morphodynamics
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India
% Performing corelation creates a folder in the directory saving the plots
% of individual cells and overall trend. Also saves in a excel file.
%Input a 0 or 1 based on the requirement to perform the computation on Before or After.


clear all;
clc


var_nm_arr = {'FRET_index'};
max_int = [300];
lvl = 1;

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';
choice = input('Before(0) or After(1)');

for diff_var = 1:length(var_nm_arr)
    clearvars -except moth_dirt max_int var_nm_arr diff_var choice lvl morph_cells
    fld = string(var_nm_arr{diff_var});
    
    if choice == 0
        sv_dirt = strcat(moth_dirt, 'Morpho_Analysis_Before\',fld,'\');
    elseif choice ==1
        sv_dirt = strcat(moth_dirt, 'Morpho_Analysis_After\',fld,'\');
    end

    mkdir (sv_dirt)
    
    lst = dir(strcat(moth_dirt,'\Cell*.*'));
    flag = [lst.isdir];
    lst = lst(flag);
    load(strcat(moth_dirt,'bef_idx.mat'));
    
    fret_morph_arr_all = [];
    tab_all = table;
    
    set(0, 'DefaultLineLineWidth', 1);
    
    for drt_idx = 1:length(lst)

        dirt = strcat(moth_dirt,lst(drt_idx).name,'\');
        l1 = load(strcat(dirt,'Cell_Intensity\',fld,'\Quad_data.mat'));
        quad_fret_data = l1.quad_fret_data;

        l2 = load(strcat(dirt,'Morph_Quant\Quad_morph_data.mat'));
        norm_tot_ar = l2.norm_tot_ar;
        if choice == 0
            quad_fret_data = quad_fret_data(1:before_idx(drt_idx),:);
            norm_tot_ar = norm_tot_ar(1:before_idx(drt_idx),:);
        elseif choice == 1
            quad_fret_data = quad_fret_data(before_idx(drt_idx)+1:end-1,:);
            norm_tot_ar = norm_tot_ar(before_idx(drt_idx)+1:end,:);
        end



        arr = [reshape(quad_fret_data(1:end,1:4),[],1),reshape(norm_tot_ar(:,1:4),[],1)];

        fret_morph_arr_all = [fret_morph_arr_all; arr];


        fret_arr = [1:1:max(ceil(arr(:,1)))]';
        morph_arr = accumarray (ceil(arr(:,1)),arr(:,2),[],@median);
        avl_fret = ~ismember( fret_arr, ceil(arr(:,1)));
        fret_arr(avl_fret) = NaN;
        morph_arr(avl_fret) = NaN;
        
        edit_arr = [ [fret_arr;nan(max_int(diff_var)-length(fret_arr),1)], [morph_arr;nan(max_int(diff_var)-length(morph_arr),1)] ];
        
        tab = array2table(edit_arr,'VariableNames',{strcat(lst(drt_idx).name,'_FRET'),strcat(lst(drt_idx).name,'_Area')});
        tab_all = [tab_all,tab];
        
        plot(fret_arr, morph_arr,'b.')
        xlim([0 max_int(diff_var)])
        ylim([0 0.1])
        xlabel('FRET Intensity');
        ylabel('Morphology');
        saveas(gcf, strcat(sv_dirt,'Graph_',lst(drt_idx).name,'.tif'))
    end
    
    
    
    writetable(tab_all,strcat(sv_dirt,'Excel_Morpho_FRET.xlsx'),'Sheet','Individual_Data');
    save(strcat(sv_dirt,'Data_Morpho_FRET.mat'), 'fret_morph_arr_all', 'tab_all');



    %% All data together
    
    e = [0:lvl:max_int];
    ds=discretize(fret_morph_arr_all(:,1),e);
    bin_data=[];

    for i=1:length(e)
        idx = find(ds==i);
        bin_data(i,1) = e(i)+(lvl)/2;
        bin_data(i,2) = median(fret_morph_arr_all(idx,2));
        bin_data(i,3) = std(fret_morph_arr_all(idx,2))/sqrt(length(idx));
    end
    idx_nan=find(isnan(bin_data(:,2))==1);
    bin_data(idx_nan,:)=[];
    figure
    errorbar(bin_data(:,1),bin_data(:,2),bin_data(:,3))
    xlabel(fld)
    ylabel('Morphodynamics')
    saveas(gcf, strcat(sv_dirt,'All_Data.tif'));
    saveas(gcf, strcat(sv_dirt,'All_Data.fig'));
    
    tab = array2table(bin_data,'VariableNames',{'FRET_T', 'Median', 'SEM'});
    writetable(tab,strcat(sv_dirt,'Excel_Morpho_FRET.xlsx'),'Sheet','All_Data_averg','WriteVariableNames',1);
    
    tab = array2table(fret_morph_arr_all, 'VariableNames',{'FRET_T', 'Morphology'});
    writetable(tab,strcat(sv_dirt, 'Excel_Morpho_FRET.xlsx'),'Sheet','All_Data','WriteVariableNames',1);
    

end
