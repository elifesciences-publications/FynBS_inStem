%% Quantification of Polarisation, Extent Of Oscillation and Pulsitility
% Code developed by Sayan Biswas(Email id: sayanbiswasjuee@gmail.com), Dr. Akash Gulyani's Lab at
% Institute for Stem Cell Science and Regenrative Medicine, Bangalore
% India
% Quantifies the parameters - Polarisation, Extent of Oscillation and
% Oscillation period and write in a excel file.


clear all;
clc

moth_dirt = 'Y:\Imag_Data\U2OS_Cell\Result\';

var_nm_arr = {'Acceptor','Donor','FRET_ch','FRET_Index','FRET Computed','Norm FRET idx','Norm FRET computed'}; % This is very important, this is based on how the files has been saved in vid files.
polr_chnl = 5;
oth_polr_ch = [1, 2, 3];

lst = dir(strcat(moth_dirt,'\Cell*.*'));
flag = [lst.isdir];
lst = lst(flag);
load(strcat(moth_dirt,'bef_idx.mat'));
time = 1; % Time gap between frames in minutes

fold = 'Cell';
polr_arr = [];
oth_ch_polr_arr = cell(length(oth_polr_ch), 1);
ext_oscc = [];
oscc_tp = [];
periodogram_cell = cell(1);


for drt_idx = 1:length(lst)
    disp(drt_idx)
    dirt = strcat(moth_dirt,lst(drt_idx).name);
    l2 = load(strcat(dirt,'\',fold,'_Intensity\',var_nm_arr{polr_chnl},'\Quad_data.mat')); % loading the required polarised chanel
    quad_fret_data = l2.quad_fret_data;
    bef_frame = before_idx(drt_idx);
    bef = quad_fret_data(1:bef_frame,:);
    aft = quad_fret_data(bef_frame+1:end,:);
    
    if isempty(aft) == 1
        quad_fret_data = [quad_fret_data; quad_fret_data];
        aft = quad_fret_data(bef_frame+1:end,:);
    end
     %% Finding the polarisation
    [frm_bf,q_bf_mx]=find(bef==max(bef(:)));
    [frm_af,q_af_mx]=find(aft==max(aft(:)));
    
    frm_af = frm_af + bef_frame; % Indexing on the total time series
    pol_bef = [max(quad_fret_data(frm_bf,:)),min(quad_fret_data(frm_bf,:))];
    q_bf_mn = find(quad_fret_data(frm_bf,:)==min(quad_fret_data(frm_bf,:)));
    q_af_mn = find(quad_fret_data(frm_af,:)==min(quad_fret_data(frm_af,:)));
    
    pol_aft = [max(quad_fret_data(frm_af,:)),min(quad_fret_data(frm_af,:))];
    polr_arr(drt_idx,:) = [frm_bf, q_bf_mx, q_bf_mn, pol_bef, pol_bef/pol_bef(2), frm_af, q_af_mx, q_af_mn, pol_aft, pol_aft/pol_aft(2)];
    
    % Storing Polarisation for other channel
    for oth_ch = 1:length(oth_polr_ch)
        l = load(strcat(dirt,'\',fold,'_Intensity\',var_nm_arr{oth_polr_ch(oth_ch)},'\Quad_data.mat')); % loading therequired polarised chanel
        data = l.quad_fret_data;
        
        polr_bef = [data(frm_bf, q_bf_mx), data(frm_bf, q_bf_mn)];
        polr_aft = [data(frm_af, q_af_mx), data(frm_af, q_af_mn)];
        
        arr = [frm_bf, q_bf_mx, q_bf_mn, polr_bef, polr_bef/polr_bef(2), frm_af, q_af_mx, q_af_mn, polr_aft, polr_aft/polr_aft(2)];
        oth_ch_polr_arr{oth_ch, 1} = [oth_ch_polr_arr{oth_ch, 1}; arr ];
    end
    %% Extent of Oscillation
    bef_oscc = std((bef));
    aft_oscc = std((aft));
    [~,idx_bef] = sort (mean((bef)), 'ascend');
    [~,idx_aft] = sort (mean((aft)), 'ascend');
    
    ext_oscc(drt_idx,:) = [bef_oscc,aft_oscc, idx_bef(end), idx_bef(1), idx_aft(end), idx_aft(1)];
    
    %% Pulsitility of Oscillation
    oscc_dirt = strcat(dirt,'\Periodogram_Data\');
    mkdir (oscc_dirt)
    
    % Before
    fs = 1./time;
    bef_oscc = detrend(bef,2,'Continuous',0);
    [bef_pxx,bef_f] = periodogram(bef_oscc,[],[],fs);
    [mx,~] = find(bef_pxx==max(bef_pxx));
    bef_tp = [1./bef_f(mx)]';
    
    % Before_Periodogram_plotting
    f = figure('visible','off');
    clf
    plot(bef_f, bef_pxx)
    xlabel('Frequency(s^{-1})')
    legend('Q1','Q2','Q3','Q4','Average')
    saveas(f,strcat(oscc_dirt,'Before_Periodogram.tif'))
    
    % Before_Proccesed_Data_plotting
    f = figure('visible','off');
    clf
    plot(bef_oscc)
    xlabel('Time')
    ylabel('Intensity')
    legend('Q1','Q2','Q3','Q4','Average')
    saveas(f,strcat(oscc_dirt,'Before_Procc_data.tif'))
    
    % After
    aft_oscc = detrend(aft,2,'Continuous',0);
    [aft_pxx,aft_f] = periodogram(aft_oscc,[],[],fs);
    [mx,~] = find(aft_pxx==max(aft_pxx));
    aft_tp = [1./aft_f(mx)]';
    
    % After_Periodogram_plotting
    f = figure('visible','off');
    clf
    plot(aft_f, aft_pxx)
    xlabel('Frequency(s^{-1})')
    legend('Q1','Q2','Q3','Q4','Average')
    saveas(f,strcat(oscc_dirt,'After_Periodogram.tif'))
    
    % After_Proccesed_Data_plotting
    f = figure('visible','off');
    clf
    plot(aft_oscc)
    xlabel('Time')
    ylabel('Intensity')
    xlabel('Frequency(s^{-1})')
    legend('Q1','Q2','Q3','Q4','Average')
    saveas(f,strcat(oscc_dirt,'After_Procc_data.tif'))
    
    oscc_tp(drt_idx,:) = [bef_tp, aft_tp, idx_bef(end), idx_bef(1), idx_aft(end), idx_aft(1)];
    periodogram_cell(drt_idx,1:6) = {bef_oscc, aft_oscc, bef_pxx, bef_f, aft_pxx, aft_f};
end

   t_polr = array2table([polr_arr, before_idx'],...
        'VariableNames',{'Frame_Bef','Quad_Bef_max','Quad_Bef_min','Before_Max_Val','Before_Min_Val','Norm_Bef_Mx','Norm_Bef_Min'...
        'Frame_Aft','Quad_Aft_max','Quad_Aft_min','After_Max_Val','After_Min_Val','Norm_Aft_Mx','Norm_Aft_Min','Before_Index'},'RowNames',{lst(:).name});
   writetable(t_polr, strcat(moth_dirt,'\Polarisation_',fold,'_',var_nm_arr{polr_chnl} ,'.xlsx'),'WriteRowNames',true, 'Sheet', 'Polarisation');
   save(strcat(moth_dirt,'\Polarisation_',fold,'_',var_nm_arr{polr_chnl},'.mat'),'polr_arr', 'oth_ch_polr_arr')
   
   
   for i = 1:length(oth_polr_ch)
        t_polr = array2table([oth_ch_polr_arr{i}, before_idx'],...
        'VariableNames',{'Frame_Bef','Quad_Bef_max','Quad_Bef_min','Before_Max_Val','Before_Min_Val','Norm_Bef_Mx','Norm_Bef_Min'...
        'Frame_Aft','Quad_Aft_max','Quad_Aft_min','After_Max_Val','After_Min_Val','Norm_Aft_Mx','Norm_Aft_Min','Before_Index'},'RowNames',{lst(:).name});
         writetable(t_polr, strcat(moth_dirt,'\Polarisation_',fold,'_',var_nm_arr{polr_chnl},'.xlsx'),'WriteRowNames', true, 'Sheet', strcat('Polarisation_',var_nm_arr{oth_polr_ch(i)})); 
   end
   
   
   t_ext_oscc = array2table([ext_oscc, before_idx'],...
       'VariableNames',{'Bef_Q1','Bef_Q2','Bef_Q3','Bef_Q4','Bef_Average',...
       'Aft_Q1','Aft_Q2','Aft_Q3','Aft_Q4','Aft_Average','Max_bef','Min_bef','Max_aft','Min_aft','Before_Index'},'RowNames',{lst(:).name});
   writetable(t_ext_oscc, strcat(moth_dirt,'\Extent_Oscillation_',fold,'_',var_nm_arr{polr_chnl},'.xlsx'),'WriteRowNames',true, 'Sheet', 'Ext_Oscc');
   save(strcat(moth_dirt,'\Extent_Oscillation_',fold,'_',var_nm_arr{polr_chnl},'.mat'),'ext_oscc');
   
   t_tp = array2table([oscc_tp, before_idx'],...
       'VariableNames',{'Bef_Q1','Bef_Q2','Bef_Q3','Bef_Q4','Bef_Average',...
       'Aft_Q1','Aft_Q2','Aft_Q3','Aft_Q4','Aft_Average','Max_bef','Min_bef','Max_aft','Min_aft','Before_Index'},'RowNames',{lst(:).name});
   writetable(t_tp, strcat(moth_dirt,'\Time_Period_',fold,'_',var_nm_arr{polr_chnl},'.xlsx'),'WriteRowNames',true, 'Sheet', 'Time_Period');
   save(strcat(moth_dirt,'\Time_Period_',fold,'_',var_nm_arr{polr_chnl},'.mat'),'oscc_tp');
   save(strcat(moth_dirt,'\Periodogram_Data_',fold,'_',var_nm_arr{polr_chnl},'.mat'),'periodogram_cell', 'oscc_tp');