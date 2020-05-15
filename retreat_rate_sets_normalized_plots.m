%Plot the normalized retreat rates for the BRIE experiments
%normalized by background RSLR rate
clc; close all; clear

%%
resultsPath = 'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_125um';

fdir = dir([resultsPath '\*.mat']);

%% Plot normalized rates for each model SLR experiment on 1 plot
m1 = figure('color','w'); 

baseline_name = 'Sea Level Experiments 9 mm per yr longshore on.mat'; %baseline name to skip

path_list = {'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\d_sh_10m';...
    'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR';...
    'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\d_sh_14m';...
    'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_100um';...
    'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_125um'};

title_list = {'d_s_h = 10m gs = 160um';...
    'd_s_h = 20m gs = 160um';...
    'd_s_h = 14m gs = 160um';...
    'd_s_h = 20m gs = 100um';...
    'd_s_h = 20m gs = 125um'};

for ii = 1:length(path_list) %make plot for every set of conditions (eg grain size, shoreface depth)
    
    resultsPath = path_list{ii};
    fdir = dir([resultsPath '\*.mat']);
    
    %save normalized mean for each retreat rate to fit a line
    shoreline_mean = [];
    shoreface_mean = [];
    slr_rate = [];
    
    for n = 1:length(fdir)
        
        results = load([resultsPath '\' fdir(n).name]); %load individual SLR model experiment
        
        if ~strcmp(fdir(n).name,baseline_name) %skip the baseline rate
            
            shoreline_mean = [shoreline_mean; results.mean_dxs];
            shoreface_mean = [shoreface_mean; results.mean_dxt];
            slr_rate = [slr_rate; results.paramValues(1,1)*1000];
            
            subplot(2,3,ii)
            hold on
            errorbar(results.paramValues(1,1)*1000,results.mean_dxs,results.stdev_dxs,'bo')
            errorbar(results.paramValues(1,1)*1000,results.mean_dxt,results.stdev_dxt,'k^')
            %paramValues(1,1) is the SLR rate in m/yr
            
            %results.mean_dxs
            %results.stdev_dxs
            
            %results.mean_dxt
            %results.stdev_dxt
            
            %          plot(results.norm_dxs_rates,'bo')
            %         hold on
            %         pause
            %
            %          plot(results.norm_dxt_rates,'k^')
            
            %         pause
            
            
        end
        
        clear results
    end
    
    %Fit a line through mean retreat rates
    shoreline_fit = polyfit(slr_rate,shoreline_mean,1);
    shoreface_fit = polyfit(slr_rate,shoreface_mean,1);
    
    %Calc r^2
    shoreline_rsq = r_squared(shoreline_mean,shoreline_fit(1)*slr_rate + shoreline_fit(2));
    shoreface_rsq = r_squared(shoreface_mean,shoreface_fit(1)*slr_rate + shoreface_fit(2));
    
    %Plot line
    plot(slr_rate,shoreline_fit(1)*slr_rate + shoreline_fit(2),'b')
    plot(slr_rate,shoreface_fit(1)*slr_rate + shoreface_fit(2),'k')
    
    text(2.5,1.75,['shoreline r^2 = ' num2str(shoreline_rsq)])
    text(2.5,1.5,['shoreface r^2 = ' num2str(shoreface_rsq)])
    
    title(title_list{ii})
    legend(['shoreline m=  ' num2str(shoreline_fit(1))],['shoreface m= ' num2str(shoreface_fit(1))],'Location','southeast')
	xlabel('SLR rate (mm/yr)')
    ylabel('Normalized (5mm/yr)Retreat Rate')
    set(gca,'ylim',[0 2],'xlim',[0 20])

end




%% clear a_star_eq_func from bstruct so the mat file stops giving a warning when it loads
%for n = 1:length(fdir)
    
%save([resultsPath '\' fdir(n).name],'b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')