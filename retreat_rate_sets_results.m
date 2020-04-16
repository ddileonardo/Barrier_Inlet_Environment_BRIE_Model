%find the median retreat rates for the BRIE experiments
clc; close all; clear

%%
% resultsPath = 'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_125um';
% 
% fdir = dir([resultsPath '\*.mat']);

            
%% Rates for every SLR condition plotted as grid of subplots
% figure('position',[-1200         -75        1137     800],'color','w') %single external monitor to the left of laptop
% for n = 1:length(fdir)
%     
%     load([resultsPath '\' fdir(n).name]);
%     
%     paramSets = fieldnames(output); %fieldnames of each parameter set
%     
%     xt_rate_median = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
%     xt_rate_mean = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
%     xs_rate_median = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
%     xs_rate_mean = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
%     
%     for jj = 1:length(paramSets)
%         
%         
%         
%         for ii = 1:length(output.(paramSets{1}).dx_tdt(1,:))
%             
%             if length(output.(paramSets{jj}).dx_tdt(1,:)) == ii %if ii is the last timestep
%                 %continue
%             elseif sum(output.(paramSets{jj}).dx_tdt(:,ii+1)) == 0 %if all the rates are zero at the next time step
%                 %the barrier has drowned, break the loop
%                 break
%             end
%             
%             xt_rate_median(ii,jj) = double(median(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
%             xt_rate_mean(ii,jj) = mean(double(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
%             
%             
%             xs_rate_median(ii,jj) = double(median(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
%             xs_rate_mean(ii,jj) = mean(double(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
%             
%             
%         end
%         
%         %subplot(1,length(paramSets),jj)
%         subplot(3,3,jj)
%         plot(xt_rate_median(:,jj),'b.')
%         hold
%         plot(xs_rate_median(:,jj),'ro')
%         grid on
%         %plot(xt_rate_mean,'co')
%         %plot(xs_rate_mean,'mo')
%         
%         plot([10 10],[0 10],'k')
%         
%         xlabel('Timestep (5 yr intervals)')
%         ylabel('Retreat Rate (m/yr)')
%         
%         set(gca,'ylim',[0 15])
%         if jj == 1
%         legend('shoreface toe','shoreline','location','northwest')
%         end
%         title(['SLR: ' num2str(paramValues(1,jj)*1000) 'mmyr^-^1  H: ' num2str(paramValues(3,jj)) 'm  T: ' num2str(paramValues(4,jj)) 's  Hb,crit:' num2str(paramValues(5,jj)) 'm'])
%         
%        
%         
%     end
%     
%     
%      pause
%      img = getframe(gcf);
%      imwrite(img.cdata, [resultsPath '\Retreat_Rates_' fdir(n).name(1:end-4), '.png']);
%      %print(gcf, '-dpng','-r150',[resultsPath '\Retreat_Rates_' fdir(n).name(1:end-3) 'png'])
%      clf
%     
%     save([resultsPath '\' fdir(n).name],'b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median')
%     
% end

%% Normalize retreat rates by the baseline rate
%Find Mean and standard deviation of normalized rates

% baseline_name = 'Sea Level Experiments 9 mm per yr longshore on';
% baseline = load([resultsPath '\' baseline_name]); %load baseline model results
% 
% for n = 1:length(fdir)
%     
%     load([resultsPath '\' fdir(n).name]); %load individual SLR model experiment
%     
%     
%     paramSets = fieldnames(output); %fieldnames of each parameter set
%     
%     
%     if ~strcmp(fdir(n).name,baseline_name) %Don't normalizze the baseline because it will add a bunch of extra 1's to the mean
%         for jj = 1:length(paramSets)
%             norm_dxt_rates = []; %holds all the rates for a single SLR experimental run (SLR rate constant)
%             norm_dxs_rates = []; %holds all the rates for a single SLR experimental run (SLR rate constant)
% 
%             %make baseline data floating points so it can hold NaN values
%             baseline.output.(paramSets{jj}).dx_tdt = double(baseline.output.(paramSets{jj}).dx_tdt);
%             baseline.output.(paramSets{jj}).dx_sdt = double(baseline.output.(paramSets{jj}).dx_sdt);
%             
%             %NaN the 0 value retreat rates of the baseline
%             baseline.output.(paramSets{jj}).dx_tdt(baseline.output.(paramSets{jj}).dx_tdt<1) = NaN;
%             baseline.output.(paramSets{jj}).dx_sdt(baseline.output.(paramSets{jj}).dx_sdt<1) = NaN;
%             
%             %%Results
%             %Did barrier drown during this model run?
%             %if column sum == 0, barrier drowns at this timestep
%             drown_time = find(sum(output.(paramSets{jj}).dx_tdt)==0,1,'first');
%             %Make sure that you aren't looking at the very early part of
%             %the results
%             if find(sum(output.(paramSets{jj}).dx_tdt)>0,1,'last')>drown_time
%                 %if the last positive value occurs after the drown time,
%                 %then you have the wrong value
%                 drown_time =[];
%             end
%                         
%             %Remove negative and 0 values from results
%             result_dxt = double(output.(paramSets{jj}).dx_tdt);
%             result_dxt(result_dxt<1) = NaN;
%             result_dxs = double(output.(paramSets{jj}).dx_sdt);
%             result_dxs(result_dxs<1) = NaN;
%             
%             %remove last column of results if barrier drowns (super high
%             %retreat rates compared to rest of model run)
%             %NaN the timestep when the barrier drowns. (one before
%             %drown_time)
%             result_dxt(:,drown_time-1) = NaN;
%             result_dxs(:,drown_time-1) = NaN;
%             
%             %normalize by the base line
%             normalized.output.(paramSets{jj}).dx_tdt = result_dxt ./ baseline.output.(paramSets{jj}).dx_tdt;
%             normalized.output.(paramSets{jj}).dx_sdt = result_dxs ./ baseline.output.(paramSets{jj}).dx_sdt;
%             
%             %Gather values for taking the mean
%             norm_dxt_rates = [norm_dxt_rates; normalized.output.(paramSets{jj}).dx_tdt(:)];
%             norm_dxs_rates = [norm_dxs_rates; normalized.output.(paramSets{jj}).dx_sdt(:)];
%             
%             clear normalized result_dxt result_dxs
%             
%         end
%     
%     
%     %alculate the mean and standard deviation
%     mean_dxt = nanmean(norm_dxt_rates);
%     mean_dxs = nanmean(norm_dxs_rates);
%     
%     stdev_dxt = nanstd(norm_dxt_rates);
%     stdev_dxs = nanstd(norm_dxs_rates);
%     
%     save([resultsPath '\' fdir(n).name],'b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')
% 
%     clear('drown_time','b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')
%     
%     end
% end
% 
% 
% return
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