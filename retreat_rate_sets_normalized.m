%find the normalized retreat rates for the BRIE experiments
%normalized by background RSLR rate
clc; close all; clear

%%
resultsPath = 'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_125um';

fdir = dir([resultsPath '\*.mat']);

%% Normalize retreat rates by the baseline rate
%Find Mean and standard deviation of normalized rates

baseline_name = 'Sea Level Experiments 9 mm per yr longshore on';
baseline = load([resultsPath '\' baseline_name]); %load baseline model results

for n = 1:length(fdir)
    
    load([resultsPath '\' fdir(n).name]); %load individual SLR model experiment
    
    
    paramSets = fieldnames(output); %fieldnames of each parameter set
    
    
    if ~strcmp(fdir(n).name,baseline_name) %Don't normalizze the baseline because it will add a bunch of extra 1's to the mean
        for jj = 1:length(paramSets)
            norm_dxt_rates = []; %holds all the rates for a single SLR experimental run (SLR rate constant)
            norm_dxs_rates = []; %holds all the rates for a single SLR experimental run (SLR rate constant)

            %make baseline data floating points so it can hold NaN values
            baseline.output.(paramSets{jj}).dx_tdt = double(baseline.output.(paramSets{jj}).dx_tdt);
            baseline.output.(paramSets{jj}).dx_sdt = double(baseline.output.(paramSets{jj}).dx_sdt);
            
            %NaN the 0 value retreat rates of the baseline
            baseline.output.(paramSets{jj}).dx_tdt(baseline.output.(paramSets{jj}).dx_tdt<1) = NaN;
            baseline.output.(paramSets{jj}).dx_sdt(baseline.output.(paramSets{jj}).dx_sdt<1) = NaN;
            
            %%Results
            %Did barrier drown during this model run?
            %if column sum == 0, barrier drowns at this timestep
            drown_time = find(sum(output.(paramSets{jj}).dx_tdt)==0,1,'first');
            %Make sure that you aren't looking at the very early part of
            %the results
            if find(sum(output.(paramSets{jj}).dx_tdt)>0,1,'last')>drown_time
                %if the last positive value occurs after the drown time,
                %then you have the wrong value
                drown_time =[];
            end
                        
            %Remove negative and 0 values from results
            result_dxt = double(output.(paramSets{jj}).dx_tdt);
            result_dxt(result_dxt<1) = NaN;
            result_dxs = double(output.(paramSets{jj}).dx_sdt);
            result_dxs(result_dxs<1) = NaN;
            
            %remove last column of results if barrier drowns (super high
            %retreat rates compared to rest of model run)
            %NaN the timestep when the barrier drowns. (one before
            %drown_time)
            result_dxt(:,drown_time-1) = NaN;
            result_dxs(:,drown_time-1) = NaN;
            
            %normalize by the base line
            normalized.output.(paramSets{jj}).dx_tdt = result_dxt ./ baseline.output.(paramSets{jj}).dx_tdt;
            normalized.output.(paramSets{jj}).dx_sdt = result_dxs ./ baseline.output.(paramSets{jj}).dx_sdt;
            
            %Gather values for taking the mean
            norm_dxt_rates = [norm_dxt_rates; normalized.output.(paramSets{jj}).dx_tdt(:)];
            norm_dxs_rates = [norm_dxs_rates; normalized.output.(paramSets{jj}).dx_sdt(:)];
            
            clear normalized result_dxt result_dxs
            
        end
    
    
    %alculate the mean and standard deviation
    mean_dxt = nanmean(norm_dxt_rates);
    mean_dxs = nanmean(norm_dxs_rates);
    
    stdev_dxt = nanstd(norm_dxt_rates);
    stdev_dxs = nanstd(norm_dxs_rates);
    
    save([resultsPath '\' fdir(n).name],'b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')

    clear('drown_time','b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')
    
    end
end





%% clear a_star_eq_func from bstruct so the mat file stops giving a warning when it loads
%for n = 1:length(fdir)
    
%save([resultsPath '\' fdir(n).name],'b_struct','output','param','paramSets','paramValues','xs_rate_mean','xs_rate_median','xt_rate_mean','xt_rate_median','norm_dxt_rates','norm_dxs_rates','mean_dxt','mean_dxs','stdev_dxt','stdev_dxs')