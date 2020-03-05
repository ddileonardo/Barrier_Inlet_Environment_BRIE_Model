%find the median retreat rates for the BRIE experiments
clc; close all

fdir = dir('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\*.mat');

for n = 1:length(fdir)
    
    load(['C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\' fdir(n).name]);
    
    paramSets = fieldnames(output); %fieldnames of each parameter set
    
    for jj = 1:length(paramSets)
        
        xt_rate_median = NaN*ones(size(output.(paramSets{1}).dx_tdt(1,:)));
        xt_rate_mean = NaN*ones(size(output.(paramSets{1}).dx_tdt(1,:)));
        xs_rate_median = NaN*ones(size(output.(paramSets{1}).dx_tdt(1,:)));
        xs_rate_mean = NaN*ones(size(output.(paramSets{1}).dx_tdt(1,:)));
        
        for ii = 1:length(output.(paramSets{1}).dx_tdt(1,:))
            
            if sum(output.(paramSets{jj}).dx_tdt(:,ii+1)) == 0 %if all the rates are zero at the next time step
                %the barrier has drowned, break the loop
                break
            end
            
            xt_rate_median(ii) = double(median(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
            xt_rate_mean(ii) = mean(double(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
            
            
            xs_rate_median(ii) = double(median(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
            xs_rate_mean(ii) = mean(double(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
            
            
        end
        
        plot(xt_rate_median,'b.')
        hold
        %plot(xt_rate_mean,'co')
        
        
        plot(xs_rate_median,'ro')
        %plot(xs_rate_mean,'mo')
        
        plot([10 10],[0 10],'k')
        
        xlabel('Timestep (5 yr intervals)')
        ylabel('Retreat Rate (m/yr)')
        
        set(gca,'ylim',[0 15])
        legend('shoreface toe','shoreline','location','southeast')
        title(['SLR: ' num2str(paramValues(1,jj)*1000) 'mmyr^-^1  H: ' num2str(paramValues(3,jj)) 'm  T: ' num2str(paramValues(4,jj)) 's  Hb,crit:' num2str(paramValues(5,jj)) 'm'])
        
        pause
        clf
        
    end
    
    save(['C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\' fdir(n).name])
    
end

%% Plot results all together

%Load output files with unique names
SLR7 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 7 mm per yr longshore on.mat');
SLR9 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 9 mm per yr longshore on.mat');
SLR11 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 11 mm per yr longshore on.mat');
SLR13 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 13 mm per yr longshore on.mat');
SLR15 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 15 mm per yr longshore on.mat');
SLR17 = load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\Sea Level Experiments 17 mm per yr longshore on.mat');

%plot medians for each timestep for each parameter set for each model
%   output
%a single plot should have the median value for 1 parameter set (at every
%timestep) for each SLR 
%Put a line at 10 time steps (50 years)
%%
h1 = figure; 
h2 = figure;
for m = 1:length(SLR9.paramValues(1,:)) %loop on columns; loop on different parameter sets for same SLR rate
    
    figure(h1)
    subplot(1,2,m)
    hold on
    
    plot(SLR7.xs_rate_median,'o','color',rgb('blue'))
    plot(SLR9.xs_rate_median,'o','color',rgb('darkred'))
    plot(SLR11.xs_rate_median,'s','color',rgb('gold'))
    plot(SLR13.xs_rate_median,'s','color',rgb('black'))
    plot(SLR15.xs_rate_median,'s','color',rgb('lime'))
    plot(SLR17.xs_rate_median,'o','color',rgb('cyan'))
    plot([10 10],[0 10],'k')
    
    set(gca,'ylim',[0 10])
    
    title(['Median Xs retreat  H: ' num2str(SLR9.paramValues(3,m)) 'm  T: ' num2str(SLR9.paramValues(4,m)) 's  Hb,crit:' num2str(SLR9.paramValues(5,m)) 'm'])

    if m == 1
    legend('7 mm/yr','9 mm/yr', '11 mm/yr','13 mm/yr','15 mm/yr','17 mm/yr','location','southeast')
    end
    
    %if m == 7 || m == 8
    xlabel('Timestep (5 yr intervals)')
    
    %end
    
    %if m == 3 || m == 7
        ylabel('Retreat Rate (m/yr)')
    %end
    
    %%%%%%%%%%%%%%%%%%%%
    figure(h2)
    subplot(1,2,m)
    hold on
    
    plot(SLR7.xt_rate_median,'o','color',rgb('blue'))
    plot(SLR9.xt_rate_median,'o','color',rgb('darkred'))
    plot(SLR11.xt_rate_median,'s','color',rgb('gold'))
    plot(SLR13.xt_rate_median,'s','color',rgb('black'))
    plot(SLR15.xt_rate_median,'s','color',rgb('lime'))
    plot(SLR17.xt_rate_median,'o','color',rgb('cyan'))
    plot([10 10],[0 10],'k')
    
    set(gca,'ylim',[0 10])
    
    title(['Median Xt retreat  H: ' num2str(SLR9.paramValues(3,m)) 'm  T: ' num2str(SLR9.paramValues(4,m)) 's  Hb,crit:' num2str(SLR9.paramValues(5,m)) 'm'])

    if m == 1
    legend('7 mm/yr','9 mm/yr', '11 mm/yr','13 mm/yr','15 mm/yr','17 mm/yr','location','southeast')
    end
    
    %if m == 7 || m == 8
    xlabel('Timestep (5 yr intervals)')
    
   % end
    
    %if m == 3 || m == 7
        ylabel('Retreat Rate (m/yr)')
    %end
    
    
end




%% for output format with cells
%clc; close all
% for jj = 1:length(output)
%     
%         
%         xt_rate_median(jj) = double(median(output{jj,1}.dx_tdt(output{jj,1}.dx_tdt>0)));
%         xt_rate_mean(jj) = mean(double(output{jj,1}.dx_tdt(output{jj,1}.dx_tdt>0)));
%         %fprintf(['shoreface toe rate is ' num2str(median(rate)) '\n'])
%         
%         xs_rate_median(jj) = double(median(output{jj,1}.dx_sdt(output{jj,1}.dx_sdt>0)));
%         xs_rate_mean(jj) = mean(double(output{jj,1}.dx_sdt(output{jj,1}.dx_sdt>0)));
%         %fprintf(['shoreline rate is ' num2str(median(rate)) '\n'])
%         
%         %     for r = 1:length(output{jj,1}.dx_sdt(1,:))
%         %     s_sf(r) = nanmean(output{jj,1}.d_sf./(double(output{jj,1}.x_s_save(output{jj,1}.dx_sdt(:,r)>0,r))- double(output{jj,1}.x_t_save(output{jj,1}.dx_sdt(:,r)>0,r))));
%         %
%         %     end
%         %fprintf(['shoreface slope is ' num2str(nanmedian(s_sf(:))) '\n'])
%         
%         
% 
%     
%     
% end
% 
%         plot(xt_rate_median,'bo')
%         hold on
%         plot(xt_rate_mean,'co')
%         
%         
%         plot(xs_rate_median,'ro')
%         plot(xs_rate_mean,'mo')