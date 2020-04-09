%find the median retreat rates for the BRIE experiments
clc; close all

resultsPath = 'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR\GS_125um';

fdir = dir([resultsPath '\*.mat']);
figure('position',[-1200         -75        1137     800],'color','w') %single
            %external monitor to the left of laptop

for n = 1:length(fdir)
    
    load([resultsPath '\' fdir(n).name]);
    
    paramSets = fieldnames(output); %fieldnames of each parameter set
    
    xt_rate_median = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
    xt_rate_mean = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
    xs_rate_median = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
    xs_rate_mean = NaN*ones(length(output.(paramSets{1}).dx_tdt(1,:)),length(paramSets));
    
    for jj = 1:length(paramSets)
        
        
        
        for ii = 1:length(output.(paramSets{1}).dx_tdt(1,:))
            
            if length(output.(paramSets{jj}).dx_tdt(1,:)) == ii %if ii is the last timestep
                %continue
            elseif sum(output.(paramSets{jj}).dx_tdt(:,ii+1)) == 0 %if all the rates are zero at the next time step
                %the barrier has drowned, break the loop
                break
            end
            
            xt_rate_median(ii,jj) = double(median(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
            xt_rate_mean(ii,jj) = mean(double(output.(paramSets{jj}).dx_tdt(output.(paramSets{jj}).dx_tdt(:,ii)>0,ii)));
            
            
            xs_rate_median(ii,jj) = double(median(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
            xs_rate_mean(ii,jj) = mean(double(output.(paramSets{jj}).dx_sdt(output.(paramSets{jj}).dx_sdt(:,ii)>0,ii)));
            
            
        end
        
        %subplot(1,length(paramSets),jj)
        subplot(3,3,jj)
        plot(xt_rate_median(:,jj),'k.')
        hold
        plot(xs_rate_median(:,jj),'o','color',rgb('dodgerblue'))
        grid on
        %plot(xt_rate_mean,'co')
        %plot(xs_rate_mean,'mo')
        
        plot([10 10],[0 10],'k')
        
        xlabel('Timestep (5 yr intervals)')
        ylabel('Retreat Rate (m/yr)')
        
        set(gca,'ylim',[0 15])
        if jj == 1
        legend('shoreface toe','shoreline','location','northwest')
        end
        title(['SLR: ' num2str(paramValues(1,jj)*1000) 'mmyr^-^1  H: ' num2str(paramValues(3,jj)) 'm  T: ' num2str(paramValues(4,jj)) 's  Hb,crit:' num2str(paramValues(5,jj)) 'm'])
        
       
        
    end
    
    
     pause
     img = getframe(gcf);
     imwrite(img.cdata, [resultsPath '\Retreat_Rates_' fdir(n).name(1:end-4), '.png']);
     %print(gcf, '-dpng','-r150',[resultsPath '\Retreat_Rates_' fdir(n).name(1:end-3) 'png'])
     clf
    
    save([resultsPath '\' fdir(n).name])
    
end

return
