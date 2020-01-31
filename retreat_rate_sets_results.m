%find the median retreat rates for the BRIE experiments
clc; close all

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

%%

%Recalc medians and measn for each model output file
fdir = dir('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR');
%Load output files with unique names
%plot medians for each timestep for each parameter set for each model
%   output
%a single plot should have the median value for 1 parameter set (at every
%timestep) for each SLR 

%Put a line at 10 time steps (50 years)

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