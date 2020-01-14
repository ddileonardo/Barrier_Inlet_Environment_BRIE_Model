%find the median retreat rates for the BRIE experiments
clc
for ii = 1:length(output)
    
    ii
    %     rate = output{ii,1}.dx_tdt(output{ii,1}.dx_tdt>0);
    %     fprintf(['shoreface toe rate is ' num2str(median(rate)) '\n'])
    %
    %     rate = output{ii,1}.dx_sdt(output{ii,1}.dx_sdt>0);
    %     fprintf(['shoreline rate is ' num2str(median(rate)) '\n'])
    %
    for jj = 1:length(output{ii,1}.dx_sdt(1,:))
    s_sf(jj) = nanmean(output{ii,1}.d_sf./(double(output{ii,1}.x_s_save(output{ii,1}.dx_sdt(:,jj)>0,jj))- double(output{ii,1}.x_t_save(output{ii,1}.dx_sdt(:,jj)>0,jj))));
    
    end
    %fprintf(['shoreface slope is ' num2str(nanmedian(s_sf(:))) '\n'])
    plot(s_sf)
    
    
    %     width = double(output{ii,1}.x_b_save(output{ii,1}.dx_sdt>0)) - double(output{ii,1}.x_s_save(output{ii,1}.dx_sdt>0));
    %     plot(width)
    %
    pause
    
end