%find the median retreat rates for the BRIE experiments
clc
for ii = 1:length(output)
    
    ii
    rate = output{ii,1}.dx_tdt(output{ii,1}.dx_tdt>0);
    fprintf(['shoreface toe rate is ' num2str(median(rate)) '\n'])
    
    rate = output{ii,1}.dx_sdt(output{ii,1}.dx_sdt>0);
    fprintf(['shoreline rate is ' num2str(median(rate)) '\n'])
    
    pause
    
end