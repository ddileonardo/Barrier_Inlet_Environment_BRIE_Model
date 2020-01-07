
wave_asym = 0.8;
wave_high = 0.2;
wave_height = 1.0;
dt = 0.05;
wave_period = 10;
g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2; %correct k from Nienhuis, Ashton, Giosan 2015 (ashton's 2006 k value is wrong)


%wave direction
wave_ang = find(wave_cdf>rand,1);

%alongshore distribution of wave energy
wave_climl = 180; %resolution of possible wave approach angles (1 per degree, easy math)
AngArray = linspace(-0.5*pi,0.5*pi,wave_climl)';
theta = 60;%made up

wave_pdf = 4*[wave_asym*wave_high*ones(wave_climl/4,1);...
    
           wave_asym*(1-wave_high)*ones(wave_climl/4,1);...
           
            (1-wave_asym)*(1-wave_high)*ones(wave_climl/4,1);...
            
            (1-wave_asym)*wave_high*ones(wave_climl/4,1)]./wave_climl;


wave_cdf = cumsum(wave_pdf);

coast_qs = wave_height.^2.4.*(wave_period.^0.2)*3600*365*24*k .*  (cos(AngArray).^1.2) .* sin(AngArray); %m3/yr



 Qs = dt.*coast_qs(min(wave_climl,max(1,round(wave_climl-wave_ang-(wave_climl./180.*theta)+1))));
 
 %% Ortiz and Ashton 2016
 
 %morphodynamically representative wave conditions
 
 %REDO STARTING WITH THE Hs_strt, etc
 
 for ii = 1:length(obsWaves.hs_in)
     %H^5*T^-5*sinh(k*z)^-5; where k is the wavenumber
     heightPart = obsWaves.hs_in(ii)/length(obsWaves.hs_in);
     %weight = (obsWaves.hs_in(ii)^5).*(obsWaves.per_in^-5)*(sinh(
     
 
 end
 
 %% 
%  1.	Shoreface depth – suggest aggregating all of the Beasley data to get 
%  retreat rates as a function of depth for all of the profiles and regions 
%  over all times, then plotting a histogram. I think that may make it possible 
%  to determine a reasonable shoreface depth in the sense that BRIE is using it. 
%  It doesn’t have to be perfect, but by comparing to your calculated values that
%  should help figure out what’s going on with your question (1) and what values 
%  we need to use to be realistic for this coast. (Eventually we may want to look 
%  at this regionally, but for now I think lumping everything together makes sense).
%  
%%%%%%
%Median retreat rate is between 3 and 10 m/yr
%%%%%%
%load('C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\Beasley\beasley.mat')

for ii = 1:12
    figure
    fieldname = ['isobath' num2str(ii)];
    data = beasley.(fieldname);
    
    data(data==0) = NaN;
    histogram(data,40,'BinLimits',[-100 100])
    set(gca,'ylim',[0,1400])
    %plot(data,'.')
    text(-75,1350,['mean = ' num2str(nanmean(data))])
    text(-75,1250,['median = ' num2str(nanmedian(data))])
    text(-75,1150,['mode = ' num2str(mode(round(data)))])
    title(['isobath ' num2str(-ii)])
    xlabel('Retreat Rate m/yr')
    ylabel('bin count')
    
    pause
    
    print(gcf,'-dpng','-r100',['C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\Beasley\retreat_rate_plots\' fieldname '_clipped'])
    
    close
    
end


%%
% 2.	Sea level rise rate. I think the value you used makes sense, but I’m 
% wondering if we can make use of the real (space/time-varying) information from 
% the local tide gauges to compare to real time/space shoreface retreat rates 
% to get a handle on the BRIE model (both how it’s working in practice and cal/val). 
% Ioannis’ input is appreciated here, especially with much more familiarity in
% this region, but I’d suggest pulling monthly time-series from a few tide gauges 
% in the area and plotting them, marking the times of Ben’s data on the plots, so
% we can look at the variation in time/space and compare to the spatial/time variation 
% in cross-shore retreat rates. There’s also the storm impacts that Ben et al. looked 
% at that aren’t explicitly/individually in BRIE, but by looking at it all together 
% I think there’s a chance we might be able to identify if there’s enough spatial/temporal 
% variability somewhere to use in calibrating BRIE. (There’s nonlinear dependencies, 
% but shoreface retreat rate scales so heavily as a function of SLR in the model that 
% I’m hoping we can see some real world variation in Ben’s cross-shore retreat 
% data that we can map back to variation in RSLR in space and/or time).

%%%%%%%
%Rates range from 4.8 mm/yr to 38 mm/yr; The most reasonable rates seem to
%be the ones that NOAA has used to calculated sea level rise - New Canal, 
%Bay Waveland, and GRand Isle. The Eugene Island is that NOAA used is very old.
%There is newer data at Eugene Island, but the record is very short.
%The rates from NOAA are: 
%Grand Isle 9.08 mm/yr
%New Canal 5.35 mm/yr
%Bay Waveland 4.64 mm/yr
%The rates calculated below are slightly different, but close.
%
%The rates from other stations are:
%Amerada Pass 18.2 mm/yr
%Belle Pass 15.9 mm/yr
%Shell Beach 13.5 mm/yr
%Pilots Stations East 38.0 mm/yr
%%%%%%%

convert = 1000*365.25; %conversion factor for m/day to mm/yr  

figure
plot(AmeradaPass.Date, AmeradaPass.MSLm,'o')
hold on
AmeradaPass_line = polyfit(AmeradaPass.Date,AmeradaPass.MSLm,1);
plot(AmeradaPass.Date,AmeradaPass.Date*AmeradaPass_line(1) + AmeradaPass_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(AmeradaPass_line(1)*convert)])
title('Amerada Pass')

figure
plot(EugeneIsland.Date,EugeneIsland.MSLm,'o')
hold on
EugeneIsland_line = polyfit(EugeneIsland.Date(~isnan(EugeneIsland.MSLm)),EugeneIsland.MSLm(~isnan(EugeneIsland.MSLm)),1);
plot(EugeneIsland.Date,EugeneIsland.Date*EugeneIsland_line(1) + EugeneIsland_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(EugeneIsland_line(1)*convert)])
title('Eugene Island')

figure
plot(BellePass.Date,BellePass.MSLm,'o')
hold on
BellePass_line = polyfit(BellePass.Date(~isnan(BellePass.MSLm)),BellePass.MSLm(~isnan(BellePass.MSLm)),1);
plot(BellePass.Date,BellePass.Date*BellePass_line(1) + BellePass_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(BellePass_line(1)*convert)])
title('Belle Pass')

figure
plot(GrandIsle.Date,GrandIsle.MSLm,'o')
hold on
GrandIsle_line = polyfit(GrandIsle.Date(~isnan(GrandIsle.MSLm)),GrandIsle.MSLm(~isnan(GrandIsle.MSLm)),1);
plot(GrandIsle.Date,GrandIsle.Date*GrandIsle_line(1) + GrandIsle_line(2))
[x,y] = ginput(1);
text(x,y,[' m = ' num2str(GrandIsle_line(1)*convert)])
title('Grand Isle')

figure
plot(BayWaveland.Date,BayWaveland.MSLm,'o')
hold on
BayWaveland_line = polyfit(BayWaveland.Date(~isnan(BayWaveland.MSLm)),BayWaveland.MSLm(~isnan(BayWaveland.MSLm)),1);
plot(BayWaveland.Date,BayWaveland.Date*BayWaveland_line(1) + BayWaveland_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(BayWaveland_line(1)*convert)])
title('Bay Waveland')

figure
plot(NewCanal.Date,NewCanal.MSLm,'o')
hold on
NewCanal_line = polyfit(NewCanal.Date(~isnan(NewCanal.MSLm)),NewCanal.MSLm(~isnan(NewCanal.MSLm)),1);
plot(NewCanal.Date,NewCanal.Date*NewCanal_line(1) + NewCanal_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(NewCanal_line(1)*convert)])
title('New Canal')

figure
plot(PilotsStationEast.Date,PilotsStationEast.MSLm,'o')
hold on
PilotsStationEast_line = polyfit(PilotsStationEast.Date(~isnan(PilotsStationEast.MSLm)),PilotsStationEast.MSLm(~isnan(PilotsStationEast.MSLm)),1);
plot(PilotsStationEast.Date,PilotsStationEast.Date*PilotsStationEast_line(1) + PilotsStationEast_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(PilotsStationEast_line(1)*convert)])
title('Pilots Station East')

figure
plot(ShellBeach.Date,ShellBeach.MSLm,'o')
hold on
ShellBeach_line = polyfit(ShellBeach.Date(~isnan(ShellBeach.MSLm)),ShellBeach.MSLm(~isnan(ShellBeach.MSLm)),1);
plot(ShellBeach.Date,ShellBeach.Date*ShellBeach_line(1) + ShellBeach_line(2))
[x,y] = ginput(1);
text(x,y,['m = ' num2str(ShellBeach_line(1)*convert)])
title('Shell Beach')


AmeradaPass_rate = AmeradaPass_line(1)*convert;
EugeneIsland_rate = EugeneIsland_line(1)*convert;
BellePass_rate = BellePass_line(1)*convert;
GrandeIsle_rate = GrandIsle_line(1)*convert;
BayWaveland_rate = BayWaveland_line(1)*convert;
NewCanal_rate = NewCanal_line(1)*convert;
PilotsStationEast_rate = PilotsStationEast_line(1)*convert;
ShellBeach_rate = ShellBeach_line(1)*convert;

%%
% 3.	Shoreface retreat rates themselves – how do they compare to what Ben is 
% seeing? Again, spatial and temporal variation might be our friend here to 
% figure out what parameters in the model may need to be tweaked and how to 
% address the discrepancies you are seeing between paper and model. Eventually 
% we will probably want to reach out to Jaap or Jorge about it, but suggest 
% we want to have a little better handle on the comparison of the model to 
% real-world first.

%%%%%%%
% With the shoreface depth at 8.9*wave height, the shoreline retreat rates are all low (< 10 m/yr)
% With the shoreface depth at Hallemeier or 14 m, the shoreline retreat rates are all very high:
%   hundreds of m/yr for Hallemeier – shoreface depths are 10 m to 43 m
%   ~150 m/yr for a 14 m shoreface depth 
%
%Shoreface retreat rates from Beasley have a median betwen 4 and 10 m/yr.
%When aggregated into 
%%%%%%%

























 