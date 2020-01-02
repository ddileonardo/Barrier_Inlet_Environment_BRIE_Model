
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
     weight = (obsWaves.hs_in(ii)^5).*(obsWaves.per_in^-5)*(sinh(
     
 
 end