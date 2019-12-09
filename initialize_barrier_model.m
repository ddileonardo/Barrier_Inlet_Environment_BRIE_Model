function b_struct = initialize_barrier_model
%#ok<*NASGU>
%which modules to run
name = 'ExampleBarrierPlot5'; 
barrier_model_on = true;
ast_model_on = true;
inlet_model_on = true;
make_gif = false;
plot_on = true;
sedstrat_on = false;

%general parameters
rho_w = 1025; %density water kg/m3
g = 9.81;

%wave climate parameters
wave_height = 1; %m
wave_period = 10; %s
wave_asym = 0.8;
wave_high = 0.2;

%barrier model parameters
slr = 2e-3; %sea level rise m/yr
s_background = 1e-3; %background slope (beta)
w_b_crit = 200; %critical barrier width m
h_b_crit = 2; %critical barrier height m
Qow_max = 20; %max overwash flux m3/m/yr

z = 10; %initial sea level (m)
bb_depth = 3; %back barrier depth (m)

grain_size = 2e-4; %(m)
R = 1.65; %relative density of sand

e_s = 0.01; %shoreface efficiency values from Aleja/Ashton
c_s = 0.01;

%alongshore stretch
dy = 100; %length of each alongshore section (m)
ny = 1000; %number of alonghsore sections

%timestepping
dt = 0.05; %years
nt = 1e5; %number of timesteps
dtsave = 1e3; %save spacing


%alongshore distribution of wave energy
wave_climl = 180; %resolution of possible wave approach angles (1 per degree, easy math)
AngArray = linspace(-0.5*pi,0.5*pi,wave_climl)';
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2; %correct k from Nienhuis, Ashton, Giosan 2015 (ashton's 2006 k value is wrong)

Jmin = 10000; %minimum inlet spacing, roos/schuttelaars (m)
a0 = 0.5; %amplitude of tide (m)
omega0 = 1.4e-4; %tidal frequency (rad s-1);
inlet_asp = sqrt(0.005); %aspect ratio inlet (gamma in swart/zimmerman)
%g = 9.81;
man_n = 0.05; %sm^-(1/3) manning n (vegetated channel)
u_e = 1; %ms-1 inlet equilibrium velocity (see swart/zimmerman)
inlet_max = 100; %maximum number of inlets in the model (mostly for debugging)
marsh_cover = 0.5; %percentage of backbarrier covered by marsh and therefore does not contribute to tidal prism

%new explicit relationship between boundary conditions and inlet area
u = @(a_star,gam,ah_star,a0) (sqrt(g.*a0).*sqrt(gam./2.*sqrt(a_star).*((-gam.*sqrt(a_star).*((a_star-ah_star).^2))+sqrt((gam.^2).*a_star.*((a_star-ah_star).^4)+4))));

%pretty function showing how the cross-sectional area varies with different back barrier configurations gamma
a_star_eq_fun = @(ah_star,gam,u_e_star) (real((2*ah_star)./3 + (2^(2/3)*((18*ah_star.*gam.^2 - 27*u_e_star.^4 - 2*ah_star.^3.*gam.^2.*u_e_star.^2 +...
    3*3^(1/2)*gam.^2.*u_e_star.^2.*(-(4*ah_star.^4.*gam.^4.*u_e_star.^4 - 4*ah_star.^3.*gam.^2.*u_e_star.^8 -...
    8*ah_star.^2.*gam.^4.*u_e_star.^2 + 36*ah_star.*gam.^2.*u_e_star.^6 + 4.*gam.^4 - 27.*u_e_star.^10)./...
    (gam.^4.*u_e_star.^6)).^(1/2))./(gam.^2.*u_e_star.^2)).^(1/3))./6 + (2^(1/3).*(ah_star.^2.*u_e_star.^2 + 3))./...
    (3.*u_e_star.^2.*((18*ah_star.*gam.^2 - 27.*u_e_star.^4 - 2*ah_star.^3.*gam.^2.*u_e_star.^2 + 3*3.^(1/2).*gam.^2.*u_e_star.^2.*...
    (-(4*ah_star.^4.*gam.^4.*u_e_star.^4 - 4.*ah_star.^3.*gam.^2.*u_e_star.^8 - 8.*ah_star.^2.*gam.^4.*u_e_star.^2 + 36*ah_star.*gam.^2.*u_e_star.^6 +...
    4.*gam.^4 - 27.*u_e_star.^10)./(gam.^4.*u_e_star.^6)).^(1/2))./(gam.^2.*u_e_star.^2)).^(1/3))));

%what are the inlet fractions (nienhuis, ashton 2016)
inlet_fraction = @(a,b,c,d,I) (a+(b./(1+c.*(I.^d))));

%1 = overtopping overwash, 2 = flood delta, 3 = dunes, 4: inlet fill
dz_strat = 0.02; %vertical resolution of the stratigraphy

%what is the dynamic equilibrium for overwash in this setting not considering inlets?
%Hde = h_b_crit-(slr*Vd_max/Qow_max);
%Wde = w_b_crit-(slr/s_background*Vd_max/Qow_max);
%Q_owde = slr/s_background*(Hde+bb_depth)+slr*Wde;
%s_sf_de = fsolve(@(s_sf_de) (((Hde+(s_background*(d_sf/s_background-Wde-d_sf/s_sf_de)))*2*(slr/s_background*(Hde+bb_depth)+slr*Wde)/(2*Hde+d_sf)/(Hde+bb_depth) - 4*k_sf*(s_sf_eq-s_sf_de)*(Hde+d_sf)/(2*Hde+d_sf)^2)-(slr/s_background)),s_sf_eq,optimoptions('fsolve','Display','none'));
%Qsf_de = k_sf*(s_sf_eq-s_sf_de);
%d_g = d_sf- d_sf*s_background/s_sf_de - s_background*Wde;
%ff = (d_g-bb_depth)/(d_g+Hde);

b_struct = v2struct([who; 'fieldNames']);