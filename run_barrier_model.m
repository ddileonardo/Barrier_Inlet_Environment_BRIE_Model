%WaveTideRun3: 
%a0=0:0.1:1, hs = 0.5:0.1:1.5, Qow_max=50, Jmin=20km,grain_size=3e-4, marsh_cover=0.9, slr=2e-3, wave_period=7s w_b_crit=300m

%WaveTideRun4:
%a0=0:0.1:1, hs = 0.5:0.1:1.5, Qow_max=30, Jmin=20km,grain_size=3e-4, marsh_cover=0.8, slr=1e-3, wave_period=7s w_b_crit=200m

%WaveTideRun5:
%a0=0.5, hs = 0.5:0.2:1.5, Qow_max=30, Jmin=20km,grain_size=3e-4, marsh_cover=0.8, slr=0:5e-4:5e-3, wave_period=7s w_b_crit=200m

%WaveTideRun6:
%a0=0:0.1:1, hs = 0.5:0.1:1.5, Qow_max=25, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=2e-3, wave_period=10s w_b_crit=200m, (depth = 3?)

%WaveTideRun7: (drowning w inlet)
%a0=0.5, hs = 0.5:0.1:1.5, Qow_max=5, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=0:1e-3:1e-2, wave_period=10s w_b_crit=200m, (depth = 5?)

%WaveTideRun8:
%a0=0.5, hs = 1, Qow_max=5:5:50, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=2e-3, wave_period=10s w_b_crit=200m, (depth = 5)

%AsymHsRun
%a0=0.5, hs = 1, Qow_max=25, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=2e-3, wave_period=10s w_b_crit=200m, depth = 3m

%WaveTideRun9: (drowning wh inlet, 1d)
%a0=0.5, hs = 0.5:0.1:1.5, Qow_max=5, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=0:1e-3:1e-2, wave_period=10s w_b_crit=200m, (depth = 5?)

%WaveTideRun10: (drowning wh inlet)
%a0=0.5, hs = 0.5:0.1:1.5, Qow_max=5, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=0:1e-3:1e-2, wave_period=10s w_b_crit=200m, (depth = 5?)

%check mass balance:
%Vbarrier = b_struct.dy.*(sum(out.d_sf.*out.d_sf./out.s_sf_save./2.*(1-(b_struct.s_background./out.s_sf_save)),1)+sum(out.h_b_save.*double(out.x_b_save-out.x_s_save),1));

%drowning run no inlet
%a0=0.5, hs = 0.5:0.25:2, Qow_max=5, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=0:1e-3:1e-2, wave_period=10s w_b_crit=200m, (depth = 5?)

%drowning run inlet
%a0=0.5, hs = 0.5:0.1:1.5, Qow_max=5, Jmin=20km,grain_size=2e-4, %marsh_cover=0.5, slr=0:1e-3:1e-2, wave_period=10s w_b_crit=200m, (depth = 5?)

%grid test run and timestep run name = 'mass_balance_inlet2';
%param = {'dt','dy'}; param1 = [0.01 0.02 0.025 0.05 0.08 0.1 0.2 0.25]; param2 = [1000 800 500 400 250 100 80 50];

%% run_barrier_model
clear;clc;close all;
savepath = 'C:\Users\ddileonardo\The Water Institute of the Gulf\TO71 - Barrier Island Modeling - General\BRIE_Tests\SLR';
slr = 17; %mm/yr
name = ['Sea Level Experiments ' num2str(slr) ' mm per yr longshore on'];
param = {'slr';'grain_size';'wave_height';'wave_period';'h_b_crit';'ast_model_on';'inlet_model_on'}; 
paramValues = [slr*[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]; ...
    1.6e-4,1.6e-4,1.6e-4,1.6e-4,1.6e-4,1.6e-4,1.6e-4,1.6e-4,1.6e-4;...
    0.75,1,1.5,0.75,1,1.5,0.75,1,1.5;...
    8,8,10,8,8,10,8,8,10;...
    1,1,1,2,2,2,3,3,3;
    1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0];

% paramValues = [7e-3,7e-3,7e-3;...
%     1.6e-4,1.6e-4,1.6e-4;...
%     0.75,1.5,1.5;...
%     8,10,10;...
%     1,2,3;
%     1,1,1;
%     0,0,0];

%Set in initialize_barrier_model
dt = 0.05;
dtsave = 1e2;


%output = cell(length(paramValues(1,:)),1);
for ii=1:length(paramValues(1,:))%loop on different parameter sets
    %parfor jj=1:length(param2)
        %ii
        %jj
        
        b_struct = initialize_barrier_model; %initialize model
        b_struct.name = name; %set name
        
        for nn = 1:length(paramValues(:,1)) %loop on diff parameters of 1 set
            %Re-set parameters that you want to change
            b_struct.(param{nn}) = paramValues(nn,ii);
        end
        
        %run model 
        b_struct = barrier_model(b_struct);
        
        %Calculate retreat rates
        b_struct = shoreline_retreat_rate(b_struct,dt,dtsave);
        
        %save to output matrix
        %output(ii) = {b_struct};
        output.(['paramSet' num2str(ii)]) = b_struct;
    
    %end
end
b_struct = initialize_barrier_model;
save([savepath '\' name],'b_struct','output','param','paramValues','-v7.3')

close all
















