function output = shoreline_retreat_rate(output,dt,dtsave)
%Calculate shoreline retreat rate from the output of barrier_model
%dt and dtsave from initialize_barrier_model

time = dt*dtsave; %years; time step times the save time
%Shoreline retreat rate
output.dx_s = diff(output.x_s_save,1,2);
output.dx_sdt = output.dx_s/time;

%shoreface toe retreat rate
output.dx_t = diff(output.x_t_save,1,2);
output.dx_tdt = output.dx_t/time;

%back barrier retreat rate
output.dx_b = diff(output.x_b_save,1,2);
output.dx_bdt = output.dx_b/time;
