function output = shoreline_retreat_rate(output,dt,dtsave)
%Calculate shoreline retreat rate from the output of barrier_model
%dt and dtsave from initialize_barrier_model

time = dt*dtsave; %years; time step times the save time
output.dx = diff(output.x_s_save,1,2);
output.dxdt = output.dx/time;