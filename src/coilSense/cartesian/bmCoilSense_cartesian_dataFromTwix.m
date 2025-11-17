% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

function y = bmCoilSense_cartesian_dataFromTwix(argFile, N_u)

myTwix      = bmTwix(argFile, false, 0, true); 
t_info      = bmTwix_cartesian_trajInfo(myTwix); 
N_u_raw     = t_info.N_u; 
dK_u_raw    = t_info.dK_u; 

y   = bmTwix_cartesian_data(myTwix); 
y   = reshape(y, [size(y, 1), N_u_raw]); 
y   = bmLowRes_cartesian(y, N_u_raw, dK_u_raw, N_u, dK_u_raw); 
y   = reshape(y, [size(y, 1), N_u]);
y   = bmPermuteToCol(y); 
y   = bmBlockReshape(y, N_u); 

end