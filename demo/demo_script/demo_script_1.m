%% Loading demonstration data

myCurrent_script_file = matlab.desktop.editor.getActiveFilename; 
d = fileparts(myCurrent_script_file); 
load([d, filesep, '..', filesep, 'demo_data', filesep, 'demo_data_1']); 

%% performing a gridded-zero-padded reconstruction

% t is the trajectory

% ve is the list of volume elements

% y is the list of complex data values for each channel

% C is the coil-sensitivity estimation for each channel

% N_u is the size of the k-space gridd for regridding, which is here

% also the size of the reconstructed image. 

% dK_u is the step size of the k-space gridd. 

x = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);
bmImage(x); 