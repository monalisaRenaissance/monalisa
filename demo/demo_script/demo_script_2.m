%% Initialisation

% Loading demonstration data
myCurrent_script_file = matlab.desktop.editor.getActiveFilename; 
d = fileparts(myCurrent_script_file); 
load([d, filesep, '..', filesep, 'demo_data', filesep, 'demo_data_2']);  


% Creating a folder for reconstruciton result and seting it as current
% directory. 
cd([d, filesep, '..', filesep, '..', filesep, '..']); 
bmCreateDir('recon_folder'); 
cd('recon_folder')

%% initial image

x0 = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);
bmImage(x0); 


%% evaluation of gridding matrices

[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);

%% Iterative-SENSE (least-square) reconstruction



file_label   = 'sensa'; 
frSize      = N_u; 
nCGD        = 4;
ve_max      = 5*prod(dK_u(:)); 
nIter       = 20;
witness_ind = 1:5:nIter; 


witnessInfo = bmWitnessInfo(file_label, witness_ind); 


x       = bmSensa(x0, y, ve, C, Gu, Gut, frSize, ...
            nCGD, ve_max, ...
            nIter, witnessInfo); 

bmImage(x); 


%% Iterative least-square regularized reconstruction 
%  with l1 norm of spatial derivative regularization

file_label   = 'steva'; 
frSize      = N_u; 
nCGD        = 4;
delta       = 0.2; 
rho         = 10*delta; 
ve_max      = 5*prod(dK_u(:)); 
nIter       = 20;
witness_ind = 1:5:nIter; 


witnessInfo = bmWitnessInfo(file_label, witness_ind); 


x       = bmSteva(x0, [], [], y, ve, C, Gu, Gut, frSize, ...
            delta, rho, nCGD, ve_max, ...
            nIter, witnessInfo); 

bmImage(x); 


%% Iterative least-square regularized reconstruction 
%  with l2 norm of spatial derivative regularization


file_label   = 'sensitiva_spatialD'; 
frSize      = N_u; 
nCGD        = 4;
delta       = 2;
regularizer = 'spatial_derivative'; % 'identity' or 'spatial_derivative'
ve_max      = 5*prod(dK_u(:)); 
nIter       = 20;
witness_ind = 1:5:nIter; 


witnessInfo = bmWitnessInfo(file_label, witness_ind); 


x = bmSensitiva(x0, y, ve, C, Gu, Gut, frSize, ...
                delta, nCGD, ve_max, ...
                regularizer, ...
                nIter, witnessInfo); 

bmImage(x); 

%% Iterative least-square regularized reconstruction 
%  with l2 norm of image regularization


file_label   = 'sensitiva_id'; 
frSize      = N_u; 
nCGD        = 4;
delta       = 2;
regularizer = 'identity'; % 'identity' or 'spatial_derivative'
ve_max      = 5*prod(dK_u(:)); 
nIter       = 20;
witness_ind = 1:5:nIter; 


witnessInfo = bmWitnessInfo(file_label, witness_ind); 


x = bmSensitiva(x0, y, ve, C, Gu, Gut, frSize, ...
                delta, nCGD, ve_max, ...
                regularizer, ...
                nIter, witnessInfo); 

bmImage(x); 