%% Initialisation

% Loading demonstration data
myCurrent_script_file = matlab.desktop.editor.getActiveFilename; 
d = fileparts(myCurrent_script_file); 
load([d, filesep, '..', filesep, 'demo_data', filesep, 'demo_data_3']);  


% Creating a folder for reconstruciton result and seting it as current
% directory. 
cd([d, filesep, '..', filesep, '..', filesep, '..']); 
bmCreateDir('recon_folder'); 
cd('recon_folder')


%% resizing C, computing ve and permuting y to correct shape

% resizing C
C       = bmImResize(C, [48, 48], N_u);


% computing ve for 2D radial
ve = cell(nFr, 1); 
ve_max  = 5*prod(dK_u(:)); 
for i = 1:nFr
    ve{i} = bmVolumeElement_voronoi_full_radial2(t{i});
    ve{i} = min(ve{i}, ve_max); 
end


% putting data in column shape for reconstruction
for i = 1:nFr
    y{i} = bmPermuteToCol(y{i}); 
end


%% initial image

x0 = cell(nFr, 1); 
for i = 1:nFr
    x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, N_u, dK_u);
end

bmImage(x0); 


%% evaluation of gridding matrices

[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);


%% TevaMorphosia_chain without deformation fields

file_label      = 'tevaMorphosia_chain'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 0.5; 
rho             = 10*delta; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmTevaMorphosia_chain(x0, ...
                          [], [], ... 
                          y, ve, C, ...
                          Gu, Gut, frSize, ...
                          [], [], ...
                          delta, rho, regul_mode, ...
                          nCGD, ve_max, ...
                          nIter, witnessInfo); 


bmImage(x)


%% deform_field evaluation with imReg Demon

maxPixDisplacement  = 10; 
nIter               = 500; 
reg_file            = []; 
reg_mask            = []; 

h = x; 
[DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 

%% Calculating deformation matrices

[Tu, Tut] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);

%% TevaMorphosia_chain with deformation fields

file_label      = 'tevaMorphosia_chain_DF'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 1.5; 
rho             = 10*delta; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmTevaMorphosia_chain(x0, ...
                          [], [], ... 
                          y, ve, C, ...
                          Gu, Gut, frSize, ...
                          Tu, Tut, ...
                          delta, rho, regul_mode, ...
                          nCGD, ve_max, ...
                          nIter, witnessInfo); 


bmImage(x)


%% TevaDuoMorphosia_chain without deformation fields



file_label      = 'tevaDuoMorphosia_chain'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 0.5; 
rho             = 10*delta; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmTevaDuoMorphosia_chain(x0, ...
                             [], [], [], [], ... 
                             y, ve, C, ...
                             Gu, Gut, frSize, ...
                             [], [], [], [], ...
                             delta, rho, regul_mode, ...
                             nCGD, ve_max, ...
                             nIter, witnessInfo); 


bmImage(x)


%% deform_field evaluation with imReg Demon

maxPixDisplacement  = 10; 
nIter               = 500; 
reg_file            = []; 
reg_mask            = []; 

h = x; 
[DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 
[DF_to_next, imReg_to_next]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_next', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


%% Calculating deformation matrices

[Tu1, Tu1t] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);
[Tu2, Tu2t] = bmImDeformField2SparseMat(DF_to_next, N_u, [], true);


%% TevaDuoMorphosia_chain with deformation fields

file_label      = 'tevaDuoMorphosia_chain_DF'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 1.5; 
rho             = 10*delta; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmTevaDuoMorphosia_chain(x0, ...
                             [], [], [], [], ... 
                             y, ve, C, ...
                             Gu, Gut, frSize, ...
                             Tu1, Tu1t, Tu2, Tu2t, ...
                             delta, rho, regul_mode, ...
                             nCGD, ve_max, ...
                             nIter, witnessInfo); 


bmImage(x)



%% SensitivaMorphosia_chain without deformation fields

file_label      = 'sensitivaMorphosia_chain'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 5; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmSensitivaMorphosia_chain(  x0, ...
                                 y, ve, C, ...
                                 Gu, Gut, frSize, ...
                                 [], [], ...
                                 delta, regul_mode, ...
                                 nCGD, ve_max, ...
                                 nIter, witnessInfo  ); 


bmImage(x)


%% deform_field evaluation with imReg Demon

maxPixDisplacement  = 10; 
nIter               = 500; 
reg_file            = []; 
reg_mask            = []; 

h = x; 
[DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 

%% Calculating deformation matrices

[Tu, Tut] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);



%% SensitivaMorphosia_chain with deformation fields

file_label      = 'sensitivaMorphosia_chain_DF'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 15; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmSensitivaMorphosia_chain(  x0, ...
                                 y, ve, C, ...
                                 Gu, Gut, frSize, ...
                                 Tu, Tut, ...
                                 delta, regul_mode, ...
                                 nCGD, ve_max, ...
                                 nIter, witnessInfo  ); 


bmImage(x)



%% SensitivaDuoMorphosia_chain without deformation fields

file_label      = 'sensitivaDuoMorphosia_chain'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 5; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmSensitivaDuoMorphosia_chain(  x0, ...
                                    y, ve, C, ...
                                    Gu, Gut, frSize, ...
                                    [], [], [], [], ...
                                    delta, regul_mode, ...
                                    nCGD, ve_max, ...
                                    nIter, witnessInfo  ); 


bmImage(x)

%% deform_field evaluation with imReg Demon

maxPixDisplacement  = 10; 
nIter               = 500; 
reg_file            = []; 
reg_mask            = []; 

h = x; 
[DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 
[DF_to_next, imReg_to_next]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_next', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


%% Calculating deformation matrices

[Tu1, Tu1t] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);
[Tu2, Tu2t] = bmImDeformField2SparseMat(DF_to_next, N_u, [], true);


%% SensitivaDuoMorphosia_chain with deformation fields

file_label      = 'sensitivaDuoMorphosia_chain_DF'; % label to name the file containing the result
frSize          = N_u; 
nCGD            = 4;
nIter           = 30;

delta           = 15; 

regul_mode      = 'normal'; 

witness_ind = 1:5:nIter; 
witnessInfo = bmWitnessInfo(file_label, witness_ind); 



x = bmSensitivaDuoMorphosia_chain(  x0, ...
                                    y, ve, C, ...
                                    Gu, Gut, frSize, ...
                                    Tu1, Tu1t, Tu2, Tu2t, ...
                                    delta, regul_mode, ...
                                    nCGD, ve_max, ...
                                    nIter, witnessInfo  ); 


bmImage(x)