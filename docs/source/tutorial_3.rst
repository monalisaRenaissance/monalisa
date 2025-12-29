Tutorial 3 : Iterative CINE reconstructions with 1 temporal dimension for non-cartesian data
============================================================================================

The code of this tutorial is written in the file `demo_script_3.m` in the demonstration folder. 
As any demonstration script, you can open it and run it right away. 

We describe here how to call the iterative CINE reconstructions of Monalisa for images with 1 temporal dimension 
for non-cartesian data. They are called **TevaMorphosia_chain**, **TevaDuoMorphosia_chain**, **SensitivaMorphosia_chain** and
**SensitivaDuoMorphosia_chain**. 

In contrast to spatially regularized reconstruction, the present temporally regularized reconstructions are quite efficient 
as demonstrated in publication ???cite publications???. They could be realistic candidates to be tested for clinical use (yet keep in mind
that the content of Monalisa cannot be used for any diagnostic purpose as it is an experimental project only).  

The script contains the following section. 

1. **Initialisation**

    The initialization section loads the demonstration data and creates a folder to store the reconstruction results. 

    .. code-block:: matlab  

        % Loading demonstration data
        myCurrent_script_file = matlab.desktop.editor.getActiveFilename; 
        d = fileparts(myCurrent_script_file); 
        load([d, filesep, '..', filesep, 'demo_data', filesep, 'demo_data_3']);  


        % Creating a folder for reconstruciton result and seting it as current
        % directory. 
        cd([d, filesep, '..', filesep, '..', filesep, '..']); 
        bmCreateDir('recon_folder'); 
        cd('recon_folder')

    Among other things, the loaded data contains a list of trajectory `t`. There is one trajectory for each frame.  
    You can plot the sampling trajectory `t` for each frame to have an idea of how good the data are sampled. 
    You may get someting as follows:      

    .. image:: ./images/data_3_t_partialSamp.gif
        :width: 60 %
        :align: center

    .. raw:: html

        <div style="margin-bottom: 30px;"></div>

    This list of trajectories contains in average 15% of the 256 original lines, still with 512 points each. It is
    obviously not fully sampled for a FoV of `600 mm`. This different trajectories (one for each frame) will be called *trajectory bins*. 

    We will often use the word *bin* to designate a pack of data that belongs to one image frame among other. The data lst `y` is also made of
    several *data bins*, one for each frame.  


2. **Resizing** `C` **, computing** `ve` **and permuting** `y`

    .. code-block:: matlab  

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


3. **Evaluating the initial image for further iterative reconstructions**

    .. code-block:: matlab  

        x0 = cell(nFr, 1); 
        for i = 1:nFr
            x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, N_u, dK_u);
        end

        bmImage(x0); 

    .. image:: ./images/data_3_mathilda.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth. 

4. **Evaluating gridding matrices**

    .. code-block:: matlab  

        [Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);


5. **Regularized least-square reconstruction with (non-motion-corrected) temporal** :math:`l_1`- **regularization for non-cartesian data (TevaMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_tevaMorph.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

6. **Estimation of deformation fields**

    .. code-block:: matlab  

        maxPixDisplacement  = 10; 
        nIter               = 500; 
        reg_file            = []; 
        reg_mask            = []; 

        h = x; 
        [DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


7. **Evaluation of deformation matrices**

    .. code-block:: matlab  

        [Tu, Tut] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);



8. **Regularized least-square reconstruction with motion-corrected temporal** :math:`l_1`- **regularization for non-cartesian data (TevaMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_tevaMorph_DF.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

9. **Regularized least-square reconstruction with backward and forward (non-motion-corrected) temporal** :math:`l_1`- **regularization for non-cartesian data (TevaDuoMorphosia_chain)**

    .. code-block:: matlab  


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

    .. image:: ./images/data_3_tevaDuoMorph.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

10. **Estimation of deformation fields**

    .. code-block:: matlab  

        maxPixDisplacement  = 10; 
        nIter               = 500; 
        reg_file            = []; 
        reg_mask            = []; 

        h = x; 
        [DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 
        [DF_to_next, imReg_to_next]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_next', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


11. **Evaluating deformation matrices**

    .. code-block:: matlab  

        [Tu1, Tu1t] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);
        [Tu2, Tu2t] = bmImDeformField2SparseMat(DF_to_next, N_u, [], true);


12. **Regularized least-square reconstruction with backward and forward motion-corrected temporal** :math:`l_1`- **regularization for non-cartesian data (TevaDuoMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_tevaDuoMorph_DF.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.


13. **Regularized least-square reconstruction with (non-motion-corrected) temporal** :math:`l_2`- **regularization for non-cartesian data (SensitivaMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_sensitivaMorph.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

14. **Estimation of deformation fields**

    .. code-block:: matlab  

        maxPixDisplacement  = 10; 
        nIter               = 500; 
        reg_file            = []; 
        reg_mask            = []; 

        h = x; 
        [DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


15. **Evalutation of deformation matrices**

    .. code-block:: matlab  

        [Tu, Tut] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);



16. **Regularized least-square reconstruction with motion-corrected temporal** :math:`l_2`- **regularization for non-cartesian data (SensitivaMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_sensitivaMorph_DF.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

17. **Regularized least-square reconstruction with backward and forward (non-motion-corrected) temporal** :math:`l_2`- **regularization for non-cartesian data (SensitivaDuoMorphosia_chain)**

    .. code-block:: matlab  

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

    .. image:: ./images/data_3_sensitivaDuoMorph.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.

18. **Estimation of deformation fields**

    .. code-block:: matlab  

        maxPixDisplacement  = 10; 
        nIter               = 500; 
        reg_file            = []; 
        reg_mask            = []; 

        h = x; 
        [DF_to_prev, imReg_to_prev]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_prev', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 
        [DF_to_next, imReg_to_next]    = bmImDeformFieldChain_imRegDemons23(  h, frSize, 'curr_to_next', nIter, 1, reg_file, reg_mask, maxPixDisplacement); 


19. **Evaluation of deformation matrices**

    .. code-block:: matlab  

        [Tu1, Tu1t] = bmImDeformField2SparseMat(DF_to_prev, N_u, [], true);
        [Tu2, Tu2t] = bmImDeformField2SparseMat(DF_to_next, N_u, [], true);


20. **Regularized least-square reconstruction with backward and ofrward motion-corrected temporal** :math:`l_2`- **regularization for non-cartesian data (SensitivaDuoMorphosia_chain)**

    .. code-block:: matlab  

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
                                            Tu1, Tu1t, Tu1, Tu1t, ...
                                            delta, regul_mode, ...
                                            nCGD, ve_max, ...
                                            nIter, witnessInfo  ); 


        bmImage(x)


    .. image:: ./images/data_3_sensitivaDuoMorph_DF.gif
        :width: 100 %
        :align: center

    Left is the reconstructed image and right the ground-truth.