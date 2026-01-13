======================
Theoretical Background
======================


As emphasized by the title of this blog, we will start from theoretical concepts in order to introduce the reconstructions functions in a second step.
If you will however, you can directly go to the tutorials and run the demonstration scripts. You can learn that way very fast how to use the reconstructions
of *Monalisa*. But if you want to grasp some understanding of what those functions are doing in order to modify them or write new reconstructions, you may need the
theoretical minimum that we present in this section. 

We will not go into the beautiful theory of convex optimization. We will therefore not (or very partially) describe the reconstruction algorithms, 
and we will not at all provide any proof of convergence of those algorithm to a minimizer of the reconstruction related *argmin*-problems. Rather, we
will write the reconstruction problem as formally as possible and define with a certain level of details each involved ingredient (linear maps, norms, 
inner-product, regularizers, ...) of the reconstruction problem. We follow the believe that a good description of the problem is already half of 
its solution. 

But we need a start point to begin that theory and since the goal of MRI reconstruction is to build an *image* from MRI data 
in order to present it to a radiologist or cardiologist for interpretation, a good start to enter in the theory may therefore be to explain 
what we mean by an *image*. Of course, an image can simply be defined by a list of numerical values, one for each voxel center. 
Or more generaly, one can also assigne :math:`m` (possibly all diferent) values to each voxel and call it an *image* (or *multimodal image*). 
But this does not say what is the meaning of the involved numerical values. This meaning can be very different depending on what we choose 
to reconstruct. We begin therefore this theory section by giving  information about the physical meaning of the numerical 
values that consitute the images reconstructed by our toolbox. This brings us necessarily by documenting the notion of *magntization*, 
which is the central physical notion in MRI. 


.. toctree::
   :maxdepth: 1

   2_1_MRI_physics
   2_2_recon_problem
   2_3_algorithms



