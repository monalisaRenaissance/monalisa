Tutorial 1 : A first gridded-zero-padded reconstrcution for non-cartesian data
==============================================================================

The present script is written in the file `demo_script_1.m` in the demonstration folder. You can open it and run it right away. 

We present here **Mathilda**: a non-iterative reconstruction for non-cartesian data. When using Monalisa, this is usually the reconstruction we start with to make a first image of the 
data we want to work on. 

The script contains the following section. 

1. **Initialization**

    The first section identify where the running script is located and deduces to location of 
    the demonstration data in order to load it. 

    .. code-block:: matlab  

        myCurrent_script_file = matlab.desktop.editor.getActiveFilename; 
        d = fileparts(myCurrent_script_file); 
        load([d, filesep, '..', filesep, 'demo_data', filesep, 'demo_data_1']); 


    You can plot the sampling trajectory `t` to have an idea of how good the data are sampled. 
    You may get someting as follows:      

    .. image:: ./images/t_radial_fullySamp.png
        :width: 100 %
        :align: center

    .. raw:: html

        <div style="margin-bottom: 30px;"></div>

    This trajectory contains 256 lines of 512 points each covering a disk homogeneously. There is no well known notion
    of "full sampling" for non-cartesian (non-uniform) trejectories for the simple reason that there exist no well known 
    sampling theorem for non-cartesian trajectory. 
    
    We will however qualify informally a non-cartesian trajectory as **fully sampled** with respect to a k-space 
    step-size :math:`\Delta k` if the trajectory is sufficiently well sampled so that the 
    interpolation of the non-cartesian data on a cartesian gridd of step-size
    :math:`\Delta k` is very close to fully sampled data. Just keep in mind that it is a somewhat abusive language. 

    This being said, we will considere that the above trajectory is fully sampled for a field of view (FoV) of
    *600 mm* (i.e. :math:`\Delta k = \frac{1}{600} \frac{1}{mm}`)  


2. **Gridded-zeropadded reconstruction without coil combination**

    The reconstruction for gridded-zero-padded reconstruction in *Monalisa* is called *Mathilda*
    (function `bmMathilda`). In abscence of coil-sensitivity estimation *C*, the function can be called as follows:  

    .. code-block:: matlab

        x = bmMathilda(y, t, ve, [], N_u, N_u, dK_u);
        bmImage(x); 

    The viewer (`bmImage`) displays then the the resulting list of images. You can brows through the list with the up- and down-arrows. The following figure shows
    a fiew of them: 

    .. image:: ./images/coil_images_demo_1.png
        :width: 100 %
        :align: center

    .. raw:: html

        <div style="margin-bottom: 30px;"></div>


    There is one image for each coil (i.e. each channel). This is why we call these images the "coil-images". In some sense, if a coil was an eye, 
    its coil-image would be the image that this coil "sees". 
    
    The function call of *Mathilda* as written in the section abbove can be described mathematically as a discrete (approximative inverse) non-uniform 
    Fourier transform for each coil. To describe it with mathematical symbols, we will call :math:`\acute{y}_c` the data vector
    of coil number :math:`c` (or channel number:math:`c`) and we will write :math:`\acute{x}_c^{\#}` the (unkown) ground-truth coil-image of 
    coil number :math:`c` sampled on the Cartisian image gridd. Image :math:`\acute{x}_c^{\#}` is the coil-image number :math:`c` that we would like to 
    ideally obtain. In the following, the word *channel* will be used interchangebly with the world *coil* and :math:`nCh` will mean *number of channels*.  

    As justified by the signal equation, the approximative link between :math:`\acute{x}_c^{\#}` and :math:`\acute{y}_c` is a discrete non-uniform Fourier transform
    that we will write :math:`\acute{F}_c`: 
    
    .. math ::

        \acute{y}_c =  \acute{F}_c \cdot \acute{x}_c^{\#}
    
    There is one linear map :math:`\acute{F}_c` for each coil and they are all identical. 
    We have however labeled them with index :math:`c` to consruct one identical copy of that Fourier transform for each coil. 
    The index :math:`c` goes from :math:`1` to :math:`nCh` (which stands for *number of channels*). 
    
    The most important thing to know about :math:`\acute{F}_c` is that it has no inverse in general. So far the author knows, the only special case
    where :math:`\acute{F}_c` has an inverse is when the sampling trajectory is uniform (Cartesian) fully sampled. Then becomes :math:`\acute{F}_c` the 
    usual uniform discrete Fourier transform with a well defined inverse. But for non-cartesian data, :math:`\acute{F}_c` has usually no inverse, 
    and even if it could have some in very special cases, we will not take time to construct one because it would be of no help for practical cases. 
    We considere that :math:`\acute{F}_c` has therefore not inverse.  
    
    Instead of building an inverse (which is anywhay not possible most of the time), we build a linear map, that is easy to implement, and that behave 
    approximately like an inverse when the data are well sampled. We will write it :math:`\acute{F}_c^{\sim 1}` and by our description it holds 
    
    .. math ::

        \acute{F}_c^{\sim 1} \cdot \acute{F}_c \approx  id

    The function *Mathilda* realizes such a map :math:`\acute{F}_c^{\sim 1}` for non-cartesian data. We give as argument a list of data vectors
    :math:`\acute{y}_1, ..., \acute{y}_{nCh}` and it returns a list of coil-images :math:`\acute{x}_1, ..., \acute{x}_{nCh}` which verify

    .. math ::

        \acute{x}_c =  \acute{F}_c^{\sim 1} \cdot \acute{y}_c \approx \acute{F}_c^{\sim 1} \cdot \acute{F}_c \cdot \acute{x}_c^{\#} \approx \acute{x}_c^{\#}
   
    
    You can combine the obtained coil-images :math:`\acute{x}_1, ..., \acute{x}_{nCh}` by a sum-of-squares operation but
    then the phase would be lost and the result would not spacially homogeneous. You need a coil-sensitivity estimation in order to combine them properly as 
    in the next section. 


3. **A look at coil-sensitivity maps**

    If a coil-sensitivity estimation *C* is present, it can be passed to *Mathilda* to perform a coil-combination. First take a look at the list of coil-sensitivity
    maps that are provided on the demonstration data. You can view it by typing

    .. code-block:: matlab

        bmImage(cat(2, real(C), imag(C))); 

    The real part of the coil-sensitivities will appear on the left, and the imaginary part on the right. Us the up- and down-arrows to brows through the different
    coil-sensitivity maps. Here are a few examples of coil-sensitivities you will see : 

    .. image:: ./images/data_2_coil_sense.png
        :width: 100 %
        :align: center

    On that figure, we displayed the real-part on the top and the imaginary part on the bottom line. 


4. **Gridded-zeropadded reconstruction with coil combination**

    You can then call *Mathilda* by passing the coil-sensitivity list `C` as argument as follows: 

    .. code-block:: matlab

        x = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);
        
        bmImage(angle(x));
        bmImage(x); 

    The obtained image is the a coil-combined image of all coil-images obtained at step 2 and it magnitude and phase looks as follows: 

    .. image:: ./images/magn_phase_demo_1.png
        :width: 100 %
        :align: center

    .. raw:: html

        <div style="margin-bottom: 30px;"></div>

    The coil-combination is performed by applying the coil-sensitivity pseudo-inverse to the list of coil-images. 

Further Explanation
-------------------

You did your first image reconstruction with Monalisa ! That is an event. 

We want to take that opportunity to introduce a few basic important definitions and some vocabulary of Monalisa. We will do that by
giving first a brief description of each argument passed to Mathilda:   

    - **y**       : List of data vectors. Each entry of that table is complex-valued. We shape it in the size `nPt x nCh` when passed to a reconstruction function, where `nPt` stands for 'number of points' and `nCh` for 'number of channels'. 
    - **t**       : List of k-space trajectory points. Each entry is real-valued. We shape it usually in the size `frDim x nPt` or a reshapable size. Here stands `frDim` for 'frame dimension' and can be 1, 2 or 3. It is the number of spatial dimensions in the image. 
    - **ve**      : List of volume element. Each entry is real-valued. There is one volume element for each point of the k-space trajectory and can be considered to be 1 divided by the k-space density compensation.  
    - **C**       : List of estimated coil-sensitivity maps. Each entry is complex-valued. There is one map for each channel. Each map has `frDim` dimension.  
    - **N_u**     : This is the grid-size of the k-space Cartesian grid used for regridding. It is a vector with `frDim` components. In our example it was `[512, 512]`. 
    - **dK_u**    : This is the grid-spacing of the k-space Cartesian grid used for regridding. It is a vector with `frDim` components. In our case, `frDim` equals 2 and we will write `dK_u` as `[dKx, dKy]`. The field of view (FoV) resulting of that grid is `[1/dKx, 1/dKy]`.  


Another quantity that is not explicitely present among that list is the *frame-size* that we write `frSize`. It is a vector with `frDim` entries like `N_u` and `dK_u`. The frame-size
epresses the spatial size of the finale reconstructed image(s). In the present case (and actually in almost all our example), we have set the frame size equal to `N_u`. 

