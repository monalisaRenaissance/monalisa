
# Welcome to Monalisa

The Monalisa toolbox is devoted to the reconstruction of both non-Cartesian and Cartesian MRI data. The present repository is  maintained by the original developpers of Monalisa. If you need some guidance, you may want to visit the [documentation](https://monalisarenaissance.github.io/monalisa/) webpage. 

## Installation

The intallation is detailed [here](https://monalisarenaissance.github.io/monalisa/1_installation.html) on the documentation webpage. For short, here are the main steps: 

1. Create a folder with a name such as "monalisaRenaissance", open a terminal and go into that folder. Then copy the URL from the green button `<> Code` here above and type in your terminal: 

    ```cmd
    git clone <the URL you just copied>
    ```

2. Open Matlab and go into the *monalisa* directory you just cloned using the Matlab console. For me, that was done by typing in Matlab the command

    ```cmd
    cd c:\main\repos\monalisaRenaissance\monalisa
    ```

3. Now check from Matlab that you have a c++ compiler installed (e.g. *gcc* on Linux and *Visual Studio C++ compiler* on Windos) by typing in the Matlab console

    ```cmd
    mex -setup C++
    ```

    You have to do the necessary until you see an answer such as 

    ```cmd
    MEX configured to use 'Microsoft Visual C++ 2022' for C++ language compilation.
    ```

    Note that you may have a c++ compiler installed but Matlab may not point to it. For that, use the command

    ```cmd
    mex -setup
    ```

    and choose the compiler in the list. 

4. In the *Environment* tab of Matlab, go to *Set Path* and *Add folder with subfolders*. Then select the monalisa folder and validate, so that the paths to all files of the toolbox are added to the Matlab list of paths. Finally, use the `bmMex` function to compile the c++ code of Monalisa. Simply type in the consol

    ```sh
    bmMex; 
    ```

    It should last a few seconds or a few minutes. If no error occured, you are done ! 
    Just be sure that `libomp` is installed on your computer to that multithreading is enabled for better performances. 




## Overview of the **Monalisa directory**

All the code, documentation, demos and tests are contained in the following sub folders :  

- `src` : this is the toolbox folder that contain the code of Monalisa. That is the only folder subject to licensing as detailed hereafter. Out policy is to include in it only: other subfolders, functions (matlab, c++ or other), classes and text files containing commands for compilation. But it should not contain data or scripts or other things. Those should be project specific and should not be shared in the toolbox. We also don't include unit tests in the `src` folder. Don't modify that folder unless you are advanced user. 

- `third_part` : The content of this subfolder is not part of Monalisa. It was written by other authors and is freely available elsewhere under some other licensing. It contains some dependencies of Monalisa that are not mandatory to use for reconstrucitons but they may be very valuable for certain tasks such as extracting data from rawdata files.   

- `demo` : Since the `src` directory does not contain data nor scripts, we built a `demo` directory to include some demonstration scripts as well as some demonstration data to run the scripts. These scripts are mainly dedicated to show how to call the reconstruction functions, which is the last step in a reconstruction. These demonstrations allows therefore see some results wihtout having to care about all preparation steps. 

- `unit_test` : Although still in construction, this subfolder is aimed to contain unit tests for some important functions of the toolbox. It will also contains some scripts that call those unit tests.

- `docs` : This is the subfolder that contains all the content of the documentation webpage and the information to build it. It has no function role in the toolbox. 

## Licensing and Contributions

This code is protected by a license written in the LICENSE file. This license was originally written by the Lausanne University Hospital (CHUV)
because it is where the toolbox was initiated. The politic for this licensing is to add as licnesor the institution where each new contributor is working while contributing to the toolbox. If the nummber of licensors become too large,  this solution may become unpractical and a more open solution will be found. 

For the moment, the toolbox can be used freely for any resonable personal or academic applications, including publications. But any redistribution or commercial usage must be done in collaboration with lincensors. 
