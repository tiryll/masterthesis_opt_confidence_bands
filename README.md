
----------------------
i) General Information
----------------------

This repository contains the files to replicate the results presented in my masterthesis: ***Optimal Simultaneous Confidence Bands for Impulse Response Functions of SVAR-Models***.

The folder *src* contains the source code, while the folder *bld* contains the raw data for the application considered in the thesis. By running the code from *src*, the subfolders of *bld*, *data* and *figures*, will be created and/or expanded to include intermediate and final results.

The code has been entirely written in MATLAB and requires at least MATLAB 2023a and the parallel computing toolbox.

---------------------
ii) Source Code & Data Overview 
---------------------

The structure of the *src* folder is as follows:

- *application*
- *MATLAB*
- *simulation*

The folder *application* contains a MATLAB-live script that creates Figure 3, discussed in section 5. The folder *MATLAB* contains all functions used in the scripts. The functions with the ending "vec" are vectorized versions of the original functions and perform the calculations for all Monte Carlo Simulations at once.
The folder *simulation* contains the source code for Figure 1, Figure 2, Table 1 and Table 2 in the respective subfolders *figure_1*, *table_1* (also includes code for Figure 1) and *table_2*.

The computations for Figure 2, Table 1 and Table 2 in *table_1* and *table_2 *are split into several numbered scripts and have to be executed in ascending order. 
