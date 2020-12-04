# AQuaSI

This repository contains the implementation of the paper

F. Schirrmacher, C. Riess and T. Köhler, "Adaptive Quantile Sparse Image (AQuaSI) Prior for Inverse Imaging Problems," in IEEE Transactions on Computational Imaging, vol. 6, pp. 503-517, 2020, doi: 10.1109/TCI.2019.2956888. [IEEE](https://ieeexplore.ieee.org/abstract/document/8931625) 

If you use this code in your work, please cite:

        @ARTICLE{8931625,
                author={F. {Schirrmacher} and C. {Riess} and T. {K\"ohler}},
                journal={IEEE Transactions on Computational Imaging}, 
                title={Adaptive Quantile Sparse Image (AQuaSI) Prior for Inverse Imaging Problems}, 
                year={2020},
                volume={6},
                pages={503--517},
                doi={10.1109/TCI.2019.2956888}}


## Getting started

To download the code, fork the repository or clone it using the following command:

```
  git clone https://github.com/franziska-schirrmacher/AQuaSI.git
```

### Code structure

- **cpp**: This folder contains the .cpp files of the AQuaSI and QuaSI implementation

- **data**: This folder contains all datasets for each of the experiments (needs to be created)

- **matlab**: The matlab folder contains the matlab source code of the experiments and the proposed algorithm
    - **algorithms**: 
    
        - **MCAQuaSI**: Contains the implementation for multi-channel (color images) inverse problems 

        - **SCAQuaSI**: Contains the implementation for single-channel (grayscale images) inverse problems 

    - **evaluation**: Contains the scripts for the experiments

    - **utility**: Contains helper functions

- **results**: Folder to store the results (needs to be created)


### Mex function

There is a pre-build mex function for Windows in the folder **matlab/utility**. If this does not work for you, build the mex function on your system:
1. Build the mex function for main_aquasi_mex.cpp and main_quasi_mex.cpp in the folder **cpp**:

        mex -R2018a main_aquasi_mex.cpp
        mex -R2018a main_quasi_mex.cpp

2. Copy the mex functions to **matlab/utility**

In case this is not working for you:
Comment the lines 21, 32, and 46 in **matlab/algorithms/SCAQuaSI/admmSC.m** and uncomment the lines 25, 36, and 49.

Comment the lines 37, 45, and 55 in **matlab/algorithms/SCAQuaSI/admmMC.m** and uncomment the lines 40, 48, and 58.

The computation with this configuration takes noticably longer.

### Datasets

- **ARRI**: Please download the ARRI dataset and store the images in the **data/ARRI** folder
- **bowls**: Please download the results from http://www.cse.cuhk.edu.hk/~leojia/projects/crossfield/ and store the bowls folder in the **data/bowls** folder


## Experiments


### IV.C. Joint filtering

In Matlab:

1. Change your current folder to **matlab/evaluation**
2. Run the script: 

        evaluateBowls.m

### IV.D. Implementation Remarks

In Matlab:

1. Change your current folder to **matlab/evaluation**
2. Run the script: 

        evaluateConvergence.m


### V.A. Comparison to Total Variation Regularization

In Matlab:

1. Change your current folder to **matlab/evaluation**
2. Run the script: 

        evaluateAQuASIvsTV.m


### V.D. RGB/NIR Cross Field Image Restoration

In Matlab:

1. Change your current folder to **matlab/evaluation**
2. Run the script: 

        evaluateMCAQuaSI.m
