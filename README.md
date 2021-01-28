# A Modified Fourier-Mellin Approach for Source Device Identification on Stabilized Videos
<img src="assets/scheme.png" width="500">

This is the official repository of **A Modified Fourier-Mellin Approach for Source Device Identification on Stabilized Videos**,
accepted to [IEEE Internation Conference on Image Processing (ICIP) 2020](https://2020.ieeeicip.org/), pp. 1266-1270, 2020 and available on [arXiv](https://arxiv.org/pdf/2005.09984.pdf).

## Code

### Requirements

- Install MATLAB Image Processing Toolbox, Global Optimization Toolbox
- Download the [Camera-fingerprint](http://dde.binghamton.edu/download/camera_fingerprint) package and run the function `compile.m` in folder `CameraFingerprint/Filter`.

### Pipeline

Given a generic video query and a reference device:
- extract the I-frames
- extract the noise residuals of the I-frames using the Camera-fingerprint package and save them as `test_noises.mat`
- extract the device PRNU using the Camera-fingerprint package, then scale and crop it as suggested in [1]. Save it as `K.mat`
- run the function `MFM_deltarho_main.m` to test the video

### Example result
Run the function `MFM_deltarho_main.m` to evaluate the proposed Modified Fourier Mellin method over a sample query video.  how the 
- [notebook showing the results](show_results.ipynb)
- You can find the complete list of results for every model [here](outputs/)

### Extract image noise residuals and device PRNU and save them 
You can extract them using the Python implementation available [here](https://github.com/polimi-ispl/prnu-python).  
For each device, create a train-validation-test split, dividing the image noise residuals in 50% training, 25% validation, 25% evaluation.  
Create 3 lists for each device reporting the paths to the noise residuals: 
- `/Noises_lists/train/list_%device_name.npy`
- `/Noises_lists/valid/list_%device_name.npy`
- `/Noises_lists/test/list_%device_name.npy`


## References
[1] S. Mandelli, P. Bestagini, L. Verdoliva, S. Tubaro, *Facing device attribution problem for stabilized video sequences*
IEEE Transactions on Information Forensics and Security 15, 14-27, 2020.
