This repository contains the MATLAB code and preliminary details/instructions for a project that investigates neural connectivity patterns in the PCC, mPFC, RIPC, and LIPC brain regions of healthy and autistic individuals, in the Default Mode Network. Using Dynamic Causal Modelling (DCM) in MATLAB and SPM, we explored the causal relationships between these regions, revealing distinct correlation patterns and excitatory/inhibitory connections, summarized in matrix form. The repo includes the project briefing, code for analysis, and results for running on MATLAB. 




The first step in the workflow involves preparing the raw resting-state fMRI data. Each subject’s functional scan, typically stored as a single 4D NIfTI file, is unzipped if necessary and split into individual 3D volumes, one per time point. This step ensures the data are compatible with SPM and ready for subsequent analyses, creating separate 3D scans in each subject’s functional folder. You can find this in directory_edit.m 

The second step uses these 3D volumes to build a General Linear Model (GLM) for each subject and extract time series from key brain regions of interest (ROIs), including the PCC, mPFC, LIPC, and RIPC. The GLM results are saved in `SPM.mat`, and each ROI’s activity over time is stored in `VOI_*.mat` files. This step provides the temporal data necessary to study neural interactions and prepare the network for causal modeling. You can find this in looping_dcm_viewed.m

Finally, the third step specifies and estimates a Dynamic Causal Model (DCM) using the extracted ROI time series. A fully connected network is defined, and resting-state interactions are modeled using stochastic, cross-spectral density analysis. The resulting `DCM_DMN.mat` files contain estimates of effective connectivity, including excitatory and inhibitory influences between regions, which can be further analyzed for group differences or network patterns. This three-step pipeline,from raw scans to DCM,forms the complete workflow for resting-state fMRI connectivity analysis. You can find this in spectral_dcm_batch_abide.m 



A few resources you might like to look into if you were interested in this project are as follows:
1) https://med.stanford.edu/content/dam/sm/scsnl/documents/Neuron_2023_Menon_20_years.pdf
2) https://direct.mit.edu/netn/article/8/1/178/117965


Also provided to you here is a folder containing the fMRI scans for healthy people and those with ASD: https://drive.google.com/drive/folders/1Mlruneli0G2UdmPuVk9IuAhF9j28QX5w?usp=sharing .


