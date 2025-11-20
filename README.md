## Target-Oriented MultiModal Fusion (TOMMF) -- SC-FC Covariation-Based Antidepressant Response Prediction

This repository contains the core algorithms for finding model coefficients using the TOMMF framework. 


### Framework Design
<img width="921" height="435" alt="image" src="https://github.com/user-attachments/assets/1f1f5f17-208d-4e7f-a461-3a1c2adf83d1" />



<!-- ### Advantages of TOMMF
-->

### Code Usage

For primary results, the code **initTOMMF**, **optimizeTOMMF_ContinuousSparsification** and **optimizeTOMMF_FineTuning** are used. These functions are called in order (Initialization -> Continuous Sparsification -> Fine Tuning), and the outputs of the previous step are the inputs of next. The inputs of initTOMMF are (normalized) structural connectivity, (normalized) functional connectivity, (normalized) prediction target (in our case, antidepressant response), and regularization parameters (one for latent dimensions and one for prediction weights).

optmizeTOMMF_L1_Counterpart is only used for results reported in Extended Data Fig. 1, using the same initialization function initTOMMF. The L1-counterpart is unable to dissect the psychopharmacological pattern into distinct components, therefore we did not included it in the main results. However, it may be a valid alternative for efficiency if only predictability is desired. We provide its code to both demonstrate the unique advantage endowed by L0-regularization and show its validity to achieve comparable predictability.

For analytical details, please refer to the accompanying paper's appendices. 

### Reference
Tong, Xiaoyu, et al. "Generalizable structure–function covariation predictive of antidepressant response revealed by target-oriented multimodal fusion." **_Nature Mental Health_** (2026).

preprint version: https://www.medrxiv.org/content/10.1101/2024.04.11.24305583v2)
