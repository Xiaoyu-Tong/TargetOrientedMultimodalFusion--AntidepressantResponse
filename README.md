## Target-Oriented MultiModal Fusion (TOMMF) -- SC-FC Covariation-Based Antidepressant Response Prediction

This repository contains the core algorithms for finding model coefficients using the TOMMF framework. In this study, structural and functional connectivity (SC and FC) are integrated to extract SC-FC covariation features, which in turn serve as predictors of antidepressant response.


### Framework Design
<img width="921" height="435" alt="image" src="https://github.com/user-attachments/assets/1f1f5f17-208d-4e7f-a461-3a1c2adf83d1" />

TOMMF is a end-to-end learning framework consisting of two steps: 1) multimodal fusion, and 2) target prediction. The multimodal fusion step takes SC and FC features and uses a CCA-inspired loss to identify SC-FC covariation dimensions that maximize the association between SC and FC dimensions. Essentially, SC-FC covariation features are weighted sums of SC and FC dimensions that represent their consensus information. Then these features are used to predict the clinical target of interest (e.g., treatment response). Notably, in addition to this forward propogation, there is also a reciprocal feedback from prediction target to multimodal fusion ensuring the clinical relevance of identified SC-FC covariation dimensions. The bidirectional propagation yields covariation features representing both SC-FC interplay and target-relevant information.

### Advantages of TOMMF
Three key advantages of TOMMF are endowed by its design: 1) Enhanced target relevance for multimodal features, 2) Ability to examine the interplay between modalities, and 3) Enhanced robustness of biomarkers. First, the end-to-end design allows backpropagation of target information to multimodal fusion that tailors the SC-FC covariation features to the target of interest. Second, the design to maximize the association between SC and FC dimension scores enables us to identify SC and FC abnormalities linking to each other, revealing interaction/relationship between SC and FC variations, especially those involving different region pairs. Furthermore, the emphasis on the information in agreement helps the model to extract more reliable features, conceptually SC-approved FC and FC-approved SC. 


<!-- ### Key Results
**0. Generalizable treatment predicion for sertraline and placebo**
xxx
-->

### Code Usage

For primary results, the code **initTOMMF**, **optimizeTOMMF_ContinuousSparsification** and **optimizeTOMMF_FineTuning** are used. These functions are called in order (Initialization -> Continuous Sparsification -> Fine Tuning), and the outputs of the previous step are the inputs of next. The inputs of initTOMMF are (normalized) structural connectivity, (normalized) functional connectivity, (normalized) prediction target (in our case, antidepressant response), and regularization parameters (one for latent dimensions and one for prediction weights).

optmizeTOMMF_L1_Counterpart is only used for results reported in Extended Data Fig. 1, using the same initialization function initTOMMF. The L1-counterpart is unable to dissect the psychopharmacological pattern into distinct components, therefore we did not included it in the main results. However, it may be a valid alternative for efficiency if only predictability is desired. We provide its code to both demonstrate the unique advantage endowed by L0-regularization and show its validity to achieve comparable predictability.

For analytical details, please refer to the accompanying paper's appendices. 

### Reference
Tong, Xiaoyu, et al. "Generalizable structure–function covariation predictive of antidepressant response revealed by target-oriented multimodal fusion." **_Nature Mental Health_** (2026).

