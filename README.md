# Target-Oriented MultiModal Fusion (TOMMF) -- Antidepressant Response

**10.6.2024 Update (Exclusive use for editors and reviewers)**

This repository contains the core algorithms for finding model coefficients using the TOMMF framework. For primary results, the code **initTOMMF**, **optimizeTOMMF_ContinuousSparsification** and **optimizeTOMMF_FineTuning** are used. These functions are called in order (Initialization -> Continuous Sparsification -> Fine Tuning), and the outputs of the previous step are the inputs of next. 

optmizeTOMMF_L1_Counterpart is only used for results reported in Extended Data Fig. 1, using the same initialization function initTOMMF. The L1-counterpart is unable to dissect the psychopharmacological pattern into distinct components, therefore we did not included it in the main results. However, it may be a valid alternative for efficiency if only predictability is desired. We provide its code to both demonstrate the unique advantage endowed by L0-regularization and show its validity to achieve comparable predictability.
