function [WS,WF,G,beta] = initTOMMF(XS,XF,Y,lambda_fusion,lambda_pred)
% % Initialize W_S, W_F, and beta for the target-oriented multimodal fusion (TOMMF) procedure % %
% 
% Input:
% XS -- data matrix for structural connectivity (N x D_S)
% XF -- data matrix for functional connectivity (N x D_F)
% Y --  data vector for prediction target (N x 1)
% lambda_pred -- hyperparameter for prediction task (L2-regularization)
% *** Note, XS and XF should be normalized so that regularization is fair
% for each features **
%
% Output:
% WS -- initial guess for WS
% WF -- initial guess for WF
% beta -- initial guess for beta
%
% by Xiauyu Tong, Lehigh, 2023-6
% xit321@lehigh.edu

% [WS,WF] = canoncorr(XS,XF);
ncomp = min([length(Y),size(XS,2),size(XF,2)])-1;
max_iter = 5; % just for initialization, the rough result is OK
rng(0);% for replicability
[WS,WF,~,~] = SparseCCA_PMD_noNorm(XS,XF,ncomp,max_iter,lambda_fusion,lambda_fusion);
G = (XS*WS+XF*WF)/2;
try
    [beta,~]=lasso(G,Y,'lambda',lambda_pred,'standardize',0);
catch
    beta = zeros(size(G,2),1);
end

end
