function [WS,WF,G,beta,L_history,iIter] = optimizeTOMMF_L1_Counterpart(rhoS,rhoF,XS,XF,G,Y,WS,WF,beta,Kpred,lambda_fusion,lambda_pred,max_iter)
% % Optimize W_S, W_F, and beta for the target-oriented multimodal fusion (TOMMF) procedure % %
% 
% Input:
% XS -- data matrix for structural connectivity (N x D_S)
% XF -- data matrix for structural connectivity (N x D_F)
% Y --  data vector for prediction target (N x 1)
% WS -- initial guess of WS (D_S x P)
% WF -- initial guess of WF (D_S x P)
% beta --  initial guess of beta (P x 1)
% Kpred -- the relative weight of prediction terms with respect to multimodal fusion
% lambda_fusion -- hyperparameter for multimodal fusion (L1-regularization)
% lambda_pred -- hyperparameter for prediction task (L1-regularization)
% max_iter -- max number of iterations
% *** Note, XS and XF should be normalized so that regularization is fair
% for each features **
%
% Output:
% WS -- initial guess for WS
% WF -- initial guess for WF
% beta -- initial guess for beta
% L_history -- record of loss function values
%
% by Xiauyu Tong, Lehigh, 2023-6
% xit321@lehigh.edu

tol = 10^-6;
iIter = 0;

% L: the loss function
L_new = norm(G - XS*WS,"fro")^2 + norm(G - XF*WF,"fro")^2 + ...
    lambda_fusion * (sum(sum(abs(WS)))+sum(sum(abs(WF)))) + ...
    Kpred * (norm(Y-G*beta,"fro")^2 + lambda_pred*sum(abs(beta)));
L = 2*L_new;

% Define learning rate
mu = 0.1;

L_history = [];
while (abs(L-L_new) > tol) && iIter < max_iter
    if L_new < tol
        break
    end
    L = L_new;
    iIter = iIter + 1;
    
    beta_new = softthresh(beta-2*mu*Kpred*G'*(G*beta-Y),mu*lambda_pred);
    for j = 1:length(beta)
        WS_new(:,j) = softthresh(WS(:,j)-2*mu*XS'*(XS*WS(:,j)-G(:,j)),mu*lambda_fusion);
        WF_new(:,j) = softthresh(WF(:,j)-2*mu*XF'*(XF*WF(:,j)-G(:,j)),mu*lambda_fusion);
    end
    rho = rhoF + rhoS;
    G_new = (rhoS*XS*WS_new + rhoF*XF*WF_new + Kpred*Y*beta_new')*(eye(length(beta_new))/rho-(Kpred*(beta_new*beta_new'))/(rho^2+rho*Kpred*(beta_new'*beta_new)));
   

    L_new = norm(G_new - XS*WS_new,"fro")^2 + norm(G_new - XF*WF_new,"fro")^2 + ...
        lambda_fusion * (sum(sum(abs(WS_new)))+sum(sum(abs(WF_new)))) + ...
        Kpred * (norm(Y-G_new*beta_new,"fro")^2 + lambda_pred*sum(abs(beta_new)));
    L_history = [L_history,L_new];
    if L_new > L && iIter > 1
%             mu_fusion = mu_fusion/2;
%             mu_pred = mu_pred/2;
        mu = mu/2;
        disp('learning rate reduced')
    else    
        WS = WS_new;
        WF = WF_new;
        beta = beta_new;
        G = G_new;
    end
end

end