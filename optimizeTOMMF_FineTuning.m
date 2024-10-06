function [WS,WF,beta,nNet,L_history] = optimizeTOMMF_FineTuning(rhoS,rhoF,XS,XF,Y,WS,WF,beta,MS,MF,Mbeta,Kpred,L_history,max_iter)
% % Optimize W_S, W_F, and beta for the target-oriented multimodal fusion (TOMMF) procedure % %
% 
% Input:
% XS -- data matrix for structural connectivity (N x D_S)
% XF -- data matrix for structural connectivity (N x D_F)
% Y --  data vector for prediction target (N x 1)
% WS -- initial guess of WS (D_S x P)
% WF -- initial guess of WF (D_S x P)
% beta --  initial guess of beta (P x 1)
% MS -- mask for structural connectivity (D_S x P)
% MF -- mask for functional connectivity (D_F x P)
% Mbeta -- mask for prediction weight (P x 1)
% Kpred -- the relative weight of prediction terms with respect to multimodal fusion
% lambda_fusion -- hyperparameter for multimodal fusion (L0-regularization)
% lambda_pred -- hyperparameter for prediction task (L0-regularization)
% max_iter -- max number of iterations
% *** Note, XS and XF should be normalized so that regularization is fair
% for each features ***
%
% Output:

%
% by Xiauyu Tong, Lehigh, 2023-7
% xit321@lehigh.edu

tol = 10^-6;
iIter = 0;
L_new = L_history(end);
L = 2*L_new;
% L_history = L_new;

H = @(x) x>=0; % Heaviside step function
hidxs = H(MS);
hidxf = H(MF);
hidxb = H(Mbeta);
if rhoS==0
    hidxb = (sum(hidxf)~=0)' & hidxb;
elseif rhoF==0
    hidxb = (sum(hidxs)~=0)' & hidxb;
else
    hidxb = (sum(hidxf)~=0)' & (sum(hidxs)~=0)' & hidxb;
end
% WS_full = WS;
% WF_full = WF;
% beta_full = beta;
WS = WS(:,hidxb);
WF = WF(:,hidxb);
beta = beta(hidxb);
hidxs = hidxs(:,hidxb);
hidxf = hidxf(:,hidxb);

mu_fusion = 0.01;
mu_pred = 0.01;
rho = rhoS + rhoF;
G = (rhoS*XS*(hidxs.*WS) + rhoF*XF*(hidxf.*WF) + Kpred*Y*beta')*(eye(length(beta))/rho-(Kpred*(beta*beta'))/(rho^2+rho*Kpred*(beta'*beta)));
while (abs(L-L_new) > tol) && iIter < max_iter
    if L_new < tol
        break
    end
    L = L_new;
    iIter = iIter + 1;
    beta_new = beta - mu_pred*2*Kpred*G'*(G*beta-Y);
    
    WS_new = WS - (2*mu_fusion*rhoS*XS'*(XS*(hidxs.*WS)-G)).*hidxs;
    WF_new = WF - (2*mu_fusion*rhoF*XF'*(XF*(hidxf.*WF)-G)).*hidxf;

    G_new = (rhoS*XS*(hidxs.*WS_new) + rhoF*XF*(hidxf.*WF_new) + Kpred*Y*beta_new')*(eye(length(beta_new))/rho-(Kpred*(beta_new*beta_new'))/(rho^2+rho*Kpred*(beta_new'*beta_new)));

    L_new = rhoS*norm(G_new - XS*(hidxs.*WS_new),"fro")^2 + rhoF*norm(G_new - XF*(hidxf.*WF_new),"fro")^2 + ...
        Kpred * norm(Y-G_new*beta_new,"fro")^2;

    L_history = [L_history,L_new];
    if L_new > L && iIter > 1
        mu_fusion = mu_fusion/2;
        mu_pred = mu_pred/2;
        disp('learning rate reduced')
    else    
        WS = WS_new;
        WF = WF_new;
        beta = beta_new;
        G = G_new;
    end 
end

WS = WS.*hidxs;
WF = WF.*hidxf;
nNet = sum(hidxb);


end









