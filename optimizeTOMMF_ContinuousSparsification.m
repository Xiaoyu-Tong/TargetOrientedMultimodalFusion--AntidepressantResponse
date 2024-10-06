function [WS,WF,G,beta,MS,MF,Mbeta,alphaW,alphaBeta,L_history,L_sub] = optimizeTOMMF_ContinuousSparsification(rhoS,rhoF,XS,XF,G,Y,WS,WF,beta,Kpred,lambda_fusion,lambda_pred,max_iter)
% % Optimize W_S, W_F, and beta for the target-oriented multimodal fusion (TOMMF) procedure % %
% 
% Input:
% XS -- data matrix for structural connectivity (N x D_S)
% XF -- data matrix for structural connectivity (N x D_F)
% Y --  data vector for prediction target (N x 1)
% WS -- initial guess of WS (D_S x P)
% WF -- initial guess of WF (D_S x P)
% G -- initial guess of G (N x P)
% beta --  initial guess of beta (P x 1)
% Kpred -- the relative weight of prediction terms with respect to multimodal fusion
% lambda_fusion -- hyperparameter for multimodal fusion (L0-regularization)
% lambda_pred -- hyperparameter for prediction task (L0-regularization)
% max_iter -- max number of iterations
% *** Note, XS and XF should be normalized so that regularization is fair
% for each features ***
%
% Output:

%
% by Xiauyu Tong, Lehigh, 2023-6
% xit321@lehigh.edu

tol = 10^-6;
KW = 2^(1/5);%1.5^(1/5)
alphaW = 1/KW; % the first alpha is 1
Kbeta = 1.5^(1/5);%2^(1/5)
alphaBeta = 1/Kbeta; % the first alpha is 1
iIter = 0;

sigm = @(x) 1./(1+exp(-x));
dsigm = @(a,m) a./(1+exp(-a*m))-a./((1+exp(-a*m)).^2);
MS = zeros(size(WS)); % Initialization of mask matrices
MF = zeros(size(WF));
Mbeta = zeros(size(beta));

% L: the loss function
L_new = rhoS*norm(G - XS*(sigm(alphaW*MS).*WS),"fro")^2 + rhoF*norm(G - XF*(sigm(alphaW*MF).*WF),"fro")^2 + ...
    lambda_fusion * (sum(sum(abs(sigm(alphaW*MS))))+sum(sum(abs(sigm(alphaW*MF))))) + ...
    Kpred * norm(Y-G*(sigm(alphaBeta*Mbeta).*beta),"fro")^2 + lambda_pred*sum(abs(sigm(alphaBeta*Mbeta)));
L = 2*L_new;

% C = eye(length(beta)); % Initialization of inverse Hessian for BFGS

% Define learning rate
% mu = 10^-7;
% mu_G = mu;
% mu_fusion = mu;
% mu_pred = mu;

% mu_G = 0.1;
mu_fusion = 0.01;%1
mu_maskB = 1;%1%0.1%0.01
mu_mask = 1;%1%0.1%0.01
mu_pred = 0.01;%0.005%0.01
% mu_ortho = 1;

% rsq=[];
L_history = L_new;
L_sub = [];
alpha_thresh = 300;%100
while alphaBeta < alpha_thresh
    alphaW = alphaW * KW;
    alphaBeta = alphaBeta * Kbeta;
%     if alpha > 1
%         Mbeta = Mbeta/K;
%         MS = MS/K;
%         MF = MF/K;
%     end
    if alphaBeta >= alpha_thresh
        max_iter = max_iter*5;%1
    end
%     alpha = 10;
    while (abs(L-L_new) > tol) && iIter < max_iter
        if L_new < tol
            break
        end
        L = L_new;
        iIter = iIter + 1;
        
%         G_new = (XS*(sigm(alpha*MS).*WS) + XF*(sigm(alpha*MF).*WF) + Kpred*Y*beta')*(eye(length(beta))/2-(Kpred*(beta*beta'))/(4+2*Kpred*(beta'*beta)));
%         G = G_new;
    
%         % Quasi-Newton BFGS algorithm for beta optimization
%         dBeta = 2*Kpred*(G'*G*beta + lambda_pred*beta - G'*Y); % Gradient
%         s = -mu_pred*C*dBeta;
%         beta_new = beta + s;
%         y = Kpred*(G'*G*beta_new + lambda_pred*beta_new) - Kpred*(G'*G*beta + lambda_pred*beta);
%         C_new = C + (s'*y+y'*C*y)*(s*s')/((s'*y)^2) - (C*y*s'+s*y'*C)/(s'*y);
%         C = C_new;
        
        beta_new = beta - mu_pred*2*Kpred*G'*(G*(sigm(alphaBeta*Mbeta).*beta)-Y).*sigm(alphaBeta*Mbeta);
        Mbeta_new = Mbeta - mu_maskB*(2*Kpred*G'*(G*(sigm(alphaBeta*Mbeta).*beta)-Y).*beta+lambda_pred).*dsigm(alphaBeta,Mbeta);
        
        WS_new = WS - (2*mu_fusion*rhoS*XS'*(XS*(sigm(alphaW*MS).*WS)-G)).*sigm(alphaW*MS);
        WF_new = WF - (2*mu_fusion*rhoF*XF'*(XF*(sigm(alphaW*MF).*WF)-G)).*sigm(alphaW*MF);
        MS_new = MS - mu_mask*((2*rhoS*XS'*(XS*(sigm(alphaW*MS).*WS)-G)).*WS+lambda_fusion).*dsigm(alphaW,MS);
        MF_new = MF - mu_mask*((2*rhoF*XF'*(XF*(sigm(alphaW*MF).*WF)-G)).*WF+lambda_fusion).*dsigm(alphaW,MF);

        gamma = sigm(alphaBeta*Mbeta_new).*beta_new;
        rho = rhoF + rhoS;
        G_new = (rhoS*XS*(sigm(alphaW*MS_new).*WS_new) + rhoF*XF*(sigm(alphaW*MF_new).*WF_new) + Kpred*Y*gamma')*(eye(length(gamma))/rho-(Kpred*(gamma*gamma'))/(rho^2+rho*Kpred*(gamma'*gamma)));
    
        L_new = rhoS*norm(G_new - XS*(sigm(alphaW*MS_new).*WS_new),"fro")^2 + rhoF*norm(G_new - XF*(sigm(alphaW*MF_new).*WF_new),"fro")^2 + ...
            lambda_fusion * (sum(sum(abs(sigm(alphaW*MS_new))))+sum(sum(abs(sigm(alphaW*MF_new))))) + ...
            Kpred * norm(Y-G_new*(sigm(alphaBeta*Mbeta_new).*beta_new),"fro")^2 + lambda_pred*sum(abs(sigm(alphaBeta*Mbeta_new)));

        L_history = [L_history,L_new];
        if L_new > L && iIter > 1
%             mu_fusion = mu_fusion/2;
%             mu_pred = mu_pred/2;
            mu_maskB = mu_maskB/2;
            mu_mask = mu_mask/2;
            disp('learning rate reduced')
        else    
            WS = WS_new;
            WF = WF_new;
            MS = MS_new;
            MF = MF_new;
            Mbeta = Mbeta_new;
            beta = beta_new;
            G = G_new;
        end
        L1 = rhoS*norm(G - XS*(sigm(alphaW*MS).*WS),"fro")^2 + rhoF*norm(G - XF*(sigm(alphaW*MF).*WF),"fro")^2;
        L2 = lambda_fusion * (sum(sum(abs(sigm(alphaW*MS))))+sum(sum(abs(sigm(alphaW*MF))))) ;
        L3 = Kpred * norm(Y-G*(sigm(alphaBeta*Mbeta).*beta),"fro")^2 ;
        L4 = lambda_pred*sum(abs(sigm(alphaBeta*Mbeta)));
        L_sub = [L_sub,[L1;L2;L3;L4]];
        
%         rsq = [rsq, evalRSQ(XS,WS,XF,WF,beta,Y,Kpred)];
    
    %     if mod(iIter,10000)==0
    %         save(['checkpoint_' num2str(iIter)],'WS','WF','G','beta','L_history')
    %         disp(['Iteration: ' num2str(iIter)])
    %     end
    
    %     if rsq > rsqThresh
    %         break
    %     end
    
    end
%     mu_mask = mu_mask * 10;
    iIter = 0;
end
% disp('for breakpoint')

end