function [net_opt,xch_opt,sm,J,H] = fluxest(simulate,net,xch,input,mea,var,iter)
%% gradient-descent method for the cumomer model to estimate optimal fluxes
% Input: simulate (function handle), net and xch (free flux vectors), input
% (input isotope substrate labeling patterns), mea and var (average
% measured labeling patterns and variances), and iter (maximum number of
% iterations)
%
% Output: net_opt and xch_opt (optimized free net and exchange fluxes)

% var can be covariance in the future

% tolerances are specific to the magnitude of fluxes
% sigma_tol=1e-10;%1e-10;
flux_tol=1e-4;%1e-6;
% sing_tol=1e-3;%singularity tolerance

net=net(:);
xch=xch(:);
net_l=length(net);
flux_l=net_l+length(xch);
net_ll=net_l+1;
free_net=net;
free_xch=xch;

%% create single vectors for measured (and simulated) labeling patterns
[obs]=svec(mea);

%% optimization steps
for ii=0:iter
%     calculate variance-weighted sum of squared residuals and simulate
%     labeling patterns with current fluxes
    [cur,score]=varssr(simulate,free_net,free_xch,input,mea,var);
    disp(['iter= ', num2str(ii) ': score= ', num2str(score)])
    
%     calculate sensitivity matrix (jacobian matrix)
    sm=numjac(simulate,free_net,free_xch,input,mea);
%     nondet=max(sm)<sing_tol;

%     calculate jacobian (jacobian vector) and (inverse) hessian
    temp=sm'/var;
    J=temp*(cur-obs);
    H=temp*sm;
%     J(nondet)=[]
%     H(nondet,:)=[];
%     H(:,nondet)=[]
%     [~,S,V]=svd(H);
%     Sd=diag(S)
%     zero_S=length(Sd(Sd<sigma_tol)); % # of close-to-zero eigenvalues
%     is the if statement necessary?
%     if Sd>0
%     v2=V(end-zero_S+1:end);
%     H_i=H+sigma_tol*v2*v2';
%     end
%     H_i=inv(H_i);
%     maybe use pinv(H) for inverse hessian?
%     H_i=pinv(H);
    
%     calculate flux step
%     dv=zeros(flux_l,1);
%     dv=-pinv(H)*J;
    dv=-H\J;
    
%     check for convergence
    if length(dv(abs(dv)<flux_tol))==flux_l
        break
    end
%     update fluxes
    free_net=free_net+dv(1:net_l);
    free_xch=free_xch+dv(net_ll:end);
end
% calculate all fluxes from the resulting free fluxes outside this function
if ii==iter
    warning('maximum iteration reached; may want to increase it')
end
net_opt=free_net;
xch_opt=free_xch;