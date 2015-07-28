function [lb,ub] = parconfest(simulate,model,net,xch,input,mea,var)
%% parallel confidence interval estimation of fluxes and for/rev ratios
% parameters for variable step size h
tol=1e-4; % maximum difference between true value and 1st order approximation
h_init=tol/2; % initial step size
h_coef1=0.998; % decreases step size by a multiplication factor (must be <= 1)
h_coef2=0.92; % decreases step size by a multiplication factor (must be <= 1)
thr=3.84; % 95%:3.84, 90%:2.71, dof=1; threshold value for confidence intervals

%% set variables and constants
% single vector for measured (and simulated) labeling patterns
[obs]=svec(mea);

% lower and upper bounds
kernel_net=model.kernel_net;
kernel_xch=model.kernel_xch;
fluxes=model.fluxes;
[opt_fr,lb_nx,N,all_free,cons_i,rm_cons]=consfree2full(kernel_net,kernel_xch,fluxes,net,xch);
ub_nx=lb_nx;

% Determine ratio of for/rev
opt_ratio=opt_fr(:,1)./opt_fr(:,2);

% constatns
N_index_l=2*size(lb_nx,1);
n_l=length(net);
nx_l=n_l+length(xch);
n_l_1=n_l+1;

% calculate optimum
[opt_cur,opt_score]=varssr(simulate,net,xch,input,mea,[],var);

% lb and ub of constrained fluxes do not need to be calculated
I=linspace(1,N_index_l,N_index_l);
I(cons_i)=[];
lb_s=zeros(length(I),1);
ub_s=lb_s;
N_I=N(I,:);

%% confidence interval: LOWERBOUND and UPPERBOUND
parfor ii=1:length(I)
    disp(['reaction ' num2str(I(ii)) ' in progress...'])
    all_free_par=all_free;
    Ni_all=N_I(ii,:);
    Ni=Ni_all(~rm_cons);
    
% matrices to store current lower bounds of all of the for and rev fluxes
% and their ratios (global updating) 
    lb_fr(:,:,ii) = opt_fr;
    lb_ratio(:,ii) = opt_ratio;
% matrices to store current upper bounds of all of the for and rev fluxes
% and their ratios (global updating) 
    ub_fr(:,:,ii) = opt_fr;
    ub_ratio(:,ii) = opt_ratio;
    
% LOWERBOUND
    lambda=0;
    h=h_init;
    d=[net; xch; lambda];
    all_free_par(~rm_cons)=d(1:nx_l);
    cur=opt_cur;
    score_new=opt_score;
    while score_new-opt_score<=thr && ~max(d(n_l_1:nx_l)<0)
        score=score_new; % need it to update h
        
        % Calculate current net or exchange lowerbound
        lb_s(ii)=Ni_all*all_free_par;
        
        % calculate sensitivity matrix (jacobian matrix)
        sm=numjac(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea);
        
        % calculate jacobian (jacobian vector) and hessian
        temp=sm'/var;
        J=temp*(cur-obs);
        H=temp*sm;
        
        % Calculate current for and rev lower bounds for all reactions
        fr_s=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:nx_l));
        
        % Update for and rev lower and upper bound matrix
        tmp=lb_fr(:,:,ii);
        I_s=fr_s<tmp;
        tmp(I_s)=fr_s(I_s);
        lb_fr(:,:,ii)=tmp;
        
        tmp=ub_fr(:,:,ii);
        I_s=fr_s>tmp;
        tmp(I_s)=fr_s(I_s);
        ub_fr(:,:,ii)=tmp;
        
        % Calculate for/rev ratios for all reactions
        ratio_s=fr_s(:,1)./fr_s(:,2);
        
        % Update if ratio is smaller or greater than current ratio
        tmp=lb_ratio(:,ii);
        I_s=ratio_s<tmp;
        tmp(I_s)=ratio_s(I_s);
        lb_ratio(:,ii)=tmp;
        
        tmp=ub_ratio(:,ii);
        I_s=ratio_s>tmp;
        tmp(I_s)=ratio_s(I_s);
        ub_ratio(:,ii)=tmp;
        
        % calculate A and b, solve for delta_d (dd), update d, get score
        A=[2*H -Ni'; -Ni 0];
        b=[-2*J; h];
        dd=A\b;
        d=d+dd;
        all_free_par(~rm_cons)=d(1:nx_l);
        [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
        
        % update step size h; if h is too large, go back to make it smaller
        diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        h=h_coef1*h*nthroot(tol/diff,3);
        while diff>tol
            d=d-dd;
            h=h_coef2*h;
            b=[-2*J; h];
            dd=A\b;
            d=d+dd;
            all_free_par(~rm_cons)=d(1:nx_l);
            [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
            diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        end
    end
    
% UPPERBOUND
    lambda=0;
    h=h_init;
    d=[net; xch; lambda];
    all_free_par(~rm_cons)=d(1:nx_l);
    cur=opt_cur;
    score_new=opt_score;
    while score_new-opt_score<=thr && ~max(d(n_l_1:nx_l)<0)
        score=score_new; % need it to update h
        
        % Calculate current net or exchange upperbound
        ub_s(ii)=Ni_all*all_free_par;
        
        % calculate sensitivity matrix (jacobian matrix)
        sm=numjac(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea);
        
        % calculate jacobian (jacobian vector) and hessian
        temp=sm'/var;
        J=temp*(cur-obs);
        H=temp*sm;
        
        % Calculate current for and rev lower bounds for all reactions
        fr_s=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:nx_l));
        
        % Update for and rev lower and upper bound matrix
        tmp=lb_fr(:,:,ii);
        I_s=fr_s<tmp;
        tmp(I_s)=fr_s(I_s);
        lb_fr(:,:,ii)=tmp;
        
        tmp=ub_fr(:,:,ii);
        I_s=fr_s>tmp;
        tmp(I_s)=fr_s(I_s);
        ub_fr(:,:,ii)=tmp;
        
        % Calculate for/rev ratios for all reactions
        ratio_s=fr_s(:,1)./fr_s(:,2);
        
        % Update if ratio is smaller or greater than current ratio
        tmp=lb_ratio(:,ii);
        I_s=ratio_s<tmp;
        tmp(I_s)=ratio_s(I_s);
        lb_ratio(:,ii)=tmp;
        
        tmp=ub_ratio(:,ii);
        I_s=ratio_s>tmp;
        tmp(I_s)=ratio_s(I_s);
        ub_ratio(:,ii)=tmp;
        
        % calculate A and b, solve for delta_d (dd), update d, get score
        A=[2*H Ni'; Ni 0];
        b=[-2*J; h];
        dd=A\b;
        d=d+dd;
        all_free_par(~rm_cons)=d(1:nx_l);
        [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
        
        % update step size h; if h is too large, go back to make it smaller
        diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        h=h_coef1*h*nthroot(tol/diff,3);
        while diff>tol
            d=d-dd;
            h=h_coef2*h;
            b=[-2*J; h];
            dd=A\b;
            d=d+dd;
            all_free_par(~rm_cons)=d(1:nx_l);
            [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
            diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        end
    end
end
lb_nx(I)=lb_s;
ub_nx(I)=ub_s;
lb_fr=min(lb_fr,[],3);
ub_fr=max(ub_fr,[],3);
lb_ratio=min(lb_ratio,[],2);
ub_ratio=max(ub_ratio,[],2);

% Concatenate and store lower and upper bounds
% (col 1=net, 2=xch, 3=for, 4=rev, 5=ratio)
lb = [lb_nx, lb_fr, lb_ratio];
ub = [ub_nx, ub_fr, ub_ratio];