function [lb,ub] = confest(simulate,model,net,xch,input,mea,var)
%% confidence interval estimation of fluxes
%lb and ub are matrices with columns 1-5 corresponding to net, xch, for, rev, and for/rev, respectively

% parameters for variable step size h
tol=1e-4; % maximum difference between true value and 1st order approximation
h_init=tol/2; % initial step size
h_coef1=0.998; % decreases step size by a multiplication factor (must be <= 1)
h_coef2=0.92; % decreases step size by a multiplication factor (must be <= 1)
thr=3.84; % 95%:3.84, 90%:2.71, dof=1; threshold value for confidence intervals
% i_tol=0;
%% set variables and constants
% single vector for measured (and simulated) labeling patterns
[obs]=svec(mea);

% lower and upper bounds; convert to for and rev and store in lb_fr
kernel_net=model.kernel_net;
kernel_xch=model.kernel_xch;
fluxes=model.fluxes;
[opt_fr,lb_nx,N,all_free,cons_i,rm_cons]=consfree2full(kernel_net,kernel_xch,fluxes,net,xch);
ub_nx=lb_nx;

% Determine ratio of for/rev
opt_ratio=opt_fr(:,1)./opt_fr(:,2);

% constants
numRxns=size(lb_nx, 1);
N_index_l=2*numRxns;
n_l=length(net);
nx_l=n_l+length(xch);
n_l_1=n_l+1;

% calculate optimum
[opt_cur,opt_score]=varssr(simulate,net,xch,input,mea,[],var);

% lb and ub of quasi/constrained fluxes do not need to be calculated
I=linspace(1,N_index_l,N_index_l);
I(cons_i)=[];

%% confidence interval: LOWERBOUND and UPPERBOUND
% matrices to store current lower bounds of all of the for and rev fluxes
% and their ratios (global updating) 
lb_fr = opt_fr;
lb_ratio = opt_ratio;
% matrices to store current upper bounds of all of the for and rev fluxes
% and their ratios (global updating) 
ub_fr = opt_fr;
ub_ratio = opt_ratio;

for ff=I
    disp(['reaction ' num2str(ff) ' in progress...'])
    Ni_all=N(ff,:);
    Ni=Ni_all(~rm_cons);
	
% LOWERBOUND
    lambda=0;
    h=h_init;
    d=[net; xch; lambda];
    all_free(~rm_cons)=d(1:nx_l);
    cur=opt_cur;
    score_new=opt_score;
    % terminate when exceeding threshold or any xch flux is less than 0
    while score_new-opt_score<=thr && ~max(d(n_l_1:nx_l)<0)
        score=score_new; % need it to update h
        
        % Calculate current net or exchange lowerbound
        lb_s=Ni_all*all_free;
        
        % calculate sensitivity matrix (jacobian matrix)
        sm=numjac(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea);

        % calculate jacobian (jacobian vector) and hessian
        temp=sm'/var;
        J=temp*(cur-obs);
        H=temp*sm;
        
        % Calculate current for and rev lower bounds for all reactions
        fr_s=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:nx_l));
        
        % Update for and rev lower and upper bound matrix
        I_s=fr_s<lb_fr;
        lb_fr(I_s)=fr_s(I_s);
        I_s=fr_s>ub_fr;
        ub_fr(I_s)=fr_s(I_s);
        
        % Calculate for/rev ratios for all reactions
        ratio_s=fr_s(:,1)./fr_s(:,2);
        
        % Update if ratio is smaller or greater than current ratio
        I_s=ratio_s<lb_ratio;
        lb_ratio(I_s)=ratio_s(I_s);
        I_s=ratio_s>ub_ratio;
        ub_ratio(I_s)=ratio_s(I_s);
        
        % calculate A and b, solve for delta_d (dd), update d, get score
        A=[2*H -Ni'; -Ni 0];
        b=[-2*J; h];
        dd=A\b;
        d=d+dd;
        all_free(~rm_cons)=d(1:nx_l);
        [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
        
%         update step size h; if h is too large, go back to make it smaller
        diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        h=h_coef1*h*nthroot(tol/diff,3);
%         if h<h_min
%             display(['h is too small: ' num2str(h) '; chi-squred: ' num2str(score_new-opt_score) '; min. flux: ' num2str(min(N*all_free)) '; diff= ' num2str(diff)]);
%             h=h/h_coef1;
%         end
        while diff>tol
%             i_tol=i_tol+1;
%             display(['warning: tol_i= ' num2str(i_tol) '; h= ' num2str(h) '; diff= ' num2str(diff)])
%             if score_new-opt_score>thr || max(N*all_free<0)
%                 break
%             end
            d=d-dd;
            h=h_coef2*h;
            b=[-2*J; h];
            dd=A\b;
            d=d+dd;
            all_free(~rm_cons)=d(1:nx_l);
            [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
            diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
              
        end
    end
    lb_nx(ff)=lb_s;
        
% UPPERBOUND
    lambda=0;
    h=h_init;
    d=[net; xch; lambda];
    all_free(~rm_cons)=d(1:nx_l);
    cur=opt_cur;
    score_new=opt_score;
    % terminate when exceeding threshold or any xch flux is less than 0
    while score_new-opt_score<=thr && ~max(d(n_l_1:nx_l)<0)
        score=score_new; % need it to update h
        
        % Calculate current net or exchange upperbound
        ub_s=Ni_all*all_free;
        
        % calculate sensitivity matrix (jacobian matrix)
        sm=numjac(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea);
        
        % calculate jacobian (jacobian vector) and hessian
        temp=sm'/var;
        J=temp*(cur-obs);
        H=temp*sm;
        
        % Calculate current for and rev upper bounds for all reactions
        fr_s=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:nx_l));
        
        % Update for and rev lower and upper bound matrix
        I_s=fr_s<lb_fr;
        lb_fr(I_s)=fr_s(I_s);
        I_s=fr_s>ub_fr;
        ub_fr(I_s)=fr_s(I_s);
        
        % Calculate for/rev ratios for all reactions
        ratio_s=fr_s(:,1)./fr_s(:,2);
        
        % Update if ratio is smaller or greater than current ratio
        I_s=ratio_s<lb_ratio;
        lb_ratio(I_s)=ratio_s(I_s);
        I_s=ratio_s>ub_ratio;
        ub_ratio(I_s)=ratio_s(I_s);
		
        % calculate A and b, solve for delta_d (dd), update d, get score
        A=[2*H Ni'; Ni 0];
        b=[-2*J; h];
        dd=A\b;
        d=d+dd;
        all_free(~rm_cons)=d(1:nx_l);
        [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
        
%         update step size h; if h is too large, go back to make it smaller
        diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        h=h_coef1*h*nthroot(tol/diff,3);
%         if h<h_min
%             display(['h is too small: ' num2str(h) '; chi-squred: ' num2str(score_new-opt_score) '; min. flux: ' num2str(min(N*all_free)) '; diff= ' num2str(diff)]);
%             h=h/h_coef1;
%         end
        while diff>tol
%             i_tol=i_tol+1;
%             display(['warning: tol_i= ' num2str(i_tol) '; h= ' num2str(h) '; diff= ' num2str(diff)])
%             if score_new-opt_score>thr || max(N*all_free<0)
%                 break
%             end
            d=d-dd;
            h=h_coef2*h;
            b=[-2*J; h];
            dd=A\b;
            d=d+dd;
            all_free(~rm_cons)=d(1:nx_l);
            [cur,score_new]=varssr(simulate,d(1:n_l),d(n_l_1:nx_l),input,mea,[],var);
            diff=abs(score_new-score-2*dd(1:nx_l)'*J-dd(1:nx_l)'*H*dd(1:nx_l));
        end
    end
    ub_nx(ff)=ub_s;
end

% Concatenate and store lower and upper bounds
% (col 1=net, 2=xch, 3=for, 4=rev, 5=ratio)
lb = [lb_nx, lb_fr, lb_ratio];
ub = [ub_nx, ub_fr, ub_ratio];