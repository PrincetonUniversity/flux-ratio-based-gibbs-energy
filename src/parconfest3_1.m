function [lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I,lb_opt_nx,ub_opt_nx] = parconfest3_1(simulate,model,net,xch,ineq,input,mea,fmea,var,lb_in,ub_in,I_in,cache_in)
%% parallel confidence interval estimation of fluxes and for/rev ratios
% Dependencies: consfree2full.m, xlsreadineq.m, varssr.m, fluxest2.m

% Input: simulate (function handle), model (MATLAB structure containing
% metabolic model), net and xch (free flux vectors), ineq (inequality
% constraints), input (input isotope
% substrate labeling patterns), mea and var (average measured labeling
% patterns and variances), fmea (structure of flux measurements), lb_in,
% ub_in, and I_in (input lower bound, upper bound, and index)

% Output: lb and ub (contain net, xch, for, rev, for/rev 95% confidence
% intervals), net_nopt, xch_nopt, nscore (new optimum found during 
% parconfest3), nflag (adds 1 everytime a new optimum is found), hs (step
% size)

%% set variables and constants
thr=3.84; % 95%:3.84, 90%:2.71, dof=1; threshold value for confidence intervals
h_rel=0.01; % minimum step size relative to the flux of interest (previously 0.025)
h_con=2; % fold increase in step size if score changes by more than thr_inc
thr_inc1=0.25; % minimum increment in the score
thr_inc2=0.75; % maximum increment in the score
i_tol=15; %10 for large models
hcoef=0.01;
% tt=datestr(clock,'YYmmDDHHMMSS');

% lower and upper bounds
kernel_net=model.kernel_net;
kernel_xch=model.kernel_xch;
fluxes=model.fluxes;
[lb_opt_fr,lb_opt_nx,N,all_free,~,rm_cons,free_i]=consfree2full(kernel_net,kernel_xch,fluxes,net,xch);

% constants
numRxns=size(lb_opt_nx, 1);
N_index_l=2*numRxns;
n_l=length(net);
n_l_1=n_l+1;
if isempty(ineq)
    [ineq]=xlsreadineq([],model,net,xch);
end

% calculate optimum
[~,opt_score]=varssr(simulate,net,xch,input,mea,fmea,var);

% find in/output fluxes and stop when they are less than 0
iofluxes=fluxes(:,2)==0 & fluxes(:,10)==1;

% if lb, ub, and I are input
cstatus=[];
if nargin>9 && nargin<13 && ~isempty(lb_in) && ~isempty(ub_in) && ~isempty(I_in)
    lb_opt_nx=lb_in(:,1:2);
    ub_opt_nx=ub_in(:,1:2);
    lb_opt_fr=lb_in(:,3:4);
    ub_opt_fr=ub_in(:,3:4);
    lb_opt_ratio=lb_in(:,5);
    ub_opt_ratio=ub_in(:,5);
    I=I_in;
    In=[];
    Ix=[];
elseif nargin>9 && nargin<13 && isempty(lb_in) && isempty(ub_in) && ~isempty(I_in)
    error('Yet to be defined: 1');
elseif nargin==13 && isempty(lb_in) && isempty(ub_in) && ~isempty(I_in)
    error('Yet to be defined: 2');
elseif nargin==13 && isempty(lb_in) && isempty(ub_in) && isempty(I_in) && ~isempty(cache_in)
    cin=load(cache_in);
    I=cin.I;
    l_I=length(I);
    cstatus=zeros(l_I,1);
    for ii=1:l_I
        if isfield(cin,['status_' num2str(I(ii))])
            eval(['cstatus(' num2str(ii) ')=cin.status_' num2str(I(ii)) ';']);
        end
    end
    ub_opt_nx=lb_opt_nx;
    ub_opt_fr=lb_opt_fr;
    lb_opt_ratio=lb_opt_fr(:,1)./lb_opt_fr(:,2);
    ub_opt_ratio=ub_opt_fr(:,1)./ub_opt_fr(:,2);
    In=I(I<=numRxns);
    Ix=I(I>numRxns)-numRxns;
else
    ub_opt_nx=lb_opt_nx;
    ub_opt_fr=lb_opt_fr;
    lb_opt_ratio=lb_opt_fr(:,1)./lb_opt_fr(:,2);
    ub_opt_ratio=ub_opt_fr(:,1)./ub_opt_fr(:,2);
%     lb and ub of quasi/constrained fluxes do not need to be calculated
    I=linspace(1,N_index_l,N_index_l);
    I(~free_i)=[];
    In=I(I<=numRxns);
    Ix=I(I>numRxns)-numRxns;
    % I(I>numRxns)=[]; % uncomment when only net confidence intervals are needed %
end

% parallel variables
lb_s=zeros(length(I),1);
ub_s=lb_s;
N_I=N(I,:);
lb_nx=inf([size(lb_opt_nx) length(I)]);
lb_fr=inf([size(lb_opt_fr) length(I)]);
lb_ratio=inf(length(lb_opt_fr), length(I));
ub_nx=-lb_nx;
ub_fr=-lb_fr;
ub_ratio=-lb_ratio;
net_nopt=zeros(length(net), length(I));
xch_nopt=zeros(length(xch), length(I));
nscore=opt_score*ones(length(I),1);
nflag=zeros(length(I),1);
hs=zeros(N_index_l,1);
hs_lb=inf(length(I),1);
hs_ub=-hs_lb;

if ~isempty(cstatus)
    for ii=1:l_I
        if cstatus(ii)==1
            eval(['lb_s(ii) = cin.lbs_' num2str(I(ii)) ';']);
            eval(['lb_nx(:,:,ii) = cin.lbnx_' num2str(I(ii)) ';']);
            eval(['ub_nx(:,:,ii) = cin.ubnx_' num2str(I(ii)) ';']);
            eval(['lb_fr(:,:,ii) = cin.lbfr_' num2str(I(ii)) ';']);
            eval(['ub_fr(:,:,ii) = cin.ubfr_' num2str(I(ii)) ';']);
            eval(['lb_ratio(:,ii) = cin.lbratio_' num2str(I(ii)) ';']);
            eval(['ub_ratio(:,ii) = cin.ubratio_' num2str(I(ii)) ';']);
        elseif cstatus(ii)==2
            eval(['lb_s(ii) = cin.lbs_' num2str(I(ii)) ';']);
            eval(['ub_s(ii) = cin.ubs_' num2str(I(ii)) ';']);
            eval(['lb_nx(:,:,ii) = cin.lbnx_' num2str(I(ii)) ';']);
            eval(['ub_nx(:,:,ii) = cin.ubnx_' num2str(I(ii)) ';']);
            eval(['lb_fr(:,:,ii) = cin.lbfr_' num2str(I(ii)) ';']);
            eval(['ub_fr(:,:,ii) = cin.ubfr_' num2str(I(ii)) ';']);
            eval(['lb_ratio(:,ii) = cin.lbratio_' num2str(I(ii)) ';']);
            eval(['ub_ratio(:,ii) = cin.ubratio_' num2str(I(ii)) ';']);
        end
    end
end
%% cache file to divide up workload
h=datestr(clock,0);
cacheFile=['cache_',h(1:11),'_',h(13:14),'_',h(16:17),'_',h(19:20)];
fullcacheFile=[pwd '/' cacheFile '.mat'];
if exist(fullcacheFile,'file')==0
    save(cacheFile,'h','I','simulate');
else
    warning('Overriding cache file');
end

%% confidence interval: LOWERBOUND and UPPERBOUND
parfor ii=1:length(I)
    if isempty(cstatus) || cstatus(ii)==0
        status=0;
    % matrices to store current lower bounds of all of the for and rev fluxes
    % and their ratios (global updating) 
        lb_fr(:,:,ii) = lb_opt_fr;
        lb_ratio(:,ii) = lb_opt_ratio;
    % matrices to store current upper bounds of all of the for and rev fluxes
    % and their ratios (global updating) 
        ub_fr(:,:,ii) = ub_opt_fr;
        ub_ratio(:,ii) = ub_opt_ratio;
    elseif cstatus(ii)==2
        disp([num2str(I(ii)) ': skip because it was cached'])
        continue
    else
        status=cstatus(ii);
    end

% matrices to store current optimum info
    net_nopt(:,ii) = net;
    xch_nopt(:,ii) = xch;

% skip duplicates
    skip=0;
    ff=I(ii);
    if ff>numRxns
        ffx=ff-numRxns;
        for fff=Ix(Ix<ffx)
            if kernel_xch(fff,:)==kernel_xch(ffx,:)
                skip=1;
                break
            end
        end
    else
        for fff=In(In<ff)
            if kernel_net(fff,:)==kernel_net(ff,:)
                skip=1;
                break
            end
        end
    end
    if skip==1
        disp([num2str(I(ii)) ': skip'])
        continue
    end
    
    disp([num2str(I(ii))  ': optimum = ' num2str(opt_score)])
    all_free_par=all_free;
    Ni_all=N_I(ii,:);
    eq=struct;
% save current results to cache
    var_lb_s=['lbs_' num2str(I(ii))];
    var_ub_s=['ubs_' num2str(I(ii))];
    var_lb_nx=['lbnx_' num2str(I(ii))];
    var_ub_nx=['ubnx_' num2str(I(ii))];
    var_lb_fr=['lbfr_' num2str(I(ii))];
    var_ub_fr=['ubfr_' num2str(I(ii))];
    var_lb_ratio=['lbratio_' num2str(I(ii))];
    var_ub_ratio=['ubratio_' num2str(I(ii))];
    var_status=['status_' num2str(I(ii))];
    
    if status==0
% LOWERBOUND
    disp([num2str(I(ii)) ': LB in progress...'])
    it=0;
    d=[net; xch];
    all_free_par(~rm_cons)=d;
    eq.A=Ni_all(~rm_cons);
    eq.b=eq.A*d;
    h=abs(h_rel*Ni_all*all_free_par);
    if h<2e-6
        h=5e-5;
    end
    info.fval=nscore(ii);
    while info.fval-nscore(ii)<=thr && ~max(d(n_l_1:end)<0) && it<i_tol && min(ineq.A*d<=ineq.b)
        % Calculate current for and rev lower bounds for all reactions
        [fr_s,nx_s]=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:end));
        
        % Stop when in/output fluxes are less than 0
        if max(fr_s(iofluxes,1)<0)
            fr_s_0I=find(fr_s(iofluxes,1)<0);
            disp([num2str(I(ii)) ': LB for or rev < 0 at: ' num2str(fr_s_0I)])
            break
        end
        
        % Update net and xch lower and upper bound matrix
        tmp=lb_nx(:,:,ii);
        I_s=nx_s<tmp;
        tmp(I_s)=nx_s(I_s);
        lb_nx(:,:,ii)=tmp;
        
        tmp=ub_nx(:,:,ii);
        I_s=nx_s>tmp;
        tmp(I_s)=nx_s(I_s);
        ub_nx(:,:,ii)=tmp;
        
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
        
        % individual lowerbound
        lb_s(ii)=Ni_all*all_free_par;
        
        % update free fluxes
%         dd=lsqlin(eq.A,h,ineq.A,ineq.b);
%         d=d-dd;
        d=d-eq.A\h;
        eq.b=eq.A*d;
        if ~min(ineq.A*d<=ineq.b) % eq.b<0 && ff>numRxns
            nhcoef=0;
            ob=find(ineq.A*d>ineq.b);
            disp([num2str(I(ii)) ': LB flux approaching ineq boundary at ' num2str(ob')])
            if sum(sum(ineq.A(ob,:)~=0,2)==1)>0
                skip=skip+1;
                while ~min(ineq.A*d<=ineq.b) && nhcoef<1/hcoef-1
%                     d=d+hcoef*dd;
                    d=d+hcoef*(eq.A\h);
                    eq.b=eq.A*d;
                    nhcoef=nhcoef+1;
                end
            end
            if skip==2
                break
            end
        end
%         disp([num2str(I(ii)) ': ineq = ', num2str(min(ineq.A*d<=ineq.b)), '; eq = ', num2str(min(eq.A*d==eq.b))])
        % bypass errors in fmincon
        try
            fprintf('%d: LB ',I(ii));
            info_fval=info.fval;
            [d(1:n_l),d(n_l_1:end),info]=fluxest2(simulate,model,d(1:n_l),d(n_l_1:end),ineq,eq,input,mea,fmea,var,info_fval);
            disp([num2str(I(ii)) ': LB optimization finishes with: ' num2str(info.fval), '; ineq = ', num2str(min(ineq.A*d<=ineq.b)), '; eq = ', num2str(min(eq.A*d==eq.b)), '; step = ', num2str(h), '; it = ', num2str(it), '; ', num2str(eq.b)])
            if nscore(ii)>info.fval && min(ineq.A*d<=ineq.b) %&& min(eq.A*d==eq.b)
                net_nopt(:,ii)=d(1:n_l);
                xch_nopt(:,ii)=d(n_l_1:end);
                nscore(ii)=info.fval;
                nflag(ii)=nflag(ii)+1;
                net_nopt(:,ii)
                xch_nopt(:,ii)
                nscore(ii)
                nflag(ii)
                disp([num2str(I(ii)) ': LB new optimum is found: ' num2str(nscore(ii)), '; ineq = ', num2str(min(ineq.A*d<=ineq.b)), '; eq = ', num2str(min(eq.A*d==eq.b))])
            end
        catch
            disp([num2str(I(ii)) ': LB has error in fmincon'])
            continue
        end
        all_free_par(~rm_cons)=d;
        it=it+1;
        hs_lb(ii)=h;
        if info.fval-info_fval<thr_inc1 && it~=i_tol
            h=h_con*h;
            disp([num2str(I(ii)) ': LB step size is larger; ', num2str(h)])
        elseif info.fval-info_fval>thr_inc2 && it~=i_tol
            h=h/h_con;
            disp([num2str(I(ii)) ': LB step size is smaller; ', num2str(h)])
        end
        if it==i_tol && info.fval-nscore(ii)<=thr
            disp([num2str(I(ii)) ': LB can be much smaller'])
        end
    end
    status=status+1;
    saveappend(cacheFile,var_lb_s,var_lb_nx,var_ub_nx,var_lb_fr,var_ub_fr,var_lb_ratio,var_ub_ratio,lb_s(ii),lb_nx(:,:,ii),ub_nx(:,:,ii),lb_fr(:,:,ii),ub_fr(:,:,ii),lb_ratio(:,ii),ub_ratio(:,ii));
    saveappend1(cacheFile,var_status,status);
    disp([num2str(I(ii)) ': LB done'])
    end
    if status==1
% UPPERBOUND
    skip=0;
    disp([num2str(I(ii)) ': UB in progress...'])
    it=0;
    d=[net; xch];
    all_free_par(~rm_cons)=d;
    eq.A=Ni_all(~rm_cons);
    eq.b=eq.A*d;
    h=abs(h_rel*Ni_all*all_free_par);
    if h<2e-6
        h=5e-5;
    end
    info.fval=nscore(ii);
    while info.fval-nscore(ii)<=thr && ~max(d(n_l_1:end)<0) && it<i_tol && min(ineq.A*d<=ineq.b)
        % Calculate current for and rev lower bounds for all reactions
        [fr_s,nx_s]=consfree2full(kernel_net,kernel_xch,fluxes,d(1:n_l),d(n_l_1:end));
        
        % Stop when in/output fluxes are less than 0
        if max(fr_s(iofluxes,1)<0)
            fr_s_0I=find(fr_s(iofluxes,1)<0);
            disp([num2str(I(ii)) ': UB for or rev < 0 at: ' num2str(fr_s_0I)])
            break
        end
        
        % Update for and rev lower and upper bound matrix
        tmp=lb_nx(:,:,ii);
        I_s=nx_s<tmp;
        tmp(I_s)=nx_s(I_s);
        lb_nx(:,:,ii)=tmp;
        
        tmp=ub_nx(:,:,ii);
        I_s=nx_s>tmp;
        tmp(I_s)=nx_s(I_s);
        ub_nx(:,:,ii)=tmp;
        
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
        
        % individual upperbound
        ub_s(ii)=Ni_all*all_free_par;
        
        % update free fluxes
%         dd=lsqlin(eq.A,h,ineq.A,ineq.b);
%         d=d+dd;
        d=d+eq.A\h;
        eq.b=eq.A*d;
        if ~min(ineq.A*d<=ineq.b)
            nhcoef=0;
            ob=find(ineq.A*d>ineq.b);
            disp([num2str(I(ii)) ': UB flux approaching ineq boundary at ' num2str(ob')])
            if sum(sum(ineq.A(ob,:)~=0,2)==1)>0
                skip=skip+1;
                while ~min(ineq.A*d<=ineq.b) && nhcoef<1/hcoef-1
%                     d=d-hcoef*dd;
                    d=d-hcoef*(eq.A\h);
                    eq.b=eq.A*d;
                    nhcoef=nhcoef+1;
                end
            end
            if skip==2
                break
            end
        end
        
        % bypass errors in fmincon
        try
            fprintf('%d: UB ',I(ii));
            info_fval=info.fval;
            [d(1:n_l),d(n_l_1:end),info]=fluxest2(simulate,model,d(1:n_l),d(n_l_1:end),ineq,eq,input,mea,fmea,var,info_fval);
            disp([num2str(I(ii)) ': UB optimization finishes with: ' num2str(info.fval), '; ineq = ', num2str(min(ineq.A*d<=ineq.b)), '; eq = ', num2str(min(eq.A*d==eq.b)), '; step = ', num2str(h), '; it = ', num2str(it), '; ', num2str(eq.b)])
            if nscore(ii)>info.fval && min(ineq.A*d<=ineq.b) %&& min(eq.A*d==eq.b)
                net_nopt(:,ii)=d(1:n_l);
                xch_nopt(:,ii)=d(n_l_1:end);
                nscore(ii)=info.fval;
                nflag(ii)=nflag(ii)+1;
                net_nopt(:,ii)
                xch_nopt(:,ii)
                nscore(ii)
                nflag(ii)
                disp([num2str(I(ii)) ': UB new optimum is found: ' num2str(nscore(ii)), '; ineq = ', num2str(min(ineq.A*d<=ineq.b)), '; eq = ', num2str(min(eq.A*d==eq.b))])
            end
        catch
            disp([num2str(I(ii)) ': UB has error in fmincon'])
            continue
        end
        all_free_par(~rm_cons)=d;
        it=it+1;
        hs_ub(ii)=h;
        if info.fval-info_fval<thr_inc1 && it~=i_tol
            h=h_con*h;
            disp([num2str(I(ii)) ': UB step size is larger; ', num2str(h)])
        elseif info.fval-info_fval>thr_inc2 && it~=i_tol
            h=h/h_con;
            disp([num2str(I(ii)) ': UB step size is smaller; ', num2str(h)])
        end
        if it==i_tol && info.fval-nscore(ii)<=thr
            disp([num2str(I(ii)) ': UB can be much larger'])
        end
    end
    status=status+1;
    saveappend(cacheFile,var_ub_s,var_lb_nx,var_ub_nx,var_lb_fr,var_ub_fr,var_lb_ratio,var_ub_ratio,ub_s(ii),lb_nx(:,:,ii),ub_nx(:,:,ii),lb_fr(:,:,ii),ub_fr(:,:,ii),lb_ratio(:,ii),ub_ratio(:,ii));
    saveappend1(cacheFile,var_status,status);
    disp([num2str(I(ii)) ': UB done'])
    end
end
lb_opt_nx(I)=lb_s;
ub_opt_nx(I)=ub_s;
lb_nx=min(lb_nx,[],3);
ub_nx=max(ub_nx,[],3);
lb_fr=min(lb_fr,[],3);
ub_fr=max(ub_fr,[],3);
lb_ratio=min(lb_ratio,[],2);
ub_ratio=max(ub_ratio,[],2);

hs(I,1)=hs_lb;
hs(I,2)=hs_ub;

% Concatenate and store lower and upper bounds
% (col 1=net, 2=xch, 3=for, 4=rev, 5=ratio)
lb = [lb_nx, lb_fr, lb_ratio];
ub = [ub_nx, ub_fr, ub_ratio];

if sum(nflag)~=0
    warning(['New optimum is found: ' num2str(min(nscore)) '; Total nflag is: ' num2str(sum(nflag))])
end

function saveappend1(fname,varstatus,status)
eval([varstatus '=status;']);
save(fname,varstatus,'-append')

function saveappend(fname,var_lub_s,var_lb_nx,var_ub_nx,var_lb_fr,var_ub_fr,var_lb_ratio,var_ub_ratio,lub_s,lb_nx,ub_nx,lb_fr,ub_fr,lb_ratio,ub_ratio)
eval([var_lub_s '=lub_s;']);
eval([var_lb_nx '=lb_nx;']);
eval([var_ub_nx '=ub_nx;']);
eval([var_lb_fr '=lb_fr;']);
eval([var_ub_fr '=ub_fr;']);
eval([var_lb_ratio '=lb_ratio;']);
eval([var_ub_ratio '=ub_ratio;']);
save(fname,var_lub_s,var_lb_fr,var_ub_fr,var_lb_ratio,var_ub_ratio,'-append')
