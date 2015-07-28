function [lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I] = cache2argout(simulate,model,net,xch,ineq,input,mea,fmea,var,lb_in,ub_in,I_in,cache_in)
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
% lower and upper bounds
kernel_net=model.kernel_net;
kernel_xch=model.kernel_xch;
fluxes=model.fluxes;
[lb_opt_fr,lb_nx,~,~,cons_i,~]=consfree2full(kernel_net,kernel_xch,fluxes,net,xch);

% constants
numRxns=size(lb_nx, 1);
N_index_l=2*numRxns;
if isempty(ineq)
    [ineq]=xlsreadineq([],model,net,xch);
end

% calculate optimum
[~,opt_score]=varssr(simulate,net,xch,input,mea,fmea,var);

% if lb, ub, and I are input
cstatus=[];
if nargin>9 && nargin<13 && ~isempty(lb_in) && ~isempty(ub_in) && ~isempty(I_in)
    lb_nx=lb_in(:,1:2);
    ub_nx=ub_in(:,1:2);
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
    ub_nx=lb_nx;
    ub_opt_fr=lb_opt_fr;
    lb_opt_ratio=lb_opt_fr(:,1)./lb_opt_fr(:,2);
    ub_opt_ratio=ub_opt_fr(:,1)./ub_opt_fr(:,2);
    In=I(I<=numRxns);
    Ix=I(I>numRxns)-numRxns;
else
    ub_nx=lb_nx;
    ub_opt_fr=lb_opt_fr;
    lb_opt_ratio=lb_opt_fr(:,1)./lb_opt_fr(:,2);
    ub_opt_ratio=ub_opt_fr(:,1)./ub_opt_fr(:,2);
%     lb and ub of quasi/constrained fluxes do not need to be calculated
    I=linspace(1,N_index_l,N_index_l);
    I(cons_i)=[];
    In=I(I<=numRxns);
    Ix=I(I>numRxns)-numRxns;
    % I(I>numRxns)=[]; % uncomment when only net confidence intervals are needed %
end

% parallel variables
lb_s=zeros(length(I),1);
ub_s=lb_s;
lb_fr=inf([size(lb_opt_fr) length(I)]);
lb_ratio=inf(length(lb_opt_fr), length(I));
ub_fr=-lb_fr;
ub_ratio=-lb_ratio;
net_nopt=zeros(length(net), length(I));
xch_nopt=zeros(length(xch), length(I));
nscore=opt_score*ones(length(I),1);
nflag=zeros(length(I),1);

if ~isempty(cstatus)
    for ii=1:l_I
        if cstatus(ii)==1
            eval(['lb_s(ii) = cin.lbs_' num2str(I(ii)) ';']);
            eval(['lb_fr(:,:,ii) = cin.lbfr_' num2str(I(ii)) ';']);
            eval(['ub_fr(:,:,ii) = cin.ubfr_' num2str(I(ii)) ';']);
            eval(['lb_ratio(:,ii) = cin.lbratio_' num2str(I(ii)) ';']);
            eval(['ub_ratio(:,ii) = cin.ubratio_' num2str(I(ii)) ';']);
        elseif cstatus(ii)==2
            eval(['lb_s(ii) = cin.lbs_' num2str(I(ii)) ';']);
            eval(['ub_s(ii) = cin.ubs_' num2str(I(ii)) ';']);
            eval(['lb_fr(:,:,ii) = cin.lbfr_' num2str(I(ii)) ';']);
            eval(['ub_fr(:,:,ii) = cin.ubfr_' num2str(I(ii)) ';']);
            eval(['lb_ratio(:,ii) = cin.lbratio_' num2str(I(ii)) ';']);
            eval(['ub_ratio(:,ii) = cin.ubratio_' num2str(I(ii)) ';']);
        end
    end
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
    
    if status==0
% LOWERBOUND
        disp([num2str(I(ii)) ': LB needs to run...'])
        status=status+1;
    end
    if status==1
% UPPERBOUND
        disp([num2str(I(ii)) ': UB needs to run...'])
        status=status+1;
    end
end
lb_nx(I)=lb_s;
ub_nx(I)=ub_s;
lb_fr=min(lb_fr,[],3);
ub_fr=max(ub_fr,[],3);
lb_ratio=min(lb_ratio,[],2);
ub_ratio=max(ub_ratio,[],2);

hs=0;

% Concatenate and store lower and upper bounds
% (col 1=net, 2=xch, 3=for, 4=rev, 5=ratio)
lb = [lb_nx, lb_fr, lb_ratio];
ub = [ub_nx, ub_fr, ub_ratio];

