function [net_opt,xch_opt,info] = fluxest2(simulate,model,net,xch,ineq,eq,input,mea,fmea,var,prevopt,CIter,keeplast)
%% gradient-descent method for the cumomer model to estimate optimal fluxes
% Dependencies: xlsreadineq.m, varssr.m
% input: simulate (function handle), model (MATLAB structure containing
% metabolic model), net and xch (free flux vectors), input (input isotope
% substrate labeling patterns), mea and var (average measured labeling
% patterns and variances), fmea (structure of flux measurements), 
% iter (maximum number of iterations)
% output: net_opt and xch_opt (optimized free net and exchange fluxes),
% info (information regarding fit)

% reshape net and xch matrices into column vectors
net=net(:);
xch=xch(:);

net_l=length(net);
% flux_l=net_l+length(xch);
net_ll=net_l+1;

F=@(x)objfunc(simulate,x(1:net_l),x(net_ll:end),input,mea,fmea,var);

%% optimization steps

% when no inequality constraints are specified
if isempty(ineq)
    [ineq]=xlsreadineq([],model,net,xch);
end

% for flux estimation
fluxest=0;
if isempty(eq)
    eq.A=zeros(1,net_l+length(xch));
    eq.b=0;
    fluxest=1;
end

% initialize
init_score=F([net;xch]);
init_fea=[min(ineq.A*[net;xch]<=ineq.b) min(eq.A*[net;xch]==eq.b)];
disp(['optimization begins with score = ', num2str(init_score), '; ineq = ', num2str(init_fea(1)), '; eq = ', num2str(init_fea(2))])

if nargin<11 || isempty(prevopt)
    prevopt=init_score-1;
end
if nargin<12 || isempty(CIter)
    CIter=500;
end
if nargin<13 || isempty(keeplast)
    keeplast=0;
end

% if init_score is less than the opt score, shorten the optimization
if sum(init_fea)==2 && init_score<=prevopt
    CIter=50;
    disp(['MaxIter = ', num2str(CIter)])
end

if fluxest
% for flux estimation
    options=optimset('Algorithm','interior-point' ...
                    ,'Display','iter-detailed' ...
                    ,'TolX',1e-10,'TolCon',1e-30,'TolFun',1e-7 ...
                    ,'MaxFunEvals',50000,'MaxIter',500);
else
% for confidence interval estimation
    options=optimset('Algorithm','interior-point' ...
                ,'Display','notify-detailed' ...
                ,'TolX',1e-10,'TolCon',1e-30,'TolFun',1e-7 ...
                ,'MaxFunEvals',100*CIter,'MaxIter',CIter);
end

% MultiStart, ga, GlobalSearch
while isnan(init_score)
    net=-2.5+5*rand(size(net));
    xch=5*rand(size(xch));
    init_score=F([net;xch]);
    init_fea=[min(ineq.A*[net;xch]<=ineq.b) min(eq.A*[net;xch]==eq.b)];
    disp(['optimization begins with score = ', num2str(init_score), '; ineq = ', num2str(init_fea(1)), '; eq = ', num2str(init_fea(2))])
end
% options=optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% [info.x,info.fval,info.exitflag,info.output]=fminsearch(F,[net;xch],options)


% turn off unnecessary warnings
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:illConditionedMatrix')
try
    [info.x,info.fval,info.exitflag,info.output]=fmincon(F,[net;xch],ineq.A,ineq.b,eq.A,eq.b,[],[],[],options);
catch
    warning('error in fmincon')
    % if there's error, try a new flux set with fewer iterations
    options=optimset('Algorithm','interior-point' ...
                    ,'Display','iter-detailed' ...
                    ,'MaxFunEvals',30000,'MaxIter',300);
    net=-1.5+3*rand(size(net));
    xch=3*rand(size(xch));
    [info.x,info.fval,info.exitflag,info.output]=fmincon(F,[net;xch],ineq.A,ineq.b,eq.A,eq.b,[],[],[],options);
end

opt_fea=[min(ineq.A*info.x<=ineq.b) min(eq.A*info.x==eq.b)];
if sum(opt_fea)<2
    warning(['infeasible opt: init = ' num2str(init_score), '; opt = ' num2str(info.fval), '; ineq = ', num2str(opt_fea(1)), '; eq = ', num2str(opt_fea(2))])
    disp(['eq.A*[net;xch] = ', num2str(eq.A*info.x), '; eq.b = ', num2str(eq.b)])
end

if sum(init_fea)==2 && init_score<info.fval && keeplast==0
    disp(['better init: init = ' num2str(init_score), '; opt = ' num2str(info.fval), '; ineq = ', num2str(opt_fea(1)), '; eq = ', num2str(opt_fea(2))])
    net_opt=net;
    xch_opt=xch;
    info.fval=init_score;
else
    net_opt=info.x(1:net_l);
    xch_opt=info.x(net_ll:end);
end

function score = objfunc(simulate,net,xch,input,mea,fmea,var)
%% objective function
[~,score]=varssr(simulate,net,xch,input,mea,fmea,var);