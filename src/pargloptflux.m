function [net_opt,xch_opt,info,net_,xch_,info_,fval_,I] = pargloptflux(simulate,model,net,xch,ineq,eq,input,mea,fmea,var,iter)
%% parallel global optimization of flux distribution
% Dependencies: fluxest2.m
% input: simulate (simulator function; model (MATLAB structure
% containing metabolic model); net (vector of free net fluxes); xch (vector
% of free exchange fluxes); ineq (inequality constraints); eq (equality
% constraints); input (input isotope substrate labeling patterns; mea
% (average measured labeling patterns), fmea (structure of flux 
% measurements), var (variances of measured labeling), iter (maximum number
% of iterations)
% output: net_opt and xch_opt (optimal net and exchange free fluxes); info
% (info regarding best fit), net_ and xch_ (all iterations of net and xch
% fluxes), info_ (info from all iterations), fval_ (function values from
% all iterations), I (index of minimal score)

% repeat fluxest2 at random starting points

% default number of iterations is 10
if nargin<11 || isempty(iter)
    iter=10;
end

% determine number of net and xch constraints; initialize cell vectors to
% store net and xch fluxes determined for each iteration
net_s=size(net);
xch_s=size(xch);
net_=cell(iter,1);
xch_=cell(iter,1);
info_=cell(iter,1);

parfor i=1:iter
    net=-1.5+3*rand(net_s);
    xch=3*rand(xch_s);
    [net_{i},xch_{i},info_{i}]=fluxest2(simulate,model,net,xch,ineq,eq,input,mea,fmea,var);
end

v=cell2mat(info_);
fval_=[v(:).fval];
[~,I]=min(fval_);

net_opt=net_{I};
xch_opt=xch_{I};
info=info_{I};