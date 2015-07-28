function [cur,score] = varssr(simulate,free_net,free_xch,input,mea,fmea,var)
%% variance weighted sum of squared residuals
% Dependencies: de2bi.m, svec.m
% input: simulate (function handle), free_net and free_xch
% (free flux vectors), input (input isotope substrate labeling patterns),
% mea and var (average measured labeling patterns and variances), 
% fmea (structure of flux measurements), 
% output: cur (simulated labeling patterns), score (score to minimize)
%% conversion vector
% isotopomers -> #C labeling fractions
% calculate nc only the first time of function call
persistent nc
if isempty(nc)
    nc=sum(de2bi(0:2^12-1),2);
end
[obs,obs_i,obs_l]=svec(mea);

%% calculate score and simulate labeling patterns
sim=simulate(free_net,free_xch,input);
meas=fieldnames(mea);
cur=zeros(obs_l,1);
for i=1:length(meas)
    simnc=zeros(2,1);
    for j=length(mea.(meas{i})):-1:1
        simnc(j)=sum(sim.(meas{i})(nc(1:length(sim.(meas{i})))==j-1));
    end
    cur(obs_i(i)+1:obs_i(i+1))=simnc;
end
if nargout>1
    diff=cur-obs;
%     score=(diff'/var)*diff;
    score=(diff'*pinv(var))*diff;
    % flux measurement
    if ~isempty(fmea)
        if size(fmea.ind,2)==1
            score=score+sum(((fmea.avg(:)-free_net(fmea.ind)).^2)./fmea.var(:));
        else
            score=score+sum(((fmea.avg(:)-fmea.ind*[free_net; free_xch]).^2)./fmea.var(:));
        end
    end
end