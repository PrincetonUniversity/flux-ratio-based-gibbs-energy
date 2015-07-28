function x = isotopomer2mid(sim,mea)
%% conversion vector
% isotopomers -> #C labeling fractions
% calculate nc only the first time of function call
persistent nc
if isempty(nc)
    nc=sum(de2bi(0:2^12-1),2);
end

%% add isotopomer fractions according to #C labeled
meas=fieldnames(mea);
for i=1:length(meas)
    simnc=0;
    try
        ll=length(sim.(meas{i}));
        l=log2(ll)+1;
    catch
        l=length(mea.(meas{i}));
    end
    for j=l:-1:1
        simnc(j)=sum(sim.(meas{i})(nc(1:length(sim.(meas{i})))==j-1));
    end
    x.(meas{i})=simnc;
end
