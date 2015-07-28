function [obs,obs_i,obs_l,meas] = svec(mea)
%% create a single vector of all observed labeling patterns
% Input: mea (structure of measurements)
% 
% Output: obs (vector of all observed labeling fractions),obs_i and obs_l
% (indices vector and length required to decipher obs)
meas=fieldnames(mea);
meas_l=length(meas);
obs_l=0;
obs_i=zeros(meas_l+1,1);
for i=1:meas_l
    len=length(mea.(meas{i}));
    obs_l=obs_l+len;
    obs_i(i+1)=obs_l;
end

obs=zeros(obs_l,1);
for i=1:meas_l
    obs(obs_i(i)+1:obs_i(i+1))=mea.(meas{i})(:);
end