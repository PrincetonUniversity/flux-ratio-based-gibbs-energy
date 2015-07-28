function J = numjac(simulate,net,xch,input,mea)
%% calculate numerical jacobian matrix
persistent F net_l xch_l
if isempty(net_l) || net_l~=length(net) || xch_l~=length(xch)
    net_l=length(net);
    xch_l=length(xch);
    F=@(x)varssr(simulate,x(1:net_l),x(net_l+1:end),input,mea,[]);
end
[~,~,~,~,~,~,J] = lsqnonlin(F,[net; xch],[],[],optimset('MaxFunEval',0,'display','off'));
% J=jacobianest(F,[net; xch]); % slow/read the license info under MATLAB/DERIVESTsuite
J=full(J);