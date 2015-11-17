function [c,ceq,rxndG] = therconsineq(metConcDfG,dg_mea,model,f0r1)

R = 8.3144621e-3; % Gas constant in kJ/(K*mol)
T = model.T; % Temperature in K
S = model.S;

nrxn=size(dg_mea,1);
l_mets=length(metConcDfG)/2;
metConc=metConcDfG(1:l_mets);
metConc(metConc<0)=0;

for i=nrxn:-1:1
    I_subs=(S(:,i)~=0);
    Q(i)=prod(metConc(I_subs,1).^S(I_subs,i));
end

if f0r1==0
    dgs_cc=dgf2dgr(model,metConcDfG(l_mets+1:end));
else % if f0r1==1
    dgs_cc=model.DrGt0;
end
rxndG=dgs_cc+R*T*log(Q(:));
c=[dg_mea(:,2)-rxndG; rxndG-dg_mea(:,3)];
ceq = [];