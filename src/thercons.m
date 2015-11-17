function [rxndG,I_inci,ninci,nmetrxn,ssr,ssr2] = thercons(metConcDfG,conc_mea,dg_mea,model,f0r1)

if nargin<5 || isempty(f0r1)
    f0r1=0;
end

nmet=size(conc_mea,1);
nrxn=size(dg_mea,1);
nmetrxn=nmet+nrxn;

[~,~,rxndG]=therconsineq(metConcDfG,dg_mea,model,f0r1);

l_mets=length(metConcDfG)/2;
metConc=metConcDfG(1:l_mets);
metDfG=metConcDfG(l_mets+1:end);

I_DrG_inci=(dg_mea(:,2)<=rxndG & rxndG<=dg_mea(:,3));
I_inci=[(conc_mea(:,2)<=metConc & metConc<=conc_mea(:,3)); I_DrG_inci];
ninci=sum(I_inci);

i_conc_mea=~(isnan(conc_mea(:,1)) | conc_mea(:,1)>0.999);
i_dg_mea=~isnan(dg_mea(:,1));

if nargout==5
    ssr=sum(((metConc(i_conc_mea)-conc_mea(i_conc_mea,1))./(conc_mea(i_conc_mea,3)-conc_mea(i_conc_mea,2))*4).^2)/sum(i_conc_mea) + ...
        sum(((metDfG-model.DfGt0)./(model.uf)).^2)/l_mets + ...
        sum(((rxndG(i_dg_mea)-dg_mea(i_dg_mea,1))./(dg_mea(i_dg_mea,3)-dg_mea(i_dg_mea,2))*4).^2)/sum(i_dg_mea);
elseif nargout==6
    ssr=0;
    ssr2=sum(((metConc(i_conc_mea)-conc_mea(i_conc_mea,1))./(conc_mea(i_conc_mea,3)-conc_mea(i_conc_mea,2))*4).^2)/sum(i_conc_mea) + ...
         sum(((rxndG(i_dg_mea)-dg_mea(i_dg_mea,1))./(dg_mea(i_dg_mea,3)-dg_mea(i_dg_mea,2))*4).^2)/sum(i_dg_mea);
end
