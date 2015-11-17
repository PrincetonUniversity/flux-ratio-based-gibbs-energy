function [rxndG,I_inci,ninci,nmetrxn,ssr] = thercons_ln(metConcDfG,conc_mea,dfg_mea,uf,dg_mea,T,S)

nmet=size(conc_mea,1);
nrxn=size(dg_mea,1);
nmetrxn=nmet+nrxn;

[~,~,rxndG]=therconsineq_ln(metConcDfG,T,S);

l_mets=length(metConcDfG)/2;
metConc=metConcDfG(1:l_mets);
metDfG=metConcDfG(l_mets+1:end);

I_DrG_inci=(dg_mea(:,2)<=rxndG & rxndG<=dg_mea(:,3));
I_inci=[(conc_mea(:,2)<=metConc & metConc<=conc_mea(:,3)); I_DrG_inci];
ninci=sum(I_inci);

i_conc_mea=~(isnan(conc_mea(:,1)) | conc_mea(:,1)>-0.001);
i_dg_mea=~isnan(dg_mea(:,1));

ssr=sum(((metConc(i_conc_mea)-conc_mea(i_conc_mea,1))./(conc_mea(i_conc_mea,3)-conc_mea(i_conc_mea,2))*4).^2)/sum(i_conc_mea) + ...
    sum(((metDfG-dfg_mea)./uf).^2)/l_mets + ...
    sum(((rxndG(i_dg_mea)-dg_mea(i_dg_mea,1))./(dg_mea(i_dg_mea,3)-dg_mea(i_dg_mea,2))*4).^2)/sum(i_dg_mea);
