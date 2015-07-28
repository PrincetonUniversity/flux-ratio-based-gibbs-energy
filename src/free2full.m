function [gross_full,net_full,Sfree,Nfree,all_free,Nfree_i,cons_i] = free2full(S,fluxes,net,xch)
%% map the full flux distribution from the free fluxes
i_free_net=logical(fluxes(:,7));
i_free_xch=logical(fluxes(:,8));
i_cons_net=logical(fluxes(:,9)); % cons net
i_cons_xch=logical(fluxes(:,10)); % cons xch

i_depd_net=logical(fluxes(:,3)+fluxes(:,5)); % depd and qcon net
% i_depd_xch=logical(fluxes(:,4)+fluxes(:,6)); % depd and qcon xch

Snet=[S(:,i_depd_net) S(:,(i_free_net | i_cons_net))];
Sgross=[Snet -S];

I=linspace(1,length(i_free_net),length(i_free_net));
% index to arrange net (and forward) fluxes back in order
index=[I(i_depd_net) I(i_free_net | i_cons_net)];

% get full net flux distribution from input "net"
u_net=fluxes(:,1);
u_net(i_free_net)=net;
u_net=u_net(i_free_net | i_cons_net);
net_full=null(Snet,'r')*u_net;
net_full(index)=net_full;

% get full forward and reverse flux distribution from input "net" and "xch"
u_xch=fluxes(:,2);
u_xch(i_free_xch)=xch;
gross_full=null(Sgross,'r')*[u_net; u_xch];
gross_full=reshape(gross_full,[],2);
gross_full(index,1)=gross_full(:,1);

%% calculate both net and gross free fluxes and corresponding null space
% need this for confidence intervals
if nargout>2
%     find xch rxns that are free or (non-zero and quasi/constrained)
    i_nzqc=(logical(fluxes(:,2)) & i_cons_xch);
    xch_logical=(i_free_xch | i_nzqc);
    Sfree=[Snet -S(:,xch_logical) ];
    Nfree=null(Sfree,'r');
    u_xch(~xch_logical)=[];
    all_free=[u_net; u_xch];
    Nfree_i=[index size(fluxes,1)+I(xch_logical)];
    cons_i=[find(i_cons_net); size(fluxes,1)+find(i_nzqc)];
end