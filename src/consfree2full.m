function [gross_full,net_xch,Nfree,all_free,cons_i,rm_cons,free_i] = consfree2full(Nnet,Nxch,fluxes,net,xch)
%% map the full flux distribution from the free fluxes and constraints
% input: Nnet(net null space matrix); Nxch (exchange null space matrix);
% fluxes (flux matrix); net (vector of free net fluxes);
% xch (vector of free exchange fluxes)
% output: gross_full (full forward and reverse fluxes); net_xch (full net
% and exchange fluxes); Nfree (null space), all_free (all free fluxes), 
% cons_i (logical indices of quasi/constraint fluxes in all net and xch 
% fluxes), rm_cons (logical indices of constraint fluxes in all_free)
i_free_net=logical(fluxes(:,7));
i_free_xch=logical(fluxes(:,8));
i_cons_net=logical(fluxes(:,9));
i_cons_xch=logical(fluxes(:,10));
free_i=[i_free_net; i_free_xch];

% get full net flux distribution from input "net"
u_net=fluxes(:,1);
u_net(i_free_net)=net;
u_net=u_net(i_free_net | i_cons_net);
net_full=Nnet*u_net;

% get full xch flux distribution from input "xch"
u_xch=fluxes(:,2);
u_xch(i_free_xch)=xch;
u_xch=u_xch(i_free_xch | i_cons_xch);
xch_full=Nxch*u_xch;

net_xch=[net_full xch_full];

% get full forward and reverse flux distribution from input "net" and "xch"
gross_full=[xch_full xch_full];
for i=1:length(net_full)
    if net_full(i)<0
        gross_full(i,2)=gross_full(i,2)-net_full(i);
    else
        gross_full(i,1)=gross_full(i,1)+net_full(i);
    end
end

if nargout>2
%     null space; Nfree=[Nnet Nxch; -Nnet Nxch];
    Nfree=[Nnet zeros(size(Nxch)); zeros(size(Nnet)) Nxch];
%     all free fluxes according to stoichiometry (including constrained ones)
    all_free=[u_net; u_xch];
%     logical indices of quasi/constraint fluxes in all net and xch fluxes
    cons_i=[(i_cons_net | logical(fluxes(:,5))); (i_cons_xch | logical(fluxes(:,6)))];
%     logical indices of constraint fluxes in all_free
    rm_cons=false(size(all_free));
    u_net_i=find(i_free_net | i_cons_net);
    net_cons_i=find(i_cons_net);
    u_net_l=length(u_net);
    for i=1:length(net_cons_i)
        rm_cons(1:u_net_l)=rm_cons(1:u_net_l) | (u_net_i==net_cons_i(i));
    end
    u_xch_i=find(i_free_xch | i_cons_xch);
    xch_cons_i=find(i_cons_xch);
    for i=1:length(xch_cons_i)
        rm_cons(u_net_l+1:end)=rm_cons(u_net_l+1:end) | (u_xch_i==xch_cons_i(i));
    end
end