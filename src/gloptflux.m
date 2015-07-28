function [net_opt,xch_opt,info] = gloptflux(simulate,model,net,xch,ineq,eq,input,mea,fmea,var,iter)
%% global optimization of flux distribution
% repeat fluxest2 at random starting points
% default number of iterations is 10

if nargin<11 || isempty(iter)
    iter=10;
end
[net_opt,xch_opt,info]=fluxest2(simulate,model,net,xch,ineq,eq,input,mea,fmea,var);
for i=2:iter
%     try
    net=-1.5+3*rand(size(net));
    xch=3*rand(size(xch));
    [net_,xch_,info_]=fluxest2(simulate,model,net,xch,ineq,eq,input,mea,fmea,var);
%     catch
%         warning([num2str(i) ': error in for fluxest2'])
%         continue
%     end
    if info_.fval<info.fval
        net_opt=net_;
        xch_opt=xch_;
        info=info_;
    end
end