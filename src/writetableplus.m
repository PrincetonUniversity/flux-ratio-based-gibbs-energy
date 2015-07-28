function [tflux,tiso,score] = writetableplus(simulate,model,net_opt,xch_opt,input,mea,fmea,var,temperature)
%% writes tables summarizing results
% Dependencies: consfree2full.m, varssr.m, svec.m
% input: simulate (function handle), model (MATLAB structure containing
% metabolic model), net_opt and xch_opt (free flux vectors), input
% (input isotope substrate labeling patterns), mea and var (average 
% measured labeling patterns and variances), fmea (structure of flux 
% measurements), temperature (intracellular temperature)
% output: tflux (summary of fluxes), tiso (summary of isotope
% distributions), score (optimum (co)variance-weighted sum of squared 
% residuals)
% not implemented yet: tconf (summary of confidence intervals)

if nargin<9 || isempty(temperature)
    temperature=310;
end

[fr,nx]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net_opt,xch_opt);
[cur,score]=varssr(simulate,net_opt,xch_opt,input,mea,fmea,var);
[obs,obs_i,obs_l,meas]=svec(mea);
diff=obs-cur;
dvar=diag(var);
fr_ratio=fr(:,1)./fr(:,2);
dg=-8.314*temperature*log(fr_ratio)/1000;

metmn=cell(obs_l,1);
for i=1:length(obs_i)-1
    k=0;
    for j=obs_i(i)+1:obs_i(i+1)
        metmn{j}=strcat(meas{i},'+',num2str(k));
        k=k+1;
    end
end

tflux=table(nx(:,1),nx(:,2),fr(:,1),fr(:,2),fr_ratio,dg,'RowNames',model.rxns(:),'VariableNames',{'net','exchange','forward','reverse','for_rev_ratio','dG_kJ_mol'});
tiso=table(obs,cur,abs(diff),dvar,diff.^2./dvar,'RowNames',metmn,'VariableNames',{'measured','simulated','abs_diff','variance','var_weighted_ssr'});

% if nargin>8
%     tconf=table
% end

% % CHECK COMPATIBLE MATLAB VERSIONS
% if verLessThan('matlab', '7.0.1')
%     % Put code to run under MATLAB older than MATLAB 7.0.1 here
% else
%     % Put code to run under MATLAB 7.0.1 and newer here
% end
% V = ver('MATLAB')
% V = 
%          Name: 'MATLAB'
%       Version: '7.11'
%       Release: '(R2010b)'
%          Date: '03-Aug-2010'