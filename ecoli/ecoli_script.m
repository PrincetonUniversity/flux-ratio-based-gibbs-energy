%% E. coli flux and confidence interval estimation script

% define input metabolites
input.CO20=isotopomervector(0,1);
input.CO20_3=isotopomervector(0,1);
input.CO20_5050=isotopomervector(0,1);

input.GLC12=isotopomervector([1 1 0 0 0 0],1);
input.GLC3=isotopomervector([0 0 1 0 0 0],1);
input.GLC5050=isotopomervector([0 0 0 0 0 0;1 1 1 1 1 1],[0.5 0.5]);

% for JAVA enabled computers (Mac, PC, Linux)
[model,free_net,free_xch]=readFluxModel('ecoli_031514.h5');
[mea,covar,raw,avgvar]=xlsreadplus('ecoli_mea.xlsx',{'12_13C', '3_13C', '5050_13C'});
fmea = xlsreadfmea('ecoli_mea.xlsx',model,free_net,free_xch);
ineq=xlsreadineq('ecoli_mea.xlsx',model,free_net,free_xch);

simulate=@ecoli_031514;

% account for extremely low variance to avoid overfitting
minvar=0.1*avgvar;
var=diag(covar);
var(var<minvar)=minvar;
var=diag(var);

%% estimate fluxes
[net_opt,xch_opt,info,net_,xch_,info_,fval_,I_]=pargloptflux(simulate,model,free_net,free_xch,ineq,[],input,mea,fmea,var,50);

%% evaluate the best fit fluxes and labeling pattern
[fr,nx]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net_opt,xch_opt);
[cur,score]=varssr(simulate,net_opt,xch_opt,input,mea,fmea,var);
[tflux,tiso,tscore]=writetableplus(simulate,model,net_opt,xch_opt,input,mea,fmea,var);

%% estimate confidence intervals
[lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I]=parconfest3_1(simulate,model,net_opt,xch_opt,ineq,input,mea,fmea,var);