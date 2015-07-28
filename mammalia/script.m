%% Mammalia flux and confidence interval estimation script
% input metabolites
input.CO2_IN=isotopomervector([0],1);

input.GLC12_IN=isotopomervector([1 1 0 0 0 0],1);
input.GLC0_IN=isotopomervector([0 0 0 0 0 0], 1);
input.GLCU_IN=isotopomervector([1 1 1 1 1 1],1);

input.GlnU_IN=isotopomervector([1 1 1 1 1], 1);

input.Ser_IN=isotopomervector([0 0 0], 1);
input.Gly_IN=isotopomervector([0 0], 1);
input.Cys_IN=isotopomervector([0 0 0], 1);
input.Val_IN=isotopomervector([0 0 0 0 0], 1);
input.Leu_IN=isotopomervector([0 0 0 0 0 0], 1);
input.Phe_IN=isotopomervector([0 0 0 0 0 0 0 0 0], 1);
input.Tyr_IN=isotopomervector([0 0 0 0 0 0 0 0 0], 1);
input.Trp_IN=isotopomervector([0 0 0 0 0 0 0 0 0 0 0], 1);
input.Thr_IN=isotopomervector([0 0 0 0], 1);
input.Gln_IN=isotopomervector([0 0 0 0 0], 1);
input.Arg_IN=isotopomervector([0 0 0 0 0 0], 1);
input.His_IN=isotopomervector([0 0 0 0 0 0], 1);
input.Ile_IN=isotopomervector([0 0 0 0 0 0], 1);
input.Met_IN=isotopomervector([0 0 0 0 0], 1);
input.Lys_IN=isotopomervector([0 0 0 0 0 0], 1);

input.AC_IN=isotopomervector([0 0], 1);

% for JAVA enabled computers (Mac, PC, Linux)
[model,free_net,free_xch]=readFluxModel('mammalia_061314.h5');
[mea,covar,raw,avgvar]=xlsreadplus('mammalia_mea.xlsx',{'12_13C','12_13C_add','UGln_13C','U_13C'});
fmea = xlsreadfmea('mammalia_mea.xlsx',model,free_net,free_xch);
ineq=xlsreadineq('mammalia_mea.xlsx',model,free_net,free_xch);
save mammalia_mea_final5

simulate=@mammalia_061314;

% account for extremely low variance to avoid overfitting
minvar=0.1*avgvar;
var=diag(covar);
var(var<minvar)=minvar;
var=diag(var);

%% evaluate the best fit fluxes and labeling pattern
[net_opt,xch_opt,info]=pargloptflux(simulate,model,free_net,free_xch,ineq,[],input,mea,fmea,var,50);

[fr,nx]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net_opt,xch_opt);
[cur,score]=varssr(simulate,net_opt,xch_opt,input,mea,fmea,var);
[tflux,tiso,tscore]=writetableplus(simulate,model,net_opt,xch_opt,input,mea,fmea,var);

%% estimate confidence intervals
[lb,ub,hs,net_nopt,xch_nopt,nscore,nflag,I]=parconfest3_1(simulate,model,net_opt,xch_opt,ineq,input,mea,fmea,var);
