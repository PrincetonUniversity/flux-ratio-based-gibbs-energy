%% integrate metabolite concentration and reaction free energy measurements
% input: Metabolic model with Flux-Ratio-based ENergy and concentrations

%% ecoli
load ../ecoli/ecoli_integrate
rxn_list={'PFK','GAPD','PGK','PGM','ENO','PYK','PDH','G6PDH2r','PGL','GND','RPI','RPE','PPC','ME1','ME2','CS','ACONTa','ACONTb','ICDHyr','AKGDH','SUCOAS','FUM','SUCDi','FRD2','FRD3','EDD','EDA'};
[ecoli_FREN2,ecoli_info2]=interconfint2(ecoli_FREN,rxn_list,1);

%% yeast
load ../scerevisiae/yeast_integrate
rxn_list={'PFK','GAPD','PGK','PGM','ENO','PYK','PDHm','G6PDH2','PGL','GND','RPI','RPE','PC','ME1m','ME2m','CSm','ACONTm','ICDHxm','SUCOASm','SUCD1m','FUMm','MDHm','AKGDam','AKGDbm','PDHcm'};
[yeast_FREN2,yeast_info2]=interconfint2(yeast_FREN,rxn_list);

%% mammalian cells
load ../mammalia/mammalia_integrate
rxn_list={'PFK','GAPD','PGK','PGM','ENO','PYK','PDHm','G6PDH2r','PGL','GND','RPI','RPE','TALA','PCm','ME1m','ME2','ME2m','CSm','ACONTm','ICDHxm','AKGDm','SUCOASm','SUCOAS1m','SUCD1m','FUMm','MDHm','PGI','FBA','TPI','TKT1','TKT2','TALA','GHMT2r','GHMT2rm','DPGM','DPGase','MTHFD','MTHFDm','MTHFD2','MTHFD2m'};
[mammalia_FREN2,mammalia_info2]=interconfint2(mammalia_FREN,rxn_list);
