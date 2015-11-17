function [m,info_] = interconfint2(m,rxns,iter)
%% 
% define constants
if nargin<3 || isempty(iter)
    iter=1;
end
l_mets=size(m.metConc,1);
R=8.3144621e-3; % Gas constant in kJ/(K*mol)
T=m.T;

%% combine divided E. coli ACONT and yeast AKGDH reactions.
l_rxns_orig=length(m.rxns);

a=strcmp(m.rxns,'ACONTa');
b=strcmp(m.rxns,'ACONTb');
if sum(a) && sum(b)
    m.S=[m.S m.S(:,a)+m.S(:,b)];
    m.rxns=[m.rxns; 'ACONT'];
    m.rev=[m.rev; min(m.rev(a),m.rev(b))];
    m.lb=[m.lb; max(m.lb(a),m.lb(b))];
    m.ub=[m.ub; min(m.ub(a),m.ub(b))];
    m.c=[m.c; min(m.c(a),m.c(b))];
    m.rules=[m.rules; ''];
    m.rxnGeneMat=[m.rxnGeneMat; m.rxnGeneMat(a,:)+m.rxnGeneMat(b,:)];
    m.grRules=[m.grRules; ''];
    m.subSystems=[m.subSystems; ''];
    m.confidenceScores=[m.confidenceScores; ''];
    m.rxnReferences=[m.rxnReferences; ''];
    m.rxnECNumbers=[m.rxnECNumbers; ''];
    m.rxnNotes=[m.rxnNotes; ''];
    m.rxnNames=[m.rxnNames; ''];
    m.DrG0=[m.DrG0; m.DrG0(a)+m.DrG0(b)];
    m.ur=[m.ur; min(m.ur(a),m.ur(b))];
    m.DrGt0=[m.DrGt0; m.DrGt0(a)+m.DrGt0(b)];
    m.DrGtMin=[m.DrGtMin; m.DrGtMin(a)+m.DrGtMin(b)];
    m.DrGtMax=[m.DrGtMin; m.DrGtMax(a)+m.DrGtMax(b)];
    m.quantDir=[m.quantDir; min(m.quantDir(a),m.quantDir(b))];
    m.rxndG_mea=[m.rxndG_mea; -2 -6 0];
    m.rxndG_meaInd=[m.rxndG_meaInd; m.rxndG_meaInd(a)];
    m.rxndG=[m.rxndG; m.rxndG(a,:)];
end

c=strcmp(m.rxns,'AKGDam');
d=strcmp(m.rxns,'AKGDbm');
e=strcmp(m.rxns,'PDHcm');
if sum(c) && sum(d) && sum(e)
    m.S=[m.S m.S(:,c)+m.S(:,d)+m.S(:,e)];
    m.rxns=[m.rxns; 'AKGDm'];
    m.rev=[m.rev; min(m.rev(c),m.rev(d))];
    m.lb=[m.lb; max(m.lb(c),m.lb(d))];
    m.ub=[m.ub; min(m.ub(c),m.ub(d))];
    m.c=[m.c; min(m.c(c),m.c(d))];
    m.rules=[m.rules; ''];
    m.rxnGeneMat=[m.rxnGeneMat; m.rxnGeneMat(c,:)+m.rxnGeneMat(d,:)+m.rxnGeneMat(e,:)];
    m.grRules=[m.grRules; ''];
    m.subSystems=[m.subSystems; ''];
    m.confidenceScores=[m.confidenceScores; ''];
    m.rxnReferences=[m.rxnReferences; ''];
    m.rxnECNumbers=[m.rxnECNumbers; ''];
    m.rxnNotes=[m.rxnNotes; ''];
    m.rxnNames=[m.rxnNames; ''];
    m.DrG0=[m.DrG0; m.DrG0(c)+m.DrG0(d)+m.DrG0(e)];
    m.ur=[m.ur; min(m.ur(c),m.ur(d))];
    m.DrGt0=[m.DrGt0; m.DrGt0(c)+m.DrGt0(d)+m.DrGt0(e)];
    m.DrGtMin=[m.DrGtMin; m.DrGtMin(c)+m.DrGtMin(d)+m.DrGtMin(e)];
    m.DrGtMax=[m.DrGtMin; m.DrGtMax(c)+m.DrGtMax(d)+m.DrGtMax(e)];
    m.quantDir=[m.quantDir; min(m.quantDir(c),m.quantDir(d))];
    m.rxndG_mea=[m.rxndG_mea; m.rxndG_mea(c,:)];
    m.rxndG_meaInd=[m.rxndG_meaInd; m.rxndG_meaInd(c)];
    m.rxndG=[m.rxndG; m.rxndG(c,:)];
end

%% find DrG that are out of bounds and which DfG needs to change
I_DrG=m.rxndG_meaInd;
if length(I_DrG)>l_rxns_orig
    I_DrG(end)=1;
end
if ~isempty(rxns)
    for i=1:length(rxns)
        I_DrG=I_DrG | strcmp(m.rxns,rxns(i));
    end
end
I_DfG=logical(full(abs(m.S))*double(I_DrG));

%% assign constraints and bounds
I_Conc=I_DfG;
F_Conc=find(I_Conc);
l_Conc=length(F_Conc);
metConc=m.metConc;
metConc(metConc<1e-10)=1e-10;
m.logMetConc=log(metConc);
lb1=m.logMetConc(I_Conc,2);
ub1=m.logMetConc(I_Conc,3);

m.DfGt0_mea=m.DfGt0;
lb2=m.DfGt0_mea(I_DfG)-2*m.uf(I_DfG);
ub2=m.DfGt0_mea(I_DfG)+2*m.uf(I_DfG);

lb=[lb1; lb2];
ub=[ub1; ub2];

I_ConcDfG=[I_Conc; I_DfG];

% save indices to the model
m.I_Conc=I_Conc;
m.I_DfG=I_DfG;
m.I_DrG=I_DrG;

% initial point
x_I_=mean(m.logMetConc(I_Conc,2:3),2);
x_I=[x_I_; m.DfGt0(I_DfG,1)];

%% optimize [met] and DfG that best match expected values
% all [met] and DrG bounds must be satified
hBool = strcmp(m.metFormulas,'H');
S=m.S;
S(hBool,:) = 0; % Set proton coefficients to 0
S=S(I_Conc,:);
S=S(:,I_DrG);
logMetConc_mea=log(m.metConc_mea(I_Conc,:));
rxndG_mea=m.rxndG_mea(I_DrG,:);
metDfG_mea=m.DfGt0_mea(I_DfG,:);
uf=m.uf(I_DfG,:);
S=full(S);
F=@(x)objfunc(x,logMetConc_mea,metDfG_mea,uf,rxndG_mea,T,S);
Fnlc=[];

maxiter=500;
options=optimset('Algorithm','interior-point' ...
                ,'Display','iter-detailed' ...
                ,'TolX',1e-10,'TolCon',1e-10,'TolFun',1e-10 ...
                ,'MaxFunEvals',300*maxiter,'MaxIter',maxiter);

info_=cell(iter,1);
for i=1:iter
    x_l=length(x_I);
    init_score=F(x_I);
    disp([num2str(i) ': number of variables = ' num2str(x_l)])
    disp([num2str(i) ': optimization begins with score = ' num2str(init_score)])
    try
        [info_{i}.x,info_{i}.fval,info_{i}.exitflag,info_{i}.output]=fmincon(F,x_I,[-R*T*S' -S'; R*T*S' S'],[-rxndG_mea(:,2); rxndG_mea(:,3)],[],[],lb,ub,Fnlc,options);
    catch
        warning('error in fmincon')
%         metConcDfG=lognrnd(-8.5,0.5,l_mets,1);
        [info_{i}.x,info_{i}.fval,info_{i}.exitflag,info_{i}.output]=fmincon(F,x_I,[-R*T*S' -S'; R*T*S' S'],[-rxndG_mea(:,2); rxndG_mea(:,3)],[],[],lb,ub,Fnlc,options);
    end
end
v=cell2mat(info_);
fval_=[v(:).fval];
[~,I]=min(fval_);
info=info_{I};

metConcDfG=[m.logMetConc(:,1); m.DfGt0(:,1)];
metConcDfG(I_ConcDfG)=info.x;
m.logMetConc(:,1)=metConcDfG(1:l_mets);
m.metConc(:,1)=exp(m.logMetConc(:,1));
m.DfGt0=metConcDfG(l_mets+1:end);
[rxndG2,m.I_inCI,m.N_inCI,m.N_mr,m.ssr,m.ssr2]=thercons([m.metConc(:,1); m.DfGt0(:,1)],m.metConc_mea,m.rxndG_mea,m);
m.rxndG2=[rxndG2 m.rxndG(:,2:3)];

%% intersection of overlapping confidence intervals 
% refined lower and upper bounds of
% metabolite concentrations and reaction free energies
% that satisfy all measurement confidence intervals

m.DrGt0_mea=m.DrGt0;
m.DrGt0=dgf2dgr(m,m.DfGt0);
lb_m=zeros(l_Conc,1);
ub_m=lb_m;

% metabolite concentrations
A=[-R*T*S'; R*T*S'];
b=[-rxndG_mea(:,2)+m.DrGt0(I_DrG); rxndG_mea(:,3)-m.DrGt0(I_DrG)];
options = optimset('Display','none');
for i=1:l_Conc
    f=zeros(l_Conc,1);
    f(i)=1;
    x=linprog(f,A,b,[],[],lb1,ub1,[],options);
    lb_m(i)=x(i);
    f=-f;
    x=linprog(f,A,b,[],[],lb1,ub1,[],options);
    ub_m(i)=x(i);
    if mod(i,10)==0
        disp([num2str(i) '/' num2str(l_Conc)])
    end
end

F_DrG=find(I_DrG);
l_DrG=length(F_DrG);
lb_g=zeros(l_DrG,1);
ub_g=lb_g;

% gibbs free energies
for i=1:l_DrG
    f=S(:,i);
    x=linprog(f,A,b,[],[],lb1,ub1,[],options);
    lb_g(i)=R*T*f'*x;
    f=-f;
    x=linprog(f,A,b,[],[],lb1,ub1,[],options);
    ub_g(i)=R*T*-f'*x;
    if mod(i,5)==0
        disp([num2str(i) '/' num2str(l_DrG)])
    end
end

% m.metConc(I_Conc,1)=exp((lb_m+ub_m)/2);
m.metConc(I_Conc,2)=exp(lb_m);
m.metConc(I_Conc,3)=exp(ub_m);

[rxndG2,m.I_inCI,m.N_inCI,m.N_mr,m.ssr,m.ssr2]=thercons([m.metConc(:,1); m.DfGt0(:,1)],m.metConc_mea,m.rxndG_mea,m);
m.rxndG2=[rxndG2 m.rxndG(:,2:3)];
m.rxndG2(I_DrG,2)=lb_g+m.DrGt0(I_DrG);
m.rxndG2(I_DrG,3)=ub_g+m.DrGt0(I_DrG);

disp('done!')

function ssr = objfunc(x,conc_mea,metDfG_mea,uf,dg_mea,T,S)
%% objective function for quadratic programming
[~,~,~,~,ssr]=thercons_ln(x,conc_mea,metDfG_mea,uf,dg_mea,T,S);
