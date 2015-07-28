function [x,free_net,free_xch] = readFluxModel(datain)
%% Read hierarchical data format (hdf5) and output as a matlab struct (x)
% Dependencies: hdf5read.m
% input:
% datain- .h5 file containing stoichiometry and cumomer equation system
% output:
% x- MATLAB structure containing metabolic model
% free_net- free net fluxes
% free_xch- free exchange fluxes

%%
% stoichiometric network data
x.S = hdf5read(datain,'/stoichiometry/matrix');
r=hdf5read(datain,'/stoichiometry/r');
c=hdf5read(datain,'/stoichiometry/c');
x.fluxes=hdf5read(datain,'/stoichiometry/fluxes');
r_fluxes=hdf5read(datain,'/stoichiometry/r_fluxes');
c_fluxes=hdf5read(datain,'/stoichiometry/c_fluxes');
x.kernel_net=hdf5read(datain,'/stoichiometry/kernel/matrix_net');
x.kernel_xch=hdf5read(datain,'/stoichiometry/kernel/matrix_xch');
r_kernel_net=hdf5read(datain,'/stoichiometry/kernel/r_net');
r_kernel_xch=hdf5read(datain,'/stoichiometry/kernel/r_xch');
c_kernel_net=hdf5read(datain,'/stoichiometry/kernel/c_net');
c_kernel_xch=hdf5read(datain,'/stoichiometry/kernel/c_xch');

% additional statistical data
x.cov_free=hdf5read(datain,'/cov_free/matrix');
r_cov_free=hdf5read(datain,'/cov_free/r');
c_cov_free=hdf5read(datain,'/cov_free/c');
x.cov_full=hdf5read(datain,'/cov/matrix');
r_cov=hdf5read(datain,'/cov/r');
c_cov=hdf5read(datain,'/cov/c');

x.jacobian=hdf5read(datain,'/jacobian/matrix');
r_jac=hdf5read(datain,'/jacobian/r');
c_jac=hdf5read(datain,'/jacobian/c');

%%
% store data in cells
for i=length(r):-1:1
    x.mets{i}=r(i).Data;
end
for i=length(c):-1:1
    x.rxns{i}=c(i).Data;
end
for i=length(r_fluxes):-1:1
    x.rxns_fluxes{i}=r_fluxes(i).Data;
end
for i=length(c_fluxes):-1:1
    x.types_fluxes{i}=c_fluxes(i).Data;
end
for i=length(r_kernel_net):-1:1
    x.rxns_kernel_net{i}=r_kernel_net(i).Data;
end
for i=length(r_kernel_xch):-1:1
    x.rxns_kernel_xch{i}=r_kernel_xch(i).Data;
end
for i=length(c_kernel_net):-1:1
    x.types_kernel_net{i}=c_kernel_net(i).Data;
end
for i=length(c_kernel_xch):-1:1
    x.types_kernel_xch{i}=c_kernel_xch(i).Data;
end
for i=length(r_cov_free):-1:1
    x.rxnsr_cov_free{i}=r_cov_free(i).Data;
end
for i=length(c_cov_free):-1:1
    x.rxnsc_cov_free{i}=c_cov_free(i).Data;
end
for i=length(r_cov):-1:1
    x.rxnsr_cov{i}=r_cov(i).Data;
end
for i=length(c_cov):-1:1
    x.rxnsc_cov{i}=c_cov(i).Data;
end
for i=length(r_jac):-1:1
    x.lab_jac{i}=r_jac(i).Data;
end
for i=length(c_jac):-1:1
    x.rxns_jac{i}=c_jac(i).Data;
end

%%
% sort reactions alphabetically
[x.rxns,I]=sort(x.rxns);
x.S=x.S(:,I);
x.fluxes=x.fluxes(I,:);
x.rxns_fluxes=x.rxns_fluxes(I);
x.kernel_net=x.kernel_net(I,:);
x.kernel_xch=x.kernel_xch(I,:);
x.rxns_kernel_net=x.rxns_kernel_net(I);
x.rxns_kernel_xch=x.rxns_kernel_xch(I);
x.types_kernel_net(1)=[];
x.types_kernel_xch(1)=[];
[x.types_kernel_net,I]=sort(x.types_kernel_net);
x.kernel_net=x.kernel_net(:,1+I);
[x.types_kernel_xch,I]=sort(x.types_kernel_xch);
x.kernel_xch=x.kernel_xch(:,1+I);

%%
% get free net and exchange reactions
free_net=x.fluxes(logical(x.fluxes(:,7)),1);
free_xch=x.fluxes(logical(x.fluxes(:,8)),2);
