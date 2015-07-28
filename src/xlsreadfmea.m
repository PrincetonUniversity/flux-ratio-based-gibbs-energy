function [fmea,textrep] = xlsreadfmea(xlsname,model,net,xch)
%% read excel file containing measured fluxes
% Dependencies: consfree2full.m
% input: xlsname (excel file name) containing reaction names and their
% measured averages and variances; left most side must contain only
% reaction names and the other side, only numbers; model (MATLAB structure
% containing metabolic model); net (vector of free net fluxes); xch (vector
% of free exchange fluxes)
% output: fmea (structure of flux measurements); textrep
% (flux measurement name(s))

[~,~,N,all_free,~,rm_cons]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net,xch);
rxns=model.rxns;

if ~isempty(xlsname)
    sheets='fmea';
    if ispc
        [~,~,alldata]=xlsread(xlsname,sheets);
    elseif isunix
        warning('off','MATLAB:xlsread:Mode')
        [~,~,alldata]=xlsread(xlsname,sheets,'','basic');
    end
    alldata(~cellfun(@ischar,alldata(:,1)),:)=[];
    alldata_l=size(alldata,1);
    fmea.ind=zeros(alldata_l,size(N,2));
    fmea.avg=zeros(alldata_l,1);
    fmea.var=zeros(alldata_l,1);
    textrep=cell(alldata_l,1);
    for i=1:alldata_l
        % remove spaces and punctuation
        alldata{i,1}(ismember(alldata{i,1},' ,.:;!')) = [];
        textrep{i}=alldata{i,1};
        psubadd=find(alldata{i,1}=='-' | alldata{i,1}=='+');
        nsubadd=length(psubadd);
        if nsubadd==0
            fmea.ind(i,:)= N(strcmp(alldata{i,1},rxns),:);
            fmea.avg(i)= alldata{i,2};
            fmea.var(i)= alldata{i,3};
        else % when there is 1 or more -/+ operators
            if psubadd(1)==1 && alldata{i,1}(1)=='-'
                fmea.ind(i,:)= -N(strcmp(alldata{i,1}(2:psubadd(2)-1),rxns),:);
                jj=2;
            else
                fmea.ind(i,:)= N(strcmp(alldata{i,1}(1:psubadd(1)-1),rxns),:);
                jj=1;
            end
            fmea.avg(i)= alldata{i,2};
            fmea.var(i)= alldata{i,3};
            for j=jj:nsubadd-1
                if alldata{i,1}(j)=='-'
                    fmea.ind(i,:)=fmea.ind(i,:)-N(strcmp(alldata{i,1}(psubadd(j)+1:psubadd(j+1)-1),rxns),:);
                elseif alldata{i,1}(j)=='+'
                    fmea.ind(i,:)=fmea.ind(i,:)+N(strcmp(alldata{i,1}(psubadd(j)+1:psubadd(j+1)-1),rxns),:);
                end
            end
            if alldata{i,1}(psubadd(nsubadd))=='-'
                fmea.ind(i,:)=fmea.ind(i,:)-N(strcmp(alldata{i,1}(psubadd(nsubadd)+1:end),rxns),:);
            elseif alldata{i,1}(j)=='+'
                fmea.ind(i,:)=fmea.ind(i,:)+N(strcmp(alldata{i,1}(psubadd(nsubadd)+1:end),rxns),:);
            end
        end
    end
end

%% matrix A must have length(net)+length(xch) number of columns
% move constrained portions to the other side
fmea.avg=fmea.avg-fmea.ind(:,rm_cons)*all_free(rm_cons);
fmea.ind=fmea.ind(:,~rm_cons);