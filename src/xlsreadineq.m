function [ineq,textrep] = xlsreadineq(xlsname,model,net,xch)
%% read excel file containing inequality constraints
% Dependencies: consfree2full.m
% input: xlsname (excel file name) containing inequalities (<=); one side
% must contain only reaction names and the other side, only numbers; model
% (MATLAB structure containing metabolic model); net (vector of free net
% fluxes); xch (vector of free exchange fluxes)
% output: ineq (structure containing A and b (such that A*free_fluxes<=b));
% textrep (text representation; currently not used)

[~,~,N,all_free,~,rm_cons]=consfree2full(model.kernel_net,model.kernel_xch,model.fluxes,net,xch);
rxns=model.rxns;
rxns_l=length(rxns);

% constraints: all xch >= 0
ineq.A=-N(rxns_l+1:end,:);
ineq.b=zeros(size(rxns(:)));
textrep=rxns(:);

% constraints: all influx efflux net >= 0
inefflux=[];
for i=1:rxns_l
    if model.S(:,i)<=0
        inefflux=[inefflux; i];
    elseif model.S(:,i)>=0
        inefflux=[inefflux; i];
    end
end
ineq.A=[ineq.A; -N(inefflux,:)];
ineq.b=[ineq.b; zeros(size(inefflux))];
textrep=[textrep; rxns(inefflux)'];

if ~isempty(xlsname)
    sheets={'net' 'xch'};
    for s=1:2
        if ispc
            [~,~,alldata]=xlsread(xlsname,sheets{s});
        elseif isunix
            warning('off','MATLAB:xlsread:Mode')
            [~,~,alldata]=xlsread(xlsname,sheets{s},'','basic');
        end
        c1=s-1;
        alldata(~(cellfun(@ischar,alldata(:,1))|cellfun(@ischar,alldata(:,3))),:)=[];
        alldata_l=size(alldata,1);
        A2.(sheets{s})=zeros(alldata_l,size(N,2));
        b2.(sheets{s})=zeros(alldata_l,1);
        textrep2.(sheets{s})=cell(alldata_l,1);
        for i=1:alldata_l
            if ischar(alldata{i,1})
                % remove spaces and punctuation
                alldata{i,1}(ismember(alldata{i,1},' ,.:;!')) = [];
                textrep2.(sheets{s}){i}=alldata{i,1};
                psub=find(alldata{i,1}=='-');
                nsub=length(psub);
                if nsub==0
                    A2.(sheets{s})(i,:)= N([false(1,c1*rxns_l) strcmp(alldata{i,1},rxns)],:);
                    b2.(sheets{s})(i)= alldata{i,3};
                else % when there is 1 or more - operators
                    A2.(sheets{s})(i,:)= N([false(1,c1*rxns_l) strcmp(alldata{i,1}(1:psub(1)-1),rxns)],:);
                    b2.(sheets{s})(i)= alldata{i,3};
                    for j=1:nsub-1
                        A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub(j)+1:psub(j+1)-1),rxns)],:);
                    end
                    A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)-N([false(1,c1*rxns_l) strcmp(alldata{i,1}(psub(nsub)+1:end),rxns)],:);
                end
            elseif ischar(alldata{i,3})
                % remove spaces and punctuation
                alldata{i,3}(ismember(alldata{i,3},' ,.:;!')) = [];
                textrep2.(sheets{s}){i}=alldata{i,3};
                psub=find(alldata{i,3}=='-');
                nsub=length(psub);
                if nsub==0
                    A2.(sheets{s})(i,:)=-N([false(1,c1*rxns_l) strcmp(alldata{i,3},rxns)],:);
                    b2.(sheets{s})(i)=-alldata{i,1};
                else % when there is 1 or more - operators
                    A2.(sheets{s})(i,:)=-N([false(1,c1*rxns_l) strcmp(alldata{i,3}(1:psub(1)-1),rxns)],:);
                    b2.(sheets{s})(i)=-alldata{i,1};
                    for j=1:nsub-1
                        A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub(j)+1:psub(j+1)-1),rxns)],:);
                    end
                    A2.(sheets{s})(i,:)=A2.(sheets{s})(i,:)+N([false(1,c1*rxns_l) strcmp(alldata{i,3}(psub(nsub)+1:end),rxns)],:);
                end
            else    
                warning([sheets{s} ' ' num2str(i) ': inequality between two rxns or nums'])
            end
        end
    end
    ineq.A=[ineq.A; A2.(sheets{1}); A2.(sheets{2})];
    ineq.b=[ineq.b; b2.(sheets{1}); b2.(sheets{2})];
    textrep=[textrep; textrep2.net; textrep2.xch];
end

%% matrix A must have length(net)+length(xch) number of columns
% move constrained portions to the other side
ineq.b=ineq.b-ineq.A(:,rm_cons)*all_free(rm_cons);
ineq.A=ineq.A(:,~rm_cons);