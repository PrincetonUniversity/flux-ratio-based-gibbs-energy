function [org_avg,covar,org_raw,avgvar] = xlsreadplus(xlsname,sheets)
%% read excel file containing labeling data and create a struct
% input: xlsname (excel file name) and sheets (sheet names)
% output: org_avg (average labeling patterns), covar (covariance matrix),
% org_raw (raw labeling patterns), and avgvar (average variance)
% % % may need to add in error checking mechanism (e.g. correct # carbons)

avgvar=0;
sheets_l=length(sheets);
nmea=zeros(sheets_l+1,1);
covmat=cell(sheets_l,1);
for s=1:sheets_l
    if ispc
        [ndata,mets]=xlsread(xlsname,sheets{s});
    elseif isunix
        warning('off','MATLAB:xlsread:Mode')
        [ndata,mets]=xlsread(xlsname,sheets{s},'','basic');
%         ndata(:,1)=[];
%     elseif ismac
%         
    end
    ndata=ndata';
    I=find(~cellfun(@isempty,mets));
    I(end+1)=1+size(ndata,2);
    
    covmat{s}=cov(ndata);
    avgvar=avgvar+sum(diag(covmat{s}));
    nmea(s+1)=nmea(s)+size(ndata,2);
    
    for i=1:length(I)-1
        org_raw.(mets{I(i)})=ndata(:,I(i):I(i+1)-1);
        org_avg.(mets{I(i)})=mean(org_raw.(mets{I(i)}),1);
    end
end
for s=1:sheets_l
    covar(nmea(s)+1:nmea(s+1),nmea(s)+1:nmea(s+1))=covmat{s};
end
% average variance if not using covariance matrix
avgvar=avgvar/nmea(end);