function x = isotopomervector(labels,frac)
%% Represents labeled input metabolite in isotopomer vector
% Dependencies: de2bi.m (The MathWorks, Inc.)
% Input:
% labels- atom activity vector of input substrate;
% frac- vector specifying mixture fraction; 1 if not a mixture
% Output:
% x- isotopomer vector
%% Check input
frac_l=length(frac);
if frac_l~=size(labels,1)
    labels=labels';
end
if frac_l~=size(labels,1)
    error('The number of labels must equal the number of fractions')
end
if sum(frac)~=1
    error('The sum of fractions must be 1')
end

%% conversion vector
% input labeling fractions -> isotopomer vector
iv_l=2^size(labels,2);
iv=de2bi(0:iv_l-1);

%% create isotopomer vector
x=zeros(1,iv_l);

for i=1:frac_l
    for j=1:iv_l
        if isequal(iv(j,:),labels(i,:))
            x(j)=frac(i);
        end
    end
end