function flag = isrxn(rowv,rxn)
flag=1;
if size(rowv,1)~=1
    error('reaction %s does not exist',rxn)
end