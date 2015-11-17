function [c,ceq,rxndG] = therconsineq_ln(metConcDfG,T,S)

R = 8.3144621e-3; % Gas constant in kJ/(K*mol)
S = full(S);

rxndG=[R*T*S' S']*metConcDfG;

c = [];
ceq = [];