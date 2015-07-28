function [ v ] = cumulate(v,dir)
% Fast-Isotopomer/Cumomer-Transformation
%
% Konvertiert zwischen Isotopomer- (=Vektor mit Summe 1) und kumulierter
% Isotopomer(=Cumomer)-Darstellung. Der zweite Parameter ist optional und
% gibt die Richtung der Transformation an.
%
% Isotopomer -> Cumomer:
%    c = cumulate(v)
%
% Cumomer -> Isotopomer:
%    v = cumulate(c,-1);
%
% Copyright: Michael Weitzel <mich@el-weitzel.de>
if nargin==1 || dir>0
    dir = +1;
else
    dir = -1;
end

% entspricht einer Radix-2-DIF-FFT ...
n = length(v);
for ldm=log2(n):-1:1
    m = 2^ldm;
    mh = m/2;
    for r=0:m:n-1
        for j=0:mh-1
            s = n-r-j;
            t = s-mh;
            % wird hier das +/- durch ein XOR ersetzt, entspricht das einer
            % Reed-Muller Transformation, deren Transformationsmatrix per
            % flipud gespiegelt wurde.
            v(t) = v(t) + dir * v(s);
        end
    end
end
