function [o] = make_translate(a, d)
%TRANSLATE Summary of this function goes here
%   Detailed explanation goes here
    v = a*ones(d,1);
    o = unifrnd(-4/5*v, 4/5*v);
end
