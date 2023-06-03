function [Q] = make_rotate(d)
%TRANSLATE Summary of this function goes here
%   Detailed explanation goes here
    m1 = normrnd(0, 1, d);
    [Q, ~] = qr(m1);
end
