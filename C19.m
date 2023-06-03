function [val] = C19(x, o1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    y = x - o1;
    if g1(y) <= 0 && g2(y) <= 0
        val = f(y);
    end
end

function [val] = f(y)
% Original objective function f
    val = sum(abs(y).^0.5 + 2*sin(y.^3));
end

function [val] = g1(y)
% Constraint function g(x)
    d = length(y);
    v = -10*exp(-0.2 * sqrt(y.^2 + circshift(y.^2, -1)));
    val = sum(v) - v(end) + (d - 1) * 10 / exp(-5);
end

function [val] = g2(y)
% Constraint function g(x)
    d = length(y);
    val = sum(sin(2*y).^2) - 0.5*d;
end
