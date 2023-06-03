function [val] = C02(x, o1, m1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = x - o1;
    y = m1 * z;
    if g(y) <= 0
        val = f(z);
    end
end

function [val] = f(z)
% Original objective function f
    val = sum(cumsum(z).^2);
end

function [val] = g(y)
% Constraint function g(x)
    val = sum(y.^2 - 5000*cos(0.1 * pi * y) - 4000);
end
