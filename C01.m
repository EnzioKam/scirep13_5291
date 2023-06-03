function [val] = C01(x, o1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = x - o1;
    if g(z) <= 0
        val = f(x);
    end
end

function [val] = f(z)
% Original objective function f
    val = sum(cumsum(z).^2);
end

function [val] = g(z)
% Constraint function g(x)
    val = sum(z.^2 - 5000*cos(0.1 * pi * z) - 4000);
end
