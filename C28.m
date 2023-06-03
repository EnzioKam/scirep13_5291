function [val] = C28(x, o1, m1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = m1 * (x - o1);
    if g1(z) <= 0 && g2(z) <= 0
        val = f(z);
    end
end

function [val] = f(z)
% Original objective function f
    val = sum(abs(z).^0.5 + 2*sin(z.^3));
end

function [val] = g1(z)
% Constraint function g(x)
    d = length(z);
    v = -10*exp(-0.2 * sqrt(z.^2 + circshift(z.^2, -1)));
    val = sum(v) - v(end) + (d - 1) * 10 / exp(-5);
end

function [val] = g2(z)
% Constraint function g(x)
    d = length(z);
    val = sum(sin(2*z).^2) - 0.5*d;
end
