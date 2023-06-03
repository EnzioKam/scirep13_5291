function [val] = C22(x, o1, m1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = m1 * (x - o1);
    if g1(z) <= 0 && g2(z) <= 0 && g3(z) <= 0
        val = f(x, z);
    end
end

function [val] = f(x, z)
% Original objective function f
    t = (z.^2 - circshift(x, -1)).^2;
    val = sum(100*(t) + (z - 1).^2);
end

function [val] = g1(z)
% Constraint function g(x)
    val = sum(z.^2 - 10 * cos(2*pi*z) + 10) - 100;
end

function [val] = g2(z)
% Constraint function g(x)
    val = sum(z) - 2*length(z);
end

function [val] = g3(z)
% Constraint function g(x)
    val = 5 - sum(z);
end
