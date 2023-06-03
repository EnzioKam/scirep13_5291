function [val] = C05(x, o1, m1, m2)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = x - o1;
    y = m1 * z;
    w = m2 * z;
    if g1(y) <= 0 && g2(w) <= 0
        val = f(z);
    end
end

function [val] = f(z)
% Original objective function f
    v = 100*(z.^2 - circshift(z, -1)).^2 + (z - 1).^2;
    val = sum(v) - v(end);
end

function [val] = g1(y)
% Constraint function g(x)
    val = sum(y.^2 - 50*cos(2*pi*y)- 40);
end

function [val] = g2(w)
% Constraint function g(x)
     val = sum(w.^2 - 50*cos(2*pi*w)- 40);
end
