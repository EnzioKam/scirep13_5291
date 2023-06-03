function [val] = C04(x, o1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    z = x - o1;
    if g1(z) <= 0 && g2(z) <= 0
        val = f(z);
    end
end

function [val] = f(z)
% Original objective function f
    val = sum(z.^2 - 10*cos(2*pi*z) + 10);
end

function [val] = g1(z)
% Constraint function g(x)
    val = -sum(z.*sin(2*z));
end

function [val] = g2(z)
% Constraint function g(x)
     val = sum(z.*sin(z));
end
