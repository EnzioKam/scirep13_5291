function [val] = C13(x, o1)
% C01 from CEC 2017
% Modified objective unction f+
    val = Inf;
    y = x - o1;
    if g1(y) <= 0 && g2(y) <= 0 && g3(y) <= 0
        val = f(y);
    end
end

function [val] = f(y)
% Original objective function f
    v = 100*(y.^2 - circshift(y, -1)).^2 + (y - 1).^2;
    val = sum(v) - v(end);
end

function [val] = g1(y)
% Constraint function g(x)
    val = sum(y.^2 - 10*cos(2*pi*y) + 10) - 100;
end

function [val] = g2(y)
% Constraint function g(x)
    val = sum(y) - 2*length(y);
end

function [val] = g3(y)
% Constraint function g(x)
    val = 5 - sum(y);
end
