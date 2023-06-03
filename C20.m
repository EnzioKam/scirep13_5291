function [val] = C20(x, o1)
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
    yplus = circshift(y, -1);
    t = sqrt(y.^2 + yplus.^2);
    num = sin(t).^2 - 0.5;
    denom = (1 + 0.001*(t)).^2;
    val = sum(0.5 + num ./ denom);
end

function [val] = g1(y)
% Constraint function g(x)
    val = cos(sum(y))^2 - 0.25*cos(sum(y)) - 0.125;
end

function [val] = g2(y)
% Constraint function g(x)
    val = exp(cos(sum(y))) - exp(0.25);
end
