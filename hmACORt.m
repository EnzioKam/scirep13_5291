function [result, glbestP, bestvaluerec, bestposrec] = hmACORt(Xini, q, zeta, sampleSize, sig, obj, iters, up, lb)
    % Runs the half-modified ACOR (Ant Colony Optimisation for continuous domains) algorithm.
    %
    % Inputs:
    % Xini: initial location of particles
    % q:
    % zeta:
    % sampleSize:
    % sig: standard deviation of pertubation noise
    % obj: objective function
    % iters: list of iteration numbers to record values
    % up: upper bound
    % lb: lower bound
    % 
    % Outputs:
    % result: position of lowest value of objective function recorded
    % glbestP: lowest values recorded for all iterations
    % bestvaluerec: lowest values recorded of iterations listed in iters
    % bestposrec: positions of lowest values recorded of iterations listed in iters
    %
    % References:
    % TBC
    %

    record_length = length(iters);
    N = iters(record_length);

    bestvaluerec = zeros(record_length, 1);
    tp = 1;

    N = N - 1;
    [d, swarmSize] = size(Xini);
    % X = zeros(d, swarmSize, N);
    % X(:, :, 1) = Xini;
    % V = Vini;
    pbest = Xini;
    bestposrec = zeros(d, record_length);

    pbestvalue = 1:1:swarmSize;
    iterNum = 1;
    for i = 1:1:swarmSize
        pbestvalue(i) = obj(pbest(:, i));
    end

    [sortedpbestvalue, sortOrder] = sort(pbestvalue);
    pbestvalue = sortedpbestvalue;
    pbest = pbest(:, sortOrder);
    glbestvalue = pbestvalue(1);
    glbest(:, iterNum) = pbest(:, 1);
    glbestP = 1:N+1;

    % w = 1/(sqrt(2*pi)*q*swarmSize)*exp(-0.5*(((1:swarmSize)-1)/(q*swarmSize)).^2);
    w = exp(-((1:swarmSize) - 1).^2 / (q^2 * swarmSize^2));
    w = w / sum(w);
    p = cumsum(w);

    while iterNum <= N
        glbestP(iterNum) = glbestvalue;
        
        % Generate the standard deviations
        stds = zeros(d, swarmSize);
        for i = 1:swarmSize
            xPos = pbest(:, i);
            total = zeros(size(xPos));
            for k = 1:swarmSize
                total = total + abs(xPos - pbest(:, k));
            end
            stds(:, i) = zeta * total / (swarmSize - 1);
        end

        % Generate new ants
        newAntsPos = zeros(d, sampleSize);
        newAntsCost = zeros(1, sampleSize);
        for t = 1:sampleSize
            for i = 1:d
                 % Choosing which distribution to generate new ant
                index = find(rand <= p, 1);
                % Create new ant and position
                newAntsPos(i, t) = pbest(i, index) + stds(i, index) * randn;
            end

            % Projection + Pertubation
            if t <= sampleSize / 2 % Apply projection + pertubation on only half of the new ants
                % noise = sig*randn(d, 1); % Random noise for pertubation
                noise = 0.01 * sqrt((sig-2)/sig) * trnd(sig, d, 1); % t-distributed noise
                tmp = max(min(up, newAntsPos(:, t)), lb) + noise;
                newAntsPos(:, t) = max(min(up, tmp), lb);
            else % Normal update without second projection and pertubation
                newAntsPos(:, t) = max(min(up, newAntsPos(:, t)), lb);
            end

            newAntsCost(t) = obj(newAntsPos(:, t));
        end

        allAntsCost = [pbestvalue newAntsCost];
        allAntsPos = [pbest newAntsPos];
        [sortedpbestvalue, sortOrder] = sort(allAntsCost);
        
        pbestvalue = sortedpbestvalue(1:swarmSize);
        pbest = allAntsPos(:, sortOrder(1:swarmSize));

        iterNum = iterNum+1;
        % Update global best value
        glbest(:, iterNum) = glbest(:, iterNum - 1);
        temp = pbestvalue(1);
        if temp < glbestvalue
            glbest(:, iterNum) = pbest(:, 1);
            glbestvalue = temp;
        end
        % Record values based on list in iters
        if iterNum == iters(tp) - 1
            bestvaluerec(tp) = glbestvalue;
            bestposrec(:, tp) = glbest(:, iterNum);
            tp = min(tp+1, record_length);
        end
    end

    glbestP(iterNum) = glbestvalue;
    result = glbest(:, iterNum);

end