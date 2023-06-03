function [result, glbestP, bestvaluerec, bestposrec] = mCSO(Xini, Vini, phi, sig, obj, iters, up, lb)
    % Runs the modified Competitive Swarm Optimiser algorithm.
    %
    % Inputs:
    % Xini: initial location of particles
    % Vini: initial velocity of particles
    % phi: social factor parameter
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
    X = zeros(d, swarmSize, N);
    X(:, :, 1) = Xini;
    V = Vini;
    pbest = Xini;
    bestposrec = zeros(d, record_length);

    pbestvalue = 1:1:swarmSize;
    iterNum = 1;
    for i = 1:1:swarmSize
        pbestvalue(i) = obj(pbest(:, i));
    end

    [glbestvalue, ind] = min(pbestvalue);
    glbest(:, iterNum) = pbest(:, ind);
    glbestP = 1:N+1;

    while iterNum <= N
        glbestP(iterNum) = glbestvalue;

        % Generate random pairs for competition
        rlist = randperm(swarmSize);
        rpairs = [rlist(1:ceil(swarmSize/2)); rlist(floor(swarmSize/2) + 1:swarmSize)]';

        % Calculate the centre position and determine the winners
        positions = X(:, :, iterNum);

        comp = (pbestvalue(rpairs(:,1)) > pbestvalue(rpairs(:,2)))';
        losers = comp.*rpairs(:,1) + ~comp.*rpairs(:,2);
        winners = ~comp.*rpairs(:,1) + comp.*rpairs(:,2);

        % Random variables for velocity update
        randco1 = rand(d, ceil(swarmSize/2));
        randco2 = rand(d, ceil(swarmSize/2));
        randco3 = rand(d, ceil(swarmSize/2));

        % Update the velocity and position of losers
        centre = mean(positions, 2) * ones(1, ceil(swarmSize/2));
        V(:,losers) = randco1.*V(:,losers) ...,
                    + randco2.*(positions(:, winners) - positions(:, losers)) ...,
                    + phi*randco3.*(centre - positions(:, losers));
        positions(:,losers) = positions(:,losers) + V(:,losers);

        % Projection + Pertubation
        for i = 1:ceil(swarmSize/2)
            tmp = min(max(positions(:, losers(i)), lb), up) + sig*randn(d, 1);
            positions(:, losers(i)) = min(max(tmp, lb), up);
        end
        X(:, :, iterNum+1) = positions;

        % Update particle best value and positions
        for i = 1:length(losers)
            j = losers(i);
            if obj(X(:, j, iterNum+1)) < pbestvalue(j)
                pbest(:, j) = X(:, j, iterNum+1);
                pbestvalue(j) = obj(pbest(:, j));
            end
        end

        iterNum = iterNum+1;
        % Update global best value
        glbest(:, iterNum) = glbest(:, iterNum - 1);
        [temp, ind] = min(pbestvalue);
        if temp < glbestvalue
            glbest(:, iterNum) = pbest(:, ind);
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