function [result, glbestP, bestvaluerec, bestposrec] = hmPSO(Xini, Vini, w, c1, c2, U1, U2, sig, obj, iters, up, lb)
    % Runs the half-modified Particle Swarm Optimisation algorithm.
    %
    % Inputs:
    % Xini: initial location of particles
    % Vini: initial velocity of particles
    % w: inertia weight
    % c1, c2: acceleration constants
    % U1, U2: distribution used for velocity update
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

        for i = 1:swarmSize
            % Velocity update
            V(:, i) = (1 - w)*V(:, i) - c1*U1(d, 1).*(X(:, i, iterNum) - pbest(:, i)) ...
                        - c2*U2(d, 1).*(X(:, i, iterNum) - glbest(:, iterNum));
            
            % Projection + Pertubation
            if i <= swarmSize/2 % Apply projection + pertubation on only half of the particles
                noise = sig*randn(d, 1); % Random noise for pertubation
                tmp = max(min(up,  X(:, i, iterNum)+V(:, i)),  lb) + noise;
                X(:, i, iterNum+1) = max(min(up, tmp), lb);
            else % Normal update without second projection and pertubation
                X(:, i, iterNum+1) = max(min(up,  X(:, i, iterNum)+V(:, i)),  lb);
            end
            
            % Update particle best value and positions
            if obj(X(:, i, iterNum+1)) < pbestvalue(i)
                pbest(:, i) = X(:, i, iterNum+1);
                pbestvalue(i) = obj(pbest(:, i));
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