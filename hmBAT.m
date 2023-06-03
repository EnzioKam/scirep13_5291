function [result, glbestP, bestvaluerec, bestposrec] = hmBAT(Xini,Vini,Qmax,Qmin,r,A,sig,obj,iters,up,lb)
    % Runs the half-modified Bat Algorithm.
    %
    % Inputs:
    % Xini: initial location of particles
    % Vini: initial velocity of particles
    % Qmax: frequency maximum parameter
    % Qmin: frequency minimum parameter
    % r: pulse emission rate parameter
    % A: loudness parameter
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

    Q = zeros(swarmSize, 1);

    while iterNum <= N
        glbestP(iterNum) = glbestvalue;

        for i = 1:swarmSize
            % Obtain current pulse frequency
            Q(i) = Qmin + (Qmin - Qmax) * rand;
            
            % Velocity update
            V(:, i) = V(:, i) + (X(:,i,iterNum) - glbest(:,iterNum)) * Q(i);
            tmp = X(:,i,iterNum) + V(:,i);

            if rand > r
                % Generate new position around the current best position
                tmp = glbest(:,iterNum) + 0.001*randn(d,1);
            end

            % Projection + Pertubation
            if i <= swarmSize/2 % Apply projection + pertubation on only half of the particles
                tmp = max(min(max(min(up, tmp), lb) + sig*randn(d,1), up), lb);
            else % Normal update without second projection and pertubation
                tmp = max(min(up, tmp), lb);
            end
            

            fNew = obj(tmp); % Evaluate new solution

            % Update if the solution improves and not too loud
            % if fNew < pbestvalue(i) && rand < A
            if fNew < obj(X(:, i, iterNum)) && rand < A
                X(:,i,iterNum+1) = tmp;
            else
                X(:,i,iterNum+1) = X(:,i,iterNum);
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