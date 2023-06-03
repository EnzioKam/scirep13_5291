function [meanEnd1, varEnd1, progP1, bestrecP1, time1, besttracj1, ...
    meanEnd2, varEnd2, progP2, bestrecP2, time2, besttracj2] ...
    = hmSwarmRunner(obj, up, lb, algo, reps, iters, sig, varargin)
    % Driver function to run comparison for half-modified swarm algorithms
    % with the original swarm algorithms
    %
    % TO BE UPDATED

    if length(up) ~= length(lb)
        error("Dimension mismatch between up (upper bound) and lb (lower bound)")
    end

    record_length = length(iters);
    N = iters(record_length);

    progP1 = zeros(1, N);
    progP2 = zeros(1, N);
    bestrecP1 = zeros(record_length, reps);
    bestrecP2 = zeros(record_length, reps);

    % For rotate and translate
    d = length(up);
    % a = up(1);

    fprintf("Modified algorithm to compare: %s\n", algo)
    time1 = zeros(1, reps);
    time2 = zeros(1, reps);

    if algo == "hmPSO"
        if length(varargin) < 5
            excp = MException('MATLAB:notEnoughInputs','Not enough input arguments.');
            throw(excp);
        end

        w = varargin{1};
        c1 = varargin{2};
        c2 = varargin{3};
        U1 = varargin{4};
        U2 = varargin{5};

        swarmSize = 32;
        if length(varargin) > 5
            swarmSize = varargin{6};
            cno = varargin{7};
        end

        A = repmat(lb, 1, swarmSize);
        B = repmat(up, 1, swarmSize);
        [b,~] = size(A);
        besttracj1 = zeros(b,record_length,reps);
        besttracj2 = zeros(b,record_length,reps);

        filename = ['data/C' num2str(cno) '_' num2str(d)];
        if cno == 1 || cno == 4 || cno == 20
            load(filename, 'tmatrix1')
        elseif cno == 2
            load(filename, 'tmatrix1', 'rmatrix1')
        elseif cno == 5
            load(filename, 'tmatrix1', 'rmatrix1', 'rmatrix2')
        else
            error("Function not implemented for hmPSO");
        end

        for i = 1:reps
            Xini = unifrnd(A, B);
            Vini = unifrnd(-(B-A), B-A);

            newobj = obj;
            if cno == 1 || cno == 4 || cno == 20
                o1 = tmatrix1(:, i);
                newobj = @(x) obj(x, o1);
            elseif cno == 2
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                newobj = @(x) obj(x, o1, r1);
            elseif cno == 5
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                r2 = rmatrix2(:, :, i);
                newobj = @(x) obj(x, o1, r1, r2);
            end

            t = cputime;
            [~,glbestP1,bestrecP1(:,i), besttracj1(:,:,i)] = hmPSO(Xini,Vini,w,c1,c2,U1,U2,sig,newobj,iters,up,lb);
            time1(i) = cputime - t;
            progP1 = progP1 + glbestP1;

            t = cputime;
            [~,glbestP2,bestrecP2(:,i), besttracj2(:,:,i)] = mPSO(Xini,Vini,w,c1,c2,U1,U2,0,newobj,iters,up,lb);
            time2(i) = cputime - t;
            progP2 = progP2 + glbestP2;

            if mod(i, 10) == 0
                fprintf("%g repetitions completed.\n", i)
            end
        end

    elseif algo == "hmBAT"
        if length(varargin) < 4
            excp = MException('MATLAB:notEnoughInputs','Not enough input arguments.');
            throw(excp);
        end

        Qmax = varargin{1};
        Qmin = varargin{2};
        r = varargin{3};
        AA = varargin{4};

        swarmSize = 32;
        if length(varargin) > 4
            swarmSize = varargin{5};
            cno = varargin{6};
        end

        A = repmat(lb, 1, swarmSize);
        B = repmat(up, 1, swarmSize);
        [b,~] = size(A);
        besttracj1 = zeros(b,record_length,reps);
        besttracj2 = zeros(b,record_length,reps);

        filename = ['data/C' num2str(cno) '_' num2str(d)];
        if cno == 1 || cno == 4 || cno == 20
            load(filename, 'tmatrix1')
        elseif cno == 2
            load(filename, 'tmatrix1', 'rmatrix1')
        elseif cno == 5
            load(filename, 'tmatrix1', 'rmatrix1', 'rmatrix2')
        else
            error("Function not implemented for hmPSO");
        end

        for i = 1:reps
            Xini = unifrnd(A, B);
            Vini = unifrnd(-(B-A), B-A);

            newobj = obj;
            if cno == 1 || cno == 4 || cno == 20
                o1 = tmatrix1(:, i);
                newobj = @(x) obj(x, o1);
            elseif cno == 2
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                newobj = @(x) obj(x, o1, r1);
            elseif cno == 5
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                r2 = rmatrix2(:, :, i);
                newobj = @(x) obj(x, o1, r1, r2);
            end

            t = cputime;
            [~,glbestP1,bestrecP1(:,i), besttracj1(:,:,i)] = hmBAT(Xini,Vini,Qmax,Qmin,r,AA,sig,newobj,iters,up,lb);
            time1(i) = cputime - t;
            progP1 = progP1 + glbestP1;

            t = cputime;
            [~,glbestP2,bestrecP2(:,i), besttracj2(:,:,i)] = mBAT(Xini,Vini,Qmax,Qmin,r,AA,0,newobj,iters,up,lb);
            time2(i) = cputime - t;
            progP2 = progP2 + glbestP2;

            if mod(i, 10) == 0
                fprintf("%g repetitions completed.\n", i)
            end
        end

    elseif algo == "hmACOR"
        if length(varargin) < 2
            excp = MException('MATLAB:notEnoughInputs','Not enough input arguments.');
            throw(excp);
        end

        q = varargin{1};
        zeta = varargin{2};

        swarmSize = 32;
        if length(varargin) > 2
            swarmSize = varargin{3};
        end

        sampleSize = swarmSize / 2;
        if length(varargin) > 3
            sampleSize = varargin{4};
            cno = varargin{5};
        end

        A = repmat(lb, 1, swarmSize);
        B = repmat(up, 1, swarmSize);
        [b,~] = size(A);
        besttracj1 = zeros(b,record_length,reps);
        besttracj2 = zeros(b,record_length,reps);

        filename = ['data/C' num2str(cno) '_' num2str(d)];
        if cno == 1 || cno == 4 || cno == 20
            load(filename, 'tmatrix1')
        elseif cno == 2
            load(filename, 'tmatrix1', 'rmatrix1')
        elseif cno == 5
            load(filename, 'tmatrix1', 'rmatrix1', 'rmatrix2')
        else
            error("Function not implemented for hmPSO");
        end

        for i = 1:reps
            Xini = unifrnd(A, B);

            newobj = obj;
            if cno == 1 || cno == 4 || cno == 20
                o1 = tmatrix1(:, i);
                newobj = @(x) obj(x, o1);
            elseif cno == 2
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                newobj = @(x) obj(x, o1, r1);
            elseif cno == 5
                o1 = tmatrix1(:, i);
                r1 = rmatrix1(:, :, i);
                r2 = rmatrix2(:, :, i);
                newobj = @(x) obj(x, o1, r1, r2);
            end

            t = cputime;
            [~,glbestP1,bestrecP1(:,i), besttracj1(:,:,i)] = hmACOR(Xini,q,zeta,sampleSize,sig,newobj,iters,up,lb);
            time1(i) = cputime - t;
            progP1 = progP1 + glbestP1;

            t = cputime;
            [~,glbestP2,bestrecP2(:,i), besttracj2(:,:,i)] = mACOR(Xini,q,zeta,sampleSize,0,newobj,iters,up,lb);
            time2(i) = cputime - t;
            progP2 = progP2 + glbestP2;

            if mod(i, 10) == 0
                fprintf("%g repetitions completed.\n", i)
            end
        end

    else
        error("Invalid algorithm specified")
    end

    % Get average over all repetitions
    progP1 = progP1 / reps;
    progP2 = progP2 / reps; 

    meanEnd1 = mean(bestrecP1(end, :));
    varEnd1 = var(bestrecP1(end, :));

    meanEnd2 = mean(bestrecP2(end, :));
    varEnd2 = var(bestrecP2(end, :));

    fprintf("Run complete. %g repetitions completed.\n", reps)
    fprintf('cno = %d\n', cno)
    fprintf('d = %d\n', d)
    fprintf('Iterations \t Avg global min-value\n')
    fprintf('-----Modified-----\n')
    fprintf('%g: \t\t\t %g\n', [iters; mean(bestrecP1, 2)'])
    fprintf('-----Original-----\n')
    fprintf('%g: \t\t\t %g\n', [iters; mean(bestrecP2, 2)'])
    fprintf('-----End of results-----\n')
end