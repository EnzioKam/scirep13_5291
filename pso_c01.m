foldername = 'results';
mkdir(foldername);

funcname = 'C01';
mkdir([foldername '/' funcname]);
obj = @C01;
cno = 1;
a = 100; % search space [-a, a]^d

reps = 100;
iters = [10,100,500,1000,2500,5000,7500,10000];
sig = 0.005;
swarmSize = 32;

% Algorithm parameters
malgo = "hmPSO";
algo = "PSO";
w = 1 - 0.729;
c1 = 1.5;
c2 = 1.5;
U1 = @rand;
U2 = @rand;

dims = [10, 30, 50, 100];
for d = dims
    up = a*ones(d, 1);
    lb = -a*ones(d, 1);

    [meanEnd1, varEnd1, progP1, bestrecP1, time1, besttracj1, ...
    meanEnd2, varEnd2, progP2, bestrecP2, time2, besttracj2] = ...
    hmSwarmRunner(obj, up, lb, malgo, reps, iters, sig, ...
                  w, c1, c2, U1, U2, swarmSize, cno);
    
    filename = [foldername '/' funcname '/' convertStringsToChars(malgo) '_' num2str(d)];
    save(filename, 'up', 'lb', 'reps', 'iters', 'sig', 'swarmSize', 'malgo', 'w', 'c1', 'c2', ...
        'U1', 'U2', 'meanEnd1', 'varEnd1', 'progP1', 'bestrecP1', 'time1', 'besttracj1')

    filename = [foldername '/' funcname '/' convertStringsToChars(algo) '_' num2str(d)];
    save(filename, 'up', 'lb', 'reps', 'iters', 'sig', 'swarmSize', 'algo', 'w', 'c1', 'c2', ...
        'U1', 'U2', 'meanEnd2', 'varEnd2', 'progP2', 'bestrecP2', 'time2', 'besttracj2')
    
end
