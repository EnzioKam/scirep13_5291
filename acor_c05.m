foldername = 'results';
mkdir(foldername);

funcname = 'C05';
mkdir([foldername '/' funcname]);
obj = @C05;
cno = 5;
a = 10; % search space [-a, a]^d

reps = 100;
iters = [10,100,500,1000,2500,5000,7500,10000];
sig = 0.005;
swarmSize = 32;

% Algorithm parameters
malgo = "hmACOR";
algo = "ACOR";
q = 0.0001;
zeta = 0.85;

dims = [10, 30, 50, 100];
for d = dims
    up = a*ones(d, 1);
    lb = -a*ones(d, 1);

    [meanEnd1, varEnd1, progP1, bestrecP1, time1, besttracj1, ...
    meanEnd2, varEnd2, progP2, bestrecP2, time2, besttracj2] = ...
    hmSwarmRunner(obj, up, lb, malgo, reps, iters, sig, ...
                  q, zeta, swarmSize, swarmSize/2, cno);
    
    filename = [foldername '/' funcname '/' convertStringsToChars(malgo) '_' num2str(d)];
    save(filename, 'up', 'lb', 'reps', 'iters', 'sig', 'swarmSize', 'malgo', 'q', 'zeta', ...
        'meanEnd1', 'varEnd1', 'progP1', 'bestrecP1', 'time1', 'besttracj1')

    filename = [foldername '/' funcname '/' convertStringsToChars(algo) '_' num2str(d)];
    save(filename, 'up', 'lb', 'reps', 'iters', 'sig', 'swarmSize', 'algo', 'q', 'zeta', ...
        'meanEnd2', 'varEnd2', 'progP2', 'bestrecP2', 'time2', 'besttracj2')
    
end
