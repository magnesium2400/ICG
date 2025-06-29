%% User inputs
pVerts = 60; % parameter for number of verts, use 60 for quick example, use 180 for 32k
nTime = 200; % length of timeseries
nModes = 50; % just used for smoothness of timeseries

doPlot = false; %#ok<*UNRCH>


%% Generate data
[v,f] = sphereMesh(pVerts, 'fib');
s = calc_geometric_eigenmode(struct('vertices', v, 'faces', f), nModes); 
n = noise('powerlaw1', struct('sz', [nModes, nTime], 'alpha', 4));
ts = s.evecs * n; 

if doPlot
videofigs(width(n), @(n) cla, @(n) plotBrain(v,f,[],ts(:,n)), ...
    @(n) clim(minmax(ts,[],'all'))); 
figure; imagesc(ts); colorbar; 
end


%%
t1 = tic;
[activity1, pairID1] = ICG(ts); 
t1 = toc(t1);

t2 = tic; 
[activity2, pairID2] = ICG_original(ts); 
t2 = toc(t2);

clc


%%
% fprintf('Time 1: %f seconds\n', t1);
% fprintf('Time 2: %f seconds\n', t2);
% fprintf('Passed check 1? %i\n', isequal(pairID1, pairID2)); 
% fprintf('Passed check 2? %i\n', isequal(activity1, activity2)); 


%%
out = ...
sprintf('Time 1: %f seconds\nTime 2: %f seconds\nPassed check 1? %i\nPassed check 2? %i\n', ...
    t1, t2, isequal(pairID1, pairID2), isequal(activity1, activity2)) %#ok<NOPTS>
writelines(out, strnow+".txt");



