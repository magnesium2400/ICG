%%% Shared variables

g = gifti(fetch_atlas('fsaverage', '10k', 'L', 'pial')); 
v = g.vertices; f = g.faces; 

m = gifti(fetch_atlas('fsaverage', '10k', 'L', '', 'desc', 'nomedialwall')).cdata; 
[v,f,~,m] = trimExcludedRois(v, f, m);

% [~,f,v] = decreaseRegularPatch(f,2,[],v);

s = calc_geometric_eigenmode(struct('vertices', v, 'faces', f), 100);


%% Run ICG

correlationFunction = @(x) x * diag(s.evals) * x';
% combinationFunction = @(x,y) x + y;

t1 = tic; 
[~,outPairID] = ICG(s.evecs, 'keepAll', true, 'correlationFunction', correlationFunction);
fprintf('Time for 1st run is %f seconds\n\n', toc(t1)); 

% t1 = tic; 
% [~,outPairID_original] = ICG_original(s.evecs);
% fprintf('Time for 2nd run is %f seconds\n\n', toc(t1));
% assert(isequaln(outPairID_original, outPairID))

clearvars t1; 


%% Plot parcels (with <500 ROIs)

figuremax;
tiledlayout('flow');

% for ii = 1:length(outPairID)
for ii = find(cellfun(@height, outPairID) < 500)

    r = outPairID{ii};

    rois = zeros(height(s.vertices),1);
    x = repmat((1:height(r))', 1, width(r));
    rois(r(~isnan(r))) = x(~isnan(r));

    ax = nexttile;
    plotBrain(s.vertices, s.faces, rois, rois, 'BoundaryMethod', 'faces', 'BoundaryColor', [1 1 1]);
    colormap( ax, distinguishable_colors(range(rois)+1) );
    title(ax, sprintf('%i parcels', height(r)))

end


%% Check interpolation between default values (for number of parcels)

parcellationsToPlot = 20:-1:10; 

figuremax;
tiledlayout('flow');
temp = [];

% for ii = 1:length(outPairID)
for ii = parcellationsToPlot

    rois = ICG2parcellation(outPairID, ii);
    temp(:,end+1) = rois; 
    % range(rois)

    ax = nexttile;
    plotBrain(s.vertices, s.faces, rois, rois, 'BoundaryMethod', 'faces');
    colormap( ax, lines(range(rois)+1) );
    title(ax, sprintf('%i parcels', ii))

end

imagesc(nexttile,diff(temp,[],2)); 

