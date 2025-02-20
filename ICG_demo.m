%%% Shared variables

[v,f] = read_vtk(fullfile(fileparts(which('ICG')), '..','..', 'DATA','D_HumanStandardSurface','fsaverage5_10k_midthickness-lh.vtk'));
v = v'; f = f';

% m = readmatrix(fullfile(fileparts(which('ICG')), '..','..', 'DATA','D_HumanStandardSurface','fsaverage5_10k_cortex-lh_mask.txt'));
% [v,f,~,m] = trimExcludedRois(v, f, m);

[~,f,v] = decreaseRegularPatch(f,2,[],v);

s = calc_geometric_eigenmode(struct('vertices', v, 'faces', f), 100);


%% Run ICG
t1 = tic; 
[~,outPairID] = ICG(s.evecs,1);
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
    plotBrain(s.vertices, s.faces, rois, rois, 'BoundaryMethod', 'faces');
    colormap( ax, lines(range(rois)+1) );
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

