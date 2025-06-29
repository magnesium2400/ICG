function [activityICG, outPairID] = ICG(allData, varargin)

%% Input Parsing
p = inputParser;
addRequired(p, 'allData', @isnumeric);
addParameter(p, 'keepAll', false, @(x) isnumeric(x) || islogical(x));
addParameter(p, 'correlationFunction', @(X) corr(X'), @(x) isa(x, 'function_handle'));
addParameter(p, 'combinationFunction', @(x,y) plus(x,y), @(x) isa(x, 'function_handle'));
parse(p, allData, varargin{:});

allData = p.Results.allData;
keepAll = p.Results.keepAll;
correlationFunction = p.Results.correlationFunction;
combinationFunction = p.Results.combinationFunction;


%% Prelims (ICG level 1 is just the original data)
% Calculate how many ICG iterations are possible
nNeurons = size(allData,1);
ICGsteps = nextpow2(nNeurons+0.5-keepAll)+keepAll; % Catch case if neurons = a power of 2

% The first level is the cellular activity
activityICG = cell(1,ICGsteps);
activityICG{1} = allData;

% Cell ids correspond to order inputted
outPairID = cell(1,ICGsteps);
outPairID{1} = (1:nNeurons)';

clearvars allData


%% Start the ICG process
for ICGlevel = 2:ICGsteps
    fprintf('========= ICG level %2i out of %2i =========\n', ICGlevel, ICGsteps);

    %% Setup for current iteration
    % Grab data
    ICGAct = activityICG{ICGlevel-1};
    nData = size(ICGAct,1);

    % How many pairs can be made
    numPairsOdd = mod(nData, 2);
    numPairsTotal = floor(nData/2) + (keepAll && numPairsOdd);

    % I also know the size of my output data = half data x time
    outdat = nan(numPairsTotal,size(ICGAct,2));

    % This is 2 times number of pairings before
    outPairID{ICGlevel} = nan(numPairsTotal,2^(ICGlevel-1));


    %% Get (sorted) correlations between neurons
    % Calculate correlation matrix
    tic
    rho = correlationFunction(ICGAct);
    fprintf('Correlation computation : %f seconds\n', toc);
    assert(all( size(rho)==nData  ))

    % Sort edges in this correlation matrix - keep track of the edge order
    tic
    upTriMask = triu(true(nData),1);
    [~,sCI] = sort(rho(upTriMask),'descend');
    invCI(sCI) = 1:numel(sCI); %#ok<AGROW>
    [allCIndx, allRIndx] = meshgrid(1:nData);
    allColIndx = allCIndx(upTriMask); % allColIndx = @(x) ceil(sqrt(2*x+0.25)+1/2);
    allRowIndx = allRIndx(upTriMask); % allRowIndx = @(x) x - 0.5*(allColIndx(x)-1)*(allColIndx(x)-2);
    fprintf('Time taken to sort edges: %f seconds\n', toc);

    clearvars  allRIndx allCIndx upTriMask rho


    %% Match the neurons (greedily)
    gdIndex = true(numel(invCI),1);
    pairUsed = false(size(outPairID{ICGlevel-1},1),1);

    tic
    for numPairCnt = 1:floor(nData/2)

        % Text counter
        if ~mod(numPairCnt,250)
            fprintf('%6i pairs out of %6i completed (%2i%%) in %f seconds\n', ...
                numPairCnt, numPairsTotal, fix(100*numPairCnt/numPairsTotal), toc);

        end

        % Find the next pairing (greedily)
        k = find(gdIndex,1,'first');

        % Grabbing data
        % top row index
        rowNew = allRowIndx(sCI(k));
        % top col index
        colNew = allColIndx(sCI(k));

        pairUsed([rowNew, colNew]) = true;

        % ICG process
        % Save data
        outdat(numPairCnt,:) = combinationFunction(ICGAct(rowNew,:), ICGAct(colNew,:));

        % Update the list of original pairs
        outPairID{ICGlevel}(numPairCnt,:) = reshape(outPairID{ICGlevel-1}([rowNew colNew],:)', 1, []);

        % Update list of available neurons to pair (optimised)
        % gdIndex = allRowIndx ~= rowNew & allRowIndx ~= colNew & allColIndx ~= colNew & allColIndx ~= rowNew & gdIndex;
        idx = getIdxToChange(nData , rowNew , colNew , invCI);
        gdIndex(idx) = false;

    end

    if keepAll && numPairsOdd
        missingEdge = find(~pairUsed);
        rowNew = missingEdge; colNew = missingEdge;
        outdat(end,:) = combinationFunction(ICGAct(rowNew,:), ICGAct(colNew,:));
        outPairID{ICGlevel}(end,1:end/2) = reshape(outPairID{ICGlevel-1}(rowNew,:), 1, []);
    end

    fprintf('Time taken to pair edges: %f seconds\n\n', toc);


    % Save all the activity
    activityICG{ICGlevel} = outdat;
    clearvars outdat invCI allRowIndx allColIndx

end



end


function idx = getIdxToChange(nData, rowNew, colNew, invCI)
% basically use the structure of allCIdx and allRIdx to quickly generate the
% indices that will need to be changed
%
% we know the locations of the values to be removed in allRowIdx and allColIdx
% given how structured they are: don't have search for where they are -->
% instead we can directly compute their locations (into a,b,c,d)

trin = @(x) (x .* (x-1))/2;  % triangular number

a = trin(colNew-1)+1:trin(colNew);
b = trin(rowNew-1)+1:trin(rowNew);

x = 0:(nData-rowNew-1);
y = trin(x);
c = trin(rowNew+1) + x.*rowNew + y;

x = 0:(nData-colNew-1);
y = trin(x);
d = trin(colNew+1) + x.*colNew + y;

% remap the indices from the (unsorted) allRIdx/allCIdx space --> into the
% sorted space
idx = invCI([a,b,c,d]);

% return these indices, and then change only those ones to false (in the main
% function sctipt)
end
