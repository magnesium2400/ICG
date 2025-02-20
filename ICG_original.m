function [activityICG,outPairID] = ICG_original(allData)

%%
%Input

%allData -  neuronal timeseries data comes in rows - neurons/ cols - time

%Ouput:

%activityICG - ICG activity for level l
%outPairID - the ids of original neurons grouped at each level

%Brandon Munn, 19/10/21


%% Prelims (ICG level 1 is just the original data)
%Calculate how many ICG iterations are possible
ICGsteps = nextpow2(size(allData,1)+0.5); % Catch case if more neurons than a power of 2

%The first level is the cellular activity
activityICG = cell(1,ICGsteps);
activityICG{1} = allData;

%Cell ids correspond to order inputted
outPairID = cell(1,ICGsteps);
outPairID{1} = (1:size(allData,1))';

clearvars allData


%% Start the ICG process
for ICGlevel = 2:ICGsteps
    fprintf('=== ICG level is %2i out of %2i ===\n', ICGlevel, ICGsteps);

    %% Setup for current iteration
    %Grab data
    ICGAct = activityICG{ICGlevel-1};
    nData = size(ICGAct,1);

    %How many pairs can be made
    numPairsTotal = floor(nData/2);

    % I also know the size of my output data = half data by time
    outdat = nan(numPairsTotal,size(ICGAct,2));

    %this is 2 times number of pairings before
    outPairID{ICGlevel} = nan(numPairsTotal,2^(ICGlevel-1));


    %% Get (sorted) correlations between neurons
    %Calculate correlation matrix
    tic
    rho = corr(ICGAct');
    fprintf('Correlation computation : %f seconds\n', toc);

    %Sort edges in this correlation matrix - keep track of the edge order
    upTriMask = triu(true(nData),1);
    [~,sCI] = sort(rho(upTriMask),'descend');
    [allCIndx, allRIndx] = meshgrid(1:nData); 
    allRIndx    = allRIndx(upTriMask);  allCIndx    = allCIndx(upTriMask); 
    allRowIndx  = allRIndx(sCI);        allColIndx  = allCIndx(sCI);
    
    clearvars allRIndx allCIndx sCI upTriMask rho


    %% Match the neurons (greedily)
    gdIndex = true(numel(allRowIndx),1);

    tic
    for numPairCnt = 1:numPairsTotal

        %Text counter
        if ~mod(numPairCnt,250)
            fprintf('%6i pairs out of %6i completed (%2i%%) in %f seconds\n', ...
                numPairCnt, numPairsTotal, fix(100*numPairCnt/numPairsTotal), toc);

        end

        %Find the next pairing (greedily)
        k = find(gdIndex,1,'first');

        %Grabbing data
        %top row index
        rowNew = allRowIndx(k);
        %top col index
        colNew = allColIndx(k);

        %ICG process
        % Save data
        outdat(numPairCnt,:) = ICGAct(rowNew,:) + ICGAct(colNew,:);

        %Update the list of original pairs
        outPairID{ICGlevel}(numPairCnt,:) = reshape(outPairID{ICGlevel-1}([rowNew colNew],:)', 1, []);

        %Update list of available neurons to pair (optimised)
        gdIndex = allRowIndx ~= rowNew & allRowIndx ~= colNew & allColIndx ~= colNew & allColIndx ~= rowNew & gdIndex;

    end
    fprintf('ICG level %2i completed in %f seconds\n\n', ICGlevel, toc);


    %Save all the activity
    activityICG{ICGlevel} = outdat;
    clearvars outdat

end



end
