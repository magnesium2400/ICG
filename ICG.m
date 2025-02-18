function [activityICG,outPairID] =  ICG(allData)

%%
%Input

%allData -  neuronal timeseries data comes in rows - neurons/ cols - time

%Ouput:

%activityICG - ICG activity for level l
%outPairID - the ids of original neurons grouped at each level

%Brandon Munn, 19/10/21 


%%
%Calculate how many ICG iterations are possible
ICGsteps = nextpow2(size(allData,1));

%Catch case if more neurons than a power of 2
if size(allData,1) == 2^ICGsteps
    ICGsteps = ICGsteps + 1;
end

%The first level is the cellular activity
activityICG = cell(1,ICGsteps);
activityICG{1} = allData;

%Cell ids correspond to order inputted
outPairID = cell(1,ICGsteps);
outPairID{1}(1:size(allData,1),1) = 1:size(allData,1);

clearvars allData
%% Start the ICG process


for ICGlevel = 2:ICGsteps
    fprintf('=== ICG level is %2i out of %2i ===\n', ICGlevel, ICGsteps); 
    
    %Grab data
    ICGAct = activityICG{ICGlevel-1};
    nData = size(ICGAct,1);
    
    
    %Calculate correlation matrix
    tic
    rho = corr(ICGAct');
    fprintf('Correlation computation : %f seconds\n', toc);
    
    C = triu(rho,1);
    clearvars rho
    
      
    %Grab just the upper triangle
    %Find the indices of each element
    upTriIndx = find(triu(true(size(C)),1)>0);
    [allRIndx, allCIndx] = ind2sub(size(C),upTriIndx);
    
    %now it is a vector of correlation values
    C = C(upTriIndx);
    clearvars upTriIndx
    
    %Sort the correlation matrix
    [~,sCI] = sort(C,'descend');
    
    
    %resort row/col indices and turn into pairings 
    allRowIndx = allRIndx(sCI);
    allColIndx = allCIndx(sCI);
    clearvars allRIndx allCIndx C sCI
    
   
    %How many pairs can be made
    numPairsTotal = floor(nData/2);
    
    
    % I also know the size of my output data = half data by time
    outdat = nan(numPairsTotal,size(ICGAct,2));
    
    %to save pairings
    outPair = nan(numPairsTotal,2);
    
    %this is 2 times number of pairings before
    outPairID{ICGlevel} = nan(numPairsTotal,2^(ICGlevel-1));

    %% Index counter 
    k = 1;
    gdIndex = true(numel(allRowIndx),1);
  
    tic
    for numPairCnt = 1:numPairsTotal
         
        %Text counter
        if ~mod(numPairCnt,250)
            fprintf('%6i pairs out of %6i completed (%3i %%) in %f seconds\n', ...
                numPairCnt, numPairsTotal, round(100*numPairCnt/numPairsTotal), toc);

        end
        
        %Grabbing data
        %top row index
        rowNew = allRowIndx(k);
        %grab corresponding data
        datRow = ICGAct(rowNew,:);
        
        %top col index
        colNew = allColIndx(k);
        %grab corresponding data
        datCol = ICGAct(colNew,:);
        
        
        %ICG process
        %Sum the data
        SpikeAct = datRow+datCol;
        
        % Save data
        outdat(numPairCnt,:) = SpikeAct;
        outPair(numPairCnt,:) = [rowNew colNew];
        
        %Update the list of original pairs
        tempPairs = outPairID{ICGlevel-1}([rowNew colNew],:);
        outPairID{ICGlevel}(numPairCnt,:) = tempPairs(:);
        
        %Update list of available neurons to pair (optimised) 
        gdIndex = allRowIndx ~= rowNew & allRowIndx ~= colNew & allColIndx ~= colNew & allColIndx ~= rowNew & gdIndex; 
        
        %Find the next pairing (greedily)
        k = find(gdIndex,1,'first');

    
        
    end
    fprintf('ICG level %2i completed in %f seconds\n\n', ICGlevel, toc); 
    
    
    %Save all the activity
    activityICG{ICGlevel} = outdat;
    clearvars outdat outPair
    
end



end
