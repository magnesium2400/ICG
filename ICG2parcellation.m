function out = ICG2parcellation(outPairID, n)

% Find the level needed 
ICGlevel = find(cellfun(@height, outPairID)<=n,1,'first');
ICGparc = outPairID{ICGlevel};
ICGn = size(ICGparc,1);

% Generate the ROI allocations based on that level
p = zerosz(ICGparc)+(1:ICGn)';

% Break apart the later parcels as needed
% fprintf('n=%i\n',n)
% 2*ICGn-n+1
% height(p)
prevOdd = all(isnan( ICGparc(end,end/2+1:end) ));
toBisect = 2*ICGn-n+1-prevOdd; 
p(toBisect:end,end/2+1:end) = p(toBisect:end,end/2+1:end) + (n-ICGn);

% Reshape and return
out = zeros(numel(outPairID{1}),1);
out(ICGparc(~isnan(ICGparc))) = p(~isnan(ICGparc));

end
