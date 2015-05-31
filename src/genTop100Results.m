load('./drugRankingResults.mat');

rIdx = [];

for d = 1:length(drugRankingResults)
    if drugRankingResults{d}.minTargetRank > 100
        rIdx = [rIdx; d];
    end
end

drugRankingResults(rIdx) = [];
drugRankingTop100Results = drugRankingResults;
save('drugRankingTop100Results', 'drugRankingTop100Results');

