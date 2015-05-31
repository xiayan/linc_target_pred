load('./drugRankingResults.mat')
load('./drugAllFeatures.mat')

iteration = 20000;
perRes = zeros(iteration + 1, 1);
numDrugs = length(drugRankingResults);

for i = 1:iteration
    success = 0;
    fprintf('Iteration %d\n', i);
    for d = 1:numDrugs
        curDrug = drugRankingResults{d};
        numTargets = length(curDrug.targets);
        numGenes = length(curDrug.geneRanking);
        randGene = randperm(numGenes);
        randIdx  = randperm(numGenes, numTargets);
        randTargetsIdx = randGene(randIdx);
        randTargets = drugFeatures{d}.genes(randTargetsIdx);

        dTi = arrayfun( @(x) find(strcmp(curDrug.geneRanking, x)), ...
            randTargets, 'UniformOutput', 0 );
        dix = cellfun(@isempty, dTi);
        dTi(dix) = [];
        dTi = cell2mat(dTi);

        minRank = min(dTi);
        if minRank <= 100
            success = success + 1;
        end
    end
    perRes(i) = success / numDrugs;
end

perRes(iteration + 1) = 0.54;
hist(perRes);

