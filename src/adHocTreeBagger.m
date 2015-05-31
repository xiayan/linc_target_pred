load('./drugAllFeatures.mat')

numTrees = 100;
numNegatives = 100;
numDrugs = length(drugFeatures);
results = zeros(numDrugs, 4);

drugRankingResults = cell(numDrugs, 1);

for d = 1:numDrugs
    XTrain = [];
    YTrain = [];
    drugCells = drugFeatures{d}.cells;

    compatibleDrugs = 0;

    fprintf('Processing %s, %d\n', drugFeatures{d}.name, d);
    for i = 1:numDrugs
        if i == d
            continue
        end

        drugInfo = drugFeatures{i};
        curCells = drugInfo.cells;

        commonCells = intersect(drugCells, curCells);
        if length(drugCells) > length(commonCells)
            continue;
        end


        % ci is the index order of cells in curCells
        ci = arrayfun( @(x) find( strcmp(curCells, x) ), drugCells );

        drugTargets = drugInfo.targets;
        cellMatch = drugInfo.ctrl;
        cellMatch = (cellMatch(ci))';
        geneIDs = drugInfo.genes;

        % ti is a list of target indices in the sorted gene id list
        ti = arrayfun( @(x) find(strcmp(geneIDs, x)), drugTargets, 'UniformOutput', 0 );
        ix = cellfun(@isempty, ti);
        ti(ix) = [];
        ti = cell2mat(ti);

        if isempty(ti)
            continue;
        end

        compatibleDrugs = compatibleDrugs + 1;

        for g = 1:length(ti)
            corrFeat = drugInfo.correlations( ti(g), ci );
            ppiFeat  = drugInfo.ppi_scores( ti(g), ci );
            meanDrugFeat = drugInfo.meanDrugPPI( ti(g), ci );
            maxDrugFeat = drugInfo.maxDrugPPI( ti(g), ci );
            feat = [corrFeat, ppiFeat, meanDrugFeat, maxDrugFeat, cellMatch];
            XTrain = [XTrain; feat];
            YTrain = [YTrain; 1.0];
        end

        numGenes = length(drugInfo.genes);
        indices = 1:numGenes;
        indices(ti) = [];
        ni = randsample(indices, numNegatives); % negative example indices
        for g = 1:length(ni)
            nCorrFeat = drugInfo.correlations( ni(g), ci );
            nPpiFeat  = drugInfo.ppi_scores( ni(g), ci );
            nMeanDrugFeat = drugInfo.meanDrugPPI( ni(g), ci );
            nMaxDrugFeat = drugInfo.maxDrugPPI( ni(g), ci );
            nFeat = [nCorrFeat, nPpiFeat, nMeanDrugFeat, nMaxDrugFeat, cellMatch];
            XTrain = [XTrain; nFeat];
            YTrain = [YTrain; 0.0];
        end
    end

    fprintf('%d cells. %d drugs compatible\n', length(drugCells), compatibleDrugs);

    disp('Training');
    B = TreeBagger(numTrees, XTrain, YTrain, 'Method', 'regression');

    % Testing
    drug = drugFeatures{d};
    dTargets = drug.targets;
    dCellMatch = (drug.ctrl)';
    dGeneIDs = drug.genes;
    dCellMatch = repmat(dCellMatch, length(dGeneIDs), 1);

    % dTi is a list of target indices in the sorted gene id list
    dTi = arrayfun( @(x) find(strcmp(dGeneIDs, x)), dTargets, 'UniformOutput', 0 );
    dix = cellfun(@isempty, dTi);
    dTi(dix) = [];
    dTi = cell2mat(dTi);

    XTest = [drug.correlations, drug.ppi_scores, drug.meanDrugPPI, ...
        drug.maxDrugPPI, dCellMatch];

    disp('Testing');
    preds = B.predict(XTest);

    [s_preds, dIdx] = sort(preds, 'descend');
    disp(s_preds(1:5));
    [~, g_rank] = sort(dIdx);
    target_rank = g_rank(dTi);
    disp(target_rank);
    min_rank = min(target_rank);
    results(d, 1) = min_rank;
    results(d, 2) = length(dGeneIDs);
    results(d, 3) = length(drugCells);
    results(d, 4) = compatibleDrugs;
    disp(min_rank);

    curDrugInfo = drugFeatures{d};
    drugRankingResults{d}.name = curDrugInfo.name;
    drugRankingResults{d}.targets = curDrugInfo.targets;
    drugRankingResults{d}.name = curDrugInfo.name;
    sorted_genes = curDrugInfo.genes(dIdx);
    drugRankingResults{d}.geneRanking = sorted_genes;
    drugRankingResults{d}.minTargetRank = min_rank;

end

save('../mat_results/drugRankingResults', 'drugRankingResults');


