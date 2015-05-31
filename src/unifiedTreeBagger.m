function [AllResult, forest] = unifiedTreeBagger
    load('./drugAllFeatures.mat');
    % Do 10-fold cross validation
    numDrugs = numel(drugFeatures);
    partition = 1 : round(numDrugs/10) : numDrugs;
    partition(end) = numDrugs;

    global cellLines;
    cellLines = {'MCF7', 'HA1E', 'A549', 'HT29', 'VCAP', 'A375', 'HCC515'};

    global numNegs;
    numNegs = 100;

    global numTrees;
    numTrees = 1;

    global numBags;
    numBags = 3500;

    AllResult = [];
    for p = 1:length(partition) - 1
        % create the index for testing set
        fprintf('Iteration %d\n', p);
        bIdx = partition(p);
        eIdx = partition(p+1) - 1;
        if p == length(partition) - 1
            eIdx = eIdx + 1;
        end

        trainingIdx = 1:numDrugs;
        trainingIdx( bIdx:eIdx ) = [];
        testingIdx = bIdx:eIdx;

        forest = genForest(trainingIdx, drugFeatures);
        result = forestPredict(forest, testingIdx, drugFeatures);
        AllResult = [AllResult; result];
    end

    for d = 1:numDrugs
        fprintf('%s\t%f\n', drugFeatures{d}.name, AllResult(d));
    end
end

function forest = genForest(trainingIdx, drugFeatures)
    global cellLines;
    global numBags;
    global numTrees;
    global numNegs;
    forest = cell(numBags, 1);
    iter = 1;
    while iter <= numBags
        if mod(iter, 50) == 0
            fprintf('Tree %d\n', iter);
        end
        cells = cellLines(randperm(7, 4));
        XTrain = [];
        YTrain = [];

        for d = 1:length(trainingIdx)
            curDrug = drugFeatures{trainingIdx(d)};
            commonCells = intersect(curDrug.cells, cells);
            if length(commonCells) < 4
                continue;
            end

            % ci is the index order of cells in curDrug
            ci = arrayfun( @(x) find( strcmp(curDrug.cells, x) ), cells );

            drugTargets = curDrug.targets;
            cellMatch   = curDrug.ctrl;
            cellMatch   = (cellMatch(ci))';
            geneIDs     = curDrug.genes;

            % ti is a list of target indices in the sorted gene id list
            ti = arrayfun( @(x) find(strcmp(geneIDs, x)), drugTargets, ...
                'UniformOutput', 0 );
            ix = cellfun(@isempty, ti);
            ti(ix) = [];
            ti = cell2mat(ti);

            if isempty(ti)
                continue;
            end

            for g = 1:length(ti)
                corrFeat = curDrug.correlations( ti(g), ci );
                ppiFeat  = curDrug.ppi_scores( ti(g), ci );
                meanDrugFeat = curDrug.meanDrugPPI( ti(g), ci );
                maxDrugFeat = curDrug.maxDrugPPI( ti(g), ci );
                feat = [corrFeat, ppiFeat, meanDrugFeat, maxDrugFeat, cellMatch];
                XTrain = [XTrain; feat];
                YTrain = [YTrain; 1.0];
            end

            numGenes = length(geneIDs);
            indices = 1:numGenes;
            indices(ti) = [];
            ni = randsample(indices, numNegs); % negative example indices
            for g = 1:length(ni)
                nCorrFeat = curDrug.correlations( ni(g), ci );
                nPpiFeat  = curDrug.ppi_scores( ni(g), ci );
                nMeanDrugFeat = curDrug.meanDrugPPI( ni(g), ci );
                nMaxDrugFeat = curDrug.maxDrugPPI( ni(g), ci );
                nFeat = [nCorrFeat, nPpiFeat, nMeanDrugFeat, nMaxDrugFeat, cellMatch];
                XTrain = [XTrain; nFeat];
                YTrain = [YTrain; 0.0];
            end
        end
        if (numel(YTrain) < 500)
            disp(size(XTrain));
            disp(size(YTrain));
        end
        B = TreeBagger(numTrees, XTrain, YTrain, 'Method', 'regression');
        forest{iter}.cells = cells;
        forest{iter}.B = B;
        iter = iter + 1;
    end
end


function result = forestPredict(forest, testingIdx, drugFeatures)
    result = zeros(length(testingIdx), 1);
    for d = 1:length(testingIdx)
        curDrug = drugFeatures{testingIdx(d)};
        curTarget = curDrug.targets;
        curCellMatch = (curDrug.ctrl)';
        curGeneIDs = curDrug.genes;
        curCellMatch = repmat(curCellMatch, length(curGeneIDs), 1);

        % dTi is a list of target indices in the sorted gene id list
        dTi = arrayfun( @(x) find(strcmp(curGeneIDs, x)), curTarget, ...
            'UniformOutput', 0 );
        dix = cellfun(@isempty, dTi);
        dTi(dix) = [];
        dTi = cell2mat(dTi);

        curRes  = zeros(length(curDrug.genes), 1);

        numCompatible = 0;
        for f = 1:numel(forest)
            cells = forest{f}.cells;
            commonCells = intersect(curDrug.cells, cells);
            if length(commonCells) < 4
                continue;
            end
            numCompatible = numCompatible + 1;

            % ci is the index order of cells in current forest
            ci = arrayfun( @(x) find( strcmp(curDrug.cells, x) ), cells );

            XTest = [curDrug.correlations(:,ci), curDrug.ppi_scores(:,ci), ...
                curDrug.meanDrugPPI(:,ci), curDrug.maxDrugPPI(:,ci), ...
                curCellMatch(:,ci)];
            preds = forest{f}.B.predict(XTest);
            curRes = curRes + preds;
        end
        fprintf('Compatible forests for drug %d: %d\n', d, numCompatible);

        [~, dIdx] = sort(curRes, 'descend');
        [~, g_rank] = sort(dIdx);
        target_ranks = g_rank(dTi);
        min_rank = min(target_ranks);
        result(d) = min_rank;
    end
end


