load('../data/drugData.mat');
load('../data/cellData.mat');

numDrugs = length(drugData);
numCells = length(cellData);

drugCorr = cell(numDrugs, 1);

for d = 1:numDrugs
    drugInfo = drugData{d};
    drugCorr{d}.name    = drugInfo.name;
    drugCorr{d}.targets = drugInfo.targets;
    drugCells = drugInfo.cells;
    drugCorr{d}.cells   = drugCells;

    fprintf('Processing %s, %d\n', drugInfo.name, d);
    % Get a intersection of common gene symbols across cell lines
    commonGenes = [];
    for c = 1:numCells
        idx = find( strcmp(cellData{c}.name, drugCells) );
        if isempty(idx)
            continue;
        end
        if isempty(commonGenes)
            commonGenes = cellData{c}.genes;
        else
            commonGenes = intersect(commonGenes, cellData{c}.genes);
        end
    end

    drugCorr{d}.genes = commonGenes;

    numCommonGenes = length(commonGenes);
    correlations = zeros(numCommonGenes, length(drugCells));
    cellCtrl     = zeros(length(drugCells), 1);

    debug = zeros(1, length(drugCells));

    for c = 1:numCells
        idx = find( strcmp(cellData{c}.name, drugCells) );
        if isempty(idx)
            continue;
        end

        debug(idx) = idx;

        find_index = @(x) find(strcmp(x, cellData{c}.genes));
        indices = arrayfun(find_index, commonGenes);

        curDrugScores = drugInfo.scores(:, idx);
        curCellScores = cellData{c}.scores(:, indices);
        curCtrlScores = cellData{c}.ctl;
        correlations(:, idx) = corr(curCellScores, curDrugScores);
        cellCtrl(idx, 1) = corr(curCtrlScores, curDrugScores);
    end

    debug = sort(debug);
    if sum( debug == 1:length(drugCells) ) ~= length(drugCells)
        disp(debug)
    end

    drugCorr{d}.correlations = correlations;
    drugCorr{d}.ctrl = cellCtrl;
    disp(size(drugCorr{d}.correlations));
    disp(size(drugCorr{d}.ctrl));
end

save('../data/drugCorr', 'drugCorr');

