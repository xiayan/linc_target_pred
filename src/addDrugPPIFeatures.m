load('../data/drugCorrPPI.mat');
drugFeatures = drugCorrPPI;
load('../data/drugData.mat');

load('../data/ppi.mat');
proteins = fieldnames(ppi);

lmSmbols = loadjson('../data/lm_symbols.json');
findIndex = @(x) find(strcmp(lmSmbols, x));

numDrugs = length(drugFeatures);
% numDrugs = 1;

for d = 1:numDrugs
    drugInfo = drugFeatures{d};
    gene_ids = drugInfo.genes;
    numGenes = length(gene_ids);
    drugCells = drugInfo.cells;
    numCells = length(drugCells);
    meanDrugPPI = zeros(numGenes, numCells);
    maxDrugPPI  = zeros(numGenes, numCells);

    fprintf('Processing %s, %d\n', drugInfo.name, d);
    for g = 1:numGenes
        if sum( strcmp(proteins, gene_ids{g}) ) == 0
            meanDrugPPI(g, :) = 0.0;
            maxDrugPPI(g, :)  = 0.0;
            continue;
        end

        partners = ppi.(gene_ids{g});
        pIdx = cell2mat( arrayfun(findIndex, partners, 'UniformOutput', false) );

        if isempty(pIdx)
            meanDrugPPI(g, :) = 0.0;
            maxDrugPPI(g, :)  = 0.0;
            continue;
        end

        for c = 1:numCells
            if ~strcmp( drugInfo.name, drugData{d}.name )
                fprintf('Error: %s\t%s\n', drugInfo.name, drugData{d}.name);
            end

            lmScores = drugData{d}.scores(:, c);
            pScores  = lmScores(pIdx);
            meanDrugPPI(g, c) = mean(pScores);
            maxDrugPPI(g, c)  = max( abs(pScores) );
        end

    end

    drugFeatures{d}.meanDrugPPI = meanDrugPPI;
    drugFeatures{d}.maxDrugPPI  = maxDrugPPI;

end

save('../data/drugAllFeatures', 'drugFeatures');

