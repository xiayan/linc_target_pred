% get landmark gene symbols
filepath = '/Volumes/Seagate/LINCS_L4/zspc_n1328098x22268.gctx';
meta = parse_gctx(filepath, 'annot_only', 1);
lm = meta.rid(1:978);

load('../data/drugs_cells_targets.mat');
drugsTargets = a;

load('../data/drugs_cells_distil.mat');
drugsDistils = a;

drugsNames = fieldnames(drugsTargets);
numDrugs   = length(drugsNames);

drugData = cell(numDrugs, 1);

for d = 1:numDrugs
    drugName = drugsNames{d};
    drugData{d}.name = drugName;
    drugData{d}.targets = drugsTargets.(drugName).('targets');

    fprintf('\nProcessing %s, %d\n', drugName, d);
    cellDistil = drugsDistils.(drugName);

    cellLines = fieldnames(cellDistil);
    drugData{d}.cells = cellLines;

    scores = zeros(978, length(cellLines));

    distilIDs = [];
    for c = 1:length(cellLines)
        cellName = cellLines{c};
        ids = cellDistil.(cellName);
        distilIDs = [distilIDs, ids];
    end

    fprintf('Retriving %d experiments', length(distilIDs));
    drugsScores = parse_gctx(filepath, 'cid', distilIDs, 'rid', lm);
    drugsScores = drugsScores.mat;

    offset = 1;
    for c = 1:length(cellLines)
        cellName = cellLines{c};
        numExp = length( cellDistil.(cellName) );
        localData = drugsScores(:, offset:offset+numExp-1);

        scores(:, c) = median(localData, 2);
        offset = offset + numExp;
    end

    if size(drugsScores, 2) + 1 ~= offset
        fprintf('Error: %d\t%d', size(drugsScores, 2), offset);
    end

    drugData{d}.scores = scores;

end

save('../data/drugData', 'drugData');

