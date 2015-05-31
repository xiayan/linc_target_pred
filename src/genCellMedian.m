% get landmark gene symbols
filepath = '/Volumes/Seagate/LINCS_L4/zspc_n1328098x22268.gctx';
meta = parse_gctx(filepath, 'annot_only', 1);
lm = meta.rid(1:978);

load('../data/cells_GS_distils.mat');
cellsDistils = a;

cellsNames = fieldnames(cellsDistils);
numCells  = length(cellsNames);
% numCells = 1;

cellData = cell(numCells, 1);

for c = 1:numCells
    cellName = cellsNames{c};
    cellData{c}.name  = cellName;

    genesDistils = cellsDistils.(cellName);
    genesNames   = fieldnames(genesDistils);

    fprintf('\nProcessing %s, %d\n', cellName, c);

    numGenes = length(genesNames);

    distilIDs = [];

    for g = 1:numGenes
        geneName = genesNames{g};
        ids = genesDistils.(geneName);
        distilIDs = [distilIDs, ids];
    end

    fprintf('Retriving %d experiments\n', length(distilIDs));
    cellScores = parse_gctx(filepath, 'cid', distilIDs, 'rid', lm);
    cellScores = cellScores.mat;

    scores = zeros(978, numGenes);
    offset = 1;
    for g = 1:numGenes
        geneName = genesNames{g};
        numExp = length( genesDistils.(geneName) );
        geneScores = cellScores(:, offset:offset+numExp-1);
        scores(:, g) = median(geneScores, 2);
        offset = offset + numExp;
    end

    if size(cellScores, 2) + 1 ~= offset
        fprintf('%d\t%d', size(cellScores, 2), offset);
    end

    ctlIdx = find( strcmp('ctl', genesNames) );
    genesNames(ctlIdx, :) = [];
    ctlVec = scores(:, ctlIdx);
    scores(:, ctlIdx) = [];

    cellData{c}.genes  = genesNames;
    cellData{c}.scores = scores;
    cellData{c}.ctl    = ctlVec;

end

save('../data/cellData', 'cellData');


