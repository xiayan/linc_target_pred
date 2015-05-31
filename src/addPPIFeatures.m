load('../data/drugCorr.mat');
load('../data/ppi.mat');

proteins = fieldnames(ppi);
num_drugs = length(drugCorr);

for d = 1:num_drugs
    drugInfo = drugCorr{d};
    ppi_scores = zeros(size(drugInfo.correlations));

    fprintf('Processing %s, %d\n', drugInfo.name, d);

    [~, indices] = sort(drugInfo.correlations, 1, 'descend');

    drugCells = drugInfo.cells;
    numCells  = length(drugCells);
    numGenes  = size(drugInfo.correlations, 1);

    for c = 1:numCells
        idx = indices(:, c);
        gene_ids = drugInfo.genes;
        sorted_gene = gene_ids(idx);
        top200 = sorted_gene(1:200);

        for g = 1:numGenes
            if sum( strcmp(proteins, gene_ids{g}) ) == 0
                ppi_scores(g, c) = 0;
                continue;
            end
            partners = ppi.(gene_ids{g});
            c_part = intersect(gene_ids, partners);
            inter  = intersect(top200, partners);
            ppi_scores(g, c) = length(inter) / (length(c_part) + 50);
        end
    end

    drugCorr{d}.ppi_scores = ppi_scores;

end

save('../data/drugCorrPPI', 'drugCorr');

