cell_lines = {'MCF7', 'HA1E', 'A549', 'HT29', 'VCAP', 'A375'};

% get landmark gene symbols
filepath = '/Volumes/Seagate/LINCS_L4/zspc_n1328098x22268.gctx';
meta = parse_gctx(filepath, 'annot_only', 1);
lm = meta.rid(1:978);

% get gene symbols for drug targets
target_ids = loadjson('common_drug_targets.json');
drug_names = fieldnames(target_ids);

num_drugs = 29;
num_cells = length(cell_lines);

% this matrix records the ranking of all 3104 common genes for each drug
% across 6 cell lines
drug_medians = zeros(978, num_drugs, num_cells);

for c = 1:num_cells
    cell_name = cell_lines{c};
    fprintf('\nProcessing %s\n', cell_name);
    cmd_str = strjoin({'python query_meta.py', cell_name}, ' ');
    fprintf('Executing %s\n', cmd_str);
    system(cmd_str);

    disp('Loading json files');
    all_drugs = loadjson('all_drugs_distil.json');
    drugs_distil = loadjson('drugs_distil_id.json');

    % get all drugs z scores
    all_drug_scores = parse_gctx(filepath, 'cid', all_drugs, 'rid', lm);
    all_drug_scores = all_drug_scores.mat;

    drug_ids = fieldnames(drugs_distil);
    num_drugs = length(drug_ids);

    % get the medians of scores of drugs and knockdown experiments
    fprintf('Taking drug medians\n');
    counter = 1;
    for i = 1:num_drugs
        num_exp = length(drugs_distil.(drug_ids{i}));
        all_scores = all_drug_scores(:, counter:counter+num_exp-1);
        di = find( strcmp(drug_names, drug_ids{i}) );
        drug_medians(:, di, c) = median(all_scores, 2);
        counter = counter + num_exp;
    end

    if counter ~= size(all_drug_scores, 2) + 1
        fprintf('Drug error: counter %d, size %d', counter, ...
            size(all_drug_scores, 2));
    end
end

save('drug_medians', 'drug_medians');

