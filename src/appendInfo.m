function drug152Results = appendInfo(molResult)
    numDrugs = numel(molResult);
    drug152Results = cell(0);

    for d = 1:numDrugs
        fprintf('Processing drug %d\n', d);
        curDrug   = molResult{d};
        drugName = curDrug.name;
        drugName = strrep(drugName, '_0x2D_', '-');
        url = strcat('http://api.lincscloud.org/a2/pertinfo?q={%22pert_id%22:%22', drugName, '%22}&user_key=c888934e11a652fcc7a23eb059244692');
        curRecord = loadjson( urlread(url) );
        curRecord = curRecord{1};
        curRecord.geneRanking = curDrug.geneRanking;
        curRecord.targets = curDrug.targets;
        curRecord.minTargetRank = curDrug.minTargetRank;

        drug152Results = [drug152Results; curRecord];
    end


