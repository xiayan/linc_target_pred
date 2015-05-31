Predicting drug targets from LINCS dataset
===

* `output_drug_cells_targets.py`: Output the cells and gene targets for each drugs in a file named `drugs_cells_targets.json`
* `output_drug_cells_distil.py`: For each drug, output the distil_ID for its cell lines.in a file named `drugs_cells_distils.json`
* `output_cells_GS_distils.json`: For each cell line, output the distil_ID for each knockdown gene symbol.
* `genDrugMedian.m`: For each drug, extract the distil_ids from LINCS, compute the median of 978 landmark genes. Results are saved in `drugData.mat`.
* `genCellMedian.m`: For each cell, extract the distil_ids for each knockdown gene symbols, and compute the median of 978 landmark genes. Results are saved in `cellData.mat`.
* `genCorrFeatures.m`: Read in `drugData.mat` and `cellData.mat`, genearate correlation features. Results are saved in `drugCorr.mat`.
* `addPPIFeatures.m`: Read in `drugCorr.mat` and `ppi.mat` (an adjancency list of PPI network), append the PPI feature. Results saved in `drugCorrPPI.mat`.
* `addDrugPPIFeature.m`: Read in `drugCorrPPI.mat`, `drugData.mat` and `ppi.mat`, append the drug response PPI feature. Results saved in `drugAllFeatures.mat`.
* `adHocTreeBagger.m`: Classification using ad hoc random forests.
* `unifiedTreeBagger.m`: Classification using unified random forests. Results saved in `drugRankingResults.mat`.