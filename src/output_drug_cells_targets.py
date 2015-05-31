""" output drug -> cell_lines -> distil_ids """

import json

def main():
    """ main function """

    # read in common drug
    f_interest = open('../data/good_drugs_gold.json', 'r')
    cell_drug = json.load(f_interest)

    # read in drug information to get the targets
    f_drug_info = open('../data/common_drugs.json', 'r')
    drug_info = json.load(f_drug_info)

    # compound_id: [set_of_cells]
    drug_cell = {}

    cell_lines = ['MCF7', 'HA1E', 'A549', 'HT29', 'VCAP', 'A375', 'HCC515']

    for cell in cell_drug:
        if cell not in cell_lines:
            continue
        drugs = cell_drug[cell]
        for drug in drugs:
            if  drug not in drug_cell:
                drug_cell[drug] = set()
            drug_cell[drug].add(cell)


    # drugs with more than 3 cell lines
    final_drug_cell = {}
    for drug in drug_cell:
        if len(drug_cell[drug]) > 3:
            final_drug_cell[drug] = {}
            final_drug_cell[drug]['cells'] = list(drug_cell[drug])
            final_drug_cell[drug]['targets'] = drug_info[drug]['gene_symbols']

    # write out a file containing distil_ids of drugs id
    with open('../data/drugs_cell_targets.json', 'w') as outfile:
        json.dump(final_drug_cell, outfile, indent = 4, \
                separators = [',', ": "], sort_keys = False)

    print len(final_drug_cell)

main()

