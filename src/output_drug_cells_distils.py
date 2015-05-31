""" output drugs -> cell_id -> distil_ids """

import json

def main():
    """ main function """
    # load drug_cells_targets
    f_drugs_info = open('../data/drugs_cells_targets.json', 'r')
    drugs_info = json.load(f_drugs_info)

    # read in inst.info
    meta_info = open('../data/inst.info', 'r')
    # get rid of header line
    meta_info.readline()

    # read in golden
    golden = json.load( open('../data/golden_distil.json', 'r') )
    golden_set = set(golden)

    drugs_cells_distil = {}
    for drug in drugs_info:
        drugs_cells_distil[drug] = {}
        for cell in drugs_info[drug]['cells']:
            drugs_cells_distil[drug][cell] = list()


    for record in meta_info:
        info = record.split('\t')
        cell_id = info[3]
        pert_id = info[1]
        pert_type = info[-10]
        distil_id = info[0]

        if (pert_type == 'trt_cp') and (distil_id in golden_set) \
                and (pert_id in drugs_cells_distil) \
                and (cell_id in drugs_cells_distil[pert_id]):
            drugs_cells_distil[pert_id][cell_id].append(distil_id)

    with open('../data/drugs_cells_distil.json', 'w') as outfile:
        json.dump(drugs_cells_distil, outfile, indent = 4, \
                separators = [',', ": "], sort_keys = False)

main()

