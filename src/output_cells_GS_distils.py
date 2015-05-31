""" output distil_number of approved drugs and
knockdowns for each cell type """

import json
import random

def main():
    """ main function """
    cell_lines = \
        set(['MCF7', 'HA1E', 'A549', 'HT29', 'VCAP', 'A375', 'HCC515'])

    # read in inst.info
    meta_info = open('../data/inst.info', 'r')
    # get rid of header line
    meta_info.readline()

    # read in golden
    golden = json.load( open('../data/golden_distil.json', 'r') )
    golden_set = set(golden)

    # cell -> KD gene symbol -> distil ids
    cell_gs_distil = {}

    for record in meta_info:
        info = record.split('\t')
        cell_id = info[3]
        pert_type = info[-10]
        pert_desc = info[2]
        distil_id = info[0]

        if cell_id not in cell_lines:
            continue


        if (pert_type == 'trt_sh') and (distil_id in golden_set):
            if cell_id not in cell_gs_distil:
                cell_gs_distil[cell_id] = {}
            if pert_desc not in cell_gs_distil[cell_id]:
                cell_gs_distil[cell_id][pert_desc] = list()
            cell_gs_distil[cell_id][pert_desc].append(distil_id)
        elif pert_type[:3] == 'ctl':
            if cell_id not in cell_gs_distil:
                cell_gs_distil[cell_id] = {}
            if 'ctl' not in cell_gs_distil[cell_id]:
                cell_gs_distil[cell_id]['ctl'] = list()
            cell_gs_distil[cell_id]['ctl'].append(distil_id)


    # random sample 100 untreated experiment
    for cell in cell_gs_distil:
        all_ctl = cell_gs_distil[cell]['ctl']

        smpl_ctl = [ all_ctl[i] \
                for i in sorted(random.sample(xrange(len(all_ctl)), 50))  ]

        cell_gs_distil[cell]['ctl'] = smpl_ctl

    # write out a file containing distil_ids of drugs id
    with open('../data/cells_GS_distils.json', 'w') as outfile:
        json.dump(cell_gs_distil, outfile, indent = 4, \
                separators = [',', ": "], sort_keys = False)

main()


