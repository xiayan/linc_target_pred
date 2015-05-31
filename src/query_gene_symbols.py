""" Read in procession ids and convert them to gene symbols """
import json, urllib2

def main():
    """ main function """
    fin = open('./pr_ids.json', 'r')
    jin = json.load(fin)

    pr_ids = jin[0]
    symbols = [''] * len(pr_ids)
    counter = 0


    ids = pr_ids[0:250]
    for pr_id in ids:
        target_url = 'http://api.lincscloud.org/a2/geneinfo?q={%22pr_id%22:%22' + pr_id + '%22}&f={%22pr_gene_symbol%22:1,%22is_l1000%22:1}&user_key=c888934e11a652fcc7a23eb059244692'
        try:
            result = json.load(urllib2.urlopen(target_url))
            symbols[counter] = result[0]['pr_gene_symbol']
            print pr_id + "\t" + symbols[counter]
            counter = counter + 1
            if not result[0]['is_l1000']:
                print "Error: " + result[0]['pr_gene_symbol']
        except urllib2.URLError:
            continue

main()

