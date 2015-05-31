import json

fin = open('../data/lm_pid_to_symbols', 'r')
result = []
for line in fin:
    record = line.split()
    result.append(record[1])


with open('../data/lm_symbols.json', 'w') as outfile:
    json.dump(result, outfile, indent = 4, \
            separators = [',', ": "], sort_keys = False)

