import sys
import os
sys.path.append('./scripts/tools')
sys.path.append('./scripts/')
import likwid
import misc

results = "./results/snellius/"

for filename in os.listdir(results):
    f = os.path.join(results, filename)
    # checking if it is a file
    if os.path.isdir(f):
        meta = misc.load_experiment_meta_data(f)
        if meta is None:
            continue
        print(f'{f}/meta.yaml')

        if meta.get("Instrumentation"):
            if meta['Instrumentation']['Collection tool'] == "likwid":
                likwid.parse_dir(f)

        if meta['General'].get("Monitor-tool"):
            if meta['General']['Monitor-tool'] == "likwid":
                likwid.parse_dir(f)
