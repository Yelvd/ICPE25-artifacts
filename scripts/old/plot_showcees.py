import sys
import glob
import os
sys.path.append('./scripts/tools')
import likwid
import misc
from enum import Enum

import pandas as pd

from bs4 import BeautifulSoup


# relative = 0
# fixed = 1


results = "./results/showcees/"

experiments = []


for filename in os.listdir(results):
    f = os.path.join(results, filename)
    # checking if it is a file
    if os.path.isdir(f):
        print(f'{f}/meta.yml')
        meta = misc.load_experiment_meta_data(f)

        if meta['Instrumentation']['Collection tool'] == "likwid":

            filename = glob.glob(f'{f}/config*.xml')[0]
            # Reading the data inside the xml file to a variable under the name  data
            with open(filename, 'r') as file:
                data = file.read()

            # Passing the stored data inside the beautifulsoup parser
            bs_data = BeautifulSoup(data, 'xml')

            # Using find() to extract attributes of the first instance of the tag
            b_name = bs_data.find('nx')
            for test in b_name.children:
                nx = int(test)

            b_name = bs_data.find('tmax')
            for test in b_name.children:
                iters = int(test)

            experiments.append((f, meta, 0 if nx == 125 else 1, iters))


data_full = {"PWR": [], "PWR DRAM": [], "Energy": [], "Energy DRAM": [], "L2_miss_ratio": [], "L3_miss_ratio": [], "FLOPS_DP": []}
for exp in experiments:

    # fixed size, full run
    if exp[2] == 1 and exp[3] > 0:
        df = likwid.load_experiment(exp[0], "ENERGY", stat=True)
        
        data_full["Energy"].append(df['Sum'][5])
        data_full["Energy DRAM"].append(df['Sum'][9])
        data_full["PWR"].append(df['Sum'][6])
        data_full["PWR DRAM"].append(df['Sum'][10])

        df = likwid.load_experiment(exp[0], "L2CACHE", stat=True)
        data_full["L2_miss_ratio"].append(df['Avg'][6])

        df = likwid.load_experiment(exp[0], "L3CACHE", stat=True)
        data_full["L3_miss_ratio"].append(df['Avg'][6])

        df = likwid.load_experiment(exp[0], "FLOPS_DP", stat=True)
        data_full["FLOPS_DP"].append(df['Sum'][4])


print(pd.DataFrame.from_dict(data_full).style.format(precision=2).to_latex())


data_full = {"PWR": [], "PWR DRAM": [], "Energy": [], "Energy DRAM": [], "L2_miss_ratio": [], "L3_miss_ratio": [], "FLOPS_DP": []}
for exp in experiments:

    # relative size, full run
    if exp[2] == 0 and exp[3] > 0:
        df = likwid.load_experiment(exp[0], "ENERGY", stat=True)
        data_full["Energy"].append(df['Sum'][5])
        data_full["Energy DRAM"].append(df['Sum'][9])
        data_full["PWR"].append(df['Sum'][6])
        data_full["PWR DRAM"].append(df['Sum'][10])

        df = likwid.load_experiment(exp[0], "L2CACHE", stat=True)
        data_full["L2_miss_ratio"].append(df['Avg'][6])

        df = likwid.load_experiment(exp[0], "L3CACHE", stat=True)
        data_full["L3_miss_ratio"].append(df['Avg'][6])

        df = likwid.load_experiment(exp[0], "FLOPS_DP", stat=True)
        data_full["FLOPS_DP"].append(df['Sum'][4])


print(pd.DataFrame.from_dict(data_full).style.format(precision=2).to_latex())
