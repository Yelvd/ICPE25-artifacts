
# OLD ========================

Metric_translation = {
        "Metric"                : "process",
        "Runtime (RDTSC) [s]"   : "RDTSC",
        "Runtime unhalted [s]"  : "runtime",
        "Clock [MHz]"           : "clock_speed",
        "Temperature [C]"       : "temperature",
        "Energy [J]"            : "energy",
        "Power [W]"             : "power",
        "Energy PP0 [J]"        : "energy_PP0",
        "Power PP0 [W]"         : "power_PP0",
        "Energy DRAM [J]"       : "energy_DRAM",
        "Power DRAM [W]"        : "power_DRAM",
        }

def parse_single_file(file, experiment_directory):
    """
    Parse a text file from the LIKWID expirment dir
    Searches for the "TABLE,Group" string, and parses the table underneath int and saves it in {experiment_directory}/{tool_name}
    """
    tables=[]
    lines=[]

    with open(file, "r") as f:

        for i, line in enumerate(f):
            if "TABLE,Group" in line:
                tables.append((i, line))
            lines.append(line)

        f.close()

    for l, table_info in tables:

        table_info = table_info.split(',')
        table_name = "_".join(([table_info[2]] + table_info[1].split(' ')[2:]))

        table = lines[l+1:l+int(table_info[3])+2]

        f = open(f'{experiment_directory}/{tool_name}/{table_name}.csv', "w")
        f.write(''.join(table))
        f.close()

        
def parse_dir(directory):
    """
    Parse every likwid output in an experiment dir
    """
    # iterate over files in
    # that directory
    if not os.path.exists(f'{directory}/LIKWID'):
        os.makedirs(f'{directory}/LIKWID')

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            parse_single_file(f, directory)


def load_experiment_2(experiment_directory, group, data_type=DATA_TYPE.METRIC,stat=False,clean_process_name=True):
    filename = f'{experiment_directory}/{tool_name}/{group}_{data_type.value}{"_STAT" if stat else ""}.csv'
    
    meta_info = misc.load_experiment_meta_data(experiment_directory)
    df = pd.read_csv(filename)

    drop_list = []
    for c in df.columns:
        if 'Unnamed:' in c:
            drop_list.append(c)
            
    df = df.drop(columns=drop_list)
    df['jobid'] = meta_info['General']['Id']
    df['system'] = meta_info['General']['Platform']

    if data_type == DATA_TYPE.METRIC:

        metric_name = []
        for m in df['Metric']:
            n = m
            if stat: 
                n = " ".join(m.split(" ")[:-1])

            metric_name.append(Metric_translation[n] if n in Metric_translation.keys() else n)

        df['metric_name'] = metric_name

    df['platform'] = [misc.translate_platform(p) for p in df['platform']]
    return df
