import os
import subprocess

def run_shell_cmd(command, cwd='./'):
    return subprocess.run(command.split(' '), stdout=subprocess.PIPE, cwd=cwd).stdout.decode('utf-8')


hemocell_dir = "/home/jelle/HemoCell-dev/"
experiment_dir = "/home/jelle/HemoCell-dev/Hemocell-Performance-Benchmarks/cube-benchmark/"
benchmarks = ["cube-fixed", "cube-fixed-large", "cube-relative"]

print("Hallo")
platform = "Ampere_Q8030" 
base_directory = "./results/eviden-perf/Ampere-aarch64-prof/" 
for b in benchmarks:
    print(f"{platform} {b}")
    d = f'{base_directory}/{b}'
    experiments = [ (f.path, f.name) for f in os.scandir(d) if f.is_dir() ]

    for e_dir, info in experiments:
        info = info.split('_')
        tool = "perf"
        tasks = info[3]
        nodes = 1
        exp_id = info[-1]

        args = f"-p {platform} -m {tool} -i {exp_id} -t {tasks} -n {nodes} {hemocell_dir} {experiment_dir} -b {b} -f {e_dir}/meta.yaml"

        run_shell_cmd(f"python3 ./scripts/old/generate-meta-data-yaml.py {args}")


platform = "AMD_Epyc" 
base_directory = "./results/eviden-perf/Amd-epyc-prof/" 
for b in benchmarks:
    print(f"{platform} {b}")
    d = f'{base_directory}/{b}'
    experiments = [ (f.path, f.name) for f in os.scandir(d) if f.is_dir() ]

    for e_dir, info in experiments:
        info = info.split('_')

        tool = "perf"
        tasks = info[3]
        nodes = 1
        exp_id = info[-1]

        args = f"-p {platform} -m {tool} -i {exp_id} -t {tasks} -n {nodes} {hemocell_dir} {experiment_dir} -b {b} -f {e_dir}/meta.yaml"

        run_shell_cmd(f"python3 ./scripts/old/generate-meta-data-yaml.py {args}")

platform = "Intel_SPR_DDR"
base_directory = "./results/eviden-perf/Intel-sapphire-rapids/ddr-prof/" 
for b in benchmarks:
    print(f"{platform} {b}")
    d = f'{base_directory}/{b}'
    experiments = [ (f.path, f.name) for f in os.scandir(d) if f.is_dir() ]

    for e_dir, info in experiments:
        info = info.split('_')

        tool = "perf"
        tasks = info[3]
        nodes = 1
        exp_id = info[-1]

        args = f"-p {platform} -m {tool} -i {exp_id} -t {tasks} -n {nodes} {hemocell_dir} {experiment_dir} -b {b} -f {e_dir}/meta.yaml"

        run_shell_cmd(f"python3 ./scripts/old/generate-meta-data-yaml.py {args}")

platform = "Intel_SPR_HBM"
base_directory = "./results/eviden-perf/Intel-sapphire-rapids/hbm-prof/" 
for b in benchmarks:
    print(f"{platform} {b}")
    d = f'{base_directory}/{b}'
    experiments = [ (f.path, f.name) for f in os.scandir(d) if f.is_dir() ]

    for e_dir, info in experiments:
        info = info.split('_')

        tool = "perf"
        tasks = info[3]
        nodes = 1
        exp_id = info[-1]

        args = f"-p {platform} -m {tool} -i {exp_id} -t {tasks} -n {nodes} {hemocell_dir} {experiment_dir} -b {b} -f {e_dir}/meta.yaml"

        run_shell_cmd(f"python3 ./scripts/old/generate-meta-data-yaml.py {args}")
