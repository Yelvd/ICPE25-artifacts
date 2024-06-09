from os import listdir
from os.path import isfile, join
from scripts.tools import misc

import sys

figsDir = "artifacts"
scriptDir = "scripts"
matplotlibStyle_ = "matplotlib-style.rc"
resultsDir_ = "results/"

PLATFORMS = ['AMD Bergamo', 'Ampere Altra', 'Intel SPR DDR', 'Intel SPR HBM', "AMD Rome 7742", "AMD Genoa", "AMD Rome 7H12"]
EXPERIMENTS = ["S1", "S2", "W1"]
METRICS = ["time", "power", "energy"]

MASCOTS24_FIGURES = ["S1-time-all.pdf", "S1-energy-all.pdf", "S1-power-all.pdf", "S1-time-AMD_Bergamo.pdf"]
MASCOTS24_FIGURES = MASCOTS24_FIGURES + [f"S1-time-{p.replace(' ','_')}.pdf" for p in PLATFORMS]
MASCOTS24_FIGURES = MASCOTS24_FIGURES + [f"S1-energy-{p.replace(' ','_')}.pdf" for p in PLATFORMS]
MASCOTS24_TABLES  = ["S1-CPI-all.tex", "S1-cache-all.tex", "S1-mem-all.tex"]

MASCOTS24_ARTIFECTS = MASCOTS24_FIGURES + MASCOTS24_TABLES
 
rule all:
    input:
        [f"{figsDir}/{e}-{m}-{p.replace(' ', '_')}.pdf"  for p in PLATFORMS for e in EXPERIMENTS for m in METRICS],
        [f"{figsDir}/{e}-{m}-all.pdf" for e in EXPERIMENTS for m in METRICS],


rule MASCOTS24:
    input:
        [f"{figsDir}/{artifect}" for artifect in MASCOTS24_ARTIFECTS]


rule plot_roofline:
    input:
        "scripts/plot-roofline.py"
    params:
        resultsDir = "./backup/results-components/",
        matplotlibStyle = matplotlibStyle_,
        platform = "AMD Genoa",
    output:
        S1_pdf = f"{figsDir}/roofline_S1_AMD_Genoa.pdf",
        S1_png = f"{figsDir}/roofline_S1_AMD_Genoa.png",
        S2_pdf = f"{figsDir}/roofline_S2_AMD_Genoa.pdf",
        S2_png = f"{figsDir}/roofline_S2_AMD_Genoa.png",
    script:
        "scripts/plot-roofline.py"


for platform in PLATFORMS:
    for metric in METRICS:
        rule:
            name: f"plot-all-{metric}-{platform.replace(' ','_')}"
            input:
                "scripts/plot-single-system-multi-bench.py"
            params:
                resultsDir = resultsDir_,
                matplotlibStyle = matplotlibStyle_,
                experiment = "all",
                platform = platform,
                metric = metric
            output:
                pdf = f"{figsDir}/all-{metric}-{platform.replace(' ','_')}.pdf",
                png = f"{figsDir}/all-{metric}-{platform.replace(' ','_')}.png"
            script:
                "scripts/plot-single-system-multi-bench.py"

        for experiment in EXPERIMENTS:
            rule:
                name: f"plot-{experiment}-{metric}-{platform.replace(' ','_')}"
                input:
                    "scripts/plot-single-system.py"
                params:
                    resultsDir = resultsDir_,
                    matplotlibStyle = matplotlibStyle_,
                    experiment = experiment,
                    platform = platform,
                    metric = metric
                output:
                    pdf = f"{figsDir}/{experiment}-{metric}-{platform.replace(' ','_')}.pdf",
                    png = f"{figsDir}/{experiment}-{metric}-{platform.replace(' ','_')}.png"
                script:
                    "scripts/plot-single-system.py"

    

for experiment in EXPERIMENTS:
    for metric in METRICS:
        rule:
            name: f"plot-{experiment}-{metric}-all"
            input:
                "scripts/plot-all-systems.py"
            params:
                resultsDir = resultsDir_,
                matplotlibStyle = matplotlibStyle_,
                experiment = experiment,
                metric = metric
            output:
                pdf = f"{figsDir}/{experiment}-{metric}-all.pdf",
                png = f"{figsDir}/{experiment}-{metric}-all.png"
            script:
                "scripts/plot-all-systems.py"

for experiment in EXPERIMENTS + ['all']:
    for platform in PLATFORMS + ['all']:
        # Instruction comparison
        rule:
            name: f"table-{experiment}-instr-comparison"
            input:
                "scripts/table-instr-comparison.py"
            params:
                resultsDir = resultsDir_,
                experiment = experiment,
            output:
                tex = f"{figsDir}/{experiment}-CPI-all.tex",
            script:
                "scripts/table-instr-comparison.py"

        rule:
            name: f"table-{experiment}-instr-component"
            input:
                "scripts/table-instr-component.py"
            params:
                resultsDir = "./results-components/",
                experiment = experiment,
            output:
                tex = f"{figsDir}/{experiment}-CPI-all-component.tex",
            script:
                "scripts/table-instr-component.py"

        rule:
            name: f"table-{experiment}-L1-component"
            input:
                "scripts/table-L1-component.py"
            params:
                resultsDir = "./results-components/",
                experiment = experiment,
            output:
                tex = f"{figsDir}/{experiment}-L1-all-component.tex",
            script:
                "scripts/table-L1-component.py"

        # Cache comparision
        rule:
            name: f"table-{experiment}-cache-{platform.replace(' ', '_')}"
            input:
                "scripts/table-cache-comparison.py"
            params:
                resultsDir = resultsDir_,
                experiment = experiment,
                platform = platform
            output:
                tex = f"{figsDir}/{experiment}-cache-{platform.replace(' ', '_')}.tex",
            script:
                "scripts/table-cache-comparison.py"

        # Cache comparision
        rule:
            name: f"table-{experiment}-cacheL3-{platform.replace(' ', '_')}"
            input:
                "scripts/table-cacheL3-comparison.py"
            params:
                resultsDir = resultsDir_,
                experiment = experiment,
                platform = platform
            output:
                tex = f"{figsDir}/{experiment}-cacheL3-{platform.replace(' ', '_')}.tex",
            script:
                "scripts/table-cacheL3-comparison.py"

        # Memory comparision
        rule:
            name: f"table-{experiment}-mem-{platform.replace(' ', '_')}"
            input:
                "scripts/table-cache-comparison.py"
            params:
                resultsDir = resultsDir_,
                experiment = experiment,
                platform = platform
            output:
                tex = f"{figsDir}/{experiment}-mem-{platform.replace(' ', '_')}.tex",
            script:
                "scripts/table-mem-comparison.py"


rule plot:
    input:
        "scripts/plot-{figname}.py"
    params:
        resultsDir = resultsDir_,
        matplotlibStyle = matplotlibStyle_
    output:
        pdf = "figs/{figname}.pdf",
        png = "figs/{figname}.png",
        tex = "figs/{figname}.tex"
    script:
        "scripts/plot-{wildcards.figname}.py"
