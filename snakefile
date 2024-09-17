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

ICPE25_FIGURES = ["S1-all.pdf", "S1-time-all.pdf", "S1-energy-all.pdf", "roofline_S1_AMD_Genoa.pdf", "roofline_S2_AMD_Genoa.pdf", "all-energy-AMD_Genoa.pdf", "all-time-AMD_Genoa.pdf", "all-power-AMD_Genoa.pdf"]
ICPE25_TABLES  = ["S1-CPI-all.tex", "S1-cache-all.tex", "S1-mem-all.tex"]

ICPE25_ARTIFECTS = ICPE25_FIGURES + ICPE25_TABLES
 
rule ICPE25:
    input:
        [f"{figsDir}/{artifect}" for artifect in ICPE25_ARTIFECTS]


rule plot_roofline:
    input:
        "scripts/plot-roofline.py"
    params:
        resultsDir = f"{resultsDir_}",
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
            name: f"plot-{experiment}-all"
            input:
                "scripts/plot-all-systems.py"
            params:
                resultsDir = resultsDir_,
                matplotlibStyle = matplotlibStyle_,
                experiment = experiment,
            output:
                pdf = f"{figsDir}/{experiment}-all.pdf",
                png = f"{figsDir}/{experiment}-all.png"
            script:
                "scripts/plot-all-systems.py"

        for experiment in EXPERIMENTS:
            rule:
                name: f"plot-{experiment}-{metric}-all"
                input:
                    "scripts/plot-single-system.py"
                params:
                    resultsDir = resultsDir_,
                    matplotlibStyle = matplotlibStyle_,
                    experiment = experiment,
                    platform = "all",
                    metric = metric
                output:
                    pdf = f"{figsDir}/{experiment}-{metric}-all.pdf",
                    png = f"{figsDir}/{experiment}-{metric}-all.png"
                script:
                    "scripts/plot-single-system.py"

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
                "scripts/table-mem-comparison.py"
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
