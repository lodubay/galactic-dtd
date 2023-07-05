# Simulations
rule simulation_powerlaw_slope11:
    output:
        directory("src/data/multizone/{migration}/{evolution}/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_powerlaw_slope14:
    output:
        directory("src/data/multizone/{migration}/{evolution}/powerlaw_slope14")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="powerlaw",
        RIa_params={"slope": -1.4},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_exponential_timescale15:
    output:
        directory("src/data/multizone/{migration}/{evolution}/exponential_timescale15")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="exponential",
        RIa_params={"timescale": 1.5},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_exponential_timescale30:
    output:
        directory("src/data/multizone/{migration}/{evolution}/exponential_timescale30")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="exponential",
        RIa_params={"timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_plateau_width03:
    output:
        directory("src/data/multizone/{migration}/{evolution}/plateau_width03")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="plateau",
        RIa_params={"width": 0.3, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_plateau_width10:
    output:
        directory("src/data/multizone/{migration}/{evolution}/plateau_width10")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="plateau",
        RIa_params={"width": 1.0, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_prompt:
    output:
        directory("src/data/multizone/{migration}/{evolution}/prompt")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="prompt",
        RIa_params={"peak": 0.05, "stdev": 0.015, "timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_triple:
    output:
        directory("src/data/multizone/{migration}/{evolution}/triple")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="{migration}",
        evolution="{evolution}",
        RIa="triple",
        RIa_params={"early_rate": 0.05, "rise_time": 0.5, "width": 0.5, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

# Figures
rule star_formation_histories:
    input:
        expand("src/data/multizone/gaussian/{evolution}/powerlaw_slope11",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        )
    output:
        "src/tex/figures/star_formation_histories.pdf"
    script:
        "src/scripts/star_formation_histories.py"

rule radial_migration:
    input:
        "src/data/multizone/diffusion/insideout/powerlaw_slope11",
        "src/data/multizone/gaussian/insideout/powerlaw_slope11"
    output:
        "src/tex/figures/radial_migration.pdf"
    script:
        "src/scripts/radial_migration.py"

# Tables
rule summary_table:
    input:
        expand("src/data/multizone/gaussian/{evolution}/{dtd}",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"],
               dtd=["powerlaw_slope11", "powerlaw_slope14", 
                    "exponential_timescale15", "exponential_timescale30", 
                    "plateau_width03", "plateau_width10", "prompt", "triple"
               ]
        )
    output:
        multiext("src/tex/output/summary_table", ".csv", ".tex")
    script:
        "src/scripts/summary_table.py"

# Variables
rule sample_size:
    input:
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/sample_size.txt"
    script:
        "src/scripts/sample_size.py"
