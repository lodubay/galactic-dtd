# Simulations
rule simulation_powerlaw_slope11:
    output:
        directory("src/data/multizone/{migration}/{evolution}/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration={migration},
        evolution={evolution},
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

# Figures
rule star_formation_histories:
    input:
        "src/data/multizone/gaussian/insideout/powerlaw_slope11",
        "src/data/multizone/gaussian/lateburst/powerlaw_slope11",
        "src/data/multizone/gaussian/earlyburst/powerlaw_slope11",
        "src/data/multizone/gaussian/twoinfall/powerlaw_slope11"
    output:
        "src/figures/star_formation_histories.pdf"
    script:
        "src/scripts/star_formation_histories.py"

# Variables
rule sample_size:
    input:
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/sample_size.txt"
    script:
        "src/scripts/sample_size.py"
