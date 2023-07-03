rule simulation_IOplaw1:
    output:
        directory("src/data/multizone/gaussian/insideout/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="powerlaw",
        RIa_params={slope=-1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"
        
rule simulations:
    output:
        directory("src/data/multizone")
    cache:
        True
    conda:
        "environment.yml"
    script:
        "src/scripts/simulations.py"

rule star_formation_histories:
    input:
        "src/data/multizone"
    output:
        "src/figures/star_formation_histories.pdf"
    script:
        "src/scripts/star_formation_histories.py"

rule sample_size:
    input:
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/sample_size.txt"
    script:
        "src/scripts/sample_size.py"
