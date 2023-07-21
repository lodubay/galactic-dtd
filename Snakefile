# Simulations
rule simulation_insideout_powerlaw_slope11:
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
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_powerlaw_slope11:
    output:
        directory("src/data/multizone/gaussian/lateburst/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_powerlaw_slope11:
    output:
        directory("src/data/multizone/gaussian/earlyburst/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_powerlaw_slope11:
    output:
        directory("src/data/multizone/gaussian/twoinfall/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_powerlaw_slope14:
    output:
        directory("src/data/multizone/gaussian/insideout/powerlaw_slope14")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="powerlaw",
        RIa_params={"slope": -1.4},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_powerlaw_slope14:
    output:
        directory("src/data/multizone/gaussian/lateburst/powerlaw_slope14")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="powerlaw",
        RIa_params={"slope": -1.4},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_powerlaw_slope14:
    output:
        directory("src/data/multizone/gaussian/earlyburst/powerlaw_slope14")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="powerlaw",
        RIa_params={"slope": -1.4},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_powerlaw_slope14:
    output:
        directory("src/data/multizone/gaussian/twoinfall/powerlaw_slope14")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="powerlaw",
        RIa_params={"slope": -1.4},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_exponential_timescale15:
    output:
        directory("src/data/multizone/gaussian/insideout/exponential_timescale15")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="exponential",
        RIa_params={"timescale": 1.5},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_exponential_timescale15:
    output:
        directory("src/data/multizone/gaussian/lateburst/exponential_timescale15")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="exponential",
        RIa_params={"timescale": 1.5},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_exponential_timescale15:
    output:
        directory("src/data/multizone/gaussian/earlyburst/exponential_timescale15")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="exponential",
        RIa_params={"timescale": 1.5},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_exponential_timescale15:
    output:
        directory("src/data/multizone/gaussian/twoinfall/exponential_timescale15")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="exponential",
        RIa_params={"timescale": 1.5},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_exponential_timescale30:
    output:
        directory("src/data/multizone/gaussian/insideout/exponential_timescale30")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="exponential",
        RIa_params={"timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_exponential_timescale30:
    output:
        directory("src/data/multizone/gaussian/lateburst/exponential_timescale30")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="exponential",
        RIa_params={"timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_exponential_timescale30:
    output:
        directory("src/data/multizone/gaussian/earlyburst/exponential_timescale30")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="exponential",
        RIa_params={"timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_exponential_timescale30:
    output:
        directory("src/data/multizone/gaussian/twoinfall/exponential_timescale30")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="exponential",
        RIa_params={"timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_plateau_width03:
    output:
        directory("src/data/multizone/gaussian/insideout/plateau_width03")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="plateau",
        RIa_params={"width": 0.3, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_plateau_width03:
    output:
        directory("src/data/multizone/gaussian/lateburst/plateau_width03")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="plateau",
        RIa_params={"width": 0.3, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_plateau_width03:
    output:
        directory("src/data/multizone/gaussian/earlyburst/plateau_width03")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="plateau",
        RIa_params={"width": 0.3, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_plateau_width03:
    output:
        directory("src/data/multizone/gaussian/twoinfall/plateau_width03")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="plateau",
        RIa_params={"width": 0.3, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_plateau_width10:
    output:
        directory("src/data/multizone/gaussian/insideout/plateau_width10")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="plateau",
        RIa_params={"width": 1.0, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_plateau_width10:
    output:
        directory("src/data/multizone/gaussian/lateburst/plateau_width10")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="plateau",
        RIa_params={"width": 1.0, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_plateau_width10:
    output:
        directory("src/data/multizone/gaussian/earlyburst/plateau_width10")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="plateau",
        RIa_params={"width": 1.0, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_plateau_width10:
    output:
        directory("src/data/multizone/gaussian/twoinfall/plateau_width10")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="plateau",
        RIa_params={"width": 1.0, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_prompt:
    output:
        directory("src/data/multizone/gaussian/insideout/prompt")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="prompt",
        RIa_params={"peak": 0.05, "stdev": 0.015, "timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_prompt:
    output:
        directory("src/data/multizone/gaussian/lateburst/prompt")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="prompt",
        RIa_params={"peak": 0.05, "stdev": 0.015, "timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_prompt:
    output:
        directory("src/data/multizone/gaussian/earlyburst/prompt")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="prompt",
        RIa_params={"peak": 0.05, "stdev": 0.015, "timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_prompt:
    output:
        directory("src/data/multizone/gaussian/twoinfall/prompt")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="prompt",
        RIa_params={"peak": 0.05, "stdev": 0.015, "timescale": 3.0},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_triple:
    output:
        directory("src/data/multizone/gaussian/insideout/triple")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="triple",
        RIa_params={"early_rate": 0.05, "rise_time": 0.5, "width": 0.5, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_lateburst_triple:
    output:
        directory("src/data/multizone/gaussian/lateburst/triple")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="lateburst",
        RIa="triple",
        RIa_params={"early_rate": 0.05, "rise_time": 0.5, "width": 0.5, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_triple:
    output:
        directory("src/data/multizone/gaussian/earlyburst/triple")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="triple",
        RIa_params={"early_rate": 0.05, "rise_time": 0.5, "width": 0.5, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_twoinfall_triple:
    output:
        directory("src/data/multizone/gaussian/twoinfall/triple")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="twoinfall",
        RIa="triple",
        RIa_params={"early_rate": 0.05, "rise_time": 0.5, "width": 0.5, "slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_insideout_greggio05_single:
    output:
        directory("src/data/multizone/gaussian/insideout/greggio05_single")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="greggio05_single",
        RIa_params={},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_fiducial_diffusion:
    output:
        directory("src/data/multizone/diffusion/insideout/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="diffusion",
        evolution="insideout",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_fiducial_delay150:
    output:
        directory("src/data/multizone/gaussian/insideout/powerlaw_slope11_delay150")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="insideout",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.15,
        yields="JW20"
    script:
        "src/scripts/snakesim.py"

rule simulation_earlyburst_C22yields:
    output:
        directory("src/data/multizone/gaussian/earlyburst_C22yields/powerlaw_slope11")
    cache:
        True
    conda:
        "environment.yml"
    params:
        migration="gaussian",
        evolution="earlyburst",
        RIa="powerlaw",
        RIa_params={"slope": -1.1},
        minimum_delay=0.04,
        yields="C22"
    script:
        "src/scripts/snakesim.py"

# Data
rule apogee_sample:
    output:
        "src/data/APOGEE/sample.csv"
    cache:
        True
    script:
        "src/scripts/apogee_tools.py"

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

rule midplane_distance:
    input:
        "src/data/multizone/diffusion/insideout/powerlaw_slope11",
        "src/data/multizone/gaussian/insideout/powerlaw_slope11"
    output:
        "src/tex/figures/midplane_distance.pdf"
    script:
        "src/scripts/midplane_distance.py"

rule feh_df_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/powerlaw_slope11",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        )
    output:
        "src/tex/figures/feh_df_sfh.pdf"
    script:
        "src/scripts/feh_df_sfh.py"

rule feh_df_dtd:
    input:
        "src/data/multizone/gaussian/insideout/powerlaw_slope14"
        "src/data/multizone/gaussian/insideout/exponential_timescale30"
    output:
        "src/tex/figures/feh_df_dtd.pdf"
    script:
        "src/scripts/feh_df_dtd.py"

rule ofe_df_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/powerlaw_slope11",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        )
    output:
        "src/tex/figures/ofe_df_sfh.pdf"
    script:
        "src/scripts/ofe_df_sfh.py"

rule ofe_df_dtd:
    input:
        expand("src/data/multizone/gaussian/insideout/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        )
    output:
        "src/tex/figures/ofe_df_dtd.pdf"
    script:
        "src/scripts/ofe_df_dtd.py"

rule ofe_feh_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/powerlaw_slope11",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        )
    output:
        "src/tex/figures/ofe_feh_sfh.pdf"
    script:
        "src/scripts/ofe_feh_sfh.py"

rule ofe_feh_dtd:
    input:
        expand("src/data/multizone/gaussian/insideout/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        )
    output:
        "src/tex/figures/ofe_feh_dtd.pdf"
    script:
        "src/scripts/ofe_feh_dtd.py"

rule age_ofe_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/powerlaw_slope11",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        )
    output:
        "src/tex/figures/age_ofe_sfh.pdf"
    script:
        "src/scripts/age_ofe_sfh.py"

rule age_ofe_dtd:
    input:
        expand("src/data/multizone/gaussian/insideout/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        )
    output:
        "src/tex/figures/age_ofe_dtd.pdf"
    script:
        "src/scripts/age_ofe_dtd.py"

# Tables
rule summary_table:
    input:
        expand("src/data/multizone/gaussian/{evolution}/{dtd}",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"],
               dtd=["powerlaw_slope11", "powerlaw_slope14", 
                    "exponential_timescale15", "exponential_timescale30", 
                    "plateau_width03", "plateau_width10", "prompt", "triple"
               ]
        ),
        "src/data/APOGEE/sample.csv"
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
