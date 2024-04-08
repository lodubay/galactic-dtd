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
        ),
        "src/data/APOGEE/sample.csv"
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
        
rule feh_df_comparison:
    input:
        "src/data/multizone/gaussian/insideout/exponential_timescale15"
        "src/data/multizone/gaussian/twoinfall/exponential_timescale15"
        "src/data/multizone/gaussian/insideout/powerlaw_slope14"
        "src/data/multizone/gaussian/insideout/exponential_timescale30"
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/feh_df_comparison.pdf"
    script:
        "src/scripts/feh_df_comparison.py"

rule ofe_df_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/exponential_timescale15",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_df_sfh.pdf"
    script:
        "src/scripts/ofe_df_sfh.py"

rule ofe_df_dtd:
    input:
        expand("src/data/multizone/gaussian/earlyburst/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_df_dtd.pdf"
    script:
        "src/scripts/ofe_df_dtd.py"

rule ofe_bimodality_summary:
    input:
        expand("src/data/multizone/gaussian/{evolution}/exponential_timescale15",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        ),
        expand("src/data/multizone/gaussian/lateburst/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_bimodality_summary.pdf"
    script:
        "src/scripts/ofe_bimodality_summary.py"

rule ofe_feh_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/exponential_timescale15",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_sfh.pdf"
    script:
        "src/scripts/ofe_feh_sfh.py"

rule ofe_feh_dtd:
    input:
        expand("src/data/multizone/gaussian/insideout/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_dtd.pdf"
    script:
        "src/scripts/ofe_feh_dtd.py"

rule ofe_feh_twoinfall:
    input:
        "src/data/multizone/gaussian/twoinfall/plateau_width10",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_twoinfall.pdf"
    script:
        "src/scripts/ofe_feh_twoinfall.py"

rule age_ofe_sfh:
    input:
        expand("src/data/multizone/gaussian/{evolution}/exponential_timescale15",
               evolution=["insideout", "lateburst", "earlyburst", "twoinfall"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/age_ofe_sfh.pdf"
    script:
        "src/scripts/age_ofe_sfh.py"

rule age_ofe_dtd:
    input:
        expand("src/data/multizone/gaussian/earlyburst/{dtd}",
               dtd=["prompt", "powerlaw_slope11", "exponential_timescale15",
                    "plateau_width10", "triple"]
        ),
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/age_ofe_dtd.pdf"
    script:
        "src/scripts/age_ofe_dtd.py"

# Tables
rule multizone_scores:
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
        "src/data/multizone/scores.csv"
    cache:
        True
    script:
        "src/scripts/score_multizone_outputs.py"

rule summary_table:
    input:
        "src/data/multizone/scores.csv",
        "src/scripts/summary_table_header.txt"
    output:
        "src/tex/output/summary_table.tex"
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

rule age_sample_size:
    input:
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/age_sample_size.txt"
    script:
        "src/scripts/age_sample_size.py"
