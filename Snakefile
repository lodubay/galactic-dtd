# Multi-zone outputs
rule multizone_output:
    input:
        "src/data/multizone.tar.gz"
    output:
        directory("src/data/multizone/")
    script:
        "src/scripts/extract_multizone_output.py"

# Figures
rule star_formation_histories:
    input:
        "src/data/multizone"
    output:
        "src/tex/figures/star_formation_histories.pdf"
    script:
        "src/scripts/star_formation_histories.py"

rule radial_migration:
    input:
        "src/data/multizone"
    output:
        "src/tex/figures/radial_migration.pdf"
    script:
        "src/scripts/radial_migration.py"

rule midplane_distance:
    input:
        "src/data/multizone"
    output:
        "src/tex/figures/midplane_distance.pdf"
    script:
        "src/scripts/midplane_distance.py"
        
rule feh_df_comparison:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/feh_df_comparison.pdf"
    script:
        "src/scripts/feh_df_comparison.py"

rule ofe_df_sfh:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_df_sfh.pdf"
    script:
        "src/scripts/ofe_df_sfh.py"

rule ofe_df_dtd:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_df_dtd.pdf"
    script:
        "src/scripts/ofe_df_dtd.py"

rule ofe_bimodality_summary:
    input:
        "src/data/multizone",
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
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_sfh.pdf"
    script:
        "src/scripts/ofe_feh_sfh.py"

rule ofe_feh_dtd:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_dtd.pdf"
    script:
        "src/scripts/ofe_feh_dtd.py"

rule ofe_feh_twoinfall:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/ofe_feh_twoinfall.pdf"
    script:
        "src/scripts/ofe_feh_twoinfall.py"

rule age_ofe_sfh:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/age_ofe_sfh.pdf"
    script:
        "src/scripts/age_ofe_sfh.py"

rule age_ofe_dtd:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/figures/age_ofe_dtd.pdf"
    script:
        "src/scripts/age_ofe_dtd.py"

# Tables
rule multizone_scores:
    input:
        "src/data/multizone",
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/scores.csv"
    script:
        "src/scripts/score_multizone_outputs.py"

rule summary_table:
    input:
        "src/tex/output/scores.csv"
    output:
        "src/tex/output/summary_table.tex"
    script:
        "src/scripts/summary_table.py"

rule scores_table:
    input:
        "src/tex/output/scores.csv"
    output:
        "src/tex/output/scores_table.tex"
    script:
        "src/scripts/scores_table.py"

rule apogee_regions_table:
    input:
        "src/data/APOGEE/sample.csv"
    output:
        "src/tex/output/apogee_regions_table.tex"
    script:
        "src/scripts/apogee_regions_table.py"

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
