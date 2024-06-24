<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/lodubay/galactic-dtd/actions/workflows/build.yml">
<img src="https://github.com/lodubay/galactic-dtd/actions/workflows/build.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/lodubay/galactic-dtd/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/lodubay/galactic-dtd/raw/main-pdf/dag.pdf">
<img src="https://img.shields.io/badge/article-dag-blue.svg?style=flat" alt="Article graph"/>
</a>
<a href="https://github.com/lodubay/galactic-dtd/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
<a href="https://doi.org/10.5281/zenodo.12207380">
<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.12207380.svg" alt="DOI">
</a>
</p>

Welcome to the repository for **Dubay et al. (2024), "Galactic Chemical 
Evolution Models Favor an Extended Type Ia Supernova Delay-Time Distribution"**,
[arXiv:2404.08059](https://arxiv.org/abs/2404.08059).

To re-build the article yourself, simply run the following from the repository
directory:
```
$ showyourwork build
```
The data and model outputs will automatically be downloaded from Zenodo deposit
[10.5281/zenodo.10961090](https://zenodo.org/doi/10.5281/zenodo.10961090).

To re-create the APOGEE sample with Leung et al. (2023) age estimates,
run the following:
```
$ python generate_sample.py
```
This will replace everything in the `src/data/APOGEE` directory.

To re-run all models yourself, run the following:
```
$ python run_all_models.py
```
This will replace everything in the `src/data/multizone` directory, including
the output files from Zenodo.

Source code for the models is located within the [src/scripts/multizone/](/src/scripts/multizone/) directory.
To run a single multi-zone model with custom parameters, run the following:
```
$ cd src/scripts
$ python -m multizone [OPTIONS...]
```
In particular, code for the stellar migration prescription is contained within the `gaussian_migration` class
located at
[src/scripts/multizone/src/migration.py](/src/scripts/multizone/src/migration.py).

The figures in the paper do not cover every possible combination of SFH + DTD. The following script
generates supplementary plots for every multi-zone model available:
```
$ cd src/scripts
$ python extra_plots.py
```
The multi-panel plots show age vs [O/Fe], [O/Fe] vs [Fe/H], MDFs, and [O/Fe] DFs for each Galactic region,
plus a plot showing the [O/Fe] bimodality. Plots are located in the [src/extra](/src/extra/) directory with the 
following structure:
`<migration scheme>/<sfh model>/<dtd model>/`.
