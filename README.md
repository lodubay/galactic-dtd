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
</p>

Welcome to the repository for **Dubay et al. (in prep), "Galactic Chemical 
Evolution Models Favor an Extended Type Ia Supernova Delay-Time Distribution"**.

To re-build the article yourself, simply run the following from the repository
directory:
```
$ showyourwork build
```
The data and model outputs will automatically be downloaded from Zenodo: 
[Zenodo DOI]

To re-run all models yourself, run the following:
```
$ python run_all_models.py
```
This will replace everything in the `src/data/multizone` directory, including
the output files from Zenodo.

To re-create the APOGEE sample with Leung et al. (2023) age estimates,
run the following:
```
$ python generate_sample.py
```
This will replace everything in the `src/data/APOGEE` directory.
