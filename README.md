

# iwa-miRNA

[![web](https://img.shields.io/badge/docker-ready-yellow.svg)](https://hub.docker.com/r/malab/iwa-mirna/) ![docker pull](https://img.shields.io/docker/pulls/malab/iwa-mirna.svg) [![web](https://img.shields.io/badge/Web server-ready-blue.svg)](https://deepngs.nwafu.edu.cn/) 

## User-friendly framework for interactive web-based annotation of plant miRNAs

- iwa-miRNA is a Galaxy-based bioinformatics platform for interactively annotating miRNAs from plant genomes. It can be used to generate a comprehensive collection of miRNA candidates, and allows users to interrogate miRNA annotation in a straightforward way, without the need for computational skills.

- The iwa-miRNA project is hosted on GitHub (https://github.com/cma2015/iwa-miRNA). A demo server is available at [https://deepngs.nwafu.edu.cn](https://deepngs.nwafu.edu.cn/). We suggest users to run iwa-miRNA locally using the Docker image (https://hub.docker.com/r/malab/iwa-mirna). This repository contains source code, and example data to help users run and test iwa-miRNA.

The following figure shows the main components of the analysis process:

<img src="assets/img/Graphical_summary.png" alt="Graphical summary of iwa-miRNA" style="zoom:25%">

## Installation

Step 1: [Docker installation](https://github.com/cma2015/PEA/blob/master/tutorial/docker_installation.md)

Step 2: Docker Image Installation

```bash
# pull image from docker hub
docker pull malab/iwa-mirna
docker run -it -p 4000:8080 malab/iwa-mirna /bin/bash
sh /home/galaxy-release_20.05/run.sh
```

## Overview of functional modules in iwa-miRNA

iwa-miRNA consists of three primary modules: 

- [**Module I: MiRNA Compilation**](Tutorials/ModuleI.md)
- [**Module II: MiRNA Selection**](Tutorials/ModuleII.md)
- [**Module III: Manual Curation**](Tutorials/ModuleIII.md)

## How to use iwa-miRNA

- Tutorials for iwa-miRNA: https://github.com/cma2015/iwa-miRNA/tree/master/Tutorials
- Test datasets referred in the tutorials for iwa-miRNA: https://github.com/cma2015/iwa-miRNA/tree/master/Test_data

## Changelog

- 2020-06-24: A demo server of iwa-miRNA was released for users running small small RNA sequencing datasets.
- 2020-06-20: Source code and Docker image of of iwa-miRNA were released for the first time.

## How to access help

- For any feedback and tool suggestions, please feel free to leave a message at Github [issues](https://github.com/cma2015/iwa-miRNA/issues). We will try our best to deal with all issues as soon as possible.
- In addition, if any suggestions are available, feel free to contact: **Zhang Ting** [zting135@gmail.com](mailto:zting135@gmail.com) or ***Ma Chuang*** [chuangma2006@gmail.com](mailto:chuangma2006@gmail.com)

## Citation

Ting Zhang, Jingjing Zhai, Xiaorong Zhang, Lei Ling, Menghan Li, Shang Xie, Minggui Song, Chuang Ma. Interactive web-based annotation of plant miRNAs with iwa-miRNA. (Submitted)