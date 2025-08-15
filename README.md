# BayesTemporalEpidemicDynamics

---

## Description

**BayesTemporalEpidemicDynamics** is a flexible **Bayesian modelling framework** for epidemic dynamics with **time-varying transmission rates**. It implements **SEMIKR** and **SIKR** mechanistic models parameterised via **B-splines**, allowing the study of interventions, behavioural changes, and stochastic effects in disease spread.

This repository contains:

* Implementations in **HAICS**, **Stan**, and **LibBi**.
* A robust **MAP-centred initialisation** strategy for efficient Hamiltonian Monte Carlo sampling.
* Scripts and data to **reproduce all experiments** from the associated paper, including synthetic and real COVID-19 case studies.

---

## Features

* Time-dependent transmission rates in SEMIKR and SIKR models.
* Differentiable B-spline parameterisation for efficient HMC inference.
* Diffusion-based stochastic modelling via LibBi.
* Under-reporting correction through a time-varying detection function.
* Ready-to-run experiments from the paper.
* Includes synthetic datasets and COVID-19 incidence data for Spanish provinces.
* Well-documented initialisation routines for stable convergence.

---

## Repository Structure

```
initialization_routine.R        # Robust initialisation procedure
run_Stan_from_initialization.R  # Stan sampling after init
run_LibBi_simulations.R         # Diffusion-based fitting via LibBi
run_HAICS.sh                    # HAICS simulation launcher
haics/simulation/src/           # C code for SEMIKR and SIKR models
Stan/model/                     # Stan model files
LibBi/models/                   # Diffusion-based epidemic models
data/                           # Synthetic & real COVID-19 incidence data
auxiliary_functions/            # Helper functions
```

---

## Installation

### HAICS

```bash
sudo apt-get install libgsl-dev
sudo apt-get install yad
# Install Sundials (to /usr/local)
# https://sundials.readthedocs.io/en/latest/Install_link.html
# Install R: https://cran.r-project.org/
```

### Stan

1. Install [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) or [Stan](https://mc-stan.org/) following the official instructions.
2. Then install the R interface:

```R
install.packages("rstan")
```

### LibBi

1. Install [LibBi](https://libbi.org/) following the platform-specific instructions.
2. Then install the R interface:

```R
install.packages("rbi")
```

---

## Usage

1. **Initialisation**

```R
source("initialization_routine.R")
```

2. **HMC Sampling**

* HAICS:

```bash
parallel -j 4 ./run_haics_examples.sh testSEI3R_02 ::: {1..10}
```

* Stan:

```R
source("run_Stan_from_initialization.R")
```

3. **Diffusion-based Sampling**

```R
source("run_LibBi_simulations.R")
```

4. **Results & Diagnostics**
   All scripts save outputs and produce visual summaries.

---

## Input / Output

* **Input:** model configuration, epidemiological parameters, initial values, and incidence data.
* **Output:** posterior samples, derived quantities (e.g., $R_0(t)$), visual diagnostics.

---

## Data

* **Synthetic:** generated with SEI3R dynamics for validation.
* **Real:** COVID-19 daily incidence data for Spanish Autonomous Communities (e.g., Basque Country), with under-reporting corrections.

---

## Authors

* Hristo Inouzhe Valdes — BCAM / UAM
* María Xosé Rodríguez Álvarez — Universidade de Vigo / CITMAga
* Lorenzo Nagar — BCAM
* Elena Akhmatskaya — BCAM / IKERBASQUE

---

## License

MIT License — see [LICENSE](LICENSE) for details.
