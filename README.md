# Strengthening disjunctive cuts
# Aleksandr M. Kazachkov

[![DOI](https://zenodo.org/badge/224919504.svg)](https://doi.org/10.5281/zenodo.14961495)

### Instructions for `Eigen`
Clone the `Eigen` library to your machine and properly link to it in the `Makefile` with the `EIG_LIB` variable.
    git clone https://gitlab.com/libeigen/eigen.git

### Instructions for `VPC` submodule
Note: The submodule for the `VPC` generator is not necessary if you have the `VPC` repository elsewhere on your machine and you properly set the `VPC_DIR` variable in the `Makefile`.

Whenever the repository is cloned on a new machine, there are several steps to get the required dependencies (from `COIN-OR` and the `VPC` cut generator). The `VPC` generator is added as a submodule, which is kept on record in the `.gitmodules` file. After cloning the main repository, the user must run `git submodule init`, then `git submodule update --remote`. Alternatively, when calling `git clone` initially, one can call `git clone --recurse-submodules`.

Afterwards, `cd` to the `vpc` directory, call `export PROJ_DIR=$PWD` from that directory, and call `./setup/install_coin.sh` to install `COIN-OR` files into `vpc/lib`.

Then, still in the `vpc` directory, modify the `makefile` to properly set the COIN and CPLEX directories and set `USE_CPLEX` to value `1`. Finally, call `make` then `make archive`.

### How to cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

Below is a BibTex entry for citing the paper.

```
@article{KazBal25+_strengthening-mpb,
    author    = {Kazachkov, Aleksandr M. and Balas, Egon},
    title     = {Monoidal Strengthening of Simple {$\mathcal{V}$}-Polyhedral Disjunctive Cuts},
    journal   = {Math. Program.},
    fjournal  = {Mathematical Programming},
    volume    = {},
    number    = {},
    pages     = {},
    month     = {},
    year      = {2025},
    doi       = {10.1007/s10107-024-02185-x},
    url       = {https://doi.org/10.1007/s10107-024-02185-x},
    note      = {Accepted and to appear},
}
```

Below is a BibTex entry for citing the latest release of the code.

```
@misc{KazBal25+_strengthening-github,
  author    = {Kazachkov, Aleksandr M. and Balas, Egon},
  title     = {Code for ``Monoidal Strengthening of Simple {$\mathcal{V}$}-Polyhedral Disjunctive Cuts''},
  year      = {2025},
  doi       = {10.5281/zenodo.14961496},
  url       = {https://github.com/akazachk/strengthening},
}
```

