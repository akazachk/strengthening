# Strengthening disjunctive cuts
# Aleksandr M. Kazachkov

### Instructions for submodule
Whenever the repository is cloned on a new machine, there are several steps to get the required dependencies (from `COIN-OR` and the `VPC` cut generator). The `VPC` generator is added as a submodule, which is kept on record in the `.gitmodules` file. After cloning the main repository, the user must run `git submodule init`, then `git submodule update`. Alternatively, when calling `git clone` initially, one can call `git clone --recurse-submodules`.

Afterwards, `cd` to the `vpc` directory, call `export PROJ_DIR=$PWD` from that directory, and call `./setup/install_coin.sh` to install `COIN-OR` files into `vpc/lib`.

Then, still in the `vpc` directory, modify the `makefile` to properly set the COIN and CPLEX directories and set `USE_CPLEX` to value `1`. Finally, call `make` then `make archive`.
