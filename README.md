# fake-factor-application
Tools to create ntuples containing fake factor weights.

## Overview
This repository contains the tools to create friend tree ntuples with extrapolation factors needed for the fake factor estimation method.
The process specific extrapolation factors are determined centrally and downloaded by the checkout script. The process fractions that are needed to calculate the inclusive extrapolation factor can be taken from centrally produced workspaces if applicable or produced with the [shape-producer](https://github.com/KIT-CMS/shape-producer) via
```bash
./fake-factor-application/produce_shapes.sh $ERA $VARIABLE
```
The latter is in particular useful if use a very specific variable in your analysis and the fractions should be binned in this variable as well. Define your binning in `fake-factor-application/config.py` using the variable as key.

The friend trees are produced with
```bash
./fake-factor-application/create_fake_factor_friends.sh $ERA $VARIABLE $CATEGORYMODE $OUTPUTDIR
```
where `$VARIABLE` may also just be a label for a composite expression like 'njets_mvis' if implemented in the code. `$CATEGORYMODE` can be set to `inclusive` if the fake factors are to be calculated for the inclusive analysis selection without any categorization applied. Any other expression will currently use the NN categorization. Note that the output directory must not exist yet as the tool refuses to overwrite an existing directory. In the shell scripts you can further define input folders and whether a fractions workspace or fractions from own production should be used.

The central functionality of the fake factor calculation is defined in some functions in `fake-factor-application/python/calculate_fake_factors.py` which can be used for integration into other frameworks as done in the [friend-tree-producer](https://github.com/KIT-CMS/friend-tree-producer).

The whole chain can be run via
```bash
./fake-factor-application/run_fake_factors.sh $ERA $VARIABLE $CATEGORYMODE $OUTPUTDIR
```

## Setting up the software
If this repository is not embedded into a cmssw environment, you need checkout a minimal cmssw setup for the fake factor utilities:
```bash
./fake-factor-application/utils/init_fake_factors.sh
```
This script also calls
```bash
./fake-factor-application/utils/get_fractions_workspaces.sh
```
which you can use to update the RooWorkspaces that contain the fake factor fractions if necessary.