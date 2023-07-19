# Bprime

In this repository we provide code to generate unified Bprime tables with positive and negative normalized mass flux.

## Dependencies
- [Mutation++](https://github.com/mutationpp/Mutationpp) (and dependencies therein). 
- Python modules `numpy` `pandas` `subprocess` `scipy`

## Instructions
- Install Mutation++ following their [installation guide](https://github.com/mutationpp/Mutationpp/blob/master/docs/installation.md#top)
- Clone the `Bprime` repo into the local directory `Bprime_cloned`
- `cd Bprime_cloned && cmake .` 
- `make generate_bprime_tables` (the file `generate_bprime_tables.cpp` is a modified version of the file `bprime.cpp` available in `$MPP_DIRECTORY/src/apps`)
- `python3 generate_table.py`

