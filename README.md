# Modified_Bprime

In this repository we provide code to generate unified Bprime tables with positive and negative normalized mass flux.

## Dependencies
- [Mutation++](https://github.com/mutationpp/Mutationpp) (and dependencies therein). 
- Python modules `numpy` `pandas` `subprocess` `scipy`

## Instructions
- Install Mutation++ following their [installation guide](https://github.com/mutationpp/Mutationpp/blob/master/docs/installation.md#top)
- Clone the `Modified_Bprime` repo into the local directory `Modified_Bprime_cloned`
- `cd Modified_Bprime_cloned && cmake .` 
- `make generate_bprime_tables` (the file `generate_bprime_tables.cpp` is a modified version of the file `bprime.cpp` available in `${MPP_DIRECTORY}/src/apps`)
- `python3 generate_table.py`

## Debugging hints
- Check that Mutation++ is properly installed by running their [tests](https://github.com/mutationpp/Mutationpp/blob/master/docs/installation.md#top)
- If running on Linux, check if the Mutation++ library way installed in `${CMAKE_INSTALL_PREFIX}/lib` or `${CMAKE_INSTALL_PREFIX}/lib64`. Then, verify that the correct path is used in the `target_link_libraries()` command in the `CMakeLists.txt` file that ships with the `Modified_Bprime` repository.
