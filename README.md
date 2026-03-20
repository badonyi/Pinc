# Pinc <a href='https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb'><img src='hexlogo.png' align="right" height="198" /></a>: a simple probabilistic AlphaFold interaction score

<!-- badges: start -->
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb)
<!-- badges: end -->

## Description
Pinc (**P**robability of **i**nterface **n**ative **c**ontacts) is a simple-to-calculate protein–protein interaction confidence metric derived from AlphaFold model coordinates and the associated predicted aligned error (PAE) matrix. It provides an interpretable probabilistic estimate of interface reliability. For example, a Pinc score of 0.8 can be read as predicting that approximately 80% of the interfacial contacts observed in the native structure are present in the model.
The method can be used on both AlphaFold2 and AlphaFold3 (including server) outputs.


## Using the Colab notebook
- [ ] Open the interactive Colab notebook [here](https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb) or click "Open in Colab" above
- [ ] Run the first cell (▶) to enable file upload
- [ ] Upload the required files: **1**) a structure file (PDB or CIF) and **2**) its matching JSON file with the PAE matrix
- [ ] Execute subsequent cells step by step (▶)
- [ ] Download the results as a zipped folder from the final cell


## Command-line usage
- [ ] Clone the repository and compile the program with `make` (requires `gcc` or `clang`, available on Linux, Mac, and Windows):

```bash
git clone https://git.mpi-cbg.de/tothpetroczylab/Pinc.git
cd Pinc
make
```


## Getting started
Calling the program empty will display its usage information:

```bash
./Pinc
```

```txt
Description:
   A program to compute the Pinc score from an AlphaFold PAE
   matrix (JSON) and the corresponding structure file (PDB/CIF).

Usage:
   Pinc <json_file> <structure_file> [options]

Options:
   None        Pinc score for all unique chain pairs (default)
   --all       Full contact probability matrix
   --pairlist  Residue pair list of non-zero contact probabilities
   --results   Optional path to output folder

Output format:
   - default:  <structure_file_name>_Pinc.csv
   - all:      <structure_file_name>_contact_probability.json
   - pairlist: <structure_file_name>_pairlist.csv
```


## Usage
The `test` directory contains example input files you can use to test the below functionalities.

**1**) Assuming you have a structure file `model_0.cif` and an associated JSON data file `model_0.json` with the PAE matrix, to calculate Pinc scores for all chain pairs you can call the program simply by:

```
./Pinc model_0.json model_0.cif
```

The program will generate a `model_0_Pinc.csv` file with columns `chain1`, `chain2`, and `Pinc`, where chain pairs are non-redundant (for A-B there is no B-A).


**2**) If you need the full contact probability matrix, add the `--all` flag:

```
./Pinc model_0.json model_0.cif --all
```

The program will additionally generate a `model_0_contact_probability.json` file with node names `token_chain_ids` and `contact_probability`, similar to the AlphaFold3 output.
Note that `token_chain_ids` define the shape of the square matrix that can be reconstructed from the numeric array. 


**3**) If you need a list of all non-zero probabilities at the token-level, call:

```
./Pinc model_0.json model_0.cif --pairlist
```

The program will additionally generate a `model_0_pairlist.csv` file with columns `token1`, `token2`, `contact_p`, and `distance` columns.
The tokens essentially represent residues for proteins and bases for nucleic acids, but can also stand for atoms in ligands or ions.
The pair list is non-redundant (for A-B there is no B-A), and represent the mean of the two contacts (A->B and B->A).
The distance is the centre-of-mass distance between the tokens in Ångströms.


Note that `--all` and `--pairlist` flags _can_ be used together.


**4**) If you would like to specify a results folder, use:

```
./Pinc model_0_full_data.json model_0.cif --results /usr/home/path_to_pinc_folder
```

**Important**: Make sure the directory has write access.


## Legacy R script
The original R implementation ([`Pinc.R`](https://git.mpi-cbg.de/tothpetroczylab/Pinc/-/blob/main/Pinc.R)) is retained in the repository for reference, as it is the version described in the preprint. It requires R (≥ 3.5.0) and the `jsonlite` package. Usage is identical to the C version, substituting `Rscript Pinc.R` for `./Pinc`:

```bash
Rscript Pinc.R model_0.json model_0.cif
```


## Contact
Please report any technical problems and questions via [Issues](https://git.mpi-cbg.de/tothpetroczylab/Pinc/-/issues) or alternatively by email: <badonyi@mpi-cbg.de>


## Reference
```
@article{Pinc,
  title={Pinc: a simple probabilistic AlphaFold interaction score},
  author={Badonyi, Mihaly and Toth-Petroczy, Agnes},
  journal={bioRxiv},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```