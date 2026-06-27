# Pinc <a href='https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb'><img src='hexlogo.png' align="right" height="198" /></a>: a simple probabilistic AlphaFold interaction score

<!-- badges: start -->
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb)
<!-- badges: end -->

## Description
Pinc (**P**robability of **i**nterface **n**ative **c**ontacts) is a simple interaction confidence metric derived from an AlphaFold model and its predicted aligned error (PAE) matrix.
A Pinc score of 0.8 means that approximately 80% of native interface contacts are predicted to be present in the model.
The method can be used on both AlphaFold2 and AlphaFold3 (including server) outputs.


## Using the Colab notebook
- [ ] Open the interactive Colab notebook [here](https://colab.research.google.com/gist/badonyi/4c976a79efbae4b64ad99d968eedd417/pinc.ipynb) or click "Open in Colab" above
- [ ] Run the first cell (▶) to enable file upload
- [ ] Upload a structure file (PDB/mmCIF) its matching PAE matrix (JSON)
- [ ] Execute subsequent cells step by step (▶)
- [ ] Download the results as a zipped folder in the final cell


## Command-line usage
Clone the repository and compile the program. A C compiler (gcc/clang) is required and is pre-installed on most systems:

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
    R script to compute the Pinc score from an AlphaFold PAE
    matrix (JSON) and the corresponding structure file (PDB/CIF).

 Usage:
    Rscript Pinc.R <json_file> <structure_file> [options]

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


Assuming you have a structure file `model_0.cif` and an associated JSON data file `model_0.json` with the PAE matrix, to calculate Pinc scores for all chain pairs you can call the program simply by:

```bash
./Pinc model_0.json model_0.cif
```

The program will generate a `model_0_Pinc.csv` file with columns `chain1`, `chain2`, `Pinc1`, `Pinc2`, and `Pinc`, where chain pairs are non-redundant (for A-B there is no B-A).
`Pinc1` and `Pinc2` are the mean contact probabilities viewed from each chain's own PAE column, and `Pinc` is their symmetric average.


If you need the full contact probability matrix, add the `--all` flag:

```bash
./Pinc model_0.json model_0.cif --all
```

In this case the program will additionally generate a `model_0_contact_probability.json` file with node names `token_chain_ids` and `contact_probability`, similar to the AlphaFold3 output.
Note that `token_chain_ids` define the shape of the square matrix that can be reconstructed from the numeric array. 


If you need a list of all non-zero probabilities at the token level, call:

```bash
./Pinc model_0.json model_0.cif --pairlist
```

The program will additionally generate a `model_0_pairlist.csv` file with columns `token1`, `token2`, `contact_p`, and `distance` columns.
The tokens represent residues for proteins and bases for nucleic acids, but can also stand for atoms in ligands or ions.
The pair list is non-redundant (for chains A-B there is no B-A), and represent the mean of the two contacts (i->j and j->i).
The distance is the centre-of-mass distance between the tokens given in Ångströms.


Note that `--all` and `--pairlist` flags _can_ be used together.


If you would like to specify a results folder, use:

```bash
./Pinc model_0_full_data.json model_0.cif --results /usr/home/path_to_pinc_folder
```

Make sure the directory has write access.


## Legacy R script
An R implementation ([`Pinc.R`](https://git.mpi-cbg.de/tothpetroczylab/Pinc/-/blob/main/Pinc.R)) is retained in the repository for reference.
It requires R (≥ 3.5.0) and the `jsonlite` package.
Note that repeated calls make it substantially slower than the C program due to interpreter overhead.
Usage is identical to the C version, substituting `Rscript Pinc.R` for `./Pinc`:

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