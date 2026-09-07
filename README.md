# Proteome modulation by opposite inotropic drugs in human engineered cardiac tissue revealed by topology-driven cross-modal integration

Data analytics pipeline combining untargeted MS-based proteomics and topological data analysis to reveal how positive and negative inotropic drugs modulate the proteome in human engineered cardiac tissues.

This repository contains two scripts:

1. [`TDArunSingleOmics_nonLog2_only24h.py`](TDArunSingleOmics_nonLog2_only24h.py) — builds the topological connectivity network (TCN) from annotated proteomics data
2. [`ProcessNodesData.py`](ProcessNodesData.py) — computes node-based correlations between proteins and metadata

## Requirements
Developed and validated with Python 3.10

## Python dependencies
Pipeline was tested with the package versions specified in `requirements.txt`.

## Installation
Install the required dependencies with:

```bash
pip install -r requirements.txt
```

## Usage

Run the scripts in order from the repository root.

### 1. Build the TCN

[`TDArunSingleOmics_nonLog2_only24h.py`](TDArunSingleOmics_nonLog2_only24h.py) builds a topological connectivity network (TCN) from annotated proteomics data.

```bash
python TDArunSingleOmics_nonLog2_only24h.py
```

- **Input:** annotated data and metadata in [`DATA`](DATA)
- **Output** ([`RESULTS`](RESULTS)): `TCN_*.html` (interactive TCN), `Nodes_*.csv` (sample composition of TCN nodes)

### 2. Compute node-based correlations

[`ProcessNodesData.py`](ProcessNodesData.py) uses the node information from script 1 to calculate node-based correlations between proteins and metadata.

```bash
python ProcessNodesData.py
```

- **Input:** node information in [`RESULTS`](RESULTS) from script 1
- **Output** ([`RESULTS`](RESULTS)): `AnnData_Nodes_*.h5ad`, `Corr_*.h5ad`, `Supp_TNE.csv`, `Supp_TNE_ProteinCorrelations.csv`

## Citation

If you use this code, please cite:

Staykova DK, Snippert D, Wessels HJCT, Passier R, Conte F. Proteome modulation by opposite inotropic drugs in human engineered cardiac tissue revealed by topology-driven cross-modal integration. bioRxiv [Preprint]. 2026 Aug 27. [https://doi.org/10.64898/2026.08.25.746955](https://doi.org/10.64898/2026.08.25.746955)