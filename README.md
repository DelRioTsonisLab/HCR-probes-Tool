# HCR PROBE CUSTOM DESIGN

## Overview

**HCR Probe Custom Design** is a streamlined tool for designing probe sets for Hybridization Chain Reaction (HCR) in-situ experiments, with the goal of **maximizing the number of probes** found for your target gene.

This tool is a **re-worked and extended version of the HCR3.0 Probe Maker** originally developed by [The Ozpolat Lab](https://bduyguozpolat.org/). 
The original software was released under the **GPL-3.0 license**, and this version maintains that license.

---

## Features

> This custom-designed pipeline increases the number of probes in a user-friendly Jupyter Notebook.

---

## Requirements

- Python 3.10.12+
- Jupyter Notebook 7.0.4
- `blastn` from **NCBI BLAST+ v2.14.1+**
- packages.txt (complete list of installed Python packages)

---

## Example Use

Run *HCR_Probe_Tool.ipynb* to design a max number of probes:
1. See the included `Vignette.html` for a full example of how the tool was used with PRDM16 cDNA and Gallus gallus as the reference species.

To perform off-target filtering using BLAST:
1. Download and install [NCBI BLAST+ 2.14.1+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/)
2. Download the cDNA reference FASTA file for your species of interest.  
   For example, for chicken (Gallus gallus): *Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa*
