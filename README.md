# Dongdong-Zhao
# MaCAST integration efficiency analysis

This package contains a Python script to quantify CAST-mediated integration from amplicon FASTQ data.

## License
MIT License (see LICENSE).

## System requirements
- OS: Linux/macOS/Windows (tested on Ubuntu 22.04)
- Python: 3.9+
- Dependencies: Biopython, XlsxWriter

## Installation
```bash
pip install -r requirements.txt
```

## Demo
Run on the included demo files:
```bash
python integration_tally.py --csv demo/sample_info.csv --fastq demo/test1.fastq.gz --excel demo_output.xlsx
```

Expected output:
- An Excel file reporting read counts, integration efficiency, and the insertion-distance distribution.

## Usage on your data
Prepare a single-sample `sample_info.csv` with columns:
- `name`
- `crrna` (spacer sequence)
- `donor_re` (donor RE sequence)
- `genome_seq` (reference amplicon sequence containing the spacer)

Provide the corresponding FASTQ(.gz) file and run the command above with your paths.

## Notes
- Integration efficiency is reported as: `integrated_reads / (integrated_reads + unedited_reads)`.
- A high-level pseudocode description is provided in `docs/pseudocode.md`.
