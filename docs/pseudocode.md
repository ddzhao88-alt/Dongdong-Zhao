## Pseudocode: integration efficiency analysis

**Inputs**
- `sample_info.csv` (single sample): `name`, `crrna` (spacer), `donor_re`, `genome_seq`
- `FASTQ(.gz)` file for the same sample

**Steps**
1. **Derive an unedited genomic marker** from `genome_seq`:
   - Locate the spacer (allowing up to one sequencing error; reverse-complement fallback).
   - Extract the 20-bp window from -25 to -5 upstream of the spacer.

2. **Process each read**:
   - Orient the read by finding the spacer (allowing up to one sequencing error; reverse-complement fallback).
   - Classify as **integrated** if a donor-end tag (first 20 bp of `donor_re`) is detected in the expected window downstream of the spacer.
   - Otherwise, classify as **unedited** if the unedited genomic marker is detected in the read.
   - Reads that cannot be oriented by spacer are counted as **waste**.

3. **Compute integration efficiency**:
   - `integration_efficiency = integrated_reads / (integrated_reads + unedited_reads)`
   - Report counts and the insertion-distance distribution.

Note: Full implementation details are provided in `integration_tally.py`.
