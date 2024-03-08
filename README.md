# FFPE Impact

## Description
This R script is designed to analyze the impact of formalin-fixed paraffin-embedded (FFPE) artefacts on mutation signatures within genomic data. It leverages the capabilities of the `signature.tools.lib` library to fit substitution (SBS) and indel (ID) signatures to mutation catalogues derived from FFPE samples. The script calculates the FFPE impact on mutation signatures and adjusts the mutation catalogue to mitigate artefact effects.

## Prerequisites
To use this script, you will need to have R and the `signature.tools.lib` library installed. The `signature.tools.lib` library is crucial as it provides the necessary functions for fitting mutation signatures and other related analyses. If you do not have `signature.tools.lib` installed, you can find installation instructions here:

https://github.com/Nik-Zainal-Group/signature.tools.lib

## Input Requirements

Before running the script, you must define the following parameters:

- `SNV_path`: The file path to the single nucleotide variant (SNV) catalogue.
- `ID_path`: The file path to the indel (ID) catalogue.
- `organ`: The organ type of the sample(s). This must be one of the organs recognized in the `FitMS` function from `signature.tools.lib`.

Make sure that the `SNV_path` and `ID_path` point to valid files containing your mutation data in a format readable by `read.table()`. The organ type should match one of the supported types in `signature.tools.lib` to ensure accurate signature fitting - it is now possible to specify ```organ="Other"``` in FitMS.

## Output

The script returns a list containing the following components:

- `SBS_exposures`: Exposure estimates for SBS signatures.
- `ID_exposures`: Exposure estimates for ID signatures.
- `FFPE_impact`: Calculated impact of FFPE artefacts on mutation signatures.
- `cleaned_catalogues`: The mutation catalogue adjusted for FFPE artefact effects.
- `newMHcount`: Proportions of microhomology-associated deletions post artefact adjustment.

These outputs facilitate further analysis of the FFPE artefact impact and the evaluation of mutation signatures in your genomic data.
