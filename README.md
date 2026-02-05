# Longitudinal Protein Stability and Individual Specificity Analysis Project

## Project Structure

### Core Analysis Scripts

- **001-PHI-different-visit.R**: Calculate Pearson correlation coefficients for proteins across 6 visit pairs (02-07, 02-12, 02-22, 07-12, 07-22, 12-22)
- **002-PHI-visualization.R**: Visualize PHI correlation results with distribution plots and scatter plots
- **003-PHI-grouped-correlation.R**: Calculate PHI correlations grouped by sex and age
- **004-PHI-grouped-visualization.R**: Visualize grouped PHI results with violin plots and statistical tests
- **005-PHI-pQTL-visualization.R**: Analyze relationships between PHI and pQTL effect sizes
- **006-GeneticVariance-PHI-visualization.R**: Analyze relationships between PHI and genetic variance
- **007-GxE-PHI-visualization.R**: Analyze relationships between PHI and GxE interaction effects
- **008-PHI-same-gene-different-seqid-scatter.R**: Compare PHI values for different protein isoforms of the same gene
- **009-protein-class-PHI-visualization.R**: Analyze PHI distributions across protein classes (intracellular, extracellular/membrane, secreted)
- **010-protein-sankey-retention.R**: Generate Sankey plots showing protein level transitions across visits with retention ratios
- **011-PHI-R2-analysis.R**: Analyze relationships between PHI and longitudinal prediction accuracy (R²)

### Visit Pairs

The project analyzes protein stability across 4 time points:

- **Dataset 1**: 2002 (02)
- **Dataset 2**: 2007 (07)
- **Dataset 3**: 2012 (12)
- **Dataset 4**: 2022 (22)

Six visit pairs are analyzed: 02-07, 02-12, 02-22, 07-12, 07-22, 12-22

### PHI Bins

Proteins are categorized into PHI bins for analysis:

- 0-0.25: Low stability
- 0.25-0.5: Moderate-low stability
- 0.5-0.75: Moderate-high stability
- 0.75-1: High stability

## Requirements

### R Packages

- dplyr
- ggplot2
- tidyr
- ggrepel
- cowplot
- ggsignif
- ggalluvial
- lme4
- lmerTest
- parallel

## Usage

1. **Data Preparation**: Ensure input data files are in the correct format and paths are updated in each script
2. **Run Analysis Scripts**: Execute scripts in numerical order  for complete analysis
3. **Review Results**: Check output directories for generated figures and CSV files

## Output

Each script generates:

- CSV files with analysis results
- PNG figures for visualizations
- Summary statistics and test results

Output directories are organized by analysis type:

- `results/`: Main results directory
- `age-sex-seperate-result/`: Grouped analysis results
- `pQTL-PHI/`: pQTL-related analyses
- `GeneticVariance-PHI/`: Genetic variance analyses
- `GxE-PHI/`: GxE interaction analyses
- `Person-Setpoint/`: Prediction accuracy analyses
- `Sankey-Plot/`: Sankey diagram visualizations

## Notes

- File paths in scripts are hardcoded and may need to be adjusted for different systems
- Some scripts require specific input files from previous analysis steps
- Batch effect removal and normalization should be performed before PHI calculations

## Citation

If you use this code, please cite the relevant publications and acknowledge the data sources.


