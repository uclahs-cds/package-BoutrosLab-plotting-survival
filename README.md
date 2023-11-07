# BoutrosLab.plotting.survival
`BoutrosLab.plotting.survival` is a collections of survival analysis visualization tools, including Kaplan-Meier curves, progonsis plots and schoenfeld residual plots.

## Installation
### Package Dependency
The following packages need to be built before installation:
- [BoutrosLab.plotting.general](https://github.com/uclahs-cds/public-R-BoutrosLab-plotting-general)
- [BoutrosLab.statistics.general](https://github.com/uclahs-cds/public-R-BoutrosLab-statistics-general)
- [BoutrosLab.statistics.survival](https://github.com/uclahs-cds/public-R-BoutrosLab-statistics-survival)
- [BoutrosLab.prognosticsignature.general](https://github.com/uclahs-cds/public-R-BoutrosLab-prognosticsignature-general)
- [BoutrosLab.utilities](https://github.com/uclahs-cds/public-R-BoutrosLab-utilities)
- [rbenchmark](https://cran.r-project.org/web/packages/rbenchmark/index.html)

### Using devtools in R:
```R
library(devtools)
install_github('https://github.com/uclahs-cds/public-R-BoutrosLab-plotting-survival')
```

### From source:
```shell script
git clone git@github.com:uclahs-cds/public-R-BoutrosLab-plotting-survival.git
R CMD INSTALL public-R-BoutrosLab-plotting-survival
```
### Getting Started
Here, we used the cancer dataset from `survival` package, which contains the basics of the survival analysis. For more details, please see the `survival` CRAN page [here](https://cran.r-project.org/web/packages/survival/).

```R
library('survival');
library('BoutrosLab.plotting.survival');

# Load the testing data
data(cancer, package='survival');
cancer.data <- data.frame(cancer);

# Plot the Kaplan-Meier plot to visualize survival curves
create.km.plot(
  survival.object = Surv(cancer.data$time, cancer.data$status),
  patient.groups = as.factor(cancer.data$sex)
);

# Plot the Competing Risk Cumlative Incidence Function
library('riskRegression');
library('cmprsk');
data(Melanoma);

create.cr.plot(
  data = Melanoma,
  formula = time ~ status,
  patient.groups = epicel,
  filename = 'test.pdf'
);

```
### Plot Gene Meta Prognosis
Example analysis can be found [here](https://github.com/uclahs-cds/BoutrosLab/blob/792c95d41d4a09bd7d42c43d660b6c240ce186a3/Biomarkers/ProstateCancer/cirDNA_Methylation/PilotStudy/compare_to_primary_data/prostate.cancer_univariate.survival.analysis.R)

### Plot multifeature meta prognosis
Example analysis can be found [here](https://github.com/uclahs-cds/BoutrosLab/blob/792c95d41d4a09bd7d42c43d660b6c240ce186a3/Collaborators/RodBremner/SB32/Hypoxia_gene_signature.R)

### Plot Schoenfeld residual plot
Example analysis can be found [here](https://github.com/uclahs-cds/BoutrosLab/blob/792c95d41d4a09bd7d42c43d660b6c240ce186a3/CPC-GENE/mitochondrial_survey/test_create_kmplot_3legend.R)


