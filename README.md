## Replication Data and Code for "Measuring Political Polarization: The Cluster-Polarization Coefficient"

# Usage

This repository contains the data and code necessary to replicate all results in both the main text and supplementary information of "Measuring Political Polarization: The Cluster-Polarization Coefficient," Isaac D. Mehlhaff (February 2021).

After downloading or cloning the repository to a local destination, ensure the working directory in the .R script is set to whatever directory contains the data files. Then, run the script. The .R script is meant to be run from top to bottom and will produce, over the course of the script, every figure reported in the main text and supplementary information as well as the information to fill in the table reported in the main text. Code chunks which produce figures and tables are commented as such.

Note that the Monte Carlo simulations are meant to be run together and in the order provided in the script. For example, to reproduce results from the two-cluster simulations, begin with the univariate simulations and then move to the bivariate simulations. This is to ensure that random seeds - which are specified at the beginning of each of the two-, three-, and four-cluster sections - are set consistently. Although different seeds in the Monte Carlo simulations will produce slightly different results, they would not change the substantive conclusions drawn from those results.

# Citation

To cite this code in publications and working papers, please use:

Mehlhaff, Isaac D. "Measuring Political Polarization: The Cluster-Polarization Coefficient," working paper (February 2021).

For BibTeX users:

```
@unpublished{Mehlhaff2021,
  title = {Measuring {{Political Polarization}}}: {{The Cluster}}-{{Polarization Coefficient}}},
  author = {Mehlhaff, Isaac D.},
  date = {2021-02},
  location = {{The University of North Carolina at Chapel Hill}}
}
```

To cite this data in publications and working papers, please refer to the relevant data owners:

Chartbook of Economic Inequality: https://www.chartbookofeconomicinequality.com

DW-NOMINATE: https://voteview.com
