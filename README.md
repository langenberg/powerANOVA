# powerANOVA: Power Analysis for Repeated Measures ANOVA.
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

`powerANOVA` is an R package intended as a companion to the article "A Tutorial on Using the Paired *t*-Test for Power Calculations in Repeated Measures ANOVA".


## 1 Installation
`powerANOVA` is not on CRAN. The development version of `powerANOVA ` can be installed directly from this GitHub repository using the additional package `devtools`. 

```
# install devtools
install.packages("devtools")

# install subgroupsem
devtools::install_github("langenberg/powerANOVA")
```

## 2 Usage

First, load the package using the following command:

```
library(powerANOVA)
```

### 2.1 Graphical User Interface (GUI)

`powerANOVA` includes an easy to use GUI. The GUI is a shiny app which can be loaded using the following command:

```
power_gui()
```

The command will open a browser tab. The GUI is very self-explaining. Simply close the browser when you want to terminate the GUI.

### 2.2 Converters

* `convert_cohens_d_petasq(cohens_d, n, population = TRUE)`
* `convert_petasq_cohens_d(p_eta_sq, n, population = TRUE)`
* `convert_petasq_f2(p_eta_sq)`
* `convert_F_petasq(Fval, df1, df2)`
* `convert_petasq_F(p_eta_sq, df1, df2)`


### 2.3 Plots

* `power_plot_cohens_d(n, cohens_d)`
* `power_plot_mu_cov_contrast(n, mu, Sigma, contrast)`
* `power_plot_p_eta_sq(n, p_eta_sq, df1)`


### 2.4 Power Calculators

* `power_cohens_d(n, cohens_d)`

	E.g.:

	```
	mu <- 58
	var <- 7200
	
	cohens_d <- mu/sqrt(var)
	
	power_cohens_d(n = 19, cohens_d = cohens_d)
	```

* `power_petasq(n, p_eta_sq, df1)`

	E.g.:

	```
	p_eta_sq <- 0.3296746
	
	power_petasq(n = 19, p_eta_sq = p_eta_sq)
	```


* `power_mu_cov_contrast(n, mu, Sigma, contrast)`

	E.g.:
	
	```
	means_dv <- matrix(c(492, 511, 483, 444), ncol = 1)
	names(means_dv) <- c("A1.B1", "A1.B2", "A2.B1", "A2.B2")
	
	# (co)variances of the dependent variables
	vcov_dv <- matrix(c(
	    9000, 7200, 7200, 7200,
	    7200, 9000, 7200, 7200,
	    7200, 7200, 9000, 7200,
	    7200, 7200, 7200, 9000
	), ncol = 4, nrow = 4)
	
	colnames(vcov_dv) <- c("A1.B1", "A1.B2", "A2.B1", "A2.B2")
	rownames(vcov_dv) <- c("A1.B1", "A1.B2", "A2.B1", "A2.B2")
	
	# contrast vector
	contrast_vec <- matrix(c(1, -1, -1, 1), nrow = 1)
	
	colnames(contrast_vec) <- c("A1.B1", "A1.B2", "A2.B1", "A2.B2")
	rownames(contrast_vec) <- c("difference")
	
	power_mu_cov_contrast(n = 19, mu = means_dv, Sigma = vcov_dv, contrast = contrast_vec)
	```
	
	
	
	
	
	