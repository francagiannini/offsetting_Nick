# Reproducing the C-sequestration case study with rCTOOL

A fully reproducible implementation of the manuscript's *Case study*
section (soil C sequestration + N2O, accounted for with the AGWP/GWP
framework of Petersen et al., 2013 and its N2O extension), built on the
[rCTOOL](https://github.com/francagiannini/rCTOOL) package.

## Files

- `R/ghg_functions.R` - all the reusable logic: a single-layer wrapper
  around rCTOOL's pool functions, the N-mineralisation/N2O sub-model, and
  the AGWP/GWP metric framework (Eq. 1-9 of the manuscript). Heavily
  commented with exactly what is original C-TOOL/rCTOOL vs. this
  contribution's own reconstruction.
- `run_case_study.R` - a plain script; run it top to bottom
  (`source("run_case_study.R")` or `Rscript run_case_study.R` from this
  directory) to regenerate every figure/table into `figures/` and
  `tables/`.
- `case_study_report.Rmd` - the same analysis as a narrative R Markdown
  report; knit it (`rmarkdown::render("case_study_report.Rmd")`) for a
  self-contained HTML report with figures and tables inline, including a
  side-by-side comparison against the numbers quoted in the manuscript
  text.

## Setup

```r
install.packages("remotes")
remotes::install_github("francagiannini/rCTOOL")
install.packages(c("dplyr", "ggplot2", "tidyr", "knitr", "rmarkdown"))
```

`rCTOOL` Depends on `easyclimate` (for pulling real climate data - not
needed for this case study, which specifies its own monthly temperatures
directly) and `data.table`; installing it will pull those in too.

Then, from this directory:

```r
source("run_case_study.R")
# or
rmarkdown::render("case_study_report.Rmd")
```

## What's reproduced, and how faithfully

- **Fig. 6, 7** (CO2 and N2O emission trajectories): reproduce the
  qualitative pattern described in the text (identical CO2 curves for the
  two plant C:N ratios; manure's initial dip and later N2O peak; delayed,
  flattened curves under the cool climate).
- **Fig. 8, 9** (instantaneous AGWP-scale warming vs. a pulse) and
  **Table 2** (net GWP effect): after one correction (see below), the
  C-sequestration offsets match the manuscript's quoted numbers to within
  about 1% for three of the four scenarios, and the N2O effective GWPs to
  within about 3-10%. See `case_study_report.Rmd`'s "Sanity checks" section
  for the full comparison table.
- A **bonus** repeated-annual-input scenario for the Discussion section's
  unfinished "Figure xx", which is under-specified in the draft - the
  implementation there is this contribution's own reasonable reading of
  that paragraph, clearly flagged as such.

## One likely typo found in Table 1

Table 1 of the draft lists `k_HUM = 0.028`. Used literally, that value
makes the HUM pool decompose ~10x faster than rCTOOL's own published
default (`0.0028`) and produces sequestration offsets 3-4x smaller than
the ones quoted in the text. Using `0.0028` instead reproduces the text's
own numbers closely (see above) - strongly suggesting `0.028` is missing a
leading zero. This is controlled by a single named constant
(`k_hum_value`) at the top of both `run_case_study.R` and the Rmd, so it's
a one-line change either way.

## Known assumptions / things to double-check against your own code

- `f_man_humification = 0.192` is not stated numerically in the draft; it
  is back-solved from the text's statement that the manure input's FOM
  C:N ratio works out to 17.9:1 (see the Rmd for the derivation), and
  happens to match rCTOOL's own real-world example value.
- The whole 1000 kg C (and, for manure, N) input is applied in a single
  month (month 1 of year 1), treating it as a discrete pulse - reasonable
  for the AGWP pulse-decomposition framework, but worth confirming against
  whatever the original simulation did.
- N mineralisation is only modelled for the FOM and HUM decomposition
  steps (matching the manuscript's Eq. 6-7); N release from ROM
  decomposition is not tracked, since ROM's decomposition rate makes it
  negligible over a 100-year horizon.
- CO2 and N2O radiative properties (Joos et al. 2013 IRF parameters,
  A_CO2 = 1.76e-15 W/m2/kg, N2O perturbation lifetime = 109 yr) are
  standard IPCC AR6-consistent values, chosen to match the manuscript's
  own stated A_N2O = 3.83e-13 (calibrated to GWP100(N2O) = 273).
