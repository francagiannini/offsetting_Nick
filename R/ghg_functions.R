## ghg_functions.R
##
## Helper functions that turn rCTOOL (the R package implementation of the
## C-TOOL soil carbon model; Taghizadeh-Toosi et al., 2014;
## https://github.com/francagiannini/rCTOOL) into the "AGWP" (Absolute
## Global Warming Potential) accounting framework for soil C sequestration
## and its associated N2O emission, as described in the manuscript
## "C sequestration accounting" (Methods: "Accounting for carbon and
## nitrogen", and "Case study": "Modelling C and N dynamics in the soil").
##
## Three things are implemented here that rCTOOL itself does NOT provide,
## because they are this manuscript's own methodological contribution:
##
##   1. A SINGLE-LAYER wrapper around rCTOOL's (two-layer, topsoil +
##      subsoil) pool functions, reproducing the paper's Methods statement:
##      "the transport parameter (tF) was set to zero and the proportion of
##      HUM scheduled for transport to the subsoil (1-fCO2-fROM) was
##      recycled to the HUM pool of the topsoil." rCTOOL's own HUM_top
##      transport term is independent of ftr (it is always
##      substrate*(1-f_romi-f_co2)), so "recycling" it is mathematically
##      equivalent to simply not removing it from HUM_top - see
##      `run_ctool_pulse()` below.
##
##   2. A nitrogen-mineralisation / N2O sub-model (paper's Eq. 6-8) that
##      rides on top of the carbon fluxes rCTOOL computes (FOM/HUM
##      decomposition, humification, romification), since rCTOOL tracks
##      carbon only.
##
##   3. The AGWP/GWP metric framework itself (paper's Eq. 1-5 for CO2,
##      Eq. 6-9 for N2O; after Petersen et al., 2013 and Aamaas et al.,
##      2013), built on the CO2 impulse-response function of Joos et al.
##      (2013) and the N2O perturbation lifetime of Forster et al. (2021).
##
## Throughout, rCTOOL's OWN exported functions are used for every carbon
## flux calculation (FOM_top_calculations(), HUM_top_calculations(),
## ROM_top_calculations(), soil_config()) - nothing about the C-TOOL
## kinetics itself is reimplemented.

suppressPackageStartupMessages({
  library(rCTOOL)
  library(dplyr)
})

# ============================================================================
# 1. SINGLE-LAYER rCTOOL PULSE WRAPPER
# ============================================================================

#' Simulate the fate of a single pulse (or an annually repeated pulse) of
#' organic carbon input, using rCTOOL's topsoil pool functions in a
#' single-layer configuration.
#'
#' @param Cin_total total C input per pulse (kg; any mass unit is fine, the
#'   model is linear - the manuscript uses kg for a single 1000 kg C
#'   "case study" input rather than rCTOOL's usual Mg/ha field-scale units)
#' @param is_manure  TRUE for a manure input (split between FOM and HUM at
#'   input, per rCTOOL's own management_config()/f_man_humification
#'   mechanism); FALSE for a plant/crop-residue input (all to FOM)
#' @param f_man_humification fraction of manure C humified immediately on
#'   input (rCTOOL parameter; default 0.192 - see notes in
#'   run_case_study.R on how this value was cross-checked against the text)
#' @param CN_input C:N ratio of the fresh organic matter input
#' @param CN_hum,CN_rom fixed C:N ratios of the HUM and ROM pools (Table 1: 11.0)
#' @param s_config a soil configuration list from rCTOOL::soil_config(),
#'   plus one extra element EF_N2O (kg N2O-N (kg N mineralised)-1)
#' @param Tavg length-12 vector of mean monthly air temperature (deg C)
#' @param n_years number of years to simulate
#' @param input_years integer vector of simulation years (1-based) in which
#'   a pulse of Cin_total is applied (default: just year 1, i.e. a single
#'   pulse; pass 1:n_years for the repeated-input scenario)
#' @param input_month month (1-12) within each input year that the pulse
#'   is applied (default 1: the whole input occurs at the start of the
#'   simulation year, consistent with treating it as a discrete pulse for
#'   the AGWP decomposition in Petersen et al., 2013)
#'
#' @return a data.frame with one row per month: step, month, year,
#'   FOM_top/HUM_top/ROM_top (pool sizes after that month), em_CO2 (CO2-C
#'   emitted that month, summed over the three pools) and N2ON_emitted
#'   (N2O-N emitted that month)
run_ctool_pulse <- function(Cin_total,
                             is_manure = FALSE,
                             f_man_humification = 0.192,
                             CN_input,
                             CN_hum = 11.0,
                             CN_rom = 11.0,
                             s_config,
                             Tavg,
                             n_years = 100,
                             input_years = 1,
                             input_month = 1) {

  stopifnot(length(Tavg) == 12)
  months <- n_years * 12

  FOM_top <- 0; HUM_top <- 0; ROM_top <- 0
  FOM_N   <- 0; HUM_N   <- 0

  out <- vector("list", months)

  for (i in seq_len(months)) {
    mon <- ((i - 1) %% 12) + 1
    yr  <- ((i - 1) %/% 12) + 1
    t_avg     <- Tavg[mon]
    amplitude <- 0   # CRAN rCTOOL (>= 3.1.0) renamed this argument from
                       # t_range to amplitude and rewrote .soil_temp() to add
                       # a genuine damped seasonal sine term on top of t_avg
                       # (run_ctool() itself auto-computes it as
                       # (max(Tavg)-min(Tavg))/2 from the monthly climatology).
                       # amplitude = 0 is used deliberately here rather than
                       # that auto-estimate: it reproduces the manuscript's
                       # quoted case-study numbers noticeably more closely
                       # (see the validation notes in run_case_study.R) than
                       # passing the "full" seasonal amplitude does - see the
                       # note in run_case_study.R for the numbers behind this
                       # call. Set amplitude <- (max(Tavg) - min(Tavg)) / 2 to
                       # use the package's own default behaviour instead.

    is_input_month <- (yr %in% input_years) && (mon == input_month)
    Cin_plant <- if (!is_manure && is_input_month) Cin_total else 0
    Cin_man   <- if ( is_manure && is_input_month) Cin_total else 0

    # ---- nitrogen entering with the C pulse (mirrors the C split) ----
    if (is_manure && is_input_month) {
      N_total      <- Cin_man / CN_input
      N_to_HUM_dir <- (Cin_man * f_man_humification) / CN_hum
      N_to_FOM_dir <- N_total - N_to_HUM_dir
    } else if (!is_manure && is_input_month) {
      N_to_FOM_dir <- Cin_plant / CN_input
      N_to_HUM_dir <- 0
    } else {
      N_to_FOM_dir <- 0
      N_to_HUM_dir <- 0
    }

    # ---- FOM: carbon via rCTOOL's own exported function ----
    FOM_top_in <- FOM_top + Cin_plant + Cin_man * (1 - f_man_humification)
    FOM_N_in   <- FOM_N + N_to_FOM_dir

    fom <- FOM_top_calculations(FOM_top_t = FOM_top_in, month = mon,
                                 t_avg = t_avg, amplitude = amplitude,
                                 s_config = s_config)

    frac_lost_FOM     <- if (FOM_top_in > 0) fom$substrate_FOM_decomp_top / FOM_top_in else 0
    N_released_FOM    <- FOM_N_in * frac_lost_FOM
    N_to_HUM_from_FOM <- fom$FOM_humified_top / CN_hum
    # Eq. 6: if there is not enough N in the decomposing FOM to satisfy the
    # C:N ratio of the newly formed HUM, the shortfall is taken from soil
    # mineral N and no N is mineralised (net) from this step.
    N_min_FOM  <- max(N_released_FOM - N_to_HUM_from_FOM, 0)
    FOM_N_out  <- FOM_N_in - N_released_FOM

    # ---- HUM: carbon via rCTOOL's own exported function ----
    HUM_top_in <- HUM_top + fom$FOM_humified_top + Cin_man * f_man_humification
    HUM_N_in   <- HUM_N + N_to_HUM_from_FOM + N_to_HUM_dir

    hum <- HUM_top_calculations(HUM_top_t = HUM_top_in, month = mon,
                                 t_avg = t_avg, amplitude = amplitude,
                                 s_config = s_config)

    # *** single-layer modification (paper's Methods) ***
    # rCTOOL's HUM_tr = substrate * (1 - f_romi - f_co2) is the amount the
    # two-layer model would ship to the subsoil, independent of ftr.
    # Recycling it straight back into topsoil HUM is equivalent to simply
    # not subtracting it - hum$HUM_top already has it subtracted, so add
    # it back:
    HUM_top_out <- hum$HUM_top + hum$HUM_tr

    frac_lost_HUM     <- if (HUM_top_in > 0) hum$substrate_HUM_decomp_top / HUM_top_in else 0
    N_released_HUM    <- HUM_N_in * frac_lost_HUM
    N_to_ROM_from_HUM <- hum$HUM_romified_top / CN_rom
    # Eq. 7: same shortfall logic as for FOM
    N_min_HUM  <- max(N_released_HUM - N_to_ROM_from_HUM, 0)
    HUM_N_out  <- HUM_N_in - N_released_HUM

    # ---- ROM: carbon via rCTOOL's own exported function ----
    # (N mineralisation from ROM decomposition is not tracked further: over
    # the century-scale horizon used here, k_ROM makes this negligible, and
    # the manuscript defines Nmin only for the FOM and HUM steps.)
    ROM_top_in <- ROM_top + hum$HUM_romified_top
    rom <- ROM_top_calculations(ROM_top_t = ROM_top_in, month = mon,
                                 t_avg = t_avg, amplitude = amplitude,
                                 s_config = s_config)

    # Eq. 8: N2O-N emission from the mineralised N
    N_N2O_N <- s_config$EF_N2O * (N_min_FOM + N_min_HUM)

    out[[i]] <- data.frame(
      step = i, month = mon, year = yr,
      FOM_top = fom$FOM_top, HUM_top = HUM_top_out, ROM_top = rom$ROM_top,
      em_CO2 = fom$em_CO2_FOM_top + hum$em_CO2_HUM_top + rom$em_CO2_ROM_top,
      N2ON_emitted = N_N2O_N
    )

    FOM_top <- fom$FOM_top; HUM_top <- HUM_top_out; ROM_top <- rom$ROM_top
    FOM_N   <- FOM_N_out;   HUM_N   <- HUM_N_out
  }

  dplyr::bind_rows(out)
}

# ============================================================================
# 2. AGWP / GWP FRAMEWORK (Eq. 1-9)
# ============================================================================

## ---- CO2: Joos et al. (2013) multi-model-mean impulse response function ----
## IRF_CO2(t) = a0 + a1*exp(-t/tau1) + a2*exp(-t/tau2) + a3*exp(-t/tau3)
## (Table 5 of Joos et al., 2013; the same parameterisation underlies the
## IPCC AR5/AR6 GWP metrics.)
irf_co2_params <- list(a0 = 0.2173, a1 = 0.2240, a2 = 0.2824, a3 = 0.2763,
                        tau1 = 394.4, tau2 = 36.54, tau3 = 4.304)

A_CO2       <- 1.76e-15  # W m-2 kg-1, radiative efficiency of CO2 (Forster et al., 2021 / IPCC AR6)
gamma_C_CO2 <- 44 / 12   # kg CO2 (kg C)-1
N2ON_to_N2O <- 44 / 28   # kg N2O (kg N2O-N)-1

#' Analytical integral of IRF_CO2(t) from 0 to H (years); this is
#' Equation 3's integral of Equation 2.
integral_IRF_CO2 <- function(H, p = irf_co2_params) {
  H <- pmax(H, 0)
  p$a0 * H +
    p$a1 * p$tau1 * (1 - exp(-H / p$tau1)) +
    p$a2 * p$tau2 * (1 - exp(-H / p$tau2)) +
    p$a3 * p$tau3 * (1 - exp(-H / p$tau3))
}

#' Eq. 3: AGWP of a pulse of CO2 (or, used generically below, of any gas
#' mass unit) integrated over time horizon H
AGWP_CO2 <- function(H) A_CO2 * integral_IRF_CO2(H)

#' Eq. 4: AGWP of the CO2 emitted from sequestered C, where beta_t is a
#' named numeric vector (names = "years after input", t = 0, 1, 2, ...)
#' giving the fraction of the original C input emitted as CO2 in year t.
#' Each year's pulse is truncated to the REMAINING (H - t) years, exactly
#' as Eq. 3 is modified in the text.
AGWP_CO2seq <- function(beta_t, H) {
  t <- as.numeric(names(beta_t))
  A_CO2 * sum(beta_t * integral_IRF_CO2(H - t))
}

#' Eq. 5: negative GWP of sequestered C, kg CO2e (kg C input)-1.
#' Bounded between 0 (all input re-emitted immediately) and -3.67 (all
#' input retained for >= H years).
GWP_Cseq <- function(beta_t, H) {
  agwp_seq   <- AGWP_CO2seq(beta_t, H)
  agwp_pulse <- AGWP_CO2(H)
  gamma_C_CO2 * (agwp_seq - agwp_pulse) / agwp_pulse
}

## ---- N2O: single-exponential IRF (Aamaas et al., 2013) ----
## tau_N2O = 109 yr perturbation lifetime (Forster et al., 2021)
tau_N2O <- 109
A_N2O   <- 3.83e-13  # W m-2 kg-1 (Table 1: chosen so pulse GWP100(N2O) = 273, IPCC AR6)

#' Eq. 6's integral: proportion-remaining IRF for N2O integrated from 0 to H
integral_IRF_N2O <- function(H) {
  H <- pmax(H, 0)
  tau_N2O * (1 - exp(-H / tau_N2O))
}

#' Eq. 7 integrated: AGWP of a pulse of N2O over horizon H
AGWP_N2O <- function(H) A_N2O * integral_IRF_N2O(H)

#' Eq. 8: AGWP of the (delayed) N2O emitted from mineralisation of the
#' organic N input, eta_t = named vector, names = "years after input",
#' giving the fraction of total N2O emitted in year t
AGWP_N2Oseq <- function(eta_t, H) {
  t <- as.numeric(names(eta_t))
  A_N2O * sum(eta_t * integral_IRF_N2O(H - t))
}

#' Eq. 9: effective GWP of the delayed N2O emission, kg CO2e (kg N2O
#' emitted)-1. Bounded above by the standard pulse GWP100 (~273 for
#' H = 100 with the A_N2O value above); delay always REDUCES it, since
#' more of the N2O's warming falls outside the H-year assessment window.
GWP_N2Oseq <- function(eta_t, H) {
  AGWP_N2Oseq(eta_t, H) / AGWP_CO2(H)
}

## ---- instantaneous (radiative-forcing-shaped) time courses, for Fig 8/9 ----

#' Instantaneous AGWP-scale warming trace from a set of CO2 pulses.
#' beta_t is the *normalised* fraction of the C input emitted each year
#' (names = years after input, as elsewhere); total_co2_mass converts
#' that fraction into an actual CO2 gas mass (kg) so the result is on the
#' same physical scale as instantaneous_warming_pulse() below (pass
#' Cin_total * gamma_C_CO2 for a pulse scenario's beta_t).
instantaneous_warming_CO2 <- function(beta_t, years_out, total_co2_mass, p = irf_co2_params) {
  t_emit <- as.numeric(names(beta_t))
  vapply(years_out, function(y) {
    dt  <- y - t_emit
    irf <- ifelse(dt < 0, 0,
                  p$a0 + p$a1 * exp(-dt / p$tau1) + p$a2 * exp(-dt / p$tau2) + p$a3 * exp(-dt / p$tau3))
    A_CO2 * total_co2_mass * sum(beta_t * irf)
  }, numeric(1))
}

#' Instantaneous AGWP-scale warming trace for a set of N2O pulses.
#' eta_t is the *normalised* fraction of the total N2O emitted each year;
#' total_n2o_mass (kg N2O gas) puts the result on an absolute scale.
instantaneous_warming_N2O <- function(eta_t, years_out, total_n2o_mass) {
  t_emit <- as.numeric(names(eta_t))
  vapply(years_out, function(y) {
    dt  <- y - t_emit
    irf <- ifelse(dt < 0, 0, exp(-dt / tau_N2O))
    A_N2O * total_n2o_mass * sum(eta_t * irf)
  }, numeric(1))
}

#' Instantaneous warming had the same total mass been emitted as an
#' immediate pulse at t = 0 (the comparison curve in Fig. 8/9)
instantaneous_warming_pulse <- function(total_mass, years_out, gas = c("CO2", "N2O")) {
  gas <- match.arg(gas)
  if (gas == "CO2") {
    p <- irf_co2_params
    A_CO2 * total_mass * (p$a0 + p$a1 * exp(-years_out / p$tau1) +
                             p$a2 * exp(-years_out / p$tau2) + p$a3 * exp(-years_out / p$tau3))
  } else {
    A_N2O * total_mass * exp(-years_out / tau_N2O)
  }
}

# ============================================================================
# 3. SCENARIO SUMMARY (reproduces the logic behind Table 2)
# ============================================================================

#' Aggregate a run_ctool_pulse() output into annual CO2/N2O totals and the
#' Table-2-style GWP summary for a *single* pulse scenario.
summarize_scenario <- function(sim, Cin_total, H = 100) {

  ann <- sim %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(em_CO2 = sum(em_CO2), N2ON = sum(N2ON_emitted), .groups = "drop") %>%
    dplyr::mutate(t = year - 1)  # "years after input": input occurs in year 1 -> t = 0

  beta_t <- stats::setNames(ann$em_CO2 / Cin_total, ann$t)
  total_N2ON <- sum(ann$N2ON)
  total_N2O  <- total_N2ON * N2ON_to_N2O
  eta_t <- if (total_N2ON > 0) {
    stats::setNames(ann$N2ON / total_N2ON, ann$t)
  } else {
    stats::setNames(rep(0, nrow(ann)), ann$t)
  }

  gwp_c      <- GWP_Cseq(beta_t, H)
  gwp_n2o_ef <- GWP_N2Oseq(eta_t, H)

  C_offset   <- gwp_c * Cin_total
  N2O_effect <- gwp_n2o_ef * total_N2O
  net        <- C_offset + N2O_effect

  list(annual = ann, beta_t = beta_t, eta_t = eta_t,
       total_N2ON_kg = total_N2ON, total_N2O_kg = total_N2O,
       gwp_c_per_kgC = gwp_c, gwp_n2o_effective = gwp_n2o_ef,
       C_offset_kgCO2e = C_offset, N2O_effect_kgCO2e = N2O_effect,
       net_kgCO2e = net)
}
