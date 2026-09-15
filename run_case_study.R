## ============================================================================
## Reproducing the "Case study" of the manuscript
## "C sequestration accounting" with the rCTOOL package
## ============================================================================
##
## This script implements, end to end and starting from a fresh R session,
## the worked example in the manuscript's Case study section:
##
##   - Fig. 6  CO2 emissions in the 10 years following a single 1000 kg C
##             input (plant material at two C:N ratios, manure, and plant
##             material under a cooler climate)
##   - Fig. 7  associated N2O(-N) emissions
##   - Fig. 8  instantaneous AGWP-scale warming from the sequestered-C CO2
##             emissions vs. an immediate (pulse) emission of the same mass
##   - Fig. 9  the same, for the N2O emitted from mineralisation of the
##             co-input organic N
##   - Table 2 the net GWP effect (C sequestration offset + N2O emission)
##             over a 100-year time horizon, in kg CO2e
##
## plus a bonus scenario (Discussion section, "Figure xx") comparing the
## conventional IPCC stock-change method against the AGWP-based method for
## a REPEATED annual C input.
##
## rCTOOL provides the carbon-turnover engine (FOM/HUM/ROM pool functions);
## everything to do with turning those carbon fluxes into the paper's own
## AGWP/GWP metric, and the N2O extension, lives in R/ghg_functions.R and
## is applied here. See that file for detailed documentation of exactly
## what is original C-TOOL/rCTOOL and what is this manuscript's own layer
## on top of it.
##
## ----------------------------------------------------------------------
## IMPORTANT NOTE ON TABLE 1 OF THE DRAFT
## ----------------------------------------------------------------------
## The draft's Table 1 lists k_HUM = 0.028 kg (kg C)-1 (confirmed by
## inspecting the .docx XML directly - it is not a text-extraction
## artefact). Taken literally, that value decomposes the HUM pool roughly
## 10x faster than rCTOOL's own published default (k_hum = 0.0028; see the
## package's basic_example / README), and produces CO2/N2O trajectories
## and C-sequestration GWP offsets that are 3-4x smaller in magnitude than
## the ones quoted in the manuscript text (the -0.167/-0.283/-0.302 kg
## CO2e (kg C)-1 sequestration offsets, and the associated N2O effective
## GWPs of 222/256/244/241 kg CO2e (kg N2O)-1).
##
## Using k_HUM = 0.0028 instead, with "years after input" measured from
## t = 0 at the start of the input year (as done in summarize_scenario()
## in R/ghg_functions.R), reproduces the manuscript's quoted C-sequestration
## offsets to within ~1% for three of the four scenarios (Plant C:N=16:
## -166.1 vs. -167; Plant, cool climate: -283.1 vs. -283; Manure: -300.3
## vs. -302; Plant C:N=50, same C dynamics as C:N=16: -166.1 vs. -167) and
## the N2O effective GWPs to within ~3-10%. This is a strong enough match
## that "0.028" in the draft table is very likely missing a leading zero.
## k_hum is set from a single named constant below (`k_hum_value`) so this
## is a one-line change if you determine the table is correct after all
## and the discrepancy lies elsewhere.
## ----------------------------------------------------------------------

## ----------------------------------------------------------------------
## NOTE ON THE INSTALLED rCTOOL VERSION (CRAN vs. GitHub)
## ----------------------------------------------------------------------
## This was originally developed against rCTOOL's GitHub "main" branch
## (v3.0.0). The CRAN release (v3.1.0, francagiannini/Serra/da Silva) is
## somewhat different: FOM_top_calculations()/HUM_top_calculations()/
## ROM_top_calculations() renamed their seasonal-temperature argument from
## t_range to amplitude, and rewrote the internal soil-temperature function
## to use proper seconds-based physics with a real damped seasonal sine
## term (the GitHub version's t_range argument was silently almost inert -
## a leftover bug CRAN's release notes describe fixing). run_ctool_pulse()
## in R/ghg_functions.R now targets the CRAN version's argument name
## (`amplitude`) and, like before, passes amplitude = 0 deliberately rather
## than the package's own auto-estimated seasonal amplitude
## ((max(Tavg)-min(Tavg))/2): using the real amplitude changes the
## C-sequestration offsets by a further ~6-17% (worst for the cool-climate
## scenario) and moves them further FROM the manuscript's quoted numbers,
## not closer, so amplitude = 0 (soil temperature = air temperature, no
## depth lag) remains the better-supported choice here. Since amplitude = 0
## makes the damping term vanish regardless of the thermal-diffusivity
## value, this choice is unaffected by the CRAN default temp_th_diff
## differing 10x from the Table 1 value (0.35e-6 vs. 0.035e-6) - see the
## soil_config() call below.
## ----------------------------------------------------------------------

## ---- 0. Setup -------------------------------------------------------------

# install.packages("remotes")
# install.packages(c("dplyr", "ggplot2", "tidyr", "knitr", "rCTOOL"))

library(rCTOOL)
library(dplyr)
library(ggplot2)
library(tidyr)

here <- function(...) file.path(getwd(), ...)
source(here("R", "ghg_functions.R"))
# 
# dir.create(here("figures"), showWarnings = FALSE)
# dir.create(here("tables"), showWarnings = FALSE)

theme_set(theme_classic(base_size = 12))

## ---- 1. Parameters from the manuscript's Table 1 --------------------------

k_hum_value <- 0.0028   # see note above; set to 0.028 to use the draft's literal value

s_config <- soil_config(
  Csoil_init = 0,        # not used by the single-layer pulse wrapper (pools start at 0)
  Cproptop   = 1,
  f_hum_top  = 0.4,       # not used (no initialize_soil_pools() call for a pulse input)
  f_rom_top  = 0.4,       # ditto
  f_hum_sub  = 0.4,       # ditto (subsoil is not used at all in the single-layer wrapper)
  f_rom_sub  = 0.4,       # ditto
  clay_top   = 0.10,      # Table 1: X = 0.10 g g-1
  clay_sub   = 0.10,      # not used
  phi        = 0.035,     # legacy parameter, unused by CRAN rCTOOL (>= 3.1.0)'s soil-temperature
                           #  function - kept only because soil_config() still accepts it
  temp_th_diff = 0.035e-6, # Table 1: "0.035 x10^-6" - this is what CRAN rCTOOL actually reads
                            #  for soil-temperature damping (see run_ctool_pulse()'s `amplitude`
                            #  argument, passed as 0 below: at amplitude = 0 this value is inert,
                            #  see the note there for why 0 was chosen)
  f_co2      = 0.628,     # Table 1: f_CO2
  f_romi     = 0.012,     # Table 1: f_ROM
  k_fom      = 0.12,      # Table 1: k_FOM
  k_hum      = k_hum_value,
  k_rom      = 3.85e-5,   # Table 1: "3.85 x10^-5"
  ftr        = 0          # single-layer configuration: transport parameter tF = 0
)
s_config$EF_N2O <- 0.01   # Table 1: EF_N2O = 0.01 kg N2O-N (kg N mineralised)-1

CN_hum <- 11.0   # Table 1
CN_rom <- 11.0   # Table 1

# Monthly mean air temperature, degC (Table 1). Row order is exactly as given
# in the draft; the annual means (14.33 and 8.08 degC) match the 14.3/8.1 degC
# quoted in the text for the standard and "cool climate" scenarios respectively,
# confirming month order/values were transcribed correctly from the table.
T_standard <- c(19, 22, 17, 16, 11, 11, 9, 8, 10, 13, 17, 19)
T_cool     <- c(14.8, 16.5, 16.0, 12.8, 9.0, 4.0, 1.8, 0.89, 1.17, 3.5, 6.0, 10.5)
stopifnot(abs(mean(T_standard) - 14.33) < 0.01, abs(mean(T_cool) - 8.08) < 0.01)

# f_man_humification (fraction of manure C humified immediately on input) is
# NOT given a numeric value anywhere in the draft's Table 1 or text. It is
# cross-checked here instead: the text states that with this parametrisation,
# the manure's FOM sub-pool ends up with a C:N ratio of "17.9:1" (not the bulk
# manure input ratio of 16:1). Given Manure C = 1000 kg, N = 1000/16 = 62.5 kg,
# and N committed directly to HUM at the fixed HUM C:N ratio of 11, the FOM
# C:N ratio of 17.9 is reproduced almost exactly by f_man_humification = 0.192
# - which is also rCTOOL's own real-world example value (Askov straw
# incorporation) - so that is the value used here.
f_man_humification <- 0.192
{
  Cman <- 1000; CNman <- 16
  f <- f_man_humification
  N_tot <- Cman / CNman
  N_hum <- f * Cman / CN_hum
  CN_fom_check <- (Cman * (1 - f)) / (N_tot - N_hum)
  message(sprintf("Check: implied manure FOM C:N ratio = %.1f (text: 17.9)", CN_fom_check))
}

H <- 100   # time horizon, years (Eq. 1: H = 100 per UNFCCC reporting convention)
n_years_sim <- 150   # simulate well beyond H so beta_t/eta_t past H (which contribute
                       # zero to the truncated AGWP integrals) are captured completely

Cin_total <- 1000  # kg C, the manuscript's single-input case study

## ---- 2. Run the four scenarios --------------------------------------------

scenario_defs <- list(
  "Plant (C:N 16:1)" = list(CN = 16, is_manure = FALSE, temps = T_standard),
  "Plant (C:N 50:1)" = list(CN = 50, is_manure = FALSE, temps = T_standard),
  "Manure"           = list(CN = 16, is_manure = TRUE,  temps = T_standard),
  "Plant, cool climate" = list(CN = 16, is_manure = FALSE, temps = T_cool)
)

sims <- lapply(scenario_defs, function(sc) {
  run_ctool_pulse(
    Cin_total = Cin_total, is_manure = sc$is_manure,
    f_man_humification = f_man_humification,
    CN_input = sc$CN, CN_hum = CN_hum, CN_rom = CN_rom,
    s_config = s_config, Tavg = sc$temps,
    n_years = n_years_sim, input_years = 1, input_month = 1
  )
})

summaries <- Map(function(sim) summarize_scenario(sim, Cin_total, H = H), sims)

## ---- 3. Figure 6: CO2 emissions in the first 10 years ---------------------

fig6_df <- bind_rows(lapply(names(sims), function(nm) {
  sims[[nm]] %>%
    mutate(scenario = nm, t_yr = (step - 1) / 12) %>%
    filter(t_yr <= 10) %>%
    select(scenario, t_yr, em_CO2)
}))

fig6 <- ggplot(fig6_df, aes(t_yr, em_CO2, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Years after C input", y = "CO2-C emission (kg / month)",
       colour = NULL,
       title = "Fig. 6 - CO2 emission in the 10 years following the input of 1000 kg C") +
  theme(legend.position = "bottom")

ggsave(here("figures", "fig6_CO2_emissions.png"), fig6, width = 7.5, height = 5, dpi = 200)

## ---- 4. Figure 7: N2O(-N) emissions in the first 10 years ------------------

fig7_df <- bind_rows(lapply(names(sims), function(nm) {
  sims[[nm]] %>%
    mutate(scenario = nm, t_yr = (step - 1) / 12) %>%
    filter(t_yr <= 10) %>%
    select(scenario, t_yr, N2ON_emitted)
}))

fig7 <- ggplot(fig7_df, aes(t_yr, N2ON_emitted, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Years after C input", y = expression(N[2] * O * "-N emission (kg / month)"),
       colour = NULL,
       title = "Fig. 7 - N2O-N emission in the 10 years following the input") +
  theme(legend.position = "bottom")

ggsave(here("figures", "fig7_N2O_emissions.png"), fig7, width = 7.5, height = 5, dpi = 200)

## ---- 5. Figure 8: instantaneous AGWP-scale warming, CO2 -------------------

years_out <- seq(0, H, by = 0.5)

pulse_CO2_mass <- Cin_total * gamma_C_CO2  # kg CO2, same total mass emitted at t=0

fig8_df <- bind_rows(lapply(names(summaries), function(nm) {
  s <- summaries[[nm]]
  data.frame(scenario = nm, year = years_out,
             warming = instantaneous_warming_CO2(s$beta_t, years_out, total_co2_mass = pulse_CO2_mass))
}))
fig8_pulse <- data.frame(
  scenario = "Pulse (immediate) emission", year = years_out,
  warming = instantaneous_warming_pulse(pulse_CO2_mass, years_out, gas = "CO2")
)

fig8 <- ggplot(bind_rows(fig8_df, fig8_pulse), aes(year, warming, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Years after input", y = expression("Instantaneous warming (W " * m^-2 * ")"),
       colour = NULL, linetype = NULL,
       title = "Fig. 8 - Instantaneous global warming from the CO2 emitted from sequestered C\nvs. an immediate pulse emission of the same mass") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))

ggsave(here("figures", "fig8_instantaneous_warming_CO2.png"), fig8, width = 8, height = 5.5, dpi = 200)

## ---- 6. Figure 9 (a-c): instantaneous AGWP-scale warming, N2O -------------

fig9_df <- bind_rows(lapply(names(summaries), function(nm) {
  s <- summaries[[nm]]
  data.frame(scenario = nm, year = years_out,
             warming = instantaneous_warming_N2O(s$eta_t, years_out, total_n2o_mass = s$total_N2O_kg))
}))

fig9_pulse <- bind_rows(lapply(names(summaries), function(nm) {
  s <- summaries[[nm]]
  data.frame(scenario = nm, year = years_out,
             warming = instantaneous_warming_pulse(s$total_N2O_kg, years_out, gas = "N2O"))
}))
fig9_df$series   <- "Sequestered (delayed) N2O"
fig9_pulse$series <- "Pulse (immediate) N2O"

fig9 <- ggplot(bind_rows(fig9_df, fig9_pulse), aes(year, warming, colour = series)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~scenario, ncol = 2, scales = "free_y") +
  labs(x = "Years after input", y = expression("Instantaneous warming (W " * m^-2 * ")"),
       colour = NULL,
       title = "Fig. 9 - Instantaneous global warming from the N2O emitted following mineralisation\nof the co-input organic N, vs. an immediate pulse emission of the same mass") +
  theme(legend.position = "bottom")

ggsave(here("figures", "fig9_instantaneous_warming_N2O.png"), fig9, width = 9, height = 7, dpi = 200)

## ---- 7. Table 2: net GWP effect over H = 100 years ------------------------

table2 <- bind_rows(lapply(names(summaries), function(nm) {
  s <- summaries[[nm]]
  data.frame(
    `Organic matter input` = nm,
    `C seq. offset (kg CO2e)` = round(s$C_offset_kgCO2e, 0),
    `N2O emission (kg CO2e)`  = round(s$N2O_effect_kgCO2e, 0),
    `Net effect (kg CO2e)`    = round(s$net_kgCO2e, 0),
    `Effective GWP of N2O`    = round(s$gwp_n2o_effective, 0),
    check.names = FALSE
  )
}))

print(table2)
write.csv(table2, here("tables", "table2_net_effect.csv"), row.names = FALSE)

## ---- 8. Bonus: repeated annual input vs. the conventional stock-change method
##        (Discussion section, "Figure xx")
## ----------------------------------------------------------------------
## The draft only sketches this scenario ("a scenario was created in which
## there is a step increase in C input that is repeated each year ... model
## parameters were the same as for the single input scenarios" - the figure
## itself is a placeholder ("Figure xx") with no numeric detail). The
## implementation below is this contribution's own reasonable reading of
## that paragraph, not a literal reproduction, and is clearly optional:
##
##   - "Conventional" (IPCC stock-change) method: each year's reported
##     negative emission is -(SOC stock this year - SOC stock last year) x
##     44/12, i.e. exactly the elemental-mass-balance approach the paper's
##     Introduction critiques.
##   - "AGWP" method: the correct, time-invariant annual offset is the one
##     already computed for a SINGLE 1000 kg C pulse (GWP_Cseq x Cin),
##     since (per the Discussion) that value "will remain constant ... even
##     if repeated each year".
##   - The over-estimation ratio is conventional / AGWP (both negative;
##     ratio > 1 means the conventional method overstates the cooling
##     effect in that year).
## ----------------------------------------------------------------------

repeat_years <- 60
sim_repeat <- run_ctool_pulse(
  Cin_total = Cin_total, is_manure = FALSE, f_man_humification = f_man_humification,
  CN_input = 16, CN_hum = CN_hum, CN_rom = CN_rom,
  s_config = s_config, Tavg = T_standard,
  n_years = repeat_years, input_years = 1:repeat_years, input_month = 1
)

soc_annual <- sim_repeat %>%
  group_by(year) %>%
  summarise(SOC_stock = last(FOM_top + HUM_top + ROM_top), .groups = "drop") %>%
  mutate(dSOC = SOC_stock - dplyr::lag(SOC_stock, default = 0),
         conventional_kgCO2e = -dSOC * gamma_C_CO2)

agwp_annual_offset <- summaries[["Plant (C:N 16:1)"]]$C_offset_kgCO2e  # constant, kg CO2e / yr

overestimate_df <- soc_annual %>%
  mutate(agwp_kgCO2e = agwp_annual_offset,
         overestimation_ratio = conventional_kgCO2e / agwp_kgCO2e)

fig_bonus <- ggplot(overestimate_df, aes(year, overestimation_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  labs(x = "Year since repeated annual input began", y = "Over-estimation ratio\n(conventional / AGWP)",
       title = "Bonus - Over-estimation of the conventional stock-change method\nrelative to the AGWP approach, for a repeated annual 1000 kg C input") +
  coord_cartesian(ylim = c(0, max(4.5, max(overestimate_df$overestimation_ratio, na.rm = TRUE))))

ggsave(here("figures", "fig_bonus_overestimation.png"), fig_bonus, width = 7.5, height = 5, dpi = 200)
write.csv(overestimate_df, here("tables", "bonus_overestimation.csv"), row.names = FALSE)

## ---- 9. Console summary ----------------------------------------------------

cat("\n================ Table 2 (net GWP effect, H = 100 yr) ================\n")
print(table2)

cat("\nFigures written to ./figures/, tables to ./tables/\n")
cat("Year-1 over-estimation ratio (bonus scenario):",
    round(overestimate_df$overestimation_ratio[1], 2), "\n")
