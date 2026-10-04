# Adsorption of Zn(II) and Mn(II) on white yam peel: reproducible analysis

Python reanalysis of the batch adsorption data from my B.Sc. research project,
*Kinetic, Isotherm and Thermodynamic Studies on the Adsorption of Zn(II) and
Mn(II) from Synthetic Aqueous Solutions Using White Yam (Dioscorea rotundata)
Peels as an Agricultural Adsorbent* (Department of Chemistry, Lagos State
University, May 2023).

Every number in the four notebooks is computed from the raw atomic absorption
readings in `data/aas_readings.csv`, which reproduce Appendix A of the thesis
exactly. Notebook 00 asserts twelve derived values against the figures printed
in the thesis tables, so a transcription or arithmetic error would fail the
run.

## What this repository is for

The point is not that three isotherms and four kinetic models can be fitted.
They can, and the code does it. The point is deciding which of those fits mean
anything on a dataset of this size, and reporting the ones that do not.

The headline results are negative, and they are the interesting part:

|Question        |Usual answer in the literature                         |What these data support                                                                                                                                     |
|----------------|-------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
|Kinetics        |Pseudo-second-order, R² > 0.95, therefore chemisorption|A constant-uptake plateau, better by ~5 AICc units. Uptake had already levelled off at the first measurement (30 min), so no rate constant is recoverable   |
|Mn(II) isotherm |Langmuir monolayer capacity                            |Langmuir degenerates to a straight line (K_L → 0). Freundlich fits best with 1/n = 1.47, which is the unfavourable regime                                   |
|Zn(II) isotherm |Whichever model has the highest R²                     |Langmuir, Freundlich and Temkin fit within 0.01 mg/g RMSE of each other. None can be preferred                                                              |
|Thermodynamics  |ΔH°, ΔS° and ΔG° to three decimals                     |Mn(II) exothermic trend, ΔH° = −25.6 ± 22.7 kJ/mol. For Zn(II) the regression is not significant (p = 0.30), so no values are reported                      |
|Effect of pH    |92 to 98% removal at pH 9 and 11                       |Above the pH where Zn(OH)₂ and Mn(OH)₂ can precipitate, so not attributable to adsorption. The defensible figures are 21.0% and 11.4%, both at pH 5         |
|Characterisation|75.52 wt% Zn and 60.96 wt% Mn by EDS                   |Zn confirmed. Mn **not** confirmed: the tallest peak matches Ca Kα, both Mn K lines are absent from the displayed range, and the composition sums to 104.97%|
|SEM micrographs |Particle and pore dimensions read from the images      |Field width rises with magnification where it must fall, and the scale bars disagree with it. No dimension is measurable; the images are qualitative only   |

Two specific checks drive those conclusions:

- **A permutation test on the linearised pseudo-second-order plot.** With five
  observations there are 120 possible assignments of the measured uptakes to
  the five times, so all of them are enumerated. The real Zn(II) data give
  R² = 0.78 on the `t/q_t` against `t` plot. The median over shuffled data is
  **0.83**. The plot cannot tell the measurements from noise, because time
  appears on both axes.
- **A rank-correlation test on the adsorbent-dose series.** Spearman
  ρ(C_e, q_e) = −0.90 for Mn(II), the opposite of what an isotherm requires, so
  that series is excluded from the isotherm fit and only the constant-dose
  concentration series is used.
- **A solubility-product screen on the pH series.** Zn(OH)₂ can begin forming
  near pH 7.1 and Mn(OH)₂ near pH 9.0 at the measured concentrations. Any
  precipitate would be caught on the filter paper and counted as removal, so
  the best-looking numbers in the dataset are set aside.
- **Matching EDS peak positions against tabulated X-ray emission lines.** The
  reported surface loadings are also 270 to 450 times the bulk loading measured
  independently by AAS, which is what identifies them as a local feature of one
  selected spot rather than a capacity.

## Notebooks

|File                             |Contents                                                                                                                                                                                                                                                                                       |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`00_data_and_controls.ipynb`     |Raw readings, the dilution-factor convention, measured versus nominal controls, a stoichiometric test of whether the unrecorded manganese salt hydration explains the gap, cross-checks against the thesis tables, and the batch-to-batch variability that sets the floor for every later claim|
|`01_isotherms.ipynb`             |Dataset selection and exclusions, non-linear Langmuir / Freundlich / Temkin fits, why AICc is undefined here, and the inflated R² of the linearised Langmuir plot                                                                                                                              |
|`02_kinetics.ipynb`              |Plateau, pseudo-first-order, pseudo-second-order, Elovich and intraparticle diffusion compared by AICc; the permutation test; sensitivity to the duplicate 120-min readings                                                                                                                    |
|`03_thermodynamics.ipynb`        |van ’t Hoff regression with 95% confidence intervals, ΔG° under both equilibrium-constant conventions, and sensitivity to the unassigned sixth temperature reading                                                                                                                             |
|`04_ph_and_precipitation.ipynb`  |Solubility-product screening for Zn(OH)₂ and Mn(OH)₂, the dilution-basis problem with the undiluted high-pH samples, and the pH range that is defensibly adsorption                                                                                                                            |
|`05_characterisation_audit.ipynb`|EDS composition sums, observed peak positions against tabulated emission lines, surface loading against the AAS uptake, which FTIR shifts exceed the instrument resolution, and whether the SEM field-width annotations are internally consistent                                              |

`src/adsorption.py` holds the model equations, the AICc calculation and the
fitting wrapper. The wrapper reports non-convergence instead of returning a
plausible-looking number, which is how the Mn(II) Langmuir failure surfaces.

## Data

`data/aas_readings.csv` is the laboratory report as reported: one row per
reading, with its dilution factor, the experimental conditions, and a note
where the laboratory record is ambiguous. Nothing is cleaned or dropped in the
file itself. Exclusions happen in the notebooks, where the reason is visible.

`data/eds_composition.csv`, `data/eds_peaks.csv` and `data/ftir_bands.csv`
reproduce Tables 6 and 7 and the peak positions read from Figures 9 to 11.
`data/xray_lines.csv` holds the reference X-ray emission energies used to check
the instrument’s element labels, and `data/sem_parameters.csv` reproduces the
Appendix C micrograph annotations.

Known limits of the dataset, carried through the analysis rather than hidden:

- No replicate flasks. Variability is estimated from three experimental series
  that happen to share one condition (1.1% RSD for Mn, 9.2% for Zn).
- Measured stock and control concentrations did not match their nominal values,
  so absolute uptake values are provisional. Removal percentages are not
  affected, because sample and control share a dilution.
- The dilution convention (`concentration = reading × dilution factor`) is
  inferred from the report layout; the laboratory’s working notes were not
  available.
- The temperature series was thermostatted for only the first 20 min of a
  140-min run.
- Some laboratory labels are undocumented (duplicate 120-min and pH 7 readings,
  a sixth temperature reading). Both alternatives are tested where they could
  change a conclusion.

## Running it

```bash
pip install -r requirements.txt
jupyter lab          # or: jupyter notebook
```

Run the notebooks in order from the repository root; they read `data/` and
`src/` by relative path. Outputs are committed so the notebooks render on
GitHub without being executed.

## Licence

Code released under the MIT Licence. The underlying experimental data is from
my own B.Sc. research project.
