# Where western trees grow vs. where they thrive

Do tree species occur most often in the places where they actually grow, survive, and reproduce best? This project uses U.S. Forest Service inventory data from California, Oregon, and Washington to compare where ~60 western tree species are **found** with where they **perform** best along a soil-water-balance (drought) gradient, and asks whether species' functional traits explain the gaps.

**Finding:** distributions and performance were often mismatched, and functional traits helped explain which species showed larger mismatches.

Published as: Bohner, T. & Diez, J. (2020). Extensive mismatches between species distributions and performance depend on functional traits. *Ecology Letters* 23(1): 33–44. [doi:10.1111/ele.13396](https://doi.org/10.1111/ele.13396)

## Where to start

If you're here to see how the modeling works, these three files are the core:

- [`Stan Models/surv_swb_temp_spp.stan`](Stan%20Models/surv_swb_temp_spp.stan): hierarchical survival model with plot-level intercepts; `generated quantities` computes each species' optimum along the drought gradient
- [`R scripts/4_All Vitals Stan Models.R`](R%20scripts/4_All%20Vitals%20Stan%20Models.R): fits the models species by species, checks convergence, and derives mismatch and demographic cost
- [`R scripts/5_Trait phylo brms.R`](R%20scripts/5_Trait%20phylo%20brms.R): phylogenetic regressions of mismatch and cost on traits, with LOO model selection

## At a glance

| | |
|---|---|
| **Data** | USDA Forest Service Forest Inventory and Analysis (FIA) plots in CA, OR, WA; survival records for ~580,000 trees across ~60 species; climate rasters; trait data from five sources |
| **Models** | Bayesian hierarchical models in Stan, one per response: occurrence (Bernoulli), abundance (negative binomial), growth (normal), survival (Bernoulli), recruitment (zero-inflated negative binomial) |
| **Key derived metrics** | For each species and vital rate: the drought-gradient optimum for performance vs. for occurrence ("mismatch"), and the performance lost at the distribution optimum ("demographic cost") |
| **Trait analysis** | Phylogenetic regressions in brms relating mismatch and cost to traits, with LOO model comparison for phylogenetic signal |
| **Stack** | R (tidyverse, data.table, rstan, brms, loo, bayesplot, tidybayes, ape, BIEN, taxize) and Stan |

## Workflow

The full annotated pipeline is in [`FIA workflow.Rmd`](FIA%20workflow.Rmd) (rendered: [`FIA-workflow.html`](FIA-workflow.html)). Scripts are numbered in run order:

1. **Data gathering & cleaning**
   - `1_Raw FIA Processing.R`: joins tree, seedling, and plot tables; attaches climate values; removes problematic records; builds occurrence and abundance datasets
   - `1_Gather All Traits.R`: merges trait data (BIEN, TRY, USDA PLANTS, World Agroforestry, Niinemets & Valladares 2006) and builds a phylogeny from NCBI
2. **Data prep**, `2_Data Prep.R`: keeps species found on more than 25 plots, uses the most recent survey for occurrence and abundance and repeat-visit plots for growth, and builds the species index table Stan needs (`spptable.csv`)
3. **Exploratory analysis**: vital-rate correlations, range coverage against Little's range maps, trait distributions
4. **Modeling**, `4_All Vitals Stan Models.R`: fits each response model, checks diagnostics, and computes optima, mismatch, and demographic cost
5. **Summaries & figures**: coefficient and mismatch summaries, phylogenetic trait regressions (`5_Trait phylo brms.R`), and final manuscript figures
6. **Sensitivity check**, `6_Pres abs with temp comparison.R`: compares occurrence models with soil water balance alone vs. soil water balance plus temperature

Stan model files are in [`Stan Models/`](Stan%20Models).

## Data

Data are not included in this repository (~8.6 GB) but are publicly available:

- **FIA:** state `TREE`, `SEEDLING`, `PLOT`, `TREE_GRM_ESTN`, and `REF_SPECIES` tables for CA, OR, and WA from the USDA Forest Service FIA DataMart
- **Soil water balance:** CGIAR-CSI Global High-Resolution Soil-Water Balance (Priestley–Taylor α, 30 arc-second)
- **Temperature:** [WorldClim v2](https://www.worldclim.org/) (30 arc-second)
- **Traits:** BIEN, TRY (request 4612), USDA PLANTS, World Agroforestry wood density, Niinemets & Valladares (2006)
- **Range maps:** Little's tree range maps (USTreeAtlas digitized shapefiles)

Plot-level climate values were matched to true (unfuzzed) plot coordinates by the FIA office. The two small species reference tables in `Processed Data/Analysis Ready Data/` show the species analyzed and the sample sizes per model.

## Notes

- This code dates from 2020 (R 3.6 era). Some spatial packages it uses (`rgdal`, `rgeos`) have since been retired from CRAN; `sf` and `terra` are the modern replacements.
- Figure output paths and saved model fits were not fully reorganized when the project was archived (see the To-Do note in the workflow file).
