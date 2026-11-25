# Sources for scientific choices

These are method and format sources, not evidence that a particular local run
has been validated. Cite the datasets, versions, models and analysis choices
actually used in each study. Original links checked 2026-11-19; CF and the added
xarray input-decoding/resource-management references checked 2026-11-25.

| Topic | Primary or official source | Use in this skill |
|---|---|---|
| CMIP6 experiment design | Eyring et al. (2016), [GMD 9, 1937–1958](https://gmd.copernicus.org/articles/9/1937/2016/), DOI 10.5194/gmd-9-1937-2016 | Historical/scenario experiment context. |
| Simulation identifiers | WCRP CMIP, [CMIP Data Access](https://www.wcrp-cmip.org/cmip-data-access/) and [CMIP6 controlled vocabularies](https://github.com/WCRP-CMIP/CMIP6_CVs) | Distinguish source, experiment, variant, variable/table and grid label. |
| Calendar, time bounds and quantity semantics | CF Community (2024), [CF Conventions 1.12](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html) | Inspect calendar, units, bounds and cell methods before averaging or integrating. |
| Input decoding and bounded reads | xarray, [open_dataset](https://docs.xarray.dev/en/latest/generated/xarray.open_dataset.html) and [Reading and writing files](https://docs.xarray.dev/en/latest/user-guide/io.html) | Decode CF fill values/scales once, subset before loading, and manage dataset lifetimes explicitly. |
| Monthly seasonal means | xarray, [Calculating Seasonal Averages from Monthly Means](https://docs.xarray.dev/en/v2025.08.0/examples/monthly-means.html) | Weight monthly means by each calendar's month lengths. |
| Quantile delta mapping | Cannon, Sobie & Murdock (2015), [Journal of Climate](https://doi.org/10.1175/JCLI-D-14-00754.1); [author repository copy](https://dspace.library.uvic.ca/bitstreams/f26cd737-49bb-42d4-b28d-f44987f6d161/download) | Preserve simulated quantile changes; distinguish this from guaranteed temporal-trend preservation. |
| GEV parameter convention | SciPy, [genextreme](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genextreme.html) | SciPy c = −ξ; use stable survival/inverse-survival operations. |
| GPD parameter convention | SciPy, [genpareto](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genpareto.html) | SciPy GPD c = ξ and the excess support. |
| Fitted-distribution goodness of fit | SciPy, [goodness_of_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.goodness_of_fit.html) | Monte Carlo calibration must re-fit unknown parameters for each sample. |
| Projection and ensemble interpretation | IPCC AR6 WGI (2021), [Chapter 10](https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-10/) | Ensemble spread is not a complete uncertainty measure; adjustment cannot repair all process biases. |

Use exact versions of software and data in the run manifest. Consult current
upstream documentation if an API, dataset endpoint or convention has changed.
