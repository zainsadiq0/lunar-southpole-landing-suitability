# Lunar South Pole Landing Site Suitability

**A multi-criteria GIS suitability and Engineering analysis of five Artemis III candidate landing regions**

By **Zain Sadiq** producing geospatial data analysis, and **Yash Agrawal** with engineering analysis.

**NOT AFFILIATED OR SUPPORTED BY NASA**

> Independent academic projectAll inputs are public NASA LRO data products.

![Overview of the five candidate landing regions near the lunar south pole](images/01_overview_map.png)

## Summary

Site07 (Peak Near Shackleton) ranks first of five candidate sites, with a mean suitability score of 5.44 on a 1 to 10 scale. Site20 (Leibnitz Beta Plateau) is second at 5.06. The ranking is a weighted overlay of four criteria (slope, illumination, a water-ice proxy and proximity to permanently shadowed regions) built in ArcGIS Pro from LRO LOLA-derived data.

![Baseline suitability ranking of the five sites](images/04_baseline_ranking.png)

| Rank | Site | Name | Mean suitability (1 to 10) |
| --- | --- | --- | --- |
| 1 | Site07 | Peak Near Shackleton | 5.44 |
| 2 | Site20 | Leibnitz Beta Plateau | 5.06 |
| 3 | Site11 | de Gerlache Rim | 4.86 |
| 4 | Site01 | Connecting Ridge | 4.75 |
| 5 | Site23 | Malapert Massif | 4.48 |

The spread between sites is small (under one point from first to last), so this is a screening result rather than a landing recommendation.

## Data

| Input | Resolution | Used for |
| --- | --- | --- |
| LOLA DEM (LRO) | 5 m | Elevation and terrain profiles |
| Slope derived from the DEM | 5 m | Slope criterion |
| Illumination | 120 m | Illumination criterion |
| Permanently shadowed regions (PSR) | 120 m | PSR proximity and water-ice proxy |

All layers use a South Polar Stereographic Moon projection. Data come from NASA GSFC LOLA/LRO products distributed through the PGDA site products for the five candidate regions.

## Methodology

Each criterion raster was reclassified to a common 1 to 10 scale (higher is better), multiplied by its weight and summed into one suitability raster per site. A site's baseline score is the mean of its suitability raster.

| Criterion | Weight |
| --- | --- |
| Slope | 40% |
| Illumination | 25% |
| Water-ice proxy | 20% |
| PSR proximity | 15% |

Slope carries the most weight because landing safety is the first constraint. The water-ice criterion is a PSR-based proxy, because the LEND neutron-detector raster was unusable for this analysis.

Tools: ArcGIS Pro (Raster Calculator, Reclassify, Zonal Statistics as Table, Stack Profile, layouts). Suitability maps use one fixed color stretch (1.5 to 8) so sites compare directly, with red for low scores and green for high.

![Baseline suitability maps for all five sites on one color scale](images/02_baseline_suitability.png)

## Site07 in detail

Site07 covers 255.68 km² centered near 88.81°S, 123.69°E, and 45.2% of its area has a slope under 10°.

| Measure | Value |
| --- | --- |
| Elevation, min / mean / max | −418.9 m / 800.6 m / 1,648.7 m |
| Slope, mean / median | 10.74° / 10.61° |
| Slope, 90th percentile / max | 16.85° / 60.48° |
| Area with slope < 5° | 11.7% |
| Area with slope < 10° | 45.2% |

Its 5.44 is carried by PSR proximity and slope, while illumination and the water-ice proxy pull it down. Each score is the zonal mean of that criterion's 1 to 10 raster over the site.

| Criterion | Weight | Mean score (1 to 10) | Weighted contribution |
| --- | --- | --- | --- |
| Slope | 40% | 6.92 | 2.77 |
| Illumination | 25% | 3.94 | 0.98 |
| Water-ice proxy | 20% | 2.18 | 0.44 |
| PSR proximity | 15% | 8.31 | 1.25 |
| **Total** | 100% | | **5.44** |

![Site07 criterion scores](images/03_site07_criterion_scores.png)

## Engineering feasibility screen (delta-v)

Yash Agrawal's MATLAB simulation estimates that landing at Site07's coordinates from a 100 km circular lunar orbit takes about **2,046 m/s** of propulsive delta-v and about **112,300 kg** of propellant, for an assumed 250,000 kg Starship HLS-inspired vehicle.

| Quantity | Result |
| --- | --- |
| Deorbit impulse | 19.46 m/s (about 1,413 kg of propellant) |
| Powered descent | about 2,026.6 m/s (110,850 kg of propellant) |
| Total propulsive delta-v | about 2,046.1 m/s |
| Total propellant | about 112,263 kg |
| Touchdown speed | 1.80 m/s (limit 2.0 m/s) |
| Landing error | 3.7 m (limit 100 m) |
| Mass above assumed dry mass at touchdown | about 17,737 kg |

**Method.** An impulsive deorbit from the 100 km orbit to a 15 km perilune, a two-body coast, then a finite-thrust, feedback-guided descent integrated with MATLAB ode45. A grid search of 80 combinations of descent start point and guidance time found 72 that met the constraints (touchdown at or below 2 m/s, error at or below 100 m, final mass at or above 125,000 kg). The reported case is the lowest-propellant one with touchdown at or below 1.85 m/s. In 60 perturbed runs (ignition state, mass, thrust and specific impulse), all met the constraints; the perturbation sizes are illustrative, so this is not a reliability estimate.

**Assumptions.** The vehicle inputs are illustrative, not verified HLS specifications: 250,000 kg initial mass, 120,000 kg dry mass, 2 MN maximum thrust, 350 s specific impulse, 10% minimum throttle. Under them the landing closes with about 17.7 t of propellant above dry mass, roughly 14% of the assumed usable propellant. The result depends directly on these numbers, so it is a first-order screen, not a statement about what HLS can do.

**What it does not use.** The model uses only Site07's latitude and longitude. It assumes a spherical, non-rotating Moon and a fixed target, and it does not use this project's slope, DEM or profile data. It leaves out hazards, navigation errors and attitude dynamics. Totals exclude Earth to Moon transfer, lunar orbit insertion, ascent and operational reserves.

## Limitations and next steps

- The water-ice criterion is a PSR-based proxy, not a direct measurement.
- Slope is at 5 m resolution while illumination and PSR are at 120 m, so those criteria are coarser.
- The weights (40/25/20/15) are a judgment call and the site scores are close together.
- The engineering screen covers delta-v only, uses assumed vehicle numbers, and is not coupled to the GIS scores.
- Next: test how sensitive the ranking is to alternative weightings, and extend the engineering screen to Site20 and Site11.

## Author

Geospatial Analysis: Zain Sadiq - Rutgers University 
Aerospace Engineering analysis: Yash Agrawal - Rutgers University
