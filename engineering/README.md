# Engineering: HLS-Inspired Lunar Descent (Delta-v Screen for Site07)

Analysis by **Yash Agrawal**. MATLAB simulation of an orbit-to-surface landing at Site07 (Peak Near Shackleton) at latitude −88.811008°, longitude 123.690068°, starting from a 100 km circular lunar orbit.

> **Read this first.** All vehicle inputs are illustrative assumptions, not verified Starship HLS specifications. This is an independent academic model, not affiliated with or endorsed by NASA or SpaceX. The full write-up is in [`HLS_Descent_Technical_Handoff.pdf`](HLS_Descent_Technical_Handoff.pdf).

## Result (simulation V17)

Estimated propulsive delta-v and propellant for the deorbit burn plus the powered descent:

| Quantity | V17 result |
| --- | --- |
| Initial mass before deorbit | 250,000 kg |
| Deorbit impulse | 19.455 m/s |
| Deorbit propellant | about 1,413 kg |
| Mass after deorbit | about 248,587 kg |
| Powered-descent propellant | 110,850 kg |
| Powered-descent delta-v | about 2,026.6 m/s |
| **Total propulsive delta-v** | **about 2,046.1 m/s** |
| **Total propellant** | **about 112,263 kg** |
| Touchdown mass | about 137,737 kg |
| Mass above assumed dry mass | about 17,737 kg |
| Selected coast fraction / guidance time | 0.928 / 660 s |
| Powered-descent duration | 665.02 s |
| Touchdown speed | 1.8028 m/s |
| Landing-site error | 3.6949 m |

Delta-v here is the accumulated propulsive increment, the integral of thrust divided by mass over time, not the difference between the initial and final velocity vectors. Totals cover deorbit and powered descent only. They exclude Earth to Moon transfer, lunar orbit insertion, ascent and operational reserves.

## Assumed parameters

| Parameter | Value |
| --- | --- |
| Lunar gravitational parameter, μ | 4.902800066 × 10¹² m³/s² |
| Lunar radius, R (spherical) | 1,737,400 m |
| Initial circular orbit altitude | 100 km |
| Transfer perilune altitude | 15 km |
| Initial mass / dry mass | 250,000 kg / 120,000 kg |
| Maximum thrust | 2.0 × 10⁶ N |
| Specific impulse | 350 s |
| Minimum commanded ON throttle | 10% of maximum |
| Maximum thrust-magnitude rate | 20,000 N/s |
| Thrust-direction slew limit (illustrative) | 2°/s |
| Touchdown speed limit / preferred target | 2.0 / 1.85 m/s |
| Maximum landing-site error | 100 m |
| Minimum final mass | 125,000 kg |

## Method

1. **Impulsive deorbit.** A retrograde impulse moves the vehicle from the 100 km circular orbit onto an ellipse with a 15 km perilune, aligned so perilune is over Site07. Propellant comes from the ideal rocket equation.
2. **Ballistic coast.** A two-body, spherical-Moon propagation (MATLAB `ode45`) up to a chosen fraction of the half-period of the transfer ellipse.
3. **Finite-thrust powered descent.** Three-dimensional equations of motion with mass depletion, driven by a feedback guidance law. Thrust is limited in magnitude and rate, and the thrust direction slews toward the commanded direction at a limited rate.
4. **Search.** V17 tested 80 combinations of coast fraction and guidance time. A feasible solution needs surface contact, touchdown at or below 2 m/s, site error at or below 100 m, and final mass at or above 125,000 kg. Among feasible solutions with touchdown at or below 1.85 m/s, the lowest-propellant case was selected. This is a discrete grid-search optimum, not a proof of global optimality.

## Validation and limitations

- **Feasibility:** 72 of the 80 candidates met the basic constraints. 58 also met the 1.85 m/s preferred speed target.
- **Solver convergence:** maximum integration steps of 2, 1 and 0.5 s gave identical displayed touchdown metrics. This supports numerical consistency at the reported precision. It does not validate the physical assumptions.
- **Perturbation study:** all 60 illustrative trials met the modeled constraints. The worst touchdown speed was 1.8054 m/s, the worst miss distance 3.696 m, and the minimum margin above dry mass 16,760.7 kg. The perturbations were Gaussian, 1 sigma: 10 m per position component, 0.1 m/s per velocity component, 100 kg of ignition mass, 1% thrust scaling and 0.5% specific-impulse scaling. These are illustrative distributions, so 60/60 successes are **not** a mission reliability estimate.
- **Not modeled:** lunar rotation, gravity harmonics, terrain and hazards, navigation errors, six-degree-of-freedom attitude dynamics, true gimbal and ignition logic, communications, structural loads, and verified Starship HLS performance data. The model uses a spherical, non-rotating Moon and a fixed inertial target.
- **Not linked to the GIS analysis:** V17 uses only Site07's latitude and longitude. It does not use the slope, DEM or terrain-profile data from the suitability analysis.

## Files

All files are MATLAB (`.m`) unless noted. The numerical results above come from **V17**.

| Role | File | Purpose |
| --- | --- | --- |
| **Primary** | `hls_integrated_descent_v17.m` | Final trajectory, fuel optimization, convergence check and robustness study |
| Support | `hls_integrated_descent_v16.m` | Thrust-direction slew diagnostics |
| Support | `hls_v15_validation.m` | Numerical validation study. Depends on the V14 workspace, so run `hls_integrated_descent_v14.m` first |
| Support | `hls_integrated_descent_v14.m` | Earlier engine-constraint baseline |
| Reference | `HLS_Descent_Technical_Handoff.pdf` | Methods, equations, assumptions, results and recommended figures |
| Illustrative only | `hls_landing_site_comparison_v18.m`, `hls_v18_site_screening.csv` | Orbital-plane accessibility screening for example sites. The non-reference coordinates are illustrative, and it does **not** compute site-specific descent fuel. Do not use it to compare landing sites |
| Archive | all other files, V1 to V13 | Development history, listed below |

### Archive (development history, not final results)

These earlier scripts led to V17. They use some different assumptions (for example a 100,000 kg dry mass in the earliest ones), and some start the descent from assumed conditions rather than from the coast simulation. Do not quote numbers from them.

| Files | What they are |
| --- | --- |
| `lunar_descent_baseline.m`, `lunar_trajectory.m` | Early two-body transfer baselines |
| `hls_powered_descent.m`, `hls_powered_descent_v2.m`, `hls_descent_optimization.m` | First two-dimensional powered-descent models and a start-altitude sweep |
| `hls_3d_lunar_orbit.m`, `hls_3d_site_comparison.m` | Orbit-plane geometry. The site comparison uses hypothetical test locations |
| `hls_3d_targeted_descent.m`, `_v2` to `_v5` | Early three-dimensional targeted-descent attempts |
| `hls_orbit_to_surface_v6.m`, `hls_site_aligned_deorbit_v7.m` | Orbit-to-surface model and site-aligned deorbit geometry |
| `hls_integrated_descent_v8.m`, `_v9.m`, `hls_feasibility_v10.m`, `_v11.m`, `_v12.m`, `_v13.m` | Integrated descent, braking feasibility and guidance searches |

## Running the code

Open `hls_integrated_descent_v17.m` in MATLAB and run it. It is standalone, prints the results and diagnostics, and draws the figures. `hls_v15_validation.m` requires the V14 workspace (run V14 first).

To export a figure at publication resolution, select the figure window and run:

```matlab
exportgraphics(gcf,'figure.png','Resolution',300)
```

## Credit

Simulation, equations and write-up: Yash Agrawal. Site coordinates are from the Site07 analysis in this repository (see the top-level README).
