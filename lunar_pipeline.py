#!/usr/bin/env python3
"""
Lunar South Pole Landing Site Suitability: reproducible pipeline.

Automates the quantitative steps of the project outside ArcGIS Pro:

  1. Weighted overlay of four reclassified (1-10) criterion rasters per site
  2. Zonal mean of every criterion and of the overlay (one number per site)
  3. Ranking table (CSV + Markdown) and the baseline ranking bar chart
  4. Delta-v / propellant screen: impulsive deorbit from the rocket equation,
     combined with the powered-descent propellant from Yash Agrawal's MATLAB
     model (hls_integrated_descent_v17.m)
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

CRITERIA = ["slope", "illumination", "water_ice", "psr"]

EXAMPLE_CONFIG = {
    "output_dir": "outputs",
    "weights": {"slope": 0.40, "illumination": 0.25, "water_ice": 0.20, "psr": 0.15},
    "highlight_site": "Site07",
    "sites": {
        "Site07": {
            "name": "Peak Near Shackleton",
            "rasters": {
                "slope": "rasters/site07_slope_1to10.tif",
                "illumination": "rasters/site07_illum_1to10.tif",
                "water_ice": "rasters/site07_water_1to10.tif",
                "psr": "rasters/site07_psr_1to10.tif",
            },
        },
        "Site20": {"name": "Leibnitz Beta Plateau", "rasters": {}},
        "Site11": {"name": "de Gerlache Rim", "rasters": {}},
        "Site01": {"name": "Connecting Ridge", "rasters": {}},
        "Site23": {"name": "Malapert Massif", "rasters": {}},
    },
    "_comment_sites": (
        "Give each site a 'rasters' block like Site07's (already 1-10), "
        "or a 'criterion_means' block {slope, illumination, water_ice, psr} "
        "if you already exported the zonal means from ArcGIS."
    ),
    "delta_v": {
        "mu_km3_s2": 4902.800066,
        "radius_km": 1737.4,
        "orbit_altitude_km": 100.0,
        "perilune_altitude_km": 15.0,
        "initial_mass_kg": 250000.0,
        "dry_mass_kg": 120000.0,
        "isp_s": 350.0,
        "powered_descent_propellant_kg": 110850.0,
        "_comment": "Powered-descent propellant comes from Yash's V17 MATLAB run.",
    },
}

G0 = 9.80665  # m/s^2


# ---------------------------------------------------------------- raster step
def reclassify(array, bins, scores):
    """Map a raw array to 1-10 scores. bins are ascending edges; scores has
    len(bins)+1 entries. NaN stays NaN. You supply the breakpoints."""
    out = np.full(array.shape, np.nan, dtype="float64")
    valid = ~np.isnan(array)
    out[valid] = np.asarray(scores, dtype="float64")[np.digitize(array[valid], bins)]
    return out


def read_raster(path):
    import rasterio  # imported here so the non-raster steps run without it

    with rasterio.open(path) as src:
        arr = src.read(1).astype("float64")
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
    return arr


def site_from_rasters(rasters, weights):
    """Return (criterion_means, overlay_mean) for one site.

    The overlay is computed only where all four criteria are valid, so its
    mean equals sum(weight * criterion mean over that same mask)."""
    stack = np.stack([read_raster(rasters[c]) for c in CRITERIA])
    if len({a.shape for a in stack}) != 1:
        raise ValueError("Criterion rasters for a site must share one grid.")
    mask = ~np.isnan(stack).any(axis=0)
    if not mask.any():
        raise ValueError("No valid pixels common to all four rasters.")
    w = np.array([weights[c] for c in CRITERIA]).reshape(-1, 1, 1)
    overlay = (stack * w).sum(axis=0)
    means = {c: float(stack[i][mask].mean()) for i, c in enumerate(CRITERIA)}
    return means, float(overlay[mask].mean())


# ------------------------------------------------------------------- ranking
def score_sites(cfg):
    weights = cfg["weights"]
    if not math.isclose(sum(weights.values()), 1.0, abs_tol=1e-9):
        raise ValueError(f"Weights must sum to 1, got {sum(weights.values())}")
    rows = []
    for sid, s in cfg["sites"].items():
        if s.get("rasters"):
            means, overlay = site_from_rasters(s["rasters"], weights)
        elif s.get("criterion_means"):
            means = s["criterion_means"]
            overlay = sum(weights[c] * means[c] for c in CRITERIA)
        else:
            print(f"  skipping {sid}: no rasters or criterion_means", file=sys.stderr)
            continue
        # Linearity check: weighted sum of criterion means == overlay mean.
        check = sum(weights[c] * means[c] for c in CRITERIA)
        if abs(check - overlay) > 1e-6:
            print(f"  WARNING {sid}: overlay mean {overlay:.4f} != "
                  f"weighted sum of means {check:.4f}", file=sys.stderr)
        rows.append({"site": sid, "name": s.get("name", ""), "mean": overlay,
                     **{f"mean_{c}": means[c] for c in CRITERIA}})
    rows.sort(key=lambda r: r["mean"], reverse=True)
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def write_ranking(rows, weights, outdir):
    fields = ["rank", "site", "name", "mean"] + [f"mean_{c}" for c in CRITERIA]
    with open(outdir / "ranking.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in r.items() if k in fields})
    md = ["| Rank | Site | Name | Mean suitability (1 to 10) |", "| --- | --- | --- | --- |"]
    md += [f"| {r['rank']} | {r['site']} | {r['name']} | {r['mean']:.2f} |" for r in rows]
    (outdir / "ranking_table.md").write_text("\n".join(md) + "\n")

    top = rows[0]
    labels = {"slope": "Slope", "illumination": "Illumination",
              "water_ice": "Water-ice proxy", "psr": "PSR proximity"}
    cm = ["| Criterion | Weight | Mean score (1 to 10) | Weighted contribution |",
          "| --- | --- | --- | --- |"]
    for c in CRITERIA:
        cm.append(f"| {labels[c]} | {weights[c]:.0%} | {top['mean_'+c]:.2f} | "
                  f"{weights[c]*top['mean_'+c]:.2f} |")
    cm.append(f"| **Total** | 100% | | **{top['mean']:.2f}** |")
    (outdir / f"{top['site'].lower()}_criterion_table.md").write_text("\n".join(cm) + "\n")


def plot_ranking(rows, highlight, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    accent, grey, ink = "#2a78d6", "#b9b8b2", "#1A2A33"
    labels = [f"{r['site']}\n{r['name']}" for r in rows]
    vals = [r["mean"] for r in rows]
    colors = [accent if r["site"] == highlight else grey for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    bars = ax.bar(labels, vals, color=colors, width=0.62)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.08, f"{v:.2f}",
                ha="center", va="bottom", fontsize=11, color=ink)
    ax.set_ylim(0, 10)
    ax.set_ylabel("Mean suitability (1-10)", color=ink)
    ax.set_title(f"{rows[0]['site']} ranks first of {len(rows)} candidate sites "
                 f"with a mean suitability of {rows[0]['mean']:.2f}",
                 fontsize=12, color=ink, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=8.5)
    fig.tight_layout()
    fig.savefig(outdir / "baseline_ranking.png", dpi=300)
    plt.close(fig)


# ------------------------------------------------------------------- delta-v
def delta_v_screen(p):
    mu, R = p["mu_km3_s2"], p["radius_km"]
    r1 = R + p["orbit_altitude_km"]
    r2 = R + p["perilune_altitude_km"]
    a = (r1 + r2) / 2
    v_circ = math.sqrt(mu / r1)
    v_apo = math.sqrt(mu * (2 / r1 - 1 / a))  # vis-viva at the 100 km end
    dv_deorbit = (v_circ - v_apo) * 1000.0     # m/s

    ve = p["isp_s"] * G0
    m0 = p["initial_mass_kg"]
    prop_deorbit = m0 * (1 - math.exp(-dv_deorbit / ve))
    m1 = m0 - prop_deorbit
    prop_pd = p["powered_descent_propellant_kg"]
    mf = m1 - prop_pd
    dv_pd = ve * math.log(m1 / mf)
    return {
        "deorbit_dv_m_s": dv_deorbit,
        "deorbit_propellant_kg": prop_deorbit,
        "mass_after_deorbit_kg": m1,
        "powered_descent_propellant_kg": prop_pd,
        "powered_descent_dv_m_s": dv_pd,
        "total_dv_m_s": dv_deorbit + dv_pd,
        "total_propellant_kg": prop_deorbit + prop_pd,
        "touchdown_mass_kg": mf,
        "margin_above_dry_mass_kg": mf - p["dry_mass_kg"],
        "note": "Deorbit burn from the rocket equation; powered-descent propellant "
                "is taken from the MATLAB V17 run (config). Illustrative vehicle "
                "assumptions, not verified HLS specifications.",
    }


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", nargs="?", help="path to config.json")
    ap.add_argument("--write-example-config", action="store_true")
    args = ap.parse_args()

    if args.write_example_config:
        print(json.dumps(EXAMPLE_CONFIG, indent=2))
        return
    if not args.config:
        ap.error("config.json required (or use --write-example-config)")

    cfg = json.loads(Path(args.config).read_text())
    base = Path(args.config).resolve().parent
    # resolve raster paths relative to the config file
    for s in cfg["sites"].values():
        s["rasters"] = {k: str((base / v)) for k, v in s.get("rasters", {}).items()}
    outdir = base / cfg.get("output_dir", "outputs")
    outdir.mkdir(parents=True, exist_ok=True)

    print("1/3 Scoring sites ...")
    rows = score_sites(cfg)
    if not rows:
        sys.exit("No sites could be scored.")
    write_ranking(rows, cfg["weights"], outdir)
    for r in rows:
        print(f"   {r['rank']}. {r['site']:7s} {r['name']:24s} {r['mean']:.2f}")

    print("2/3 Plotting ranking chart ...")
    plot_ranking(rows, cfg.get("highlight_site", rows[0]["site"]), outdir)

    print("3/3 Delta-v screen ...")
    dv = delta_v_screen(cfg["delta_v"])
    (outdir / "delta_v_screen.json").write_text(json.dumps(dv, indent=2))
    print(f"   deorbit {dv['deorbit_dv_m_s']:.3f} m/s, "
          f"total {dv['total_dv_m_s']:.1f} m/s, "
          f"propellant {dv['total_propellant_kg']:.0f} kg, "
          f"margin above dry {dv['margin_above_dry_mass_kg']:.0f} kg")
    print(f"Done. Outputs in {outdir}")


if __name__ == "__main__":
    main()
