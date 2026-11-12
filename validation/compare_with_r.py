"""Compare PyQUEFTS with the R-QUEFTS reference outputs produced by generate_r_reference.R.

Usage:  python compare_with_r.py [reference_dir] [--report report.md]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pyquefts as pq  # noqa: E402
from pyquefts.core.result import R_OUTPUTS  # noqa: E402

RTOL, ATOL = 1e-9, 1e-7


def compare_model(ref: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    crops = {name: pq.get_crop(name) for name in ref["crop"].unique()}
    # run per crop (parameters differ) but vectorized within crop
    py = pd.DataFrame(index=ref.index, columns=list(R_OUTPUTS), dtype=float)
    for name, g in ref.groupby("crop", sort=False):
        soil = pq.SoilParameters(N_base_supply=g.N_base_supply.values, P_base_supply=g.P_base_supply.values,
                                 K_base_supply=g.K_base_supply.values, N_recovery=g.N_recovery.values,
                                 P_recovery=g.P_recovery.values, K_recovery=g.K_recovery.values)
        fert = pq.Fertilizer(N=g.N.values, P=g.P.values, K=g.K.values)
        biom = pq.Biomass(leaf_att=g.leaf_att.values, stem_att=g.stem_att.values, store_att=g.store_att.values,
                          SeasonLength=g.SeasonLength.values)
        res = pq.run_quefts(soil, crops[name], fert, biom, keep_intermediate=False)
        for k in R_OUTPUTS:
            py.loc[g.index, k] = res[k]
    rows = []
    for k in R_OUTPUTS:
        r, p = ref[k].to_numpy(float), py[k].to_numpy(float)
        both_nan = np.isnan(r) & np.isnan(p)
        d = np.abs(r - p)
        rel = d / np.maximum(np.abs(r), 1e-12)
        ok = both_nan | np.isclose(r, p, rtol=RTOL, atol=ATOL)
        rows.append({"variable": k, "n": len(r), "n_mismatch": int((~ok).sum()),
                     "max_abs_diff": float(np.nanmax(d)) if np.isfinite(d).any() else np.nan,
                     "max_rel_diff": float(np.nanmax(rel)) if np.isfinite(rel).any() else np.nan,
                     "n_nan_r": int(np.isnan(r).sum()), "n_nan_py": int(np.isnan(p).sum())})
    summary = pd.DataFrame(rows)
    mism = ref.copy()
    for k in R_OUTPUTS:
        mism[f"py_{k}"] = py[k]
    bad = ~np.all([np.isclose(ref[k], py[k], rtol=RTOL, atol=ATOL) | (np.isnan(ref[k]) & np.isnan(py[k]))
                   for k in R_OUTPUTS], axis=0)
    return summary, mism[bad]


def compare_supply(ref: pd.DataFrame) -> pd.DataFrame:
    s1 = pq.nut_supply1(ref.pH, ref.SOC, ref.Kex, ref.Polsen)
    s1t = pq.nut_supply1(ref.pH, ref.SOC, ref.Kex, ref.Polsen, ptotal=ref.Ptotal)
    s2 = pq.nut_supply2(ref.temp, ref.pH, ref.SOC, ref.Kex, ref.Polsen, ptotal=ref.Ptotal)
    rows = []
    for prefix, sup in (("s1_", s1), ("s1tot_", s1t), ("s2_", s2)):
        for nut in ("N", "P", "K"):
            col = f"{prefix}{nut}_base_supply"
            r, p = ref[col].to_numpy(float), np.asarray(getattr(sup, f"{nut}_base_supply"), float)
            d = np.abs(r - p)
            rows.append({"variable": col, "n": len(r), "n_mismatch": int((~np.isclose(r, p, rtol=RTOL, atol=ATOL)).sum()),
                         "max_abs_diff": float(d.max()), "max_rel_diff": float((d / np.maximum(np.abs(r), 1e-12)).max())})
    return pd.DataFrame(rows)


def compare_examples(ref: pd.DataFrame) -> pd.DataFrame:
    cases = {}
    cases["default"] = pq.Quefts().run()
    cases["barley_N100"] = pq.Quefts(crop="Barley", fert=dict(N=100, P=0, K=0),
                                     biom=dict(leaf_att=2200, stem_att=2700, store_att=4800, SeasonLength=110)).run()
    cases["barley_long"] = pq.Quefts(crop="Barley", fert=dict(N=100, P=50, K=50),
                                     biom=dict(leaf_att=2651, stem_att=5053, store_att=8208, SeasonLength=162)).run()
    for N in (1, 51, 101, 191):
        cases[f"rice_N{N}"] = pq.Quefts(crop="Rice", fert=dict(N=N, P=0, K=0),
                                        biom=dict(leaf_att=2651, stem_att=5053, store_att=8208)).run()
    rows = []
    for _, r in ref.iterrows():
        py = cases[r["case"]]
        for k in R_OUTPUTS:
            rows.append({"case": r["case"], "variable": k, "r": r[k], "py": py.scalar(k),
                         "abs_diff": abs(r[k] - py.scalar(k))})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("refdir", nargs="?", default=str(Path(__file__).parent))
    ap.add_argument("--report", default=None)
    a = ap.parse_args()
    refdir = Path(a.refdir)

    model_ref = pd.read_csv(refdir / "r_reference_model.csv", keep_default_na=True)
    model_ref["crop"] = model_ref["crop"].fillna("")
    summ, bad = compare_model(model_ref)
    sup_ref = pd.read_csv(refdir / "r_reference_supply.csv")
    sup_summ = compare_supply(sup_ref)
    ex_ref = pd.read_csv(refdir / "r_reference_examples.csv")
    ex = compare_examples(ex_ref)
    info = (refdir / "r_session_info.txt").read_text().strip().splitlines()

    pd.set_option("display.width", 200)
    print("== model runs ==")
    print(summ.to_string(index=False))
    print(f"mismatching scenarios: {len(bad)} of {len(model_ref)}")
    print("== supply ==")
    print(sup_summ.to_string(index=False))
    print("== documented examples ==")
    print(f"max abs diff: {ex['abs_diff'].max():.3e}")

    if a.report:
        lines = ["# PyQUEFTS vs R-QUEFTS validation report", "",
                 f"Reference: {'; '.join(info)}", f"Tolerance: rtol={RTOL}, atol={ATOL}", "",
                 f"## Model runs ({len(model_ref)} randomized scenarios, {model_ref['crop'].nunique()} crop parameter sets)", "",
                 summ.to_markdown(index=False), "", f"Scenarios with any mismatch: **{len(bad)}**", "",
                 "Scenario design: random base supplies (incl. zeros and extremes), recoveries, fertilizer rates "
                 "(30 % zeros), attainable yields 500-15000 kg/ha with 8 % tiny storage organs (exercising the "
                 "leaf/stem organ branches), season lengths 60-400 days, all 28 crops plus the default set.", "",
                 f"## Supply functions ({len(sup_ref)} random soils, pH 3.5-9.5, SOC 0-80, Kex 0-40, P-Olsen 0-40, Ptotal 0-2000, T 10-35)", "",
                 sup_summ.to_markdown(index=False), "",
                 "## Documented examples (reagro.org tutorial / Rquefts manual)", "",
                 ex.pivot(index="case", columns="variable", values="abs_diff").max(axis=1).rename("max_abs_diff").to_frame().to_markdown(), ""]
        if len(bad):
            lines += ["## Mismatching scenarios (first 20)", "", bad.head(20).to_markdown(index=False)]
        Path(a.report).write_text("\n".join(lines), encoding="utf-8")
        print("report written to", a.report)
    return 0 if len(bad) == 0 and sup_summ.n_mismatch.sum() == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
