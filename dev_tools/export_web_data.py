"""Compute every number shown on the web page with ion_track_mc and write ../data.js.

Nothing on the page is simulated in the browser: the tracks, curves and tables
come from this script, which runs the real code. Re-run it after changing the
code or its data tables:

    cd dev_tools
    ../../ion-track-mc/.venv/bin/python export_web_data.py        # ~10 min (the beam)
    ../../ion-track-mc/.venv/bin/python export_web_data.py --quick  # 5-proton beam, for testing
    ../../ion-track-mc/.venv/bin/python export_web_data.py --tracks-only  # keep the beam of data.js

It also copies the tutorial PDFs from the repository into ../docs/ and makes their
cover thumbnails (pdftoppm).

The repository is located through REPO below (sibling folder ion-track-mc/).
"""
import copy
import json
import sys
import time
import shutil
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APP = HERE.parent
REPO = APP.parent / "ion-track-mc"
sys.path.insert(0, str(REPO / "docs"))

from _tables import PSTAR_KEV, PSTAR_MEV_CM2_G, edep, eemean, xs  # noqa: E402

from ion_track_mc.analysis.beam import run_beam  # noqa: E402
from ion_track_mc.analysis.stopping_power import (  # noqa: E402
    analytic_equilibrium_stopping_power, analytic_stopping_power, equilibrium_charge_fractions)
from ion_track_mc.physics.channels import Channel  # noqa: E402
from ion_track_mc.physics.medium import WATER  # noqa: E402
from ion_track_mc.scoring.dose import RadialDoseScorer, make_bragg_curve_scorer  # noqa: E402
from ion_track_mc.transport.engine import TrackEngine  # noqa: E402
from ion_track_mc.transport.particle import Particle  # noqa: E402

QUICK = "--quick" in sys.argv
TRACKS_ONLY = "--tracks-only" in sys.argv
CH = {Channel.IONIZATION: "ion", Channel.EXCITATION: "exc", Channel.CAPTURE: "cap",
      Channel.LOSS: "los", Channel.ELASTIC: "ela"}


def r(a, nd=4):
    """Round to `nd` significant digits (keeps data.js small)."""
    return [float(f"{v:.{nd}g}") for v in np.asarray(a, float)]


def track(E_eV, z_max_nm=None, seed=1):
    """One proton track: per collision z, x [nm], process, charge before, energy after [keV], ΔE [eV]."""
    eng = TrackEngine(WATER, xs, edep_tables=edep, eemean_tables=eemean, seed=seed)
    p = Particle(kind="proton", position=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]), energy=E_eV, q=1)
    eng.transport(p)
    dep = [d for d in eng.deposits if d.kind == "proton"]
    dE = np.array([c[1] for c in eng.collision_log])
    E_after = (E_eV - np.cumsum(dE)) / 1e3
    n = len(dep) if z_max_nm is None else sum(1 for d in dep if d.position[2] <= z_max_nm)
    # A deposit records the charge state the projectile had when it arrived
    # (capture/loss change q after the deposit), so d.q is the state on the
    # segment that ends at that collision.
    q_in = [d.q for d in dep]
    return {"E0_keV": E_eV / 1e3, "n_total": len(dep),
            "range_nm": round(float(max(d.position[2] for d in dep)), 2),
            "z": [round(float(d.position[2]), 3) for d in dep[:n]],
            "x": [round(float(d.position[0]), 4) for d in dep[:n]],
            "c": [CH[d.channel] for d in dep[:n]],
            "q": q_in[:n],
            "E": [round(float(e), 3) for e in E_after[:n]],
            "dE": [round(float(e), 2) for e in dE[:n]]}


(APP / "docs").mkdir(exist_ok=True)
for pdf in ("tutorial.pdf", "tutorial_es.pdf"):
    shutil.copy2(REPO / "docs" / "tutorial" / pdf, APP / "docs" / pdf)
    # cover thumbnail for the Tutorial section (first slide; needs pdftoppm from poppler-utils)
    subprocess.run(["pdftoppm", "-png", "-f", "1", "-l", "1", "-scale-to-x", "720", "-scale-to-y", "-1",
                    "-singlefile", str(APP / "docs" / pdf), str(APP / "docs" / pdf.replace(".pdf", "_cover"))], check=True)

t0 = time.time()
out = {"generated": time.strftime("%Y-%m-%d"), "quick": QUICK}

# ---------------------------------------------------------------- tracks
out["tracks"] = {
    "t15_25": track(15e3, 25.0),
    "t15_full": track(15e3),
    "t500_130": track(500e3, 131.0),
}
# Whole 1 MeV track, decimated to (z, x, process) for the zoomable view.
full = track(1e6)
out["track1MeV"] = {"z": full["z"], "x": full["x"], "c": full["c"], "E": full["E"],
                    "n_total": full["n_total"], "range_nm": full["range_nm"]}
print(f"tracks done ({time.time() - t0:.0f} s)")

# ---------------------------------------------------------------- cross sections
xsd = {}
for q in (1, 0):
    e_min, e_max = xs.energy_range(q)
    E = np.geomspace(max(e_min, 8e3), min(e_max, 1e7), 160)
    d = {"E_keV": r(E / 1e3)}
    for ch in xs.channels(q):
        d[CH[ch]] = r([xs.sigma(ch, e, q) for e in E])
    d["mfp_nm"] = r([1.0 / (WATER.N * xs.sigma_total(e, q)) for e in E])
    xsd[str(q)] = d
out["xs"] = xsd

# ---------------------------------------------------------------- charge fractions
E = np.geomspace(8e3, 1e6, 160)
f = np.array([equilibrium_charge_fractions(xs, e) for e in E])
out["charge"] = {"E_keV": r(E / 1e3), "fp": r(f[:, 0]), "f0": r(f[:, 1])}

# ---------------------------------------------------------------- stopping power
Ek = np.geomspace(8, 1e5, 160)
out["sp"] = {
    "model": {"E_keV": r(Ek),
              "S": r([analytic_equilibrium_stopping_power(WATER, xs, e * 1e3, edep, eemean) for e in Ek]),
              "Sp": r([analytic_stopping_power(WATER, xs, e * 1e3, 1, edep, eemean) for e in Ek]),
              "S0": r([analytic_stopping_power(WATER, xs, e * 1e3, 0, edep, eemean) for e in Ek])},
    "pstar": {"E_keV": PSTAR_KEV.tolist(), "S": r(PSTAR_MEV_CM2_G * WATER.density_g_cm3 * 0.1)},
    "exp": [],
}
DC = REPO / "data_comparison" / "stopping_power_water"
for name, fname, phase, doi in (
        ("Phillips 1953", "Exp_SP_Phillips_1953.dat", "vapour", "10.1103/PhysRev.90.532"),
        ("Reynolds et al. 1953", "Exp_SP_Reynolds_1953.dat", "vapour", "10.1103/PhysRev.92.742"),
        ("Mitterschiffthaler & Bauer 1990", "Exp_SP_Mitterschiffthaler_1990.dat", "vapour", "10.1016/0168-583X(90)90073-4"),
        ("Baek et al. 2006", "Exp_SP_Baek_2006.dat", "vapour", "10.1093/rpd/ncl514"),
        ("Siiskonen et al. 2011", "Exp_SP_Siiskonen_2011.dat", "liquid", "10.1088/0031-9155/56/8/003")):
    a = np.loadtxt(DC / fname)
    out["sp"]["exp"].append({"name": name, "phase": phase, "doi": doi,
                             "E_keV": r(a[:, 0]), "S": r(a[:, 1]), "err": r(a[:, 2])})
table = []
for e in (10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 100000):
    m = analytic_equilibrium_stopping_power(WATER, xs, e * 1e3, edep, eemean)
    p = float(np.interp(np.log(e), np.log(PSTAR_KEV), PSTAR_MEV_CM2_G)) * WATER.density_g_cm3 * 0.1
    table.append([e, round(m, 1), round(p, 1), round(m / p, 2)])
out["sp"]["table"] = table

# ---------------------------------------------------------------- beam: Bragg curve, radial dose, LET
if TRACKS_ONLY:   # reuse the (slow) beam results already in data.js
    old = json.loads((APP / "data.js").read_text(encoding="utf-8").split("=", 1)[1].rstrip().rstrip(";"))
    out["beam"], out["quick"] = old["beam"], old["quick"]
    (APP / "data.js").write_text("// Generated by dev_tools/export_web_data.py from ion_track_mc. Do not edit by hand.\n"
                                 "window.ITMC_DATA = " + json.dumps(out, separators=(",", ":"), ensure_ascii=False) + ";\n",
                                 encoding="utf-8")
    print("tracks/curves updated, beam kept")
    sys.exit(0)
N_P = 5 if QUICK else 50
z_edges = np.linspace(0, 26e3, 105)
bragg = make_bragg_curve_scorer(WATER, z_edges=z_edges, transverse_half_width_nm=1e6)
r_edges = np.geomspace(0.05, 200, 31)
radial = RadialDoseScorer(WATER, r_edges=r_edges, length_nm=23e3)
res = run_beam(WATER, xs, edep, eemean, n_particles=N_P, E_beam_eV=1e6,
               scorers=[bragg, radial], seed=0, let_segment_length_nm=1000.0)
D = bragg.dose_Gy()[0, 0, :]
Dr = radial.dose_Gy() / N_P
m = Dr > 0
out["beam"] = {"n": N_P, "E_keV": 1000,
               "z_um": r(0.5 * (z_edges[1:] + z_edges[:-1]) / 1e3), "D": r(D / D.max()),
               "r_nm": r(np.sqrt(r_edges[1:] * r_edges[:-1])[m]), "Dr": r(Dr[m]),
               "LET_T": round(res.track_averaged_LET(), 2), "LET_D": round(res.dose_averaged_LET(), 2)}
print(f"beam done ({time.time() - t0:.0f} s): LET_T {out['beam']['LET_T']}, LET_D {out['beam']['LET_D']}")

js = ("// Generated by dev_tools/export_web_data.py from ion_track_mc. Do not edit by hand.\n"
      "window.ITMC_DATA = " + json.dumps(out, separators=(",", ":"), ensure_ascii=False) + ";\n")
(APP / "data.js").write_text(js, encoding="utf-8")
print(f"wrote data.js ({len(js) / 1e3:.0f} kB) in {time.time() - t0:.0f} s")
