"""Boundary-condition construction for the benchmark cases (CGS units: cm, g, s, dyn).

Inflow waveforms are analytic, documented shapes scaled to a target mean flow, so every case is
fully reproducible without patient data. Outlet RCR / resistance values are derived from a target
mean pressure and a flow split (prescribed or Murray-law-by-area), the way a first-pass
SimVascular tuning is done. Replace with VMR `rcrt.dat` / `inflow.flow` when the VMR project
is available: `from_vmr_project()` does that.

Credits: Inflow/RCR/Resistance types from SlicerSimVascular svromsetup; Windkessel after Westerhof et al.
2009 (doi:10.1007/s11517-008-0359-2); Murray 1926 (doi:10.1073/pnas.12.3.207). See REFERENCES.md.
"""
import math
import numpy as np

from svromsetup.boundary_conditions import Inflow, RCR, Resistance, read_rcrt, read_flow_file

MMHG = 1333.22  # dyn/cm^2 per mmHg


def aortic_waveform(mean_flow, period=1.0, n=101, systole_fraction=0.35, backflow_fraction=0.08):
    """Smooth adult-aortic-type flow waveform (cm^3/s): half-sine systolic ejection, a short
    negative dip at valve closure, near-zero diastolic flow; scaled to the requested mean."""
    t = np.linspace(0.0, period, n)
    ts = systole_fraction * period
    tb = backflow_fraction * period
    q = np.zeros_like(t)
    sys_mask = t <= ts
    q[sys_mask] = np.sin(np.pi * t[sys_mask] / ts)
    back = (t > ts) & (t <= ts + tb)
    q[back] = -0.12 * np.sin(np.pi * (t[back] - ts) / tb)
    q += 0.02  # small diastolic forward flow
    q[-1] = q[0]
    mean = np.trapezoid(q, t) / period
    q *= mean_flow / mean
    return Inflow(tuple(float(x) for x in t), tuple(float(x) for x in q), source="analytic_aortic")


def venous_waveform(mean_flow, period=1.0, n=61, pulsatility=0.15, phase=0.0):
    """Low-pulsatility caval inflow for Fontan / Glenn cases: mean plus a 15% respiratory-cardiac
    modulation. Steady if pulsatility == 0."""
    t = np.linspace(0.0, period, n)
    q = mean_flow * (1.0 + pulsatility * np.sin(2 * np.pi * t / period + phase))
    q[-1] = q[0]
    return Inflow(tuple(float(x) for x in t), tuple(float(x) for x in q), source="analytic_caval")


def murray_split(areas, exponent=1.0):
    """Flow fraction per outlet proportional to area**exponent (exponent 1 == area, 1.5 == Murray's
    law in diameter^3)."""
    a = np.array([areas[k] for k in areas], float) ** exponent
    frac = a / a.sum()
    return dict(zip(areas.keys(), frac.tolist()))


def rcr_from_targets(outlet_areas, mean_flow, mean_pressure_mmHg, *, flow_split=None,
                     rp_fraction=0.06, total_compliance=1.0e-3, distal_pressure_mmHg=0.0,
                     area_exponent=1.0):
    """Three-element Windkessel per outlet.

    R_i = (P_mean - P_d) / Q_i with Q_i = f_i * Q_mean; Rp = rp_fraction * R, Rd = R - Rp;
    compliance distributed C_i = C_total * (f_i) (proportional to flow, i.e. inversely to R).
    `flow_split` may fix some outlets' fractions; the remainder is shared by area.
    """
    names = list(outlet_areas)
    frac = {}
    fixed = dict(flow_split or {})
    fixed_total = sum(fixed.values())
    if fixed_total > 1.0 + 1e-9:
        raise ValueError("prescribed flow split exceeds 1")
    free = {k: outlet_areas[k] for k in names if k not in fixed}
    if free:
        for k, f in murray_split(free, area_exponent).items():
            frac[k] = f * (1.0 - fixed_total)
    frac.update(fixed)
    P = (mean_pressure_mmHg - distal_pressure_mmHg) * MMHG
    out = {}
    for k in names:
        q = frac[k] * mean_flow
        R = P / q
        out[k] = RCR(Rp=float(rp_fraction * R), C=float(total_compliance * frac[k]),
                     Rd=float((1 - rp_fraction) * R), Pd=float(distal_pressure_mmHg * MMHG))
    return out, frac


def resistance_from_targets(outlet_areas, total_flow, pvr_wood_units=None, mean_pressure_mmHg=None,
                            distal_pressure_mmHg=0.0, flow_split=None, area_exponent=1.0):
    """Pure resistance outlets (Fontan/Glenn pulmonary beds). Total resistance from PVR in Wood
    units (mmHg·min/L) or from a mean-pressure target; split by prescribed fractions / area."""
    names = list(outlet_areas)
    fixed = dict(flow_split or {})
    free = {k: outlet_areas[k] for k in names if k not in fixed}
    frac = {}
    if free:
        for k, f in murray_split(free, area_exponent).items():
            frac[k] = f * (1.0 - sum(fixed.values()))
    frac.update(fixed)
    if pvr_wood_units is not None:
        # Wood unit = mmHg per L/min -> dyn s / cm^5 : 1 WU = 1333.22 / (1000/60) = 79.99 dyn s/cm^5
        R_total = pvr_wood_units * 79.993
    else:
        R_total = (mean_pressure_mmHg - distal_pressure_mmHg) * MMHG / total_flow
    out = {k: Resistance(R=float(R_total / frac[k])) for k in names}
    return out, frac


def describe(conditions, frac=None):
    rows = []
    for name, c in conditions.items():
        if isinstance(c, RCR):
            rows.append(f"| {name} | RCR | Rp={c.Rp:.5g} C={c.C:.3g} Rd={c.Rd:.5g} Pd={c.Pd:.0f} (dyn·s/cm⁵, cm⁵/dyn, dyn/cm²) | {100*frac[name]:.1f}% |" if frac else f"| {name} | RCR | Rp={c.Rp:.5g} C={c.C:.3g} Rd={c.Rd:.5g} | |")
        elif isinstance(c, Resistance):
            rows.append(f"| {name} | R | R={c.R:.5g} dyn·s/cm⁵ | {100*frac[name]:.1f}% |" if frac else f"| {name} | R | R={c.R:.5g} | |")
        elif isinstance(c, Inflow):
            rows.append(f"| {name} | inflow | mean {c.mean:.3g} cm³/s, period {c.period} s, {len(c.time)} pts ({c.source}) | |")
    return "\n".join(["| face | type | values | flow share |", "|---|---|---|---|"] + rows)


def from_vmr_project(project_dir):
    """Pull svSolver boundary conditions out of a downloaded VMR SimVascular project:
    `flow-files/inflow_3d.flow` (negative = into the model in svSolver) and
    `Simulations/*/rcrt.dat`. Returns (Inflow, {face name: RCR}) or raises if absent."""
    import glob, os
    flows = glob.glob(os.path.join(project_dir, "flow-files", "*.flow")) + \
        glob.glob(os.path.join(project_dir, "Simulations", "*", "*.flow"))
    rcrts = glob.glob(os.path.join(project_dir, "Simulations", "*", "rcrt.dat")) + \
        glob.glob(os.path.join(project_dir, "Simulations", "*", "*", "rcrt.dat"))
    inflow = read_flow_file(flows[0]) if flows else None
    if inflow is not None and inflow.mean < 0:
        inflow = inflow.flipped()
    rcr = read_rcrt(rcrts[0]) if rcrts else {}
    return inflow, rcr, {"flow_files": flows, "rcrt_files": rcrts}
