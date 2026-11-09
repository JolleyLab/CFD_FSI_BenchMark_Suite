#!/usr/bin/env python3
"""Physics and boundary-condition sanity checks for a built case (called by verify_case.py).

Checks, each a PASS / WARN / FAIL line:
  * BC completeness: every cap in solver.xml has exactly one condition; exactly the inflow caps
    are Dirichlet; outlets are RCR or resistance with positive values.
  * Inflow sign and period: the .flow file drives flow INTO the domain (negative in the
    svMultiPhysics convention), its header point count matches its rows, period > 0.
  * 0D pressure estimate: with every outlet's total resistance (Rp + Rd, or R) in parallel to
    its distal pressure, mean pressure = Q_mean * R_total + Pd. Systemic cases are expected in
    50-130 mmHg, Fontan/pulmonary cases in 5-25 mmHg (WARN outside).
  * Flow regime: peak inlet Reynolds number Re = 4 rho Q / (pi mu D) and Womersley number
    alpha = R sqrt(omega rho / mu), reported; WARN if peak Re > 4000 (laminar assumption).
  * Time step: convective CFL = u_peak dt / h_min_inlet, reported with the dt that gives CFL 1.
  * Caps: each cap is planar (out-of-plane deviation < 2 % of its diameter) and its
    triangles face out of the domain.
  * MFEM smoke solve (PyMFEM, levels up to 2 M elements): a Laplace problem with u = 1 on the
    inflow cap(s), u = 0 on the outlets and natural (no-flux) walls converges with PCG, and the
    flux balance (inflow = outflow) closes to 1 %. This proves the mesh, its attributes and the
    BC binding work in an MFEM solve, not only that the file loads.

References: Celik et al. 2008, doi:10.1115/1.2960953 (discretization uncertainty);
Westerhof et al. 2009, doi:10.1007/s11517-008-0359-2 (Windkessel). See REFERENCES.md.
"""
import json, math, os, re
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

MMHG = 1333.22


def _read(path):
    r = (vtk.vtkXMLPolyDataReader if path.endswith(".vtp") else vtk.vtkXMLUnstructuredGridReader)(); r.SetFileName(path); r.Update(); return r.GetOutput()


def _bcs(xml):
    out = {}
    for m in re.finditer(r'<Add_BC name="([^"]+)"\s*>(.*?)</Add_BC>', xml, flags=re.S):
        name, body = m.group(1), m.group(2)
        tag = lambda t: (re.search(rf"<{t}>\s*([^<]+?)\s*</{t}>", body) or [None, None])[1]
        out[name] = {"type": tag("Type"), "time": tag("Time_dependence"), "file": tag("Temporal_values_file_path"),
                     "Rp": tag("Proximal_resistance"), "Rd": tag("Distal_resistance"), "C": tag("Capacitance"),
                     "Pd": tag("Distal_pressure"), "R": tag("Value") if tag("Time_dependence") == "Resistance" else None}
    return out


def run(case_dir, check, warn, note, mfem_max_elements=2.0e6):
    xml_path = os.path.join(case_dir, "solver.xml")
    if not os.path.exists(xml_path):
        return
    xml = open(xml_path).read()
    if "<Add_equation type=\"fluid\"" not in xml.replace("'", '"').replace("type=\"fluid\" ", "type=\"fluid\""):
        if 'type="FSI"' in xml:
            note("physics checks skipped (FSI case)"); return
    bcs = _bcs(xml)
    faces = re.findall(r'<Add_face name="([^"]+)"', xml)
    caps = [f for f in faces if f.startswith("cap_") or f in ("inlet", "outlet")]
    for c in caps:
        check(c in bcs, f"BC: cap {c} has a boundary condition")
    inflows = [c for c in caps if bcs.get(c, {}).get("type") == "Dir"]
    outlets = [c for c in caps if c not in inflows]
    check(len(inflows) >= 1, f"BC: {len(inflows)} Dirichlet inflow cap(s): {inflows}")
    rd = os.path.join(case_dir, "mesh", "mesh-surfaces")
    rho = float((re.search(r"<Density>\s*([0-9.eE+-]+)", xml) or [0, 1.06])[1])
    mu = float((re.search(r"<Viscosity[^>]*>\s*<Value>\s*([0-9.eE+-]+)", xml) or [0, 0.04])[1])
    # inflow
    Q_mean = 0.0; Q_peak = 0.0; period = None; d_in = None; h_in = None
    for c in inflows:
        b = bcs[c]
        if b["file"]:
            lines = [l.split() for l in open(os.path.join(case_dir, b["file"])).read().strip().splitlines()]
            npts = int(lines[0][0]); data = np.array([[float(x) for x in l[:2]] for l in lines[1:]])
            check(npts == len(data), f"BC: {b['file']} header lists {npts} points, file has {len(data)}")
            q = -data[:, 1]                       # svMultiPhysics: negative = into the domain
            period = float(data[-1, 0] - data[0, 0])
            check(period > 0, f"BC: {b['file']} period {period:.4g} s")
            qm = np.trapezoid(q, data[:, 0]) / period
            check(qm > 0, f"BC: {c} mean inflow {qm:.4g} cm^3/s into the domain (sign convention)")
            Q_mean += qm; Q_peak += q.max()
        else:
            v = float((re.search(rf'<Add_BC name="{c}"\s*>.*?<Value>\s*([0-9.eE+-]+)', xml, flags=re.S) or [0, 0])[1])
            Q_mean += -v; Q_peak += -v
        face = _read(os.path.join(rd, f"{c}.vtp"))
        tri = vtk.vtkTriangleFilter(); tri.SetInputData(face); tri.Update()
        mp = vtk.vtkMassProperties(); mp.SetInputData(tri.GetOutput()); mp.Update()
        d = 2 * math.sqrt(mp.GetSurfaceArea() / math.pi); d_in = d if d_in is None else max(d_in, d)
        ed = vtk.vtkExtractEdges(); ed.SetInputData(face); ed.Update(); e = ed.GetOutput()
        P = vtk_to_numpy(e.GetPoints().GetData()); L = vtk_to_numpy(e.GetLines().GetConnectivityArray()).reshape(-1, 2)
        hl = np.linalg.norm(P[L[:, 0]] - P[L[:, 1]], axis=1); h_in = hl.min() if h_in is None else min(h_in, hl.min())
    # outlets, 0D estimate
    conductance = 0.0; pd_w = 0.0
    for c in outlets:
        b = bcs.get(c, {})
        if b.get("time") == "RCR":
            Rp, Rd, C = float(b["Rp"]), float(b["Rd"]), float(b["C"])
            check(Rp >= 0 and Rd > 0 and C > 0, f"BC: {c} RCR Rp {Rp:.4g}, C {C:.3g}, Rd {Rd:.4g} positive")
            R = Rp + Rd; Pd = float(b["Pd"] or 0)
        elif b.get("time") == "Resistance":
            R = float(b["R"]); Pd = 0.0
            check(R > 0, f"BC: {c} resistance {R:.4g} positive")
        else:
            check(False, f"BC: {c} outlet type {b.get('time')} is RCR or Resistance"); continue
        conductance += 1.0 / R; pd_w += Pd / R
    if conductance > 0 and Q_mean > 0:
        P_mean = (Q_mean + pd_w) / conductance / MMHG
        pulmonary = any(k in case_dir.lower() for k in ("pulm", "tcpc", "fontan"))
        lo, hi = (5, 25) if pulmonary else (50, 130)
        warn(lo <= P_mean <= hi, f"0D: mean pressure {P_mean:.1f} mmHg from Q {Q_mean:.3g} cm^3/s and outlet resistances ({'pulmonary' if pulmonary else 'systemic'} range {lo}-{hi})")
    if d_in:
        Re = 4 * rho * Q_peak / (math.pi * mu * d_in)
        omega = 2 * math.pi / period if period else None
        alpha = (d_in / 2) * math.sqrt(omega * rho / mu) if omega else float("nan")
        u_peak = 2 * Q_peak / (math.pi * (d_in / 2) ** 2)     # parabolic centreline
        dt = float((re.search(r"<Time_step_size>\s*([0-9.eE+-]+)", xml) or [0, 1e-3])[1])
        cfl = u_peak * dt / h_in
        warn(Re < 4000, f"flow: peak inlet Re {Re:.0f}, Womersley {alpha:.1f}, inlet d {d_in:.2f} cm (laminar assumption warns above 4000)")
        note(f"time step: dt {dt:g} s gives inlet CFL {cfl:.2f} (u_peak {u_peak:.1f} cm/s, h_min {h_in:.3f} cm); CFL 1 at dt {h_in / u_peak:.2e} s")
    # caps planar and outward
    vol = _read(os.path.join(case_dir, "mesh", "mesh-complete.mesh.vtu"))
    centre = np.array(vol.GetCenter())
    for c in caps:
        f = os.path.join(rd, f"{c}.vtp")
        if not os.path.exists(f):
            continue
        face = _read(f); P = vtk_to_numpy(face.GetPoints().GetData())
        cen = P.mean(axis=0); u, s, vt = np.linalg.svd(P - cen, full_matrices=False); n = vt[2]
        dev = np.abs((P - cen) @ n).max(); diam = 2 * math.sqrt(max(s[0], 1e-30) ** 2 / len(P) * 2)
        span = np.ptp((P - cen) @ vt[0])
        warn(dev < 0.02 * span, f"cap {c}: planar (max out-of-plane {dev:.2e} cm, {100*dev/span:.2f} % of span)")
        nrm = vtk.vtkPolyDataNormals(); nrm.SetInputData(face); nrm.ComputeCellNormalsOn(); nrm.ComputePointNormalsOff(); nrm.AutoOrientNormalsOff(); nrm.ConsistencyOff(); nrm.Update()
        cn = vtk_to_numpy(nrm.GetOutput().GetCellData().GetNormals()).mean(axis=0)
        # outward = pointing away from the nearest volume interior: probe a small step inside along -cn
        loc = vtk.vtkCellLocator(); loc.SetDataSet(vol); loc.BuildLocator()
        step = 0.25 * (span / max(4, math.sqrt(face.GetNumberOfCells())))
        inside_back = loc.FindCell(list(cen - step * cn / np.linalg.norm(cn))) >= 0
        warn(inside_back, f"cap {c}: triangles face out of the domain")
    # MFEM smoke solve
    _mfem_solve(case_dir, inflows, outlets, check, note, mfem_max_elements)


def _mfem_solve(case_dir, inflows, outlets, check, note, max_el):
    try:
        import mfem.ser as mfem
    except ImportError:
        note("MFEM smoke solve skipped (PyMFEM not installed)"); return
    d = os.path.join(case_dir, "mfem")
    files = sorted(os.listdir(d)) if os.path.isdir(d) else []
    meshes = [m for m in files if m.endswith(".msh")] or [m for m in files if m.endswith(".mesh")]   # release ships .msh
    if not meshes:
        return
    path = os.path.join(d, meshes[0])
    attrs = json.load(open(path.rsplit(".", 1)[0] + "_attributes.json"))["boundary_attributes"]
    by_name = {v: int(k) for k, v in attrs.items()}
    in_ids = [by_name[n] for n in inflows if n in by_name]; out_ids = [by_name[n] for n in outlets if n in by_name]
    if not in_ids or not out_ids:
        note("MFEM smoke solve skipped (inflow/outlet attributes not found)"); return
    m = mfem.Mesh(path, 1, 1)
    if m.GetNE() > max_el:
        note(f"MFEM smoke solve skipped ({m.GetNE():,} elements > {max_el:,.0f}; run on a smaller level)"); return
    fes = mfem.FiniteElementSpace(m, mfem.H1_FECollection(1, 3))
    nb = m.bdr_attributes.Max()
    ess = mfem.intArray(nb); ess.Assign(0)
    for a in in_ids + out_ids:
        ess[a - 1] = 1
    tdofs = mfem.intArray(); fes.GetEssentialTrueDofs(ess, tdofs)
    x = mfem.GridFunction(fes); x.Assign(0.0)
    mark_in = mfem.intArray(nb); mark_in.Assign(0)
    for a in in_ids:
        mark_in[a - 1] = 1
    x.ProjectBdrCoefficient(mfem.ConstantCoefficient(1.0), mark_in)
    a = mfem.BilinearForm(fes); a.AddDomainIntegrator(mfem.DiffusionIntegrator()); a.Assemble()
    b = mfem.LinearForm(fes); b.Assemble()
    A = mfem.OperatorPtr(); B = mfem.Vector(); X = mfem.Vector()
    a.FormLinearSystem(tdofs, x, b, A, X, B)
    M = mfem.GSSmoother(mfem.Opr2SparseMat(A.Ptr()))
    solver = mfem.CGSolver(); solver.SetPreconditioner(M); solver.SetOperator(A.Ptr())
    solver.SetRelTol(1e-8); solver.SetMaxIter(4000); solver.SetPrintLevel(-1)
    solver.Mult(B, X)
    converged = solver.GetConverged(); iters = solver.GetNumIterations()
    a.RecoverFEMSolution(X, b, x)
    check(bool(converged), f"MFEM smoke solve: Laplace u=1 on inflow, u=0 on outlets converged in {iters} PCG iterations ({m.GetNE():,} elements)")
    # flux balance through boundary: residual r = A x (full, unconstrained) summed on inflow vs outlet dofs
    a2 = mfem.BilinearForm(fes); a2.AddDomainIntegrator(mfem.DiffusionIntegrator()); a2.Assemble(); a2.Finalize()
    r = mfem.Vector(fes.GetVSize()); a2.Mult(x, r); r = r.GetDataArray()
    def dofs_on(ids):
        mk = mfem.intArray(nb); mk.Assign(0)
        for i in ids:
            mk[i - 1] = 1
        lst = mfem.intArray(); fes.GetEssentialTrueDofs(mk, lst); return set(lst.ToList())
    din, dout = dofs_on(in_ids), dofs_on(out_ids)
    shared = din & dout
    fin = r[list(din - shared)].sum(); fout = r[list(dout - shared)].sum()
    bal = abs(fin + fout) / max(abs(fin), 1e-300)
    check(bal < 0.01, f"MFEM smoke solve: flux balance inflow {fin:.4g} vs outflow {-fout:.4g} ({100*bal:.3f} %)")
