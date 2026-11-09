"""Mesh quality report for the volume elements of a VTU (tets and/or hexes), via vtkMeshQuality.

Reports, per element type: count, element volume (min / median / max, negative count),
aspect ratio, radius ratio (tets), edge ratio, minimum and maximum dihedral angle (tets),
scaled Jacobian, and the fraction of elements beyond common "poor element" thresholds.

Credits: element metrics from VTK vtkMeshQuality (Verdict); max dihedral angle computed here. See REFERENCES.md.
"""
import json
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

_T = vtk.vtkMeshQuality.QualityMeasureTypes
TET_MEASURES = {
    "aspect_ratio": _T.ASPECT_RATIO,
    "radius_ratio": _T.RADIUS_RATIO,
    "edge_ratio": _T.EDGE_RATIO,
    "min_dihedral_angle_deg": _T.MIN_ANGLE,
    "scaled_jacobian": _T.SCALED_JACOBIAN,
    "shape": _T.SHAPE,
    "volume": _T.VOLUME,
}
HEX_MEASURES = {
    "aspect_ratio": _T.MAX_ASPECT_FROBENIUS,
    "edge_ratio": _T.EDGE_RATIO,
    "skew": _T.SKEW,
    "scaled_jacobian": _T.SCALED_JACOBIAN,
    "shape": _T.SHAPE,
    "volume": _T.VOLUME,
}


def _stats(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    return {
        "min": float(v.min()), "p01": float(np.percentile(v, 1)), "median": float(np.median(v)),
        "mean": float(v.mean()), "p99": float(np.percentile(v, 99)), "max": float(v.max()),
    }


def _extract_type(mesh, vtk_type):
    types = vtk_to_numpy(mesh.GetDistinctCellTypesArray()) if False else np.array([mesh.GetCellType(i) for i in range(mesh.GetNumberOfCells())], dtype=np.uint8) if mesh.GetNumberOfCells() < 200000 else vtk_to_numpy(mesh.GetCellTypesArray())
    idx = np.flatnonzero(types == vtk_type)
    if idx.size == 0:
        return None, 0
    sel_ids = vtk.vtkIdTypeArray(); sel_ids.SetNumberOfTuples(idx.size)
    for i, c in enumerate(idx):
        sel_ids.SetValue(i, int(c))
    node = vtk.vtkSelectionNode(); node.SetFieldType(vtk.vtkSelectionNode.CELL)
    node.SetContentType(vtk.vtkSelectionNode.INDICES); node.SetSelectionList(sel_ids)
    sel = vtk.vtkSelection(); sel.AddNode(node)
    ex = vtk.vtkExtractSelection(); ex.SetInputData(0, mesh); ex.SetInputData(1, sel); ex.Update()
    return ex.GetOutput(), int(idx.size)


def _max_dihedral_tets(ug):
    """Max dihedral angle per tet, computed directly (vtkMeshQuality only offers the minimum)."""
    pts = vtk_to_numpy(ug.GetPoints().GetData())
    conn = vtk_to_numpy(ug.GetCells().GetConnectivityArray()).reshape(-1, 4)
    p = pts[conn]  # (n,4,3)
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    normals = []
    for a, b, c in faces:
        n = np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a])
        n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-300
        # orient outward from the opposite vertex
        d = 6 - (a + b + c)
        s = np.sign(np.einsum("ij,ij->i", n, p[:, a] - p[:, d]))
        s[s == 0] = 1
        normals.append(n * s[:, None])
    normals = np.stack(normals, axis=1)  # (n,4,3)
    maxang = np.zeros(len(p))
    for i in range(4):
        for j in range(i + 1, 4):
            cos = np.einsum("ij,ij->i", normals[:, i], normals[:, j])
            ang = np.degrees(np.pi - np.arccos(np.clip(cos, -1, 1)))
            maxang = np.maximum(maxang, ang)
    return maxang


def quality_report(mesh, units="cm"):
    report = {"units": units, "n_points": int(mesh.GetNumberOfPoints()), "element_types": {}}
    for name, vtk_type, measures in (("tet4", vtk.VTK_TETRA, TET_MEASURES),
                                     ("hex8", vtk.VTK_HEXAHEDRON, HEX_MEASURES),
                                     ("wedge6", vtk.VTK_WEDGE, {"volume": _T.VOLUME})):
        sub, n = _extract_type(mesh, vtk_type)
        if n == 0:
            continue
        entry = {"count": n, "metrics": {}}
        for mname, measure in measures.items():
            mq = vtk.vtkMeshQuality(); mq.SetInputData(sub)
            mq.SaveCellQualityOn()
            if vtk_type == vtk.VTK_TETRA:
                mq.SetTetQualityMeasure(int(measure))
            elif vtk_type == vtk.VTK_HEXAHEDRON:
                mq.SetHexQualityMeasure(int(measure))
            else:
                mq.SetWedgeQualityMeasure(int(measure))
            mq.Update()
            vals = vtk_to_numpy(mq.GetOutput().GetCellData().GetArray("Quality"))
            if mname == "volume":
                entry["negative_volume_count"] = int((vals <= 0).sum())
                entry["total_volume"] = float(vals.sum())
            entry["metrics"][mname] = _stats(vals)
            if mname == "aspect_ratio":
                entry["fraction_aspect_ratio_gt_5"] = float((vals > 5).mean())
                entry["fraction_aspect_ratio_gt_10"] = float((vals > 10).mean())
            if mname == "scaled_jacobian":
                entry["fraction_scaled_jacobian_lt_0p2"] = float((vals < 0.2).mean())
            if mname == "min_dihedral_angle_deg":
                entry["fraction_min_dihedral_lt_10deg"] = float((vals < 10).mean())
        if vtk_type == vtk.VTK_TETRA:
            mx = _max_dihedral_tets(sub)
            entry["metrics"]["max_dihedral_angle_deg"] = _stats(mx)
            entry["fraction_max_dihedral_gt_160deg"] = float((mx > 160).mean())
        report["element_types"][name] = entry
    return report


def report_markdown(report, title=""):
    lines = []
    if title:
        lines.append(f"### {title}")
    lines.append(f"Points: {report['n_points']:,}")
    for etype, e in report["element_types"].items():
        lines.append("")
        lines.append(f"**{etype}**: {e['count']:,} elements, total volume {e.get('total_volume', float('nan')):.4g} {report['units']}^3, "
                     f"negative-volume elements: {e.get('negative_volume_count', 0)}")
        lines.append("")
        lines.append("| metric | min | p01 | median | mean | p99 | max |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for m, s in e["metrics"].items():
            if s is None:
                continue
            lines.append(f"| {m} | {s['min']:.4g} | {s['p01']:.4g} | {s['median']:.4g} | {s['mean']:.4g} | {s['p99']:.4g} | {s['max']:.4g} |")
        flags = {k: v for k, v in e.items() if k.startswith("fraction_")}
        if flags:
            lines.append("")
            lines.append("Poor-element fractions: " + ", ".join(f"{k.replace('fraction_', '')} = {v:.3%}" for k, v in flags.items()))
    return "\n".join(lines)


def read_ug(path):
    r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(path); r.Update()
    return r.GetOutput()


if __name__ == "__main__":
    path = sys.argv[1]
    rep = quality_report(read_ug(path))
    print(report_markdown(rep, path))
    if len(sys.argv) > 2:
        json.dump(rep, open(sys.argv[2], "w"), indent=2)
