"""Headless wrapper around the CFD Mesh Generator pipeline (SlicerExtension-VMTK) and
SimVascular Mesh Prep (svmeshcomplete / svmpsetup), so a labelled surface becomes a
mesh-complete case from a terminal with the same code the Slicer panels run.

Credits: wraps CfdMeshGeneratorLib.MeshingPipeline from SlicerExtension-VMTK (VMTK: Antiga et al. 2008,
doi:10.1007/s11517-008-0420-1) and Netgen (Schoberl 1997, doi:10.1007/s007910050004). See REFERENCES.md.
"""
import os, sys, time, json, logging
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

VMTK_LIB = os.environ.get("CFDMESHGEN_LIB",
    os.path.join(os.path.dirname(__file__), "..", "ext", "SlicerExtension-VMTK", "CfdMeshGenerator"))
for cand in (VMTK_LIB, "/home/claude/ext/SlicerExtension-VMTK/CfdMeshGenerator"):
    if os.path.isdir(cand) and cand not in sys.path:
        sys.path.insert(0, cand)
# The Slicer extension imports VMTK's wrapped libraries by their bare module names (as Slicer
# exposes them); the pip "vmtk" wheel puts them inside the vmtk package. Alias them.
import importlib  # noqa: E402
for _name in ("vtkvmtkCommonPython", "vtkvmtkMiscPython", "vtkvmtkComputationalGeometryPython",
              "vtkvmtkDifferentialGeometryPython", "vtkvmtkContribPython", "vtkvmtkIOPython",
              "vtkvmtkSegmentationPython", "vtkvmtkITKPython"):
    if _name not in sys.modules:
        try:
            sys.modules[_name] = importlib.import_module("vmtk." + _name)
        except Exception as _e:  # pragma: no cover
            log_import_err = _e
from CfdMeshGeneratorLib.MeshingPipeline import MeshingPipeline, Mesher  # noqa: E402

log = logging.getLogger("cvbench")


def read_polydata(path):
    ext = os.path.splitext(path)[1].lower()
    r = {".vtp": vtk.vtkXMLPolyDataReader, ".stl": vtk.vtkSTLReader, ".ply": vtk.vtkPLYReader,
         ".vtk": vtk.vtkPolyDataReader}[ext]()
    r.SetFileName(path); r.Update()
    return r.GetOutput()


def write_dataset(ds, path):
    if isinstance(ds, vtk.vtkPolyData):
        w = vtk.vtkXMLPolyDataWriter()
    else:
        w = vtk.vtkXMLUnstructuredGridWriter()
    w.SetFileName(path); w.SetInputData(ds); w.SetCompressorTypeToZLib(); w.SetDataModeToAppended(); w.Write()
    return path


def face_summary(surface, array="ModelFaceID"):
    ids = vtk_to_numpy(surface.GetCellData().GetArray(array))
    out = {}
    for fid in np.unique(ids):
        sel = np.flatnonzero(ids == fid)
        thr = vtk.vtkThreshold(); thr.SetInputData(surface)
        thr.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, array)
        thr.SetLowerThreshold(fid); thr.SetUpperThreshold(fid); thr.Update()
        geo = vtk.vtkGeometryFilter(); geo.SetInputData(thr.GetOutput()); geo.Update()
        pd = geo.GetOutput()
        mp = vtk.vtkMassProperties(); mp.SetInputData(pd); mp.Update()
        area = mp.GetSurfaceArea()
        com = vtk.vtkCenterOfMass(); com.SetInputData(pd); com.Update()
        out[int(fid)] = dict(ncells=int(len(sel)), area=float(area),
                             eq_diameter=float(2*np.sqrt(area/np.pi)), centroid=list(com.GetCenter()))
    return out


def generate(surface, edge_length, *, cell_entity_ids="ModelFaceID", cap_ids="CapID",
             mesher="netgen", volume_scale=0.8, boundary_layer=False, n_sublayers=2,
             bl_thickness_factor=0.25, skip_capping=True, log_cb=None):
    """Run the CFD Mesh Generator pipeline on a closed, face-labelled surface."""
    p = MeshingPipeline()
    if log_cb: p.stepCallback = log_cb
    t0 = time.time()
    mesh = p.generateMesh(surface,
                          targetEdgeLength=edge_length,
                          cellEntityIdsArrayName=cell_entity_ids,
                          capIdsArrayName=cap_ids,
                          skipCapping=skip_capping,
                          mesher=mesher,
                          volumeElementScaleFactor=volume_scale,
                          boundaryLayer=boundary_layer,
                          numberOfSubLayers=n_sublayers,
                          boundaryLayerThicknessFactor=bl_thickness_factor,
                          tetrahedralize=boundary_layer)
    if isinstance(mesh, tuple):
        mesh, remeshed_surface = mesh
    if p.lastTetrahedralizationFailed:
        raise RuntimeError("tetrahedralization failed")
    return mesh, time.time() - t0


def count_cells(mesh):
    types = vtk_to_numpy(mesh.GetCellTypesArray())
    return {"tet4": int((types == vtk.VTK_TETRA).sum()), "hex8": int((types == vtk.VTK_HEXAHEDRON).sum()),
            "wedge6": int((types == vtk.VTK_WEDGE).sum()), "tri3": int((types == vtk.VTK_TRIANGLE).sum()),
            "quad4": int((types == vtk.VTK_QUAD).sum()), "points": int(mesh.GetNumberOfPoints())}
