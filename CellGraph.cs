using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace VolumeGeneratorBasedOnGraph
{
    public class CellGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CellGraph class.
        /// </summary>
        public CellGraph()
          : base("CellGraph", "CellGraph",
              "Cells to be represnted as graph nodes, whose adjacencies to be represented as graph links",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("DualConvexPolygon", "DualConvexPolygon", "Cells to be represnted as graph nodes, whose adjacencies to be represented as graph links", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("DualConvexConnectivityGraph", "Graph", "A graph representation as an adjacency list datatree", GH_ParamAccess.list);
            pManager.AddPointParameter("PointsRespresentDualConvex", "GraphVertices", "Graph vertices as [edge] centrpoids of cells", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> dualConvexPolygonCurve = new List<Curve>();

            // 判断输入的DualConvexPolygon是不是闭合的
            if (DA.GetDataList<Curve>("DualConvexPolygon", dualConvexPolygonCurve))
            {
                int count = 0;
                foreach (Curve curve in dualConvexPolygonCurve)
                {
                    if (!curve.IsClosed)
                    {
                        count++;
                    }
                }
                if (count > 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All curves must be closed!");
                }
            }


            // 用对偶ConvexPolygon的中心点来代表这个ConvexPolygon
            DataTree<int> dualConvexPolygonConnectivityDT = new DataTree<int>();
            List<Point3d> dualConvexCenterPoints = new List<Point3d>();

            for (int i = 0; i < dualConvexPolygonCurve.Count; i++)
            {
                dualConvexPolygonConnectivityDT.EnsurePath(i);

                Polyline polyline = null;
                if (dualConvexPolygonCurve[i].TryGetPolyline(out polyline))
                {
                    dualConvexCenterPoints.Add(polyline.CenterPoint());
                }

                for (int j = 0; j < dualConvexPolygonCurve.Count; j++)
                {
                    if (j != i)
                    {
                        // 判断两条Polyline是否相交，如果是，那么它俩有邻接关系
                        RegionContainment regionContainment = Curve.PlanarClosedCurveRelationship(dualConvexPolygonCurve[i], dualConvexPolygonCurve[j], Plane.WorldXY, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                        if (regionContainment == RegionContainment.MutualIntersection)
                        {
                            dualConvexPolygonConnectivityDT.Branch(i).Add(j);
                        }
                    }
                }
            }

            DA.SetDataTree(0, dualConvexPolygonConnectivityDT);
            DA.SetDataList("PointsRespresentDualConvex", dualConvexCenterPoints);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("a9a2f070-d7e0-4391-816a-752daa2e9e1b"); }
        }
    }
}