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
              "Description",
              "Constructs a graph from a set of adjacent polygoncells", "Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Cell", "Cell", "Cells to be represnted as graph nodes, whose adjacencies to be represented as graph links", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Graph", "Graph", "A graph representation as an adjacency list datatree", GH_ParamAccess.list);
            pManager.AddPointParameter("GraphVertices", "GraphVertices", "Graph vertices as [edge] centrpoids of cells", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> list = new List<Curve>();

            if (DA.GetDataList<Curve>("Cell", list))
            {
                int count = 0;
                foreach (Curve curve in list)
                {
                    if (!curve.IsClosed)
                    {
                        count++;
                    }
                }

                bool flag = count > 0;
                if (flag)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All curves must be closed!");
                }
            }

            DataTree<int> dataTree = new DataTree<int>();
            List<Point3d> list2 = new List<Point3d>();

            for (int i = 0; i < list.Count; i++)
            {
                dataTree.EnsurePath(new int[]
                {
                    i
                });

                Polyline polyline = null;
                if (list[i].TryGetPolyline(out polyline))
                {
                    list2.Add(polyline.CenterPoint());
                }

                for (int j = 0; j < list.Count; j++)
                {
                    if (j != i)
                    {
                        RegionContainment regionContainment = Curve.PlanarClosedCurveRelationship(list[i], list[j], Plane.WorldXY, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                        if (regionContainment == RegionContainment.MutualIntersection)
                        {
                            dataTree.Branch(i).Add(j);
                        }
                    }
                }
            }

            DA.SetDataTree(0, dataTree);
            DA.SetDataList("GraphVertices", list2);
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