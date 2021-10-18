using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class WeightedVoroni : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the WeightedVoroni class.
        /// </summary>
        public WeightedVoroni()
          : base("WeightedVoroni", "WeightedVoroni",
              "Draws a Voroni diagram for points with various radii on a base plane; this diagram will be sort of a bubble diagram",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPlaneParameter("BasePlane", "BP", "The plane on which the diagram is to be drawn", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddPointParameter("Points", "P", "The centre points of bubbles/Voronoi cells", GH_ParamAccess.list);
            pManager.AddNumberParameter("Radii", "R", "The list of radii of cells correspoding to the list of points", GH_ParamAccess.list);
            pManager.AddCurveParameter("BoundingCurve", "BC", "A polyline confining the area of the diagram", GH_ParamAccess.item);
            pManager[3].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("AdjacencyGraph", "AG", "A graph (adjacency lists stored in a datatree) describing how the cells are adjacent to each other", GH_ParamAccess.tree);
            pManager.AddCurveParameter("VoronoiCells", "VC", "Bubbles in the form of Voronoi cells", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Point3d> list = new List<Point3d>();
            List<double> list2 = new List<double>();
            Curve bc = null;
            bool flag = DA.GetDataList<Point3d>("Points", list) & DA.GetDataList<double>("Radii", list2);
            if (flag)
            {
                Plane bp = Plane.WorldXY;
                DA.GetData<Plane>("BasePlane", ref bp);
                DA.GetData<Curve>("BoundingCurve", ref bc);

                DataTree<int> dataTree = null;
                List<Curve> list3 = 

            }
        }

        public List<Curve> AlphaVoronoi(
            Plane bp, 
            List<Point3d> v, 
            List<double> r, 
            Curve bc, 
            ref DataTree<int> adjacencyGraph)
        {
            bool flag = !bp.IsValid || v == null || r == null;
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
            get { return new Guid("18dd9201-de43-4ae4-88c4-be8a5db150ec"); }
        }
    }
}