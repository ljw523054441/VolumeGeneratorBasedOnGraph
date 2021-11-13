using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDebugGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDebugGraph class.
        /// </summary>
        public GhcDebugGraph()
          : base("DebugGraph", "DebugGraph",
              "显示图结构",
              "VolumeGeneratorBasedOnGraph", "Debug")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Graph", "G", "要显示的图结构", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("DebugGraphNode", "DebugGN", "", GH_ParamAccess.list);

            pManager.AddGenericParameter("DebugGraphLoL", "DebugGL", "", GH_ParamAccess.list);

            pManager.AddGenericParameter("DebugGraphVertex", "DebugGV", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Graph graph = new Graph();
            if (DA.GetData<Graph>("Graph", ref graph))
            {
                List<string> printGraphNode = UtilityFunctions.PrintGraphNode(graph);
                List<string> printGraphLoL = UtilityFunctions.PrintGraphLoL(graph);
                List<string> printGraphVertex = UtilityFunctions.PrintGraphVertex(graph);

                DA.SetDataList("DebugGraphNode", printGraphNode);
                DA.SetDataList("DebugGraphLoL", printGraphLoL);
                DA.SetDataList("DebugGraphVertex", printGraphVertex);
            }

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
            get { return new Guid("bfd48d0e-5fbf-42d9-b487-220032061370"); }
        }
    }
}