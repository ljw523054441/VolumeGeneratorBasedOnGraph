using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class ConstructGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ConstructGraph class.
        /// </summary>
        public ConstructGraph()
          : base("ConstructGraph", "ConstructGraph",
              "构建图结构",
              "VolumeGeneratorBasedOnGraph", "Construct Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddPointParameter("VolumeNodePoints", "VolumeNodePoints", "能够代表volumeNode节点的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePoints", "BoundaryNodePoints", "能够代表除了NEWS点外的其他BoundaryNode节点的抽象点", GH_ParamAccess.list);
            pManager.AddIntegerParameter("VolumeConnectivityTree", "ConnectivityTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "AdjacencyTree", "体量邻近关系树", GH_ParamAccess.tree);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
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
            get { return new Guid("95f16532-4629-4afd-8227-75f14332553e"); }
        }
    }
}