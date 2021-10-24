using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class GlobalParameters : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Globals class.
        /// </summary>
        public GlobalParameters()
          : base("Globals", "SetGlobalVaribles",
              "设置全局变量",
              "VolumeGeneratorBasedOnGraph", "GlobalParameter")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("VolumeNodePoints", "VNPoints", "能够代表volumeNode节点的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePoints", "BNPoints", "能够代表BoundaryNode节点的抽象点（point）", GH_ParamAccess.list);
            // pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GP", "构造后的全局参数", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Point3d> volumeNodes = new List<Point3d>();
            List<Point3d> boundaryNodes = new List<Point3d>();

            GlobalParameter globalParameter;

            if (DA.GetDataList<Point3d>("VolumeNodePoints", volumeNodes)
                & DA.GetDataList<Point3d>("BoundaryNodePoints", boundaryNodes))
            {
                globalParameter = new GlobalParameter(volumeNodes.Count, boundaryNodes.Count);
                globalParameter.VolumeNodePointLocations = volumeNodes.ToArray();
                globalParameter.BoundaryNodePointLocations = boundaryNodes.ToArray();

                DA.SetData("GlobalParameter", globalParameter);
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
            get { return new Guid("280c58e4-bb79-45ea-bfe8-5e2f553d453b"); }
        }
    }
}