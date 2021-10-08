using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class Globals : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Globals class.
        /// </summary>
        public Globals()
          : base("Globals", "SetGlobalVaribles",
              "Set global varibles",
              "VolumeGeneratorBasedOnGraph", "CSV Import")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("VolumeNodePoints", "", "", GH_ParamAccess.list);
            pManager.AddPointParameter("AdditionalBoundaryNodePoints", "", "", GH_ParamAccess.list);
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
            get { return new Guid("280c58e4-bb79-45ea-bfe8-5e2f553d453b"); }
        }
    }
}