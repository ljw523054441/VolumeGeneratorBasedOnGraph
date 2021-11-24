using Grasshopper.Kernel;
using Plankton;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcSplitTest : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcSplitTest class.
        /// </summary>
        public GhcSplitTest()
          : base("GhcSplitTest", "Nickname",
              "Description",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
            pManager.AddIntegerParameter("first", "", "", GH_ParamAccess.item);
            pManager.AddIntegerParameter("second", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
        }


        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int first = -1;
            int second = -1;

            DA.GetData("first", ref first);
            DA.GetData("second", ref second);
            
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM))
            {
                DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);

                dualGraphWithHMDP.DualPlanktonMesh.Vertices.SplitVertex(first, second);

                List<string> debug = UtilityFunctions.PrintFacesVertices(dualGraphWithHMDP.DualPlanktonMesh);

                DA.SetData("DualGraphWithHM", dualGraphWithHMDP);
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
            get { return new Guid("c68dae24-7b58-4019-9662-f0010603cce0"); }
        }
    }
}