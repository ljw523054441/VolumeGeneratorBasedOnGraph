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
    public class GhcSiteDivision : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcSiteDivision class.
        /// </summary>
        public GhcSiteDivision()
          : base("GhcSiteDivision", "SiteDivision",
              "进行场地划分",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
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
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM))
            {
                DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);

                PlanktonMesh D = dualGraphWithHMDP.DualPlanktonMesh;

                // 找到边缘的面
                List<int> boundaryFaceIndex = new List<int>();
                for (int i = 0; i < D.Faces.Count; i++)
                {
                    if (D.Faces.NakedEdgeCount(i) >= 3)
                    {
                        boundaryFaceIndex.Add(i);
                    }
                }

                for (int i = 0; i < boundaryFaceIndex.Count; i++)
                {

                }
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
            get { return new Guid("3b884dcf-666e-4908-9a87-5e719ab579d8"); }
        }
    }
}