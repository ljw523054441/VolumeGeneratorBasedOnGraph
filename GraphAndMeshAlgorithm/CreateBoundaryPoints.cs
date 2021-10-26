using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class CreateBoundaryPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateBoundaryPoints class.
        /// </summary>
        public CreateBoundaryPoints()
          : base("CreateBoundaryPoints", "CreateBoundaryPoints",
              "依据对偶图的关系创建边界点",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualHalfedgeMesh", "DHM", "生成的对偶图（半边数据结构）", GH_ParamAccess.item);
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
            // 在boundary多段线上，对各个角点的index重新排序，依据对偶图中对应点的序号
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
            get { return new Guid("0224cb05-f2b2-46ea-acfe-ce1ab488c585"); }
        }
    }
}