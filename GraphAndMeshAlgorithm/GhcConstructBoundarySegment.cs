using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcConstructBoundarySegment : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcConstructBoundarySegment class.
        /// </summary>
        public GhcConstructBoundarySegment()
          : base("ConstructBoundarySegment", "ConstructBoundarySegment",
              "根据自己的设定，构造边界polyline的segment",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            base.Message = "开发中";
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("SelectedBoundaryPolylineSegment", "S", "被选的边界多边形的Segment", GH_ParamAccess.item);
            pManager.AddTextParameter("SegmentLabel", "L", "被选的边界多边形的Segment的标签名", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("BoundarySegment", "S", "构造生成的BoundarySegment对象", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Line line = new Line();
            string label = null;
            
            if (DA.GetData<Line>("SelectedBoundaryPolylineSegment", ref line)
                && DA.GetData<string>("SegmentLabel", ref label))
            {
                BoundarySegment boundarySegment = new BoundarySegment(line, label);

                DA.SetData("BoundarySegment", boundarySegment);
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
            get { return new Guid("ef42b4dc-b65f-4806-8786-a6d80cd5db5f"); }
        }
    }
}