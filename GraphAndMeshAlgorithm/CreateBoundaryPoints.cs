using Grasshopper.Kernel;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Linq;

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
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddCurveParameter("BoundaryPolyline", "BPolyline", "场地边界的多段线", GH_ParamAccess.item);
            pManager.AddGenericParameter("DualHalfedgeMesh", "DHM", "生成的对偶图（半边数据结构）", GH_ParamAccess.item);
            pManager.AddIntegerParameter("faceIndexsFromOuterNodes", "", "", GH_ParamAccess.tree);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("OuterCorner", "OC", "在边界上的角点", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int volumeNodeCount = globalParameter.VolumeNodeCount;
            int boundaryNodeCount = globalParameter.BoundaryNodeCount;

            // todo: 在boundary多段线上，对各个角点的index重新排序，依据对偶图中对应点的序号

            Curve boundaryPolylineCurve = null;
            Polyline boundaryPolyline = null;

            if (DA.GetData("BoundaryPolyline", ref boundaryPolylineCurve))
            {
                if (!boundaryPolylineCurve.IsPlanar())
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的边界多段线不是平面的");
                    return;
                }
                // 对输入的Curve类型的boundaryPolyline进行类型转换，转换成Curve类的子类Polyline
                boundaryPolylineCurve.TryGetPolyline(out boundaryPolyline);
                

                if (boundaryPolyline.Count < 3)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "边界多段线至少要有三个顶点");
                    return;
                }
                
                Line[] boundaryPolylineSegments = boundaryPolyline.GetSegments();
                if (boundaryPolylineSegments.Length != boundaryNodeCount)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "边界多段线的段数一定要跟BoundaryNode的数量一样");
                    return;
                }

                Vector3d vector01 = new Vector3d(boundaryPolyline[1] - boundaryPolyline[0]);
                Vector3d vector02 = new Vector3d(boundaryPolyline[2] - boundaryPolyline[0]);
                Vector3d vectorZ = Vector3d.CrossProduct(vector01, vector02);
                /* 右手定则
                 * 如果叉积向外，则原来Polyline中点的排序是逆时针的，我们需要逆时针的排序
                 * 如果叉积向内，则原来Polyline中点的排序是顺时针的
                 */
                if (vectorZ.Z < 0)
                {
                    boundaryPolyline.Reverse();
                }

                List<Vector3d> segmentVerticalDirections = new List<Vector3d>();
                for (int i = 0; i < boundaryPolylineSegments.Length; i++)
                {
                    segmentVerticalDirections.Add(Vector3d.CrossProduct(boundaryPolylineSegments[i].Direction, Vector3d.ZAxis));
                }

                List<int> northSegmentIndexes = new List<int>();
                List<int> westSegmentIndexes = new List<int>();
                List<int> southSegmentIndexes = new List<int>();
                List<int> eastSegmentIndexes = new List<int>();
                for (int i = 0; i < boundaryPolylineSegments.Length; i++)
                {
                    if (segmentVerticalDirections[i] * Vector3d.YAxis > 0)
                    {
                        northSegmentIndexes.Add(i);
                    }
                    else if (segmentVerticalDirections[i] * (-1 * Vector3d.XAxis) > 0)
                    {
                        westSegmentIndexes.Add(i);
                    }
                    else if (segmentVerticalDirections[i] * (-1 * Vector3d.YAxis) > 0)
                    {
                        southSegmentIndexes.Add(i);
                    }
                    else
                    {
                        eastSegmentIndexes.Add(i);
                    }
                }

                List<Line> northSegments = new List<Line>();
                List<Line> westSegments = new List<Line>();
                List<Line> southSegments = new List<Line>();
                List<Line> eastSegments = new List<Line>();
                for (int i = 0; i < northSegmentIndexes.Count; i++)
                {
                    northSegments.Add(boundaryPolylineSegments[northSegmentIndexes[i]]);
                }
                for (int i = 0; i < westSegmentIndexes.Count; i++)
                {
                    westSegments.Add(boundaryPolylineSegments[westSegmentIndexes[i]]);
                }
                for (int i = 0; i < southSegmentIndexes.Count; i++)
                {
                    southSegments.Add(boundaryPolylineSegments[southSegmentIndexes[i]]);
                }
                for (int i = 0; i < eastSegmentIndexes.Count; i++)
                {
                    eastSegments.Add(boundaryPolylineSegments[eastSegmentIndexes[i]]);
                }

                List<Line> sortedNorthSegments = UtilityFunctions.SortSameOrientationSegments(northSegments, "N");
                List<Line> sortedWestSegments = UtilityFunctions.SortSameOrientationSegments(westSegments, "W");
                List<Line> sortedSouthSegments = UtilityFunctions.SortSameOrientationSegments(southSegments, "S");
                List<Line> sortedEastSegments = UtilityFunctions.SortSameOrientationSegments(eastSegments, "E");


                /* todo：
                 * 把按照逆时针规则排序好的线段的端点，附上根据faceIndexsFromOuterNodes输出的dual图顶点序号，形成新的列表。
                 */


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
            get { return new Guid("0224cb05-f2b2-46ea-acfe-ce1ab488c585"); }
        }
    }
}