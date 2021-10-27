using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcReorganizeBoundary : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateBoundaryPoints class.
        /// </summary>
        public GhcReorganizeBoundary()
          : base("ReorganizeBoundary", "ReorganizeBoundary",
              "依据对偶图的关系重新组织边界polyline",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            base.Message = "开发中 Todo:对偶图中，度为2的顶点布置方式";
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            // pManager.AddCurveParameter("BoundaryPolyline", "BPolyline", "场地边界的多段线", GH_ParamAccess.item);
            pManager.AddGenericParameter("BoundarySegments", "S", "构造生成的BoundarySegment对象列表", GH_ParamAccess.list);
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

            pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);
            pManager.AddIntegerParameter("ReorganizedCornerIndex", "RCIndex", "经过重新组织的边界polyline的角点对应的对偶图中角点的index", GH_ParamAccess.list);
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
            int innerNodeCount = globalParameter.VolumeNodeCount;
            int outerNodeCount = globalParameter.BoundaryNodeCount;

            // Curve boundaryPolylineCurve = null;
            List<BoundarySegment> boundarySegments = new List<BoundarySegment>();

            

            GH_Structure<GH_Integer> gh_structure_FaceIndexs = new GH_Structure<GH_Integer>();
            DataTree<int> faceIndexsFromOuterNodes = new DataTree<int>();

            List<Node> nodes = new List<Node>();

            if (DA.GetDataList<BoundarySegment>("BoundarySegments", boundarySegments)
                && DA.GetDataTree<GH_Integer>("faceIndexsFromOuterNodes", out gh_structure_FaceIndexs)
                && DA.GetDataList<Node>("GraphNode", nodes))
            {
                //if (!boundaryPolylineCurve.IsPlanar())
                //{
                //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的边界多段线不是平面的");
                //    return;
                //}
                //// 对输入的Curve类型的boundaryPolyline进行类型转换，转换成Curve类的子类Polyline
                //boundaryPolylineCurve.TryGetPolyline(out boundaryPolyline);

                
                
                /* 先对无序的boundarySegment进行排序
                 * 按照Label进行
                 * 排完后，调整From和To的顺序
                 */
                List<string> labels = new List<string>();
                for (int i = 0; i < nodes.Count; i++)
                {
                    labels.Add(nodes[i].NodeAttribute.NodeLabel);
                }

                // 按照Label进行排序
                List<BoundarySegment> sortedBoundarySegments = new List<BoundarySegment>();
                for (int i = 0; i < labels.Count; i++)
                {
                    for (int j = 0; j < boundarySegments.Count; j++)
                    {
                        if (boundarySegments[j].Label == labels[i])
                        {
                            sortedBoundarySegments.Add(boundarySegments[j]);
                        }
                    }
                }

                List<Point3d> boundaryCorner = new List<Point3d>();
                bool flag = UtilityFunctions.PointCompare(sortedBoundarySegments[0].From,
                                                          sortedBoundarySegments[0].To,
                                                          new Point3d((sortedBoundarySegments[0].From + sortedBoundarySegments[0].To) / 2));
                // 如果[0]在[1]的逆时针方向，那么首尾交换
                if (flag == false)
                {
                    sortedBoundarySegments[0].Reverse();
                }
                boundaryCorner.Add(sortedBoundarySegments[0].From);
                // 从第二个开始，调整from和to
                for (int i = 1; i < sortedBoundarySegments.Count; i++)
                {
                    if (sortedBoundarySegments[i].From != sortedBoundarySegments[i-1].To)
                    {
                        sortedBoundarySegments[i].Reverse();
                    }
                    boundaryCorner.Add(sortedBoundarySegments[i].From);
                }
                if (sortedBoundarySegments.Last().To != sortedBoundarySegments[0].From)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Boundary不闭合");
                }

                // 输出闭合的多段线
                List<Point3d> closedBoundaryPoints = boundaryCorner;
                closedBoundaryPoints.Add(boundaryCorner[0]);
                Polyline reorganizedBoundary = new Polyline(closedBoundaryPoints);
                DA.SetData("ReorganizedBoundary", reorganizedBoundary);



                // 计算ReorganizedCornerIndex的部分
                // 从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_structure_FaceIndexs, ref faceIndexsFromOuterNodes);
                // 去除-1的情况
                for (int i = 0; i < faceIndexsFromOuterNodes.BranchCount; i++)
                {
                    faceIndexsFromOuterNodes.Branch(i).RemoveAt(0);
                }

                List<int> boundaryCornerIndexs = new List<int>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    boundaryCornerIndexs.Add(faceIndexsFromOuterNodes.Branch(i).First());
                }
                DA.SetDataList("ReorganizedCornerIndex", boundaryCornerIndexs);





                //if (boundaryPolyline.Count < 3)
                //{
                //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "边界多段线至少要有三个顶点");
                //    return;
                //}

                //// 判断boundaryPolyline的绘制方向
                //Vector3d vector01 = new Vector3d(boundaryPolyline[1] - boundaryPolyline[0]);
                //Vector3d vector02 = new Vector3d(boundaryPolyline[2] - boundaryPolyline[0]);
                //Vector3d vectorZ = Vector3d.CrossProduct(vector01, vector02);
                ///* 叉积右手定则
                // * 如果叉积向外，则原来Polyline中点的排序是逆时针的，我们需要逆时针的排序
                // * 如果叉积向内，则原来Polyline中点的排序是顺时针的
                // */
                //if (vectorZ.Z < 0)
                //{
                //    boundaryPolyline.Reverse();
                //}

                //Line[] boundaryPolylineSegments = boundaryPolyline.GetSegments();
                //if (boundaryPolylineSegments.Length != outerNodeCount)
                //{
                //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "边界多段线的段数一定要跟BoundaryNode的数量一样");
                //    return;
                //}

                //List<Vector3d> segmentVerticalDirections = new List<Vector3d>();
                //for (int i = 0; i < boundaryPolylineSegments.Length; i++)
                //{
                //    segmentVerticalDirections.Add(Vector3d.CrossProduct(boundaryPolylineSegments[i].Direction, Vector3d.ZAxis));
                //}

                //List<int> northSegmentIndexes = new List<int>();
                //List<int> westSegmentIndexes = new List<int>();
                //List<int> southSegmentIndexes = new List<int>();
                //List<int> eastSegmentIndexes = new List<int>();
                //// int lastIndex = 0;
                //for (int i = 0; i < boundaryPolylineSegments.Length; i++)
                //{
                //    // 求segmentVerticalDirections与NEWS方向的夹角
                //    List<double> radians = new List<double>();
                //    // 与N方向的夹角
                //    radians.Add(Vector3d.VectorAngle(segmentVerticalDirections[i], Vector3d.YAxis, Plane.WorldXY));
                //    // 与W方向的夹角
                //    radians.Add(Vector3d.VectorAngle(segmentVerticalDirections[i], (-1 * Vector3d.XAxis), Plane.WorldXY));
                //    // 与S方向的夹角
                //    radians.Add(Vector3d.VectorAngle(segmentVerticalDirections[i], (-1 * Vector3d.YAxis), Plane.WorldXY));
                //    // 与E方向的夹角
                //    radians.Add(Vector3d.VectorAngle(segmentVerticalDirections[i], Vector3d.XAxis, Plane.WorldXY));

                //    // 找出最小的那个夹角对应的NEWS方向
                //    //if (lastIndex != )
                //    //{

                //    //}
                //    int index = radians.IndexOf(radians.Min());
                //    switch (index)
                //    {
                //        case 0:
                //            northSegmentIndexes.Add(i);
                //            break;
                //        case 1:
                //            westSegmentIndexes.Add(i);
                //            break;
                //        case 2:
                //            southSegmentIndexes.Add(i);
                //            break;
                //        case 3:
                //            eastSegmentIndexes.Add(i);
                //            break;
                //    }
                //    // lastIndex = index;

                //    //if (segmentVerticalDirections[i] * Vector3d.YAxis > 0)
                //    //{
                //    //    northSegmentIndexes.Add(i);
                //    //}
                //    //else if (segmentVerticalDirections[i] * (-1 * Vector3d.XAxis) > 0)
                //    //{
                //    //    westSegmentIndexes.Add(i);
                //    //}
                //    //else if (segmentVerticalDirections[i] * (-1 * Vector3d.YAxis) > 0)
                //    //{
                //    //    southSegmentIndexes.Add(i);
                //    //}
                //    //else
                //    //{
                //    //    eastSegmentIndexes.Add(i);
                //    //}
                //}

                //List<Line> northSegments = new List<Line>();
                //List<Line> westSegments = new List<Line>();
                //List<Line> southSegments = new List<Line>();
                //List<Line> eastSegments = new List<Line>();
                //for (int i = 0; i < northSegmentIndexes.Count; i++)
                //{
                //    northSegments.Add(boundaryPolylineSegments[northSegmentIndexes[i]]);
                //}
                //for (int i = 0; i < westSegmentIndexes.Count; i++)
                //{
                //    westSegments.Add(boundaryPolylineSegments[westSegmentIndexes[i]]);
                //}
                //for (int i = 0; i < southSegmentIndexes.Count; i++)
                //{
                //    southSegments.Add(boundaryPolylineSegments[southSegmentIndexes[i]]);
                //}
                //for (int i = 0; i < eastSegmentIndexes.Count; i++)
                //{
                //    eastSegments.Add(boundaryPolylineSegments[eastSegmentIndexes[i]]);
                //}

                //List<Line> sortedNorthSegments = UtilityFunctions.SortSameOrientationSegments(northSegments, "N");
                //List<Line> sortedWestSegments = UtilityFunctions.SortSameOrientationSegments(westSegments, "W");
                //List<Line> sortedSouthSegments = UtilityFunctions.SortSameOrientationSegments(southSegments, "S");
                //List<Line> sortedEastSegments = UtilityFunctions.SortSameOrientationSegments(eastSegments, "E");


                ///* todo：
                // * 把按照逆时针规则排序好的线段的端点，附上根据faceIndexsFromOuterNodes输出的dual图顶点序号，形成新的列表。
                // */
                //// todo: 在boundary多段线上，对各个角点的index重新排序，依据对偶图中对应点的序号


                //// 从GH_Structure<GH_Integer>转化为DataTree<int>
                //UtilityFunctions.GH_StructureToDataTree_Int(gh_structure_FaceIndexs, ref faceIndexsFromOuterNodes);
                //// 去除-1的情况
                //for (int i = 0; i < faceIndexsFromOuterNodes.BranchCount; i++)
                //{
                //    faceIndexsFromOuterNodes.Branch(i).RemoveAt(0);
                //}

                //List<string> outerNodeLabels = new List<string>();
                //for (int i = 0; i < outerNodeCount; i++)
                //{
                //    outerNodeLabels.Add(nodes[i + innerNodeCount].NodeAttribute.NodeLabel);
                //}

                //// 将outerNode中的label按照NEWS进行分类
                //List<string> nLabels = new List<string>();
                //List<string> wLabels = new List<string>();
                //List<string> sLabels = new List<string>();
                //List<string> eLabels = new List<string>();
                //List<int> nLabelIndexs = new List<int>();
                //List<int> wLabelIndexs = new List<int>();
                //List<int> sLabelIndexs = new List<int>();
                //List<int> eLabelIndexs = new List<int>();
                //for (int i = 0; i < outerNodeLabels.Count; i++)
                //{
                //    if (outerNodeLabels[i].Contains('N'))
                //    {
                //        nLabels.Add(outerNodeLabels[i]);
                //        nLabelIndexs.Add(i);
                //    }
                //    else if (outerNodeLabels[i].Contains('W'))
                //    {
                //        wLabels.Add(outerNodeLabels[i]);
                //        wLabelIndexs.Add(i);
                //    }
                //    else if (outerNodeLabels[i].Contains('S'))
                //    {
                //        sLabels.Add(outerNodeLabels[i]);
                //        sLabelIndexs.Add(i);
                //    }
                //    else if (outerNodeLabels[i].Contains('E'))
                //    {
                //        eLabels.Add(outerNodeLabels[i]);
                //        eLabelIndexs.Add(i);
                //    }
                //    else
                //    {
                //        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "OuterNode的标签中有除了NEWS之外的其他字符");
                //    }
                //}


                //List<Point3d[]> nBoundarySegmentStartAndEndPoints = new List<Point3d[]>();
                //List<int[]> nBoundarySegmentStartAndEndPointsIndexs = new List<int[]>();
                //// 对于N相关的边
                //for (int i = 0; i < sortedNorthSegments.Count; i++)
                //{
                //    Point3d[] pointsTuple = new Point3d[] { sortedNorthSegments[i].From, sortedNorthSegments[i].To };
                //    int[] indexsTuple = new int[] { faceIndexsFromOuterNodes.Branch(nLabelIndexs[i]).First(), 
                //                                    faceIndexsFromOuterNodes.Branch(nLabelIndexs[i]).Last() };

                //    nBoundarySegmentStartAndEndPoints.Add(pointsTuple);
                //    nBoundarySegmentStartAndEndPointsIndexs.Add(indexsTuple);
                //}

                //List<Point3d[]> wBoundarySegmentStartAndEndPoints = new List<Point3d[]>();
                //List<int[]> wBoundarySegmentStartAndEndPointsIndexs = new List<int[]>();
                //// 对于W相关的边
                //for (int i = 0; i < sortedWestSegments.Count; i++)
                //{
                //    Point3d[] pointsTuple = new Point3d[] { sortedWestSegments[i].From, sortedWestSegments[i].To };
                //    int[] indexsTuple = new int[] { faceIndexsFromOuterNodes.Branch(wLabelIndexs[i]).First(),
                //                                    faceIndexsFromOuterNodes.Branch(wLabelIndexs[i]).Last() };

                //    wBoundarySegmentStartAndEndPoints.Add(pointsTuple);
                //    wBoundarySegmentStartAndEndPointsIndexs.Add(indexsTuple);
                //}

                //List<Point3d[]> sBoundarySegmentStartAndEndPoints = new List<Point3d[]>();
                //List<int[]> sBoundarySegmentStartAndEndPointsIndexs = new List<int[]>();
                //// 对于S相关的边
                //for (int i = 0; i < sortedSouthSegments.Count; i++)
                //{
                //    Point3d[] pointsTuple = new Point3d[] { sortedSouthSegments[i].From, sortedSouthSegments[i].To };
                //    int[] indexsTuple = new int[] { faceIndexsFromOuterNodes.Branch(sLabelIndexs[i]).First(),
                //                                    faceIndexsFromOuterNodes.Branch(sLabelIndexs[i]).Last() };

                //    sBoundarySegmentStartAndEndPoints.Add(pointsTuple);
                //    sBoundarySegmentStartAndEndPointsIndexs.Add(indexsTuple);
                //}

                //List<Point3d[]> eBoundarySegmentStartAndEndPoints = new List<Point3d[]>();
                //List<int[]> eBoundarySegmentStartAndEndPointsIndexs = new List<int[]>();
                //// 对于E相关的边
                //for (int i = 0; i < sortedEastSegments.Count; i++)
                //{
                //    Point3d[] pointsTuple = new Point3d[] { sortedEastSegments[i].From, sortedEastSegments[i].To };
                //    int[] indexsTuple = new int[] { faceIndexsFromOuterNodes.Branch(eLabelIndexs[i]).First(),
                //                                    faceIndexsFromOuterNodes.Branch(eLabelIndexs[i]).Last() };

                //    eBoundarySegmentStartAndEndPoints.Add(pointsTuple);
                //    eBoundarySegmentStartAndEndPointsIndexs.Add(indexsTuple);
                //}

                //List<Point3d> boundaryPoints = new List<Point3d>();
                //List<int> boundaryPointIndexs = new List<int>();

                ///* 从W开始，W1,W2...
                // * 把W1的起始点放入boundaryPoints
                // * 把W1的起始点对应的在faceIndexsFromOuterNodes中的序号放入boundaryPointIndexs
                // * ...
                // */
                //for (int i = 0; i < wBoundarySegmentStartAndEndPoints.Count; i++)
                //{
                //    boundaryPoints.Add(wBoundarySegmentStartAndEndPoints[i][0]);
                //    boundaryPointIndexs.Add(wBoundarySegmentStartAndEndPointsIndexs[i][0]);
                //}
                ///* 继续添加S
                // * 把S1的起始点放入boundaryPoints
                // * 把S1的起始点对应的在faceIndexsFromOuterNodes中的序号放入boundaryPointIndexs
                // * ...
                // */
                //for (int i = 0; i < sBoundarySegmentStartAndEndPoints.Count; i++)
                //{
                //    boundaryPoints.Add(sBoundarySegmentStartAndEndPoints[i][0]);
                //    boundaryPointIndexs.Add(sBoundarySegmentStartAndEndPointsIndexs[i][0]);
                //}
                ///* 继续添加E
                // * 把E1的起始点放入boundaryPoints
                // * 把E1的起始点对应的在faceIndexsFromOuterNodes中的序号放入boundaryPointIndexs
                // * ...
                // */
                //for (int i = 0; i < eBoundarySegmentStartAndEndPoints.Count; i++)
                //{
                //    boundaryPoints.Add(eBoundarySegmentStartAndEndPoints[i][0]);
                //    boundaryPointIndexs.Add(eBoundarySegmentStartAndEndPointsIndexs[i][0]);
                //}
                ///* 继续添加N
                // * 把N1的起始点放入boundaryPoints
                // * 把N1的起始点对应的在faceIndexsFromOuterNodes中的序号放入boundaryPointIndexs
                // * ...
                // */
                //for (int i = 0; i < nBoundarySegmentStartAndEndPoints.Count; i++)
                //{
                //    boundaryPoints.Add(nBoundarySegmentStartAndEndPoints[i][0]);
                //    boundaryPointIndexs.Add(nBoundarySegmentStartAndEndPointsIndexs[i][0]);
                //}

                //List<Point3d> closedBoundaryPoints = boundaryPoints;
                //closedBoundaryPoints.Add(boundaryPoints.First());

                //Polyline reorganizedBoundary = new Polyline(closedBoundaryPoints);
                //DA.SetData("ReorganizedBoundary", reorganizedBoundary);

                //DA.SetDataList("ReorganizedCornerIndex", boundaryPointIndexs);
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