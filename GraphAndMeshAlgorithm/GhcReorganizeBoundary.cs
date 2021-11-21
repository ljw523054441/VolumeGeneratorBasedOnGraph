using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

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
            base.Message = "开发中 Todo:对偶图中，边界顶点布置方式";

            BoundaryPolylinePoints = new List<Point3d>();
            BoundaryCornerTextDots = new List<TextDot>();
            BoundarySegmentTextDots = new List<TextDot>();

            VolumeJunctionTextDot = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> BoundaryPolylinePoints;

        private List<TextDot> BoundaryCornerTextDots;
        private List<TextDot> BoundarySegmentTextDots;

        private List<TextDot> VolumeJunctionTextDot;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);

            // pManager.AddCurveParameter("BoundaryPolyline", "BPolyline", "场地边界的多段线", GH_ParamAccess.item);
            pManager.AddGenericParameter("BoundarySegments", "S", "构造生成的BoundarySegment对象列表", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DualHalfedgeMesh", "DHM", "生成的对偶图（半边数据结构）", GH_ParamAccess.item);
            // pManager.AddIntegerParameter("faceIndexsFromOuterNodes", "", "", GH_ParamAccess.tree);

            // pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("BoundaryCorner", "BC", "在边界上的角点", GH_ParamAccess.list);

            pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);
            pManager.AddIntegerParameter("ReorganizedCornerIndex", "RCIndex", "经过重新组织的边界polyline的角点对应的对偶图中角点的index", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 局部变量初始化
            Thickness = 2;
            List<BoundarySegment> boundarySegments = new List<BoundarySegment>();

            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();

            #endregion

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetDataList<BoundarySegment>("BoundarySegments", boundarySegments))
            {
                PlanktonMesh D = dualGraphWithHM.DualPlanktonMesh;
                PlanktonMesh I = dualGraphWithHM.IntegrateDualPlanktonMesh;
                List<GraphNode> decomposedPGraphNodes = dualGraphWithHM.GraphNodes;

                SortedDictionary<int, GraphNode> volume_volumeNode = dualGraphWithHM.Volume_VolumeNode;

                List<GraphNode> undividedPGraphNodes = dualGraphWithHM.UndividedGraphNodes;
                List<List<int>> undividedPGraphTables = dualGraphWithHM.UndividedGraphTable;



                int innerNodeCount = dualGraphWithHM.InnerNodeCount;
                int outerNodeCount = dualGraphWithHM.OuterNodeCount;

                /* 先对无序的boundarySegment进行排序
                 * 按照Label进行
                 * 排完后，调整From和To的顺序
                 */
                #region 获取outerNodeLabel列表
                List<string> outerNodeLabels = new List<string>();
                //for (int i = 0; i < nodes.Count - innerNodeCount; i++)
                //{
                //    outerLabels.Add(nodes[i + innerNodeCount].NodeAttribute.NodeLabel);
                //}
                for (int i = 0; i < decomposedPGraphNodes.Count; i++)
                {
                    if (!decomposedPGraphNodes[i].IsInner)
                    {
                        outerNodeLabels.Add(decomposedPGraphNodes[i].NodeAttribute.NodeLabel);
                    }
                }
                #endregion
                #region 按照Label进行排序：把输入的boundarySegement按照outerNodeLabels的顺序来排序
                List<BoundarySegment> sortedBoundarySegments = new List<BoundarySegment>();
                for (int i = 0; i < outerNodeLabels.Count; i++)
                {
                    for (int j = 0; j < boundarySegments.Count; j++)
                    {
                        if (boundarySegments[j].Label == outerNodeLabels[i])
                        {
                            sortedBoundarySegments.Add(boundarySegments[j]);
                        }
                    }
                }
                #endregion

                #region 按照从W开始，逆时针的顺序，输出场地边界的角点
                List<Point3d> boundaryCorner = new List<Point3d>();

                // 利用叉积来判断第一个Segment的方向时候需要反转
                if (IsFirstSegmentOpppsite(sortedBoundarySegments[0].Line.Direction,
                                           new Vector3d((sortedBoundarySegments[1].From + sortedBoundarySegments[1].To) / 2 - sortedBoundarySegments[0].From)))
                {
                    sortedBoundarySegments[0].Reverse();
                }

                boundaryCorner.Add(sortedBoundarySegments[0].From);
                // 从第二个开始，调整from和to
                for (int i = 1; i < sortedBoundarySegments.Count; i++)
                {
                    if (sortedBoundarySegments[i].From != sortedBoundarySegments[i - 1].To)
                    {
                        sortedBoundarySegments[i].Reverse();
                    }
                    boundaryCorner.Add(sortedBoundarySegments[i].From);
                }
                // 判断生成的polyline是否闭合
                if (sortedBoundarySegments.Last().To != sortedBoundarySegments[0].From)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Boundary不闭合");
                }
                #endregion
                DA.SetDataList("BoundaryCorner", boundaryCorner);

                #region 输出闭合的多段线
                List<Point3d> closedBoundaryPoints = new List<Point3d>();
                closedBoundaryPoints.AddRange(boundaryCorner);
                closedBoundaryPoints.Add(boundaryCorner[0]);
                Polyline reorganizedBoundary = new Polyline(closedBoundaryPoints);
                #endregion
                DA.SetData("ReorganizedBoundary", reorganizedBoundary);



                #region 计算场地边界与角点所对应的对偶图中顶点的序号

                #region 输出从W开始的逆时针方向的，场地边界每个角点所对应的对偶图顶点的序号
                // 从GH_Structure<GH_Integer>转化为DataTree<int>
                // UtilityFunctions.GH_StructureToDataTree_Int(gh_structure_FaceIndexs, ref faceIndexsAroundOuterNodes);
                List<List<int>> faceIndexsAroundOuterNodes = dualGraphWithHM.DFaceIndexsAroundOuterNodes;
                // 去除-1的情况，第0个永远是-1
                for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
                {
                    faceIndexsAroundOuterNodes[i].RemoveAt(0);
                }
                List<int> sortedBoundaryCornerIndexs = new List<int>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    sortedBoundaryCornerIndexs.Add(faceIndexsAroundOuterNodes[i].First());
                }
                #endregion
                DA.SetDataList("ReorganizedCornerIndex", sortedBoundaryCornerIndexs);

                #region 输出从W开始的逆时针方向的，场地边界上所应该布置的对偶图顶点序号
                List<List<int>> verticesIndexForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    verticesIndexForEachBoundarySegment.Add(new List<int>());
                    verticesIndexForEachBoundarySegment[i].AddRange(faceIndexsAroundOuterNodes[i]);
                }
                #endregion

                #endregion

                List<List<int>> faceIndexsAroundVertex = new List<List<int>>();
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(D.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(0);
                }

                //#region 计算每个segment上的中间点的t值
                //List<List<double>> tOnSegment = new List<List<double>>();
                //for (int i = 0; i < verticesIndexForEachBoundarySegment.Count; i++)
                //{
                //    if (verticesIndexForEachBoundarySegment[i].Count == 2)
                //    {
                //        tOnSegment.Add(new List<double>());
                //    }
                //    else
                //    {
                //        List<int> centralVertexIndex = verticesIndexForEachBoundarySegment[i].GetRange(1, verticesIndexForEachBoundarySegment[i].Count - 2);


                //    }
                //}

                //#endregion

                SortedDictionary<int, List<int>> volumeJunction = new SortedDictionary<int, List<int>>();
                foreach (var item in dualGraphWithHM.DVertice_Volume)
                {
                    if (item.Value.Count > 1 && I.Vertices.GetVertexFaces(item.Key).Contains(-1))
                    {
                        volumeJunction.Add(item.Key, item.Value);
                    }
                }


                // 
                SortedDictionary<int, int> segment_volumeJunction = new SortedDictionary<int, int>();
                // 求volumeJunction点的t值
                SortedDictionary<int, double> volumeJunction_t = CalVolumeJunction_T(volumeJunction, verticesIndexForEachBoundarySegment, dualGraphWithHM, out segment_volumeJunction);

                // SortedDictionary<int, Point3d> segment_volumeJunctionPoint = new SortedDictionary<int, Point3d>();
                List<int> volumeJunctionIndex = new List<int>();
                List<Point3d> volumeJunctionPoint = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    foreach (var item in segment_volumeJunction)
                    {
                        if (i == item.Key)
                        {
                            volumeJunctionIndex.Add(item.Value);
                            double t = volumeJunction_t[item.Value];
                            volumeJunctionPoint.Add(sortedBoundarySegments[i].Line.PointAt(t));
                        }
                    }
                }




                #region 可视化部分
                BoundaryPolylinePoints.Clear();
                BoundaryCornerTextDots.Clear();
                BoundarySegmentTextDots.Clear();

                VolumeJunctionTextDot.Clear();

                #region 设置表示场地边界的Polyline
                BoundaryPolylinePoints = closedBoundaryPoints;
                #endregion

                #region 设置场地角点的TextDot
                List<string> args = new List<string>();
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    args.Add(string.Join<int>(";", D.Vertices.GetVertexFaces(i)));
                }
                for (int i = 0; i < boundaryCorner.Count; i++)
                {
                    TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], args[sortedBoundaryCornerIndexs[i]]), boundaryCorner[i]);
                    BoundaryCornerTextDots.Add(boundaryCornerTextDot);
                }
                #endregion

                #region 设置volumeJunction
                // List<string> args2 = new List<string>();
                for (int i = 0; i < volumeJunctionIndex.Count; i++)
                {
                    List<int> volumeIndex = volumeJunction[volumeJunctionIndex[i]];
                    string str = string.Join<int>(";", volumeIndex);
                    TextDot volumeJunctionTextDot = new TextDot(string.Format("{0} | {1}", volumeJunctionIndex[i], str), volumeJunctionPoint[i]);
                    VolumeJunctionTextDot.Add(volumeJunctionTextDot);
                }
                #endregion



                #region 设置场地边界Segment的TextDot
                List<Point3d> boundarySegmentTextDotLocations = new List<Point3d>();
                for (int i = 0; i < boundarySegments.Count; i++)
                {
                    Vector3d verticalVector = Vector3d.CrossProduct(boundarySegments[i].Line.Direction, Vector3d.ZAxis);
                    verticalVector.Unitize();
                    Point3d boundarySegmentTextDotLocation = new Point3d((boundarySegments[i].From + boundarySegments[i].To) / 2);
                    boundarySegmentTextDotLocation += verticalVector * 2;
                    boundarySegmentTextDotLocations.Add(boundarySegmentTextDotLocation);
                }
                for (int i = 0; i < boundarySegments.Count; i++)
                {
                    TextDot boundarySegmentTextDot = new TextDot(string.Format("{0} | {1}", i + innerNodeCount, outerNodeLabels[i]), boundarySegmentTextDotLocations[i]);
                    BoundarySegmentTextDots.Add(boundarySegmentTextDot);
                }
                #endregion

                #endregion



                
            }

        }

        /// <summary>
        /// 用叉积来判断第一段Segment是否反了
        /// </summary>
        /// <param name="lineDirection"></param>
        /// <param name="vector02"></param>
        /// <returns></returns>
        private bool IsFirstSegmentOpppsite(Vector3d lineDirection, Vector3d vector02)
        {
            Vector3d Z = Vector3d.CrossProduct(lineDirection, vector02);
            if (Z.Z > 0)
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        public SortedDictionary<int, double> CalVolumeJunction_T(SortedDictionary<int, List<int>> volumeJunction, 
                                                                 List<List<int>> verticesIndexForEachBoundarySegment,
                                                                 DualGraphWithHM dualGraphWithHM,
                                                                 out SortedDictionary<int,int> segment_volumeJunction)
        {
            SortedDictionary<int, double> volumeJunction_t = new SortedDictionary<int, double>();
            segment_volumeJunction = new SortedDictionary<int, int>();
            foreach (var item in volumeJunction)
            {
                List<int> verticesIndexForCurrentBoundarySegment = new List<int>();
                for (int i = 0; i < verticesIndexForEachBoundarySegment.Count; i++)
                {
                    if (verticesIndexForEachBoundarySegment[i].Contains(item.Key))
                    {
                        verticesIndexForCurrentBoundarySegment.AddRange(verticesIndexForEachBoundarySegment[i]);
                        segment_volumeJunction.Add(i, item.Key);
                        break;
                    }
                }

                // 对于分割两个或几个Volume的边缘点来说
                switch (item.Value.Count)
                {
                    // 如果分割2个volume
                    case 2:
                        int vface1 = -1;
                        int vface2 = -1;
                        for (int i = 0; i < dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.Count; i++)
                        {
                            if (dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[i].StartVertex == item.Key
                                && dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.EndVertex(i) == verticesIndexForCurrentBoundarySegment[verticesIndexForCurrentBoundarySegment.IndexOf(item.Key) - 1])
                            {
                                vface1 = dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[i].AdjacentFace;
                            }
                            if (dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[i].StartVertex == verticesIndexForCurrentBoundarySegment[verticesIndexForCurrentBoundarySegment.IndexOf(item.Key) + 1]
                                && dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.EndVertex(i) == item.Key)
                            {
                                vface2 = dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[i].AdjacentFace;
                            }
                        }

                        int volume1 = -1;
                        int volume2 = -1;
                        foreach (var item1 in dualGraphWithHM.Volume_VFace)
                        {
                            if (vface1 == item1.Value)
                            {
                                volume1 = item1.Key;
                            }
                            if (vface2 == item1.Value)
                            {
                                volume2 = item1.Key;
                            }
                        }

                        double volume1AreaProportion = dualGraphWithHM.UndividedGraphNodes[volume1].NodeAttribute.NodeAreaProportion;
                        double volume2AreaProportion = dualGraphWithHM.UndividedGraphNodes[volume2].NodeAttribute.NodeAreaProportion;

                        double t = volume1AreaProportion / (volume1AreaProportion + volume2AreaProportion);

                        volumeJunction_t.Add(item.Key, t);

                        break;
                    // 如果分割3个volume
                    case 3:
                        break;
                    default:
                        break;
                }
            }

            return volumeJunction_t;
        }




        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);

            args.Display.DrawPolyline(BoundaryPolylinePoints, Color.DarkOrange, Thickness);

            for (int i = 0; i < BoundaryCornerTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(BoundaryCornerTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < BoundarySegmentTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(BoundarySegmentTextDots[i], Color.Gray, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < VolumeJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // base.DrawViewportMeshes(args);

            args.Display.DrawPolyline(BoundaryPolylinePoints, Color.DarkOrange, Thickness);

            for (int i = 0; i < BoundaryCornerTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(BoundaryCornerTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < BoundarySegmentTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(BoundarySegmentTextDots[i], Color.Gray, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < VolumeJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
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