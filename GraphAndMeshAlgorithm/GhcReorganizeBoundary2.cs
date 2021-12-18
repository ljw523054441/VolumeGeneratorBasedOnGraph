using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcReorganizeBoundary2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcReorganizeBoundary2 class.
        /// </summary>
        public GhcReorganizeBoundary2()
          : base("ReorganizeBoundary2", "ReorganizeBoundary2",
              "依据对偶图的关系重新组织边界polyline",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            BoundaryPolylinePoints = new List<Point3d>();
            BoundaryCornerTextDots = new List<TextDot>();
            BoundarySegmentTextDots = new List<TextDot>();

        }

        private int Thickness;

        private List<Point3d> BoundaryPolylinePoints;

        private List<TextDot> BoundaryCornerTextDots;
        private List<TextDot> BoundarySegmentTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
            pManager.AddGenericParameter("BoundarySegments", "S", "构造生成的BoundarySegment对象列表", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);
            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);
            pManager.AddIntegerParameter("BSIndexContainVolumeJunctions", "BSI", "包含边界分裂点的BoundarySegment序号", GH_ParamAccess.list);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);

            pManager.AddBooleanParameter("NeedToConnect", "NTC", "的这一对分裂点是否作为分界线", GH_ParamAccess.list);
            
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
                DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);
                List<GraphNode> decomposedPGraphNodes = dualGraphWithHMDP.GraphNodes;

                int innerNodeCount = dualGraphWithHMDP.InnerNodeCount;
                int outerNodeCount = dualGraphWithHMDP.OuterNodeCount;

                /* 先对无序的boundarySegment进行排序
                 * 按照Label进行
                 * 排完后，调整From和To的顺序
                 */
                #region 获取outerNodeLabel列表
                List<string> outerNodeLabels = new List<string>();
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
                DA.SetDataList("SortedBoundarySegments", sortedBoundarySegments);

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

                #region 输出经过重新组织后，闭合的多段线
                List<Point3d> closedBoundaryPoints = new List<Point3d>();
                closedBoundaryPoints.AddRange(boundaryCorner);
                closedBoundaryPoints.Add(boundaryCorner[0]);
                Polyline reorganizedBoundary = new Polyline(closedBoundaryPoints);
                #endregion
                DA.SetData("ReorganizedBoundary", reorganizedBoundary);

                #region 计算场地边界与角点所对应的对偶图中顶点的序号

                #region 输出从W开始的逆时针方向的，场地边界每个角点所对应的对偶图顶点的序号
                List<List<int>> faceIndexsAroundOuterNodes = dualGraphWithHMDP.DFaceIndexsAroundOuterNodes;
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
                // DA.SetDataList("ReorganizedCornerIndex", sortedBoundaryCornerIndexs);
                #endregion

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
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(-1);
                }

                #region 找到每个segment上的volumeJunction
                List<List<int>> volumeJunctionForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    volumeJunctionForEachBoundarySegment.Add(new List<int>());
                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        if (faceIndexsAroundVertex[verticesIndexForEachBoundarySegment[i][j]].Count > 1
                            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBoundarySegment[i][j]).Contains(-1))
                        {
                            volumeJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                        }
                    }
                }
                #endregion

                #region 找到包含边界分裂点的BoundarySegment序号
                List<int> bSIndexContainVolumeJunctions = new List<int>();
                List<int> volumeJunctionIndex = new List<int>();
                for (int i = 0; i < volumeJunctionForEachBoundarySegment.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionForEachBoundarySegment[i].Count; j++)
                    {
                        bSIndexContainVolumeJunctions.Add(i);
                        volumeJunctionIndex.Add(volumeJunctionForEachBoundarySegment[i][j]);
                    }
                }
                #endregion
                DA.SetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions);

                #region 输出表示边界分裂点的Text
                List<string> volumeJunctionText = new List<string>();
                for (int i = 0; i < volumeJunctionForEachBoundarySegment.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionForEachBoundarySegment[i].Count; j++)
                    {
                        string arg = string.Format("{0} | {1}", volumeJunctionForEachBoundarySegment[i][j], string.Join<int>(";", faceIndexsAroundVertex[verticesIndexForEachBoundarySegment[i][j]]));
                        volumeJunctionText.Add(arg);
                    }
                }
                #endregion
                DA.SetDataList("VolumeJunctionsTexts", volumeJunctionText);

                #region 得到相邻的两个VolumeJunction是否构成分界线的判断
                List<int[]> pairVolumeJunctionsIndex = new List<int[]>();

                List<int[]> indexPair = new List<int[]>();
                for (int i = 0; i < volumeJunctionIndex.Count; i++)
                {
                    int[] pair = new int[2] { volumeJunctionIndex[i], volumeJunctionIndex[(i + 1) % volumeJunctionIndex.Count] };
                    indexPair.Add(pair);
                }

                List<bool> needToConnect = new List<bool>();
                foreach (int[] pair in indexPair)
                {
                    int faceIndex = -1;
                    int[] hs = null;
                    for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Faces.Count; i++)
                    {
                        int[] faceVertex = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetFaceVertices(i);
                        int[] intersect = faceVertex.Intersect(pair).ToArray();
                        if (intersect.Length == pair.Length)
                        {
                            faceIndex = i;
                            hs = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetHalfedges(i);
                        }
                    }

                    if (faceIndex == -1)
                    {
                        needToConnect.Add(false);
                        continue;
                    }

                    int hIndex = -1;
                    for (int i = 0; i < hs.Length; i++)
                    {
                        if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[hs[i]].StartVertex == pair[0])
                        {
                            hIndex = hs[i];
                        }
                    }

                    int currH = hIndex;
                    List<int> path = new List<int>();
                    while (dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(currH) != pair[1])
                    {
                        path.Add(currH);
                        currH = dualGraphWithHMDP.DualPlanktonMesh.Halfedges[currH].NextHalfedge;
                    }

                    HashSet<int> adjacentFaces = new HashSet<int>();
                    for (int i = 0; i < path.Count; i++)
                    {
                        adjacentFaces.Add(dualGraphWithHMDP.DualPlanktonMesh.Halfedges[path[i]].AdjacentFace);
                    }

                    if (adjacentFaces.Contains(-1))
                    {
                        needToConnect.Add(false);
                    }
                    else
                    {
                        needToConnect.Add(true);
                        pairVolumeJunctionsIndex.Add(pair);
                    }
                }
                #endregion
                DA.SetDataList("NeedToConnect", needToConnect);




                #region 可视化部分
                BoundaryPolylinePoints.Clear();
                BoundaryCornerTextDots.Clear();
                BoundarySegmentTextDots.Clear();

                #region 设置表示场地边界的Polyline
                BoundaryPolylinePoints = closedBoundaryPoints;
                #endregion

                #region 设置场地角点的TextDot
                List<string> args = new List<string>();
                //for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                //{
                //    List<int> dFaces = dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i).ToList();
                //    args.Add(string.Join<int>(";", faceIndexsAroundVertex))
                //}
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        args.Add(string.Join<int>(";", faceIndexsAroundVertex[verticesIndexForEachBoundarySegment[i][j]]));
                    }
                    TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], args[i]), sortedBoundarySegments[i].To);
                    BoundaryCornerTextDots.Add(boundaryCornerTextDot);
                }
                //for (int i = 0; i < boundaryCorner.Count; i++)
                //{
                //    TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], args[sortedBoundaryCornerIndexs[i]]), boundaryCorner[i]);
                //    BoundaryCornerTextDots.Add(boundaryCornerTextDot);
                //}
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

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            //base.DrawViewportWires(args);

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
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            //base.DrawViewportMeshes(args);

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
            get { return new Guid("50b86078-83e3-46a8-85b5-3bee8c71c1ea"); }
        }
    }
}