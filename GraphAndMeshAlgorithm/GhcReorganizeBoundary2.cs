using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Attributes;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;
using Grasshopper.GUI.Canvas;

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
        /// 场地划分时所需要的tX的数量
        /// </summary>
        private int TXCount;

        /// <summary>
        /// 场地划分时所需要的tY的数量
        /// </summary>
        private int TYCount;

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
            // pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);
            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);
            // pManager.AddIntegerParameter("BSIndexContainVolumeJunctions", "BSI", "包含边界分裂点的BoundarySegment序号", GH_ParamAccess.list);

            pManager.AddIntegerParameter("VerticesIndexForEachBS", "VIFBS", "每个BS上点对应的对偶图中的index", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("indexOnEachBS", "IOBS", "每个BS上的VolumeJunctionIndex", GH_ParamAccess.tree);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);


            // pManager.AddBooleanParameter("NeedToConnect", "NTC", "的这一对分裂点是否作为分界线", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectVolumeJunctionsIndex", "N_VJI", "需要构成分界线的两个volumeJunction的Index", GH_ParamAccess.list);
            pManager.AddGenericParameter("NeedToConnectBSIndex", "N_BSI", "这两个volumeJunction所对应的bsIndex", GH_ParamAccess.list);

            pManager.AddIntegerParameter("tXCount", "tXCount", "tX的数量", GH_ParamAccess.item);
            pManager.AddIntegerParameter("tYCount", "tYCount", "tY的数量", GH_ParamAccess.item);
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
                #region 对无序的boundarySegment进行排序
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
                #endregion
                DA.SetDataList("SortedBoundarySegments", sortedBoundarySegments);

                #region 按照从W开始，逆时针的顺序，输出场地边界的角点
                List<Point3d> boundaryCorner = new List<Point3d>();

                // 利用叉积来判断第一个Segment的方向时候需要反转
                if (IsFirstSegmentOpppsite(sortedBoundarySegments[0].Lines[0].Direction,
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
                //DA.SetData("ReorganizedBoundary", reorganizedBoundary);

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
                DataTree<int> verticesIndexForEachBSDT = UtilityFunctions.LoLToDataTree<int>(verticesIndexForEachBoundarySegment);
                DA.SetDataTree(1, verticesIndexForEachBSDT);

                #endregion

                #region 找到每个segment上的volumeJunction
                // 找到每个vertex所对应的faceIndex
                List<List<int>> faceIndexsAroundVertex = new List<List<int>>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(-1);
                }

                DataTree<int> indexOnEachBS = new DataTree<int>();
                List<List<int>> volumeJunctionForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    indexOnEachBS.EnsurePath(i);
                    volumeJunctionForEachBoundarySegment.Add(new List<int>());
                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        if (faceIndexsAroundVertex[verticesIndexForEachBoundarySegment[i][j]].Count > 1
                            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBoundarySegment[i][j]).Contains(-1))
                        {
                            indexOnEachBS.Branch(i).Add(verticesIndexForEachBoundarySegment[i][j]);
                            volumeJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                        }
                    }
                }
                #endregion

                #region DataTree<int> indexOnEachBS
                DA.SetDataTree(2, indexOnEachBS);
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
                //DA.SetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions);

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

                #region 得到相邻的两个VolumeJunction是否构成分界线的判断，以及需要构成分界线的两个volumeJunction的Index，这两个volumeJunction所对应的bsIndex
                List<int[]> needToConnectVolumeJunctionsIndex = new List<int[]>();
                List<int[]> needToConnectBSIndex = new List<int[]>();

                List<int[]> indexPair = new List<int[]>();
                List<int[]> bsIndexPair = new List<int[]>();
                for (int i = 0; i < volumeJunctionIndex.Count; i++)
                {
                    int[] pair = new int[2] { volumeJunctionIndex[i], volumeJunctionIndex[(i + 1) % volumeJunctionIndex.Count] };
                    indexPair.Add(pair);

                    int[] bsPair = new int[2] { bSIndexContainVolumeJunctions[i], bSIndexContainVolumeJunctions[(i + 1) % bSIndexContainVolumeJunctions.Count] };
                    bsIndexPair.Add(bsPair);
                }

                List<bool> needToConnect = new List<bool>();
                for (int i = 0; i < indexPair.Count; i++)
                {
                    int faceIndex = -1;
                    int[] hs = null;
                    for (int j = 0; j < dualGraphWithHMDP.DualPlanktonMesh.Faces.Count; j++)
                    {
                        int[] faceVertex = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetFaceVertices(j);
                        int[] intersect = faceVertex.Intersect(indexPair[i]).ToArray();
                        if (intersect.Length == indexPair[i].Length)
                        {
                            faceIndex = j;
                            hs = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetHalfedges(j);
                        }
                    }

                    if (faceIndex == -1)
                    {
                        needToConnect.Add(false);
                        continue;
                    }

                    int hIndex = -1;
                    for (int j = 0; j < hs.Length; j++)
                    {
                        if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[hs[j]].StartVertex == indexPair[i][0])
                        {
                            hIndex = hs[j];
                        }
                    }

                    int currH = hIndex;
                    List<int> path = new List<int>();
                    while (dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(currH) != indexPair[i][1])
                    {
                        path.Add(currH);
                        currH = dualGraphWithHMDP.DualPlanktonMesh.Halfedges[currH].NextHalfedge;
                    }

                    HashSet<int> adjacentFaces = new HashSet<int>();
                    for (int j = 0; j < path.Count; j++)
                    {
                        adjacentFaces.Add(dualGraphWithHMDP.DualPlanktonMesh.Halfedges[path[j]].AdjacentFace);
                    }

                    if (adjacentFaces.Contains(-1))
                    {
                        needToConnect.Add(false);
                    }
                    else
                    {
                        needToConnect.Add(true);
                        needToConnectVolumeJunctionsIndex.Add(indexPair[i]);
                        needToConnectBSIndex.Add(bsIndexPair[i]);
                    }
                }

                #endregion
                // DA.SetDataList("NeedToConnect", needToConnect);
                DA.SetDataList("NeedToConnectVolumeJunctionsIndex", needToConnectVolumeJunctionsIndex);
                DA.SetDataList("NeedToConnectBSIndex", needToConnectBSIndex);

                #region 计算tXDT和tYDT的数量
                int txCount = 0;
                int tyCount = 0;
                for (int i = 0; i < needToConnectVolumeJunctionsIndex.Count; i++)
                {
                    int interval = CalInterval(needToConnectBSIndex[i], 4);

                    int[] tCount = CalTCount(interval,
                                                 sortedBoundarySegments[needToConnectBSIndex[i][0]],
                                                 sortedBoundarySegments[needToConnectBSIndex[i][1]]);

                    //tXDT.EnsurePath(i);
                    //tXDT.Branch(i).AddRange(tXList.Skip(txCount).Take(tCount[0]));
                    //tYDT.EnsurePath(i);
                    //tYDT.Branch(i).AddRange(tYList.Skip(tyCount).Take(tCount[1]));

                    txCount += tCount[0];
                    tyCount += tCount[1];
                }
                #endregion
                TXCount = txCount;
                TYCount = tyCount;
                DA.SetData("tXCount", txCount);
                DA.SetData("tYCount", tyCount);


                #region 可视化部分
                BoundaryPolylinePoints.Clear();
                BoundaryCornerTextDots.Clear();
                BoundarySegmentTextDots.Clear();

                #region 设置表示场地边界的Polyline
                BoundaryPolylinePoints = closedBoundaryPoints;
                #endregion

                #region 设置场地角点的TextDot
                //List<string> args = new List<string>();
                //for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                //{
                //    List<int> dFaces = dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i).ToList();
                //    args.Add(string.Join<int>(";", faceIndexsAroundVertex))
                //}
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        string arg = string.Join<int>(";", faceIndexsAroundVertex[verticesIndexForEachBoundarySegment[i][0]]);
                        TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], arg), sortedBoundarySegments[i].From);
                        BoundaryCornerTextDots.Add(boundaryCornerTextDot);
                    }
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
                    Vector3d verticalVector = Vector3d.CrossProduct(boundarySegments[i].Lines[0].Direction, Vector3d.ZAxis);
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
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="bsIndexPair"></param>
        /// <param name="bsCount"></param>
        /// <returns></returns>
        private int CalInterval(int[] bsIndexPair, int bsCount)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            int index0 = bsIndexPair[0];
            int index1 = bsIndexPair[1];

            if (index0 < index1)
            {
                counterClockwiseCount = index1 - index0;
                clockwiseCount = bsCount - index1 + index0;
            }
            else
            {
                counterClockwiseCount = bsCount - index0 + index1;
                clockwiseCount = index0 - index1;
            }

            if (counterClockwiseCount < clockwiseCount)
            {
                return counterClockwiseCount;
            }
            else
            {
                return clockwiseCount;
            }
        }

        private int[] CalTCount(int interval,
                                BoundarySegment bsForPoint1,
                                BoundarySegment bsForPoint2)
        {
            #region 基础计算
            #region 求point1和point2所对应的BS的法线（逆时针90度）
            // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
            Vector3d vectorForBS1 = new Vector3d(bsForPoint1.To - bsForPoint1.From);
            Vector3d vectorForBS2 = new Vector3d(bsForPoint2.To - bsForPoint2.From);
            Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.From, bsForPoint1.To);
            Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.From, bsForPoint2.To);

            // 求在x轴正方向或y轴正方向的投影
            Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
            Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
            // 投影后向量的normal，为正X正Y方向
            Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
            Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
            #endregion
            #endregion

            int[] tCount;
            if (interval == 0)
            {
                // interval == 0
                if (nVectorForBS1.X == 0)
                {
                    // 法向量为y轴
                    tCount = new int[2] { 0, 1 };
                }
                else
                {
                    // 法向量为x轴
                    tCount = new int[2] { 1, 0 };
                }
            }
            else if (interval == 1)
            {
                // interval == 1
                tCount = new int[2] { 0, 0 };
            }
            else
            {
                // interval == 2
                if (nVectorForBS1.X == 0)
                {
                    // 法向量为y轴
                    tCount = new int[2] { 0, 1 };
                }
                else
                {
                    // 法向量为x轴
                    tCount = new int[2] { 1, 0 };
                }
            }

            return tCount;
        }

        private int WhichQuadrant(Point3d point1, Point3d point2)
        {
            double dx = point2.X - point1.X;
            double dy = point2.Y - point1.Y;

            if (dx == 0)
            {
                // y轴
                if (dy > 0)
                {
                    // y轴正方向，并入第一象限的情况计算
                    return 1;
                }
                else
                {
                    // y轴反方向，并入第三象限的情况计算
                    return 3;
                }
            }
            if (dy == 0)
            {
                // x轴
                if (dx > 0)
                {
                    // x轴正方向，并入第四象限计算
                    return 4;
                }
                else
                {
                    // x轴反方向，并入第二象限计算
                    return 2;
                }
            }

            if (dy / dx > 0 && dx > 0)
            {
                return 1;
            }
            else if (dy / dx > 0 && dx < 0)
            {
                return 3;
            }
            else if (dy / dx < 0 && dx < 0)
            {
                return 2;
            }
            else
            {
                return 4;
            }
        }


        private Vector3d CalProjectVector(Point3d point1, Point3d point2)
        {
            Vector3d projectVector;

            // 要考虑原来的向量是属于哪个象限的问题
            int quadrant = WhichQuadrant(point1, point2);

            if (point2.X - point1.X != 0)
            {
                double k = (point2.Y - point1.Y) / (point2.X - point1.X);
                if (k > 0 && k <= 1)
                {
                    // 斜率小于45度时
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k > 1)
                {
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                }
                else if (k == 0)
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k >= -1 && k < 0)
                {
                    // 斜率小于45度时
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                }
            }
            else
            {
                // 即y轴正方向
                if (quadrant == 1)
                {
                    projectVector = new Vector3d(0, 1, 0);
                }
                else//quadrant == 3
                {
                    projectVector = new Vector3d(0, -1, 0);
                }
            }

            return projectVector;
        }

        /// <summary>
        /// 计算向量逆时针旋转九十度的单位法向量
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        private Vector3d Normal(Vector3d vector)
        {
            double L = Math.Sqrt(vector * vector);
            return new Vector3d(-vector.Y / L, vector.X / L, vector.Z);
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


        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new CompLabelAttribute(this);
        }

        public class CompLabelAttribute : GH_ComponentAttributes
        {
            public CompLabelAttribute(GhcReorganizeBoundary2 component) : base(component) { }

            protected override void Layout()
            {
                base.Layout();
                /* 先执行base.Layout()，可以按GH电池默认方式计算电池的出/入口需要的高度，我们在下面基于这个高度进行更改 */
                Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 42.0f);
            }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);/* 执行基本的电池渲染 */

                /* 额外的电池渲染，仅在“Objects”这个渲染轨道绘制 */
                if (channel == GH_CanvasChannel.Objects)
                {
                    RectangleF buttonRect1 = /* 按钮的位置 */ new RectangleF(Bounds.X, Bounds.Bottom - 42, Bounds.Width, 20.0f);

                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect1.Inflate(-2.0f, -2.0f);

                    using(GH_Capsule capsule1 = GH_Capsule.CreateCapsule(buttonRect1, GH_Palette.Normal))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule1.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("tXCount:{0}", ((GhcReorganizeBoundary2)Owner).TXCount),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.Black,
                                        buttonRect1,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });

                    RectangleF buttonRect2 = new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);
                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect2.Inflate(-2.0f, -2.0f);

                    using (GH_Capsule capsule2 = GH_Capsule.CreateCapsule(buttonRect2, GH_Palette.Normal))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule2.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("tYCount:{0}", ((GhcReorganizeBoundary2)Owner).TYCount),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.Black,
                                        buttonRect2,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });
                }
            }
        }
    }
}