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
            // base.Message = "开发中 Todo:对偶图中，边界顶点布置方式";

            BoundaryPolylinePoints = new List<Point3d>();
            BoundaryCornerTextDots = new List<TextDot>();
            BoundarySegmentTextDots = new List<TextDot>();

            JunctionTextDot = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> BoundaryPolylinePoints;

        private List<TextDot> BoundaryCornerTextDots;
        private List<TextDot> BoundarySegmentTextDots;

        private List<TextDot> JunctionTextDot;

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
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "经过顶点分裂的对偶图", GH_ParamAccess.item);

            pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);

            pManager.AddGenericParameter("SubBoundarySegments", "SS", "分割过后形成的子BoundarySegment对象", GH_ParamAccess.list);
            // pManager.AddIntegerParameter("HalfedgeIndexOfSubBoundarySegment", "hI", "分割过后形成的子BoundarySegment对象所对应的Halfedge的序号", GH_ParamAccess.list);

            // pManager.AddPointParameter("BoundaryVerticesPoints", "BVP", "", GH_ParamAccess.list);

            // pManager.AddPointParameter("BoundaryCorner", "BC", "在边界上的角点", GH_ParamAccess.list);
            // pManager.AddIntegerParameter("ReorganizedCornerIndex", "RCIndex", "经过重新组织的边界polyline的角点对应的对偶图中角点的index", GH_ParamAccess.list);
            //pManager.AddPointParameter("DivisionPoints", "DP", "要生成分割线的边界点", GH_ParamAccess.tree);
            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);
            pManager.AddIntegerParameter("BSIndexContainVolumeJunctions", "BSI", "包含边界分裂点的BoundarySegment序号", GH_ParamAccess.list);
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
                // DA.SetDataList("BoundaryCorner", boundaryCorner);

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

                //#region 为BoundarySegment的IncludedDVertice属性赋值
                //for (int i = 0; i < sortedBoundarySegments.Count; i++)
                //{
                //    sortedBoundarySegments[i].IncludedDVertice = new List<int>();
                //    sortedBoundarySegments[i].IncludedDVertice.AddRange(verticesIndexForEachBoundarySegment[i]);
                //}
                //#endregion

                #endregion

                List<List<int>> faceIndexsAroundVertex = new List<List<int>>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(-1);
                }


                #region 区分出每个Segment上的volumeJunction
                List<List<int>> volumeJunctionForEachBoundarySegment = new List<List<int>>();

                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    volumeJunctionForEachBoundarySegment.Add(new List<int>());

                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        if (dualGraphWithHMDP.DVertice_Inner[verticesIndexForEachBoundarySegment[i][j]].Count > 1
                            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBoundarySegment[i][j]).Contains(-1))
                        {
                            if (dualGraphWithHMDP.DVertice_Volume[verticesIndexForEachBoundarySegment[i][j]].Count != 1)
                            {
                                volumeJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                            }
                        }
                    }
                }
                #endregion

                #region 找到包含边界分裂点的BoundarySegment序号
                List<int> bSIndexContainVolumeJunctions = new List<int>();
                for (int i = 0; i < volumeJunctionForEachBoundarySegment.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionForEachBoundarySegment[i].Count; j++)
                    {
                        bSIndexContainVolumeJunctions.Add(i);
                    }
                }
                #endregion
                DA.SetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions);

                #region 每个segment上的junction列表
                // 每个segment上，除去首尾两点，剩下的都是junction
                List<List<int>> juntionsIndexForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    // 深拷贝
                    juntionsIndexForEachBoundarySegment.Add(new List<int>());
                    juntionsIndexForEachBoundarySegment[i].AddRange(verticesIndexForEachBoundarySegment[i]);
                    // 去除首尾顶点
                    juntionsIndexForEachBoundarySegment[i].RemoveAt(0);
                    juntionsIndexForEachBoundarySegment[i].RemoveAt(juntionsIndexForEachBoundarySegment[i].Count - 1);
                }
                #endregion

                #region 对于对偶图边界上所有的Junction点进行分裂，并生成新的vertices和junctions列表
                DualGraphWithHM newDualGraphWithHM = new DualGraphWithHM(dualGraphWithHMDP);
                List<List<int>> newJunctionForAllSegment = new List<List<int>>();
                List<List<int>> newVerticesForAllSegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    List<int> newJunctionForCurrentSegment;
                    newDualGraphWithHM = DecomposeJunction(newDualGraphWithHM, 
                                                           juntionsIndexForEachBoundarySegment[i], 
                                                           out newJunctionForCurrentSegment);

                    // List<string> debug1 = UtilityFunctions.PrintFacesVertices(newDualGraphWithHM.DualPlanktonMesh);

                    List<int> newVerticesForCurrentSegment = new List<int>();
                    newVerticesForCurrentSegment.AddRange(newJunctionForCurrentSegment);
                    newVerticesForCurrentSegment.Insert(0, verticesIndexForEachBoundarySegment[i].First());
                    newVerticesForCurrentSegment.Add(verticesIndexForEachBoundarySegment[i].Last());

                    newJunctionForAllSegment.Add(new List<int>());
                    newVerticesForAllSegment.Add(new List<int>());
                    newJunctionForAllSegment[i].AddRange(newJunctionForCurrentSegment);
                    newVerticesForAllSegment[i].AddRange(newVerticesForCurrentSegment);
                }
                #endregion

                // List<string> debug = UtilityFunctions.PrintFacesVertices(newDualGraphWithHM.DualPlanktonMesh);
                DA.SetData("DualGraphWithHM", newDualGraphWithHM);

                #region 用junction点分割boundarysegment
                List<List<FaceEdgeSegment>> sortedSubBoundarySegment = new List<List<FaceEdgeSegment>>();
                List<List<List<int>>> junctionCorrespondingWhichInnerForAllSegment = new List<List<List<int>>>();
                List<List<Point3d>> junctionPointForAllSegment = new List<List<Point3d>>();
                // List<List<int>> hIndexForAllSegment = new List<List<int>>();

                // 此时边界上的点，全部都是juntion点（除去每条边的首尾）
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    // 添加分支
                    sortedSubBoundarySegment.Add(new List<FaceEdgeSegment>());
                    junctionCorrespondingWhichInnerForAllSegment.Add(new List<List<int>>());
                    junctionPointForAllSegment.Add(new List<Point3d>());
                    // hIndexForAllSegment.Add(new List<int>());

                    // 初始化out参数
                    List<List<int>> junctionCorrespondingWhichInnerForCurrentSegment;
                    List<Point3d> junctionPointForCurrentSegment;
                    // List<int> hIndexForCurrentSegment;

                    List<FaceEdgeSegment> subBoundarySegment = SplitBoundarySegment(sortedBoundarySegments[i],
                                                                                    // verticesIndexForEachBoundarySegment[i],
                                                                                    newJunctionForAllSegment[i],
                                                                                    newVerticesForAllSegment[i],
                                                                                    newDualGraphWithHM,
                                                                                    out junctionCorrespondingWhichInnerForCurrentSegment,
                                                                                    out junctionPointForCurrentSegment);
                                                                                    //out hIndexForCurrentSegment);

                    // 将结果添加到对应分支上
                    sortedSubBoundarySegment[i].AddRange(subBoundarySegment);
                    junctionCorrespondingWhichInnerForCurrentSegment.ForEach(item =>
                    {
                        junctionCorrespondingWhichInnerForAllSegment[i].Add(new List<int>());
                        junctionCorrespondingWhichInnerForAllSegment[i].Last().AddRange(item);
                    });
                    junctionPointForAllSegment[i].AddRange(junctionPointForCurrentSegment);
                    // hIndexForAllSegment[i].AddRange(hIndexForCurrentSegment);
                }
                #endregion
                DA.SetDataList("SubBoundarySegments", sortedSubBoundarySegment);


                #region 计算DivisionPoints
                //Dictionary<int, Point3d> index_JunctionPoint = new Dictionary<int, Point3d>();
                //for (int i = 0; i < sortedSubBoundarySegment.Count; i++)
                //{
                //    for (int j = 0; j < sortedSubBoundarySegment[i].Count; j++)
                //    {
                //        index_JunctionPoint.TryAdd(sortedSubBoundarySegment[i][j].IncludedDVertice[0], sortedSubBoundarySegment[i][j].To);
                //        index_JunctionPoint.TryAdd(sortedSubBoundarySegment[i][j].IncludedDVertice[1], sortedSubBoundarySegment[i][j].From);
                //    }
                //}

                //List<List<Point3d>> volumeJunctionPointLoL = new List<List<Point3d>>();
                //for (int i = 0; i < volumeJunctionForEachBoundarySegment.Count; i++)
                //{
                //    volumeJunctionPointLoL.Add(new List<Point3d>());
                //    if (volumeJunctionForEachBoundarySegment[i].Count != 0)
                //    {
                //        for (int j = 0; j < volumeJunctionForEachBoundarySegment[i].Count; j++)
                //        {
                //            volumeJunctionPointLoL[i].Add(index_JunctionPoint[volumeJunctionForEachBoundarySegment[i][j]]);
                //        }
                //    }
                //    else
                //    {
                //        Point3d point = Point3d.Unset;
                //        volumeJunctionPointLoL[i].Add(point);
                //    }
                //}

                //DataTree<Point3d> volumeJunctionPointDT = UtilityFunctions.LoLToDataTree<Point3d>(volumeJunctionPointLoL);
                //DA.SetDataTree(3, volumeJunctionPointDT);
                #endregion

                #region 可视化部分
                BoundaryPolylinePoints.Clear();
                BoundaryCornerTextDots.Clear();
                BoundarySegmentTextDots.Clear();

                JunctionTextDot.Clear();

                #region 设置表示场地边界的Polyline
                BoundaryPolylinePoints = closedBoundaryPoints;
                #endregion

                #region 设置场地角点的TextDot
                List<string> args = new List<string>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    List<int> dFaces = dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i).ToList();
                    // 去除-1
                    dFaces.RemoveAt(0);
                    args.Add(string.Join<int>(";", dFaces));
                }
                for (int i = 0; i < boundaryCorner.Count; i++)
                {
                    TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], args[sortedBoundaryCornerIndexs[i]]), boundaryCorner[i]);
                    BoundaryCornerTextDots.Add(boundaryCornerTextDot);
                }
                #endregion

                #region 设置Junction
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    for (int j = 0; j < newJunctionForAllSegment[i].Count; j++)
                    {
                        string str = string.Join<int>(";", junctionCorrespondingWhichInnerForAllSegment[i][j]);
                        TextDot junctionTextDot = new TextDot(string.Format("{0}|{1}", newJunctionForAllSegment[i][j], str), junctionPointForAllSegment[i][j]);
                        JunctionTextDot.Add(junctionTextDot);
                    }
                }
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

        private DualGraphWithHM DecomposeJunction(//INode<BoundarySegment> INodeToSplit,
                                      DualGraphWithHM dualGraphWithHM,
                                      List<int> junctionForCurrentSegment,
                                      out List<int> newJunctionForCurrentSegment)
                                      // out List<int> newVerticesForCurrentSegment)
                                      // SortedDictionary<int, List<int>> volumeJunction_itsSplit)
        {
            PlanktonMesh newPlanktonMesh = new PlanktonMesh(dualGraphWithHM.DualPlanktonMesh);

            // 存储分裂后，该segment上新的junction列表
            newJunctionForCurrentSegment = new List<int>();

            // 判断是否需要分裂
            for (int i = 0; i < junctionForCurrentSegment.Count; i++)
            {
                List<int> adjacentFace = newPlanktonMesh.Vertices.GetVertexFaces(junctionForCurrentSegment[i]).ToList();
                // adjacentFace.Remove(-1);
                adjacentFace.Reverse();

                // 与adjacentFace对应的该顶点的出半边的列表
                List<int> correspondingH = new List<int>();
                List<int> hs = newPlanktonMesh.Vertices.GetHalfedges(junctionForCurrentSegment[i]).ToList();
                for (int j = 0; j < adjacentFace.Count; j++)
                {
                    for (int k = 0; k < hs.Count; k++)
                    {
                        if (newPlanktonMesh.Halfedges[hs[k]].AdjacentFace == adjacentFace[j])
                        {
                            correspondingH.Add(hs[k]);
                            break;
                        }
                    }
                }

                List<int> correspondingHEnd = new List<int>();
                for (int j = 0; j < adjacentFace.Count; j++)
                {
                    for (int k = 0; k < hs.Count; k++)
                    {
                        if (newPlanktonMesh.Halfedges[hs[k]].AdjacentFace == adjacentFace[j])
                        {
                            correspondingHEnd.Add(newPlanktonMesh.Halfedges.EndVertex(hs[k]));
                            break;
                        }
                    }
                }

                if (correspondingHEnd.Count > 3)
                {
                    int iter = 0;
                    int newVertexIndex = -1;
                    int newAddedHalfedgeEndStartFromOrigin = correspondingHEnd[0];
                    // int count = 0;
                    // List<int> newVertice = new List<int>();
                    while (iter < correspondingHEnd.Count - 3)
                    {
                        int[,] currentHPairIndexs = new int[2, 2] {
                                /* 初始化索引号为0的行 */
                                { junctionForCurrentSegment[i],
                                  newAddedHalfedgeEndStartFromOrigin},
                                /* 初始化索引号为1的行 */
                                { junctionForCurrentSegment[i],
                                  correspondingHEnd[1 + iter]}};

                        newPlanktonMesh = SplitVertex(newPlanktonMesh,
                                                      currentHPairIndexs,
                                                      junctionForCurrentSegment[i],
                                                      out newVertexIndex,
                                                      out newAddedHalfedgeEndStartFromOrigin);

                        List<string> debugPrint = UtilityFunctions.PrintFacesVertices(newPlanktonMesh);
                        // newVertice.Add(newVertexIndex);

                        newJunctionForCurrentSegment.Add(newVertexIndex);

                        iter++;
                    }
                    // volumeJunction_itsSplit.Add(volumeJunctionForCurrentSegment[i], newVertice);

                    newJunctionForCurrentSegment.Add(junctionForCurrentSegment[i]);
                }
                else
                {
                    newPlanktonMesh = new PlanktonMesh(newPlanktonMesh);

                    newJunctionForCurrentSegment.Add(junctionForCurrentSegment[i]);
                    continue;
                }
            }

            // 构造新的DualPlanktonMesh对象
            DualGraphWithHM newDualGraphWithHM = new DualGraphWithHM(newPlanktonMesh,
                                                                     dualGraphWithHM.GraphNodes,
                                                                     dualGraphWithHM.GraphTables,
                                                                     dualGraphWithHM.Volume_VolumeNode,
                                                                     dualGraphWithHM.Volume_Inner,
                                                                     dualGraphWithHM.Inner_DFace,
                                                                     dualGraphWithHM.DFaceIndexsAroundOuterNodes);
            newDualGraphWithHM.UndividedGraphNodes = dualGraphWithHM.UndividedGraphNodes;
            newDualGraphWithHM.UndividedGraphTable = dualGraphWithHM.UndividedGraphTable;
            return newDualGraphWithHM;
        }

        /// <summary>
        /// 分裂生成新的顶点，并且将边的连接关系转移到新的顶点上
        /// </summary>
        /// <param name="P">要修改的PlanktonMesh</param>
        /// <param name="HStartAndEndIndexs">一对相邻的半边的首尾顶点序号</param>
        /// <param name="currentVolumeJunction">当前要分裂的VolumeJunction点</param>
        /// <param name="newVertexIndex">新生成的顶点的序号</param>
        /// <param name="newAddedHEndStartFromOrigin">新生成的半边的序号（由原来的volumeJunction点发出的）</param>
        /// <returns></returns>
        private PlanktonMesh SplitVertex(PlanktonMesh P,
                                        int[,] HStartAndEndIndexs,
                                        int currentVolumeJunction,
                                        out int newVertexIndex,
                                        out int newAddedHEndStartFromOrigin)
        {
            PlanktonMesh pDeepCopy = new PlanktonMesh(P);

            #region 分裂，生成新的点
            int halfedgeStartVertexIndex = HStartAndEndIndexs[0, 0];
            int halfedgeEndVertexIndex = HStartAndEndIndexs[0, 1];

            // 计算要分裂的边两侧的Face，准备好两个Face的顶点列表，新增midVertex作为分裂产生的新顶点

            // 找到halfedgeStartVertexIndex和halfedgeEndVertexIndex所对应的那条半边
            int halfedgeIndex = -1;
            foreach (int hfIndex in pDeepCopy.Vertices.GetHalfedges(halfedgeStartVertexIndex))
            {
                if (pDeepCopy.Halfedges.EndVertex(hfIndex) == halfedgeEndVertexIndex)
                {
                    halfedgeIndex = hfIndex;
                    break;
                }
            }
            if (halfedgeIndex == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ResetHalfedgeStart函数出错，找不到输入起点和终点所对应的半边。");
                newVertexIndex = -1;
                newAddedHEndStartFromOrigin = -1;
                return null;
            }

            int adjacentFaceIndex = pDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = pDeepCopy.Halfedges[pDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            List<int> viAroundAdjacentFace = new List<int>();
            // 包围AdjacentFace的顶点
            viAroundAdjacentFace = pDeepCopy.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();
            List<int> viAroundPairAdjacentFace = new List<int>();
            if (pairAdjacentFaceIndex == -1)
            {
                viAroundPairAdjacentFace = new List<int>();
            }
            else
            {
                // 包围PairAdjcentFace的顶点
                viAroundPairAdjacentFace = pDeepCopy.Faces.GetFaceVertices(pairAdjacentFaceIndex).ToList<int>();
            }

            int startVertexIndex = halfedgeStartVertexIndex;
            int endVertexIndex = halfedgeEndVertexIndex;

            // 向P中增加顶点
            Point3d startVertex = pDeepCopy.Vertices[startVertexIndex].ToPoint3d();
            Point3d endVertex = pDeepCopy.Vertices[endVertexIndex].ToPoint3d();
            Point3d midVertex = (startVertex + endVertex) / 2;
            pDeepCopy.Vertices.Add(midVertex);
            int midVertexIndex = pDeepCopy.Vertices.Count - 1;
            newVertexIndex = midVertexIndex;

            #region 构造添加顶点后的新Face
            for (int i = 0; i < viAroundAdjacentFace.Count; i++)
            {
                if (viAroundAdjacentFace[i] == startVertexIndex)
                {
                    viAroundAdjacentFace.Insert(i + 1, midVertexIndex);
                    break;
                }
            }
            for (int i = 0; i < viAroundPairAdjacentFace.Count; i++)
            {
                if (viAroundPairAdjacentFace[i] == endVertexIndex)
                {
                    viAroundPairAdjacentFace.Insert(i + 1, midVertexIndex);
                    break;
                }
            }
            #endregion
            #endregion

            #region 将原来的表示连接关系的一对半边移动到一个新分裂生成的顶点上 
            int anotherHalfedgeStartVertexIndex = HStartAndEndIndexs[1, 0];
            int anotherHalfedgeEndVertexIndex = HStartAndEndIndexs[1, 1];

            // 找到anotherHalfedgeStartVertexIndex和anotherHalfedgeEndVertexIndex所对应的那条半边
            int anotherHalfedgeIndex = -1;
            foreach (int hfIndex in pDeepCopy.Vertices.GetHalfedges(anotherHalfedgeStartVertexIndex))
            {
                if (pDeepCopy.Halfedges.EndVertex(hfIndex) == anotherHalfedgeEndVertexIndex)
                {
                    anotherHalfedgeIndex = hfIndex;
                    break;
                }
            }
            if (anotherHalfedgeIndex == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ResetHalfedgeStart函数出错，找不到输入起点和终点所对应的半边。");
                newAddedHEndStartFromOrigin = -1;
                return null;
            }

            // int pairAnotherHalfedgeIndex = pDeepCopy.Halfedges.GetPairHalfedge(anotherHalfedgeIndex);

            int anotherAdjacentFaceIndex = pDeepCopy.Halfedges[anotherHalfedgeIndex].AdjacentFace;
            int pairAnotherAdjacentFaceIndex = pDeepCopy.Halfedges[pDeepCopy.Halfedges.GetPairHalfedge(anotherHalfedgeIndex)].AdjacentFace;

            #region 包围anotherAdjacentFace的顶点
            List<int> viAroundAnotherAdjacentFace = new List<int>();
            int currentHalfedgeIndex = anotherHalfedgeIndex;
            do
            {
                viAroundAnotherAdjacentFace.Add(pDeepCopy.Halfedges[currentHalfedgeIndex].StartVertex);
                currentHalfedgeIndex = pDeepCopy.Halfedges[currentHalfedgeIndex].NextHalfedge;
            } while (currentHalfedgeIndex != anotherHalfedgeIndex);
            #endregion

            #region 包围anotherPairAdjacentFace的顶点，此时anotherPairAdjacentFace就是上面的AdjacentFace
            List<int> viAroundAnotherPairAdjacentFace = viAroundAdjacentFace;
            #endregion

            #region 修改包围AdjacentFace的顶点列表以及包围PairAdjacentFace的顶点列表
            // 用目标顶点属于哪个面，来判断该怎么修改viAroundAdjacentFace和viAroundPairAdjacentFace列表
            // 获取newVertex在Face顶点列表的序号
            int indexInAnotherFace = viAroundAnotherAdjacentFace.IndexOf(newVertexIndex);
            // 获取newVertex在PairFace顶点列表的序号
            int indexInAnotherPairFace = viAroundAnotherPairAdjacentFace.IndexOf(newVertexIndex);

            // 当目标顶点newVertex同时不属于这两个Face时报错
            if (indexInAnotherFace == -1 && indexInAnotherPairFace == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的targetVertex " + currentVolumeJunction.ToString() + " 不属于要改变的面");
                newAddedHEndStartFromOrigin = -1;
                return null;
            }
            // 当目标顶点newVertex属于PairFace时
            else if (indexInAnotherFace == -1 && indexInAnotherPairFace != -1)
            {
                int currentVolumeJunctionIndex = viAroundAnotherAdjacentFace.IndexOf(currentVolumeJunction);
                viAroundAnotherAdjacentFace.Insert(currentVolumeJunctionIndex + 1, newVertexIndex);

                // 从PairFace顶点列表中删除第一个（即vertexToSplit）
                viAroundAnotherPairAdjacentFace.Remove(currentVolumeJunction);

            }
            // 当目标顶点newVertex属于AdjacentFace时
            else
            {
                int currentVolumeJunctionIndex = viAroundAnotherPairAdjacentFace.IndexOf(currentVolumeJunction);
                viAroundAnotherPairAdjacentFace.Insert(currentVolumeJunctionIndex + 1, newVertexIndex);

                // 从Face顶点列表中删除第一个（即vertexToSplit）
                viAroundAnotherAdjacentFace.Remove(currentVolumeJunction);
            }
            #endregion
            #endregion

            #region 转移，更新，生成新的PlanktonMesh
            // 建立需要进行置换的面的index列表和VertexIndex列表，方便转移P的Face属性
            List<int> needChangeFaceIndexs = new List<int>();
            // 先anotherPair，再another
            needChangeFaceIndexs.Add(pairAnotherAdjacentFaceIndex);
            needChangeFaceIndexs.Add(anotherAdjacentFaceIndex);
            List<List<int>> needChangeViAround = new List<List<int>>();
            // 先anotherPair，再another
            needChangeViAround.Add(viAroundAnotherPairAdjacentFace);
            needChangeViAround.Add(viAroundAnotherAdjacentFace);

            // 转移P的Vertex属性
            List<PlanktonXYZ> pPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < pDeepCopy.Vertices.Count; i++)
            {
                pPlanktonVertex.Add(pDeepCopy.Vertices[i].ToXYZ());
            }

            // 转移P的Face属性，同时替换掉两个应该删除的面
            List<List<int>> pFaceVertexOrder = new List<List<int>>();

            for (int i = 0; i < pDeepCopy.Faces.Count; i++)
            {
                if (needChangeFaceIndexs.Contains(i))
                {
                    pFaceVertexOrder.Add(new List<int>());
                    pFaceVertexOrder[i].AddRange(needChangeViAround[needChangeFaceIndexs.IndexOf(i)]);
                }
                else
                {
                    pFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = pDeepCopy.Faces.GetFaceVertices(i);
                    pFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }

            // 用转移的P的Vertex属性和P的Face属性来构造新的PlanktonMesh newP
            PlanktonMesh newP = new PlanktonMesh();

            newP.Vertices.AddVertices(pPlanktonVertex);
            newP.Faces.AddFaces(pFaceVertexOrder);
            #endregion

            newP.Compact();

            #region 输出新生成的半边的序号（由原来的volumeJunction点发出的）
            newAddedHEndStartFromOrigin = -1;
            for (int i = 0; i < newP.Vertices.GetHalfedges(currentVolumeJunction).Length; i++)
            {
                if (newP.Halfedges.EndVertex(newP.Vertices.GetHalfedges(currentVolumeJunction)[i]) == newVertexIndex)
                {
                    newAddedHEndStartFromOrigin = newP.Halfedges.EndVertex(newP.Vertices.GetHalfedges(currentVolumeJunction)[i]);
                }
            }

            return newP;
            #endregion
        }

        private List<FaceEdgeSegment> SplitBoundarySegment(BoundarySegment currentBoundarySegment,
                                          // List<int> currentVerticesIndexForBoundarySegment,
                                          List<int> junctionsOnCurrentBoundary,
                                          List<int> verticeOnCurrentBoundary,
                                          DualGraphWithHM dualGraphWithHMDP,
                                          out List<List<int>> junctionCorrespondingWhichInnerForCurrentSegment,
                                          out List<Point3d> junctionPointForCurrentSegment)
                                          // out List<int> hIndexForCurrentSegment)
        {
            List<FaceEdgeSegment> subSegments = new List<FaceEdgeSegment>();

            // 没有junction在该segment上时
            if (junctionsOnCurrentBoundary.Count == 0)
            {
                junctionCorrespondingWhichInnerForCurrentSegment = new List<List<int>>();
                junctionPointForCurrentSegment = new List<Point3d>();

                int hIndex = -1;
                // hIndexForCurrentSegment = new List<int>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Halfedges.Count; i++)
                {
                    if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[i].StartVertex == verticeOnCurrentBoundary.Last() 
                        && dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(i) == verticeOnCurrentBoundary.First())
                    {
                        //hIndexForCurrentSegment.Add(i);
                        hIndex = i;
                    }
                }

                Line segment = currentBoundarySegment.Lines[0];
                string label = currentBoundarySegment.Label;
                List<int> includedDVertice = verticeOnCurrentBoundary;
                // includedDVertice中的顺序（即includedDVertex的顺序），应该跟halfedge一样，是逆时针的
                // 而Line中From，To的顺序，跟boundarySegment一样，是顺时针的
                includedDVertice.Reverse();
                
                FaceEdgeSegment subSegment = new FaceEdgeSegment(segment, label, includedDVertice, hIndex);

                subSegments.Add(subSegment);
                return subSegments;
            }

            #region 计算分割后每段子segment上所包含的顶点
            List<int> junctionIndexsInList = new List<int>();
            for (int i = 0; i < junctionsOnCurrentBoundary.Count; i++)
            {
                junctionIndexsInList.Add(verticeOnCurrentBoundary.IndexOf(junctionsOnCurrentBoundary[i]));
            }
            // 构造[0,junctionIndex,...,verticeOnCurrentBoundary.Count - 1]的列表，用来在下面截取verticeOnCurrentBoundary列表
            junctionIndexsInList.Insert(0, 0);
            junctionIndexsInList.Add(verticeOnCurrentBoundary.Count - 1);
            List<int[]> pairIndex = new List<int[]>();
            for (int i = 0; i < junctionIndexsInList.Count - 1; i++)
            {
                int[] pair = new int[] { junctionIndexsInList[i], junctionIndexsInList[i + 1] };
                pairIndex.Add(pair);
            }
            // 截取verticeOnCurrentBoundary列表
            List<List<int>> dVerticesForEachSubSegment = new List<List<int>>();
            for (int i = 0; i < pairIndex.Count; i++)
            {
                List<int> subList = verticeOnCurrentBoundary.Skip(pairIndex[i][0]).Take(pairIndex[i][1] - pairIndex[i][0] + 1).ToList();
                // dVerticesForEachSubSegment中的顺序（即includedDVertex的顺序），应该跟halfedge一样，是逆时针的
                // 而Line中From，To的顺序，跟boundarySegment一样，是顺时针的
                subList.Reverse();
                dVerticesForEachSubSegment.Add(subList);
            }
            #endregion

            #region 计算这段BoundarySegment上所有的Junction点所对应的t，及Point点
            // 计算分割点的t
            junctionCorrespondingWhichInnerForCurrentSegment = new List<List<int>>();

            List<double> innerAreaProportionForCurrentSegment = new List<double>();
            List<int> hIndexForCurrentSegment = new List<int>();
            for (int i = 0; i < junctionsOnCurrentBoundary.Count; i++)
            {
                // List<int> value = dualGraphWithHMDP.DVertice_Inner[junctionsOnCurrentBoundary[i]];

                junctionCorrespondingWhichInnerForCurrentSegment.Add(new List<int>());

                int dface1 = -1;
                int dface2 = -1;

                int h1 = -1;
                int h2 = -1;

                for (int j = 0; j < dualGraphWithHMDP.DualPlanktonMesh.Halfedges.Count; j++)
                {
                    if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[j].StartVertex == junctionsOnCurrentBoundary[i]
                        && dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(j) == verticeOnCurrentBoundary[verticeOnCurrentBoundary.IndexOf(junctionsOnCurrentBoundary[i]) - 1])
                    {
                        dface1 = dualGraphWithHMDP.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                        h1 = j;
                        continue;
                    }
                    if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[j].StartVertex == verticeOnCurrentBoundary[verticeOnCurrentBoundary.IndexOf(junctionsOnCurrentBoundary[i]) + 1]
                        && dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(j) == junctionsOnCurrentBoundary[i])
                    {
                        dface2 = dualGraphWithHMDP.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                        h2 = j;
                        continue;
                    }
                }

                int inner1 = -1;
                int inner2 = -1;
                foreach (var item in dualGraphWithHMDP.Inner_DFace)
                {
                    if (dface1 == item.Value)
                    {
                        inner1 = item.Key;
                        continue;
                    }
                    if (dface2 == item.Value)
                    {
                        inner2 = item.Key;
                        continue;
                    }
                }

                junctionCorrespondingWhichInnerForCurrentSegment[i].Add(inner1);
                junctionCorrespondingWhichInnerForCurrentSegment[i].Add(inner2);

                double inner1AreaProportion = dualGraphWithHMDP.GraphNodes[inner1].NodeAttribute.NodeAreaProportion;
                double inner2AreaProportion = dualGraphWithHMDP.GraphNodes[inner2].NodeAttribute.NodeAreaProportion;

                if (i == 0)
                {
                    innerAreaProportionForCurrentSegment.Add(inner1AreaProportion);
                    innerAreaProportionForCurrentSegment.Add(inner2AreaProportion);
                    hIndexForCurrentSegment.Add(h1);
                    hIndexForCurrentSegment.Add(h2);
                }
                else
                {
                    innerAreaProportionForCurrentSegment.Add(inner2AreaProportion);
                    hIndexForCurrentSegment.Add(h2);
                }
            }



            // 计算t及Point点
            List<double> tList = new List<double>();
            junctionPointForCurrentSegment = new List<Point3d>();
            double sum = innerAreaProportionForCurrentSegment.Sum();
            // 注意t值是需要累加的
            double previousT = 0;
            for (int i = 0; i < junctionsOnCurrentBoundary.Count; i++)
            {
                double currentProportion = innerAreaProportionForCurrentSegment[i] / sum;
                double t = previousT + currentProportion;
                previousT = t;
                Point3d tPoint = currentBoundarySegment.Lines[0].PointAt(t);
                tList.Add(t);
                junctionPointForCurrentSegment.Add(tPoint);
                
            }
            #endregion

            #region 分割生成子Sement
            Curve segmentLineCurve = currentBoundarySegment.Lines[0].ToNurbsCurve();

            Curve[] subSegmentLineCurve = segmentLineCurve.Split(tList);

            List<Line> subLines = new List<Line>();
            foreach (var item in subSegmentLineCurve)
            {
                Point3d start = item.PointAtStart;
                Point3d end = item.PointAtEnd;
                subLines.Add(new Line(start, end));
            }

            for (int i = 0; i < subLines.Count; i++)
            {
                
                string label = string.Format("{0}{1}", currentBoundarySegment.Label, i);
                // BoundarySegment subSegment = new BoundarySegment(subLines[i], label);


                FaceEdgeSegment subSegment = new FaceEdgeSegment(subLines[i], label, dVerticesForEachSubSegment[i], hIndexForCurrentSegment[i]);

                //subSegment.IncludedDVertice = new List<int>();
                //subSegment.IncludedDVertice.AddRange(dVerticesForEachSubSegment[i]);
                subSegments.Add(subSegment);
            }
            #endregion

            return subSegments;
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

            for (int i = 0; i < JunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(JunctionTextDot[i], Color.Red, Color.White, Color.White);
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

            for (int i = 0; i < JunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(JunctionTextDot[i], Color.Red, Color.White, Color.White);
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