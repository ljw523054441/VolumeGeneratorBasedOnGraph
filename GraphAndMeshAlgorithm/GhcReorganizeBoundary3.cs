using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
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
    public class GhcReorganizeBoundary3 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcReorganizeBoundary3 class.
        /// </summary>
        public GhcReorganizeBoundary3()
          : base("ReorganizeBoundary3", "ReorganizeBoundary3",
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
        /// 场地划分时所需要的BoundaryT的数量
        /// </summary>
        private string AllLayerBoundaryTCountsString;

        /// <summary>
        /// 场地划分时所需要的InnerT的数量
        /// </summary>
        private string AllLayerInnerTMaxCountString;

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

            pManager.AddGenericParameter("VerticesIndexForEachBS", "VIFBS", "每个BS上点对应的对偶图中的index", GH_ParamAccess.list);

            pManager.AddGenericParameter("IndexOnEachBS", "IOBS", "每个BS上的VolumeJunctionIndex", GH_ParamAccess.list);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);


            // pManager.AddBooleanParameter("NeedToConnect", "NTC", "的这一对分裂点是否作为分界线", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectVolumeJunctionsIndex", "N_VJI", "需要构成分界线的两个volumeJunction的Index", GH_ParamAccess.list);
            pManager.AddGenericParameter("NeedToConnectBSIndex", "N_BSI", "这两个volumeJunction所对应的bsIndex", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectPairCorrespondingFaceIndex", "N_CFI", "需要构成分界线的两个volumeJunction所对对应的FaceIndex", GH_ParamAccess.list);

            pManager.AddIntegerParameter("BoundaryTCounts", "boundaryTCounts", "在BS上的点所需要的t值的数量", GH_ParamAccess.list);
            pManager.AddIntegerParameter("InnerTMaxCounts", "InnerTMaxCounts", "由BS上的点生成的转折所需要的t值的最大数量", GH_ParamAccess.list);
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

                List<DataTree<int>> allLayerVerticesIndexForEachBS = new List<DataTree<int>>();
                List<DataTree<int>> allLayerVolumeJunctionsIndexOnEachBS = new List<DataTree<int>>();
                List<List<int[]>> allLayerVertexIndexPairs = new List<List<int[]>>();
                List<List<int[]>> allLayerBsIndexPairs = new List<List<int[]>>();

                List<List<int>> allLayerPairCorrespondingFaceIndex = new List<List<int>>();

                // List<int> allLayerTXCount = new List<int>();
                // List<int> allLayerTYCount = new List<int>();

                List<int> allLayerBoundaryTCounts = new List<int>();
                List<int> allLayerInnerTMaxCounts = new List<int>();

                // 记录目前有哪几个Face已经有了vertexPairIndex，即处理过，correspondingFaceIndexs与allFaceIndexs比较，可以判断是否需要再计算下一层
                List<int> correspondingFaceIndexs = new List<int>();
                List<int> allFaceIndexs = new List<int>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Faces.Count; i++)
                {
                    allFaceIndexs.Add(i);
                }

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
                List<BoundarySegment> sortedBSFromW = new List<BoundarySegment>();
                for (int i = 0; i < outerNodeLabels.Count; i++)
                {
                    for (int j = 0; j < boundarySegments.Count; j++)
                    {
                        if (boundarySegments[j].Label == outerNodeLabels[i])
                        {
                            sortedBSFromW.Add(boundarySegments[j]);
                        }
                    }
                }
                #endregion
                #endregion
                //DA.SetDataList("SortedBoundarySegments", sortedBSFromW);

                #region 按照从W开始，逆时针的顺序，输出场地边界的角点
                List<Point3d> boundaryCorner = new List<Point3d>();

                // 利用叉积来判断第一个Segment的方向时候需要反转
                if (IsFirstSegmentOpppsite(sortedBSFromW[0].Lines[0].Direction,
                                           new Vector3d((sortedBSFromW[1].From + sortedBSFromW[1].To) / 2 - sortedBSFromW[0].From)))
                {
                    sortedBSFromW[0].Reverse();
                }

                boundaryCorner.Add(sortedBSFromW[0].From);
                // 从第二个开始，调整from和to
                for (int i = 1; i < sortedBSFromW.Count; i++)
                {
                    if (sortedBSFromW[i].From != sortedBSFromW[i - 1].To)
                    {
                        sortedBSFromW[i].Reverse();
                    }
                    boundaryCorner.Add(sortedBSFromW[i].From);
                }
                // 判断生成的polyline是否闭合
                if (sortedBSFromW.Last().To != sortedBSFromW[0].From)
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

                #region 还没有开始对bs进行shift时

                #region 输出从W开始的逆时针方向的，场地边界每个角点所对应的对偶图顶点的序号
                List<List<int>> faceIndexsAroundOuterNodes = dualGraphWithHMDP.DFaceIndexsAroundOuterNodes;
                // 去除-1的情况，第0个永远是-1
                for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
                {
                    faceIndexsAroundOuterNodes[i].RemoveAt(0);
                }
                List<int> sortedBoundaryCornerIndexs = new List<int>();
                for (int i = 0; i < sortedBSFromW.Count; i++)
                {
                    sortedBoundaryCornerIndexs.Add(faceIndexsAroundOuterNodes[i].First());
                }
                #endregion

                #region 输出从W开始的逆时针方向的，场地边界上所应该布置的对偶图顶点序号(最外圈的verticesIndexForEachBSDT)
                List<List<int>> verticesIndexForEachBSLoL = new List<List<int>>();
                for (int i = 0; i < sortedBSFromW.Count; i++)
                {
                    verticesIndexForEachBSLoL.Add(new List<int>());
                    verticesIndexForEachBSLoL[i].AddRange(faceIndexsAroundOuterNodes[i]);
                }
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

                //DataTree<int> volumeJunctionsIndexOnEachBS = new DataTree<int>();
                List<List<int>> volumeJunctionForEachBSLoL = new List<List<int>>();
                for (int i = 0; i < sortedBSFromW.Count; i++)
                {
                    //volumeJunctionsIndexOnEachBS.EnsurePath(i);
                    volumeJunctionForEachBSLoL.Add(new List<int>());
                    for (int j = 0; j < verticesIndexForEachBSLoL[i].Count; j++)
                    {
                        if (faceIndexsAroundVertex[verticesIndexForEachBSLoL[i][j]].Count > 1
                            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBSLoL[i][j]).Contains(-1))
                        {
                            //volumeJunctionsIndexOnEachBS.Branch(i).Add(verticesIndexForEachBoundarySegment[i][j]);
                            volumeJunctionForEachBSLoL[i].Add(verticesIndexForEachBSLoL[i][j]);
                        }
                    }
                }
                #endregion

                #region 输出表示边界分裂点的Text
                List<string> volumeJunctionText = new List<string>();
                for (int i = 0; i < volumeJunctionForEachBSLoL.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionForEachBSLoL[i].Count; j++)
                    {
                        string arg = string.Format("{0} | {1}", volumeJunctionForEachBSLoL[i][j], string.Join<int>(";", faceIndexsAroundVertex[verticesIndexForEachBSLoL[i][j]]));
                        volumeJunctionText.Add(arg);
                    }
                }
                #endregion
                DA.SetDataList("VolumeJunctionsTexts", volumeJunctionText);

                #endregion

                int currentLayer = 0;

                #region shift的部分
                #region 找到包含边界分裂点的BoundarySegment序号
                List<int> bSIndexContainVolumeJunctions = new List<int>();
                List<int> volumeJunctionIndex = new List<int>();
                for (int i = 0; i < volumeJunctionForEachBSLoL.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionForEachBSLoL[i].Count; j++)
                    {
                        bSIndexContainVolumeJunctions.Add(i);
                        volumeJunctionIndex.Add(volumeJunctionForEachBSLoL[i][j]);
                    }
                }
                #endregion

                #region 求shift的值
                int firstVolumeJunctionCorrespondingBSIndex = bSIndexContainVolumeJunctions[0];
                int shift = firstVolumeJunctionCorrespondingBSIndex;
                #endregion

                #region 对BS进行shift
                // 把每个volumeJuntion所对应的bs的Index都向前减去bsShift个
                for (int i = 0; i < bSIndexContainVolumeJunctions.Count; i++)
                {
                    bSIndexContainVolumeJunctions[i] = bSIndexContainVolumeJunctions[i] - shift;
                }
                // 对BS的列表进行Shift
                List<BoundarySegment> sortedBS = Shift<BoundarySegment>(sortedBSFromW, shift);
                #endregion
                DA.SetDataList("SortedBoundarySegments", sortedBS);

                #region 对verticesIndexForEachBS进行shift
                List<List<int>> shiftedVerticesIndexForEachBSLoL = Shift<List<int>>(verticesIndexForEachBSLoL, shift);
                DataTree<int> verticesIndexForEachBSDT = UtilityFunctions.LoLToDataTree<int>(shiftedVerticesIndexForEachBSLoL);
                allLayerVerticesIndexForEachBS.Add(verticesIndexForEachBSDT);
                #endregion

                #region 对volumeJunctionsIndexOnEachBS进行shift
                List<List<int>> shiftedVolumeJunctionsIndexOnEachBSLoL = Shift<List<int>>(volumeJunctionForEachBSLoL, shift);
                DataTree<int> volumeJunctionsIndexOnEachBS = UtilityFunctions.LoLToDataTree<int>(shiftedVolumeJunctionsIndexOnEachBSLoL);
                allLayerVolumeJunctionsIndexOnEachBS.Add(volumeJunctionsIndexOnEachBS);
                #endregion

                #endregion

                #region 计算最外层的BoundaryTCount
                int currentLayerBoundaryTCount = volumeJunctionsIndexOnEachBS.DataCount;
                allLayerBoundaryTCounts.Add(currentLayerBoundaryTCount);
                #endregion

                #region 得到相邻的两个VolumeJunction是否构成分界线的判断，以及需要构成分界线的两个volumeJunction的Index，这两个volumeJunction所对应的bsIndex
                List<int[]> needToConnectVolumeJunctionsIndex = new List<int[]>();
                List<int[]> needToConnectBSIndex = new List<int[]>();

                List<int[]> indexPairs = new List<int[]>();
                List<int[]> bsIndexPairs = new List<int[]>();
                for (int i = 0; i < volumeJunctionIndex.Count; i++)
                {
                    int[] pair = new int[2] { volumeJunctionIndex[i], volumeJunctionIndex[(i + 1) % volumeJunctionIndex.Count] };
                    indexPairs.Add(pair);

                    int[] bsPair = new int[2] { bSIndexContainVolumeJunctions[i], bSIndexContainVolumeJunctions[(i + 1) % bSIndexContainVolumeJunctions.Count] };
                    bsIndexPairs.Add(bsPair);
                }

                List<int> pairsCorrespondingFaceIndexs = new List<int>();
                for (int i = 0; i < indexPairs.Count; i++)
                {
                    int faceIndex = -1;
                    int[] hs = null;
                    for (int j = 0; j < dualGraphWithHMDP.DualPlanktonMesh.Faces.Count; j++)
                    {
                        int[] faceVertex = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetFaceVertices(j);
                        int[] intersect = faceVertex.Intersect(indexPairs[i]).ToArray();
                        if (intersect.Length == indexPairs[i].Length)
                        {
                            faceIndex = j;
                            hs = dualGraphWithHMDP.DualPlanktonMesh.Faces.GetHalfedges(j);
                        }
                    }

                    if (faceIndex == -1)
                    {
                        continue;
                    }

                    pairsCorrespondingFaceIndexs.Add(faceIndex);

                    int hIndex = -1;
                    for (int j = 0; j < hs.Length; j++)
                    {
                        if (dualGraphWithHMDP.DualPlanktonMesh.Halfedges[hs[j]].StartVertex == indexPairs[i][0])
                        {
                            hIndex = hs[j];
                        }
                    }

                    int currH = hIndex;
                    List<int> path = new List<int>();
                    while (dualGraphWithHMDP.DualPlanktonMesh.Halfedges.EndVertex(currH) != indexPairs[i][1])
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
                        continue;
                    }
                    else
                    {
                        needToConnectVolumeJunctionsIndex.Add(indexPairs[i]);
                        needToConnectBSIndex.Add(bsIndexPairs[i]);
                    }
                }

                #endregion
                allLayerVertexIndexPairs.Add(needToConnectVolumeJunctionsIndex);
                allLayerBsIndexPairs.Add(needToConnectBSIndex);
                allLayerPairCorrespondingFaceIndex.Add(pairsCorrespondingFaceIndexs);

                #region 原来的tXDT和tYDT合并为innerTForEachBS
                // 首层是能够确定有几个InnerT的，所以变量名不叫TMax
                int currentLayerInnerTCount = 0;
                for (int i = 0; i < needToConnectVolumeJunctionsIndex.Count; i++)
                {
                    int interval = CalFirstBSInterval(needToConnectBSIndex[i], sortedBS.Count);

                    //// 不论interval为0,1,2哪个值，tCount要么是0，要么是1，取最大值tMaxCount
                    //int tMaxCount = 1;
                    //currentLayerInnerTMaxCountForEachBS.Add(tMaxCount);
                    int tCount = CalTCount(interval,
                                           sortedBS[allLayerBsIndexPairs[currentLayer][i][0]],
                                           sortedBS[allLayerBsIndexPairs[currentLayer][i][1]],
                                           0,
                                           0);
                    //currentLayerInnerTCountForEachBS.Add(tCount);
                    currentLayerInnerTCount += tCount;
                }
                allLayerInnerTMaxCounts.Add(currentLayerInnerTCount);
                #endregion

                bool needToCalNextLayer = false;
                correspondingFaceIndexs.AddRange(pairsCorrespondingFaceIndexs);
                List<int> except = allFaceIndexs.Except(correspondingFaceIndexs).ToList();
                // 判断是否要计算inner
                if (except.Count == 0)
                {
                    // 不需要计算下一层
                    needToCalNextLayer = false;
                }
                else
                {
                    // 需要计算下一层
                    needToCalNextLayer = true;
                }

                while (needToCalNextLayer)
                {
                    DataTree<int> innerVerticesIndexForEachBS;
                    DataTree<int> innerVolumeJunctionsIndexOnEachBS;
                    List<int[]> innerNeedToConnectVolumeJunctionsIndex;
                    List<int[]> innerNeedToConnectBSIndex;
                    List<int> innerIndexPairCorrespondingFaceIndexs;

                    GetInnerPart(dualGraphWithHMDP.DualPlanktonMesh,
                                 pairsCorrespondingFaceIndexs,
                                 needToConnectVolumeJunctionsIndex,
                                 out innerVerticesIndexForEachBS,
                                 out innerVolumeJunctionsIndexOnEachBS,
                                 out innerNeedToConnectVolumeJunctionsIndex,
                                 out innerNeedToConnectBSIndex,
                                 out innerIndexPairCorrespondingFaceIndexs);

                    int currentLayerInnerBoundaryTCount = innerVolumeJunctionsIndexOnEachBS.DataCount;
                    allLayerBoundaryTCounts.Add(currentLayerInnerBoundaryTCount);

                    allLayerVerticesIndexForEachBS.Add(innerVerticesIndexForEachBS);
                    allLayerVolumeJunctionsIndexOnEachBS.Add(innerVolumeJunctionsIndexOnEachBS);
                    allLayerVertexIndexPairs.Add(innerNeedToConnectVolumeJunctionsIndex);
                    allLayerBsIndexPairs.Add(innerNeedToConnectBSIndex);
                    allLayerPairCorrespondingFaceIndex.Add(innerIndexPairCorrespondingFaceIndexs);

                    #region 计算这一层的tXDT和tYDT的数量
                    //int innerTXCount = 0;
                    //int innerTYCount = 0;
                    //for (int i = 0; i < innerNeedToConnectVolumeJunctionsIndex.Count; i++)
                    //{
                    //    int interval = CalInterval(innerNeedToConnectBSIndex[i], indexPairsCorrespondingFaceIndexs.Count);

                    //    int[] tCount = CalTCount(interval,)
                    //}
                    #endregion
                    #region 原来的tXDT和tYDT合并为innerTForEachBS
                    // 其他内部层不能够确定有几个InnerT的，所以变量名叫TMax
                    //List<int> currentLayerInnerTMaxCountForEachBS = new List<int>();
                    int currentLayerInnerTMaxCount = 0;
                    for (int i = 0; i < innerNeedToConnectVolumeJunctionsIndex.Count; i++)
                    {
                        // 不论interval为0,1,2哪个值，tCount要么是0，要么是1，取最大值tMaxCount
                        int tMaxCount = 1;
                        //currentLayerInnerTMaxCountForEachBS.Add(tMaxCount);
                        currentLayerInnerTMaxCount += tMaxCount;
                    }
                    allLayerInnerTMaxCounts.Add(currentLayerInnerTMaxCount);
                    #endregion

                    correspondingFaceIndexs.AddRange(innerIndexPairCorrespondingFaceIndexs);
                    except = allFaceIndexs.Except(correspondingFaceIndexs).ToList();
                    // 判断是否要计算下一层
                    if (except.Count == 0)
                    {
                        // 不需要计算下一层
                        needToCalNextLayer = false;
                    }
                    else
                    {
                        // 需要计算下一层
                        needToCalNextLayer = true;
                    }
                }

                //// List<DataTree<int>> 转换为 List<List<List<int>>>
                //List<List<List<int>>> allLayerVerticesIndexForEachBSLoL = new List<List<List<int>>>();
                //for (int i = 0; i < allLayerVerticesIndexForEachBS.Count; i++)
                //{
                //    List<List<int>> currentVerticesIndexForEachBSLoL = UtilityFunctions.DataTreeToLoL<int>(allLayerVerticesIndexForEachBS[i]);
                //    allLayerVerticesIndexForEachBSLoL.Add(currentVerticesIndexForEachBSLoL);
                //}
                //List<List<List<int>>> allLayerVolumeJunctionsIndexOnEachBSLoL = new List<List<List<int>>>();
                //for (int i = 0; i < allLayerVolumeJunctionsIndexOnEachBS.Count; i++)
                //{
                //    List<List<int>> currentVolumeJunctionsIndexOnEachBS = UtilityFunctions.DataTreeToLoL<int>(allLayerVolumeJunctionsIndexOnEachBS[i]);
                //    allLayerVolumeJunctionsIndexOnEachBSLoL.Add(currentVolumeJunctionsIndexOnEachBS);
                //}

                //List<int> abc = allLayerVolumeJunctionsIndexOnEachBS[0].Branch(1);


                //List<GH_Structure<GH_Integer>> allLayerVerticesIndexForEachBS_GHS = new List<GH_Structure<GH_Integer>>();
                //for (int i = 0; i < allLayerVerticesIndexForEachBS.Count; i++)
                //{
                //    GH_Structure<GH_Integer> ghstructure = UtilityFunctions.DataTreeToGH_Structure_Int(allLayerVerticesIndexForEachBS[i]);
                //    allLayerVerticesIndexForEachBS_GHS.Add(ghstructure);
                //}

                //DataTree<DataTree<int>> dt = new DataTree<DataTree<int>>();
                //dt.Add(allLayerVerticesIndexForEachBS[0]);


                //DA.SetData("VerticesIndexForEachBS", dt);
                DA.SetDataList("VerticesIndexForEachBS", allLayerVerticesIndexForEachBS);
                DA.SetDataList("IndexOnEachBS", allLayerVolumeJunctionsIndexOnEachBS);
                DA.SetDataList("NeedToConnectVolumeJunctionsIndex", allLayerVertexIndexPairs);
                DA.SetDataList("NeedToConnectBSIndex", allLayerBsIndexPairs);
                DA.SetDataList("NeedToConnectPairCorrespondingFaceIndex", allLayerPairCorrespondingFaceIndex);
                DA.SetDataList("BoundaryTCounts", allLayerBoundaryTCounts);
                DA.SetDataList("InnerTMaxCounts", allLayerInnerTMaxCounts);

                string string1 = string.Join(";", allLayerBoundaryTCounts);
                string string2 = string.Join(";", allLayerInnerTMaxCounts);
                AllLayerBoundaryTCountsString = string1;
                AllLayerInnerTMaxCountString = string2;


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
                for (int i = 0; i < sortedBSFromW.Count; i++)
                {
                    for (int j = 0; j < verticesIndexForEachBSLoL[i].Count; j++)
                    {
                        string arg = string.Join<int>(";", faceIndexsAroundVertex[verticesIndexForEachBSLoL[i][0]]);
                        TextDot boundaryCornerTextDot = new TextDot(string.Format("{0} | {1}", sortedBoundaryCornerIndexs[i], arg), sortedBSFromW[i].From);
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

        private List<T> Shift<T>(List<T> originList, int shift)
        {
            // 深拷贝
            List<T> newList = new List<T>();
            newList.AddRange(originList);

            int iter = 0;
            while (iter < shift)
            {
                T item = newList[0];
                newList.RemoveAt(0);
                newList.Add(item);
                iter++;
            }

            return newList;
        }

        private void GetInnerPart(PlanktonMesh D,
                                  List<int> indexPairCorrespondingFaceIndexs,
                                  List<int[]> indexPairs,
                                  //out PlanktonMesh subD,
                                  out DataTree<int> innerVerticesIndexForEachBS,
                                  out DataTree<int> innerVolumeJunctionsIndexOnEachBS,
                                  out List<int[]> innerNeedToConnectVolumeJunctionsIndex,
                                  out List<int[]> innerNeedToConnectBSIndex,
                                  out List<int> innerIndexPairCorrespondingFaceIndexs)
        {
            PlanktonMesh DDeepCopy = new PlanktonMesh(D);

            #region 找到每个vertex所对应的FaceIndex
            List<List<int>> faceIndexsAroundVertex = new List<List<int>>();
            for (int i = 0; i < DDeepCopy.Vertices.Count; i++)
            {
                faceIndexsAroundVertex.Add(new List<int>());
                faceIndexsAroundVertex[i].AddRange(DDeepCopy.Vertices.GetVertexFaces(i));
                // 去除-1
                faceIndexsAroundVertex[i].Remove(-1);
            }
            #endregion

            #region 求verticesIndexForEachBS
            innerVerticesIndexForEachBS = new DataTree<int>();
            for (int i = 0; i < indexPairCorrespondingFaceIndexs.Count; i++)
            {
                innerVerticesIndexForEachBS.EnsurePath(i);

                List<int> faceVertices = DDeepCopy.Faces.GetFaceVertices(indexPairCorrespondingFaceIndexs[i]).ToList();

                // int vertex0InFaceVertices = faceVertices.IndexOf(indexPairs[i][0]);
                int vertex1InFaceVertices = faceVertices.IndexOf(indexPairs[i][1]);
                // 找到从由vertex1到vertex0的切片
                List<int> slice = new List<int>();
                int iter1 = vertex1InFaceVertices;
                while (faceVertices[(iter1 + faceVertices.Count) % faceVertices.Count] != indexPairs[i][0])
                {
                    iter1++;
                    slice.Add(faceVertices[(iter1 + faceVertices.Count) % faceVertices.Count]);
                }
                // 移除indexPairs[i][0]
                slice.RemoveAt(slice.Count - 1);

                innerVerticesIndexForEachBS.Branch(i).AddRange(slice);
            }
            #endregion

            #region 求volumeJunctionsIndexOnEachBS
            innerVolumeJunctionsIndexOnEachBS = new DataTree<int>();
            for (int i = 0; i < indexPairCorrespondingFaceIndexs.Count; i++)
            {
                innerVolumeJunctionsIndexOnEachBS.EnsurePath(i);

                List<int> volumeJunctions = new List<int>();
                for (int j = 0; j < innerVerticesIndexForEachBS.Branch(i).Count; j++)
                {
                    if (faceIndexsAroundVertex[innerVerticesIndexForEachBS.Branch(i)[j]].Count > 1
                        && DDeepCopy.Vertices.GetVertexFaces(innerVerticesIndexForEachBS.Branch(i)[j]).Contains(indexPairCorrespondingFaceIndexs[i]))
                    {
                        volumeJunctions.Add(innerVerticesIndexForEachBS.Branch(i)[j]);
                    }
                }

                innerVolumeJunctionsIndexOnEachBS.Branch(i).AddRange(volumeJunctions);
            }
            #endregion

            #region 找到包含边界分裂点的BoundarySegment序号
            List<int> bsIndexContainVolumeJunctions = new List<int>();
            List<int> volumeJunctionIndex = new List<int>();
            for (int i = 0; i < innerVolumeJunctionsIndexOnEachBS.BranchCount; i++)
            {
                for (int j = 0; j < innerVolumeJunctionsIndexOnEachBS.Branch(i).Count; j++)
                {
                    bsIndexContainVolumeJunctions.Add(i);
                    volumeJunctionIndex.Add(innerVolumeJunctionsIndexOnEachBS.Branch(i)[j]);
                }
            }
            #endregion

            #region 得到相邻的两个volumeJunction是否构成分界线的判断，以及需要构成分界线的两个volumeJunction的Index，这两个volumeJunction所对应的bsIndex
            innerNeedToConnectVolumeJunctionsIndex = new List<int[]>();
            innerNeedToConnectBSIndex = new List<int[]>();

            List<int[]> innerIndexPairs = new List<int[]>();
            List<int[]> innerBSIndexPairs = new List<int[]>();
            for (int i = 0; i < volumeJunctionIndex.Count; i++)
            {
                int[] pair = new int[2] { volumeJunctionIndex[i], volumeJunctionIndex[(i + 1) % volumeJunctionIndex.Count] };
                innerIndexPairs.Add(pair);

                int[] bsPair = new int[2] { bsIndexContainVolumeJunctions[i], bsIndexContainVolumeJunctions[(i + 1) % bsIndexContainVolumeJunctions.Count] };
                innerBSIndexPairs.Add(bsPair);
            }

            innerIndexPairCorrespondingFaceIndexs = new List<int>();
            for (int i = 0; i < innerIndexPairs.Count; i++)
            {
                int faceIndex = -1;
                int[] hs = null;
                for (int j = 0; j < DDeepCopy.Faces.Count; j++)
                {
                    int[] faceVertex = DDeepCopy.Faces.GetFaceVertices(j);
                    int[] intersect = faceVertex.Intersect(innerIndexPairs[i]).ToArray();
                    if (intersect.Length == innerIndexPairs[i].Length)
                    {
                        faceIndex = j;
                        hs = DDeepCopy.Faces.GetHalfedges(j);
                    }
                }

                if (faceIndex == -1)
                {
                    continue;
                }

                innerIndexPairCorrespondingFaceIndexs.Add(faceIndex);

                int hIndex = -1;
                for (int j = 0; j < hs.Length; j++)
                {
                    if (DDeepCopy.Halfedges[hs[j]].StartVertex == innerIndexPairs[i][0])
                    {
                        hIndex = hs[j];
                    }
                }

                int currH = hIndex;
                List<int> path = new List<int>();
                while (DDeepCopy.Halfedges.EndVertex(currH) != innerIndexPairs[i][1])
                {
                    path.Add(currH);
                    currH = DDeepCopy.Halfedges[currH].NextHalfedge;
                }

                HashSet<int> adjacentFaces = new HashSet<int>();
                for (int j = 0; j < path.Count; j++)
                {
                    adjacentFaces.Add(DDeepCopy.Halfedges[path[j]].AdjacentFace);
                }

                if (adjacentFaces.Contains(-1))
                {
                    continue;
                }
                else
                {
                    innerNeedToConnectVolumeJunctionsIndex.Add(innerIndexPairs[i]);
                    innerNeedToConnectBSIndex.Add(innerBSIndexPairs[i]);
                }
            }
            #endregion
        }

        /// <summary>
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="bsIndexPair"></param>
        /// <param name="bsCount"></param>
        /// <returns></returns>
        private int CalFirstBSInterval(int[] bsIndexPair, int bsCount)
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

        /// <summary>
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="pointOnWhichLineIndexPairs"></param>
        /// <param name="turningIndex"></param>
        /// <returns></returns>
        private int CalInnerBSInterval(int[] pointOnWhichLineIndexPairs, List<int> turningIndex)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            //int index0;
            //int index1;

            List<int[]> turningRange = new List<int[]>();
            for (int i = 0; i < turningIndex.Count - 1; i++)
            {
                int[] pair = new int[2] { turningIndex[i], turningIndex[(i + 1) % turningIndex.Count] };
                turningRange.Add(pair);
            }

            // 需要注意等于的位置
            List<int> indexs = new List<int>();
            for (int i = 0; i < pointOnWhichLineIndexPairs.Length; i++)
            {
                if (pointOnWhichLineIndexPairs[i] <= turningIndex[0])
                {
                    indexs.Add(turningRange.Count);
                }
                else if (pointOnWhichLineIndexPairs[i] > turningIndex.Last())
                {
                    indexs.Add(turningRange.Count);
                }
                else
                {
                    for (int j = 0; j < turningRange.Count; j++)
                    {
                        if (pointOnWhichLineIndexPairs[i] > turningRange[j][0] && pointOnWhichLineIndexPairs[i] <= turningRange[j][1])
                        {
                            indexs.Add(j);
                        }
                    }
                }
            }

            if (indexs[0] < indexs[1])
            {
                counterClockwiseCount = indexs[1] - indexs[0];
                clockwiseCount = (turningRange.Count + 1) - indexs[1] + indexs[0];
            }
            else
            {
                counterClockwiseCount = (turningRange.Count + 1) - indexs[0] + indexs[1];
                clockwiseCount = indexs[0] - indexs[1];
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

        private int CalTCount(int interval,
                              BoundarySegment bsForPoint1,
                              BoundarySegment bsForPoint2,
                              int point1OnWhichSegment,
                              int point2OnWhichSegment)
        {
            #region 基础计算
            #region 求point1和point2所对应的BS的法线（逆时针90度）
            // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
            Vector3d vectorForBS1 = new Vector3d(bsForPoint1.Lines[point1OnWhichSegment].To - bsForPoint1.Lines[point1OnWhichSegment].From);
            Vector3d vectorForBS2 = new Vector3d(bsForPoint2.Lines[point2OnWhichSegment].To - bsForPoint2.Lines[point2OnWhichSegment].From);
            Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.Lines[point1OnWhichSegment].From, bsForPoint1.Lines[point1OnWhichSegment].To);
            Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.Lines[point2OnWhichSegment].From, bsForPoint2.Lines[point2OnWhichSegment].To);

            // 求在x轴正方向或y轴正方向的投影
            Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
            Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
            // 投影后向量的normal，为正X正Y方向
            Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
            Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
            #endregion
            #endregion

            int tCount;
            if (interval == 0)
            {
                tCount = 1;
            }
            else if (interval == 1)
            {
                tCount = 0;
            }
            else
            {
                tCount = 1;
            }
            //int[] tCount;
            //if (interval == 0)
            //{
            //    // interval == 0
            //    if (nVectorForBS1.X == 0)
            //    {
            //        // 法向量为y轴
            //        tCount = new int[2] { 0, 1 };
            //    }
            //    else
            //    {
            //        // 法向量为x轴
            //        tCount = new int[2] { 1, 0 };
            //    }
            //}
            //else if (interval == 1)
            //{
            //    // interval == 1
            //    tCount = new int[2] { 0, 0 };
            //}
            //else
            //{
            //    // interval == 2
            //    if (nVectorForBS1.X == 0)
            //    {
            //        // 法向量为y轴
            //        tCount = new int[2] { 0, 1 };
            //    }
            //    else
            //    {
            //        // 法向量为x轴
            //        tCount = new int[2] { 1, 0 };
            //    }
            //}

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
            get { return new Guid("2b753b60-fa89-48f2-a3ba-9475cebd0117"); }
        }

        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new CompLabelAttribute(this);
        }

        public class CompLabelAttribute : GH_ComponentAttributes
        {
            public CompLabelAttribute(GhcReorganizeBoundary3 component) : base(component) { }

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

                    using (GH_Capsule capsule1 = GH_Capsule.CreateCapsule(buttonRect1, GH_Palette.Normal))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule1.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("BoundaryTCounts:{0}", ((GhcReorganizeBoundary3)Owner).AllLayerBoundaryTCountsString),
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

                    graphics.DrawString(string.Format("InnerTMaxCounts:{0}", ((GhcReorganizeBoundary3)Owner).AllLayerInnerTMaxCountString),
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