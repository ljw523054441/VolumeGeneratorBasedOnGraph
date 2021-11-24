using Grasshopper.Kernel;
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
            base.Message = "开发中 Todo:对偶图中，边界顶点布置方式";

            BoundaryPolylinePoints = new List<Point3d>();
            BoundaryCornerTextDots = new List<TextDot>();
            BoundarySegmentTextDots = new List<TextDot>();

            VolumeJunctionTextDot = new List<TextDot>();
            InnerJunctionTextDot = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> BoundaryPolylinePoints;

        private List<TextDot> BoundaryCornerTextDots;
        private List<TextDot> BoundarySegmentTextDots;

        private List<TextDot> VolumeJunctionTextDot;
        private List<TextDot> InnerJunctionTextDot;

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
            pManager.AddPointParameter("BoundaryVerticesPoints", "BVP", "", GH_ParamAccess.list);
            
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
                DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);
                // PlanktonMesh D = dualGraphWithHMDP.DualPlanktonMesh;
                // PlanktonMesh I = dualGraphWithHM.IntegrateDualPlanktonMesh;
                List<GraphNode> decomposedPGraphNodes = dualGraphWithHMDP.GraphNodes;

                //SortedDictionary<int, GraphNode> volume_volumeNode = dualGraphWithHM.Volume_VolumeNode;

                //List<GraphNode> undividedPGraphNodes = dualGraphWithHM.UndividedGraphNodes;
                //List<List<int>> undividedPGraphTables = dualGraphWithHM.UndividedGraphTable;



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

                #region 为BoundarySegment的IncludedDVertice属性赋值
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    sortedBoundarySegments[i].IncludedDVertice = new List<int>();
                    sortedBoundarySegments[i].IncludedDVertice.AddRange(verticesIndexForEachBoundarySegment[i]);
                }
                #endregion

                #endregion

                List<List<int>> faceIndexsAroundVertex = new List<List<int>>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(0);
                }

                //#region 区分出每个Segment上的volumeJunction和innerJunction
                //List<int> volumeJunctionList = new List<int>();
                //List<int> innerJunctionList = new List<int>();

                //List<List<int>> volumeJunctionForEachBoundarySegment = new List<List<int>>();
                //List<List<int>> innerJunctionForEachBoundarySegment = new List<List<int>>();
                //for (int i = 0; i < sortedBoundarySegments.Count; i++)
                //{
                //    volumeJunctionForEachBoundarySegment.Add(new List<int>());
                //    innerJunctionForEachBoundarySegment.Add(new List<int>());

                //    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                //    {
                //        if (dualGraphWithHMDP.DVertice_Inner[verticesIndexForEachBoundarySegment[i][j]].Count > 1
                //            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBoundarySegment[i][j]).Contains(-1))
                //        {
                //            if (dualGraphWithHMDP.DVertice_Volume[verticesIndexForEachBoundarySegment[i][j]].Count == 1)
                //            {
                //                innerJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                //                innerJunctionList.Add(verticesIndexForEachBoundarySegment[i][j]);
                //            }
                //            else
                //            {
                //                volumeJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                //                volumeJunctionList.Add(verticesIndexForEachBoundarySegment[i][j]);
                //            }
                //        }
                //    }
                //}
                //#endregion

                List<List<int>> juntionsIndexForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    // 深拷贝
                    juntionsIndexForEachBoundarySegment.Add(new List<int>());
                    juntionsIndexForEachBoundarySegment[i].AddRange(verticesIndexForEachBoundarySegment[i]);
                    // 去除首尾顶点
                    juntionsIndexForEachBoundarySegment[i].Remove(0);
                    juntionsIndexForEachBoundarySegment[i].Remove(juntionsIndexForEachBoundarySegment[i].Count - 1);
                }

                #region 对于对偶图边界上所有的Junction点进行分裂
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    // 分裂后输入的dualGraphWithHMDP改变
                    DecomposeJunction(dualGraphWithHMDP, juntionsIndexForEachBoundarySegment[i]);
                }
                #endregion

                // 此时边界上的点，全部都是juntion点（除去每条边的首尾）



                #region 分割boundarySegment
                // 分裂每个Segment
                List<List<List<int>>> innerJunctionForEachSubSegmentForAllSegment = new List<List<List<int>>>();

                List<List<Point3d>> volumeJunctionPointForAllSegment = new List<List<Point3d>>();
                List<List<List<int>>> volumeJunctionCorrespondingWhichVolumeForAllSegment = new List<List<List<int>>>();
                List<List<List<Point3d>>> innerJunctionPointForAllSegment = new List<List<List<Point3d>>>();
                List<List<List<List<int>>>> innerJunctionCorrespondingWhichInnerForAllSegment = new List<List<List<List<int>>>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    ITree<BoundarySegment> currentBoundary = NodeTree<BoundarySegment>.NewTree();
                    INode<BoundarySegment> top = currentBoundary.AddChild(sortedBoundarySegments[i]);

                    // SortedDictionary<int, List<int>> volumeJunction_itsSplit = new SortedDictionary<int, List<int>>();
                    DecomposeVolumeJunction(dualGraphWithHMDP, volumeJunctionForEachBoundarySegment[i]);

                    // List<string> debugPrint = UtilityFunctions.PrintFacesVertices(dualGraphWithHMDP.DualPlanktonMesh);

                    List<List<int>> innerJunctionForEachSubSegment;

                    List<Point3d> volumeJunctionPointForCurrentSegment;
                    List<List<int>> volumeJunctionCorrespondingWhichInnerForCurrentSegment;
                    List<List<int>> volumeJunctionCorrespondingWhichVolumeForCurrentSegment;
                    List<List<Point3d>> innerJunctionPointForCurrentSegmentAllSubSegment;
                    List<List<List<int>>> innerJunctionCorrespondingWhichInnerForCurrentSegmentAllSubSegment;

                    DecomposeSegmentWithVolumeJunctions(top, 
                                                        dualGraphWithHMDP, 
                                                        volumeJunctionForEachBoundarySegment[i],
                                                        innerJunctionForEachBoundarySegment[i],
                                                        out volumeJunctionPointForCurrentSegment,
                                                        out volumeJunctionCorrespondingWhichInnerForCurrentSegment,
                                                        out volumeJunctionCorrespondingWhichVolumeForCurrentSegment,
                                                        out innerJunctionForEachSubSegment,
                                                        out innerJunctionPointForCurrentSegmentAllSubSegment,
                                                        out innerJunctionCorrespondingWhichInnerForCurrentSegmentAllSubSegment);
                    #region 添加新分支
                    innerJunctionForEachSubSegmentForAllSegment.Add(new List<List<int>>());

                    volumeJunctionPointForAllSegment.Add(new List<Point3d>());
                    volumeJunctionCorrespondingWhichVolumeForAllSegment.Add(new List<List<int>>());
                    innerJunctionPointForAllSegment.Add(new List<List<Point3d>>());
                    innerJunctionCorrespondingWhichInnerForAllSegment.Add(new List<List<List<int>>>());
                    #endregion

                    #region 向新分支中添加值
                    for (int j = 0; j < innerJunctionForEachSubSegment.Count; j++)
                    {
                        innerJunctionForEachSubSegmentForAllSegment[i].Add(new List<int>());
                        innerJunctionForEachSubSegmentForAllSegment[i][j].AddRange(innerJunctionForEachSubSegment[j]);
                    }

                    volumeJunctionPointForAllSegment[i].AddRange(volumeJunctionPointForCurrentSegment);
                    for (int j = 0; j < volumeJunctionCorrespondingWhichInnerForCurrentSegment.Count; j++)
                    {
                        volumeJunctionCorrespondingWhichVolumeForAllSegment[i].Add(new List<int>());
                        volumeJunctionCorrespondingWhichVolumeForAllSegment[i][j].AddRange(volumeJunctionCorrespondingWhichInnerForCurrentSegment[j]);
                    }
                    for (int j = 0; j < innerJunctionPointForCurrentSegmentAllSubSegment.Count; j++)
                    {
                        innerJunctionPointForAllSegment[i].Add(new List<Point3d>());
                        innerJunctionPointForAllSegment[i][j].AddRange(innerJunctionPointForCurrentSegmentAllSubSegment[j]);
                    }
                    for (int j = 0; j < innerJunctionCorrespondingWhichInnerForCurrentSegmentAllSubSegment.Count; j++)
                    {
                        innerJunctionCorrespondingWhichInnerForAllSegment[i].Add(new List<List<int>>());
                        for (int k = 0; k < innerJunctionCorrespondingWhichInnerForCurrentSegmentAllSubSegment[j].Count; k++)
                        {
                            innerJunctionCorrespondingWhichInnerForAllSegment[i][j].Add(new List<int>());
                            innerJunctionCorrespondingWhichInnerForAllSegment[i][j][k].AddRange(innerJunctionCorrespondingWhichInnerForCurrentSegmentAllSubSegment[j][k]);
                        }
                    }
                    #endregion
                }
                #endregion


                List<int> boundaryVerticesOnDual = new List<int>();
                for (int i = 0; i < dualGraphWithHMDP.DualPlanktonMesh.Vertices.Count; i++)
                {
                    if (dualGraphWithHMDP.DualPlanktonMesh.Vertices.NakedEdgeCount(i) >= 2)
                    {
                        boundaryVerticesOnDual.Add(i);
                    }
                }

                List<Point3d> boundaryVerticesPoints = new List<Point3d>();
                for (int i = 0; i < boundaryVerticesOnDual.Count; i++)
                {
                    int viIndex;
                    int vjIndex;
                    int iiIndex;
                    int ijIndex;
                    int ikIndex;
                    // 如果是VolumeJunction点
                    if (IsInVolumeJunctions(volumeJunctionForEachBoundarySegment,
                                            boundaryVerticesOnDual[i],
                                            out viIndex,
                                            out vjIndex))
                    {
                        boundaryVerticesPoints.Add(volumeJunctionPointForAllSegment[viIndex][vjIndex]);
                    }
                    // 如果是InnerJunction点
                    else if (IsInInnerJunctions(innerJunctionForEachSubSegmentForAllSegment,
                                                boundaryVerticesOnDual[i],
                                                out iiIndex,
                                                out ijIndex,
                                                out ikIndex))
                    {
                        boundaryVerticesPoints.Add(innerJunctionPointForAllSegment[viIndex][vjIndex][ikIndex]);
                    }
                    // 如果是边界上的角点
                    else
                    {
                        boundaryVerticesPoints.Add(boundaryCorner[sortedBoundaryCornerIndexs.IndexOf(boundaryVerticesOnDual[i])]);
                    }
                }
                DA.SetDataList("BoundaryVerticesPoints", boundaryVerticesPoints);



                #region 可视化部分
                BoundaryPolylinePoints.Clear();
                BoundaryCornerTextDots.Clear();
                BoundarySegmentTextDots.Clear();

                VolumeJunctionTextDot.Clear();
                InnerJunctionTextDot.Clear();

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

                #region 设置volumeJunction
                for (int i = 0; i < volumeJunctionPointForAllSegment.Count; i++)
                {
                    for (int j = 0; j < volumeJunctionPointForAllSegment[i].Count; j++)
                    {
                        int currentVolumeJunctionIndex = volumeJunctionForEachBoundarySegment[i][j];
                        List<int> currentVolumeJunctionCorrespondingWhichVolume = volumeJunctionCorrespondingWhichVolumeForAllSegment[i][j];
                        string str = string.Join<int>(";", currentVolumeJunctionCorrespondingWhichVolume);
                        TextDot volumeJunctionTextDot = new TextDot(string.Format("{0}|{1}", currentVolumeJunctionIndex, str), volumeJunctionPointForAllSegment[i][j]);
                        VolumeJunctionTextDot.Add(volumeJunctionTextDot);
                    }
                }
                #endregion

                #region 设置innerJunction
                for (int i = 0; i < innerJunctionPointForAllSegment.Count; i++)
                {
                    for (int j = 0; j < innerJunctionPointForAllSegment[i].Count; j++)
                    {
                        for (int k = 0; k < innerJunctionPointForAllSegment[i][j].Count; k++)
                        {
                            int currentInnerJunctionIndex = innerJunctionForEachSubSegmentForAllSegment[i][j][k];
                            List<int> currentInnerJunctionCorrespondingWhichInner = innerJunctionCorrespondingWhichInnerForAllSegment[i][j][k];
                            string str = string.Join<int>(";", currentInnerJunctionCorrespondingWhichInner);
                            TextDot innerJunctionTextDot = new TextDot(string.Format("{0}|{1}", currentInnerJunctionIndex, str), innerJunctionPointForAllSegment[i][j][k]);
                            InnerJunctionTextDot.Add(innerJunctionTextDot);
                        }
                    }
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

        private bool IsInVolumeJunctions(List<List<int>> volumeJunctionForEachBoundarySegment, 
                                         int boundaryVertex,
                                         out int iIndex,
                                         out int jIndex)
        {
            bool flag = false;
            iIndex = -1;
            jIndex = -1;
            for (int i = 0; i < volumeJunctionForEachBoundarySegment.Count; i++)
            {
                if (volumeJunctionForEachBoundarySegment[i].Contains(boundaryVertex))
                {
                    flag = true;
                    iIndex = i;
                    jIndex = volumeJunctionForEachBoundarySegment[i].IndexOf(boundaryVertex);
                    break;
                }
            }
            return flag;
        }

        private bool IsInInnerJunctions(List<List<List<int>>> innerJunctionForEachBoundarySegment,
                                        int boundaryVertex,
                                        out int iIndex,
                                        out int jIndex,
                                        out int kIndex)
        {
            bool flag = false;
            iIndex = -1;
            jIndex = -1;
            kIndex = -1;
            for (int i = 0; i < innerJunctionForEachBoundarySegment.Count; i++)
            {
                for (int j = 0; j < innerJunctionForEachBoundarySegment[i].Count; j++)
                {
                    if (innerJunctionForEachBoundarySegment[i][j].Contains(boundaryVertex))
                    {
                        flag = true;
                        iIndex = i;
                        jIndex = j;
                        kIndex = innerJunctionForEachBoundarySegment[i][j].IndexOf(boundaryVertex);
                    }
                }
                
            }
            return flag;
        }

        private void DecomposeJunction(//INode<BoundarySegment> INodeToSplit,
                                      DualGraphWithHM dualGraphWithHM,
                                      List<int> junctionForCurrentSegment)
                                      // SortedDictionary<int, List<int>> volumeJunction_itsSplit)
        {
            // 实际上经过此函数后，dualGraphWithHM就已经被修改了

            // 判断是否需要分裂
            for (int i = 0; i < junctionForCurrentSegment.Count; i++)
            {
                List<int> adjacentFace = dualGraphWithHM.DualPlanktonMesh.Vertices.GetVertexFaces(junctionForCurrentSegment[i]).ToList();
                // adjacentFace.Remove(-1);
                adjacentFace.Reverse();

                // 与adjacentFace对应的该顶点的出半边的列表
                List<int> correspondingH = new List<int>();
                List<int> hs = dualGraphWithHM.DualPlanktonMesh.Vertices.GetHalfedges(junctionForCurrentSegment[i]).ToList();
                for (int j = 0; j < adjacentFace.Count; j++)
                {
                    for (int k = 0; k < hs.Count; k++)
                    {
                        if (dualGraphWithHM.DualPlanktonMesh.Halfedges[hs[k]].AdjacentFace == adjacentFace[j])
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
                        if (dualGraphWithHM.DualPlanktonMesh.Halfedges[hs[k]].AdjacentFace == adjacentFace[j])
                        {
                            correspondingHEnd.Add(dualGraphWithHM.DualPlanktonMesh.Halfedges.EndVertex(hs[k]));
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
                    List<int> newVertice = new List<int>();
                    while (iter < correspondingHEnd.Count - 3)
                    {
                        int[,] currentHPairIndexs = new int[2, 2] {
                                /* 初始化索引号为0的行 */
                                { junctionForCurrentSegment[i],
                                  newAddedHalfedgeEndStartFromOrigin},
                                /* 初始化索引号为1的行 */
                                { junctionForCurrentSegment[i],
                                  correspondingHEnd[1 + iter]}};

                        dualGraphWithHM.DualPlanktonMesh = SplitVertex(dualGraphWithHM.DualPlanktonMesh,
                                                                           currentHPairIndexs,
                                                                           junctionForCurrentSegment[i],
                                                                           out newVertexIndex,
                                                                           out newAddedHalfedgeEndStartFromOrigin);

                        List<string> debugPrint = UtilityFunctions.PrintFacesVertices(dualGraphWithHM.DualPlanktonMesh);
                        newVertice.Add(newVertexIndex);

                        iter++;
                    }
                    // volumeJunction_itsSplit.Add(volumeJunctionForCurrentSegment[i], newVertice);
                }
                else
                {
                    continue;
                }
            }

            // 将分裂的结果也同样应用到IntegratedPlanktonMesh上

            // 更新dualGraphWithHM的DVertice_Inner属性
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

        private void SplitBoundarySegment(BoundarySegment currentBoundarySegment,
                                          List<int> junctionsOnCurrentBoundary,
                                          List<int> verticeOnCurrentBoundary,
                                          DualGraphWithHM dualGraphWithHMDP)
        {
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
                dVerticesForEachSubSegment.Add(subList);
            }
            #endregion

            #region 计算这段BoundarySegment上所有的Junction点所对应的t，及Point点
            // 计算分割点的t
            List<List<int>> JunctionCorrespondingWhichInnerForCurrentSubSegment = new List<List<int>>();

            List<double> innerAreaProportionForCurrentSegment = new List<double>();
            for (int i = 0; i < junctionsOnCurrentBoundary.Count; i++)
            {
                List<int> value = dualGraphWithHMDP.DVertice_Inner[junctionsOnCurrentBoundary[i]];

                JunctionCorrespondingWhichInnerForCurrentSubSegment.Add(new List<int>());


            }
            #endregion
        }


        public void DecomposeSegmentWithVolumeJunctions(INode<BoundarySegment> INodeToSplit, 
                              DualGraphWithHM dualGraphWithHM, 
                              List<int> volumeJunctionForCurrentSegment,
                              List<int> innerJunctionForCurrentSegment,
                              out List<Point3d> volumeJunctionPointForCurrentSegment,
                              out List<List<int>> volumeJunctionCorrespondingWhichInner,
                              out List<List<int>> volumeJunctionCorrespondingWhichVolume,
                              out List<List<int>> innerJunctionForEachSubSegment,
                              out List<List<Point3d>> innerJunctionPointForAllSubSegment,
                              out List<List<List<int>>> innerJunctionCorrespondingWhichInner) 
        {
            if (volumeJunctionForCurrentSegment.Count > 0)
            {
                List<int> volumeJunctionForCurrentSegmentDP = new List<int>();
                volumeJunctionForCurrentSegmentDP.AddRange(volumeJunctionForCurrentSegment);
                List<int> verticesIndexForCurrentSegmentDP = new List<int>();
                verticesIndexForCurrentSegmentDP.AddRange(INodeToSplit.Data.IncludedDVertice);

                #region 计算分割后每段子Segment上所包含的顶点
                // 计算分割后每段的顶点
                List<int> volumeJunctionIndexs = new List<int>();
                for (int i = 0; i < volumeJunctionForCurrentSegmentDP.Count; i++)
                {
                    int debug = volumeJunctionForCurrentSegmentDP[i];
                    volumeJunctionIndexs.Add(verticesIndexForCurrentSegmentDP.IndexOf(debug));
                }
                // 构造[0,volumeJunctionIndex,...,volumeJunctionForCurrentSegmentDP.Count - 1]的列表，用来在下面截取volumeJunctionForCurrentSegmentDP列表
                volumeJunctionIndexs.Insert(0, 0);
                volumeJunctionIndexs.Add(verticesIndexForCurrentSegmentDP.Count - 1);
                List<int[]> pairIndex = new List<int[]>();
                for (int i = 0; i < volumeJunctionIndexs.Count - 1; i++)
                {
                    int[] pair = new int[] { volumeJunctionIndexs[i], volumeJunctionIndexs[i + 1] };
                    pairIndex.Add(pair);
                }
                // 截取volumeJunctionForCurrentSegmentDP列表
                List<List<int>> dVerticesForEachSubSegment = new List<List<int>>();
                for (int i = 0; i < pairIndex.Count; i++)
                {
                    List<int> subList = verticesIndexForCurrentSegmentDP.Skip(pairIndex[i][0]).Take(pairIndex[i][1] - pairIndex[i][0] + 1).ToList();
                    dVerticesForEachSubSegment.Add(subList);
                }
                #endregion

                #region 每段子Segment上所包含的InnerJunction
                innerJunctionForEachSubSegment = new List<List<int>>();
                for (int i = 0; i < dVerticesForEachSubSegment.Count; i++)
                {
                    innerJunctionForEachSubSegment.Add(new List<int>());
                    for (int j = 0; j < dVerticesForEachSubSegment[i].Count; j++)
                    {
                        if (dualGraphWithHM.DVertice_Inner[dVerticesForEachSubSegment[i][j]].Count > 1
                            && !volumeJunctionForCurrentSegmentDP.Contains(dVerticesForEachSubSegment[i][j]))
                        {
                            innerJunctionForEachSubSegment[i].Add(dVerticesForEachSubSegment[i][j]);
                        }
                    }
                }
                #endregion


                #region 计算这段Segment上VolumeJunction点所对应的t，及Point点
                // 计算分割点的t
                volumeJunctionCorrespondingWhichInner = new List<List<int>>();
                volumeJunctionCorrespondingWhichVolume = new List<List<int>>();

                List<double> volumeAreaProportionForCurrentSegment = new List<double>();
                for (int i = 0; i < volumeJunctionForCurrentSegmentDP.Count; i++)
                {
                    List<int> correspondingVolumeIndex = dualGraphWithHM.DVertice_Volume[volumeJunctionForCurrentSegmentDP[i]];

                    volumeJunctionCorrespondingWhichInner.Add(new List<int>());
                    volumeJunctionCorrespondingWhichVolume.Add(new List<int>());
                    // 对于分割两个或几个Volume的边缘点来说
                    switch (correspondingVolumeIndex.Count)
                    {
                        // 如果分割2个volume
                        case 2:
                            int vface1 = -1;
                            int vface2 = -1;

                            int dface1 = -1;
                            int dface2 = -1;
                            for (int j = 0; j < dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.Count; j++)
                            {
                                if (dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[j].StartVertex == volumeJunctionForCurrentSegmentDP[i]
                                    && dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.EndVertex(j) == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(volumeJunctionForCurrentSegmentDP[i]) - 1])
                                {
                                    vface1 = dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[j].AdjacentFace;
                                }
                                if (dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[j].StartVertex == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(volumeJunctionForCurrentSegmentDP[i]) + 1]
                                    && dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges.EndVertex(j) == volumeJunctionForCurrentSegmentDP[i])
                                {
                                    vface2 = dualGraphWithHM.IntegrateDualPlanktonMesh.Halfedges[j].AdjacentFace;
                                }
                            }
                            for (int j = 0; j < dualGraphWithHM.DualPlanktonMesh.Halfedges.Count; j++)
                            {
                                if (dualGraphWithHM.DualPlanktonMesh.Halfedges[j].StartVertex == volumeJunctionForCurrentSegmentDP[i]
                                    && dualGraphWithHM.DualPlanktonMesh.Halfedges.EndVertex(j) == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(volumeJunctionForCurrentSegmentDP[i]) - 1])
                                {
                                    dface1 = dualGraphWithHM.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                                }
                                if (dualGraphWithHM.DualPlanktonMesh.Halfedges[j].StartVertex == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(volumeJunctionForCurrentSegmentDP[i]) + 1]
                                    && dualGraphWithHM.DualPlanktonMesh.Halfedges.EndVertex(j) == volumeJunctionForCurrentSegmentDP[i])
                                {
                                    dface2 = dualGraphWithHM.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                                }
                            }

                            int volume1 = -1;
                            int volume2 = -1;
                            foreach (var item in dualGraphWithHM.Volume_VFace)
                            {
                                if (vface1 == item.Value)
                                {
                                    volume1 = item.Key;
                                }
                                if (vface2 == item.Value)
                                {
                                    volume2 = item.Key;
                                }
                            }

                            volumeJunctionCorrespondingWhichVolume[i].Add(volume1);
                            volumeJunctionCorrespondingWhichVolume[i].Add(volume2);

                            int inner1 = -1;
                            int inner2 = -1;
                            foreach (var item in dualGraphWithHM.Inner_DFace)
                            {
                                if (dface1 == item.Value)
                                {
                                    inner1 = item.Key;
                                }
                                if (dface2 == item.Value)
                                {
                                    inner2 = item.Key;
                                }
                            }

                            volumeJunctionCorrespondingWhichInner[i].Add(inner1);
                            volumeJunctionCorrespondingWhichInner[i].Add(inner2);

                            double volume1AreaProportion = dualGraphWithHM.UndividedGraphNodes[volume1].NodeAttribute.NodeAreaProportion;
                            double volume2AreaProportion = dualGraphWithHM.UndividedGraphNodes[volume2].NodeAttribute.NodeAreaProportion;

                            if (i == 0)
                            {
                                volumeAreaProportionForCurrentSegment.Add(volume1AreaProportion);
                                volumeAreaProportionForCurrentSegment.Add(volume2AreaProportion);
                            }
                            else
                            {
                                volumeAreaProportionForCurrentSegment.Add(volume2AreaProportion);
                            }

                            break;
                        // 如果分割3个volume
                        case 3:
                            // 那么就涉及到分裂这个volumeJunction点的问题了
                            break;
                        default:
                            break;

                    }
                }

                // 计算t及Point点
                List<double> tList = new List<double>();
                volumeJunctionPointForCurrentSegment = new List<Point3d>();
                double sum = volumeAreaProportionForCurrentSegment.Sum();
                for (int i = 0; i < volumeJunctionForCurrentSegmentDP.Count; i++)
                {
                    double t = volumeAreaProportionForCurrentSegment[i] / sum;
                    Point3d tPoint = INodeToSplit.Data.Line.PointAt(t);
                    tList.Add(t);
                     volumeJunctionPointForCurrentSegment.Add(tPoint);
                }
                #endregion

                #region 分割生成子Segment，并生成叶子节点
                Curve segmentLineCurve = INodeToSplit.Data.Line.ToNurbsCurve();

                Curve[] subSegmentLineCurve = segmentLineCurve.Split(tList);

                List<Line> subLines = new List<Line>();
                foreach (var item in subSegmentLineCurve)
                {
                    Point3d start = item.PointAtStart;
                    Point3d end = item.PointAtEnd;
                    subLines.Add(new Line(start, end));
                }

                innerJunctionPointForAllSubSegment = new List<List<Point3d>>();
                innerJunctionCorrespondingWhichInner = new List<List<List<int>>>();
                for (int i = 0; i < subLines.Count; i++)
                {
                    string label = string.Format("{0}{1}", INodeToSplit.Data.Label, i);
                    BoundarySegment subSegment = new BoundarySegment(subLines[i], label);
                    subSegment.IncludedDVertice = new List<int>();
                    subSegment.IncludedDVertice.AddRange(dVerticesForEachSubSegment[i]);
                    INode<BoundarySegment> childBoundarySegment = INodeToSplit.AddChild(subSegment);

                    innerJunctionPointForAllSubSegment.Add(new List<Point3d>());
                    innerJunctionCorrespondingWhichInner.Add(new List<List<int>>());
                    List<Point3d> innerJunctionPointForCurrentSubSegment;
                    List<List<int>> innerJunctionCorrespondingWhichInnerForCurrentSubSegment;
                    DecomposeSubSegmentWithInnerJunctions(childBoundarySegment,
                                                          dualGraphWithHM,
                                                          innerJunctionForEachSubSegment[i],
                                                          out innerJunctionPointForCurrentSubSegment,
                                                          out innerJunctionCorrespondingWhichInnerForCurrentSubSegment);
                    innerJunctionPointForAllSubSegment[i].AddRange(innerJunctionPointForCurrentSubSegment);
                    for (int j = 0; j < innerJunctionCorrespondingWhichInnerForCurrentSubSegment.Count; j++)
                    {
                        innerJunctionCorrespondingWhichInner[i].Add(new List<int>());
                        innerJunctionCorrespondingWhichInner[i][j].AddRange(innerJunctionCorrespondingWhichInnerForCurrentSubSegment[j]);
                    }
                }
            }
            else if (innerJunctionForCurrentSegment.Count > 0)
            {
                // out
                volumeJunctionPointForCurrentSegment = new List<Point3d>();
                volumeJunctionCorrespondingWhichInner = new List<List<int>>();
                volumeJunctionCorrespondingWhichVolume = new List<List<int>>();
                innerJunctionForEachSubSegment = new List<List<int>>();

                innerJunctionForEachSubSegment.Add(new List<int>());
                innerJunctionForEachSubSegment.Last().AddRange(innerJunctionForCurrentSegment);

                innerJunctionPointForAllSubSegment = new List<List<Point3d>>();
                innerJunctionCorrespondingWhichInner = new List<List<List<int>>>();

                List<Point3d> innerJunctionPointForCurrentSegment;
                List<List<int>> innerJunctionCorrespondingWhichInnerForCurrentSegment;
                DecomposeSubSegmentWithInnerJunctions(INodeToSplit,
                                                      dualGraphWithHM,
                                                      innerJunctionForCurrentSegment,
                                                      out innerJunctionPointForCurrentSegment,
                                                      out innerJunctionCorrespondingWhichInnerForCurrentSegment);

                innerJunctionPointForAllSubSegment.Add(new List<Point3d>());
                innerJunctionCorrespondingWhichInner.Add(new List<List<int>>());
                innerJunctionPointForAllSubSegment.Last().AddRange(innerJunctionPointForCurrentSegment);
                for (int i = 0; i < innerJunctionCorrespondingWhichInnerForCurrentSegment.Count; i++)
                {
                    innerJunctionCorrespondingWhichInner.Last().Add(new List<int>());
                    innerJunctionCorrespondingWhichInner.Last()[i].AddRange(innerJunctionCorrespondingWhichInnerForCurrentSegment[i]);
                }
            }
            else
            {
                volumeJunctionPointForCurrentSegment = new List<Point3d>();
                volumeJunctionCorrespondingWhichInner = new List<List<int>>();
                volumeJunctionCorrespondingWhichVolume = new List<List<int>>();
                innerJunctionForEachSubSegment = new List<List<int>>();
                innerJunctionPointForAllSubSegment = new List<List<Point3d>>();
                innerJunctionCorrespondingWhichInner = new List<List<List<int>>>();
            }
            #endregion
            return;
        }

        public void DecomposeSubSegmentWithInnerJunctions(INode<BoundarySegment> INodeToSplit,
                                            DualGraphWithHM dualGraphWithHM,
                                            List<int> innerJunctionForCurrentSubSegment,
                                            out List<Point3d> innerJunctionPointForCurrentSubSegment,
                                            out List<List<int>> innerJunctionCorrespondingWhichInnerForCurrentSubSegment)
        {
            List<int> innerJunctionForCurrentSegmentDP = new List<int>();
            innerJunctionForCurrentSegmentDP.AddRange(innerJunctionForCurrentSubSegment);
            List<int> verticesIndexForCurrentSegmentDP = new List<int>();
            verticesIndexForCurrentSegmentDP.AddRange(INodeToSplit.Data.IncludedDVertice);

            #region 计算分割后每段子Segment上所包含的顶点
            // 计算分割后每段子Segment上所包含的顶点
            List<int> innerJunctionIndexs = new List<int>();
            for (int i = 0; i < innerJunctionForCurrentSegmentDP.Count; i++)
            {
                innerJunctionIndexs.Add(verticesIndexForCurrentSegmentDP.IndexOf(innerJunctionForCurrentSegmentDP[i]));
            }
            // 构造[0,innerJunctionIndex,...,innerJunctionForCurrentSegmentDP.Count - 1]的列表，用来在下面截取innerJunctionForCurrentSegmentDP列表
            innerJunctionIndexs.Insert(0, 0);
            innerJunctionIndexs.Add(verticesIndexForCurrentSegmentDP.Count - 1);
            List<int[]> pairIndex = new List<int[]>();
            for (int i = 0; i < innerJunctionIndexs.Count - 1; i++)
            {
                int[] pair = new int[] { innerJunctionIndexs[i], innerJunctionIndexs[i + 1] };
                pairIndex.Add(pair);
            }
            // 截取innerJunctionForCurrentSegmentDP列表
            List<List<int>> dVerticesForEachSubSegment = new List<List<int>>();
            for (int i = 0; i < pairIndex.Count; i++)
            {
                List<int> subList = verticesIndexForCurrentSegmentDP.Skip(pairIndex[i][0]).Take(pairIndex[i][1] - pairIndex[i][0] + 1).ToList();
                dVerticesForEachSubSegment.Add(subList);
            }
            #endregion

            #region 计算这段Segment上innerJunction点所对应的t，及Point点
            // 计算分割点的t
            innerJunctionCorrespondingWhichInnerForCurrentSubSegment = new List<List<int>>();

            List<double> innerAreaProportionForCurrentSegment = new List<double>();
            for (int i = 0; i < innerJunctionForCurrentSegmentDP.Count; i++)
            {
                List<int> value = dualGraphWithHM.DVertice_Inner[innerJunctionForCurrentSegmentDP[i]];

                innerJunctionCorrespondingWhichInnerForCurrentSubSegment.Add(new List<int>());
                // 对于分割两个或几个Inner的边缘点来说
                switch (value.Count)
                {
                    case 2:
                        int dface1 = -1;
                        int dface2 = -1;

                        for (int j = 0; j < dualGraphWithHM.DualPlanktonMesh.Halfedges.Count; j++)
                        {
                            if (dualGraphWithHM.DualPlanktonMesh.Halfedges[j].StartVertex == innerJunctionForCurrentSegmentDP[i] 
                                && dualGraphWithHM.DualPlanktonMesh.Halfedges.EndVertex(j) == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(innerJunctionForCurrentSegmentDP[i]) - 1])
                            {
                                dface1 = dualGraphWithHM.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                            }
                            if (dualGraphWithHM.DualPlanktonMesh.Halfedges[j].StartVertex == verticesIndexForCurrentSegmentDP[verticesIndexForCurrentSegmentDP.IndexOf(innerJunctionForCurrentSegmentDP[i]) + 1]
                                && dualGraphWithHM.DualPlanktonMesh.Halfedges.EndVertex(j) == innerJunctionForCurrentSegmentDP[i])
                            {
                                dface2 = dualGraphWithHM.DualPlanktonMesh.Halfedges[j].AdjacentFace;
                            }
                        }

                        int inner1 = -1;
                        int inner2 = -1;
                        foreach (var item in dualGraphWithHM.Inner_DFace)
                        {
                            if (dface1 == item.Value)
                            {
                                inner1 = item.Key;
                            }
                            if (dface2 == item.Value)
                            {
                                inner2 = item.Key;
                            }
                        }

                        innerJunctionCorrespondingWhichInnerForCurrentSubSegment[i].Add(inner1);
                        innerJunctionCorrespondingWhichInnerForCurrentSubSegment[i].Add(inner2);

                        double inner1AreaProportion = dualGraphWithHM.GraphNodes[inner1].NodeAttribute.NodeAreaProportion;
                        double inner2AreaProportion = dualGraphWithHM.GraphNodes[inner2].NodeAttribute.NodeAreaProportion;

                        if (i == 0)
                        {
                            innerAreaProportionForCurrentSegment.Add(inner1AreaProportion);
                            innerAreaProportionForCurrentSegment.Add(inner2AreaProportion);
                        }
                        else
                        {
                            innerAreaProportionForCurrentSegment.Add(inner2AreaProportion);
                        }

                        break;
                    // 如果分割3个inner
                    case 3:
                        // 要思考怎么办
                        break;
                    default:
                        break;
                }
            }

            // 计算t及Point点
            List<double> tList = new List<double>();
            innerJunctionPointForCurrentSubSegment = new List<Point3d>();
            double sum = innerAreaProportionForCurrentSegment.Sum();
            for (int i = 0; i < innerJunctionForCurrentSegmentDP.Count; i++)
            {
                double t = innerAreaProportionForCurrentSegment[i] / sum;
                Point3d tPoint = INodeToSplit.Data.Line.PointAt(t);
                tList.Add(t);
                innerJunctionPointForCurrentSubSegment.Add(tPoint);
            }
            #endregion

            #region 分割生成子Segment，并生成叶子节点
            Curve segmentLineCurve = INodeToSplit.Data.Line.ToNurbsCurve();

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
                string label = string.Format("{0}{1}", INodeToSplit.Data.Label, i);
                BoundarySegment subSegment = new BoundarySegment(subLines[i], label);
                subSegment.IncludedDVertice = new List<int>();
                subSegment.IncludedDVertice.AddRange(dVerticesForEachSubSegment[i]);
                INode<BoundarySegment> childBoundarySegment = INodeToSplit.AddChild(subSegment);
            }
            #endregion

            return;
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

            for (int i = 0; i < VolumeJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.Red, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < InnerJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerJunctionTextDot[i], Color.DarkOrange, Color.White, Color.White);
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
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.Red, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < InnerJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerJunctionTextDot[i], Color.DarkOrange, Color.White, Color.White);
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