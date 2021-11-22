using Grasshopper.Kernel;
using Plankton;
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
                PlanktonMesh D = dualGraphWithHMDP.DualPlanktonMesh;
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
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    faceIndexsAroundVertex.Add(new List<int>());
                    faceIndexsAroundVertex[i].AddRange(D.Vertices.GetVertexFaces(i));
                    // 去除-1
                    faceIndexsAroundVertex[i].Remove(0);
                }

                #region 区分出每个Segment上的volumeJunction和innerJunction
                List<List<int>> volumeJunctionForEachBoundarySegment = new List<List<int>>();
                List<List<int>> innerJunctionForEachBoundarySegment = new List<List<int>>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    volumeJunctionForEachBoundarySegment.Add(new List<int>());
                    innerJunctionForEachBoundarySegment.Add(new List<int>());

                    for (int j = 0; j < verticesIndexForEachBoundarySegment[i].Count; j++)
                    {
                        if (dualGraphWithHMDP.DVertice_Inner[verticesIndexForEachBoundarySegment[i][j]].Count > 1
                            && dualGraphWithHMDP.DualPlanktonMesh.Vertices.GetVertexFaces(verticesIndexForEachBoundarySegment[i][j]).Contains(-1))
                        {
                            if (dualGraphWithHMDP.DVertice_Volume[verticesIndexForEachBoundarySegment[i][j]].Count == 1)
                            {
                                innerJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                            }
                            else
                            {
                                volumeJunctionForEachBoundarySegment[i].Add(verticesIndexForEachBoundarySegment[i][j]);
                            }
                        }
                    }
                }
                #endregion

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

                    #region 像新分支中添加值
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
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    List<int> dFaces = D.Vertices.GetVertexFaces(i).ToList();
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