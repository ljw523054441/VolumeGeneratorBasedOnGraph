using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDecomposeVertex : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDecomposeVertex class.
        /// </summary>
        public GhcDecomposeVertex()
          : base("DecomposeVertex", "DecomposeVertex",
              "拆解度大于4的顶点",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            PRhinoMesh = new Mesh();
            PHalfedgeDottedCurves = new List<Curve>();

            GraphEdges = new List<Line>();

            InnerNodeTextDot = new List<TextDot>();
            OuterNodeTextDot = new List<TextDot>();


        }

        private int Thickness;

        /// <summary>
        /// 用RhinoMesh来表示的结果PlanktonMesh
        /// </summary>
        private Mesh PRhinoMesh;
        /// <summary>
        /// 用虚线来表示结果PlanktonMesh的边
        /// </summary>
        private List<Curve> PHalfedgeDottedCurves;

        /// <summary>
        /// 图结构的边
        /// </summary>
        private List<Line> GraphEdges;

        private List<TextDot> InnerNodeTextDot;
        private List<TextDot> OuterNodeTextDot;



        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);

            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddIntegerParameter("IndexOfSplittedHalfedgeMesh", "I", "某一种顶点分裂后的结果", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("NewTriangleHalfedgeMesh", "NTHMesh", "新生成的三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddGenericParameter("NewGraph", "NG", "描述VolumeNode和BoundaryNode的所有连接关系的图结构(包括新分裂产生的点)", GH_ParamAccess.tree);
            pManager.AddGenericParameter("NewGraphNode", "NGN", "更新后的图结构中的节点", GH_ParamAccess.list);

            pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);

            pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);
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

            PlanktonMesh P = new PlanktonMesh();

            List<Node> pNodes = new List<Node>();

            PlanktonMesh edgeSplitedP = new PlanktonMesh();
            List<PlanktonMesh> edgeResetStartP = new List<PlanktonMesh>();

            GH_Structure<GH_Integer> gh_Structure_graph = null;
            DataTree<int> graph = new DataTree<int>();
            List<List<int>> pGraphLoL = new List<List<int>>();

            int index = 0;

            Thickness = 2;

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh", ref P)
                && DA.GetDataList<Node>("GraphNode", pNodes)
                && DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_graph))
            {
                // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>，再转化为LoL
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_graph, ref graph);
                pGraphLoL = UtilityFunctions.DataTreeToLoL<int>(graph);

                // 利用PlanktonMesh为复制准备的构造函数，进行深拷贝
                PlanktonMesh PDeepCopy = new PlanktonMesh(P);

                #region 存储所有结果的列表，包括图结构，图节点，PlanktonMesh
                List<List<List<int>>> allPossiblePGraphLoL = new List<List<List<int>>>();
                List<List<Node>> allPossiblePNodes = new List<List<Node>>();
                List<List<PlanktonMesh>> allPossiblePMeshForCurrentVertexSplit = new List<List<PlanktonMesh>>();
                #endregion

                // 对于每个InnerNode来说，判断它的度是否大于4
                //for (int i = 0; i < PDeepCopy.Vertices.Count - outerNodeCount; i++)
                //{
                //    //int degree = PDeepCopy.Vertices.GetValence(i);
                //}

                /* 这里先只分割0点 */
                // 计算0点的度
                int degree = PDeepCopy.Vertices.GetValence(0);

                switch (degree)
                {
                    case 5:
                        #region 对于度为5的顶点

                        #region 获取5条半边中的相邻2条半边，以及它们的起点和终点
                        // 对于一个顶点来说的所有可能的一对半边的index（相邻的半边）
                        List<int[,]> allPossibleHalfedgeVertexIndexs = new List<int[,]>();

                        /* 这里先只分割0点 */
                        // 点0周围的halfedgeindex列表，顺时针
                        int[] halfEdgesIndexStartFromVertex = PDeepCopy.Vertices.GetHalfedges(0);

                        for (int i = 0; i < halfEdgesIndexStartFromVertex.Length; i++)
                        {

                            /* 不能用边的index来作为查找的根据，因为边的index是会变的
                             * 所以把找到这条边变为找到这个边的两个端点
                             */
                            int[,] possibleHalfedgeVertexIndexs = new int[2, 2] {
                                /* 初始化索引号为0的行 */
                                { PDeepCopy.Halfedges[halfEdgesIndexStartFromVertex[i]].StartVertex,
                                  PDeepCopy.Halfedges.EndVertex(halfEdgesIndexStartFromVertex[i])},
                                /* 初始化索引号为1的行 */
                                { PDeepCopy.Halfedges[halfEdgesIndexStartFromVertex[(i + 1) % halfEdgesIndexStartFromVertex.Length]].StartVertex,
                                  PDeepCopy.Halfedges.EndVertex(halfEdgesIndexStartFromVertex[(i + 1) % halfEdgesIndexStartFromVertex.Length])}};

                            allPossibleHalfedgeVertexIndexs.Add(possibleHalfedgeVertexIndexs);
                        }
                        #endregion

                        // debug Print
                        List<string> printFacesEdgeSplitedP = new List<string>();
                        List<string> printFacesHalfedgeSplitedP = new List<string>();
                        List<List<string>> printFacesEdgeResetP = new List<List<string>>();

                        #region 循环计算所有半边情况的所有分割可能性
                        // 对于每一种相邻2条半边的情况
                        for (int i = 0; i < allPossibleHalfedgeVertexIndexs.Count; i++)
                        {
                            // 每循环一次，都需要把List<string> printFacesEdgeResetP清理一次
                            printFacesEdgeResetP.Clear();

                            allPossiblePGraphLoL.Add(new List<List<int>>());
                            allPossiblePNodes.Add(new List<Node>());
                            allPossiblePMeshForCurrentVertexSplit.Add(new List<PlanktonMesh>());


                            
                            // 向存储当前顶点所有可能的分裂结果的列表中添加分裂的可能性
                            List<List<int>> newPGraphLoL;
                            List<Node> newPNodes;
                            int newVertexIndex;
                            edgeSplitedP = SplitEdgeIntoTwo(PDeepCopy,
                                                            allPossibleHalfedgeVertexIndexs[i][0, 0],
                                                            allPossibleHalfedgeVertexIndexs[i][0, 1],
                                                            2,
                                                            pGraphLoL,
                                                            pNodes,
                                                            out newPGraphLoL,
                                                            out newPNodes,
                                                            out newVertexIndex);

                            allPossiblePGraphLoL[i].AddRange(newPGraphLoL);
                            allPossiblePNodes[i].AddRange(newPNodes);

                            // debug printFacesEdgeSplitedP
                            printFacesEdgeSplitedP = UtilityFunctions.PrintFacesVertices(edgeSplitedP);
                            printFacesHalfedgeSplitedP = UtilityFunctions.PrintFacesHalfedges(edgeSplitedP);

                            // 得到这种相邻2条半边情况下的所有可能的分裂情况
                            edgeResetStartP = ResetHalfedgeStart(edgeSplitedP,
                                                                 allPossibleHalfedgeVertexIndexs[i][1, 0],
                                                                 allPossibleHalfedgeVertexIndexs[i][1, 1],
                                                                 0,
                                                                 newVertexIndex);

                            // debug printFacesEdgeResetP
                            for (int j = 0; j < edgeResetStartP.Count; j++)
                            {
                                printFacesEdgeResetP.Add(new List<string>());
                                printFacesEdgeResetP[j] = UtilityFunctions.PrintFacesVertices(edgeResetStartP[j]);
                            }

                            // 将得到的所有的可能的分裂情况，对应添加到代表对应的相邻2条半边情况的分支中
                            allPossiblePMeshForCurrentVertexSplit[i].AddRange(edgeResetStartP);

                            //int exceptionIndex;
                            //try
                            //{

                            //}
                            //catch (Exception)
                            //{
                            //    allPossiblePMeshForCurrentVertexSplit[i].AddRange(new List<PlanktonMesh>());
                            //    exceptionIndex = i;
                            //}
                        }
                        #endregion

                        #endregion
                        break;

                    default:
                        break;
                }


                // 获取要选择的序号
                DA.GetData<int>("IndexOfSplittedHalfedgeMesh", ref index);
                if (index >= allPossiblePMeshForCurrentVertexSplit.Count)
                {
                    index = index % allPossiblePMeshForCurrentVertexSplit.Count;
                }

                #region 电池结果输出
                #region 输出所选的HalfedgeMesh结果
                /* 注意这里第二个[]中，暂时填0 */
                PlanktonMesh selectedHalfedgeMesh = allPossiblePMeshForCurrentVertexSplit[index][0];
                DA.SetData("NewTriangleHalfedgeMesh", selectedHalfedgeMesh);
                #endregion
                #region 输出所选的HalfedgeMesh对应的NewGraphLoL
                List<List<int>> newGraphLoL = allPossiblePGraphLoL[index];
                DA.SetDataTree(1, UtilityFunctions.LoLToDataTree<int>(newGraphLoL));
                #endregion
                #region 输出所选的HalfedgeMesh对应的NewGraphNode
                List<Node> newGraphNode = allPossiblePNodes[index];
                DA.SetDataList("NewGraphNode", newGraphNode);
                #endregion
                #endregion

                #region 用于输出Debug数据

                #region HalfedgeMesh的顶点数据
                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(selectedHalfedgeMesh);
                DA.SetDataList("DebugVerticesOutput", printVertices);
                #endregion

                #region HalfedgeMesh的半边数据
                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(selectedHalfedgeMesh);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);
                #endregion

                #region HalfedgeMesh的每个面由哪些顶点构成
                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFacesVertices(selectedHalfedgeMesh);
                DA.SetDataList("DebugFacesOutput", printFaces);
                #endregion

                #region HalfedgeMesh的每个面由哪些半边构成
                List<string> printFacesHalfedge = new List<string>();
                printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(selectedHalfedgeMesh);
                DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);
                #endregion

                #endregion

                #region 可视化部分

                #region 设置用于可视化的TextDot
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();

                List<List<int>> selectedHalfedgeMeshPGraphLoL = allPossiblePGraphLoL[index];
                List<Node> selectedHalfedgeMeshPNodes = allPossiblePNodes[index];


                for (int i = 0; i < selectedHalfedgeMeshPNodes.Count; i++)
                {
                    if (selectedHalfedgeMeshPNodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, selectedHalfedgeMeshPNodes[i].NodeAttribute.NodeLabel), selectedHalfedgeMeshPNodes[i].NodeVertex);
                        InnerNodeTextDot.Add(textDot);
                    }
                    else
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, selectedHalfedgeMeshPNodes[i].NodeAttribute.NodeLabel), selectedHalfedgeMeshPNodes[i].NodeVertex);
                        OuterNodeTextDot.Add(textDot);
                    }
                }
                #endregion

                #region 转换为RhinoMesh，在DrawViewportMeshes绘制选中的mesh
                PRhinoMesh = RhinoSupport.ToRhinoMesh(selectedHalfedgeMesh);
                #endregion

                #region 在DrawViewportWires绘制mesh的edge（虚线）
                List<Line> selectedHalfedgeMeshPEdges = new List<Line>();
                for (int i = 0; i < PRhinoMesh.TopologyEdges.Count; i++)
                {
                    selectedHalfedgeMeshPEdges.Add(PRhinoMesh.TopologyEdges.EdgeLine(i));
                }

                PHalfedgeDottedCurves.Clear();

                double[] pattern = { 1.0 };
                for (int i = 0; i < selectedHalfedgeMeshPEdges.Count; i++)
                {
                    IEnumerable<Curve> segments = UtilityFunctions.ApplyDashPattern(selectedHalfedgeMeshPEdges[i].ToNurbsCurve(), pattern);
                    foreach (Curve segment in segments)
                    {
                        PHalfedgeDottedCurves.Add(segment);
                    }
                }
                #endregion

                #region 绘制表示图结构关系的实线
                GraphEdges.Clear();
                GraphEdges.AddRange(GraphEdgeLine(selectedHalfedgeMeshPGraphLoL, selectedHalfedgeMeshPNodes));
                #endregion

                #endregion
            }
        }

        /// <summary>
        /// 将选中的halfedgeIndex的半边进行分割（即这条半边的起点，分裂成了两个点）
        /// </summary>
        /// <param name="P"></param>
        /// <param name="halfedgeIndex"></param>
        /// <param name="pGraphLoL"></param>
        /// <param name="pNodes"></param>
        /// <param name="splitParts"></param>
        /// <returns></returns>
        public PlanktonMesh SplitEdgeIntoTwo(PlanktonMesh P, 
                                             int halfedgeStartVertexIndex, 
                                             int halfedgeEndVertexIndex, 
                                             int splitParts, 
                                             List<List<int>> pGraphLoL, 
                                             List<Node> pNodes, 
                                             out List<List<int>> newPGraphLoL, 
                                             out List<Node> newPNodes, 
                                             out int newVertexIndex)
        {
            #region 深拷贝
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);

            //List<List<int>> newPDeepCopyGraphLoL = new List<List<int>>();
            //pGraphLoL.ForEach(i => newPDeepCopyGraphLoL.Add(i));
            //List<Node> newPDeepCopyNodes = new List<Node>();
            //pNodes.ForEach(i => newPDeepCopyNodes.Add(i));

            //List<List<int>> newPDeepCopyGraphLoL = new List<List<int>>();
            //for (int i = 0; i < pGraphLoL.Count; i++)
            //{
            //    newPDeepCopyGraphLoL.Add(Clone<List<int>>(pGraphLoL[i]));
            //}

            //List<Node> newPDeepCopyNodes = new List<Node>();
            //for (int i = 0; i < pNodes.Count; i++)
            //{
            //    newPDeepCopyNodes.Add(Clone<Node>(pNodes[i]));
            //}

            List<List<int>> newPDeepCopyGraphLoL = new List<List<int>>();
            List<Node> newPDeepCopyNodes = new List<Node>();

            for (int i = 0; i < pGraphLoL.Count; i++)
            {
                newPDeepCopyGraphLoL.Add(new List<int>());
                newPDeepCopyGraphLoL[i].AddRange(pGraphLoL[i].Select(num => num));
            }

            for (int i = 0; i < pNodes.Count; i++)
            {
                newPDeepCopyNodes.Add(new Node(pNodes[i]));
            }

            #endregion


            #region 计算要分裂的边两侧的Face，准备好两个Face的顶点列表，新增midVertex作为分裂产生的新顶点

            // 找到halfedgeStartVertexIndex和halfedgeEndVertexIndex所对应的那条半边
            int halfedgeIndex = -1;
            foreach (int hfIndex in PDeepCopy.Vertices.GetHalfedges(halfedgeStartVertexIndex))
            {
                if (PDeepCopy.Halfedges.EndVertex(hfIndex) == halfedgeEndVertexIndex)
                {
                    halfedgeIndex = hfIndex;
                    break;
                }
            }
            if (halfedgeIndex == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ResetHalfedgeStart函数出错，找不到输入起点和终点所对应的半边。");
                newPGraphLoL = null;
                newPNodes = null;
                newVertexIndex = -1;
                return null;
            }

            int adjacentFaceIndex = PDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = PDeepCopy.Halfedges[PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            //// 包围AdjacentFace的边
            //List<int> hiAroundAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[adjacentFaceIndex].FirstHalfedge).ToList<int>();
            //// 包围PairAdjacentFace的边
            //List<int> hiAroundPairAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[pairAdjacentFaceIndex].FirstHalfedge).ToList<int>();

            // 包围AdjacentFace的顶点
            List<int> viAroundAdjacentFace = PDeepCopy.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();
            // 包围PairAdjcentFace的顶点
            List<int> viAroundPairAdjacentFace = PDeepCopy.Faces.GetFaceVertices(pairAdjacentFaceIndex).ToList<int>();

            //int startVertexIndex = PDeepCopy.Halfedges[halfedgeIndex].StartVertex;
            //int endVertexIndex = PDeepCopy.Halfedges.EndVertex(halfedgeIndex);
            int startVertexIndex = halfedgeStartVertexIndex;
            int endVertexIndex = halfedgeEndVertexIndex;

            Point3d startVertex = PDeepCopy.Vertices[startVertexIndex].ToPoint3d();
            Point3d endVertex = PDeepCopy.Vertices[endVertexIndex].ToPoint3d();
            Point3d midVertex = (startVertex + endVertex) / 2;

            // 向P中增加顶点
            PDeepCopy.Vertices.Add(midVertex);
            int midVertexIndex = PDeepCopy.Vertices.Count - 1;
            newVertexIndex = midVertexIndex;
            #endregion

            #region 为splitEdge的操作更新图结构和Node对象
            // 更新图结构
            for (int i = 0; i < newPDeepCopyGraphLoL.Count; i++)
            {
                if (i == startVertexIndex)
                {
                    newPDeepCopyGraphLoL[i].Add(midVertexIndex);
                }
                if (i == endVertexIndex)
                {
                    newPDeepCopyGraphLoL[i].Add(midVertexIndex);
                }
            }
            newPDeepCopyGraphLoL.Add((new int[2] { startVertexIndex, endVertexIndex }).ToList<int>());

            // 对startVertex的Node的相关属性进行修正，后续需要根据规则修改补充
            string startNodeLabel = newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeLabel;
            double startNodeArea = newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeArea;
            // double startNodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion;

            newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeLabel = newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeLabel + "0";
            newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeArea = newPDeepCopyNodes[startVertexIndex].NodeAttribute.NodeArea / splitParts;
            // nodes[startVertexIndex].NodeAttribute.NodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion / splitParts;

            // 对新生成的顶点，构造新的node
            NodeAttribute nodeAttribute = new NodeAttribute((startNodeLabel + (midVertexIndex - 8).ToString()),
                                                             startNodeArea / splitParts);
            nodeAttribute.ConnectivityTable = new int[] { startVertexIndex, endVertexIndex };
            nodeAttribute.AdjacencyTable = new int[] { };
            newPDeepCopyNodes.Add(new Node(midVertex, nodeAttribute, true));

            // 输出out参数
            newPGraphLoL = newPDeepCopyGraphLoL;
            newPNodes = newPDeepCopyNodes;
            #endregion

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

            #region 转移，更新，生成新的PlanktonMesh
            // 建立需要进行置换的面的index列表和VertexIndex列表，方便转移P的Face属性
            List<int> needChangeFaceIndexs = new List<int>();
            needChangeFaceIndexs.Add(adjacentFaceIndex);
            needChangeFaceIndexs.Add(pairAdjacentFaceIndex);
            List<List<int>> needChangeViAround = new List<List<int>>();
            needChangeViAround.Add(viAroundAdjacentFace);
            needChangeViAround.Add(viAroundPairAdjacentFace);

            // 转移P的Vertex属性
            List<PlanktonXYZ> pPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
            {
                pPlanktonVertex.Add(PDeepCopy.Vertices[i].ToXYZ());
            }

            // 转移P的Face属性，同时替换掉两个应该删除的面
            List<List<int>> pFaceVertexOrder = new List<List<int>>();

            for (int i = 0; i < PDeepCopy.Faces.Count; i++)
            {
                if (needChangeFaceIndexs.Contains(i))
                {
                    pFaceVertexOrder.Add(new List<int>());
                    pFaceVertexOrder[i].AddRange(needChangeViAround[needChangeFaceIndexs.IndexOf(i)]);
                }
                else
                {
                    pFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = PDeepCopy.Faces.GetFaceVertices(i);
                    pFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }

            // 用转移的P的Vertex属性和P的Face属性来构造新的PlanktonMesh newP
            PlanktonMesh newP = new PlanktonMesh();

            newP.Vertices.AddVertices(pPlanktonVertex);
            newP.Faces.AddFaces(pFaceVertexOrder);
            #endregion

            return newP;
        }

        /// <summary>
        /// 将一对半边（其中一条起点是vertexToSplitIndex）移动到一个新的顶点（targetVertexIndex）上，并且将这个顶点的度，补满到4
        /// </summary>
        /// <param name="P"></param>
        /// <param name="halfedgeIndex">要移动的一对半边中的一个，它的起点是要改变的</param>
        /// <param name="vertexToSplitIndex">要移动的一对半边中的一个，该半边的起点</param>
        /// <param name="targetVertexIndex">移动到哪一个顶点上，起点改为这个顶点</param>
        /// <returns></returns>
        public List<PlanktonMesh> ResetHalfedgeStart(PlanktonMesh P, int halfedgeStartVertexIndex, int halfedgeEndVertexIndex, int vertexToSplitIndex, int targetVertexIndex)
        {
            #region 深拷贝
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);
            #endregion

            #region 将一对半边移动到一个新的顶点上 

            // 找到halfedgeStartVertexIndex和halfedgeEndVertexIndex所对应的那条半边
            int halfedgeIndex = -1;
            foreach (int hfIndex in PDeepCopy.Vertices.GetHalfedges(halfedgeStartVertexIndex))
            {
                if (PDeepCopy.Halfedges.EndVertex(hfIndex) == halfedgeEndVertexIndex)
                {
                    halfedgeIndex = hfIndex;
                    break;
                }
            }
            if (halfedgeIndex == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ResetHalfedgeStart函数出错，找不到输入起点和终点所对应的半边。");
                return null;
            }

            int pairHalfedgeIndex = PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex);

            int adjacentFaceIndex = PDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = PDeepCopy.Halfedges[PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            #region 包围AdjacentFace的顶点
            List<int> viAroundAdjacentFace = new List<int>();
            int currentHalfedgeIndex = halfedgeIndex;
            do
            {
                viAroundAdjacentFace.Add(PDeepCopy.Halfedges[currentHalfedgeIndex].StartVertex);
                currentHalfedgeIndex = PDeepCopy.Halfedges[currentHalfedgeIndex].NextHalfedge;
            } while (currentHalfedgeIndex != halfedgeIndex);
            #endregion

            #region 包围PairAdjacentFace的顶点
            List<int> viAroundPairAdjacentFace = new List<int>();
            currentHalfedgeIndex = pairHalfedgeIndex;
            do
            {
                viAroundPairAdjacentFace.Add(PDeepCopy.Halfedges.EndVertex(currentHalfedgeIndex));
                currentHalfedgeIndex = PDeepCopy.Halfedges[currentHalfedgeIndex].NextHalfedge;
            } while (currentHalfedgeIndex != pairHalfedgeIndex);
            #endregion

            #region 修改包围AdjacentFace的顶点列表以及包围PairAdjacentFace的顶点列表
            // 用目标顶点属于哪个面，来判断该怎么修改viAroundAdjacentFace和viAroundPairAdjacentFace列表
            //bool targetVertexBelongsToAdjacent = viAroundAdjacentFace.Contains(targetVertexIndex);
            //bool targetVertexBelongsToPairAdjacent = viAroundPairAdjacentFace.Contains(targetVertexIndex);

            // 获取targetVertex在Face顶点列表的序号
            int indexInFace = viAroundAdjacentFace.IndexOf(targetVertexIndex);
            // 获取targetVertex在PairFace顶点列表的序号
            int indexInPairFace = viAroundPairAdjacentFace.IndexOf(targetVertexIndex);

            // 当目标顶点targetVertex同时不属于这两个Face时报错
            if (indexInFace == -1 && indexInPairFace == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的targetVertex " + vertexToSplitIndex.ToString() + " 不属于要改变的面");
                return null;
            }
            // 当目标顶点targetVertex属于PairFace时
            else if (indexInFace == -1 && indexInPairFace != -1)
            {
                // 当targetVertex在PairFace中时
                // 从PairFace顶点列表中删除第一个（即vertexToSplit）
                viAroundPairAdjacentFace.RemoveAt(0);
                // 把targetVertex按照它在PairFace顶点列表中的序号，放在Face顶点列表中
                // 需要判断是放在队尾，还是插入
                if (indexInPairFace >= viAroundAdjacentFace.Count)
                {
                    viAroundAdjacentFace.Add(targetVertexIndex);
                }
                else
                {
                    viAroundAdjacentFace.Insert(indexInPairFace, targetVertexIndex);
                }
            }
            // 当目标顶点targetVertex属于Face时
            else
            {
                // 当targetVertex在Face中时
                // 从Face顶点列表中删除第一个（即vertexToSplit）
                viAroundAdjacentFace.RemoveAt(0);
                // 把targetVertex按照它在Face顶点列表中的序号，放在PairFace顶点列表中
                // 需要判断是放在队尾，还是插入
                if (indexInFace >= viAroundPairAdjacentFace.Count)
                {
                    viAroundPairAdjacentFace.Add(targetVertexIndex);
                }
                else
                {
                    viAroundPairAdjacentFace.Insert(indexInFace, targetVertexIndex);
                }
            }

            //if (!targetVertexBelongsToAdjacent && !targetVertexBelongsToPairAdjacent)
            //{
            //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的targetVertex " + vertexToSplitIndex.ToString() + " 不属于要改变的面");
            //    return null;
            //}
            //// 当目标顶点targetVertex属于pairAdjacentFace时
            //else if (!targetVertexBelongsToAdjacent && targetVertexBelongsToPairAdjacent)
            //{
            //    //viAroundAdjacentFace.Insert(1, targetVertexIndex);
            //    //viAroundPairAdjacentFace.RemoveAt(0);
            //}
            //// 当目标顶点targetVertex属于AdjacentFace时
            //else
            //{
            //    //viAroundAdjacentFace.RemoveAt(0);
            //    //viAroundPairAdjacentFace.Insert(1, targetVertexIndex);
            //}
            #endregion

            #region 用变量whichFaceNeedSplitIndex记录哪一个面在下面将度补满的过程中，需要再进行split
            List<bool> whichFaceNeedSplit = new List<bool>();
            //whichFaceNeedSplit.Add(!targetVertexBelongsToAdjacent);
            //whichFaceNeedSplit.Add(!targetVertexBelongsToPairAdjacent);
            whichFaceNeedSplit.Add(indexInFace < 0);
            whichFaceNeedSplit.Add(indexInPairFace < 0);
            int whichFaceNeedSplitLaterIndex = -1;
            #endregion
            #endregion

            #region 构建resetHalfEdgeStartP，用来存储将一对半边移动后形成的新mesh,同时将未改动的PDeepCopy中的其他顶点和面的数据，转移过来

            #region 建立需要进行置换的面的index列表和VertexIndex列表，方便转移P的Face属性
            List<int> needChangeFaceIndexs = new List<int>();
            needChangeFaceIndexs.Add(adjacentFaceIndex);
            needChangeFaceIndexs.Add(pairAdjacentFaceIndex);
            List<List<int>> needChangeViAround = new List<List<int>>();
            needChangeViAround.Add(viAroundAdjacentFace);
            needChangeViAround.Add(viAroundPairAdjacentFace);
            #endregion

            #region 转移P的Vertex属性
            List<PlanktonXYZ> resetPPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
            {
                resetPPlanktonVertex.Add(PDeepCopy.Vertices[i].ToXYZ());
            }
            #endregion

            #region 转移P的Face属性，同时记录下，哪个面需要在下面新生成的顶点度补满的过程中，进行分裂
            List<List<int>> resetPFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < PDeepCopy.Faces.Count; i++)
            {
                if (needChangeFaceIndexs.Contains(i))
                {
                    resetPFaceVertexOrder.Add(new List<int>());
                    resetPFaceVertexOrder[i].AddRange(needChangeViAround[needChangeFaceIndexs.IndexOf(i)]);
                    // 如果这个面是被加过点的，那么记下它的i
                    if (whichFaceNeedSplit[needChangeFaceIndexs.IndexOf(i)])
                    {
                        whichFaceNeedSplitLaterIndex = i;
                    }
                }
                else
                {
                    resetPFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = PDeepCopy.Faces.GetFaceVertices(i);
                    resetPFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            #endregion

            #region 用转移的P的Vertex属性和P的Face属性来构造新的PlanktonMesh resetHalfEdgeStartP
            PlanktonMesh resetHalfEdgeStartP = new PlanktonMesh();
            resetHalfEdgeStartP.Vertices.AddVertices(resetPPlanktonVertex);
            resetHalfEdgeStartP.Faces.AddFaces(resetPFaceVertexOrder);
            #endregion

            // debug
            List<string> printFaces1 = new List<string>();
            printFaces1 = UtilityFunctions.PrintFacesVertices(resetHalfEdgeStartP);
            #endregion

            #region 为了将新得到的顶点的度补满，需要对一个面进行新的分割，形成新的一对半边
            List<List<int>> splitedFaceWithOriginIndexLoL = new List<List<int>>();
            List<List<int>> splitedFaceWithNewIndexLoL = new List<List<int>>();

            #region 用vertexToSplit属于哪个面，来判断该怎么修改经过半边移动后新的viAroundAdjacentFace和viAroundPairAdjacentFace列表
            bool vertexToSplitBelongsToAdjacent = viAroundAdjacentFace.Contains(vertexToSplitIndex);
            bool vertexToSplitBelongsToPairAdjacent = viAroundPairAdjacentFace.Contains(vertexToSplitIndex);

            #region 如果分裂顶点vertexToSplit，属于AdjacentFace时
            if (vertexToSplitBelongsToAdjacent && !vertexToSplitBelongsToPairAdjacent)
            {
                #region 判断VertexToSplit和TargetVertex谁在谁左边
                bool flag = RelationOfTwoVerticesInFaceVertexList(vertexToSplitIndex, targetVertexIndex, viAroundAdjacentFace);
                #endregion

                SplitFaceVertexList(flag, vertexToSplitIndex, targetVertexIndex, viAroundAdjacentFace, out splitedFaceWithOriginIndexLoL, out splitedFaceWithNewIndexLoL);
            }
            #endregion
            #region 如果分裂顶点vertexToSplit，属于PairFace时
            else
            {
                #region 判断VertexToSplit和TargetVertex谁在谁左边
                bool flag = RelationOfTwoVerticesInFaceVertexList(vertexToSplitIndex, targetVertexIndex, viAroundPairAdjacentFace);
                #endregion

                SplitFaceVertexList(flag, vertexToSplitIndex, targetVertexIndex, viAroundPairAdjacentFace, out splitedFaceWithOriginIndexLoL, out splitedFaceWithNewIndexLoL);
            }
            #endregion
            #endregion

            #region 转移resetP的Vertex属性
            List<PlanktonXYZ> rebuildPPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < resetHalfEdgeStartP.Vertices.Count; i++)
            {
                rebuildPPlanktonVertex.Add(resetHalfEdgeStartP.Vertices[i].ToXYZ());
            }
            #endregion

            #region 转移resetP的Face属性
            /* 对于每一种分割面的情况（即一对splitedFaceWithOriginIndexLoL[i]和splitedFaceWithNewIndexLoL[i]）作为List<List<List<int>>>的第一层List
             * 第二层List存储PlanktonMesh上每个Face的情况
             * 第三层List存储每个Face的顶点
             */
            List<List<List<int>>> rebuildPFaceVertexOrder = new List<List<List<int>>>();
            for (int i = 0; i < splitedFaceWithOriginIndexLoL.Count; i++)
            {
                rebuildPFaceVertexOrder.Add(new List<List<int>>());
                for (int j = 0; j < resetHalfEdgeStartP.Faces.Count; j++)
                {
                    // rebuildPFaceVertexOrder[i].Add(new List<int>());
                    if (j == whichFaceNeedSplitLaterIndex)
                    {
                        rebuildPFaceVertexOrder[i].Add(splitedFaceWithOriginIndexLoL[i]);
                    }
                    else
                    {
                        int[] faceVertexOrder = resetHalfEdgeStartP.Faces.GetFaceVertices(j);
                        rebuildPFaceVertexOrder[i].Add(faceVertexOrder.ToList<int>());
                    }
                }
                // 注意这里还要额外加上split后新生成的面
                rebuildPFaceVertexOrder[i].Add(splitedFaceWithNewIndexLoL[i]);
            }
            #endregion

            #region 用转移的resetP的Vertex属性和P的Face属性来构造新的PlanktonMesh rebulidNewHalfedgeP
            /* 对于每一种分割面的情况，构造一个新的PlanktonMesh
             * 每个mesh有自己的顶点，自己的面
             */
            List<PlanktonMesh> rebulidNewHalfedgeP = new List<PlanktonMesh>();
            for (int i = 0; i < rebuildPFaceVertexOrder.Count; i++)
            {
                rebulidNewHalfedgeP.Add(new PlanktonMesh());
                rebulidNewHalfedgeP[i].Vertices.AddVertices(rebuildPPlanktonVertex);
            }
            for (int i = 0; i < rebuildPFaceVertexOrder.Count; i++)
            {
                for (int j = 0; j < rebuildPFaceVertexOrder[i].Count; j++)
                {
                    rebulidNewHalfedgeP[i].Faces.AddFace(rebuildPFaceVertexOrder[i][j]);
                }
            }
            #endregion

            // debug
            List<List<string>> printFaces2 = new List<List<string>>();
            for (int i = 0; i < rebulidNewHalfedgeP.Count; i++)
            {
                printFaces2.Add(new List<string>());
                printFaces2[i] = UtilityFunctions.PrintFacesVertices(rebulidNewHalfedgeP[i]);
            }
            #endregion
            return rebulidNewHalfedgeP;
        }

        public bool RelationOfTwoVerticesInFaceVertexList(int vertexToSplitIndex, int targetVertexIndex, List<int> viAroundFaceList)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            int index0 = 0;
            int index1 = 0;
            bool flag = false;

            for (int i = 0; i < viAroundFaceList.Count; i++)
            {

                if (viAroundFaceList[i] == vertexToSplitIndex)
                {
                    index0 = i;
                }
                if (viAroundFaceList[i] == targetVertexIndex)
                {
                    index1 = i;
                }
            }

            if (index0 < index1)
            {
                counterClockwiseCount = index1 - index0;
                clockwiseCount = viAroundFaceList.Count - index1 + index0;
            }
            else
            {
                counterClockwiseCount = viAroundFaceList.Count - index0 + index1;
                clockwiseCount = index0 - index1;
            }

            // 如果VertexToSplit和TargetVertex逆时针数的间隔大于顺时针数的间隔，那么VertexToSplit在TargetVertex的左边
            if (counterClockwiseCount < clockwiseCount)
            {
                flag = true;
            }
            else
            {
                flag = false;
            }

            return flag;
        }

        public void SplitFaceVertexList(bool flag,
                                        int vertexToSplitIndex, 
                                        int targetVertexIndex, 
                                        List<int> faceVertexList, 
                                        out List<List<int>> splitedFaceWithOriginIndexLoL, 
                                        out List<List<int>> splitedFaceWithNewIndexLoL)
        {
            splitedFaceWithNewIndexLoL = new List<List<int>>();
            splitedFaceWithOriginIndexLoL = new List<List<int>>();

            if (flag)
            {
                #region usedVertexIndex用来存储目前构成过子面的顶点
                List<int> usedVertexIndexs = new List<int>();
                // 添加vertexToSplit
                usedVertexIndexs.Add(vertexToSplitIndex);
                // 添加targetVertex
                usedVertexIndexs.Add(targetVertexIndex);
                // 添加targetVertex的下一个顶点，同时为了防止数组越界，用取余的方式来做index
                // 添加这个点进入usedVertex的目的是，避免vertexToSplit产生新的连接（即构成了vertexToSplit，targetVertex，targerVertexNext的三角形），使得它的度超过了4
                usedVertexIndexs.Add(faceVertexList[(faceVertexList.IndexOf(targetVertexIndex) + 1) % faceVertexList.Count]);
                #endregion
                // for循环计数
                int iteration = 0;

                for (int i = (faceVertexList.IndexOf(targetVertexIndex) + 2) % faceVertexList.Count; i < faceVertexList.Count - 1; i++)
                {
                    #region 求差集，得到除了usedVertex外的其他的点
                    List<int> unusedVertexIndexs = new List<int>();
                    unusedVertexIndexs = faceVertexList.Except(usedVertexIndexs).ToList<int>();
                    #endregion

                    #region 从targetVertex的下下一个顶点开始构造新分割产生的三角形newFace
                    // 此时usedVertexIndexs中不包含当前的faceVertexList[i]
                    splitedFaceWithNewIndexLoL.Add(new List<int>());
                    splitedFaceWithNewIndexLoL[iteration].Add(vertexToSplitIndex);
                    splitedFaceWithNewIndexLoL[iteration].Add(targetVertexIndex);
                    splitedFaceWithNewIndexLoL[iteration].AddRange(unusedVertexIndexs);
                    #endregion

                    #region 构造完新的三角形后剩余的顶点，组成originFace
                    splitedFaceWithOriginIndexLoL.Add(new List<int>());
                    // 求差集即可
                    List<int> remainingVertex = faceVertexList.Except(splitedFaceWithNewIndexLoL[iteration]).ToList<int>();
                    splitedFaceWithOriginIndexLoL[iteration].AddRange(remainingVertex);
                    splitedFaceWithOriginIndexLoL[iteration].Add(targetVertexIndex);
                    // 注意这里是add了splitedFaceWithNewIndexLoL[iteration].Last()
                    splitedFaceWithOriginIndexLoL[iteration].Add(splitedFaceWithNewIndexLoL[iteration].Last());
                    #endregion

                    // 持续更新usedVertexIndex和iteration
                    usedVertexIndexs.Add(faceVertexList[i]);
                    iteration++;
                }
            }
            else
            {
                // 把列表反一下，就可以复用上面的计算过程，相当于正好把vertexToSplit跟targetVertex两个顺序反过来
                faceVertexList.Reverse();

                #region usedVertexIndex用来存储目前构成过子面的顶点
                List<int> usedVertexIndexs = new List<int>();
                // 添加vertexToSplit
                usedVertexIndexs.Add(vertexToSplitIndex);
                // 添加targetVertex
                usedVertexIndexs.Add(targetVertexIndex);
                // 添加targetVertex的下一个顶点，同时为了防止数组越界，用取余的方式来做index
                // 添加这个点进入usedVertex的目的是，避免vertexToSplit产生新的连接（即构成了vertexToSplit，targetVertex，targerVertexNext的三角形），使得它的度超过了4
                usedVertexIndexs.Add(faceVertexList[(faceVertexList.IndexOf(targetVertexIndex) + 1) % faceVertexList.Count]);
                #endregion
                // for循环计数
                int iteration = 0;

                for (int i = (faceVertexList.IndexOf(targetVertexIndex) + 2) % faceVertexList.Count; i < faceVertexList.Count - 1; i++)
                {
                    #region 求差集，得到除了usedVertex外的其他点
                    List<int> unusedVertexIndexs = new List<int>();
                    unusedVertexIndexs = faceVertexList.Except(usedVertexIndexs).ToList<int>();
                    #endregion

                    #region 从targetVertex的下下一个顶点开始构造新分割产生的三角形newFace
                    // 此时usedVertexIndexs中不包含当前的faceVertexList[i]
                    splitedFaceWithOriginIndexLoL.Add(new List<int>());
                    splitedFaceWithOriginIndexLoL[iteration].Add(vertexToSplitIndex);
                    splitedFaceWithOriginIndexLoL[iteration].Add(targetVertexIndex);
                    splitedFaceWithOriginIndexLoL[iteration].AddRange(unusedVertexIndexs);

                    splitedFaceWithOriginIndexLoL[iteration].Reverse();
                    #endregion

                    #region 构造完新的三角形后剩余的顶点，组成originFace
                    splitedFaceWithNewIndexLoL.Add(new List<int>());

                    // 求差集即可
                    List<int> remainingVertex = faceVertexList.Except(splitedFaceWithOriginIndexLoL[iteration]).ToList<int>();
                    splitedFaceWithNewIndexLoL[iteration].AddRange(remainingVertex);
                    splitedFaceWithNewIndexLoL[iteration].Add(targetVertexIndex);
                    // 注意这里是add了splitedFaceWithNewIndexLoL[iteration].First()
                    splitedFaceWithNewIndexLoL[iteration].Add(splitedFaceWithOriginIndexLoL[iteration].First());
                    // splitedFaceWithNewIndexLoL[iteration].Reverse();
                    #endregion

                    // 持续更新usedVertexIndex和iteration
                    usedVertexIndexs.Add(faceVertexList[i]);
                    iteration++;
                }

                faceVertexList.Reverse();

                //int iteration = 0;
                //List<int> usedVertexIndexs = new List<int>();
                //// 添加targetVertex的上一个顶点，同时为了防止数组越界，用取余的方式来做index
                //// usedVertexIndexs.Add(faceVertexList[(faceVertexList.IndexOf(targetVertexIndex) - 1 + faceVertexList.Count) % faceVertexList.Count]);
                //// 添加targetVertex
                //usedVertexIndexs.Add(targetVertexIndex);
                //// 添加vertexToSplit
                //usedVertexIndexs.Add(vertexToSplitIndex);


                //List<int> unusedVertexIndexs = new List<int>();

                //for (int i = (faceVertexList.IndexOf(vertexToSplitIndex) + 1) % faceVertexList.Count; i < faceVertexList.Count; i++)
                //{
                //    usedVertexIndexs.Add(faceVertexList[i]);

                //    // 求差集
                //    unusedVertexIndexs = (List<int>)faceVertexList.Except(usedVertexIndexs);

                //    #region 从vertexToSplit的下一个顶点开始构造新分割产生的三角形newFace

                //    #endregion
                //    // 此时usedVertexIndexs中不包含当前的faceVertexList[i]
                //    splitedFaceWithNewIndexLoL.Add(new List<int>());
                //    splitedFaceWithNewIndexLoL[iteration].AddRange(usedVertexIndexs);




                //    // 持续更新iteration
                //    iteration++;

                //}
            }
        }

        /// <summary>
        /// 绘制表示图结构关系的Line
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="graphVertices"></param>
        /// <returns></returns>
        public List<Line> GraphEdgeLine(List<List<int>> graph, List<Node> graphNodes)
        {
            List<Point3d> graphVertices = new List<Point3d>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                graphVertices.Add(graphNodes[i].NodeVertex);
            }

            List<Line> list = new List<Line>();
            for (int i = 0; i < graph.Count; i++)
            {
                for (int j = 0; j < graph[i].Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graph[i][j]]);
                    list.Add(item);
                }
            }
            return list;
        }

        /// <summary>
        /// 利用System.Runtime.Serialization序列化与反序列化完成引用对象的复制
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="RealObject"></param>
        /// <returns></returns>
        //public static T Clone<T>(T RealObject)
        //{
        //    using (Stream objectStream = new MemoryStream())
        //    {
        //        // 利用System.Runtime.Serialization序列化与反序列化完成引用对象的复制
        //        IFormatter formatter = new BinaryFormatter();
        //        formatter.Serialize(objectStream, RealObject);
        //        objectStream.Seek(0, SeekOrigin.Begin);
        //        return (T)formatter.Deserialize(objectStream);
        //    }
        //}

        /// <summary>
        /// 预览模式为WireFrame模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);

            for (int i = 0; i < PHalfedgeDottedCurves.Count; i++)
            {
                args.Display.DrawCurve(PHalfedgeDottedCurves[i], Color.DarkGreen, Thickness);
            }

            // 后画实线
            args.Display.DrawLines(GraphEdges, Color.BlueViolet, Thickness);

            for (int i = 0; i < InnerNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerNodeTextDot[i], Color.Black, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
            for (int i = 0; i < OuterNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(OuterNodeTextDot[i], Color.Gray, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        /// <summary>
        /// 预览模式为Shaded模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportMeshes(args);

            args.Display.DrawMeshShaded(PRhinoMesh, new Rhino.Display.DisplayMaterial(Color.White, 0));

            for (int i = 0; i < PHalfedgeDottedCurves.Count; i++)
            {
                args.Display.DrawCurve(PHalfedgeDottedCurves[i], Color.DarkGreen, Thickness);
            }

            // 后画实线
            args.Display.DrawLines(GraphEdges, Color.BlueViolet, Thickness);

            for (int i = 0; i < InnerNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerNodeTextDot[i], Color.Black, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
            for (int i = 0; i < OuterNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(OuterNodeTextDot[i], Color.Gray, Color.White, Color.White);
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
            get { return new Guid("8505b065-3ef0-4f90-bcb8-9cbf18874dca"); }
        }
    }
}