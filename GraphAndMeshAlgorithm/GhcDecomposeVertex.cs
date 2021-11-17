using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Attributes;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;
using Grasshopper.GUI.Canvas;

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

            Tree = NodeTree<GraphWithHM>.NewTree();

            ExceptionReports = new List<string>();
            CorrectReports = new List<string>();
            DefaultReports = new List<string>();

            ParentINode = null;
            CurrentINode = null;


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

        private ITree<GraphWithHM> Tree;

        private List<string> ExceptionReports;
        private List<string> CorrectReports;
        private List<string> DefaultReports;

        private INode<GraphWithHM> ParentINode;
        private INode<GraphWithHM> CurrentINode;

        private int CountOfAllLeafNodes = 0;
        private int CurrentLeafNodeIndex = 0;

        private int SelectedBranchDepth = 0;
        private int CurrentDepth = 0;

        private int RecursionDepth = -1;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphWithHM", "GHM", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.item);
            // pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);

            // pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddIntegerParameter("IndexOfGraphWithHM", "I", "某一种顶点分裂后的结果", GH_ParamAccess.item, 0);

            pManager.AddIntegerParameter("IndexInItsGenealogy", "IG", "", GH_ParamAccess.item, 0);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // pManager.AddGenericParameter("NewTriangleHalfedgeMesh", "NTHMesh", "新生成的三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddIntegerParameter("CountOfAllPossibleGraphWithHM", "C", "所有可能的结果的总数量", GH_ParamAccess.item);

            // pManager.AddGenericParameter("OriginGraphWithHM", "OGHM", "描述VolumeNode和BoundaryNode的所有连接关系的图结构(包括新分裂产生的点)", GH_ParamAccess.item);

            pManager.AddGenericParameter("SelectedGraphWithHM", "SGHM", "描述VolumeNode和BoundaryNode的所有连接关系的图结构(包括新分裂产生的点)", GH_ParamAccess.item);
            // pManager.AddGenericParameter("NewGraphNode", "NGN", "更新后的图结构中的节点", GH_ParamAccess.list);

            pManager.AddGenericParameter("SelectedAllPath", "SP", "", GH_ParamAccess.list);

            // pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugVertexConnection", "DebugVC", "Debug结果顶点连接关系", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);

            pManager.AddGenericParameter("SelectedTreeNodeString", "S", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("SelectedTreeNodeHistory", "SH", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Tree.Clear();
            #region 局部变量初始化

            GraphWithHM graphWithHM = new GraphWithHM();

            int indexOfLeafNodes = 0;
            int indexInItsGenealogy = 0;

            ExceptionReports.Clear();
            CorrectReports.Clear();
            DefaultReports.Clear();

            Thickness = 2;
            #endregion

            if (DA.GetData<GraphWithHM>("GraphWithHM", ref graphWithHM)
                && DA.GetData<int>("IndexOfGraphWithHM", ref indexOfLeafNodes))
            {
                DA.GetData<int>("IndexInItsGenealogy", ref indexInItsGenealogy);
                
                #region 深拷贝GraphWithHM

                GraphWithHM graphWithHMDeepCopy = new GraphWithHM(graphWithHM);
                graphWithHMDeepCopy.TreeNodeLabel = "TreeTop";

                #endregion


                /* 考虑用树形结构来处理和存储每一次分裂 */
                // 创建一个树
                // 对输入的GraphWithHM，形成树的根节点
                INode<GraphWithHM> top = Tree.AddChild(graphWithHMDeepCopy);

                // 新分裂产生的Vertex的Index集合
                List<int> viOfNewDecomposed = new List<int>();
                // 已经分裂过的Vertex的Index集合
                List<int> viHasBeenDecomposed = new List<int>();

                #region 原来能正确运行的代码
                ///* 这里先只分割0点 */
                //// 计算0点的度
                //int degree = PDeepCopy.Vertices.GetValence(0);

                //switch (degree)
                //{
                //    case 5:
                //        #region 对于度为5的顶点

                //        #region 获取5条半边中的相邻2条半边，以及它们的起点和终点
                //        /* 这里先只分割0点 */
                //        List<int[,]> allPossibleHalfedgeVertexIndexs = GetAllPossibleHStartAndEndIndexs(PDeepCopy, 0);
                //        #endregion

                //        #region 循环计算所有半边情况的所有分割可能性
                //        // 对于每一种相邻2条半边的情况
                //        for (int i = 0; i < allPossibleHalfedgeVertexIndexs.Count; i++)
                //        {
                //            //allPossiblePGraphLoL.Add(new List<List<int>>());
                //            //allPossiblePNodes.Add(new List<Node>());
                //            //allPossiblePMeshForCurrentVertexSplit.Add(new List<PlanktonMesh>());

                //            // 向存储当前顶点所有可能的分裂结果的列表中添加分裂的可能性
                //            List<List<int>> newPGraphLoL;
                //            List<Node> newPNodes;
                //            int newVertexIndex;
                //            edgeSplitedP = SplitEdgeIntoTwo(PDeepCopy,
                //                                            allPossibleHalfedgeVertexIndexs[i][0, 0],
                //                                            allPossibleHalfedgeVertexIndexs[i][0, 1],
                //                                            2,
                //                                            pGraphLoL,
                //                                            pNodes,
                //                                            out newPGraphLoL,
                //                                            out newPNodes,
                //                                            out newVertexIndex);

                //            //allPossiblePGraphLoL[i].AddRange(newPGraphLoL);
                //            //allPossiblePNodes[i].AddRange(newPNodes);

                //            #region debug printFacesEdgeSplitedP
                //            List<string> printFacesEdgeSplitedP = new List<string>();
                //            List<string> printFacesHalfedgeSplitedP = new List<string>();
                //            printFacesEdgeSplitedP = UtilityFunctions.PrintFacesVertices(edgeSplitedP);
                //            printFacesHalfedgeSplitedP = UtilityFunctions.PrintFacesHalfedges(edgeSplitedP);
                //            #endregion


                //            // 得到这种相邻2条半边情况下的所有可能的分裂情况
                //            edgeResetStartP = ResetHalfedgeStart(edgeSplitedP,
                //                                                 allPossibleHalfedgeVertexIndexs[i][1, 0],
                //                                                 allPossibleHalfedgeVertexIndexs[i][1, 1],
                //                                                 0,
                //                                                 newVertexIndex);

                //            #region debug printFacesEdgeResetP
                //            List<List<string>> printFacesEdgeResetP = new List<List<string>>();
                //            // 每循环一次，都需要把List<string> printFacesEdgeResetP清理一次
                //            // printFacesEdgeResetP.Clear();
                //            for (int j = 0; j < edgeResetStartP.Count; j++)
                //            {
                //                printFacesEdgeResetP.Add(new List<string>());
                //                printFacesEdgeResetP[j] = UtilityFunctions.PrintFacesVertices(edgeResetStartP[j]);
                //            }
                //            #endregion

                //            //// 将得到的所有的可能的分裂情况，对应添加到代表对应的相邻2条半边情况的分支中
                //            //allPossiblePMeshForCurrentVertexSplit[i].AddRange(edgeResetStartP);

                //            #region 构造对于每一种半边选择可能性所产生的DecomposedHM
                //            List<DecomposedHM> currentIterDHM = new List<DecomposedHM>();
                //            for (int j = 0; j < edgeResetStartP.Count; j++)
                //            {
                //                currentIterDHM.Add(new DecomposedHM(edgeResetStartP[i], newPGraphLoL, newPNodes));
                //            }
                //            #endregion
                //        }
                //        #endregion

                //        #endregion
                //        break;

                //    default:
                //        break;
                //}


                //// 获取要选择的序号
                //DA.GetData<int>("IndexOfSplittedHalfedgeMesh", ref index);
                //if (index >= allPossiblePMeshForCurrentVertexSplit.Count)
                //{
                //    index = index % allPossiblePMeshForCurrentVertexSplit.Count;
                //}
                #endregion

                #region 测试递归

                List<int> currentInnerNodesNeedToSplit = new List<int>();
                for (int i = 0; i < graphWithHMDeepCopy.GraphNodes.Count; i++)
                {
                    if (graphWithHMDeepCopy.GraphNodes[i].IsInner)
                    {
                        currentInnerNodesNeedToSplit.Add(i);
                    }
                }

                // GenerateDecomposedHMs(top, viOfNewDecomposed, viHasBeenDecomposed);
                Decompose(top, currentInnerNodesNeedToSplit);
                // 结束递归运算后，深度回归-1
                RecursionDepth = -1;

                List<INode<GraphWithHM>> allLeafNodes = new List<INode<GraphWithHM>>();
                foreach (INode<GraphWithHM> treeNodeGraphWithHM in Tree.AllChildren.Nodes)
                {
                    if (treeNodeGraphWithHM.HasChild == false)
                    {
                        //if (treeNodeGraphWithHM.Data.InnerNodeCount != 0)
                        //{
                        //    allLeafNodes.Add(treeNodeGraphWithHM.Data);
                        //}
                        allLeafNodes.Add(treeNodeGraphWithHM);
                    }
                }


                // TreeDebug
                INode<GraphWithHM> treeTop = Tree.Root.Child;

                List<INode<GraphWithHM>> allChildren = new List<INode<GraphWithHM>>();
                foreach (var item in Tree.AllChildren.Nodes)
                {
                    allChildren.Add(item);
                }


                INode<GraphWithHM> selectedINode = allLeafNodes[indexOfLeafNodes];

                // -1 是为了减去Root的影响
                SelectedBranchDepth = selectedINode.Depth - 1;

                List<INode<GraphWithHM>> selectedGenealogy = new List<INode<GraphWithHM>>();
                INode<GraphWithHM> current = selectedINode;
                while (current.HasParent)
                {
                    selectedGenealogy.Add(current);
                    current = current.Parent;
                }

                selectedGenealogy.Reverse();

                GraphWithHM seletedGHM = selectedGenealogy[indexInItsGenealogy].Data;
                DA.SetData("SelectedTreeNodeString", seletedGHM.TreeNodeLabel);

                GraphWithHM embededGraphWithHM = new GraphWithHM();
                if (seletedGHM.VolumeContainsWhichInnerNode != null)
                {
                    embededGraphWithHM = GraphWithHM.ReTutteEmbeding(seletedGHM);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "当前节点为空叶子节点，现在显示的是其父节点");
                    embededGraphWithHM = GraphWithHM.ReTutteEmbeding(selectedGenealogy[indexInItsGenealogy].Parent.Data);
                }


                #endregion

                


                List<string> selectedTreeNodeHistory = new List<string>();
                foreach (var keyValuePair in embededGraphWithHM.VolumeContainsWhichInnerNode)
                {
                    selectedTreeNodeHistory.Add(string.Format("{0}:{1}", keyValuePair.Key, string.Join(",", keyValuePair.Value)));
                }
                // string str = string.Join(";", selectedTreeNodeHistory);
                DA.SetDataList("SelectedTreeNodeHistory", selectedTreeNodeHistory);

                #region 电池结果输出

                // -1 是消除Root节点的影响
                CountOfAllLeafNodes = allLeafNodes.Count - 1;
                CurrentLeafNodeIndex = indexOfLeafNodes;

                DA.SetData("CountOfAllPossibleGraphWithHM", allLeafNodes.Count - 1);
                DA.SetData("SelectedGraphWithHM", embededGraphWithHM);
                #endregion

                #region 用于输出Debug数据

                //#region HalfedgeMesh的顶点数据
                //List<string> printVertices = new List<string>();
                //printVertices = UtilityFunctions.PrintVertices(embededGraphWithHM.PlanktonMesh);
                //DA.SetDataList("DebugVerticesOutput", printVertices);
                //#endregion

                //#region HalfedgeMesh的顶点连接数据
                //List<string> printVertexConnection = new List<string>();
                //printVertexConnection = UtilityFunctions.PrintVertexConnection(embededGraphWithHM.PlanktonMesh);
                //DA.SetDataList("DebugVertexConnection", printVertexConnection);
                //#endregion


                //#region HalfedgeMesh的半边数据
                //List<string> printHalfedges = new List<string>();
                //printHalfedges = UtilityFunctions.PrintHalfedges(embededGraphWithHM.PlanktonMesh);
                //DA.SetDataList("DebugHalfedgesOutput", printHalfedges);
                //#endregion

                //#region HalfedgeMesh的每个面由哪些顶点构成
                //List<string> printFaces = new List<string>();
                //printFaces = UtilityFunctions.PrintFacesVertices(embededGraphWithHM.PlanktonMesh);
                //DA.SetDataList("DebugFacesOutput", printFaces);
                //#endregion

                //#region HalfedgeMesh的每个面由哪些半边构成
                //List<string> printFacesHalfedge = new List<string>();
                //printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(embededGraphWithHM.PlanktonMesh);
                //DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);
                //#endregion

                #endregion

                #region 可视化部分

                #region 设置用于可视化的TextDot
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();

                for (int i = 0; i < embededGraphWithHM.GraphNodes.Count; i++)
                {
                    if (embededGraphWithHM.GraphNodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, embededGraphWithHM.GraphNodes[i].NodeAttribute.NodeLabel), embededGraphWithHM.GraphNodes[i].NodeVertex);
                        InnerNodeTextDot.Add(textDot);
                    }
                    else
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, embededGraphWithHM.GraphNodes[i].NodeAttribute.NodeLabel), embededGraphWithHM.GraphNodes[i].NodeVertex);
                        OuterNodeTextDot.Add(textDot);
                    }
                }
                #endregion

                #region 转换为RhinoMesh，在DrawViewportMeshes绘制选中的mesh
                // PRhinoMesh = RhinoSupport.ToRhinoMesh(embededGraphWithHM.PlanktonMesh);

                PlanktonMesh embededPlanktonMesh = embededGraphWithHM.PlanktonMesh;

                List<Line> edges = new List<Line>();
                List<int> usedHFIndex = new List<int>();
                for (int i = 0; i < embededPlanktonMesh.Halfedges.Count; i++)
                {
                    if (!usedHFIndex.Contains(i))
                    {
                        Point3d start = embededPlanktonMesh.Vertices[embededPlanktonMesh.Halfedges[i].StartVertex].ToPoint3d();
                        Point3d end = embededPlanktonMesh.Vertices[embededPlanktonMesh.Halfedges.EndVertex(i)].ToPoint3d();
                        edges.Add(new Line(start, end));
                        int pairHFIndex = embededPlanktonMesh.Halfedges.GetPairHalfedge(i);
                        usedHFIndex.Add(i);
                        usedHFIndex.Add(pairHFIndex);
                    }
                }
                #endregion

                #region 在DrawViewportWires绘制mesh的edge（虚线）

                PHalfedgeDottedCurves.Clear();

                double[] pattern = { 1.0 };
                for (int i = 0; i < edges.Count; i++)
                {
                    IEnumerable<Curve> segments = UtilityFunctions.ApplyDashPattern(edges[i].ToNurbsCurve(), pattern);
                    foreach (Curve segment in segments)
                    {
                        PHalfedgeDottedCurves.Add(segment);
                    }
                }
                #endregion

                #region 绘制表示图结构关系的实线
                GraphEdges.Clear();
                GraphEdges.AddRange(UtilityFunctions.GraphEdgeLine(embededGraphWithHM.GraphTables, embededGraphWithHM.GraphNodes));
                #endregion

                #endregion
            }
        }

        public void Decompose(INode<GraphWithHM> INodeToSplit, List<int> verticesToSplit)
        {
            //List<int> allInnerNodeNeedToSplit = new List<int>();
            //List<int> allInnerNodeNeedToSplitDegrees = new List<int>();
            //for (int i = 0; i < INodeToSplit.Data.GraphNodes.Count; i++)
            //{
            //    if (INodeToSplit.Data.GraphNodes[i].IsInner)
            //    {
            //        if (!INodeToSplit.Data.GraphNodes[i].IsNewDecomposed)
            //        {
            //            int degree = INodeToSplit.Data.PlanktonMesh.Vertices.GetValence(i);
            //            if (degree != 4)
            //            {
            //                allInnerNodeNeedToSplit.Add(i);
            //                allInnerNodeNeedToSplitDegrees.Add(degree);
            //            }
            //        }
            //    }
            //}

            // 开始执行当前层级的递归，递归深度++
            RecursionDepth++;


            int currentVertexToSplit = -1;
            // 当递归深度大于数组的长度时，直接结束执行当前递归，并递归深度--
            if (RecursionDepth < verticesToSplit.Count)
            {
                currentVertexToSplit = verticesToSplit[RecursionDepth];
            }
            else
            {
                // 结束执行当前递归，递归深度--
                RecursionDepth--;
                return;
            }

            GraphWithHM gHMDeepCopy = new GraphWithHM(INodeToSplit.Data);

            ParentINode = INodeToSplit;

            // int currentLevel = allInnerNodeNeedToSplit[i];

            int degree = gHMDeepCopy.PlanktonMesh.Vertices.GetValence(currentVertexToSplit);

            switch (degree)
            {
                case 5:
                    #region 获取5条半边中的相邻2条半边，以及它们的起点和终点
                    List<int[,]> allPossibleHalfedgeVertexIndexs = GetAllPossibleHStartAndEndIndexs(gHMDeepCopy.PlanktonMesh, currentVertexToSplit);
                    #endregion

                    #region 循环计算所有半边情况的所有分割可能性
                    // 对于每一种相邻2条半边的情况

                    for (int i = 0; i < allPossibleHalfedgeVertexIndexs.Count; i++)
                    {

                        // 向存储当前顶点所有可能的分裂结果的列表中添加分裂的可能性
                        // List<List<int>> newPGraphLoL;
                        List<GraphNode> newPNodes;
                        int newVertexIndex;
                        PlanktonMesh edgeSplitedP = SplitEdgeIntoTwo(gHMDeepCopy.PlanktonMesh,
                                                        allPossibleHalfedgeVertexIndexs[i],
                                                        2,
                                                        gHMDeepCopy.GraphNodes,
                                                        out newPNodes,
                                                        out newVertexIndex);

                        #region 当SplitEdgeIntoTwo函数无法正确分裂时，抛出异常
                        if (newVertexIndex == -1)
                        {
                            //string exceptionReport = INodeToSplit.Data.TreeNodeLabel + "->" + "VertexIndex:" +
                            //                         allInnerNodeNeedToSplit[i].ToString() +
                            //                         "-" + j.ToString() + " " +
                            //                         "cant split edge into two. Halfedge StartVertex is " +
                            //                         allPossibleHalfedgeVertexIndexs[j][0, 0].ToString() +
                            //                         "Halfedge EndVertex is " +
                            //                         allPossibleHalfedgeVertexIndexs[j][0, 1].ToString();
                            //ExceptionReports.Add(exceptionReport);

                            // continue继续下一次对allPossibleHalfedgeVertexIndexs的循环，本次的循环没有产生新的GraphWithHM（新的叶子节点）
                            // 为了debug，这里添加一个空的叶子节点
                            GraphWithHM newChildGraphWithHM = new GraphWithHM();
                            newChildGraphWithHM.TreeNodeLabel = "分裂了第" + currentVertexToSplit.ToString() + "个Vertex" +
                                                                "，所选的两个半边的起点终点是" +
                                                                allPossibleHalfedgeVertexIndexs[i][0, 0].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][0, 1].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][1, 0].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][1, 1].ToString() +
                                                                "，该子节点为空";
                            INode<GraphWithHM> childDHMToSplit = INodeToSplit.AddChild(newChildGraphWithHM);
                            CurrentINode = childDHMToSplit;
                            // 还要继续分裂剩余需要分裂的点

                            continue;
                        }
                        #endregion

                        List<int> newInnerNodeIndexs = new List<int>();
                        for (int k = 0; k < newPNodes.Count; k++)
                        {
                            if (newPNodes[k].IsInner)
                            {
                                newInnerNodeIndexs.Add(k);
                            }
                        }

                        // 得到这种相邻2条半边情况下的所有可能的分裂情况
                        List<int> newEndVertexAfterResetHalfedgeList;
                        List<PlanktonMesh> edgeResetStartP = ResetHalfedgeStart(edgeSplitedP,
                                                             allPossibleHalfedgeVertexIndexs[i],
                                                             currentVertexToSplit,
                                                             newVertexIndex,
                                                             newInnerNodeIndexs,
                                                             out newEndVertexAfterResetHalfedgeList);

                        #region 当ResetHalfedgeStart函数无法正确产生结果时，直接执行下一次的循环
                        // 如果新生成的半边的端点，是innerNode（包括初始的和新分裂产生的），那么返回null，然后结束这次循环（即不产生新的叶子节点）
                        if (edgeResetStartP == null)
                        {
                            // continue继续下一次对allPossibleHalfedgeVertexIndexs的循环，本次的循环没有产生新的GraphWithHM（新的叶子节点）
                            // 为了debug，这里添加一个空的叶子节点
                            GraphWithHM newChildGraphWithHM = new GraphWithHM();
                            newChildGraphWithHM.TreeNodeLabel = "分裂了第" + currentVertexToSplit.ToString() + "个Vertex" +
                                                                "，所选的两个半边的起点终点是" +
                                                                allPossibleHalfedgeVertexIndexs[i][0, 0].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][0, 1].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][1, 0].ToString() + ";" +
                                                                allPossibleHalfedgeVertexIndexs[i][1, 1].ToString() +
                                                                "，该子节点为空";
                            INode<GraphWithHM> childDHMToSplit = INodeToSplit.AddChild(newChildGraphWithHM);
                            CurrentINode = childDHMToSplit;
                            continue;
                        }
                        #endregion

                        List<int> newValue = new int[] { currentVertexToSplit, newVertexIndex }.ToList();
                        Dictionary<int, List<int>> volumeContainsWhichInnerNode = INodeToSplit.Data.VolumeContainsWhichInnerNode;
                        if (volumeContainsWhichInnerNode.ContainsKey(currentVertexToSplit))
                        {
                            List<int> subtraction = volumeContainsWhichInnerNode[currentVertexToSplit].Except(newValue).ToList();
                            if (subtraction.Count != 0)
                            {
                                volumeContainsWhichInnerNode[currentVertexToSplit].Add(newVertexIndex);
                            }
                        }
                        else
                        {
                            volumeContainsWhichInnerNode.Add(currentVertexToSplit, new int[] { currentVertexToSplit, newVertexIndex }.ToList());
                        }

                        #region 构造对于每一种半边选择可能性所产生的DecomposedHM，并且将它作为INode<DecomposedHM> dHMToSplit的叶子节点，并且递归分裂叶子节点
                        for (int j = 0; j < edgeResetStartP.Count; j++)
                        {
                            List<List<int>> newGraphTables = RenewGraphTable(gHMDeepCopy.GraphTables, allPossibleHalfedgeVertexIndexs[i], newVertexIndex, newEndVertexAfterResetHalfedgeList[j]);

                            string label = "分裂了第" + currentVertexToSplit.ToString() + "个Vertex" +
                                           "，所选的两个半边的起点终点是" +
                                           allPossibleHalfedgeVertexIndexs[i][0, 0].ToString() + ";" +
                                           allPossibleHalfedgeVertexIndexs[i][0, 1].ToString() + ";" +
                                           allPossibleHalfedgeVertexIndexs[i][1, 0].ToString() + ";" +
                                           allPossibleHalfedgeVertexIndexs[i][1, 1].ToString() +
                                           "，并且是第" + j.ToString() + "分裂可能性";



                            GraphWithHM currentDHM = new GraphWithHM(edgeResetStartP[j], newPNodes, newGraphTables, label, volumeContainsWhichInnerNode);

                            INode<GraphWithHM> childDHMToSplit = INodeToSplit.AddChild(currentDHM);
                            CurrentINode = childDHMToSplit;
                            //if (!viOfNewDecomposed.Contains(newVertexIndex))
                            //{
                            //    viOfNewDecomposed.Add(newVertexIndex);
                            //}
                            //if (!viHasBeenDecomposed.Contains(innerNodeToSplitIndexList[i]))
                            //{
                            //    viHasBeenDecomposed.Add(innerNodeToSplitIndexList[i]);
                            //}

                            // 递归分裂叶子节点
                            Decompose(childDHMToSplit, verticesToSplit);
                            // GenerateDecomposedHMs(childDHMToSplit, viOfNewDecomposed, viHasBeenDecomposed);
                        }
                        #endregion


                    }


                    #endregion

                    break;

                default:
                    break;
            }

            // 结束执行当前递归，递归深度--
            RecursionDepth--;
        }

        ///// <summary>
        ///// 对顶点进行递归分裂
        ///// </summary>
        ///// <param name="gHMToSplit"></param>
        ///// <param name="viOfNewDecomposed"></param>
        ///// <param name="viHasBeenDecomposed"></param>
        //public void GenerateDecomposedHMs(INode<GraphWithHM> gHMToSplit, List<int> viOfNewDecomposed, List<int> viHasBeenDecomposed)
        //{
        //    // 深拷贝DecomposedHM
        //    GraphWithHM gHMDeepCopy = new GraphWithHM(gHMToSplit.Data);
        //    ParentINode = gHMDeepCopy;

        //    List<int> innerNodeToSplitIndexList = new List<int>();
        //    for (int i = 0; i < gHMDeepCopy.PlanktonMesh.Vertices.Count; i++)
        //    {
        //        if (gHMDeepCopy.GraphNodes[i].IsInner)
        //        {
        //            innerNodeToSplitIndexList.Add(i);
        //        }
        //    }

        //    innerNodeToSplitIndexList = innerNodeToSplitIndexList.Except<int>(viOfNewDecomposed).ToList();
        //    innerNodeToSplitIndexList = innerNodeToSplitIndexList.Except<int>(viHasBeenDecomposed).ToList();

        //    if (innerNodeToSplitIndexList.Count == 0)
        //    {
                
        //        return;
        //    }
        //    else
        //    {
        //        for (int i = 0; i < innerNodeToSplitIndexList.Count; i++)
        //        {
        //            int degree = gHMDeepCopy.PlanktonMesh.Vertices.GetValence(innerNodeToSplitIndexList[i]);

                    

        //            switch (degree)
        //            {
        //                case 5:
        //                    #region 对于度为5的顶点

        //                    #region 获取5条半边中的相邻2条半边，以及它们的起点和终点
        //                    List<int[,]> allPossibleHalfedgeVertexIndexs = GetAllPossibleHStartAndEndIndexs(gHMDeepCopy.PlanktonMesh, innerNodeToSplitIndexList[i]);
        //                    #endregion

        //                    #region 循环计算所有半边情况的所有分割可能性
        //                    // 对于每一种相邻2条半边的情况
        //                    for (int j = 0; j < allPossibleHalfedgeVertexIndexs.Count; j++)
        //                    {

        //                        // 向存储当前顶点所有可能的分裂结果的列表中添加分裂的可能性
        //                        List<GraphNode> newPNodes;
        //                        int newVertexIndex;
        //                        PlanktonMesh edgeSplitedP = SplitEdgeIntoTwo(gHMDeepCopy.PlanktonMesh,
        //                                                        allPossibleHalfedgeVertexIndexs[j],
        //                                                        2,
        //                                                        gHMDeepCopy.GraphNodes,
        //                                                        out newPNodes,
        //                                                        out newVertexIndex);

        //                        #region 当SplitEdgeIntoTwo函数无法正确分裂时，抛出异常
        //                        if (newVertexIndex == -1)
        //                        {
        //                            string exceptionReport = gHMToSplit.Data.TreeNodeLabel + "->" + "VertexIndex:" +
        //                                                     innerNodeToSplitIndexList[i].ToString() +
        //                                                     "-" + j.ToString() + " " +
        //                                                     "cant split edge into two. Halfedge StartVertex is " +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 0].ToString() +
        //                                                     "Halfedge EndVertex is " +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 1].ToString();
        //                            ExceptionReports.Add(exceptionReport);

        //                            // continue继续下一次对allPossibleHalfedgeVertexIndexs的循环，本次的循环没有产生新的GraphWithHM（新的叶子节点）
        //                            // 为了debug，这里添加一个空的叶子节点
        //                            GraphWithHM newChildGraphWithHM = new GraphWithHM();
        //                            newChildGraphWithHM.TreeNodeLabel = "分裂了第" + innerNodeToSplitIndexList[i].ToString() + "个Vertex" +
        //                                                     "，所选的两个半边的起点终点是" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 0].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 1].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][1, 0].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][1, 1].ToString() +
        //                                                     "，该子节点为空";
        //                            INode<GraphWithHM> childDHMToSplit = gHMToSplit.AddChild(newChildGraphWithHM);
        //                            continue;
        //                        }
        //                        #endregion

        //                        // 得到这种相邻2条半边情况下的所有可能的分裂情况
        //                        List<int> newEndVertexAfterResetHalfedgeList;
        //                        List<PlanktonMesh> edgeResetStartP = ResetHalfedgeStart(edgeSplitedP,
        //                                                             allPossibleHalfedgeVertexIndexs[j],
        //                                                             innerNodeToSplitIndexList[i],
        //                                                             newVertexIndex,
        //                                                             innerNodeToSplitIndexList,
        //                                                             out newEndVertexAfterResetHalfedgeList);

        //                        #region 当ResetHalfedgeStart函数无法正确产生结果时，直接执行下一次的循环
        //                        // 如果新生成的半边的端点，是innerNode（包括初始的和新分裂产生的），那么返回null，然后结束这次循环（即不产生新的叶子节点）
        //                        if (edgeResetStartP == null)
        //                        {
        //                            // continue继续下一次对allPossibleHalfedgeVertexIndexs的循环，本次的循环没有产生新的GraphWithHM（新的叶子节点）
        //                            // 为了debug，这里添加一个空的叶子节点
        //                            GraphWithHM newChildGraphWithHM = new GraphWithHM();
        //                            newChildGraphWithHM.TreeNodeLabel = "分裂了第" + innerNodeToSplitIndexList[i].ToString() + "个Vertex" +
        //                                                     "，所选的两个半边的起点终点是" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 0].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][0, 1].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][1, 0].ToString() + ";" +
        //                                                     allPossibleHalfedgeVertexIndexs[j][1, 1].ToString() +
        //                                                     "，该子节点为空";
        //                            INode<GraphWithHM> childDHMToSplit = gHMToSplit.AddChild(newChildGraphWithHM);
        //                            continue;
        //                        }
        //                        #endregion

        //                        List<int> newValue = new int[] { innerNodeToSplitIndexList[i], newVertexIndex }.ToList();

        //                        Dictionary<int, List<int>> volumeContainsWhichInnerNode = new Dictionary<int, List<int>>();
        //                        foreach (KeyValuePair<int, List<int>> pair in ParentINode.VolumeContainsWhichInnerNode)
        //                        {
        //                            volumeContainsWhichInnerNode.Add(pair.Key, pair.Value);
        //                        }

        //                        if (volumeContainsWhichInnerNode.ContainsKey(innerNodeToSplitIndexList[i]))
        //                        {
        //                            List<int> subtraction = volumeContainsWhichInnerNode[innerNodeToSplitIndexList[i]].Except(newValue).ToList();
        //                            if (subtraction .Count != 0 )
        //                            {
        //                                volumeContainsWhichInnerNode[innerNodeToSplitIndexList[i]].Add(newVertexIndex);
        //                            }
        //                        }
        //                        else
        //                        {
        //                            volumeContainsWhichInnerNode.Add(innerNodeToSplitIndexList[i], new int[] { innerNodeToSplitIndexList[i], newVertexIndex }.ToList());
        //                        }

        //                        #region 构造对于每一种半边选择可能性所产生的DecomposedHM，并且将它作为INode<DecomposedHM> dHMToSplit的叶子节点，并且递归分裂叶子节点
        //                        for (int k = 0; k < edgeResetStartP.Count; k++)
        //                        {
        //                            List<List<int>> newGraphTables = RenewGraphTable(gHMDeepCopy.GraphTables, allPossibleHalfedgeVertexIndexs[j], newVertexIndex, newEndVertexAfterResetHalfedgeList[k]);

        //                            string label = "分裂了第" + innerNodeToSplitIndexList[i].ToString() + "个Vertex" +
        //                                           "，所选的两个半边的起点终点是" +
        //                                           allPossibleHalfedgeVertexIndexs[j][0, 0].ToString() + ";" +
        //                                           allPossibleHalfedgeVertexIndexs[j][0, 1].ToString() + ";" +
        //                                           allPossibleHalfedgeVertexIndexs[j][1, 0].ToString() + ";" +
        //                                           allPossibleHalfedgeVertexIndexs[j][1, 1].ToString() +
        //                                           "，并且是第" + k.ToString() + "分裂可能性";

        //                            GraphWithHM currentDHM = new GraphWithHM(edgeResetStartP[k], newPNodes, newGraphTables, label, volumeContainsWhichInnerNode);
        //                            CurrentINode = currentDHM;

        //                            INode<GraphWithHM> childDHMToSplit = gHMToSplit.AddChild(currentDHM);

        //                            if (!viOfNewDecomposed.Contains(newVertexIndex))
        //                            {
        //                                viOfNewDecomposed.Add(newVertexIndex);
        //                            }
        //                            if (!viHasBeenDecomposed.Contains(innerNodeToSplitIndexList[i]))
        //                            {
        //                                viHasBeenDecomposed.Add(innerNodeToSplitIndexList[i]);
        //                            }
                                    
        //                            // 递归分裂叶子节点
        //                            GenerateDecomposedHMs(childDHMToSplit, viOfNewDecomposed, viHasBeenDecomposed);
        //                        }
        //                        #endregion
        //                    }
        //                    #endregion

        //                    #endregion
        //                    break;
        //                case 4:
        //                    // 记录度为4的正确的情况
        //                    string correctReport = gHMDeepCopy.TreeNodeLabel + " -> " + "VertexIndex: " +
        //                                           innerNodeToSplitIndexList[i].ToString() + " " +
        //                                           "degree is " + degree.ToString();
        //                    CorrectReports.Add(correctReport);
        //                    break;

        //                default:
        //                    // 暂时先记录一下度为其他值的情况
        //                    string defaultReport = gHMDeepCopy.TreeNodeLabel + " -> " + "VertexIndex: " +
        //                                           innerNodeToSplitIndexList[i].ToString() + " " +
        //                                           "degree is " + degree.ToString();
        //                    DefaultReports.Add(defaultReport);
        //                    break;
        //            }
        //        }
        //    }

            
        //    return;
        //}


        /// <summary>
        /// 找到要分割的Vertex发出的两条相邻的半边，更进一步得到这两条边4个端点的Index，作为后面查找的根据
        /// </summary>
        /// <param name="PDeepCopy"></param>
        /// <param name="VertexToSplitIndex"></param>
        /// <returns></returns>
        private static List<int[,]> GetAllPossibleHStartAndEndIndexs(PlanktonMesh PDeepCopy, int VertexToSplitIndex)
        {
            List<int[,]> allPossibleHStartAndEndIndexs = new List<int[,]>();

            // 对于一个顶点来说的所有可能的一对半边的index（相邻的半边）
            // 点0周围的halfedgeindex列表，顺时针
            int[] halfEdgesIndexStartFromVertex = PDeepCopy.Vertices.GetHalfedges(VertexToSplitIndex);

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

                allPossibleHStartAndEndIndexs.Add(possibleHalfedgeVertexIndexs);
            }

            return allPossibleHStartAndEndIndexs;
        }

        /// <summary>
        /// 将选中的halfedgeIndex的半边进行分割（即这条半边的起点，分裂成了两个点）
        /// </summary>
        /// <param name="P">要进行操作的PlanktonMesh</param>
        /// <param name="possibleHStartAndEndIndexs">对应的那种可能的两个相邻的半边</param>
        /// <param name="splitParts">将半边分成几份</param>
        /// <param name="graphNodes">原来的graphNodes</param>
        /// <param name="newPNodes">更新后的graphNodes</param>
        /// <param name="newVertexIndex">新生成的点的Index</param>
        /// <returns></returns>
        public PlanktonMesh SplitEdgeIntoTwo(PlanktonMesh P, 
                                             int[,] possibleHStartAndEndIndexs,
                                             int splitParts,  
                                             List<GraphNode> graphNodes, 
                                             out List<GraphNode> newPNodes, 
                                             out int newVertexIndex)
        {
            #region 深拷贝
            PlanktonMesh pDeepCopy = new PlanktonMesh(P);

            // 深拷贝graphNodes
            List<GraphNode> graphNodesDeepCopy = new List<GraphNode>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                graphNodesDeepCopy.Add(new GraphNode(graphNodes[i]));
            }

            #endregion

            int halfedgeStartVertexIndex = possibleHStartAndEndIndexs[0, 0];
            int halfedgeEndVertexIndex = possibleHStartAndEndIndexs[0, 1];
            int anotherHalfedgeStartVertexIndex = possibleHStartAndEndIndexs[1, 0];
            int anotherHalfedgeEndVertexIndex = possibleHStartAndEndIndexs[1, 1];

            #region 计算要分裂的边两侧的Face，准备好两个Face的顶点列表，新增midVertex作为分裂产生的新顶点

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
                newPNodes = null;
                newVertexIndex = -1;
                return null;
            }

            int adjacentFaceIndex = pDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = pDeepCopy.Halfedges[pDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            //// 包围AdjacentFace的边
            //List<int> hiAroundAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[adjacentFaceIndex].FirstHalfedge).ToList<int>();
            //// 包围PairAdjacentFace的边
            //List<int> hiAroundPairAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[pairAdjacentFaceIndex].FirstHalfedge).ToList<int>();

            List<int> viAroundAdjacentFace = new List<int>();
            try
            {
                // 包围AdjacentFace的顶点
                viAroundAdjacentFace = pDeepCopy.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();
            }
            catch (ArgumentOutOfRangeException)
            {
                newPNodes = graphNodesDeepCopy;
                newVertexIndex = -1;
                return pDeepCopy;
            }

            List<int> viAroundPairAdjacentFace = new List<int>();
            try
            {
                // 包围PairAdjcentFace的顶点
                viAroundPairAdjacentFace = pDeepCopy.Faces.GetFaceVertices(pairAdjacentFaceIndex).ToList<int>();
            }
            catch (ArgumentOutOfRangeException)
            {
                newPNodes = graphNodesDeepCopy;
                newVertexIndex = -1;
                return pDeepCopy;
            }

            int startVertexIndex = halfedgeStartVertexIndex;
            int endVertexIndex = halfedgeEndVertexIndex;

            Point3d startVertex = pDeepCopy.Vertices[startVertexIndex].ToPoint3d();
            Point3d endVertex = pDeepCopy.Vertices[endVertexIndex].ToPoint3d();
            Point3d midVertex = (startVertex + endVertex) / 2;

            // 向P中增加顶点
            pDeepCopy.Vertices.Add(midVertex);
            int midVertexIndex = pDeepCopy.Vertices.Count - 1;
            newVertexIndex = midVertexIndex;
            #endregion

            #region 为splitEdge的操作更新图结构和Node对象

            /* 对startVertex的Node的相关属性进行修正，后续需要根据规则修改补充 */
            string startNodeLabel = graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeLabel;
            double startNodeArea = graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeArea;
            // double startNodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion;

            graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeLabel = graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeLabel + "0";
            graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeArea = graphNodesDeepCopy[startVertexIndex].NodeAttribute.NodeArea / splitParts;
            // nodes[startVertexIndex].NodeAttribute.NodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion / splitParts;

            // 对新生成的顶点，构造新的node
            GraphNodeAttribute nodeAttribute = new GraphNodeAttribute((startNodeLabel + (midVertexIndex - (graphNodes.Count - 1)).ToString()),
                                                             startNodeArea / splitParts);
            nodeAttribute.ConnectivityTable = new int[] { startVertexIndex, endVertexIndex };
            nodeAttribute.AdjacencyTable = new int[] { };
            GraphNode newNode = new GraphNode(midVertex, nodeAttribute, true);
            newNode.IsNewDecomposed = true;
            graphNodesDeepCopy.Add(newNode);

            // 输出out参数
            newPNodes = graphNodesDeepCopy;
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

            return newP;
        }

        /// <summary>
        /// 将一对半边（其中一条起点是vertexToSplitIndex）移动到一个新的顶点（targetVertexIndex）上，并且将这个顶点的度，补满到4
        /// </summary>
        /// <param name="P"></param>
        /// <param name="possibleHStartAndEndIndexs"></param>
        /// <param name="vertexToSplitIndex">要移动的一对半边中的一个，该半边的起点</param>
        /// <param name="targetVertexIndex">移动到哪一个顶点上，起点改为这个顶点</param>
        /// <param name="innerNodeToSplitIndexList">剩余的还未被分割的innerNodes</param>
        /// <param name="newEndVertexAfterResetHalfedgeList"></param>
        /// <returns></returns>
        public List<PlanktonMesh> ResetHalfedgeStart(PlanktonMesh P, 
                                                     int[,] possibleHStartAndEndIndexs,
                                                     int vertexToSplitIndex, 
                                                     int targetVertexIndex, 
                                                     List<int> innerNodeIndexs,
                                                     out List<int> newEndVertexAfterResetHalfedgeList)
        {
            #region 深拷贝
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);
            #endregion

            //// debug 显示每条边的起点与终点index
            //List<string> printHalfedgeStartAndEnd = new List<string>();
            //printHalfedgeStartAndEnd = UtilityFunctions.PrintHalfedgeStartAndEnd(PDeepCopy);

            int halfedgeStartVertexIndex = possibleHStartAndEndIndexs[0, 0];
            int halfedgeEndVertexIndex = possibleHStartAndEndIndexs[0, 1];
            int anotherHalfedgeStartVertexIndex = possibleHStartAndEndIndexs[1, 0];
            int anotherHalfedgeEndVertexIndex = possibleHStartAndEndIndexs[1, 1];


            #region 将一对半边移动到一个新的顶点上 

            // 找到anotherHalfedgeStartVertexIndex和anotherHalfedgeEndVertexIndex所对应的那条半边
            int halfedgeIndex = -1;
            foreach (int hfIndex in PDeepCopy.Vertices.GetHalfedges(anotherHalfedgeStartVertexIndex))
            {
                if (PDeepCopy.Halfedges.EndVertex(hfIndex) == anotherHalfedgeEndVertexIndex)
                {
                    halfedgeIndex = hfIndex;
                    break;
                }
            }
            if (halfedgeIndex == -1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ResetHalfedgeStart函数出错，找不到输入起点和终点所对应的半边。");
                newEndVertexAfterResetHalfedgeList = null;
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
                newEndVertexAfterResetHalfedgeList = null;
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
            //List<string> printFaces1 = new List<string>();
            //printFaces1 = UtilityFunctions.PrintFacesVertices(resetHalfEdgeStartP);
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




            // 如果新生成的半边的端点，是innerNode（包括初始的和新分裂产生的），那么返回null
            // 同时输出anotherEndVertexAfterResetHalfedgeList的列表，对应每种splitedFace的情况
            newEndVertexAfterResetHalfedgeList = new List<int>();
            for (int i = 0; i < splitedFaceWithOriginIndexLoL.Count; i++)
            {
                List<int> intersect = splitedFaceWithOriginIndexLoL[i].Intersect<int>(splitedFaceWithNewIndexLoL[i]).ToList<int>();
                intersect.Remove(targetVertexIndex);

                //if (innerNodeToSplitIndexList.Intersect<int>(intersect).Count() != 0)
                //{
                //    newEndVertexAfterResetHalfedgeList = null;
                //    return null;
                //}
                if (innerNodeIndexs.Contains(intersect[0]))
                {
                    newEndVertexAfterResetHalfedgeList = null;
                    return null;
                }


                newEndVertexAfterResetHalfedgeList.Add(intersect[0]);
            }






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
                    rebuildPFaceVertexOrder[i].Add(new List<int>());
                    if (j == whichFaceNeedSplitLaterIndex)
                    {
                        // rebuildPFaceVertexOrder[i].Add(splitedFaceWithOriginIndexLoL[i]);

                        rebuildPFaceVertexOrder[i][j].AddRange(splitedFaceWithOriginIndexLoL[i]);
                    }
                    else
                    {
                        int[] faceVertexOrder = resetHalfEdgeStartP.Faces.GetFaceVertices(j);
                        // rebuildPFaceVertexOrder[i].Add(faceVertexOrder.ToList<int>());

                        rebuildPFaceVertexOrder[i][j].AddRange(faceVertexOrder.ToList<int>());
                    }
                }
                // 注意这里还要额外加上split后新生成的面
                rebuildPFaceVertexOrder[i].Add(new List<int>());
                rebuildPFaceVertexOrder[i][resetHalfEdgeStartP.Faces.Count].AddRange(splitedFaceWithNewIndexLoL[i]);
            }

            //int iIndex = -1;
            //int jIndex = -1;
            //for (int i = 0; i < resetHalfEdgeStartP.Faces.Count; i++)
            //{
            //    rebuildPFaceVertexOrder.Add(new List<List<int>>());
            //    iIndex = i;
            //    for (int j = 0; j < splitedFaceWithOriginIndexLoL.Count; j++)
            //    {
            //        if (i == whichFaceNeedSplitLaterIndex)
            //        {
            //            rebuildPFaceVertexOrder[i].Add(splitedFaceWithOriginIndexLoL[j]);
            //            jIndex = j;
            //        }
            //        else
            //        {
            //            int[] faceVertexOrder = resetHalfEdgeStartP.Faces.GetFaceVertices(i);
            //            rebuildPFaceVertexOrder[i].Add(faceVertexOrder.ToList<int>());
            //        }
            //    }
                
            //}
            //// 注意这里还要额外加上split后新生成的面
            //rebuildPFaceVertexOrder.Add(new List<List<int>>());
            //rebuildPFaceVertexOrder[iIndex + 1].Add(splitedFaceWithNewIndexLoL[jIndex]);

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

            //// debug
            //List<List<string>> printFaces2 = new List<List<string>>();
            //for (int i = 0; i < rebulidNewHalfedgeP.Count; i++)
            //{
            //    printFaces2.Add(new List<string>());
            //    printFaces2[i] = UtilityFunctions.PrintFacesVertices(rebulidNewHalfedgeP[i]);
            //}
            #endregion

            for (int i = 0; i < rebulidNewHalfedgeP.Count; i++)
            {
                rebulidNewHalfedgeP[i].Compact();

            }

            for (int i = 0; i < rebulidNewHalfedgeP.Count; i++)
            {

            }

            return rebulidNewHalfedgeP;
        }


        public List<List<int>> RenewGraphTable(List<List<int>> graphTables,int[,] possibleHStartAndEndIndexs,int newVertexIndex, int newEndVertexAfterResetHalfedge)
        {
            List<List<int>> graphTablesDeepCopy = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                graphTablesDeepCopy.Add(new List<int>());
                graphTablesDeepCopy[i].AddRange(graphTables[i]);
            }

            int halfedgeStartVertexIndex = possibleHStartAndEndIndexs[0, 0];
            int halfedgeEndVertexIndex = possibleHStartAndEndIndexs[0, 1];
            int anotherHalfedgeStartVertexIndex = possibleHStartAndEndIndexs[1, 0];
            int anotherHalfedgeEndVertexIndex = possibleHStartAndEndIndexs[1, 1];

            // 更新图结构
            // 如果要分裂的这条边，是在原来的GraphLoL里的（注意不是PlanktonMeshGraphLoL）
            // 那么首尾点的邻接表上都要加上这个新增点的序号
            if (graphTablesDeepCopy[halfedgeStartVertexIndex].Contains(halfedgeEndVertexIndex))
            {
                // 所有的添加
                graphTablesDeepCopy[halfedgeStartVertexIndex].Add(newVertexIndex);
                graphTablesDeepCopy[halfedgeEndVertexIndex].Add(newVertexIndex);
                graphTablesDeepCopy.Add((new int[2] { halfedgeStartVertexIndex, halfedgeEndVertexIndex }).ToList<int>());

                // 如果另一条半边在原来的GraphLoL里，那么把anotherHalfedgeEndVertexIndex添加到新增点的邻接表中
                bool flag = false;
                if (graphTablesDeepCopy[halfedgeStartVertexIndex].Contains(anotherHalfedgeEndVertexIndex))
                {
                    graphTablesDeepCopy[newVertexIndex].Add(anotherHalfedgeEndVertexIndex);
                    graphTablesDeepCopy[anotherHalfedgeEndVertexIndex].Add(newVertexIndex);
                    flag = true;
                }

                // 所有的减去
                graphTablesDeepCopy[halfedgeStartVertexIndex].Remove(halfedgeEndVertexIndex);
                graphTablesDeepCopy[halfedgeEndVertexIndex].Remove(halfedgeStartVertexIndex);

                // 如果另一条半边在原来的GraphLoL里，那么把anotherHalfedgeEndVertexIndex的邻接表中应该移除halfedgeStartVertexIndex
                if (flag)
                {
                    graphTablesDeepCopy[halfedgeStartVertexIndex].Remove(anotherHalfedgeEndVertexIndex);
                    graphTablesDeepCopy[anotherHalfedgeEndVertexIndex].Remove(halfedgeStartVertexIndex);
                }

            }
            // 如果要分裂的这条边，不在原来的GraphLoL里的（注意不是PlanktonMeshGraphLoL）
            // 那么还需要判断另一条半边anotherHalfedge是不是在原来的GraphLoL里
            // 如果不在，那么什么都不用做
            // 如果在，就先在首点的邻接表上要加上这个新增点的序号，并且移除末点的序号；然后在末点的邻接表，移除首点
            // 然后在GraphLoL末尾，这个新增点的邻接表上，先加上首点的序号，再加上这个新增点连接的另一半（由ResetHalfedgeStart函数输出的anotherEndVertexAfterResetHalfedge）
            else
            {
                if (graphTablesDeepCopy[halfedgeStartVertexIndex].Contains(anotherHalfedgeEndVertexIndex))
                {
                    // 所有的添加
                    graphTablesDeepCopy[halfedgeStartVertexIndex].Add(newVertexIndex);
                    graphTablesDeepCopy[newEndVertexAfterResetHalfedge].Add(newVertexIndex);
                    graphTablesDeepCopy.Add((new int[2] { halfedgeStartVertexIndex, newEndVertexAfterResetHalfedge }).ToList<int>());

                    // 所有的减去
                    // graphTablesDeepCopy[halfedgeStartVertexIndex].Remove(halfedgeEndVertexIndex);
                    // graphTablesDeepCopy[halfedgeEndVertexIndex].Remove(halfedgeStartVertexIndex);
                    graphTablesDeepCopy[halfedgeStartVertexIndex].Remove(anotherHalfedgeEndVertexIndex);
                    graphTablesDeepCopy[anotherHalfedgeEndVertexIndex].Remove(halfedgeStartVertexIndex);
                }
            }
            return graphTablesDeepCopy;
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




        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new CompLabelAttribute(this);
        }

        public class CompLabelAttribute : GH_ComponentAttributes
        {
            public CompLabelAttribute(GhcDecomposeVertex component) : base(component) { }

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

                    using (GH_Capsule capsule1 = GH_Capsule.CreateCapsule(buttonRect1, GH_Palette.Black))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule1.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("CurrentLeafNode: {0} / {1}", ((GhcDecomposeVertex)Owner).CurrentLeafNodeIndex.ToString(), ((GhcDecomposeVertex)Owner).CountOfAllLeafNodes.ToString()),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.White,
                                        buttonRect1,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });

                    RectangleF buttonRect2 = new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);
                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect2.Inflate(-2.0f, -2.0f);

                    using (GH_Capsule capsule2 = GH_Capsule.CreateCapsule(buttonRect2, GH_Palette.Black))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule2.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("CurrentTreeBranch: {0} / {1}", ((GhcDecomposeVertex)Owner).CurrentDepth.ToString(), ((GhcDecomposeVertex)Owner).SelectedBranchDepth.ToString()),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.White,
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