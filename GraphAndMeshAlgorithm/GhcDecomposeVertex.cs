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
using System.Reflection;

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
            RhinoMesh = new Mesh();
            
            ConvexPolylinesPoints = new List<List<Point3d>>();
            SelectedTriangleMeshEdges = new List<Line>();

            InnerNodeTextDot = new List<TextDot>();
            OuterNodeTextDot = new List<TextDot>();

            DottedCurve = new List<Curve>();
        }

        private int Thickness;

        private Mesh RhinoMesh;
        // private int IndexOfTriangularMesh;

        private List<List<Point3d>> ConvexPolylinesPoints;
        private List<Line> SelectedTriangleMeshEdges;

        private List<TextDot> InnerNodeTextDot;
        private List<TextDot> OuterNodeTextDot;

        private List<Curve> DottedCurve;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);

            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("NewTriangleHalfedgeMesh", "NTHMesh", "新生成的三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddGenericParameter("NewGraph", "NG", "描述VolumeNode和BoundaryNode的所有连接关系的图结构(包括新分裂产生的点)", GH_ParamAccess.tree);

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

            //List<int> innerNodeIndexs = new List<int>();
            //for (int i = 0; i < innerNodeCount; i++)
            //{
            //    innerNodeIndexs.Add(i);
            //}
            //List<int> outerNodeIndexs = new List<int>();
            //for (int i = 0; i < outerNodeCount; i++)
            //{
            //    outerNodeIndexs.Add(i + innerNodeCount);
            //}

            PlanktonMesh P = new PlanktonMesh();

            List<Node> nodes = new List<Node>();

            PlanktonMesh edgeSplitedP = new PlanktonMesh();
            PlanktonMesh edgeResetStartP = new PlanktonMesh();

            GH_Structure<GH_Integer> gh_Structure_graph = null;
            DataTree<int> graph = new DataTree<int>();
            List<List<int>> graphLoL = new List<List<int>>();

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh", ref P)
                && DA.GetDataList<Node>("GraphNode", nodes)
                && DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_graph))
            {
                ConvexPolylinesPoints.Clear();
                SelectedTriangleMeshEdges.Clear();
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();

                

                // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>，再转化为LoL
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_graph, ref graph);
                graphLoL = UtilityFunctions.DataTreeToLoL<int>(graph);

                // 利用PlanktonMesh为复制准备的构造函数，进行深拷贝
                PlanktonMesh PDeepCopy = new PlanktonMesh(P);
                

                // 对于每个InnerNode来说，判断它的度是否大于4
                for (int i = 0; i < PDeepCopy.Vertices.Count - outerNodeCount; i++)
                {
                    
                }

                int degree = PDeepCopy.Vertices.GetValence(1);

                // 点1周围的halfedgeindex列表，顺时针
                int[] halfEdgesIndexStartFromVertex = PDeepCopy.Vertices.GetHalfedges(1);

                switch (degree)
                {
                    case 5:
                        // SplitVertex(PDeepCopy, 1, innerNodeIndexs);
                        // P.Vertices.SetVertex(P.Vertices.Count-1, P.Vertices[0].ToXYZ().X, P.Vertices[0].ToXYZ().Y + (float)0.5, P.Vertices[0].ToXYZ().Z);

                        //Point3d startPoint = PDeepCopy.Vertices[1].ToPoint3d();
                        //Point3d endPoint = PDeepCopy.Vertices[PDeepCopy.Halfedges.EndVertex(15)].ToPoint3d();

                        //int newHalfEdgeIndex = PDeepCopy.Vertices.SplitVertex(15, 16);
                        // PDeepCopy.Vertices.SetVertex(PDeepCopy.Halfedges[15].StartVertex, (startPoint + endPoint) / 2);

                        edgeSplitedP = SplitEdgeIntoTwo(PDeepCopy, 15, graphLoL, nodes, 2);
                        edgeResetStartP = ResetHalfedgeStart(edgeSplitedP, 35, 1, 9);

                        break;

                    default:
                        break;
                }

                DA.SetData("NewTriangleHalfedgeMesh", edgeResetStartP);


                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(edgeResetStartP);
                DA.SetDataList("DebugVerticesOutput", printVertices);

                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(edgeResetStartP);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);

                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFaces(edgeResetStartP);
                DA.SetDataList("DebugFacesOutput", printFaces);

                List<string> printFacesHalfedge = new List<string>();
                printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(edgeResetStartP);
                DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);


                #region 可视化部分
                // 设置用于可视化的TextDot
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, nodes[i].NodeAttribute.NodeLabel), nodes[i].NodeVertex);
                        InnerNodeTextDot.Add(textDot);
                    }
                    else
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, nodes[i].NodeAttribute.NodeLabel), nodes[i].NodeVertex);
                        OuterNodeTextDot.Add(textDot);
                    }
                }


                // 转换为RhinoMesh，在DrawViewportMeshes绘制选中的mesh
                RhinoMesh = RhinoSupport.ToRhinoMesh(edgeResetStartP);
                // 在DrawViewportWires绘制mesh的edge
                for (int i = 0; i < RhinoMesh.TopologyEdges.Count; i++)
                {
                    SelectedTriangleMeshEdges.Add(RhinoMesh.TopologyEdges.EdgeLine(i));
                }

                DottedCurve.Clear();

                double[] pattern = { 1.0 };
                for (int i = 0; i < SelectedTriangleMeshEdges.Count; i++)
                {
                    IEnumerable<Curve> segments = UtilityFunctions.ApplyDashPattern(SelectedTriangleMeshEdges[i].ToNurbsCurve(), pattern);
                    foreach (Curve segment in segments)
                    {
                        DottedCurve.Add(segment);
                    }
                }
                #endregion
            }
        }


        /// <summary>
        /// 将选中的halfedgeIndex的半边进行分割（即这条半边的起点，分裂成了两个点）
        /// </summary>
        /// <param name="P"></param>
        /// <param name="halfedgeIndex"></param>
        /// <param name="graphLoL"></param>
        /// <param name="nodes"></param>
        /// <param name="splitParts"></param>
        /// <returns></returns>
        public PlanktonMesh SplitEdgeIntoTwo(PlanktonMesh P, int halfedgeIndex, List<List<int>> graphLoL, List<Node> nodes, int splitParts)
        {
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);

            #region 计算要分裂的边两侧的Face，准备好两个Face的顶点列表，新增midVertex作为分裂产生的新顶点
            //PlanktonVertexList vertices = PDeepCopy.Vertices;
            //PlanktonHalfEdgeList halfedges = PDeepCopy.Halfedges;
            //PlanktonFaceList faces = PDeepCopy.Faces;

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

            int startVertexIndex = PDeepCopy.Halfedges[halfedgeIndex].StartVertex;
            int endVertexIndex = PDeepCopy.Halfedges.EndVertex(halfedgeIndex);

            Point3d startVertex = PDeepCopy.Vertices[startVertexIndex].ToPoint3d();
            Point3d endVertex = PDeepCopy.Vertices[endVertexIndex].ToPoint3d();
            Point3d midVertex = (startVertex + endVertex) / 2;

            // 向P中增加顶点
            PDeepCopy.Vertices.Add(midVertex);
            int midVertexIndex = PDeepCopy.Vertices.Count - 1;
            #endregion

            #region 为splitEdge的操作更新图结构和Node对象
            // 更新图结构
            for (int i = 0; i < graphLoL.Count; i++)
            {
                if (i == startVertexIndex)
                {
                    graphLoL[i].Add(midVertexIndex);
                }
                if (i == endVertexIndex)
                {
                    graphLoL[i].Add(midVertexIndex);
                }
            }
            graphLoL.Add((new int[2] { startVertexIndex, endVertexIndex }).ToList<int>());

            // 对startVertex的Node的相关属性进行修正，后续需要根据规则修改补充
            string startNodeLabel = nodes[startVertexIndex].NodeAttribute.NodeLabel;
            double startNodeArea = nodes[startVertexIndex].NodeAttribute.NodeArea;
            // double startNodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion;

            nodes[startVertexIndex].NodeAttribute.NodeLabel = nodes[startVertexIndex].NodeAttribute.NodeLabel + "0";
            nodes[startVertexIndex].NodeAttribute.NodeArea = nodes[startVertexIndex].NodeAttribute.NodeArea / splitParts;
            // nodes[startVertexIndex].NodeAttribute.NodeAreaProportion = nodes[startVertexIndex].NodeAttribute.NodeAreaProportion / splitParts;

            // 对新生成的顶点，构造新的node
            NodeAttribute nodeAttribute = new NodeAttribute((startNodeLabel + (midVertexIndex - 8).ToString()),
                                                            startNodeArea / splitParts);
            nodes.Add(new Node(midVertex, nodeAttribute, true));
            #endregion

            #region 构造添加顶点后的新Face
            for (int i = 0; i < viAroundAdjacentFace.Count; i++)
            {
                if (viAroundAdjacentFace[i] == startVertexIndex)
                {
                    viAroundAdjacentFace.Insert(i + 1, midVertexIndex);
                }
            }
            for (int i = 0; i < viAroundPairAdjacentFace.Count; i++)
            {
                if (viAroundPairAdjacentFace[i] == endVertexIndex)
                {
                    viAroundPairAdjacentFace.Insert(i + 1, midVertexIndex);
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
        public PlanktonMesh ResetHalfedgeStart(PlanktonMesh P, int halfedgeIndex, int vertexToSplitIndex, int targetVertexIndex)
        {
            /* 深拷贝 */
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);


            /* 将一对半边移动到一个新的顶点上 
             *
             */
            int pairHalfedgeIndex = PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex);

            int adjacentFaceIndex = PDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = PDeepCopy.Halfedges[PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            // 包围AdjacentFace的顶点
            List<int> viAroundAdjacentFace = new List<int>();
            int currentHalfedgeIndex = halfedgeIndex;
            do
            {
                viAroundAdjacentFace.Add(PDeepCopy.Halfedges[currentHalfedgeIndex].StartVertex);
                currentHalfedgeIndex = PDeepCopy.Halfedges[currentHalfedgeIndex].NextHalfedge;
            } while (currentHalfedgeIndex != halfedgeIndex);
            // 包围PairAdjacentFace的顶点
            List<int> viAroundPairAdjacentFace = new List<int>();
            currentHalfedgeIndex = pairHalfedgeIndex;
            do
            {
                viAroundPairAdjacentFace.Add(PDeepCopy.Halfedges.EndVertex(currentHalfedgeIndex));
                currentHalfedgeIndex = PDeepCopy.Halfedges[currentHalfedgeIndex].NextHalfedge;
            } while (currentHalfedgeIndex != pairHalfedgeIndex);

            // 修改包围AdjacentFace的顶点列表以及包围PairAdjacentFace的顶点列表
            // 用目标顶点属于哪个面，来判断该怎么修改viAroundAdjacentFace和viAroundPairAdjacentFace列表
            bool targetVertexBelongsToAdjacent = viAroundAdjacentFace.Contains(targetVertexIndex);
            bool targetVertexBelongsToPairAdjacent = viAroundPairAdjacentFace.Contains(targetVertexIndex);

            // 当目标顶点targetVertex同时不属于这两个Face时报错
            if (!targetVertexBelongsToAdjacent && !targetVertexBelongsToPairAdjacent)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的targetVertex " + vertexToSplitIndex.ToString() + " 不属于要改变的面");
                return null;
            }
            // 当目标顶点targetVertex属于pairAdjacentFace时
            else if (!targetVertexBelongsToAdjacent && targetVertexBelongsToPairAdjacent)
            {
                viAroundAdjacentFace.Insert(1, targetVertexIndex);
                viAroundPairAdjacentFace.RemoveAt(0);
            }
            // 当目标顶点targetVertex属于AdjacentFace时
            else
            {
                viAroundAdjacentFace.RemoveAt(0);
                viAroundPairAdjacentFace.Insert(1, targetVertexIndex);
            }

            // 为了记录哪一个面在下面的，将度补满的过程中，需要再进行split
            List<bool> whichFaceNeedSplit = new List<bool>();
            whichFaceNeedSplit.Add(!targetVertexBelongsToAdjacent);
            whichFaceNeedSplit.Add(!targetVertexBelongsToPairAdjacent);
            int whichFaceNeedSplitIndex = -1;

            /* 构建resetHalfEdgeStartP，用来存储将一对半边移动后形成的新mesh
             * 同时将未改动的PDeepCopy中的其他顶点和面的数据，转移过来
             */

            // 建立需要进行置换的面的index列表和VertexIndex列表，方便转移P的Face属性
            List<int> needChangeFaceIndexs = new List<int>();
            needChangeFaceIndexs.Add(adjacentFaceIndex);
            needChangeFaceIndexs.Add(pairAdjacentFaceIndex);
            List<List<int>> needChangeViAround = new List<List<int>>();
            needChangeViAround.Add(viAroundAdjacentFace);
            needChangeViAround.Add(viAroundPairAdjacentFace);

            // 转移P的Vertex属性
            List<PlanktonXYZ> resetPPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
            {
                resetPPlanktonVertex.Add(PDeepCopy.Vertices[i].ToXYZ());
            }
            // 转移P的Face属性
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
                        whichFaceNeedSplitIndex = i;
                    }
                }
                else
                {
                    resetPFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = PDeepCopy.Faces.GetFaceVertices(i);
                    resetPFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }

            // 用转移的P的Vertex属性和P的Face属性来构造新的PlanktonMesh newP
            PlanktonMesh resetHalfEdgeStartP = new PlanktonMesh();
            resetHalfEdgeStartP.Vertices.AddVertices(resetPPlanktonVertex);
            resetHalfEdgeStartP.Faces.AddFaces(resetPFaceVertexOrder);

            // debug
            List<string> printFaces1 = new List<string>();
            printFaces1 = UtilityFunctions.PrintFaces(resetHalfEdgeStartP);


            /* 为了将新得到的顶点的度补满，需要对一个面进行新的分割，形成新的一对半边
             * 
             */
            List<int> splitedFaceWithOriginIndex = new List<int>();
            List<int> splitedFaceWithNewIndex = new List<int>();

            // 用vertexToSplit属于哪个面，来判断该怎么修改经过半边移动后新的viAroundAdjacentFace和viAroundPairAdjacentFace列表
            bool vertexToSplitBelongsToAdjacent = viAroundAdjacentFace.Contains(vertexToSplitIndex);
            bool vertexToSplitBelongsToPairAdjacent = viAroundPairAdjacentFace.Contains(vertexToSplitIndex);

            // 如果分裂顶点vertexToSplit，属于pairAdjacentFace时
            if (vertexToSplitBelongsToAdjacent && !vertexToSplitBelongsToPairAdjacent)
            {
                int currentIndex = viAroundAdjacentFace.IndexOf(targetVertexIndex);
                int iteration = 0;
                do
                {
                    splitedFaceWithOriginIndex.Add(viAroundAdjacentFace[currentIndex]);
                    if (currentIndex + 1 > viAroundAdjacentFace.Count - 1)
                    {
                        currentIndex = 0;
                    }
                    else
                    {
                        currentIndex++;
                    }
                    iteration++;
                } while (iteration < 3);

                // 把最后多做的一次currentIndex++撤回掉
                if (currentIndex == 0)
                {
                    currentIndex = viAroundAdjacentFace.Count - 1;
                }
                else
                {
                    currentIndex -= 1;
                }

                do
                {
                    splitedFaceWithNewIndex.Add(viAroundAdjacentFace[currentIndex]);
                    if (currentIndex + 1 > viAroundAdjacentFace.Count - 1)
                    {
                        currentIndex = 0;
                    }
                    else
                    {
                        currentIndex++;
                    }
                    iteration++;
                } while (iteration < 6);
            }
            // 如果分裂顶点vertexToSplit，属于pairAdjacentFace时
            else
            {
                int currentIndex = viAroundAdjacentFace.IndexOf(targetVertexIndex);
                int iteration = 0;
                do
                {
                    splitedFaceWithOriginIndex.Add(viAroundPairAdjacentFace[currentIndex]);
                    if (currentIndex + 1 > viAroundPairAdjacentFace.Count - 1)
                    {
                        currentIndex = 0;
                    }
                    else
                    {
                        currentIndex++;
                    }
                    iteration++;
                } while (iteration < 3);

                // 把最后多做的一次currentIndex++撤回掉
                if (currentIndex == 0)
                {
                    currentIndex = viAroundAdjacentFace.Count - 1;
                }
                else
                {
                    currentIndex -= 1;
                }

                do
                {
                    splitedFaceWithNewIndex.Add(viAroundPairAdjacentFace[currentIndex]);
                    if (currentIndex + 1 > viAroundPairAdjacentFace.Count - 1)
                    {
                        currentIndex = 0;
                    }
                    else
                    {
                        currentIndex++;
                    }
                    iteration++;
                } while (iteration < 6);
            }

            // 转移P的Vertex属性
            List<PlanktonXYZ> rebuildPPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < resetHalfEdgeStartP.Vertices.Count; i++)
            {
                rebuildPPlanktonVertex.Add(resetHalfEdgeStartP.Vertices[i].ToXYZ());
            }
            // 转移P的Face属性
            List<List<int>> rebuildPFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < resetHalfEdgeStartP.Faces.Count; i++)
            {
                if (i == whichFaceNeedSplitIndex)
                {
                    rebuildPFaceVertexOrder.Add(new List<int>());
                    rebuildPFaceVertexOrder[i].AddRange(splitedFaceWithOriginIndex);
                }
                else
                {
                    rebuildPFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = resetHalfEdgeStartP.Faces.GetFaceVertices(i);
                    rebuildPFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            rebuildPFaceVertexOrder.Add(splitedFaceWithNewIndex);

            PlanktonMesh rebulidNewHalfedgeP = new PlanktonMesh();
            rebulidNewHalfedgeP.Vertices.AddVertices(rebuildPPlanktonVertex);
            // rebulidNewHalfedgeP.Faces.AddFaces(rpFaceVertexOrder);
            for (int i = 0; i < rebuildPFaceVertexOrder.Count; i++)
            {
                rebulidNewHalfedgeP.Faces.AddFace(rebuildPFaceVertexOrder[i]);
            }

            List<string> printFaces2 = new List<string>();
            printFaces2 = UtilityFunctions.PrintFaces(rebulidNewHalfedgeP);

            return rebulidNewHalfedgeP;
        }

        /// <summary>
        /// 预览模式为WireFrame模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);

            for (int i = 0; i < DottedCurve.Count; i++)
            {
                args.Display.DrawCurve(DottedCurve[i], Color.DarkGreen, Thickness);
            }
            //args.Display.EnableDepthTesting(true);

            // 后画实线
            for (int i = 0; i < ConvexPolylinesPoints.Count; i++)
            {
                args.Display.DrawPolyline(ConvexPolylinesPoints[i], Color.BlueViolet, Thickness);
            }

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

            args.Display.DrawMeshShaded(RhinoMesh, new Rhino.Display.DisplayMaterial(Color.White, 0));

            for (int i = 0; i < DottedCurve.Count; i++)
            {
                args.Display.DrawCurve(DottedCurve[i], Color.DarkGreen, Thickness);
            }

            // 后画实线
            for (int i = 0; i < ConvexPolylinesPoints.Count; i++)
            {
                args.Display.DrawPolyline(ConvexPolylinesPoints[i], Color.BlueViolet, Thickness);
            }

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