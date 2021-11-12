using Grasshopper;
using Grasshopper.GUI.Gradient;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcTutteConvexDrawing : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the TutteConvexDrawing class.
        /// </summary>
        public GhcTutteConvexDrawing()
          : base("TutteConvexDrawing", "TutteDrawing",
              "Finds a unique convex drawing of a bi-connected graph. A bi-connected graph is a graph in which every vertex(node) is connected to others at least through two edges (links)",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddPlaneParameter("BasePlane", "BP", "生成Tutte嵌入的平面", GH_ParamAccess.item, Plane.WorldXY);
            // pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);
            // pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);

            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // pManager.AddGenericParameter("GraphNode", "GNode", "经过Tutte嵌入后得到的节点", GH_ParamAccess.list);
            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);

            pManager.AddCurveParameter("ConvexFaceBorders", "CFBorders", "Tutte嵌入后得到的TutteConvex列表", GH_ParamAccess.list);

            pManager.AddGenericParameter("SubGraph", "SG", "只包含InnerNode的子图", GH_ParamAccess.item);

            //pManager.AddIntegerParameter("SubGraph", "SubG", "InnerNode内部点之间的连接关系", GH_ParamAccess.tree);
            //pManager.AddGenericParameter("SubGraphNode", "SubGNode", "InnerNode内部节点", GH_ParamAccess.list);
            pManager.AddLineParameter("GraphEdge", "GEdge", "用直线来代表图结构中的连接关系", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 局部变量初始化
            Graph graph = new Graph();
            Graph subGraph = new Graph();

            Plane worldXY = Plane.WorldXY;
            #endregion

            if (DA.GetData("Graph", ref graph))
            {
                DA.GetData<Plane>("BasePlane", ref worldXY);
                #region 计算相关矩阵的前置条件

                int outerNodeCount = graph.OuterNodeCount;
                int innerNodeCount = graph.InnerNodeCount;
                List<int> outerNodeIndexList = graph.OuterNodeIndexList;
                List<int> innerNodeIndexList = graph.InnerNodeIndexList;

                List<List<int>> graphLoL = graph.GraphTables;
                // List<Node> graphNode = graph.GraphNodes;

                // 每个volume点，其连接的数量，即该节点的度
                List<int> innerNodeDegreeList = new List<int>();
                for (int i = 0; i < innerNodeCount; i++)
                {
                    innerNodeDegreeList.Add(graphLoL[i].Count);
                }
                #endregion

                #region 计算矩阵 P(outer)
                // 回0,0,0附近进行计算
                Point3d boundaryCenter = new Point3d(0.0, 0.0, 0.0);
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (!graph.GraphNodes[i].IsInner)
                    {
                        boundaryCenter += graph.GraphNodes[i].NodeVertex;
                    }
                }
                boundaryCenter = boundaryCenter / outerNodeCount;

                //矩阵 P(outer)，回0,0,0附近进行计算
                Matrix P_outer = new Matrix(outerNodeCount, 2);
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (!graph.GraphNodes[i].IsInner)
                    {
                        P_outer[outerNodeIndexList.IndexOf(i), 0] = graph.GraphNodes[i].NodeVertex.X - boundaryCenter.X;
                        P_outer[outerNodeIndexList.IndexOf(i), 1] = graph.GraphNodes[i].NodeVertex.Y - boundaryCenter.Y;
                    }
                }

                #endregion

                #region 整个大的矩阵，包含四个部分  inner行-inner列  inner行-outer列：Q  outer行-inner列：Q(上标T)  outer行-outer列

                Matrix inner_innerM;
                // 矩阵Q：inner_OuterAdjacencyMatrix
                Matrix inner_outerM;
                Matrix outer_innerM;
                Matrix outer_outerM;
                Matrix wholeM = GraphLoLToMatrix(graphLoL, 
                                 innerNodeIndexList, 
                                 outerNodeIndexList, 
                                 out inner_innerM, 
                                 out inner_outerM, 
                                 out outer_innerM, 
                                 out outer_outerM);
                #endregion

                List<string> Mstring = PrintMatrix(wholeM);

                #region 矩阵L(下标1)：inner_InnerLaplacianMatrix
                //DataTree<int> inner_InnerAdjacencyDataTree = SubAdjacencyMatrix(wholeAdjacencyDataTree, volumeNodeIndexList);
                //Matrix inner_InnerAdjacencyMatrix = ToGHMatrix(inner_InnerAdjacencyDataTree);
                // 拉普拉斯矩阵
                //DataTree<int> inner_InnerLaplacianDataTree = LaplacianMatrix(inner_InnerAdjacencyDataTree, volumeNodeDegreeList);
                //Matrix inner_InnerLaplacianMatrix = ToGHMatrix(inner_InnerLaplacianDataTree);

                Matrix inner_InnerLaplacianMatrix = LaplacianMatrix(inner_innerM, innerNodeDegreeList);
                #endregion

                List<string> LMstring = PrintMatrix(inner_InnerLaplacianMatrix);

                #region 拉普拉斯矩阵的逆
                inner_InnerLaplacianMatrix.Invert(0.0);
                // 拉普拉斯矩阵的逆 L(下标1)(上标-1)：inverse_Inner_InnerLaplacianMatrix
                Matrix inverse_Inner_InnerLaplacianMatrix = inner_InnerLaplacianMatrix;
                #endregion

                #region 求解矩阵 P(inner)
                Matrix P_inner = new Matrix(innerNodeCount, 2);
                P_inner = inverse_Inner_InnerLaplacianMatrix * inner_outerM * P_outer;
                // 原插件这里没乘负号
                // P_inner.Scale(-1.0);
                #endregion

                #region 从P_inner矩阵中生成新的inner点的坐标
                List<Point3d> newInnerPoints = new List<Point3d>();
                Point3d innerPoint;
                for (int i = 0; i < innerNodeCount; i++)
                {
                    innerPoint = new Point3d(P_inner[i, 0] + boundaryCenter.X, P_inner[i, 1] + boundaryCenter.Y, 0.0);
                    newInnerPoints.Add(innerPoint);
                }
                #endregion

                #region Relocate重定位，形成Graph的Node
                // 将新的inner点与原来的outer点合并成一个新的列表
                List<Point3d> newNodePoints = new List<Point3d>();
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (graph.GraphNodes[i].IsInner)
                    {
                        newNodePoints.Add(newInnerPoints[innerNodeIndexList.IndexOf(i)]);
                    }
                    else
                    {
                        newNodePoints.Add(graph.GraphNodes[i].NodeVertex);
                    }
                }


                // 将新的inner+outer列表重新定位在基于worldXY变量的另一个地方
                List<Point3d> newNodePointsRelocated = UtilityFunctions.Relocate(newNodePoints, Plane.WorldXY, worldXY);
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    graph.GraphNodes[i].NodeVertex = newNodePointsRelocated[i];
                }
                #endregion
                // DA.SetDataList("GraphNode", nodes);
                DA.SetData("Graph", graph);




                #region 绘制Graph的edge
                // Graph的所有Edge，在新的位置上画Edge
                List<Line> graphEdge = GraphEdgeList(graph.GraphTables, newNodePointsRelocated);
                DA.SetDataList("GraphEdge", graphEdge);
                #endregion

                #region 绘制Graph形成的Convex
                // Graph形成的每个convex，在新的位置上画Convex
                List<Point3d> sortedBoundaryNodePointsRelocated = new List<Point3d>();
                for (int i = 0; i < outerNodeCount; i++)
                {
                    sortedBoundaryNodePointsRelocated.Add(newNodePointsRelocated[outerNodeIndexList[i]]);
                }
                List<Polyline> allSplitFaceBoundarys = GMeshFaceBoundaries(sortedBoundaryNodePointsRelocated, graphEdge, worldXY);
                #endregion
                DA.SetDataList("ConvexFaceBorders", allSplitFaceBoundarys);

                #region 输出subGraph

                List<List<int>> subGraphLoL = new List<List<int>>();
                List<Node> subGraphNodes = new List<Node>();
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    
                    if (graph.GraphNodes[i].IsInner)
                    {
                        subGraphLoL.Add(new List<int>());
                        int index = subGraphLoL.Count - 1;
                        for (int j = 0; j < graph.GraphNodes.Count; j++)
                        {
                            if (graph.GraphNodes[j].IsInner)
                            {
                                if (wholeM[i, j] == 1)
                                {
                                    subGraphLoL[index].Add(j);
                                }
                            }
                        }
                        subGraphNodes.Add(graph.GraphNodes[i]);
                    }

                    
                }
                #endregion

                subGraph = new Graph(subGraphNodes, subGraphLoL);
                DA.SetData("SubGraph", subGraph);
                #region 输出subGraphNode
                // DA.SetDataTree(2, subGraph);
                //for (int i = 0; i < innerNodeCount; i++)
                //{
                //    subNodes.Add(nodes[i]);
                //}

                
                #endregion
                // DA.SetDataList("SubGraphNode", subNodes);
            }
        }

        public Matrix GraphLoLToMatrix(List<List<int>> graphLoL,
                                     List<int> innerNodeIndexList,
                                     List<int> outerNodeIndexList,
                                     out Matrix inner_innerM, 
                                     out Matrix inner_outerM, 
                                     out Matrix outer_innerM, 
                                     out Matrix outer_outerM)
        {
            Matrix wholeMatrix = new Matrix(graphLoL.Count, graphLoL.Count);

            List<List<int>> sortedGraphLoL = new List<List<int>>();
            for (int i = 0; i < innerNodeIndexList.Count; i++)
            {
                sortedGraphLoL.Add(new List<int>(graphLoL[innerNodeIndexList[i]]));
            }
            for (int i = 0; i < outerNodeIndexList.Count; i++)
            {
                sortedGraphLoL.Add(new List<int>(graphLoL[outerNodeIndexList[i]]));
            }

            for (int i = 0; i < sortedGraphLoL.Count; i++)
            {
                for (int j = 0; j < sortedGraphLoL.Count; j++)
                {
                    if (sortedGraphLoL[i].Contains(j))
                    {
                        wholeMatrix[i, j] = 1;
                    }
                    else
                    {
                        wholeMatrix[i, j] = 0;
                    }

                }
            }

            inner_innerM = new Matrix(innerNodeIndexList.Count, innerNodeIndexList.Count);
            for (int i = 0; i < innerNodeIndexList.Count; i++)
            {
                for (int j = 0; j < innerNodeIndexList.Count; j++)
                {
                    inner_innerM[i, j] = wholeMatrix[i, j];
                }
            }

            inner_outerM = new Matrix(innerNodeIndexList.Count, outerNodeIndexList.Count);
            for (int i = 0; i < innerNodeIndexList.Count; i++)
            {
                for (int j = innerNodeIndexList.Count; j < graphLoL.Count; j++)
                {
                    inner_outerM[i, j - innerNodeIndexList.Count] = wholeMatrix[i, j];
                }
            }

            outer_innerM = new Matrix(outerNodeIndexList.Count, innerNodeIndexList.Count);
            for (int i = innerNodeIndexList.Count; i < graphLoL.Count; i++)
            {
                for (int j = 0; j < innerNodeIndexList.Count; j++)
                {
                    outer_innerM[i - innerNodeIndexList.Count, j] = wholeMatrix[i, j];
                }
            }

            outer_outerM = new Matrix(outerNodeIndexList.Count, outerNodeIndexList.Count);
            for (int i = innerNodeIndexList.Count; i < graphLoL.Count; i++)
            {
                for (int j = innerNodeIndexList.Count; j < graphLoL.Count; j++)
                {
                    outer_outerM[i - innerNodeIndexList.Count, j - innerNodeIndexList.Count] = wholeMatrix[i, j];
                }
            }

            return wholeMatrix;
        }




        /// <summary>
        /// 调和矩阵L Harmonic Matrix，又称拉普拉斯矩阵 Laplacian Matrix。是图的矩阵表示
        /// 在图论和计算机科学中，邻接矩阵A（英语：adjacency matrix）是一种方阵，用来表示有限图。它的每个元素代表各点之间是否有边相连。
        /// 在数学领域图论中，度数矩阵D是一个对角矩阵 ，其中包含的信息为的每一个顶点的度数，也就是说，每个顶点相邻的边数[1] 它可以和邻接矩阵一起使用以构造图的拉普拉斯算子矩阵。
        /// L = D - A
        /// </summary>
        /// <param name="M">邻接矩阵 Adjacency Matrix，对角线都是0，其他部分有数据0或1</param>
        /// <param name="degrees">度数矩阵 Degrees Matrix，因为是只有对角线有数据，其他部分都是零，所以可以用list存储</param>
        /// <returns></returns>
        public Matrix LaplacianMatrix(Matrix M, List<int> degrees)
        {
            if (!M.IsSquare)
            {
                return null;
            }
            if (M.RowCount != degrees.Count)
            {
                return null;
            }
            
            Matrix LaplacianM = new Matrix(M.RowCount, M.ColumnCount);
            for (int i = 0; i < M.RowCount; i++)
            {
                for (int j = 0; j < M.ColumnCount; j++)
                {
                    if (i == j)
                    {
                        LaplacianM[i, j] = degrees[j] - 0;
                    }
                    else
                    {
                        LaplacianM[i, j] = 0 - degrees[j];
                    }
                }
            }

            return LaplacianM;
        }

        public List<string> PrintMatrix(Matrix M)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < M.RowCount; i++)
            {
                string str = string.Format("{0}行:", i);
                for (int j = 0; j < M.ColumnCount; j++)
                {
                    str += string.Format("{0},", M[i, j]);
                }
                output.Add(str);
            }

            return output;
        }


        /// <summary>
        /// 将排好序的点连接成多段线（NurbsCurve类型）
        /// </summary>
        /// <param name="sortedBoundaryPoints"></param>
        /// <param name="basePlane"></param>
        /// <returns></returns>
        public Polyline SortedPointsToPolyline(List<Point3d> sortedBoundaryPoints, Plane basePlane)
        {
            // 把startPoint放到最后作为endPoint，这样polyline才会闭合
            sortedBoundaryPoints.Add(sortedBoundaryPoints[0]);

            Polyline convex = new Polyline(sortedBoundaryPoints);
            return convex;
        }


        public List<Polyline> GMeshFaceBoundaries(List<Point3d> sortedBoundaryPoints, List<Line> edges, Plane basePlane)
        {
            List<Curve> splitLines = new List<Curve>();
            foreach (Line line in edges)
            {
                splitLines.Add(line.ToNurbsCurve());
            }

            Polyline polyline = SortedPointsToPolyline(sortedBoundaryPoints, basePlane);
            Curve curve = polyline.ToNurbsCurve();

            Brep[] PlanarBrepArray = Brep.CreatePlanarBreps(curve, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            Brep planarBreps = PlanarBrepArray[0];
            BrepFace planarBrepFace = planarBreps.Faces[0];
            Brep splittedPlanarBrepFaces = planarBrepFace.Split(splitLines, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            
            List<Polyline> allConvexPolyline = new List<Polyline>();
            foreach (BrepFace brepFace in splittedPlanarBrepFaces.Faces)
            {
                Brep singleSplitFaceBrep = brepFace.DuplicateFace(true);
                Point3d[] singleSplitFaceVertices = singleSplitFaceBrep.DuplicateVertices();
                Curve splitFaceConvexCure = Curve.JoinCurves(singleSplitFaceBrep.Faces[0].DuplicateFace(true).DuplicateEdgeCurves())[0];
                Polyline splitFaceConvexPolyline = null;
                splitFaceConvexCure.TryGetPolyline(out splitFaceConvexPolyline);
                allConvexPolyline.Add(splitFaceConvexPolyline);
            }

            return allConvexPolyline;
        }

        public List<Line> GraphEdgeList(List<List<int>> graphLoL, List<Point3d> graphVertices)
        {
            List<Line> list = new List<Line>();

            for (int i = 0; i < graphLoL.Count; i++)
            {
                for (int j = 0; j < graphLoL[i].Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graphLoL[i][j]]);
                    list.Add(item);
                }
            }

            return list;
        }


        /// <summary>
        /// 预览模式为WireFrame模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);
        }

        /// <summary>
        /// 预览模式为Shaded模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            base.DrawViewportMeshes(args);
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
            get { return new Guid("97181a6b-fd35-4eb5-9938-4e7a7e7b91ab"); }
        }
    }
}