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
            // 0
            pManager.AddPlaneParameter("BasePlane", "BP", "生成Tutte嵌入的平面", GH_ParamAccess.item, Plane.WorldXY);
            // 1
            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // 0
            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);
            // 1
            pManager.AddCurveParameter("ConvexFaceBorders", "CFBorders", "Tutte嵌入后得到的TutteConvex列表", GH_ParamAccess.list);
            // 2
            pManager.AddGenericParameter("SubGraph", "SG", "只包含InnerNode的子图", GH_ParamAccess.item);
            // 3
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
            Plane targetPlane = Plane.WorldXY;
            #endregion

            if (DA.GetData<Graph>("Graph", ref graph))
            {
                DA.GetData<Plane>("BasePlane", ref targetPlane);

                List<List<int>> graphLoL = graph.GraphTables;
                List<Node> graphNodes = graph.GraphNodes;

                int innerNodeCount = graph.InnerNodeCount;
                int outerNodeCount = graph.OuterNodeCount;
                List<int> outerNodeIndexList = graph.OuterNodeIndexList;
                List<int> innerNodeIndexList = graph.InnerNodeIndexList;

                DataTree<int> graphDT = UtilityFunctions.LoLToDataTree<int>(graphLoL);

                #region 计算相关矩阵的前置条件
                // 每个volume点，其连接的数量，即该节点的度
                List<int> innerNodeDegreeList = new List<int>();
                for (int i = 0; i < innerNodeCount; i++)
                {
                    innerNodeDegreeList.Add(graphDT.Branch(i).Count);
                }

                // 每个volume点与其他点的连接关系，包括其他volume点和boundary点
                DataTree<int> volumeNodeGraph = new DataTree<int>();
                for (int i = 0; i < innerNodeCount; i++)
                {
                    volumeNodeGraph.EnsurePath(i);
                    volumeNodeGraph.Branch(i).AddRange(graphDT.Branch(i));
                }


                #endregion

                #region 计算矩阵 P(outer)
                Matrix P_outer = new Matrix(outerNodeCount, 2);
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (!graph.GraphNodes[i].IsInner)
                    {
                        P_outer[outerNodeIndexList.IndexOf(i), 0] = graph.GraphNodes[i].NodeVertex.X;
                        P_outer[outerNodeIndexList.IndexOf(i), 1] = graph.GraphNodes[i].NodeVertex.Y;
                    }
                }
                #endregion

                #region 计算inner，outer相关矩阵
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

                Matrix inner_InnerLaplacianM = LaplacianMatrix(inner_innerM, innerNodeDegreeList);

                bool flag = inner_InnerLaplacianM.Invert(0.0);
                Matrix inverse_Inner_InnerLaplacianM = inner_InnerLaplacianM;

                #endregion



                #region 求解矩阵 P(inner)
                // 求解矩阵 P(inner)
                Matrix P_inner = new Matrix(innerNodeCount, 2);
                P_inner = inverse_Inner_InnerLaplacianM * inner_outerM * P_outer;
                // P_inner.Scale(-1.0);
                #endregion


                #region 从P_inner矩阵中生成新的inner点的坐标
                List<Point3d> newInnerPoints = new List<Point3d>();
                Point3d innerPoint;
                for (int i = 0; i < innerNodeCount; i++)
                {
                    innerPoint = new Point3d(P_inner[i, 0], P_inner[i, 1], 0.0);
                    newInnerPoints.Add(innerPoint);
                }
                #endregion

                #region Relocate重定位，形成Graph的Node

                #region 新的NodePoint列表
                // 将新的inner点与原来的outer点合并成一个新的列表
                List<Point3d> newNodePoints = new List<Point3d>();
                newNodePoints.AddRange(newInnerPoints);
                List<Point3d> sortedBoundaryNodePoints = new List<Point3d>();
                for (int i = 0; i < outerNodeCount; i++)
                {
                    sortedBoundaryNodePoints.Add(graphNodes[innerNodeCount + i].NodeVertex);
                }
                newNodePoints.AddRange(sortedBoundaryNodePoints);
                #endregion

                Point3d currentCenter = new Point3d();
                for (int i = 0; i < newNodePoints.Count; i++)
                {
                    currentCenter += newNodePoints[i];
                }
                currentCenter /= newNodePoints.Count;

                Plane currentPlane = Plane.WorldXY;
                currentPlane.Origin = currentCenter;

                // 将新的inner+outer列表重新定位在基于worldXY变量的另一个地方
                List<Point3d> newNodePointsRelocated = UtilityFunctions.Relocate(newNodePoints, currentPlane, targetPlane);
                for (int i = 0; i < graphNodes.Count; i++)
                {
                    graphNodes[i].NodeVertex = newNodePointsRelocated[i];
                }

                #endregion

                Graph newGraph = new Graph(graphNodes, graphLoL);
                DA.SetData("Graph", newGraph);


                #region 绘制Graph的edge
                // Graph的所有Edge，在新的位置上画Edge
                List<Line> graphEdge = GraphEdgeList(graphDT, newNodePointsRelocated);
                DA.SetDataList("GraphEdge", graphEdge);
                #endregion

                #region 绘制Graph形成的Convex
                // Graph形成的每个convex，在新的位置上画Convex
                List<Point3d> newOuterNodePointsRelocated = new List<Point3d>();
                //for (int i = innerNodeCount; i < newNodePointsRelocated.Count; i++)
                //{
                //    sortedBoundaryNodePointsRelocated.Add(newNodePointsRelocated[i]);
                //}
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (!graph.GraphNodes[i].IsInner)
                    {
                        newOuterNodePointsRelocated.Add(newNodePointsRelocated[i]);
                    }
                }

                List<Polyline> allSplitFaceBoundarys = GMeshFaceBoundaries(newOuterNodePointsRelocated, graphEdge, targetPlane);
                #endregion
                DA.SetDataList("ConvexFaceBorders", allSplitFaceBoundarys);



                /* 这里要把inner_InnerAdjacencyDataTree找回来*/



                #region 输出subGraph
                // subGraph的每一支表示一个volumeNode与哪些其他的volumeNode相连
                List<List<int>> subGraphLoL = new List<List<int>>();
                for (int i = 0; i < inner_innerM.RowCount; i++)
                {
                    subGraphLoL.Add(new List<int>());
                    for (int j = 0; j < inner_innerM.ColumnCount; j++)
                    {
                        if (inner_innerM[i,j] == 1)
                        {
                            subGraphLoL[i].Add(j);
                        }
                    }
                }
                #endregion

                #region 输出subGraphNode
                List<Node> subGraphNodes = new List<Node>();
                for (int i = 0; i < graph.GraphNodes.Count; i++)
                {
                    if (graph.GraphNodes[i].IsInner)
                    {
                        subGraphNodes.Add(graph.GraphNodes[i]);
                    }
                }
                #endregion

                Graph subGraph = new Graph(subGraphNodes, subGraphLoL);
                DA.SetData("SubGraph", subGraph);
            }
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

        public List<Line> GraphEdgeList(DataTree<int> graph, List<Point3d> graphVertices)
        {
            List<Line> list = new List<Line>();

            for (int i = 0; i < graph.BranchCount; i++)
            {
                for (int j = 0; j < graph.Branch(i).Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graph.Branch(i)[j]]);
                    list.Add(item);
                }
            }

            return list;
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
                        LaplacianM[i, j] = 0 - M[i,j];
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