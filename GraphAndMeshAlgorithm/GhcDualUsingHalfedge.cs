﻿using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDualUsingHalfedge : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDualUsingHalfedge class.
        /// </summary>
        public GhcDualUsingHalfedge()
          : base("DualUsingHalfedge", "DualUsingHalfedge",
              "生成对应的对偶图（基于半边数据结构）",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            DualVertices = new List<Point3d>();
            DualPolylines = new List<List<Point3d>>();
            DualVertexTextDots = new List<TextDot>();

            NodePoints = new List<Point3d>();
            NodePointIndex = new List<int>();
            GraphEdges = new List<Line>();
            NodeTextDots = new List<TextDot>();
        }


        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        private List<Point3d> NodePoints;
        private List<int> NodePointIndex;
        private List<Line> GraphEdges;
        private List<TextDot> NodeTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果。", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // pManager.AddPointParameter("DualVertices", "DV", "所有构成对偶多边形的顶点", GH_ParamAccess.list);
            // pManager.AddTextParameter("DualFacesDescriptions", "DFD", "生成的对偶面的描述，包括对于第几个innerNode，周围由哪几个三角形的中心点包围形成", GH_ParamAccess.list);
            // pManager.AddCurveParameter("DualConvexPolygons", "DFB", "生成的对偶多边形", GH_ParamAccess.list);

            // pManager.AddIntegerParameter("DualConvexConnectivityGraph", "DCCGraph", "A graph representation as an adjacency list datatree", GH_ParamAccess.list);
            // pManager.AddPointParameter("DualConvexCenterPointsAsGraphNodes", "DCGraphNode", "Graph vertices as [edge] centrpoids of cells", GH_ParamAccess.list);

            pManager.AddGenericParameter("DualHalfedgeMesh", "DHM", "生成的对偶图（半边数据结构）", GH_ParamAccess.item);

            pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);

            pManager.AddIntegerParameter("faceIndexsFromOuterNodes", "", "", GH_ParamAccess.tree);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int volumeNodeCount = globalParameter.VolumeNodeCount;
            int boundaryNodeCount = globalParameter.BoundaryNodeCount;
            #endregion

            #region 局部变量初始化
            Thickness = 2;

            PlanktonMesh P = new PlanktonMesh();

            List<Node> nodes = new List<Node>();
            #endregion

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh",ref P))
            {
                // 获取节点
                DA.GetDataList<Node>("GraphNode", nodes);

                #region 获得对偶图
                // 利用半边数据结构求出对偶
                PlanktonMesh D = P.Dual();
                DA.SetData("DualHalfedgeMesh", D);
                #endregion

                #region DebugPrint

                #region HalfedgeMesh的顶点数据
                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(D);
                DA.SetDataList("DebugVerticesOutput", printVertices);
                #endregion

                #region HalfedgeMesh的半边数据
                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(D);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);
                #endregion

                #region HalfedgeMesh的每个面由哪些顶点构成
                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFacesVertices(D);
                DA.SetDataList("DebugFacesOutput", printFaces);
                #endregion

                #region HalfedgeMesh的每个面由哪些半边构成
                List<string> printFacesHalfedge = new List<string>();
                printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(D);
                DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);
                #endregion

                #endregion

                #region 得到所有的innerNode的序号和outerNode的序号
                List<int> innerNodeIndexs = new List<int>();
                List<int> outerNodeIndexs = new List<int>();
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].IsInner)
                    {
                        innerNodeIndexs.Add(i);
                    }
                    else
                    {
                        outerNodeIndexs.Add(i);
                    }
                }
                #endregion

                #region 找到每个outerNode相关的面和发出的半边
                // 对于每个outerNode找到由它发出的半边
                List<List<int>> halfedgeIndexsStartFromOuterNodes = new List<List<int>>();
                // 对于每个outerNode找到与它相关的面
                List<List<int>> faceIndexsAroundOuterNodes = new List<List<int>>();

                for (int i = 0; i < outerNodeIndexs.Count; i++)
                {
                    // 找到每个outerNode所发出的halfedge的index
                    int[] halfedgeIndexsStartFromOuterNode = P.Vertices.GetHalfedges(outerNodeIndexs[i]);
                    halfedgeIndexsStartFromOuterNodes.Add(new List<int>());
                    halfedgeIndexsStartFromOuterNodes[i].AddRange(halfedgeIndexsStartFromOuterNode);
                    // 找到每个outerNode所邻接的Face的index，-1表示邻接外界
                    int[] faceIndexsAroundOuterNode = P.Vertices.GetVertexFaces(outerNodeIndexs[i]);
                    faceIndexsAroundOuterNodes.Add(new List<int>());
                    faceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNode);
                }
                #endregion
                DA.SetDataTree(5, UtilityFunctions.LoLToDataTree<int>(faceIndexsAroundOuterNodes));

                #region 将每个outerNode相关的面的index的列表，分割成两个一组两个一组
                List<List<int[]>> pairFaceIndexsAroundOuterNodes = new List<List<int[]>>();
                /* 对于每个outerNode，找到跟它相关的面中，每两个相邻面的index，这个面的index就是Dual中顶点的index
                 * int[2]中表示围绕顶点顺时针排列的两个相邻面的index，如[0][1],[1][2],[2][3]...
                 */
                for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
                {
                    pairFaceIndexsAroundOuterNodes.Add(new List<int[]>());
                    for (int j = 0; j < faceIndexsAroundOuterNodes[i].Count - 1; j++)
                    {
                        int[] pairFaceIndexs = new int[2] { faceIndexsAroundOuterNodes[i][j], faceIndexsAroundOuterNodes[i][j + 1] };
                        pairFaceIndexsAroundOuterNodes[i].Add(pairFaceIndexs);
                    }
                }
                #endregion

                #region 可视化部分

                #region clear
                // 对偶图
                DualVertices.Clear();
                DualPolylines.Clear();
                DualVertexTextDots.Clear();
                // 原来的图
                NodePoints.Clear();
                NodePointIndex.Clear();
                GraphEdges.Clear();
                NodeTextDots.Clear();
                #endregion

                #region 对偶图
                DualVertices = PlanktonGh.RhinoSupport.GetPositions(D).ToList();

                Polyline[] dualPolylines = PlanktonGh.RhinoSupport.ToPolylines(D);
                for (int i = 0; i < dualPolylines.Length; i++)
                {
                    DualPolylines.Add(new List<Point3d>());
                    DualPolylines[i].AddRange(dualPolylines[i]);
                }

                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    string arg = string.Join<int>(";", D.Vertices.GetVertexFaces(i));
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), DualVertices[i]);
                    DualVertexTextDots.Add(textDot);
                }
                #endregion

                #region 原来的图

                #region 找到对偶图中的面所对应的InnerNode的序号，即NodePointIndex
                List<List<int>> pFaceIndexAroundInnerNode = new List<List<int>>();
                List<int> correspondingInnerNodeIndex = new List<int>();
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].IsInner)
                    {
                        //pFaceIndexAroundInnerNode.Add(new List<int>());
                        //pFaceIndexAroundInnerNode[i].AddRange(P.Vertices.GetVertexFaces(i));

                        pFaceIndexAroundInnerNode.Add(P.Vertices.GetVertexFaces(i).ToList<int>());
                        correspondingInnerNodeIndex.Add(i);
                    }
                    
                }

                for (int i = 0; i < D.Faces.Count; i++)
                {
                    for (int j = 0; j < pFaceIndexAroundInnerNode.Count; j++)
                    {
                        /* D.Faces.GetFaceVertices(i):对偶图中，每个面的顶点集
                         * pFaceIndexAroundInnerNode[j]:原图中，每个inner顶点所邻接的面的index集
                         * 如果这两个集没有差集（即完全相同的时候），那么把这个inner顶点对应的序号给对偶图中对应的面
                         */
                        if (D.Faces.GetFaceVertices(i).Except(pFaceIndexAroundInnerNode[j]).ToArray().Length == 0)
                        {
                            NodePointIndex.Add(correspondingInnerNodeIndex[j]);
                        }
                    }
                }
                #endregion

                #region 添加对偶图所有面的中心点作为NodePoint的坐标位置
                for (int i = 0; i < D.Faces.Count; i++)
                {
                    NodePoints.Add(PlanktonGh.RhinoSupport.ToPoint3d(D.Faces.GetFaceCenter(i)));
                }
                #endregion

                #region 根据innerNode的序号，写入对应的NodeTextDots
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].IsInner)
                    {
                        NodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, nodes[i].NodeAttribute.NodeLabel), this.NodePoints[NodePointIndex.IndexOf(i)]));
                    }
                }
                #endregion

                #region 代表innerNode的点之间连线
                List<List<int>> faceAdjacency = UtilityFunctions.GetAdjacencyFaceIndexs(D);
                for (int i = 0; i < NodePoints.Count; i++)
                {
                    for (int j = 0; j < faceAdjacency[i].Count; j++)
                    {
                        GraphEdges.Add(new Line(NodePoints[i], NodePoints[faceAdjacency[i][j]]));
                    }
                }
                #endregion
                #endregion

                #endregion

            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);

            for (int i = 0; i < DualPolylines.Count; i++)
            {
                args.Display.DrawPolyline(DualPolylines[i], Color.DarkOrange, Thickness);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < NodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(NodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);

            for (int i = 0; i < DualPolylines.Count; i++)
            {
                args.Display.DrawPolyline(DualPolylines[i], Color.DarkOrange, Thickness);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < NodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(NodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
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
            get { return new Guid("f73e455b-814f-4409-84da-b82a6df8f8f6"); }
        }
    }
}