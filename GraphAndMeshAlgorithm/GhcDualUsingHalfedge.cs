using Grasshopper;
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
            GraphEdges = new List<Line>();
            NodeTextDots = new List<TextDot>();
        }


        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        private List<Point3d> NodePoints;
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

            pManager.AddIntegerParameter("faceIndexsFromOuterNodes", "", "", GH_ParamAccess.tree);

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
            int volumeNodeCount = globalParameter.VolumeNodeCount;
            int boundaryNodeCount = globalParameter.BoundaryNodeCount;

            Thickness = 2;

            PlanktonMesh P = new PlanktonMesh();

            List<Node> nodes = new List<Node>();

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh",ref P))
            {
                // 获取节点
                DA.GetDataList<Node>("GraphNode", nodes);



                // 利用半边数据结构求出对偶
                PlanktonMesh D = P.Dual();
                DA.SetData("DualHalfedgeMesh", D);

                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(D);
                DA.SetDataList("DebugVerticesOutput", printVertices);

                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(D);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);

                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFacesVertices(D);
                DA.SetDataList("DebugFacesOutput", printFaces);


                //for (int i = 0; i < P.Faces.Count; i++)
                //{
                //    int[] pFaceHalfedges = P.Faces.GetHalfedges(i);

                //    for (int j = 0; j < pFaceHalfedges.Length; j++)
                //    {
                //        PlanktonHalfedge pFaceHalfedge = P.Halfedges[pFaceHalfedges[j]];
                //        PlanktonHalfedge pFacePairHalfedge = P.Halfedges[P.Halfedges.GetPairHalfedge(pFaceHalfedges[j])];
                //    }
                //}

                // 得到所有的innerNode的序号和outerNode的序号
                List<int> innerNodeIndexs = new List<int>();
                List<int> outerNodeIndexs = new List<int>();
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        innerNodeIndexs.Add(i);
                    }
                    else
                    {
                        outerNodeIndexs.Add(i);
                    }
                }

                // 对于每个outerNode找到由它发出的半边
                List<List<int>> halfedgeIndexsFromOuterNodes = new List<List<int>>();
                // 对于每个outerNode找到与它相关的面
                List<List<int>> faceIndexsFromOuterNodes = new List<List<int>>();

                for (int i = 0; i < outerNodeIndexs.Count; i++)
                {
                    int[] halfedgeIndexsFromOuterNode = P.Vertices.GetHalfedges(outerNodeIndexs[i]);
                    halfedgeIndexsFromOuterNodes.Add(new List<int>());
                    halfedgeIndexsFromOuterNodes[i].AddRange(halfedgeIndexsFromOuterNode);
                    int[] faceIndexsFromOuterNode = P.Vertices.GetVertexFaces(outerNodeIndexs[i]);
                    faceIndexsFromOuterNodes.Add(new List<int>());
                    faceIndexsFromOuterNodes[i].AddRange(faceIndexsFromOuterNode);
                }

                DA.SetDataTree(4, UtilityFunctions.LoLToDataTree<int>(faceIndexsFromOuterNodes));





                List<List<int[]>> pairFaceIndexsFromOuterNodes = new List<List<int[]>>();
                /* 对于每个outerNode，找到跟它相关的面中，每两个相邻面的index，这个面的index就是Dual中顶点的index
                 * int[2]中表示围绕顶点顺时针排列的两个相邻面的index，如[0][1],[1][2],[2][3]...
                 */
                for (int i = 0; i < faceIndexsFromOuterNodes.Count; i++)
                {
                    pairFaceIndexsFromOuterNodes.Add(new List<int[]>());
                    for (int j = 0; j < faceIndexsFromOuterNodes[i].Count - 1; j++)
                    {
                        int[] pairFaceIndexs = new int[2] { faceIndexsFromOuterNodes[i][j], faceIndexsFromOuterNodes[i][j + 1] };
                        pairFaceIndexsFromOuterNodes[i].Add(pairFaceIndexs);
                    }
                }


                // 绘图部分
                // 对偶图
                DualVertices.Clear();
                DualPolylines.Clear();
                DualVertexTextDots.Clear();

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

                // 原来的图
                NodePoints.Clear();
                GraphEdges.Clear();
                NodeTextDots.Clear();

                for (int i = 0; i < D.Faces.Count; i++)
                {
                    NodePoints.Add(PlanktonGh.RhinoSupport.ToPoint3d(D.Faces.GetFaceCenter(i)));
                }

                List<List<int>> faceAdjacency = UtilityFunctions.GetAdjacencyFaceIndexs(D);
                for (int i = 0; i < NodePoints.Count; i++)
                {
                    for (int j = 0; j < faceAdjacency[i].Count; j++)
                    {
                        GraphEdges.Add(new Line(NodePoints[i], NodePoints[faceAdjacency[i][j]]));
                    }
                }

                for (int i = 0; i < volumeNodeCount; i++)
                {
                    NodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, nodes[i].NodeAttribute.NodeLabel), this.NodePoints[i]));
                }
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