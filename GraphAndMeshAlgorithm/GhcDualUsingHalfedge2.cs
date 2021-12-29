using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Attributes;
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
using System.Windows.Forms;
using VolumeGeneratorBasedOnGraph.Class;
using Grasshopper.GUI.Canvas;
using Grasshopper.GUI;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDualUsingHalfedge2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDualUsingHalfedge2 class.
        /// </summary>
        public GhcDualUsingHalfedge2()
          : base("DualUsingHalfedge2", "DualUsingHalfedge2",
              "生成对应的对偶图（基于半边数据结构）",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            DualVertices = new List<Point3d>();
            DualPolylines = new List<List<Point3d>>();
            DualVertexTextDots = new List<TextDot>();

            GraphNodePoints = new List<Point3d>();
            GraphEdges = new List<Line>();
            GraphNodeTextDots = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        private List<Point3d> GraphNodePoints;
        private List<Line> GraphEdges;
        private List<TextDot> GraphNodeTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GraphWithHM", "GHM", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // 0
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GraphWithHM graphWithHM = new GraphWithHM();

            if (DA.GetData<GraphWithHM>("GraphWithHM", ref graphWithHM))
            {
                Thickness = 2;
                
                // 获取相关属性
                List<GraphNode> pGraphNodes = graphWithHM.GraphNodes;
                List<List<int>> pGraphTables = graphWithHM.GraphTables;
                PlanktonMesh pGraphWithHM = graphWithHM.PlanktonMesh;


                #region 获得对偶图
                // 利用半边数据结构求出对偶
                PlanktonMesh D = pGraphWithHM.Dual();
                #endregion

                #region 得到所有的innerNode的序号和outerNode的序号
                List<int> innerNodeIndexs = graphWithHM.InnerNodeIndexList;
                List<int> outerNodeIndexs = graphWithHM.OuterNodeIndexList;
                int innerNodeCount = graphWithHM.InnerNodeCount;
                int outerNodeCount = graphWithHM.OuterNodeCount;
                #endregion

                #region 计算faceIndexsAroundOuterNodes，找到每个outerNode相关的面和发出的半边
                // 对于每个outerNode找到由它发出的半边
                List<List<int>> halfedgeIndexsStartFromOuterNodes = new List<List<int>>();
                // 对于每个outerNode找到与它相关的面
                List<List<int>> faceIndexsAroundOuterNodes = new List<List<int>>();

                for (int i = 0; i < outerNodeIndexs.Count; i++)
                {
                    // 找到每个outerNode所发出的halfedge的index
                    int[] halfedgeIndexsStartFromOuterNode = pGraphWithHM.Vertices.GetHalfedges(outerNodeIndexs[i]);
                    halfedgeIndexsStartFromOuterNodes.Add(new List<int>());
                    halfedgeIndexsStartFromOuterNodes[i].AddRange(halfedgeIndexsStartFromOuterNode);
                    // 找到每个outerNode所邻接的Face的index，-1表示邻接外界
                    int[] faceIndexsAroundOuterNode = pGraphWithHM.Vertices.GetVertexFaces(outerNodeIndexs[i]);
                    faceIndexsAroundOuterNodes.Add(new List<int>());
                    faceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNode);
                }
                #endregion

                //#region 将每个outerNode相关的面的index的列表，分割成两个一组两个一组
                //List<List<int[]>> pairFaceIndexsAroundOuterNodes = new List<List<int[]>>();
                ///* 对于每个outerNode，找到跟它相关的面中，每两个相邻面的index，这个面的index就是Dual中顶点的index
                // * int[2]中表示围绕顶点顺时针排列的两个相邻面的index，如[0][1],[1][2],[2][3]...
                // */
                //for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
                //{
                //    pairFaceIndexsAroundOuterNodes.Add(new List<int[]>());
                //    for (int j = 0; j < faceIndexsAroundOuterNodes[i].Count - 1; j++)
                //    {
                //        int[] pairFaceIndexs = new int[2] { faceIndexsAroundOuterNodes[i][j], faceIndexsAroundOuterNodes[i][j + 1] };
                //        pairFaceIndexsAroundOuterNodes[i].Add(pairFaceIndexs);
                //    }
                //}
                //#endregion

                //SortedDictionary<int, int> inner_dFace = new SortedDictionary<int, int>();
                //#region 找到对偶图中的面所对应的InnerNode的序号，即NodePointIndex，NodePointIndex实际上是对偶图每个面所对应的原来图中的InnerNode的序号
                //List<List<int>> pFaceIndexAroundInnerNode = new List<List<int>>();
                //List<int> correspondingInnerNodeIndex = new List<int>();
                //for (int i = 0; i < decomposedPGraphNodes.Count; i++)
                //{
                //    if (decomposedPGraphNodes[i].IsInner)
                //    {
                //        pFaceIndexAroundInnerNode.Add(pGraphWithHM.Vertices.GetVertexFaces(i).ToList<int>());
                //        correspondingInnerNodeIndex.Add(i);
                //    }
                //}

                //for (int i = 0; i < D.Faces.Count; i++)
                //{
                //    for (int j = 0; j < pFaceIndexAroundInnerNode.Count; j++)
                //    {
                //        /* D.Faces.GetFaceVertices(i):对偶图中，每个面的顶点集
                //         * pFaceIndexAroundInnerNode[j]:原图中，每个inner顶点所邻接的面的index集
                //         * 如果这两个集没有差集（即完全相同的时候），那么把这个inner顶点对应的序号给对偶图中对应的面
                //         */
                //        if (D.Faces.GetFaceVertices(i).Except(pFaceIndexAroundInnerNode[j]).ToArray().Length == 0)
                //        {
                //            inner_dFace.Add(correspondingInnerNodeIndex[j], i);
                //        }
                //    }
                //}

                //#endregion

                DualGraphWithHM dualGraphWithHM = new DualGraphWithHM(D,
                                                                      pGraphNodes,
                                                                      pGraphTables,
                                                                      faceIndexsAroundOuterNodes);
                

                DA.SetData("DualGraphWithHM", dualGraphWithHM);

                #region 可视化部分
                PlanktonMesh MeshForVisualize = new PlanktonMesh();

                MeshForVisualize = dualGraphWithHM.DualPlanktonMesh;

                #region clear
                // 对偶图
                DualVertices.Clear();
                DualPolylines.Clear();
                DualVertexTextDots.Clear();
                // 原来的图
                GraphNodePoints.Clear();
                GraphEdges.Clear();
                GraphNodeTextDots.Clear();
                #endregion

                #region 对偶图
                DualVertices = PlanktonGh.RhinoSupport.GetPositions(MeshForVisualize).ToList();

                Polyline[] dualPolylines = PlanktonGh.RhinoSupport.ToPolylines(MeshForVisualize);
                for (int i = 0; i < dualPolylines.Length; i++)
                {
                    DualPolylines.Add(new List<Point3d>());
                    DualPolylines[i].AddRange(dualPolylines[i]);
                }

                for (int i = 0; i < MeshForVisualize.Vertices.Count; i++)
                {
                    string arg = string.Join<int>(";", D.Vertices.GetVertexFaces(i));
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), DualVertices[i]);
                    DualVertexTextDots.Add(textDot);
                }
                #endregion

                #region 原来的图
                #region 添加对偶图所有面的中心点作为GraphNodePoints的坐标位置
                for (int i = 0; i < MeshForVisualize.Faces.Count; i++)
                {
                    GraphNodePoints.Add(PlanktonGh.RhinoSupport.ToPoint3d(MeshForVisualize.Faces.GetFaceCenter(i)));
                }
                #endregion

                #region 根据innerNode的序号，写入对应的GraphNodeTextDots
                List<GraphNode> graphNodes = new List<GraphNode>();
                for (int i = 0; i < pGraphNodes.Count; i++)
                {
                    if (pGraphNodes[i].IsInner)
                    {
                        graphNodes.Add(pGraphNodes[i]);
                    }
                }
                for (int i = 0; i < graphNodes.Count; i++)
                {
                    GraphNodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, graphNodes[i].NodeAttribute.NodeLabel), GraphNodePoints[i]));
                }
                #endregion

                #region 代表innernode之间的连线
                List<List<int>> faceAdjacency = UtilityFunctions.GetAdjacencyFaceIndexs(MeshForVisualize);
                for (int i = 0; i < GraphNodePoints.Count; i++)
                {
                    for (int j = 0; j < faceAdjacency[i].Count; j++)
                    {
                        GraphEdges.Add(new Line(GraphNodePoints[i], GraphNodePoints[faceAdjacency[i][j]]));
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

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < GraphNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(GraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // base.DrawViewportMeshes(args);
            for (int i = 0; i < DualPolylines.Count; i++)
            {
                args.Display.DrawPolyline(DualPolylines[i], Color.DarkOrange, Thickness);
            }

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < GraphNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(GraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
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
            get { return new Guid("27b89768-b836-4b2d-aae8-39a6020ca8ae"); }
        }
    }
}