﻿using Grasshopper;
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

            OriginGraphNodePoints = new List<Point3d>();
            DualFaceIndexCorrespondingInnerNodeIndex = new List<int>();
            OriginGraphEdges = new List<Line>();
            OriginGraphNodeTextDots = new List<TextDot>();
        }


        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        private List<Point3d> OriginGraphNodePoints;
        private List<int> DualFaceIndexCorrespondingInnerNodeIndex;
        private List<int> DualFaceIndexCorrespondingVolumeNodeIndex;
        private List<Line> OriginGraphEdges;
        private List<TextDot> OriginGraphNodeTextDots;

        public enum CalMode { Daul, IntegrateDual}
        public CalMode CompWorkMode { get; set; } = CalMode.Daul;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // 0
            pManager.AddGenericParameter("GraphWithHM", "GHM", "图结构", GH_ParamAccess.item);
            // 1
            // pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果。", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // 0
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
            // 1
            pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);
            // 2
            // pManager.AddIntegerParameter("faceIndexsFromOuterNodes", "", "", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 局部变量初始化
            Thickness = 2;

            GraphWithHM graphWithHM = new GraphWithHM();
            
            #endregion

            if (DA.GetData<GraphWithHM>("GraphWithHM", ref graphWithHM))
            {
                if (graphWithHM.VolumeContainsWhichInnerNode == null)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "当前传入的GraphWithHM是来自一个空的叶子节点");
                    return;
                }
                
                // 获取相关属性
                List<GraphNode> originGraphNodes = graphWithHM.GraphNodes;
                List<List<int>> originGraphTables = graphWithHM.GraphTables;
                PlanktonMesh pGraphWithHM = graphWithHM.PlanktonMesh;

                List<int> volumeNodesIndex = graphWithHM.VolumeNodesIndex;
                List<GraphNode> volumeNodes = graphWithHM.VolumeNodes;
                Dictionary<int, List<int>> volumeContainsWhichInnerNode = graphWithHM.VolumeContainsWhichInnerNode;

                #region 获得对偶图
                // 利用半边数据结构求出对偶
                PlanktonMesh D = pGraphWithHM.Dual();
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
                // DA.SetDataTree(5, UtilityFunctions.LoLToDataTree<int>(faceIndexsAroundOuterNodes));





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


                #region 计算FaceCorresponding相关
                List<int> dualFaceCorrespondingInnerNodeIndex = new List<int>();
                #region 找到对偶图中的面所对应的InnerNode的序号，即NodePointIndex，NodePointIndex实际上是对偶图每个面所对应的原来图中的InnerNode的序号
                List<List<int>> pFaceIndexAroundInnerNode = new List<List<int>>();
                List<int> correspondingInnerNodeIndex = new List<int>();
                for (int i = 0; i < originGraphNodes.Count; i++)
                {
                    if (originGraphNodes[i].IsInner)
                    {
                        //pFaceIndexAroundInnerNode.Add(new List<int>());
                        //pFaceIndexAroundInnerNode[i].AddRange(P.Vertices.GetVertexFaces(i));

                        pFaceIndexAroundInnerNode.Add(pGraphWithHM.Vertices.GetVertexFaces(i).ToList<int>());
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
                            dualFaceCorrespondingInnerNodeIndex.Add(correspondingInnerNodeIndex[j]);
                        }
                    }
                }

                List<GraphNode> dualFaceCorrespondingInnerNodes = new List<GraphNode>();
                for (int i = 0; i < dualFaceCorrespondingInnerNodeIndex.Count; i++)
                {
                    dualFaceCorrespondingInnerNodes.Add(originGraphNodes[dualFaceCorrespondingInnerNodeIndex[i]]);
                }


                DualFaceIndexCorrespondingInnerNodeIndex.Clear();
                DualFaceIndexCorrespondingInnerNodeIndex = dualFaceCorrespondingInnerNodeIndex;

                #endregion
                #endregion

                




                #region 找到对偶图中每个顶点属于哪个InnerNode
                // 找到对偶图中每个顶点邻接的Face的Index，要去掉-1
                List<List<int>> dFaceIndexAroundVertice = new List<List<int>>();
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    dFaceIndexAroundVertice.Add(new List<int>());
                    List<int> currentDFaceIndex = D.Vertices.GetVertexFaces(i).ToList<int>();
                    // 去掉-1
                    currentDFaceIndex.RemoveAt(0);
                    dFaceIndexAroundVertice[i].AddRange(currentDFaceIndex);
                }

                // 对于每个对偶图上的vertex，它属于哪个InnerNode
                List<List<int>> dVertexBelongsToWhichInnerNode = new List<List<int>>();
                for (int i = 0; i < D.Vertices.Count; i++)
                {
                    dVertexBelongsToWhichInnerNode.Add(new List<int>());
                    foreach (int dFaceIndex in dFaceIndexAroundVertice[i])
                    {
                        // NodePointIndex实际上是对偶图每个面所对应的原来图中的InnerNode的序号
                        dVertexBelongsToWhichInnerNode[i].Add(DualFaceIndexCorrespondingInnerNodeIndex[dFaceIndex]);
                    }
                }
                #endregion


                List<List<int>> innerNodeContainsWhichDVertex = new List<List<int>>();
                for (int i = 0; i < graphWithHM.InnerNodeIndexList.Count; i++)
                {
                    innerNodeContainsWhichDVertex.Add(new List<int>());
                }




                #region 找到对偶图中每个顶点属于哪个Volume
                Dictionary<int, List<int>> volumeCorrespondingDualVertices = new Dictionary<int, List<int>>();
                foreach (KeyValuePair<int, List<int>> pair in volumeContainsWhichInnerNode)
                {
                    List<int> dualVertices = new List<int>();
                    for (int i = 0; i < pair.Value.Count; i++)
                    {
                        for (int j = 0; j < dVertexBelongsToWhichInnerNode[pair.Value[i]].Count; j++)
                        {
                            if (!dualVertices.Contains(dVertexBelongsToWhichInnerNode[pair.Value[i]][j]))
                            {
                                dualVertices.Add(dVertexBelongsToWhichInnerNode[pair.Value[i]][j]);
                            }
                        }
                    }
                    volumeCorrespondingDualVertices.Add(pair.Key, dualVertices);
                }





                #endregion

                


                DualGraphWithHM dualGraphWithHM = new DualGraphWithHM(D, 
                                                                      originGraphNodes, 
                                                                      originGraphTables, 
                                                                      dualFaceCorrespondingInnerNodeIndex,
                                                                      volumeNodesIndex,
                                                                      volumeNodes,
                                                                      volumeContainsWhichInnerNode,
                                                                      faceIndexsAroundOuterNodes);


                DA.SetData("DualGraphWithHM", dualGraphWithHM);


                #region 可视化部分
                PlanktonMesh MeshForVisualize = new PlanktonMesh();

                if (CompWorkMode == CalMode.Daul)
                {
                    MeshForVisualize = dualGraphWithHM.PlanktonMesh;

                    #region clear
                    // 对偶图
                    DualVertices.Clear();
                    DualPolylines.Clear();
                    DualVertexTextDots.Clear();
                    // 原来的图
                    OriginGraphNodePoints.Clear();
                    OriginGraphEdges.Clear();
                    OriginGraphNodeTextDots.Clear();
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
                        // 这里不能是pMeshForVisualize.Vertices.GetVertexFaces(i)，这样才能正确显示原来图中的InnerNode的序号
                        // string arg = string.Join<int>(";", D.Vertices.GetVertexFaces(i));
                        string arg = string.Join<int>(";", dualGraphWithHM.DVertexBelongsToWhichInnerNode()[i]);
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), DualVertices[i]);
                        DualVertexTextDots.Add(textDot);
                    }
                    #endregion

                    #region 原来的图

                    #region 添加对偶图所有面的中心点作为NodePoint的坐标位置
                    for (int i = 0; i < MeshForVisualize.Faces.Count; i++)
                    {
                        OriginGraphNodePoints.Add(PlanktonGh.RhinoSupport.ToPoint3d(MeshForVisualize.Faces.GetFaceCenter(i)));
                    }
                    #endregion

                    #region 根据innerNode的序号，写入对应的NodeTextDots
                    for (int i = 0; i < originGraphNodes.Count; i++)
                    {
                        if (originGraphNodes[i].IsInner)
                        {
                            OriginGraphNodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, originGraphNodes[i].NodeAttribute.NodeLabel), this.OriginGraphNodePoints[DualFaceIndexCorrespondingInnerNodeIndex.IndexOf(i)]));
                        }
                    }
                    #endregion

                    #region 代表innerNode的点之间连线
                    List<List<int>> faceAdjacency = UtilityFunctions.GetAdjacencyFaceIndexs(MeshForVisualize);
                    for (int i = 0; i < OriginGraphNodePoints.Count; i++)
                    {
                        for (int j = 0; j < faceAdjacency[i].Count; j++)
                        {
                            OriginGraphEdges.Add(new Line(OriginGraphNodePoints[i], OriginGraphNodePoints[faceAdjacency[i][j]]));
                        }
                    }
                    #endregion
                    #endregion
                }
                else
                {

                    MeshForVisualize = dualGraphWithHM.IntegratePlanktonMesh;

                    DualFaceIndexCorrespondingVolumeNodeIndex = dualGraphWithHM.DFaceBelongsToWhichVolumeNode();

                    #region clear
                    // 对偶图
                    DualVertices.Clear();
                    DualPolylines.Clear();
                    DualVertexTextDots.Clear();
                    // 原来的图
                    OriginGraphNodePoints.Clear();
                    OriginGraphEdges.Clear();
                    OriginGraphNodeTextDots.Clear();
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
                        // 这里不能是pMeshForVisualize.Vertices.GetVertexFaces(i)，这样才能正确显示原来图中的InnerNode的序号
                        // string arg = string.Join<int>(";", D.Vertices.GetVertexFaces(i));
                        string arg = string.Join<int>(";", dualGraphWithHM.DVertexBelongsToWhichVolumeNode()[i]);
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), DualVertices[i]);
                        DualVertexTextDots.Add(textDot);
                    }
                    #endregion

                    #region 原来的图

                    //#region 添加对偶图所有面的中心点作为NodePoint的坐标位置
                    //for (int i = 0; i < MeshForVisualize.Faces.Count; i++)
                    //{
                    //    OriginGraphNodePoints.Add(PlanktonGh.RhinoSupport.ToPoint3d(MeshForVisualize.Faces.GetFaceCenter(i)));
                    //}
                    //#endregion

                    //#region 根据innerNode的序号，写入对应的NodeTextDots
                    //for (int i = 0; i < originGraphNodes.Count; i++)
                    //{
                    //    if (originGraphNodes[i].IsInner)
                    //    {
                    //        OriginGraphNodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, originGraphNodes[i].NodeAttribute.NodeLabel), this.OriginGraphNodePoints[DualFaceIndexCorrespondingInnerNodeIndex.IndexOf(i)]));
                    //    }
                    //}
                    //#endregion

                    //#region 代表innerNode的点之间连线
                    //List<List<int>> faceAdjacency = UtilityFunctions.GetAdjacencyFaceIndexs(MeshForVisualize);
                    //for (int i = 0; i < OriginGraphNodePoints.Count; i++)
                    //{
                    //    for (int j = 0; j < faceAdjacency[i].Count; j++)
                    //    {
                    //        OriginGraphEdges.Add(new Line(OriginGraphNodePoints[i], OriginGraphNodePoints[faceAdjacency[i][j]]));
                    //    }
                    //}
                    //#endregion
                    #endregion
                }



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

            args.Display.DrawLines(OriginGraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < OriginGraphNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(OriginGraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
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

            args.Display.DrawLines(OriginGraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < OriginGraphNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(OriginGraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
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

        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new GhcAttribute(this);
        }
    }

    public class GhcAttribute : GH_ComponentAttributes
    {
        public GhcAttribute(GhcDualUsingHalfedge component) : base(component) { }

        protected override void Layout()
        {
            base.Layout();
            /* 先执行base.Layout()，可以按GH电池默认方式计算电池的出/入口需要的高度，我们在下面基于这个高度进行更改 */
            Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 20.0f);
        }

        protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
        {
            base.Render(canvas, graphics, channel);/* 执行基本的电池渲染 */

            /* 额外的电池渲染，仅在“Objects”这个渲染轨道绘制 */
            if (channel == GH_CanvasChannel.Objects)
            {
                RectangleF buttonRect = /* 按钮的位置 */ new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);

                /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                buttonRect.Inflate(-2.0f, -2.0f);

                using (GH_Capsule capsule = GH_Capsule.CreateCapsule(buttonRect,GH_Palette.Black))
                {
                    /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                    /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                    /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                    capsule.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                }

                graphics.DrawString(((GhcDualUsingHalfedge)Owner).CompWorkMode.ToString(),
                                    new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                    Brushes.White,
                                    buttonRect,
                                    new StringFormat()
                                    {
                                        Alignment = StringAlignment.Center,
                                        LineAlignment = StringAlignment.Center,
                                    });
            }
        }

        public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
        {
            RectangleF buttonRect = /* -重新计算按钮的区域大小- */ new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);

            if (e.Button == MouseButtons.Left && buttonRect.Contains(e.CanvasLocation))
            {
                GhcDualUsingHalfedge comp = (GhcDualUsingHalfedge)Owner;/* 通过Owner属性来获得电池本身 */

                /* 依照电池当前工作状态来改变电池 */
                if (comp.CompWorkMode == GhcDualUsingHalfedge.CalMode.Daul)
                {
                    comp.CompWorkMode = GhcDualUsingHalfedge.CalMode.IntegrateDual;
                }
                else
                {
                    comp.CompWorkMode = GhcDualUsingHalfedge.CalMode.Daul;
                }

                /* 改变完电池后，重启计算 */
                comp.ExpireSolution(true);

                /* 结束鼠标事件处理，通知GH已经处理完毕 */
                return GH_ObjectResponse.Handled;
            }

            /* 若上述条件未满足，则直接返回“未处理” */
            return GH_ObjectResponse.Ignore;
        }
    }
}