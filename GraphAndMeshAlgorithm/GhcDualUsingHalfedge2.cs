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

        public enum ShowMode { ShowTopology, NotShowTopology };/* 定义一个enum类型 */
        public ShowMode CompWorkMode { get; set; } = ShowMode.ShowTopology;/* 使用这个enum类型来定义一个代表电池工作状态的变量 */


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

                DualGraphWithHM dualGraphWithHM = new DualGraphWithHM(D,
                                                                      pGraphNodes,
                                                                      pGraphTables,
                                                                      faceIndexsAroundOuterNodes);

                List<int> faceCorrespondingGraphNodeIndex = new List<int>();
                for (int i = 0; i < pGraphNodes.Count; i++)
                {
                    if (pGraphNodes[i].IsInner)
                    {
                        faceCorrespondingGraphNodeIndex.Add(i);
                    }
                }
                dualGraphWithHM.FaceCorrespondingGraphNodeIndex = faceCorrespondingGraphNodeIndex;

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

            if (this.CompWorkMode == ShowMode.ShowTopology)
            {
                args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

                for (int i = 0; i < GraphNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(GraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
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

            if (this.CompWorkMode == ShowMode.ShowTopology)
            {
                args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

                for (int i = 0; i < GraphNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(GraphNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
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

        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new ShowAttribute(this);
        }

        public class ShowAttribute : GH_ComponentAttributes
        {
            public ShowAttribute(GhcDualUsingHalfedge2 component) : base(component) { }

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

                    using (GH_Capsule capsule = GH_Capsule.CreateCapsule(buttonRect, GH_Palette.Black))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(((GhcDualUsingHalfedge2)Owner).CompWorkMode.ToString(),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.White,
                                        buttonRect,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });
                }
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                RectangleF buttonRect = /* -重新计算按钮的区域大小- */ new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);

                if (e.Button == MouseButtons.Left && buttonRect.Contains(e.CanvasLocation))
                {
                    GhcDualUsingHalfedge2 comp = (GhcDualUsingHalfedge2)Owner; /* 通过Owner属性来获得电池本身 */

                    /* 依照电池当前工作状态来改变电池 */
                    if (comp.CompWorkMode == GhcDualUsingHalfedge2.ShowMode.NotShowTopology)
                        comp.CompWorkMode = GhcDualUsingHalfedge2.ShowMode.ShowTopology;
                    else
                        comp.CompWorkMode = GhcDualUsingHalfedge2.ShowMode.NotShowTopology;

                    /* 改变完电池后，重启计算 */
                    comp.ExpireSolution(true);

                    /* 结束鼠标事件处理，通知GH已经处理完毕 */
                    return GH_ObjectResponse.Handled;
                }


                return GH_ObjectResponse.Ignore;
            }
        }
    }
}