using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcLaplaceSmooth : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcLaplaceSmooth class.
        /// </summary>
        public GhcLaplaceSmooth()
          : base("LaplaceSmooth", "LaplaceSmooth",
              "对生成的半边网格进行拉普拉斯平滑",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            PRhinoMesh = new Mesh();
            PHalfedgeDottedCurves = new List<Curve>();

            GraphEdges = new List<Line>();

            InnerNodeTextDot = new List<TextDot>();
            OuterNodeTextDot = new List<TextDot>();
        }


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //0
            pManager.AddGenericParameter("PlanktonMesh", "PM", "输入半边网格", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);
            pManager.AddGenericParameter("GraphNode", "GN", "图结构中的节点", GH_ParamAccess.list);


            //4
            pManager.AddIntegerParameter("Flip", "Flip", "Criterion used to decide when to flip edges (0 for valence based, 1 for angle based)", GH_ParamAccess.item, 1);

            //5
            pManager.AddNumberParameter("PullStrength", "Pull", "Strength of pull to target geometry (between 0 and 1). Set to 0 for minimal surfaces", GH_ParamAccess.item, 0.8);

            //6
            pManager.AddIntegerParameter("Iterations", "Iter", "Number of steps between outputs", GH_ParamAccess.item, 1);

            //7
            pManager.AddBooleanParameter("Reset", "Reset", "True to initialize, false to run remeshing. Connect a timer for continuous remeshing", GH_ParamAccess.item, true);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "M", "Remeshed result as Plankton Mesh", GH_ParamAccess.item);
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



        private PlanktonMesh PDeepCopy = new PlanktonMesh();
        private Mesh M = new Mesh();
        private List<int> AnchorV = new List<int>();
        // private List<int> FeatureV = new List<int>();
        // private List<int> FeatureE = new List<int>();
        private bool initialized;

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Thickness = 2;

            // ITargetLength TargetLength = null;
            bool reset = false;
            int Flip = 0;
            // List<Curve> FC = new List<Curve>();
            List<Point3d> FV = new List<Point3d>();
            double FixT = 0.01;
            double PullStrength = 0.8;
            double SmoothStrength = 0.8;
            // double LengthTol = 0.15;
            bool Minim = false;
            int Iters = 1;

            List<GraphNode> nodes = new List<GraphNode>();

            PlanktonMesh p = new PlanktonMesh();


            GH_ObjectWrapper Geo = new GH_ObjectWrapper();
            DA.GetData<GH_ObjectWrapper>("PlanktonMesh", ref Geo);

            // GH_ObjectWrapper Obj = null;
            // DA.GetData<GH_ObjectWrapper>(1, ref Obj);
            // TargetLength = Obj.Value as ITargetLength;

            //DA.GetDataList<Curve>(2, FC);
            //DA.GetDataList<Point3d>(3, FV);
            DA.GetData<int>("Flip", ref Flip);
            DA.GetData<double>("PullStrength", ref PullStrength);
            DA.GetData<int>("Iterations", ref Iters);
            DA.GetData<bool>("Reset", ref reset);




            if (PullStrength == 0) { Minim = true; }

            if (Geo.Value is PlanktonMesh)
            {
                DA.GetData<PlanktonMesh>("PlanktonMesh", ref p);
                PDeepCopy = new PlanktonMesh(p);
                M = PDeepCopy.ToRhinoMesh();
            }

            DA.GetDataList<GraphNode>("GraphNode", nodes);
            List<GraphNode> nodesDeepCopy = new List<GraphNode>();
            for (int i = 0; i < nodes.Count; i++)
            {
                nodesDeepCopy.Add(new GraphNode(nodes[i]));
            }

            GH_Structure<GH_Integer> gh_Structure_graph = null;
            DataTree<int> graph = new DataTree<int>();
            List<List<int>> pGraphLoL = new List<List<int>>();
            DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_graph);
            // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>，再转化为LoL
            UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_graph, ref graph);
            pGraphLoL = UtilityFunctions.DataTreeToLoL<int>(graph);


            for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
            {
                if (!nodesDeepCopy[i].IsInner)
                {
                    FV.Add(nodesDeepCopy[i].NodeVertex);
                }
            }

            if (reset || initialized == false)
            {
                #region reset

                initialized = true;

                AnchorV.Clear();

                // 标记OuterNode为AnchorV
                for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
                {
                    Point3d Pt = p.Vertices[i].ToPoint3d();
                    AnchorV.Add(-1);
                    for (int j = 0; j < FV.Count; j++)
                    {
                        if (Pt.DistanceTo(FV[j]) < FixT)
                        {
                            AnchorV[AnchorV.Count - 1] = j;
                        }
                    }
                }
                #endregion
            }
            else
            {
                for (int iter = 0; iter < Iters; iter++)
                {
                    int EdgeCount = PDeepCopy.Halfedges.Count / 2;
                    double[] EdgeLength = PDeepCopy.Halfedges.GetLengths();
                    List<bool> Visited = new List<bool>();
                    Vector3d[] Normals = new Vector3d[PDeepCopy.Vertices.Count];

                    for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
                    {
                        Visited.Add(false);
                        Normals[i] = Normal(PDeepCopy, i);
                    }

                    // double t = LengthTol;     //a tolerance for when to split/collapse edges
                    double smooth = SmoothStrength;  //smoothing strength
                    double pull = PullStrength;  //pull to target mesh strength




                    EdgeCount = PDeepCopy.Halfedges.Count / 2;

                    if ((Flip == 0) && (PullStrength > 0))
                    {
                        //Flip edges to reduce valence error
                        for (int i = 0; i < EdgeCount; i++)
                        {
                            if (!PDeepCopy.Halfedges[2 * i].IsUnused
                              && (PDeepCopy.Halfedges[2 * i].AdjacentFace != -1)
                              && (PDeepCopy.Halfedges[2 * i + 1].AdjacentFace != -1))
                            {
                                int Vert1 = PDeepCopy.Halfedges[2 * i].StartVertex;
                                int Vert2 = PDeepCopy.Halfedges[2 * i + 1].StartVertex;
                                int Vert3 = PDeepCopy.Halfedges[PDeepCopy.Halfedges[PDeepCopy.Halfedges[2 * i].NextHalfedge].NextHalfedge].StartVertex;
                                int Vert4 = PDeepCopy.Halfedges[PDeepCopy.Halfedges[PDeepCopy.Halfedges[2 * i + 1].NextHalfedge].NextHalfedge].StartVertex;

                                int Valence1 = PDeepCopy.Vertices.GetValence(Vert1);
                                int Valence2 = PDeepCopy.Vertices.GetValence(Vert2);
                                int Valence3 = PDeepCopy.Vertices.GetValence(Vert3);
                                int Valence4 = PDeepCopy.Vertices.GetValence(Vert4);

                                if (PDeepCopy.Vertices.NakedEdgeCount(Vert1) > 0) { Valence1 += 2; }
                                if (PDeepCopy.Vertices.NakedEdgeCount(Vert2) > 0) { Valence2 += 2; }
                                if (PDeepCopy.Vertices.NakedEdgeCount(Vert3) > 0) { Valence3 += 2; }
                                if (PDeepCopy.Vertices.NakedEdgeCount(Vert4) > 0) { Valence4 += 2; }

                                int CurrentError =
                                  Math.Abs(Valence1 - 6) +
                                  Math.Abs(Valence2 - 6) +
                                  Math.Abs(Valence3 - 6) +
                                  Math.Abs(Valence4 - 6);
                                int FlippedError =
                                  Math.Abs(Valence1 - 7) +
                                  Math.Abs(Valence2 - 7) +
                                  Math.Abs(Valence3 - 5) +
                                  Math.Abs(Valence4 - 5);
                                if (CurrentError > FlippedError)
                                {
                                    PDeepCopy.Halfedges.FlipEdge(2 * i);
                                }
                            }
                        }
                    }
                    else
                    {
                        //Flip edges based on angle
                        for (int i = 0; i < EdgeCount; i++)
                        {
                            if (!PDeepCopy.Halfedges[2 * i].IsUnused
                              && (PDeepCopy.Halfedges[2 * i].AdjacentFace != -1)
                              && (PDeepCopy.Halfedges[2 * i + 1].AdjacentFace != -1))
                            {
                                int Vert1 = PDeepCopy.Halfedges[2 * i].StartVertex;
                                int Vert2 = PDeepCopy.Halfedges[2 * i + 1].StartVertex;
                                int Vert3 = PDeepCopy.Halfedges[PDeepCopy.Halfedges[PDeepCopy.Halfedges[2 * i].NextHalfedge].NextHalfedge].StartVertex;
                                int Vert4 = PDeepCopy.Halfedges[PDeepCopy.Halfedges[PDeepCopy.Halfedges[2 * i + 1].NextHalfedge].NextHalfedge].StartVertex;

                                Point3d P1 = PDeepCopy.Vertices[Vert1].ToPoint3d();
                                Point3d P2 = PDeepCopy.Vertices[Vert2].ToPoint3d();
                                Point3d P3 = PDeepCopy.Vertices[Vert3].ToPoint3d();
                                Point3d P4 = PDeepCopy.Vertices[Vert4].ToPoint3d();

                                double A1 = Vector3d.VectorAngle(new Vector3d(P3 - P1), new Vector3d(P4 - P1))
                                  + Vector3d.VectorAngle(new Vector3d(P4 - P2), new Vector3d(P3 - P2));

                                double A2 = Vector3d.VectorAngle(new Vector3d(P1 - P4), new Vector3d(P2 - P4))
                                  + Vector3d.VectorAngle(new Vector3d(P2 - P3), new Vector3d(P1 - P3));

                                if (A2 > A1)
                                {
                                    PDeepCopy.Halfedges.FlipEdge(2 * i);
                                }
                            }
                        }
                    }

                    if (Minim)
                    {
                        Vector3d[] SmoothC = LaplacianSmooth(PDeepCopy, 1, smooth);

                        for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
                        {
                            if (AnchorV[i] == -1) // don't smooth feature vertices
                            {
                                PDeepCopy.Vertices.MoveVertex(i, 0.5 * SmoothC[i]);
                            }
                        }
                    }

                    Vector3d[] Smooth = LaplacianSmooth(PDeepCopy, 0, smooth);

                    for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
                    {
                        if (AnchorV[i] == -1) // don't smooth feature vertices
                        {
                            // make it tangential only
                            Vector3d VNormal = Normal(PDeepCopy, i);
                            double ProjLength = Smooth[i] * VNormal;
                            Smooth[i] = Smooth[i] - (VNormal * ProjLength);

                            PDeepCopy.Vertices.MoveVertex(i, Smooth[i]);

                            if (PDeepCopy.Vertices.NakedEdgeCount(i) != 0)//special smoothing for feature edges
                            {
                                int[] Neighbours = PDeepCopy.Vertices.GetVertexNeighbours(i);
                                int ncount = 0;
                                Point3d Avg = new Point3d();

                                for (int j = 0; j < Neighbours.Length; j++)
                                {
                                    if (PDeepCopy.Vertices.NakedEdgeCount(Neighbours[j]) != 0)
                                    {
                                        ncount++;
                                        Avg = Avg + PDeepCopy.Vertices[Neighbours[j]].ToPoint3d();
                                    }
                                }
                                Avg = Avg * (1.0 / ncount);
                                Vector3d move = Avg - PDeepCopy.Vertices[i].ToPoint3d();
                                move = move * smooth;
                                PDeepCopy.Vertices.MoveVertex(i, move);
                            }


                            //projecting points onto the target along their normals

                            if (pull > 0)
                            {
                                Point3d Point = PDeepCopy.Vertices[i].ToPoint3d();
                                Vector3d normal = Normal(PDeepCopy, i);
                                Ray3d Ray1 = new Ray3d(Point, normal);
                                Ray3d Ray2 = new Ray3d(Point, -normal);
                                double RayPt1 = Rhino.Geometry.Intersect.Intersection.MeshRay(M, Ray1);
                                double RayPt2 = Rhino.Geometry.Intersect.Intersection.MeshRay(M, Ray2);
                                Point3d ProjectedPt;

                                if ((RayPt1 < RayPt2) && (RayPt1 > 0) && (RayPt1 < 1.0))
                                {
                                    ProjectedPt = Point * (1 - pull) + pull * Ray1.PointAt(RayPt1);
                                }
                                else if ((RayPt2 < RayPt1) && (RayPt2 > 0) && (RayPt2 < 1.0))
                                {
                                    ProjectedPt = Point * (1 - pull) + pull * Ray2.PointAt(RayPt2);
                                }
                                else
                                {
                                    ProjectedPt = Point * (1 - pull) + pull * M.ClosestPoint(Point);
                                }

                                PDeepCopy.Vertices.SetVertex(i, ProjectedPt);
                            }

                        }
                        else
                        {
                            PDeepCopy.Vertices.SetVertex(i, FV[AnchorV[i]]); //pull anchor vertices onto their points
                        }
                    }


                    //end new

                    AnchorV = CompactByVertex(PDeepCopy, AnchorV); //compact the fixed points along with the vertices
                    // FeatureV = CompactByVertex(P, FeatureV);
                    // FeatureE = CompactByEdge(P, FeatureE);

                    PDeepCopy.Compact(); //this cleans the mesh data structure of unused elements
                }
            }
            DA.SetData(0, PDeepCopy);



            #region 更新nodes
            for (int i = 0; i < nodesDeepCopy.Count; i++)
            {
                nodesDeepCopy[i].NodeVertex = PDeepCopy.Vertices.GetPositions()[i].ToPoint3d();
            }
            #endregion


            #region 设置用于可视化的TextDot
            InnerNodeTextDot.Clear();
            OuterNodeTextDot.Clear();

            // List<List<int>> selectedHalfedgeMeshPGraphLoL = allPossiblePGraphLoL[index];
            // List<Node> selectedHalfedgeMeshPNodes = allPossiblePNodes[index];


            for (int i = 0; i < nodesDeepCopy.Count; i++)
            {
                if (nodesDeepCopy[i].IsInner)
                {
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, nodesDeepCopy[i].NodeAttribute.NodeLabel), nodesDeepCopy[i].NodeVertex);
                    InnerNodeTextDot.Add(textDot);
                }
                else
                {
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, nodesDeepCopy[i].NodeAttribute.NodeLabel), nodesDeepCopy[i].NodeVertex);
                    OuterNodeTextDot.Add(textDot);
                }
            }
            #endregion

            #region 转换为RhinoMesh，在DrawViewportMeshes绘制选中的mesh
            PRhinoMesh = RhinoSupport.ToRhinoMesh(PDeepCopy);
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
            GraphEdges.AddRange(GraphEdgeLine(pGraphLoL, nodesDeepCopy));
            #endregion

        }

        private Vector3d Normal(PlanktonMesh P, int V)
        {
            Point3d Vertex = P.Vertices[V].ToPoint3d();
            Vector3d Norm = new Vector3d();

            int[] OutEdges = P.Vertices.GetHalfedges(V);
            int[] Neighbours = P.Vertices.GetVertexNeighbours(V);
            Vector3d[] OutVectors = new Vector3d[Neighbours.Length];
            int Valence = P.Vertices.GetValence(V);

            for (int j = 0; j < Valence; j++)
            {
                OutVectors[j] = P.Vertices[Neighbours[j]].ToPoint3d() - Vertex;
            }

            for (int j = 0; j < Valence; j++)
            {
                if (P.Halfedges[OutEdges[(j + 1) % Valence]].AdjacentFace != -1)
                {
                    Norm += (Vector3d.CrossProduct(OutVectors[(j + 1) % Valence], OutVectors[j]));
                }
            }

            Norm.Unitize();
            return Norm;
        }

        private Point3d MidPt(PlanktonMesh P, int E)
        {
            Point3d Pos1 = P.Vertices[P.Halfedges[2 * E].StartVertex].ToPoint3d();
            Point3d Pos2 = P.Vertices[P.Halfedges[2 * E + 1].StartVertex].ToPoint3d();
            return (Pos1 + Pos2) * 0.5;
        }

        private List<int> CompactByVertex(PlanktonMesh P, List<int> L)
        {
            List<int> L2 = new List<int>();

            for (int i = 0; i < P.Vertices.Count; i++)
            {
                if (P.Vertices[i].IsUnused == false)
                {
                    L2.Add(L[i]);
                }
            }
            return L2;
        }

        private static Vector3d[] LaplacianSmooth(PlanktonMesh P, int W, double Strength)
        {
            int VertCount = P.Vertices.Count;
            Vector3d[] Smooth = new Vector3d[VertCount];

            for (int i = 0; i < VertCount; i++)
            {
                if ((P.Vertices[i].IsUnused == false) && (P.Vertices.IsBoundary(i) == false))
                {
                    int[] Neighbours = P.Vertices.GetVertexNeighbours(i);
                    Point3d Vertex = P.Vertices[i].ToPoint3d();
                    Point3d Centroid = new Point3d();
                    if (W == 0)
                    {
                        for (int j = 0; j < Neighbours.Length; j++)
                        { Centroid = Centroid + P.Vertices[Neighbours[j]].ToPoint3d(); }
                        Smooth[i] = ((Centroid * (1.0 / P.Vertices.GetValence(i))) - Vertex) * Strength;
                    }
                    if (W == 1)
                    {
                        //get the radial vectors of the 1-ring
                        //get the vectors around the 1-ring
                        //get the cotangent weights for each edge

                        int valence = Neighbours.Length;

                        Point3d[] NeighbourPts = new Point3d[valence];
                        Vector3d[] Radial = new Vector3d[valence];
                        Vector3d[] Around = new Vector3d[valence];
                        double[] CotWeight = new double[valence];
                        double WeightSum = 0;

                        for (int j = 0; j < valence; j++)
                        {
                            NeighbourPts[j] = P.Vertices[Neighbours[j]].ToPoint3d();
                            Radial[j] = NeighbourPts[j] - Vertex;
                        }

                        for (int j = 0; j < valence; j++)
                        {
                            Around[j] = NeighbourPts[(j + 1) % valence] - NeighbourPts[j];
                        }

                        for (int j = 0; j < Neighbours.Length; j++)
                        {
                            //get the cotangent weights
                            int previous = (j + valence - 1) % valence;
                            Vector3d Cross1 = Vector3d.CrossProduct(Radial[previous], Around[previous]);
                            double Cross1Length = Cross1.Length;
                            double Dot1 = Radial[previous] * Around[previous];

                            int next = (j + 1) % valence;
                            Vector3d Cross2 = Vector3d.CrossProduct(Radial[next], Around[j]);
                            double Cross2Length = Cross2.Length;
                            double Dot2 = Radial[next] * Around[j];

                            CotWeight[j] = Math.Abs(Dot1 / Cross1Length) + Math.Abs(Dot2 / Cross2Length);
                            WeightSum += CotWeight[j];
                        }

                        double InvWeightSum = 1.0 / WeightSum;

                        Vector3d ThisSmooth = new Vector3d();

                        for (int j = 0; j < Neighbours.Length; j++)
                        {
                            ThisSmooth = ThisSmooth + Radial[j] * CotWeight[j];
                        }

                        Smooth[i] = ThisSmooth * InvWeightSum * Strength;
                    }

                }
            }
            return Smooth;
        }


        /// <summary>
        /// 绘制表示图结构关系的Line
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="graphVertices"></param>
        /// <returns></returns>
        public List<Line> GraphEdgeLine(List<List<int>> graph, List<GraphNode> graphNodes)
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
            get { return new Guid("0f550e45-ec2a-4958-bf39-f051b1f901c7"); }
        }
    }
}