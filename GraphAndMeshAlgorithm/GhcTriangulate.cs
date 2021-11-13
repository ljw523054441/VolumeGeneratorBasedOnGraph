using Grasshopper.Kernel;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcTriangulate : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Triangulate class.
        /// </summary>
        public GhcTriangulate()
          : base("Triangulate", "Triangulate",
              "Finds all possible triangulations of a [convex] mesh",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            SelectedIsomorphismTriangleMesh = new Mesh();

            ConvexPolylinesPoints = new List<List<Point3d>>();
            SelectedTriangleMeshEdges = new List<Line>();

            InnerNodeTextDot = new List<TextDot>();
            OuterNodeTextDot = new List<TextDot>();

            DottedCurve = new List<Curve>();
        }

        private int Thickness;

        private Mesh SelectedIsomorphismTriangleMesh;
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
            // 0
            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);
            // 1
            pManager.AddCurveParameter("ConvexFaceBorders", "CFBorders", "Convex face borders of the Tutte algorithm as a list of polyline curves.", GH_ParamAccess.list);
            // 2
            pManager.AddIntegerParameter("IndexOfTriangularMesh", "I", "Index of a triangulation to be visualized from the list of all triangulations", GH_ParamAccess.item);
            // 3
            pManager.AddBooleanParameter("ExcludeDegenerateTINS", "ExT", "排除那些产生三角形房间的三角剖分 Exclude those triangulations that give rise to triangular rooms", GH_ParamAccess.item);
            pManager[2].Optional = true;
            pManager[3].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // 0
            pManager.AddMeshParameter("AllTriangularMeshes", "AllTMesh", "所有可能的三角形剖分结果。All Computed triangulations; these triangulations describe all possible planar topologies for your plan diagram, the added links are those of adjacencies not connectivities", GH_ParamAccess.list);
            // 1
            // pManager.AddMeshParameter("TheChosenTriangularMesh", "TMesh", "所选择的那个三角形剖分结果。The one you have chosen with index I", GH_ParamAccess.item);
            // 2
            pManager.AddGenericParameter("TheChosenGraphWithHM", "THFMesh", "所选择的那个三角形剖分结果(GraphWithHFMesh对象)", GH_ParamAccess.item);
            // 3
            // pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 局部变量初始化
            Graph graph = new Graph();
            List<Curve> convexFaceBorderCurves = new List<Curve>();

            Thickness = 2;
            #endregion

            if (DA.GetData<Graph>("Graph", ref graph)
                && DA.GetDataList<Curve>("ConvexFaceBorders", convexFaceBorderCurves))
            {
                int innerNodeCount = graph.InnerNodeCount;
                // int outerNodeCount = graph.OuterNodeCount;
                // List<int> innerNodeIndexList = graph.InnerNodeIndexList;
                // List<int> outerNodeIndexList = graph.OuterNodeIndexList;

                #region 对输入的Curve类型的ConvexFaceBorder进行类型转换，转换成Curve类的子类Polyline
                List<Polyline> convexFaceBorderPolylines = new List<Polyline>();
                for (int i = 0; i < convexFaceBorderCurves.Count; i++)
                {
                    Polyline polyline = null;
                    if (convexFaceBorderCurves[i].TryGetPolyline(out polyline))
                    {
                        if (polyline.IsClosed)
                        {
                            convexFaceBorderPolylines.Add(polyline);

                            #region 设置可视化中的实线部分
                            ConvexPolylinesPoints.Add(new List<Point3d>());
                            ConvexPolylinesPoints[i].AddRange(polyline);
                            #endregion
                        }
                    }
                }
                #region 可能的报错
                if (convexFaceBorderCurves.Count != convexFaceBorderPolylines.Count)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the ConvexFaceBorders are supposed to be closed polylines or polygons, but apparently at least one is not closed or is not a polyline!");
                    return;
                }
                #endregion
                #endregion



                #region 获取Node的坐标和构造Point3d

                List<Node> nodes = graph.GraphNodes;

                List<Point3d> nodePoints = new List<Point3d>();
                List<Point3d> outerPoints = new List<Point3d>();

                for (int i = 0; i < nodes.Count; i++)
                {
                    nodePoints.Add(nodes[i].NodeVertex);
                }
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (!nodes[i].IsInner)
                    {
                        outerPoints.Add(nodes[i].NodeVertex);
                    }
                }
                #endregion

                int index = 0;
                bool flag = false;

                DA.GetData<int>("IndexOfTriangularMesh", ref index);
                DA.GetData<bool>("ExcludeDegenerateTINS", ref flag);

                #region 得到这样一个图结构下，所有的同形异构的整体的三角形网格
                List<Mesh> allPolylineCorrespondIsomorphismTriangleMeshes = GetAllIsomorphismTriangleMeshes(ClosedPolylineToTriangleMesh(convexFaceBorderPolylines));
                #endregion

                #region 剔除细分后的三角形全部由outerNode构成的情况
                List<Mesh> allPolylineCorrespondIsomorphismTriangleMeshesExceptOuterNode = RemoveTriangleMeshWithAllOuterNode(allPolylineCorrespondIsomorphismTriangleMeshes, outerPoints);
                #endregion

                #region 对于新生成的三角网格，其中的顶点顺序是跟原来的nodePoints列表中的点的顺序不同的，需要重新匹配
                List<Mesh> regeneratedIsomorphismMeshes = new List<Mesh>();
                foreach (Mesh IsomorphismMesh in allPolylineCorrespondIsomorphismTriangleMeshesExceptOuterNode)
                {
                    regeneratedIsomorphismMeshes.Add(MatchVerticesIndex(nodePoints, IsomorphismMesh));
                }
                // 是否删除存在顶点的度为3的剖分
                List<Mesh> result;
                if (flag)
                {
                    result = RemoveTriangleMeshWithInnerNodeDegree3(regeneratedIsomorphismMeshes, innerNodeCount);
                }
                else
                {
                    result = regeneratedIsomorphismMeshes;
                }
                #endregion
                DA.SetDataList("AllTriangularMeshes", result);

                #region 选择生成的三角剖分
                if (index >= result.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的序号超过了所有三角形网格的总数");
                    return;
                }
                #endregion
                // DA.SetData("TheChosenTriangularMesh", result[index]);

                #region 构造半边数据结构
                PlanktonMesh planktonMesh = new PlanktonMesh();

                for (int i = 0; i < nodePoints.Count; i++)
                {
                    planktonMesh.Vertices.Add(nodePoints[i].X, nodePoints[i].Y, nodePoints[i].Z);
                }

                List<List<int>> faceVertexOrder = new List<List<int>>();
                for (int i = 0; i < result[index].Faces.Count; i++)
                {
                    int[] faceTopologicalVerticesArray = result[index].Faces.GetTopologicalVertices(i);
                    List<int> faceTopologicalVerticesList = faceTopologicalVerticesArray.ToList<int>();
                    faceTopologicalVerticesList.RemoveAt(faceTopologicalVerticesList.Count - 1);

                    faceVertexOrder.Add(new List<int>());
                    faceVertexOrder[i].AddRange(faceTopologicalVerticesList);
                }
                planktonMesh.Faces.AddFaces(faceVertexOrder);
                #endregion

                GraphWithHM theChosenOne = new GraphWithHM(planktonMesh, graph.GraphNodes, graph.GraphTables);

                DA.SetData("TheChosenGraphWithHM", theChosenOne);

                //#region Debug显示

                //#region HalfedgeMesh的顶点数据
                //List<string> printVertices = new List<string>();
                //printVertices = UtilityFunctions.PrintVertices(planktonMesh);
                //DA.SetDataList("DebugVerticesOutput", printVertices);
                //#endregion

                //#region HalfedgeMesh的半边数据
                //List<string> printHalfedges = new List<string>();
                //printHalfedges = UtilityFunctions.PrintHalfedges(planktonMesh);
                //DA.SetDataList("DebugHalfedgesOutput", printHalfedges);
                //#endregion

                //#region HalfedgeMesh的每个面由哪些顶点构成
                //List<string> printFaces = new List<string>();
                //printFaces = UtilityFunctions.PrintFacesVertices(planktonMesh);
                //DA.SetDataList("DebugFacesOutput", printFaces);
                //#endregion

                //#region HalfedgeMesh的每个面由哪些半边构成
                //List<string> printFacesHalfedge = new List<string>();
                //printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(planktonMesh);
                //DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);
                //#endregion

                //#endregion

                #region 可视化部分
                // ConvexPolylinesPoints.Clear();
                SelectedTriangleMeshEdges.Clear();
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();

                #region 设置用于可视化的Mesh
                // 在DrawViewportMeshes绘制选中的mesh
                SelectedIsomorphismTriangleMesh = result[index];
                #endregion

                #region 设置用于可视化的虚线
                // 在DrawViewportWires绘制mesh的edge
                for (int i = 0; i < SelectedIsomorphismTriangleMesh.TopologyEdges.Count; i++)
                {
                    SelectedTriangleMeshEdges.Add(SelectedIsomorphismTriangleMesh.TopologyEdges.EdgeLine(i));
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

                #region 设置用于可视化的TextDot
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
                #endregion
                #endregion
            }
        }

        /// <summary>
        /// 将每个polyline所对应的一列表代表所有可能的连接方式的小的Convex形状的三角网格，穷尽所有的排列组合，合成一列表大的三角网格
        /// </summary>
        /// <param name="allPolylineCorrespondIsomorphismTriangleMeshes"></param>
        /// <returns></returns>
        public List<Mesh> GetAllIsomorphismTriangleMeshes(List<List<Mesh>> allPolylineCorrespondIsomorphismTriangleMeshes)
        {
            List<Mesh> result;

            if (allPolylineCorrespondIsomorphismTriangleMeshes.Count <= 1)
            {
                result = null;
            }
            else
            {
                while (allPolylineCorrespondIsomorphismTriangleMeshes.Count > 1)
                {
                    List<Mesh> combinedMeshTypes = new List<Mesh>();
                    List<Mesh> polyline_0_CorrespondIsomorphismTriangleMeshes = allPolylineCorrespondIsomorphismTriangleMeshes[0];
                    List<Mesh> polyline_1_CorrespondIsomorphismTriangleMeshes = allPolylineCorrespondIsomorphismTriangleMeshes[1];

                    // 穷尽所有组合，将每个convex的同形异构的三角形网格，组合形成一堆同形异构的整体的三角形网格
                    foreach (Mesh polyline_0_IsomorphismTriangleMeshes in polyline_0_CorrespondIsomorphismTriangleMeshes)
                    {
                        foreach (Mesh polyline_1_IsomorphismTriangleMeshes in polyline_1_CorrespondIsomorphismTriangleMeshes)
                        {
                            Mesh mesh = new Mesh();
                            mesh.Append(polyline_0_IsomorphismTriangleMeshes);
                            mesh.Append(polyline_1_IsomorphismTriangleMeshes);
                            mesh.Vertices.CombineIdentical(true, true);
                            combinedMeshTypes.Add(mesh);
                        }
                    }
                    allPolylineCorrespondIsomorphismTriangleMeshes.RemoveAt(0);
                    allPolylineCorrespondIsomorphismTriangleMeshes.RemoveAt(0);
                    allPolylineCorrespondIsomorphismTriangleMeshes.Insert(0, combinedMeshTypes);
                }
                
                result = allPolylineCorrespondIsomorphismTriangleMeshes[0];
            }
            return result;
        }

        /// <summary>
        /// 由封闭的polyline生成树形结构的三角形网格
        /// </summary>
        /// <param name="closedPolylines"></param>
        /// <returns></returns>
        public List<List<Mesh>> ClosedPolylineToTriangleMesh(List<Polyline> closedPolylines)
        {
            if (closedPolylines.Count <= 1)
            {
                return null;
            }
            else
            {
                List<List<Mesh>> allPolylineCorrespondTriangleMesh = new List<List<Mesh>>();
                foreach (Polyline polyline in closedPolylines)
                {
                    // Polyline类型转化为ConvexPolygon
                    ConvexPolygon eachConvexPolygon = new ConvexPolygon(polyline);
                    // 将每个polyline对应的ConvexPolygon进行三角剖分后的的LoLMeshFace，添加上vertice信息，形成真正的mesh。
                    // 列表中每一个元素，是由这个polyline形成的对应的一种三角形网格（即三角剖分方式）
                    List<Mesh> polylineCorrespondTrianglationType = GenerateMeshOfTriangleConvexPolygon(eachConvexPolygon.TriangulateConvexPolygon(), polyline.GetRange(0, polyline.Count - 1));
                    List<Mesh> cleanedPolylineCorrespondTrianglationType = RemoveInvaildTriangleMesh(polylineCorrespondTrianglationType);
                    allPolylineCorrespondTriangleMesh.Add(cleanedPolylineCorrespondTrianglationType);
                }
                return allPolylineCorrespondTriangleMesh;
            }
        }

        /// <summary>
        /// 根据输入的三角形MeshFace的树形结构，每支上的每个MeshFace只是顶点的排序，需要配合实际的Point3d组合形成一个Mesh，最后输出Mesh列表
        /// </summary>
        /// <param name="LoLConvexPolygonMeshFace">mesh网格中每个三角面之间点的序号顺序</param>
        /// <param name="meshVertices">组成mesh网格的三角面的顶点</param>
        /// <returns></returns>
        public List<Mesh> GenerateMeshOfTriangleConvexPolygon(List<List<MeshFace>> LoLConvexPolygonMeshFace, IEnumerable<Point3d> meshVertices)
        {
            // 生成的三角形mesh的列表
            List<Mesh> generatedTriangleMesh = new List<Mesh>();

            // 对于整个大的ConvexPolygon进行三角化后，形成的树形结构的每一支，每支上有好多子MeshFace，把它们按支组合成一个Mesh
            // LoLConvexPolygonMeshFace中的每支上表示不同的对ConvexPolygon的三角剖分方式，每一支上的每个meshFace表示这种剖分方式下的一个三角形
            foreach (List<MeshFace> subTriangleConvexMeshFaceList in LoLConvexPolygonMeshFace)
            {
                Mesh mesh = new Mesh();
                mesh.Vertices.AddVertices(meshVertices);
                mesh.Faces.AddFaces(subTriangleConvexMeshFaceList);
                generatedTriangleMesh.Add(mesh);
            }

            return generatedTriangleMesh;
        }

        /// <summary>
        /// 删除无效的三角形网格(face中心点到边的距离小于0.001的)
        /// </summary>
        /// <param name="triangleMeshList"></param>
        /// <returns></returns>
        public List<Mesh> RemoveInvaildTriangleMesh(List<Mesh> triangleMeshList)
        {
            List<Mesh> resultTriangleMeshList = new List<Mesh>();

            foreach (Mesh mesh in triangleMeshList)
            {
                MeshVertexList vertices = mesh.Vertices;
                MeshTopologyEdgeList topologyEdges = mesh.TopologyEdges;
                MeshFaceList faces = mesh.Faces;

                List<int> unqualifiedEdgeIndexList = new List<int>();
                // 对于mesh的每一个面
                for (int i = 0; i < faces.Count; i++)
                {
                    // 找到构成第i个face的边的序号列表
                    int[] edgesForFace = topologyEdges.GetEdgesForFace(i);
                    // 对于每个构成第i个面的边
                    for (int j = 0; j < edgesForFace.Length; j++)
                    {
                        // 如果这个face的中心点到这个边的距离，注意不一定是垂直距离，小于0.001的话，就把这条边的序号加入列表
                        if ((topologyEdges.EdgeLine(edgesForFace[j])).DistanceTo(faces.GetFaceCenter(i),true) < 0.001)
                        {
                            unqualifiedEdgeIndexList.Add(edgesForFace[j]);
                        }
                    }
                }

                // 如果这个mesh的每个面的每条边都合规，那么就允许这个mesh输出，否则这个mesh就不输出（即删掉）
                if (unqualifiedEdgeIndexList.Count == 0)
                {
                    resultTriangleMeshList.Add(mesh);
                }
            }
            return resultTriangleMeshList;
        }

        /// <summary>
        /// 排除细分后的三角形全部由outerNode构成的情况
        /// </summary>
        /// <param name="triangleMeshList"></param>
        /// <param name="outerPoints"></param>
        /// <returns></returns>
        public List<Mesh> RemoveTriangleMeshWithAllOuterNode(List<Mesh> triangleMeshList, List<Point3d> outerPoints)
        {
            List<Mesh> resultTriangleMeshList = new List<Mesh>();

            foreach (Mesh mesh in triangleMeshList)
            {
                MeshFaceList faces = mesh.Faces;
                int faceCount = 0;

                for (int i = 0; i < faces.Count; i++)
                {
                    int vertexContainsCount = 0;

                    List<Point3d> faceVertices = new List<Point3d>();
                    Point3f a;
                    Point3f b;
                    Point3f c;
                    Point3f d;
                    faces.GetFaceVertices(i, out a, out b, out c, out d);
                    faceVertices.Add(a);
                    faceVertices.Add(b);
                    faceVertices.Add(c);
                    
                    foreach (Point3d vertex in faceVertices)
                    {
                        Point3dList point3dList = new Point3dList(outerPoints);
                        if (EpsilonContains(point3dList, vertex, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance))
                        {
                            vertexContainsCount++;
                        }
                    }
                    if (vertexContainsCount == 3)
                    {
                        continue;
                    }
                    else
                    {
                        faceCount++;
                    }

                }
                if (faceCount == faces.Count)
                {
                    resultTriangleMeshList.Add(mesh);
                }
                else
                {
                    continue;
                }
            }
            return resultTriangleMeshList;
        }

        /// <summary>
        /// 让新生成的同构异形三角网格顶点序号，与原来的nodePoints的序号做匹配。点的坐标是一样的，只是序号要重新匹配
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="isomorphismMesh"></param>
        /// <returns></returns>
        public Mesh MatchVerticesIndex(List<Point3d> vertices, Mesh isomorphismMesh)
        {
            Point3d[] isomorphismMeshVerticesArray = isomorphismMesh.Vertices.ToPoint3dArray();
            Point3dList nodeVertices = new Point3dList(vertices);
            List<int> indexList = new List<int>();
            foreach (Point3d testPoint in isomorphismMeshVerticesArray)
            {
                indexList.Add(nodeVertices.ClosestIndex(testPoint));
            }

            List<MeshFace> regeneratedMeshFaces = new List<MeshFace>();
            foreach (MeshFace meshFace in isomorphismMesh.Faces)
            {
                MeshFace regeneratedMeshFace = new MeshFace(indexList[meshFace.A], indexList[meshFace.B], indexList[meshFace.C]);
                regeneratedMeshFaces.Add(regeneratedMeshFace);
            }

            Mesh regeneratedMesh = new Mesh();
            regeneratedMesh.Vertices.AddVertices(vertices);
            regeneratedMesh.Faces.AddFaces(regeneratedMeshFaces);
            return regeneratedMesh;
        }

        /// <summary>
        /// 通过顶点的度与生成的对偶图多边形边数的对应关系，来删除对偶图中出现三角形的情况（即该顶点的度为3的情况）
        /// </summary>
        /// <param name="triangleMeshes"></param>
        /// <returns></returns>
        public List<Mesh> RemoveTriangleMeshWithInnerNodeDegree3(List<Mesh> triangleMeshes, int innerNodeCount)
        {
            List<Mesh> removedTriangleMesh = new List<Mesh>();

            foreach (Mesh triangleMesh in triangleMeshes)
            {
                List<int> verticeDegrees = new List<int>();
                for (int i = 0; i < innerNodeCount; i++)
                {
                    int verticeDegree = triangleMesh.TopologyVertices.ConnectedTopologyVertices(i).Length;

                    // list2.Add(Enumerable.Count<int>(triangleMesh.TopologyVertices.ConnectedTopologyVertices(i)));
                    verticeDegrees.Add(verticeDegree);

                }


                  
                if (!verticeDegrees.Contains(3))
                {
                    removedTriangleMesh.Add(triangleMesh);
                }
            }
            return removedTriangleMesh;
        }

        /// <summary>
        /// 判断vertices中的点跟对应的最近的testpoint中的点，距离是否小于公差
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="testPoint"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public bool EpsilonContains(Point3dList vertices, Point3d testPoint, double tolerance)
        {
            Point3d point3d = vertices[vertices.ClosestIndex(testPoint)];
            bool bool1 = Math.Abs(point3d.X - testPoint.X) < tolerance;
            bool bool2 = Math.Abs(point3d.Y - testPoint.Y) < tolerance;
            bool bool3 = Math.Abs(point3d.Z - testPoint.Z) < tolerance;
            return bool1 & bool2 & bool3;
        }

        /// <summary>
        /// 预览模式为WireFrame模式时，调用此函数
        /// </summary>
        /// <param name="args"></param>
        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);

            // args.Display.DrawMeshShaded(SelectedIsomorphismTriangleMesh, new Rhino.Display.DisplayMaterial(Color.White, 0));

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

            args.Display.DrawMeshShaded(SelectedIsomorphismTriangleMesh, new Rhino.Display.DisplayMaterial(Color.White, 0));

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
            get { return new Guid("88a4820f-ef86-4c32-ba34-74e80870373f"); }
        }
    }
}