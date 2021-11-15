using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using System;
using System.Collections.Generic;
using System.Drawing;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class DualGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DualGraph class.
        /// </summary>
        public DualGraph()
          : base("DualGraph", "DualGraph",
              "生成对应的对偶图",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            DualConvexVertices = new List<Point3d>();
            DualConvexPolygon = new List<List<Point3d>>();
            DualConvexVertexTextDots = new List<TextDot>();

            DualConvexCenterPoints = new List<Point3d>();
            GraphEdges = new List<Line>();
            NodeTextDots = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> DualConvexVertices;
        private List<List<Point3d>> DualConvexPolygon;
        private List<TextDot> DualConvexVertexTextDots;

        private List<Point3d> DualConvexCenterPoints;
        private List<Line> GraphEdges;
        private List<TextDot> NodeTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddMeshParameter("TriangleMesh", "TMesh", "所选择的那个三角形剖分结果。", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);

            // pManager.AddGenericParameter("SubAttributes", "SubAttributes", "(Optional)Put labels, areas and colors as an attribute LOL, if you want to get a legned for your analysis.", GH_ParamAccess.list);
            // pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("DualVertices", "DV", "所有构成对偶多边形的顶点", GH_ParamAccess.list);
            pManager.AddTextParameter("DualFacesDescriptions", "DFD", "生成的对偶面的描述，包括对于第几个innerNode，周围由哪几个三角形的中心点包围形成", GH_ParamAccess.list);
            pManager.AddCurveParameter("DualConvexPolygons", "DFB", "生成的对偶多边形", GH_ParamAccess.list);

            pManager.AddIntegerParameter("DualConvexConnectivityGraph", "DCCGraph", "A graph representation as an adjacency list datatree", GH_ParamAccess.list);
            pManager.AddPointParameter("DualConvexCenterPointsAsGraphNodes", "DCGraphNode", "Graph vertices as [edge] centrpoids of cells", GH_ParamAccess.list);
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

            Mesh mesh = new Mesh();
            // List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();

            List<GraphNode> nodes = new List<GraphNode>();

            if (DA.GetData<Mesh>("TriangleMesh", ref mesh))
            {
                // 获取节点
                DA.GetDataList<GraphNode>("GraphNode", nodes);






                // 获取网格上每个面的闭合polyline，组成列表List<Polyline> facePolylines，注意这里获取的一切都不能保持原来node的顺序
                MeshFaceList meshFaces = mesh.Faces;
                MeshVertexList meshVertices = mesh.Vertices;
                List<Polyline> facePolylines = new List<Polyline>();

                for (int i = 0; i < meshFaces.Count; i++)
                {
                    Point3f pointA;
                    Point3f pointB;
                    Point3f pointC;
                    Point3f pointD;

                    meshFaces.GetFaceVertices(i, out pointA, out pointB, out pointC, out pointD);
                    Point3dList faceVerticesPointList = new Point3dList();
                    faceVerticesPointList.AddRange(new Point3d[]
                    {
                        pointA,
                        pointB,
                        pointC,
                        pointA
                    });
                    Polyline facePolyline = new Polyline(faceVerticesPointList);
                    facePolylines.Add(facePolyline);
                }



                // meshVerticesPointsList mesh中所有顶点的集合
                List<Point3d> innerNodePoints = new List<Point3d>();

                // meshVerticesPointList中，是分innerNode和outerNode的，innerNode部分的顺序跟nodePoint原来的顺序一样
                /*
                 * 这里要进行确认
                 */

                // 所以在计算是要把outerNode去掉进行计算
                for (int i = 0; i < meshVertices.Count - boundaryNodeCount; i++)
                {
                    innerNodePoints.Add(meshVertices[i]);
                }

                List<Point3d> triangleCenterPoints = new List<Point3d>();
                List<List<int>> sortedCenterPointIndexLoL = new List<List<int>>();
                List<Polyline> dualConvexPolygons = new List<Polyline>();
                List<List<int>> polylineCenterBelongsToWhichNodeIndexLoL = new List<List<int>>();
                List<string> faceDescriptions = new List<string>();

                ConstructDualConvex(
                    innerNodePoints,
                    facePolylines,
                    ref triangleCenterPoints,
                    ref sortedCenterPointIndexLoL,
                    ref dualConvexPolygons,
                    ref polylineCenterBelongsToWhichNodeIndexLoL,
                    ref faceDescriptions);

                // 输出所有的对偶多边形
                DA.SetDataList("DualConvexPolygons", dualConvexPolygons);
                // 输出所有的对偶多边形的顶点，即每个三角形的中心点
                DA.SetDataList("DualVertices", triangleCenterPoints);
                // 输出生成的对偶面的描述，包括对于第几个innerNode，周围由哪几个三角形的中心点包围形成
                DA.SetDataList("DualFacesDescriptions", faceDescriptions);



                // 将DualConvex由polygon转为curve，方便调用函数
                List<Curve> dualConvexPolygonCurve = new List<Curve>();
                for (int i = 0; i < dualConvexPolygons.Count; i++)
                {
                    dualConvexPolygonCurve.Add(dualConvexPolygons[i].ToNurbsCurve());
                }
                //// 对没有按照node顺序的dualConvexPolygon进行排序
                //List<Curve> sortedDualConvexPolygonCurve = new List<Curve>();
                //for (int i = 0; i < nodes.Count; i++)
                //{
                //    for (int j = 0; j < dualConvexPolygonCurve.Count; j++)
                //    {
                //        // 如果这个对偶多边形包含了某个
                //        PointContainment pointContainment = dualConvexPolygonCurve[j].Contains(nodes[i].NodeVertex, Plane.WorldXY, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                //        if (pointContainment == PointContainment.Inside || pointContainment == PointContainment.Coincident)
                //        {
                //            sortedDualConvexPolygonCurve.Add(dualConvexPolygonCurve[j]);
                //        }
                //    }
                //}

                // 判断对偶的ConvexPolygon是否有没闭合的
                int count = 0;
                foreach (Curve curve in dualConvexPolygonCurve)
                {
                    if (!curve.IsClosed)
                    {
                        count++;
                    }
                }
                if (count > 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All dual Convex must be closed!");
                }

                // 用对偶ConvexPolygon的中心点来代表这个ConvexPolygon
                DataTree<int> DualConvexPolygonConnectivityDT = new DataTree<int>();
                List<Point3d> DualConvexCenterPoints = new List<Point3d>();

                for (int i = 0; i < dualConvexPolygonCurve.Count; i++)
                {
                    DualConvexPolygonConnectivityDT.EnsurePath(i);

                    Polyline polyline = null;
                    if (dualConvexPolygonCurve[i].TryGetPolyline(out polyline))
                    {
                        DualConvexCenterPoints.Add(polyline.CenterPoint());
                    }

                    for (int j = 0; j < dualConvexPolygonCurve.Count; j++)
                    {
                        if (j != i)
                        {
                            // 判断两条Polyline是否相交，如果是，那么它俩有邻接关系
                            RegionContainment regionContainment = Curve.PlanarClosedCurveRelationship(dualConvexPolygonCurve[i], dualConvexPolygonCurve[j], Plane.WorldXY, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                            if (regionContainment == RegionContainment.MutualIntersection)
                            {
                                DualConvexPolygonConnectivityDT.Branch(i).Add(j);
                            }
                        }
                    }
                }
                // 输出对偶多边形之间的连接关系
                DA.SetDataTree(3, DualConvexPolygonConnectivityDT);
                // 输出对偶多边形的中心作为对偶多边形之间的图结构的node
                DA.SetDataList("DualConvexCenterPointsAsGraphNodes", DualConvexCenterPoints);
                





                // 绘图部分
                // 对偶图
                DualConvexVertices.Clear();
                DualConvexPolygon.Clear();
                DualConvexVertexTextDots.Clear();

                DualConvexVertices = triangleCenterPoints;
                for (int i = 0; i < dualConvexPolygons.Count; i++)
                {
                    DualConvexPolygon.Add(new List<Point3d>());
                    DualConvexPolygon[i].AddRange(dualConvexPolygons[i]);
                }
                
                for (int i = 0; i < polylineCenterBelongsToWhichNodeIndexLoL.Count; i++)
                {
                    string arg = string.Join<int>(";", polylineCenterBelongsToWhichNodeIndexLoL[i]);
                    
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), triangleCenterPoints[i]);
                    DualConvexVertexTextDots.Add(textDot);
                }

                // 原来的Graph
                this.DualConvexCenterPoints.Clear();
                this.DualConvexCenterPoints.AddRange(DualConvexCenterPoints);

                GraphEdges.Clear();
                List<Line> edges = new List<Line>();
                for (int i = 0; i < DualConvexPolygonConnectivityDT.BranchCount; i++)
                {
                    for (int j = 0; j < DualConvexPolygonConnectivityDT.Branch(i).Count; j++)
                    {
                        edges.Add(new Line(this.DualConvexCenterPoints[i], this.DualConvexCenterPoints[DualConvexPolygonConnectivityDT.Branch(i)[j]]));
                    }
                }
                GraphEdges.AddRange(edges);

                NodeTextDots.Clear();
                List<TextDot> nodeTextDots = new List<TextDot>();
                for (int i = 0; i < nodes.Count - boundaryNodeCount; i++)
                {
                    nodeTextDots.Add(new TextDot(string.Format("{0} | {1}", i, nodes[i].NodeAttribute.NodeLabel), this.DualConvexCenterPoints[i]));
                }
                NodeTextDots.AddRange(nodeTextDots);
            }
        }

        /// <summary>
        /// 构造对偶多边形，并且对于每个对偶多边形的顶点，即每个三角面的centerPoints，按照每个centerPoint向量与x轴正方向的夹角的弧度值大小，对每个对偶多边形的顶点所邻接的三角面的中心进行排序
        /// </summary>
        /// <param name="innerNodePoints">innerNode的列表</param>
        /// <param name="facePolylines">mesh中每个三角面对应的polyline</param>
        /// <param name="facePolylineCenterPoints">mesh中每个三角面对应的中心点</param>
        /// <param name="sortedCenterPointIndexLoL">innerNode列表中，每个innerNode所对应的经过夹角排序的的centerPoints序号的列表</param>
        /// <param name="dualConvexPolygons">生成的对偶多边形，顺序与innerNodePoints对应</param>
        /// <param name="polylineCenterBelongsToWhichNodeIndexLoL">每个polyline属于哪几个Node，即哪几个Convex</param>
        /// <param name="faceDescriptions">生成的ConvexFace的描述用于debug，包括对于第几个innerNode，周围由哪几个三角形的中心点包围形成</param>
        public void ConstructDualConvex(
            List<Point3d> innerNodePoints,
            List<Polyline> facePolylines,
            ref List<Point3d> facePolylineCenterPoints,
            ref List<List<int>> sortedCenterPointIndexLoL,
            ref List<Polyline> dualConvexPolygons,
            ref List<List<int>> polylineCenterBelongsToWhichNodeIndexLoL,
            ref List<string> faceDescriptions)
        {
            // 获取mesh中每个三角面对应的中心点
            facePolylineCenterPoints = PolygonCenters(facePolylines);

            // 获取innerNodePoints中每个元素，分别属于哪几个三角面对应的polyline
            List<List<int>> innerNodeContainInWhichFacePolylineIndexLoL = PointsInWhichPolylines(innerNodePoints, facePolylines, Plane.WorldXY);
            // 根据上面获取的polyline序号列表，找到innerNodePoints中每个元素，对应的邻接哪些polyline的中心点
            List<List<Point3d>> innerNodeAdjacentWhichPolylineCenterLoL = FetchPolylineCentersForDualVerticesLoL(facePolylineCenterPoints, innerNodeContainInWhichFacePolylineIndexLoL);



            // 对于每个polylineCenter，它是属于那个innerNode
            for (int i = 0; i < facePolylineCenterPoints.Count; i++)
            {
                polylineCenterBelongsToWhichNodeIndexLoL.Add(new List<int>());
                for (int j = 0; j < innerNodeAdjacentWhichPolylineCenterLoL.Count; j++)
                {
                    if (innerNodeAdjacentWhichPolylineCenterLoL[j].Contains(facePolylineCenterPoints[i]))
                    {
                        polylineCenterBelongsToWhichNodeIndexLoL[i].Add(j);
                    }
                }
            }


            if (innerNodeAdjacentWhichPolylineCenterLoL == null | innerNodePoints == null | innerNodePoints.Count != innerNodeAdjacentWhichPolylineCenterLoL.Count)
            {
                return;
            }
            else
            {
                List<List<Point3d>> sortedVertexAdjacentCenterPointLoL = new List<List<Point3d>>();
                List<double> radianList = new List<double>();

                // 对于mesh上每一个顶点
                for (int i = 0; i < innerNodePoints.Count; i++)
                {
                    radianList.Clear();
                    sortedVertexAdjacentCenterPointLoL.Add(new List<Point3d>());
                    sortedCenterPointIndexLoL.Add(new List<int>());
                    

                    // 每个顶点会对应几个PolygonCenter
                    int[] centerPointIndexArray = new int[innerNodeAdjacentWhichPolylineCenterLoL[i].Count];

                    // 对于每一支上的每个元素（即每个PolylineCenter）
                    for (int j = 0; j < innerNodeAdjacentWhichPolylineCenterLoL[i].Count; j++)
                    {
                        // 在centerPoint列表中找到与vertex邻接的对应的centerpoint的序号
                        centerPointIndexArray[j] = facePolylineCenterPoints.IndexOf(innerNodeAdjacentWhichPolylineCenterLoL[i][j]);
                        /* atan2(y, x)是4象限反正切，它的取值不仅取决于正切值y/x，还取决于点 (x, y) 落入哪个象限：
                         * 当点(x, y) 落入第一象限时，atan2(y, x)的范围是 0 ~ pi/2;
                         * 当点(x, y) 落入第二象限时，atan2(y, x)的范围是 pi/2 ~ pi;
                         * 当点(x, y) 落入第三象限时，atan2(y, x)的范围是 －pi～－pi/2;
                         * 当点(x, y) 落入第四象限时，atan2(y, x)的范围是 -pi/2～0. */
                        double radian = Math.Atan2(innerNodeAdjacentWhichPolylineCenterLoL[i][j].Y - innerNodePoints[i].Y, innerNodeAdjacentWhichPolylineCenterLoL[i][j].X - innerNodePoints[i].X);
                        // 得到每个centerPoint向量与x轴正方向的夹角的弧度值大小
                        if (radian >= 0.0)
                        {
                            radianList.Add(radian);
                        }
                        else
                        {
                            radianList.Add(2 * Math.PI + radian);
                        }
                    }

                    // 按照与X轴的夹角大小，对CenterPoint进行排序
                    Point3d[] AdjacentPolylineCenterArray = innerNodeAdjacentWhichPolylineCenterLoL[i].ToArray();
                    Array.Sort<double, Point3d>(radianList.ToArray(), AdjacentPolylineCenterArray);
                    Array.Sort<double, int>(radianList.ToArray(), centerPointIndexArray);

                    sortedVertexAdjacentCenterPointLoL[i].AddRange(AdjacentPolylineCenterArray);
                    sortedCenterPointIndexLoL[i].AddRange(centerPointIndexArray);

                    // 构造对偶的convex，并输出
                    List<Point3d> pointsForDualPolyline = new List<Point3d>();
                    pointsForDualPolyline.AddRange(AdjacentPolylineCenterArray);
                    // 补上起始点，使首尾闭合
                    pointsForDualPolyline.Add(AdjacentPolylineCenterArray[0]);
                    Polyline dualConvexPolyline = new Polyline(pointsForDualPolyline);
                    dualConvexPolygons.Add(dualConvexPolyline);

                    // 用；连接sortedCenterPointIndexLoL中的所有序号
                    string arg = string.Join<int>(";", sortedCenterPointIndexLoL[i]);
                    string description = string.Format("{0}|{1}", i, arg);
                    // 生成face的描述用于debug
                    faceDescriptions.Add(description);
                }
            }
        }

        /// <summary>
        /// 根据mesh上每个顶点属于哪几个polygon的关系列表，找到这个顶点对应的那几个polygon的中心的关系列表
        /// </summary>
        /// <param name="centerPoints"></param>
        /// <param name="vertexBelongtoWhichPolygon"></param>
        /// <returns></returns>
        public List<List<Point3d>> FetchPolylineCentersForDualVerticesLoL(List<Point3d> centerPoints, List<List<int>> vertexBelongtoWhichPolygon)
        {
            List<List<Point3d>> vertexBelongtoWhichPolygonCenter = new List<List<Point3d>>();
            for (int i = 0; i < vertexBelongtoWhichPolygon.Count; i++)
            {
                vertexBelongtoWhichPolygonCenter.Add(new List<Point3d>());
                foreach (int index in vertexBelongtoWhichPolygon[i])
                {
                    vertexBelongtoWhichPolygonCenter[i].Add(centerPoints[index]);
                }
            }

            return vertexBelongtoWhichPolygonCenter;
        }

        /// <summary>
        /// 获取每个多边形的中心
        /// </summary>
        /// <param name="convexPolygons"></param>
        /// <returns></returns>
        public List<Point3d> PolygonCenters(List<Polyline> convexPolygons)
        {

            if (convexPolygons == null)
            {
                return null;
            }
            else
            {
                List<Point3d> polygonCenterPoints = new List<Point3d>();
                foreach (Polyline polyline in convexPolygons)
                {
                    polygonCenterPoints.Add(polyline.CenterPoint());
                }

                return polygonCenterPoints;
            }
        }

        /// <summary>
        /// 每个点在哪些polygon里，输出每个点对应的polygon序号
        /// 每个vertex属于哪些三角形polygon
        /// </summary>
        /// <param name="points"></param>
        /// <param name="polygons"></param>
        /// <param name="basePlane"></param>
        /// <returns></returns>
        public List<List<int>> PointsInWhichPolylines(List<Point3d> points, List<Polyline> polygons, Plane basePlane)
        {
            // double num = 0.01 * RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
            // 每一支上的每个元素代表，这个点在哪个polygon里（polygon的序号）。
            List<List<int>> pointContainInWhichPolygon = new List<List<int>>();

            for (int i = 0; i < points.Count; i++)
            {
                pointContainInWhichPolygon.Add(new List<int>());

                for (int j = 0; j < polygons.Count; j++)
                {
                    Curve curve = polygons[j].ToNurbsCurve();
                    PointContainment pointContainment = curve.Contains(points[i], basePlane, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                    if (pointContainment == PointContainment.Coincident | pointContainment == PointContainment.Inside)
                    {
                        pointContainInWhichPolygon[i].Add(j);
                    }
                }
            }
            return pointContainInWhichPolygon;
        }


        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);

            for (int i = 0; i < DualConvexPolygon.Count; i++)
            {
                args.Display.DrawPolyline(DualConvexPolygon[i], Color.DarkOrange, Thickness);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualConvexVertices.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualConvexVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
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
            // 屏蔽掉电池原本的预览
            // base.DrawViewportMeshes(args);

            for (int i = 0; i < DualConvexPolygon.Count; i++)
            {
                args.Display.DrawPolyline(DualConvexPolygon[i], Color.DarkOrange, Thickness);
            }

            args.Display.DrawLines(GraphEdges, Color.ForestGreen, Thickness);

            for (int i = 0; i < DualConvexVertices.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualConvexVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
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
            get { return new Guid("72a758d7-2dcc-4b72-bd48-5bb5ff8bfde9"); }
        }
    }
}