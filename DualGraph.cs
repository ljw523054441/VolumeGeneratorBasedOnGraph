using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph
{
    public class DualGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DualGraph class.
        /// </summary>
        public DualGraph()
          : base("DualGraph", "DualGraph",
              "Find corresponding dual graph",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            DualVertices = new List<Point3d>();
            DualConvexPolygon = new List<List<Point3d>>();
            DualVertexTextDots = new List<TextDot>();
    }

        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualConvexPolygon;
        private List<TextDot> DualVertexTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddMeshParameter("TriangleMesh", "TriangleMesh", "A completely trianngulated mesh based on the NEWS graph convex drawing", GH_ParamAccess.item);
            pManager.AddGenericParameter("SubAttributes", "SubAttributes", "(Optional)Put labels, areas and colors as an attribute LOL, if you want to get a legned for your analysis.", GH_ParamAccess.list);
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("DualVertices", "DV", "Vertices of the dual mesh", GH_ParamAccess.list);
            pManager.AddTextParameter("DualFaces", "DF", "Faces of the dual mesh as a list of General_Graph meshface", GH_ParamAccess.list);
            pManager.AddCurveParameter("DualConvexPolygon", "DFB", "Face borders of the dual mesh as a list of closed polyline curves", GH_ParamAccess.list);
            // pManager.AddGenericParameter("DualFaceColors", "DFC", "Face colors as a list of materials based on the colors ", GH_ParamAccess.list);
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
            List<NodeAttribute> list = new List<NodeAttribute>();

            if (DA.GetData<Mesh>("TriangleMesh",ref mesh))
            {

                // 获取网格上每个面的闭合polyline，组成列表List<Polyline> facePolylines
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

                List<Point3d> meshVerticesPointList = new List<Point3d>();

                // meshVerticesPointList中，是分innerNode和outerNode的，顺序跟nodePoint原来的顺序一样
                // 所以在计算是要把outerNode去掉进行计算
                for (int i = 0; i < meshVertices.Count - boundaryNodeCount; i++)
                {
                    meshVerticesPointList.Add(meshVertices[i]);
                }

                List<Point3d> centerPoints = new List<Point3d>();
                // List<List<Point3d>> list4 = new List<List<Point3d>>();
                List<List<int>> sortedCenterPointIndexLoL = new List<List<int>>();
                List<Polyline> dualConvexPolygon = new List<Polyline>();
                List<string> faceDescriptions = new List<string>();

                SortAdjacentCenterPoints(meshVerticesPointList, facePolylines, ref centerPoints, ref sortedCenterPointIndexLoL, ref dualConvexPolygon, ref faceDescriptions);

                DA.SetDataList("DualVertices", centerPoints);
                DA.SetDataList("DualFaces", faceDescriptions);
                DA.SetDataList("DualFaceBorders", dualConvexPolygon);

                DualVertices.Clear();
                DualConvexPolygon.Clear();
                DualVertexTextDots.Clear();

                DualVertices = centerPoints;
                for (int i = 0; i < dualConvexPolygon.Count; i++)
                {
                    DualConvexPolygon.Add(new List<Point3d>());
                    DualConvexPolygon[i].AddRange(dualConvexPolygon[i]);
                }

                for (int i = 0; i < DualVertices.Count; i++)
                {
                    TextDot textDot = new TextDot(string.Format("DualVertex {0}", i), DualVertices[i]);
                    DualVertexTextDots.Add(textDot);
                }



                if (DA.GetDataList<NodeAttribute>("SubAttributes", list))
                {
                    // List<Color> list7 = (List<Color>)list[2].Value;
                }


            }
        }

        /// <summary>
        /// 按照每个centerPoint向量与x轴正方向的夹角的弧度值大小，对每个顶点所邻接的dualConvexPolygon的中心进行排序
        /// </summary>
        /// <param name="vertices"></param>
        /// <param name="facePolylines"></param>
        /// <param name="facePolylineCenterPoints"></param>
        /// <param name="sortedCenterPointIndexLoL"></param>
        /// <param name="dualConvexPolygons"></param>
        /// <param name="faceDescriptions"></param>
        public void SortAdjacentCenterPoints(List<Point3d> vertices, List<Polyline> facePolylines, ref List<Point3d> facePolylineCenterPoints, ref List<List<int>> sortedCenterPointIndexLoL, ref List<Polyline> dualConvexPolygons, ref List<string> faceDescriptions)
        {
            facePolylineCenterPoints = PolygonCenters(facePolylines);
            List<List<int>> verticeContainInWhichFacePolyline = PointsInPolygons(vertices, facePolylines, Plane.WorldXY);

            List<List<Point3d>> vertexAdjacentWhichPolygonCenter = FetchPolygonCentersForDualVerticesLoL(facePolylineCenterPoints, verticeContainInWhichFacePolyline);

            

            if (vertexAdjacentWhichPolygonCenter == null | vertices == null |vertices.Count != vertexAdjacentWhichPolygonCenter.Count)
            {
                return;
            }
            else
            {
                List<List<Point3d>> sortedAdjacentCenterLoL = new List<List<Point3d>>();
                List<double> radianList = new List<double>();

                // 对于mesh上每一个顶点
                for (int i = 0; i < vertexAdjacentWhichPolygonCenter.Count; i++)
                {
                    radianList.Clear();
                    sortedAdjacentCenterLoL.Add(new List<Point3d>());
                    sortedCenterPointIndexLoL.Add(new List<int>());

                    // 每个顶点会对应几个PolygonCenter
                    int[] centerPointIndexArray = new int[vertexAdjacentWhichPolygonCenter[i].Count];

                    // 对于每一支上的每个元素（即PolygonCenter）
                    for (int j = 0; j < vertexAdjacentWhichPolygonCenter[i].Count; j++)
                    {
                        // 在centerPoint列表中找到与vertex邻接的对应的centerpoint的序号
                        centerPointIndexArray[j] = facePolylineCenterPoints.IndexOf(vertexAdjacentWhichPolygonCenter[i][j]);
                        /* atan2(y, x)是4象限反正切，它的取值不仅取决于正切值y/x，还取决于点 (x, y) 落入哪个象限：
                         * 当点(x, y) 落入第一象限时，atan2(y, x)的范围是 0 ~ pi/2;
                         * 当点(x, y) 落入第二象限时，atan2(y, x)的范围是 pi/2 ~ pi;
                         * 当点(x, y) 落入第三象限时，atan2(y, x)的范围是 －pi～－pi/2;
                         * 当点(x, y) 落入第四象限时，atan2(y, x)的范围是 -pi/2～0. */
                        double radian = Math.Atan2(vertexAdjacentWhichPolygonCenter[i][j].Y - vertices[i].Y, vertexAdjacentWhichPolygonCenter[i][j].X - vertices[i].X);
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
                    Point3d[] AdjacentCenterArray = vertexAdjacentWhichPolygonCenter[i].ToArray();
                    Array.Sort<double, Point3d>(radianList.ToArray(), AdjacentCenterArray);
                    Array.Sort<double, int>(radianList.ToArray(), centerPointIndexArray);

                    sortedAdjacentCenterLoL[i].AddRange(AdjacentCenterArray);
                    sortedCenterPointIndexLoL[i].AddRange(centerPointIndexArray);

                    List<Point3d> pointsForPolyline = new List<Point3d>();
                    pointsForPolyline.AddRange(AdjacentCenterArray);
                    pointsForPolyline.Add(AdjacentCenterArray[0]);
                    Polyline dualConvexPolyline = new Polyline(pointsForPolyline);
                    dualConvexPolygons.Add(dualConvexPolyline);

                    string arg = string.Join<int>(";", sortedCenterPointIndexLoL[i]);
                    string item2 = string.Format("{0}|{1}", sortedCenterPointIndexLoL[i].Count, arg);
                    faceDescriptions.Add(item2);
                }
            }
        }

        /// <summary>
        /// 根据mesh上每个顶点属于哪几个polygon的关系列表，找到这个顶点对应的那几个polygon的中心的关系列表
        /// </summary>
        /// <param name="centerPoints"></param>
        /// <param name="vertexBelongtoWhichPolygon"></param>
        /// <returns></returns>
        public List<List<Point3d>> FetchPolygonCentersForDualVerticesLoL(List<Point3d> centerPoints, List<List<int>> vertexBelongtoWhichPolygon)
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
        public List<List<int>> PointsInPolygons(List<Point3d> points, List<Polyline> polygons, Plane basePlane)
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

            for (int i = 0; i < DualVertices.Count; i++)
            {
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
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

            for (int i = 0; i < DualVertices.Count; i++)
            {
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
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