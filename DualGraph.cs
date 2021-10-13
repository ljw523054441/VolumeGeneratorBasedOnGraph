using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using System;
using System.Collections.Generic;
using System.Drawing;

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
        }

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
            pManager.AddCurveParameter("DualFaceBorders", "DFB", "Face borders of the dual mesh as a list of closed polyline curves", GH_ParamAccess.list);
            pManager.AddGenericParameter("DualFaceColors", "DFC", "Face colors as a list of materials based on the colors ", GH_ParamAccess.list);
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
                // 注意这个减4
                for (int i = 0; i < meshVertices.Count - 4; i++)
                {
                    meshVerticesPointList.Add(meshVertices[i]);
                }

                List<Point3d> centerPoints = new List<Point3d>();
                List<List<Point3d>> list4 = new List<List<Point3d>>();
                List<List<int>> sortedIndices = new List<List<int>>();
                List<Polyline> polygons = new List<Polyline>();
                List<string> faceDescriptions = new List<string>();

                PolarSortNeighbors(meshVerticesPointList, facePolylines, ref centerPoints, ref sortedIndices, ref polygons, ref faceDescriptions);

                DA.SetDataList("DualVertices", centerPoints);
                DA.SetDataList("DualFaces", faceDescriptions);
                DA.SetDataList("DualFaceBorders", polygons);

                if (DA.GetDataList<NodeAttribute>("SubAttributes", list))
                {
                    // List<Color> list7 = (List<Color>)list[2].Value;
                }


            }
        }

        public void PolarSortNeighbors(List<Point3d> vertices, List<Polyline> facePolylines, ref List<Point3d> facePolylineCenterPoints, ref List<List<int>> sortedIndices, ref List<Polyline> polygons, ref List<string> faceDescriptions)
        {
            facePolylineCenterPoints = PolygonCenters(facePolylines);
            List<List<int>> verticeContainInWhichFacePolyline = PointsInPolygons(vertices, facePolylines, Plane.WorldXY);

            List<List<Point3d>> vertexBelongtoWhichPolygonCenter = FetchPolygonCentersForDualVerticesLoL(facePolylineCenterPoints, verticeContainInWhichFacePolyline);

            

            if (vertexBelongtoWhichPolygonCenter == null | vertices == null |vertices.Count != vertexBelongtoWhichPolygonCenter.Count)
            {
                return;
            }
            else
            {
                List<List<Point3d>> list2 = new List<List<Point3d>>();
                List<double> list3 = new List<double>();

                // 对于mesh上每一个顶点
                for (int i = 0; i < vertexBelongtoWhichPolygonCenter.Count; i++)
                {
                    list3.Clear();
                    list2.Add(new List<Point3d>());
                    sortedIndices.Add(new List<int>());

                    // 每个顶点会对应几个PolygonCenter
                    int[] array = new int[vertexBelongtoWhichPolygonCenter[i].Count];

                    // 对于每一支上的每个元素（即PolygonCenter）
                    for (int j = 0; j < vertexBelongtoWhichPolygonCenter[i].Count; j++)
                    {
                        array[j] = facePolylineCenterPoints.IndexOf(vertexBelongtoWhichPolygonCenter[i][j]);

                        double num10 = Math.Atan2(vertexBelongtoWhichPolygonCenter[i][j].Y - vertices[i].Y, vertexBelongtoWhichPolygonCenter[i][j].X - vertices[i].X);
                        if (num10 >= 0.0)
                        {
                            list3.Add(num10);
                        }
                        else
                        {
                            list3.Add(2 * Math.PI + num10);
                        }
                    }

                    Point3d[] array2 = vertexBelongtoWhichPolygonCenter[i].ToArray();
                    Array.Sort<double, Point3d>(list3.ToArray(), array2);
                    Array.Sort<double, int>(list3.ToArray(), array);
                    list2[i].AddRange(array2);
                    sortedIndices[i].AddRange(array);
                    List<Point3d> list4 = new List<Point3d>();
                    list4.AddRange(array2);
                    list4.Add(array2[0]);
                    Polyline item = new Polyline(list4);
                    polygons.Add(item);
                    string arg = string.Join<int>(";", sortedIndices[i]);
                    string item2 = string.Format("{0}|{1}", sortedIndices[i].Count, arg);
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
        /// <param name="polygon"></param>
        /// <returns></returns>
        public List<Point3d> PolygonCenters(List<Polyline> polygon)
        {
            
            if (polygon == null)
            {
                return null;
            }
            else
            {
                List<Point3d> polygonCenterPoints = new List<Point3d>();
                foreach (Polyline polyline in polygon)
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