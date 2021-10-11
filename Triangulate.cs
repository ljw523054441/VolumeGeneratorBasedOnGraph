using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph
{
    public class Triangulate : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Triangulate class.
        /// </summary>
        public Triangulate()
          : base("Triangulate", "Triangulate",
              "Finds all possible triangulations of a [convex] mesh",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddPointParameter("TutteOutputVertices", "NGV", "The set of locations of vertices as resulted from Tutte algorithm", GH_ParamAccess.list);
            pManager.AddCurveParameter("ConvexFaceBorders", "CFB", "Convex face borders of the Tutte algorithm as a list of polyline curves.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("IndexOfTriangulation", "I", "Index of a triangulation to be visualized from the list of all triangulations", GH_ParamAccess.item);
            pManager.AddBooleanParameter("ExcludeDegenerateTINS", "ExT", "排除那些产生三角形房间的三角剖分 Exclude those triangulations that give rise to triangular rooms", GH_ParamAccess.item);
            pManager[3].Optional = true;
            pManager[4].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("AllTriangularMeshes", "TMS", "All Computed triangulations; these triangulations describe all possible planar topologies for your plan diagram, the added links are those of adjacencies not connectivities", GH_ParamAccess.list);
            pManager.AddMeshParameter("TheChosenTriangularMesh", "TMI", "The one you have chosen with index I", GH_ParamAccess.item);
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

            List<Point3d> nodePoints = new List<Point3d>();           // list -> nodePoints
            List<Curve> convexFaceBorderCurves = new List<Curve>();              // list2 -> convexFaceBorderCurves
            List<Polyline> convexFaceBorderPolylines = new List<Polyline>();        // list3 -> convexFaceBorderPolylines

            int index = 0;
            bool flag = false;

            if (DA.GetDataList<Point3d>("TutteOutputVertices", nodePoints)
                & DA.GetDataList<Curve>("ConvexFaceBorders", convexFaceBorderCurves))
            {
                // 对输入的Curve类型的ConvexFaceBorder进行类型转换，转换成Curve类的子类Polyline
                foreach (Curve curve in convexFaceBorderCurves)
                {
                    Polyline polyline = null;
                    if (curve.TryGetPolyline(out polyline))
                    {
                        if (polyline.IsClosed)
                        {
                            convexFaceBorderPolylines.Add(polyline);
                        }
                    }
                }

                if (convexFaceBorderCurves.Count != convexFaceBorderPolylines.Count)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the ConvexFaceBorders are supposed to be closed polylines or polygons, but apparently at least one is not closed or is not a polyline!");
                }

                DA.GetData<int>("IndexOfTriangulation", ref index);
                DA.GetData<bool>("ExcludeDegenerateTINS", ref flag);

                List<Mesh> list4 = EnumerateAllTINs()
            }
        }

        public List<Mesh> EnumerateAllTINs(List<List<Mesh>> M)
        {
            List<Mesh> result;

            if (M.Count <= 1)
            {
                result = null;
            }
            else
            {
                while (M.Count > 1)
                {
                    List<Mesh> list = new List<Mesh>();
                    List<Mesh> list2 = M[0];
                    List<Mesh> list3 = M[1];

                    foreach (Mesh other in list2)
                    {
                        foreach (Mesh other2 in list3)
                        {
                            Mesh mesh = new Mesh();
                            mesh.Append(other);
                            mesh.Append(other2);
                            mesh.Vertices.CombineIdentical(true, true);
                            list.Add(mesh);
                        }
                    }
                    M.RemoveAt(0);
                    M.RemoveAt(0);
                    M.Insert(0, list);
                }
                
                result = M[0];
            }
            return result;
        }

        public List<List<Mesh>> ToTINs(List<Polyline> Polygons)
        {
            List<List<Mesh>> result;

            if (Polygons.Count <= 1)
            {
                result = null;
            }
            else
            {
                List<List<Mesh>> list = new List<List<Mesh>>();
                foreach (Polyline polyline in Polygons)
                {
                    Trangulate
                }
            }
        }


        public class PolygonFace : List<int>
        {
            // 用一个index的列表来存储和表示Polygon
            private List<int> indices;

            public string Description
            {
                get
                {
                    return string.Join<int>(",", indices);
                }
            }

            public bool IsTriangle
            {
                get
                {
                    return indices.Count == 3;
                }
            }

            public bool IsQuadrangle
            {
                get
                {
                    return indices.Count == 4;
                }
            }

            public PolygonFace(List<int> vertex_indices)
            {
                indices = vertex_indices;
            }

            public PolygonFace(Polyline PolygonalPolyline)
            {
                if (PolygonalPolyline.IsClosed)
                {
                    for (int i = 0; i < PolygonalPolyline.Count - 1; i++)
                    {
                        indices.Add(i);
                    }
                }
            }

            public Polyline ToPolyline(IEnumerable<Point3d> Vertices)
            {
                List<Point3d> list = new List<Point3d>();

                foreach (int num in indices)
                {
                    list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, num));
                }

                // 保证Polygon闭合
                list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, indices[0]));

                return new Polyline(list);
            }

            public List<List<MeshFace>> TriangulatePolygonFace()
            {
                //List<List<MeshFace>> result;
                if (this == null)
                {
                    return null;
                }

                List<List<MeshFace>> list = new List<List<MeshFace>>();
                if (this.IsTriangle)
                {
                    // List<List<MeshFace>> list2 = new List<List<MeshFace>>();
                    list.Add(new List<MeshFace>());
                    List<MeshFace> list3 = list[0];
                    MeshFace item = new MeshFace(indices[0], indices[1], indices[2]);
                    list3.Add(item);
                    //result = list;
                    return list;
                }

                for (int i = 0; i < indices.Count; i++)
                {
                    List<PolygonFace> list4 = 
                }
            }

            public List<PolygonFace> DividePolygonFace(int ear)
            {
                if (ear <= 1 | ear > indices.Count - 1)
                {
                    return null;
                }

                List<int> list = new List<int>();
                List<int> list2 = new List<int>();
                for (int i = ear; i < indices.Count; i++)
                {
                    list.Add(indices[i]);
                }
                list.Add(indices[0]);
                for (int i = 1; i <= ear; i++)
                {
                    list2.Add(indices[i]);
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
            get { return new Guid("88a4820f-ef86-4c32-ba34-74e80870373f"); }
        }
    }
}