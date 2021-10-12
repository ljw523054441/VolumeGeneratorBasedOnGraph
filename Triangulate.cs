using Grasshopper.Kernel;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
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

                List<Mesh> triangleMeshList = EnumerateAllTriangleMeshList(ClosedPolylineToTriangleMesh(convexFaceBorderPolylines));
                List<Mesh> list5 = new List<Mesh>();

                foreach (Mesh firstMesh in triangleMeshList)
                {
                    list5.Add(RegenerateMesh(nodePoints, firstMesh));
                }

                List<Mesh> list6;
                if (flag)
                {
                    list6 = RemoveTriangularDuals(list5);
                }
                else
                {
                    list6 = list5;
                }

                DA.SetDataList("AllTriangularMeshes", list6);
                DA.SetData("TheChosenTriangularMesh", list6[index]);
            }
        }

        public List<Mesh> EnumerateAllTriangleMeshList(List<List<Mesh>> allPolylineCorrespondTriangleMesh)
        {
            List<Mesh> result;

            if (allPolylineCorrespondTriangleMesh.Count <= 1)
            {
                result = null;
            }
            else
            {
                while (allPolylineCorrespondTriangleMesh.Count > 1)
                {
                    List<Mesh> list = new List<Mesh>();
                    List<Mesh> polyline_0_CorrespondTriangleMesh = allPolylineCorrespondTriangleMesh[0];
                    List<Mesh> polyline_1_CorrespondTriangleMesh = allPolylineCorrespondTriangleMesh[1];

                    foreach (Mesh polyline_0_TriangleMeshType in polyline_0_CorrespondTriangleMesh)
                    {
                        foreach (Mesh polyline_1_TriangleMeshType in polyline_1_CorrespondTriangleMesh)
                        {
                            Mesh mesh = new Mesh();
                            mesh.Append(polyline_0_TriangleMeshType);
                            mesh.Append(polyline_1_TriangleMeshType);
                            mesh.Vertices.CombineIdentical(true, true);
                            list.Add(mesh);
                        }
                    }
                    allPolylineCorrespondTriangleMesh.RemoveAt(0);
                    allPolylineCorrespondTriangleMesh.RemoveAt(0);
                    allPolylineCorrespondTriangleMesh.Insert(0, list);
                }
                
                result = allPolylineCorrespondTriangleMesh[0];
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
        /// 删除无效的三角形网格
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

        public Mesh RegenerateMesh(List<Point3d> vertices, Mesh firstMesh)
        {
            Point3d[] array = firstMesh.Vertices.ToPoint3dArray();
            Point3dList point3dList = new Point3dList(vertices);
            List<int> list = new List<int>();
            foreach (Point3d testPoint in array)
            {
                list.Add(point3dList.ClosestIndex(testPoint));
            }

            List<MeshFace> list2 = new List<MeshFace>();
            foreach (MeshFace meshFace in firstMesh.Faces)
            {
                MeshFace item = new MeshFace(list[meshFace.A], list[meshFace.B], list[meshFace.C]);
                list2.Add(item);
            }

            Mesh mesh = new Mesh();
            mesh.Vertices.AddVertices(vertices);
            mesh.Faces.AddFaces(list2);
            return mesh;
        }

        public List<Mesh> RemoveTriangularDuals(List<Mesh> TINs)
        {
            List<Mesh> list = new List<Mesh>();

            foreach (Mesh mesh in TINs)
            {
                List<int> list2 = new List<int>();
                for (int i = 0; i < mesh.TopologyVertices.Count; i++)
                {
                    list2.Add(Enumerable.Count<int>(mesh.TopologyVertices.ConnectedTopologyVertices(i)));
                }

                if (!list2.Contains(3))
                {
                    list.Add(mesh);
                }
            }
            return list;
        }

        /// <summary>
        /// 凸包类
        /// </summary>
        public class ConvexPolygon : List<int>
        {
            /// <summary>
            /// 用一个index的列表来存储和表示凸包ConvexPolygon，它不像Polyline那样需要在末尾再存上[0]
            /// </summary>
            private List<int> indices;

            /// <summary>
            /// Debug用
            /// </summary>
            public string Description
            {
                get
                {
                    return string.Join<int>(",", indices);
                }
            }

            /// <summary>
            /// 判断凸包ConvexPolygon是否已经是三角形了（三角形也是一种凸包）
            /// </summary>
            public bool IsTriangle
            {
                get
                {
                    return indices.Count == 3;
                }
            }

            /// <summary>
            /// 判断凸包ConvexPolygon是否已经是四边形了（四边形也是一种凸包）
            /// </summary>
            public bool IsQuadrangle
            {
                get
                {
                    return indices.Count == 4;
                }
            }

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="vertex_indices"></param>
            public ConvexPolygon(List<int> vertex_indices)
            {
                indices = new List<int>();
                indices = vertex_indices;
            }

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="PolygonalPolyline"></param>
            public ConvexPolygon(Polyline PolygonalPolyline)
            {
                if (PolygonalPolyline.IsClosed)
                {
                    indices = new List<int>();
                    // 因为闭合的polyline会有首尾重复的一个点
                    for (int i = 0; i < PolygonalPolyline.Count - 1; i++)
                    {
                        indices.Add(i);
                    }
                }
            }

            /// <summary>
            /// 将凸包ConvexPolygon转化为Polyline类，凸包ConvexPolygon只有点序号，通过输入实际的点，来形成闭合的Polyline
            /// </summary>
            /// <param name="Vertices"></param>
            /// <returns></returns>
            public Polyline ToPolyline(IEnumerable<Point3d> Vertices)
            {
                List<Point3d> list = new List<Point3d>();

                foreach (int num in indices)
                {
                    list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, num));
                }

                // 保证生成的Polyline的首尾闭合，将起始点又放在队尾
                list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, indices[0]));

                return new Polyline(list);
            }

            /// <summary>
            /// 将凸包ConvexPolygon三角形化
            /// </summary>
            /// <returns></returns>
            public List<List<MeshFace>> TriangulateConvexPolygon()
            {
                if (this == null)
                {
                    return null;
                }


                List<List<MeshFace>> LoLConvexPolygonMeshFace = new List<List<MeshFace>>();

                // 如果对象是三角化的
                if (IsTriangle)
                {
                    LoLConvexPolygonMeshFace.Add(new List<MeshFace>());

                    MeshFace currentConvexPolygon = new MeshFace(indices[0], indices[1], indices[2]);
                    LoLConvexPolygonMeshFace[0].Add(currentConvexPolygon);

                    return LoLConvexPolygonMeshFace;
                }
                // 如果对象不是三角化的
                else
                {
                    // 就从第二个顶点开始，对对象进行分割
                    for (int i = 2; i < indices.Count; i++)
                    {
                        List<ConvexPolygon> subConvexPolygonList = DivideConvexPolygon(i);

                        // 对于分割后生成的子对象列表来说，如果只有一个元素
                        if (subConvexPolygonList.Count == 1)
                        {
                            ConvexPolygon subConvexPolygon = subConvexPolygonList[0];
                            List<List<MeshFace>> LoLSubConvexPolygonMeshFace = subConvexPolygon.TriangulateConvexPolygon();
                            
                            foreach (List<MeshFace> meshFace in LoLSubConvexPolygonMeshFace)
                            {
                                List<MeshFace> list = new List<MeshFace>();
                                MeshFace currentConvexPolygonMeshFace = new MeshFace(indices[0], indices[1], indices[i]);

                                list.Add(currentConvexPolygonMeshFace);
                                list.AddRange(meshFace);

                                LoLConvexPolygonMeshFace.Add(list);
                            }
                        }
                        else
                        {
                            ConvexPolygon behindDivideConvexPolygon = subConvexPolygonList[0];
                            ConvexPolygon frontDivideConvexPolygon = subConvexPolygonList[1];

                            List<List<MeshFace>> LoLBehindDivideConvexPolygonMeshFace = behindDivideConvexPolygon.TriangulateConvexPolygon();
                            List<List<MeshFace>> LoLFrontDivideConvexPolygonMeshFace = frontDivideConvexPolygon.TriangulateConvexPolygon();

                            foreach (List<MeshFace> behindMeshFace in LoLBehindDivideConvexPolygonMeshFace)
                            {
                                foreach (List<MeshFace> frontMeshFace in LoLFrontDivideConvexPolygonMeshFace)
                                {
                                    List<MeshFace> list = new List<MeshFace>();
                                    MeshFace currentConvexPolygonMeshFace = new MeshFace(indices[0], indices[1], indices[i]);
                                    list.Add(currentConvexPolygonMeshFace);
                                    list.AddRange(behindMeshFace);
                                    list.AddRange(frontMeshFace);
                                    LoLConvexPolygonMeshFace.Add(list);
                                }
                            }
                        }
                    }
                }
                return LoLConvexPolygonMeshFace;
            }

            /// <summary>
            /// 在要分割位置的点处，分割整个大的Polygon，形成一个或两个子Polygon(behindConvexPolygon和frontConvexPolygon)
            /// </summary>
            /// <param name="divideIndex">要分割位置的点的序号</param>
            /// <returns></returns>
            public List<ConvexPolygon> DivideConvexPolygon(int divideIndex)
            {
                // 分割后得到的两个ConvexPolygon的列表
                List<ConvexPolygon> dividedConvexPolygonList;

                // 当分割位置的点的序号是0（表示polyline的第一个点），1（polyline的第二个点）或者indices.Count - 1（表示polyline的最后一个点）时，不能分割出两个新的ConvexPolygon
                if (divideIndex <= 1 | divideIndex > indices.Count - 1)
                {
                    dividedConvexPolygonList = null;
                }
                else
                {
                    // 从分割点往后，包括分割点和整个大的Polygon起点的这些点，所形成的子Polygon
                    List<int> behindDividedIndexList = new List<int>();
                    // 从分割点往前，包括分割点
                    List<int> frontDividedIndexList = new List<int>();

                    // 包括分割点并且往后
                    for (int i = divideIndex; i < indices.Count; i++)
                    {
                        behindDividedIndexList.Add(indices[i]);
                    }
                    // 包括整个大的Polygon起点
                    behindDividedIndexList.Add(indices[0]);

                    // 从第二个点开始，到分割点（包括分割点），所形成的的子Polygon
                    for (int i = 1; i <= divideIndex; i++)
                    {
                        frontDividedIndexList.Add(indices[i]);
                    }

                    ConvexPolygon behindDivideConvexPolygon = null;
                    ConvexPolygon frontDivideConvexPolygon = null;

                    // 大于三个点，才形成ConvexPolygon
                    if (behindDividedIndexList.Count >= 3)
                    {
                        // polygonFace -> behindDivideConvexPolygon
                        behindDivideConvexPolygon = new ConvexPolygon(behindDividedIndexList);
                    }
                    if (frontDividedIndexList.Count >= 3)
                    {
                        // polygonFace2 -> frontDivideConvexPolygon
                        frontDivideConvexPolygon = new ConvexPolygon(frontDividedIndexList);
                    }

                    // 输出最后分割出来符合顶点数>=3的子ConvexPolygon
                    if (behindDivideConvexPolygon != null & frontDivideConvexPolygon != null)
                    {
                        dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                        {
                            behindDivideConvexPolygon,
                            frontDivideConvexPolygon
                        });
                    }
                    else
                    {
                        if (behindDivideConvexPolygon == null)
                        {
                            dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                            {
                                frontDivideConvexPolygon
                            });
                        }
                        else if (frontDivideConvexPolygon == null)
                        {
                            dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                            {
                                behindDivideConvexPolygon
                            });
                        }
                        else
                        {
                            dividedConvexPolygonList = null;
                        }
                    }
                }
                return dividedConvexPolygonList;
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