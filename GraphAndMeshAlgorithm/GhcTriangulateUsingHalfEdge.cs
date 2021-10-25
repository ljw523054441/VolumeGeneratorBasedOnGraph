using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Plankton;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcTriangulateUsingHalfEdge : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcTriangulateUsingHalfEdge class.
        /// </summary>
        public GhcTriangulateUsingHalfEdge()
          : base("TriangulateUsingHalfEdge", "Nickname",
              "Description",
              "Category", "Subcategory")
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
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);

            // pManager.AddPointParameter("TutteOutputVertices", "NGV", "The set of locations of vertices as resulted from Tutte algorithm", GH_ParamAccess.list);
            pManager.AddCurveParameter("ConvexFaceBorders", "CFBorders", "Convex face borders of the Tutte algorithm as a list of polyline curves.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("IndexOfTriangularMesh", "I", "Index of a triangulation to be visualized from the list of all triangulations", GH_ParamAccess.item);
            pManager.AddBooleanParameter("ExcludeDegenerateTINS", "ExT", "排除那些产生三角形房间的三角剖分 Exclude those triangulations that give rise to triangular rooms", GH_ParamAccess.item);
            pManager[3].Optional = true;
            pManager[4].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("AllTriangularMeshes", "AllTMesh", "所有可能的三角形剖分结果。All Computed triangulations; these triangulations describe all possible planar topologies for your plan diagram, the added links are those of adjacencies not connectivities", GH_ParamAccess.list);
            pManager.AddMeshParameter("TheChosenTriangularMesh", "TMesh", "所选择的那个三角形剖分结果。The one you have chosen with index I", GH_ParamAccess.item);
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
            int innerNodeCount = globalParameter.VolumeNodeCount;
            int outerNodeCount = globalParameter.BoundaryNodeCount;

            List<Node> nodes = new List<Node>();
            List<Point3d> nodePoints = new List<Point3d>();
            List<Point3d> outerPoints = new List<Point3d>();

            List<Curve> convexFaceBorderCurves = new List<Curve>();
            List<Polyline> convexFaceBorderPolylines = new List<Polyline>();

            int index = 0;
            bool flag = false;

            Thickness = 2;


            if (DA.GetDataList<Node>("GraphNode", nodes)
                && DA.GetDataList<Curve>("ConvexFaceBorders", convexFaceBorderCurves))
            {
                // 从Node类中获取点的坐标
                for (int i = 0; i < nodes.Count; i++)
                {
                    nodePoints.Add(nodes[i].NodeVertex);
                }
                for (int i = 0; i < nodes.Count - innerNodeCount; i++)
                {
                    outerPoints.Add(nodes[i + innerNodeCount].NodeVertex);
                }


                ConvexPolylinesPoints.Clear();
                SelectedTriangleMeshEdges.Clear();
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();

                // 设置用于可视化的TextDot
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (i < innerNodeCount)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} {1}", i, nodes[i].NodeAttribute.NodeLabel), nodes[i].NodeVertex);
                        InnerNodeTextDot.Add(textDot);
                    }
                    else
                    {
                        TextDot textDot = new TextDot(string.Format("{0} {1}", i, nodes[i].NodeAttribute.NodeLabel), nodes[i].NodeVertex);
                        OuterNodeTextDot.Add(textDot);
                    }
                }

                // 对输入的Curve类型的ConvexFaceBorder进行类型转换，转换成Curve类的子类Polyline
                for (int i = 0; i < convexFaceBorderCurves.Count; i++)
                {
                    Polyline polyline = null;
                    if (convexFaceBorderCurves[i].TryGetPolyline(out polyline))
                    {
                        if (polyline.IsClosed)
                        {
                            convexFaceBorderPolylines.Add(polyline);

                            ConvexPolylinesPoints.Add(new List<Point3d>());
                            ConvexPolylinesPoints[i].AddRange(polyline);
                            // ConvexPolylinesPoints[i].RemoveAt(ConvexPolylinesPoints[i].Count - 1);
                        }
                    }
                }

                if (convexFaceBorderCurves.Count != convexFaceBorderPolylines.Count)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "ConvexFaceBorders应该都是闭合的polyline");
                }

                DA.GetData<int>("IndexOfTriangularMesh", ref index);
                DA.GetData<bool>("ExcludeDegenerateTINS", ref flag);




            }


        }


        List<List<PlanktonMesh>> ClosedPolylineTriangleMesh(List<Polyline> closedPolylines)
        {
            if (closedPolylines.Count <= 1)
            {
                return null;
            }
            else
            {
                List<List<PlanktonMesh>> allPolylineCorrespondTriangleMesh = new List<List<PlanktonMesh>>();
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
            get { return new Guid("b790ac7d-ff2e-49a7-91fb-20bf3f7b4751"); }
        }
    }
}