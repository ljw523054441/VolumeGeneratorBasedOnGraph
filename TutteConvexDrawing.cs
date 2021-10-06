using Grasshopper;
using Grasshopper.GUI.Gradient;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;

namespace VolumeGeneratorBasedOnGraph
{
    public class TutteConvexDrawing : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the TutteConvexDrawing class.
        /// </summary>
        public TutteConvexDrawing()
          : base("TutteConvexDrawing", "TutteDrawing",
              "Finds a unique convex drawing of a bi-connected graph. A bi-connected graph is a graph in which every vertex(node) is connected to others at least through two edges (links)",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPlaneParameter("BasePlane", "BasePlane", "A plane on which you want to get the Tutte convex drawing", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddIntegerParameter("Graph", "Graph", "A graph describing which node is connected to which other nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("NodePoints", "Vertices", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            pManager.AddGenericParameter("NodeAttributes", "Attributes", "图结构中节点所包含的属性", GH_ParamAccess.list);
            pManager[3].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("New Graph Vertices", "New Graph Vertices", "A set of new vertices on which a neat coin/kissing disk drawing can be drawn", GH_ParamAccess.list);
            pManager.AddCurveParameter("ConvexFaceBorders", "ConvexFaceBorders", "A list of convex polylines represening the boundaries of Tutte convex faces.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("SubGraph", "SubGraph", "A graph describing the relations of main nodes only, without NEWS nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("SubVertices", "SubVertices", "A list of points containing inner vertices only, excluding NEWS vertices", GH_ParamAccess.list);
            pManager.AddGenericParameter("SubAttributes", "SubAttributes", "A list of lists of attributes pertaining to the main vertices only", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Integer> gh_Structure_graph = null;
            DataTree<int> graph = new DataTree<int>();

            List<Point3d> nodePointList = new List<Point3d>();                   // list -> nodePointList
            List<GH_ObjectWrapper> list2 = new List<GH_ObjectWrapper>();// list2 ->
            Plane worldXY = Plane.WorldXY;
            List<string> list3 = new List<string>();                    // list3 ->
            List<double> list4 = new List<double>();                    // list4 ->
            List<Color> list5 = new List<Color>();                      // list5 ->

            // bool flag = DA.GetDataTree<GH_Integer>("Graph", out gh_Structure) & DA.GetDataList<Point3d>("Vertices", list);
            if (DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_graph) & DA.GetDataList<Point3d>("NodePoints", nodePointList))
            {
                DA.GetData<Plane>("BasePlane", ref worldXY);

                // List<List<int>> list6 = new List<List<int>>();          // list6 -> graph

                //foreach (GH_Path gh_Path in gh_Structure_graph.Paths)
                //{
                //    List<int> list7 = new List<int>();
                //    foreach (object obj in gh_Structure_graph.get_Branch(gh_Path))
                //    {
                //        GH_Integer gH_Integer = 
                //    }
                //}

                // // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_graph, ref graph);

                List<int> volumeNodeIndexList = new List<int>();                      // list8 -> volumeNodeIndexList
                int num = nodePointList.Count - 5;                               // num ->
                
                // int num3 = nodePointList.Count - 5;
                // int num4 = 0;

                for (int i = 0; i < nodePointList.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    volumeNodeIndexList.Add(i);
                }

                List<int> list9 = new List<int>();                      // list9 ->
                List<int> graphBranchCountList = new List<int>();                     // list10 -> graphBranchCountList
                // int num7 = 0;
                // int num8 = graph.DataCount - 1;
                // int num9 = 0;

                for (int i = 0; i < graph.BranchCount; i++)
                {
                    graphBranchCountList.Add(graph.Branch(i).Count);
                }

                List<int> list11 = new List<int>();                     // list11 ->
                // int num11 = 0;
                int num12 = num;
                int num13 = 0;
                for (int i = 0; i < nodePointList.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    list11.Add(graphBranchCountList[volumeNodeIndexList[i]]);
                }

                List<List<int>> list12 = new List<List<int>>();         // list12 ->
                int num15 = 0;
                int num16 = num;
                num13 = num15;
                for (int i = 0; i < nodePointList.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    list12.Add(new List<int>());
                    list12[i].AddRange(graph.Branch(volumeNodeIndexList[i]));
                }

                List<int> list13 = new List<int>();                     // list13 ->
                int num18 = nodePointList.Count - 4;
                int num19 = nodePointList.Count - 1;
                int num20 = num18;
                for (int i = nodePointList.Count - NodeAttribute.BoundaryNodeCount; i < nodePointList.Count; i++)
                {
                    list13.Add(i);
                }

                List<int> list14 = new List<int>();                     // list14 ->
                List<List<int>> list15 = new List<List<int>>();         // list15 ->
                int num22 = list12.Count - 1;
                int num23 = 0;
                int num24 = num22;
                //num9 = num23;
                for (int i = 0; i < list12.Count; i++)
                {
                    list15.Add(new List<int>());
                    int num26 = 0;
                    int num27 = list12[i].Count - 1;
                    for (int j = 0; j < list12[j].Count; j++)
                    {
                        if (list13.Contains(list12[i][j]))
                        {
                            list15[i].Add(list12[i][j]);
                        }
                    }
                }

                List<List<double>> list16 = new List<List<double>>();   // list16 ->
                List<List<double>> list17 = new List<List<double>>();   // list17 ->

                int num29 = 0;
                int num30 = num22;
                // num9 = num29;
                Point3d item;
                for (int i = 0; i < list12.Count; i++)
                {
                    list16.Add(new List<double>());
                    list17.Add(new List<double>());
                    int num32 = 0;
                    int num33 = list15[i].Count - 1;
                    num13 = num32;
                    for (int j = 0; j < list12.Count; j++)
                    {
                        List<double> list18 = list16[i];
                        item = nodePointList[list15[i][j]];
                        list18.Add(item.X);
                        List<double> list19 = list17[i];
                        item = nodePointList[list15[i][j]];
                        list19.Add(item.Y);
                    }
                }

                List<double> list20 = new List<double>();
                List<double> list21 = new List<double>();
                Matrix matrix = new Matrix(list12.Count, 1);
                Matrix matrix2 = new Matrix(list12.Count, 1);
                int num35 = 0;
                int num36 = num22;
                num13 = num35;
                for (int i = 0; i < list12.Count; i++)
                {
                    list20.Add(Sum(list16[i]));
                    matrix[i, 0] = Sum(list16[i]);
                    list21.Add(Sum(list17[i]));
                    matrix2[i, 0] = Sum(list17[i]);
                }

                DataTree<int> dataTree = new DataTree<int>();
                int num38 = 0;
                int num39 = nodePointList.Count - 5;
                // num9 = num38;
                for (int i = 0; i < nodePointList.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    dataTree.EnsurePath(new int[]
                    {
                        i
                    });
                    int num41 = 0;
                    int num42 = nodePointList.Count - 5;
                    num13 = num41;
                    for (int j = 0; j < nodePointList.Count - NodeAttribute.BoundaryNodeCount; j++)
                    {
                        if (SubGraph(GraphMatrix(graph), volumeNodeIndexList)[i][j] == 1)
                        {
                            dataTree.Branch(i).Add(j);
                        }
                    }
                }

                List<Point3d> list22 = new List<Point3d>();
                DA.SetDataTree(2, dataTree);

                Matrix matrix3 = ToGHMatrix(In_LaplacianM(SubGraph(GraphMatrix(graph), volumeNodeIndexList), list11));
                matrix3.Invert(0.0);
                Matrix matrix4 = matrix3;
                Matrix matrix5 = matrix4 * matrix;
                Matrix matrix6 = matrix4 * matrix2;
                List<Point3d> list23 = new List<Point3d>();

                int num44 = 0;
                int num45 = list11.Count - 1;
                num13 = num44;
                for (int i = 0; i < list11.Count; i++)
                {
                    List<Point3d> list24 = list23;
                    item = new Point3d(matrix5[i, 0], matrix6[i, 0], 0.0);
                    list24.Add(item);
                }

                int num47 = 0;
                int num48 = volumeNodeIndexList.Count - 1;
                num13 = num47;
                for (int i = 0; i < volumeNodeIndexList.Count; i++)
                {
                    nodePointList[volumeNodeIndexList[i]] = list23[i];
                }

                List<Point3d> list25 = UtilityFunctions.Relocate(nodePointList, Plane.WorldXY, worldXY);
                DA.SetDataList("New Graph Vertices", list25);
                DA.SetDataList("ConvexFaceBorders", GMeshFaceBoundaries(list25, GraphEdgeList(graph, list25), worldXY));


                for (int i = 0; i < nodePointList.Count; i++)
                {
                    list22.Add(nodePointList[i]);
                }
                DA.SetDataList("SubVertices", list22);
                if (DA.GetDataList("NodeAttributes", list2))
                {
                    List<NodeAttribute> list26 = ;
                    DA.SetDataList("SubAttributes", list26);
                }
            }

        }

        public double Sum(List<double> x)
        {
            //int num = 0;
            //int num2 = x.Count - 1;
            //int num3 = num;
            double sum = 0;
            for (int i = 0; i < x.Count - 1; i++)
            {
                sum += x[i];
            }
            return sum;
        }

        public List<NodeAttribute> SubAttributes(List<NodeAttribute> attributes, List<int> subIndices)
        {
            List<NodeAttribute> list = new List<NodeAttribute>();





            return list;
        }

        public List<List<int>> SubGraph(List<List<int>> graph, List<int> indices)
        {
            List<List<int>> list = new List<List<int>>();
            List<List<int>> list2 = new List<List<int>>();

            // int graphDataCount = graph.Count - 1;          // num -> graphDataCount
            // int num2 = 0;
            // int indexListCount = indices.Count - 1;                   // num3 -> indexListCount
            // int num4 = num2;
            for (int i = 0; i < indices.Count; i++)
            {
                list.Add(new List<int>());
                list[i].AddRange(graph[indices[i]]);
            }

            // int num7 = 0;
            // int num8 = indices.Count - 1;
            // num4 = num7;
            for (int i = 0; i < indices.Count; i++)
            {
                list2.Add(new List<int>());
                // int num10 = 0;
                // int num11 = indices.Count - 1;
                // int num12 = num10;
                for (int j = 0; j < indices.Count; j++)
                {
                    list2[i].Add(list[i][indices[j]]);
                }
            }

            return list2;
        }

        public List<List<int>> GraphMatrix(List<List<int>> graph)
        {
            List<List<int>> list = new List<List<int>>();
            
            //int num = graph.Count - 1;
            //int num2 = 0;
            //int num3 = num;
            //int num4 = num2;
            for (int i = 0; i < graph.Count; i++)
            {
                list.Add(new List<int>());
                //int num7 = 0;
                //int num8 = num;
                //int num9 = num7;
                for (int j = 0; j < graph.Count; j++)
                {
                    if (graph[i].Contains(j))
                    {
                        list[i].Add(1);
                    }
                    else
                    {
                        list[i].Add(0);
                    }
                }
            }

            return list;
        }

        public Matrix ToGHMatrix(List<List<int>> x)
        {
            Matrix result;
            if (x == null)
            {
                result = null;
            }
            else
            {
                int num = 1;
                int num2 = x.Count - 1;
                int num3 = num;
                for (int i = 0; i < x.Count; i++)
                {
                    // goto Block_4;

                    if (x[i].Count != x[0].Count)
                    {
                        break;
                    }

                }
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Matrix must be rectagular! Your input graph is not valid");
                // return null;
            // Block_4:
                Matrix matrix = new Matrix(x[0].Count, x.Count);
                int num6 = 0;
                int num7 = x[0].Count - 1;
                int num8 = num6;
                for (int i = 0; i < x[0].Count; i++)
                {
                    for (int j = 0; j < x[0].Count; j++)
                    {
                        matrix[i, j] = (double)x[j][i];
                    }
                }
                result = matrix;
            }

            return result;
        }

        public List<List<int>> In_LaplacianM(List<List<int>> graph, List<int> Degrees)
        {
            List<List<int>> list = new List<List<int>>();

            int num = graph.Count - 1;
            int num2 = 0;
            int num3 = num;
            int num4 = num2;
            for (int i = 0; i < graph.Count; i++)
            {
                list.Add(new List<int>());
                int num7 = 0;
                int num8 = num;
                int num9 = num7;
                for (int j = 0; j < graph.Count; j++)
                {
                    if (j == i)
                    {
                        list[i].Add(Degrees[j]);
                    }
                    else
                    {
                        list[i].Add(0 - graph[i][j]);
                    }
                }
            }

            return list;
        }

        //public List<Point3d> Relocate(List<Point3d> p, Plane location1, Plane location2)
        //{
        //    Point3d point3d = new Point3d(0.0, 0.0, 0.0);
        //    List<Point3d> list = new List<Point3d>();
        //    foreach (Point3d point3d2 in p)
        //    {
        //        point3d += point3d2;
        //    }
        //    point3d /= (double)p.Count;
        //    location1.Origin = point3d;
        //    foreach (Point3d item in p)
        //    {
        //        item.Transform(Transform.PlaneToPlane(location1, location2));
        //        list.Add(item);
        //    }
        //    return list;
        //}

        public Curve ConvexBorder(List<Point3d> vertices, Plane basePlane)
        {
            Point3dList point3dList = new Point3dList(vertices);
            point3dList.Sort();
            Plane plane = basePlane;
            // 旋转 PI / 4
            plane.Rotate(0.7853981633974483, plane.ZAxis);
            Rectangle3d rectangle3d = new Rectangle3d(plane, point3dList.First, point3dList.Last);
            return rectangle3d.ToNurbsCurve();
        }

        public List<Polyline> GMeshFaceBoundaries(List<Point3d> vertices, List<Line> edges, Plane basePlane)
        {
            List<Curve> list = new List<Curve>();
            foreach (Line line in edges)
            {
                list.Add(line.ToNurbsCurve());
            }

            Curve curve = ConvexBorder(vertices, basePlane);
            Brep brep = Brep.CreatePlanarBreps(curve, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)[0];
            BrepFace brepFace = brep.Faces[0];
            Brep brep2 = brepFace.Split(list, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            List<Polyline> list2 = new List<Polyline>();

            foreach (BrepFace brepFace2 in brep2.Faces)
            {
                Brep brep3 = brepFace2.DuplicateFace(true);
                Point3d[] array = brep3.DuplicateVertices();
                Curve curve2 = Curve.JoinCurves(brep3.Faces[0].DuplicateFace(true).DuplicateEdgeCurves())[0];
                Polyline item = null;
                curve2.TryGetPolyline(out item);
                list2.Add(item);
            }

            return list2;
        }

        public List<Line> GraphEdgeList(List<List<int>> graph, List<Point3d> graphVertices)
        {
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
            get { return new Guid("97181a6b-fd35-4eb5-9938-4e7a7e7b91ab"); }
        }
    }
}