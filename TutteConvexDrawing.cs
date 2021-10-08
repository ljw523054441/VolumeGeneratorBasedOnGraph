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
            // pManager.AddGenericParameter("SubAttributes", "SubAttributes", "A list of lists of attributes pertaining to the main vertices only", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Integer> gh_Structure_graph = null;
            DataTree<int> graph = new DataTree<int>();

            List<Point3d> nodePoints = new List<Point3d>();                     // list -> nodePointList

            List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();
            // List<GH_ObjectWrapper> list2 = new List<GH_ObjectWrapper>();        // list2 -> nodeAttributes

            Plane worldXY = Plane.WorldXY;
            List<string> nodeLabels = new List<string>();                    // list3 -> nodeLabels
            List<double> nodeAreas = new List<double>();                    // list4 -> nodeAreas
            List<Color> nodeMatColor = new List<Color>();                      // list5 -> nodeMatColor

            // bool flag = DA.GetDataTree<GH_Integer>("Graph", out gh_Structure) & DA.GetDataList<Point3d>("Vertices", list);
            if (DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_graph) & DA.GetDataList<Point3d>("NodePoints", nodePoints))
            {
                DA.GetData<Plane>("BasePlane", ref worldXY);

                // List<List<int>> list6 = new List<List<int>>();          // list6 -> graph

                // // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_graph, ref graph);

                List<int> volumeNodeIndexList = new List<int>();                      // list8 -> volumeNodeIndexList

                for (int i = 0; i < nodePoints.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    volumeNodeIndexList.Add(i);
                }

                // List<int> list9 = new List<int>();                                      // list9 ->
                List<int> graphBranchCountList = new List<int>();                       // list10 -> graphBranchCountList
                for (int i = 0; i < graph.BranchCount; i++)
                {
                    graphBranchCountList.Add(graph.Branch(i).Count);
                }

                List<int> volumeNodeBranchCountList = new List<int>();                     // list11 -> volumeNodeBranchCountList
                for (int i = 0; i < nodePoints.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    volumeNodeBranchCountList.Add(graphBranchCountList[volumeNodeIndexList[i]]);
                }

                // List<List<int>> volumeNodeGraph = new List<List<int>>();         // list12 -> volumeNodeGraph
                DataTree<int> volumeNodeGraph = new DataTree<int>();
                for (int i = 0; i < nodePoints.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    volumeNodeGraph.EnsurePath(i);
                    volumeNodeGraph.Branch(i).AddRange(graph.Branch(volumeNodeIndexList[i]));
                }

                List<int> boundaryNodeIndexList = new List<int>();                     // list13 -> boundaryNodeIndexList
                //for (int i = nodePoints.Count - NodeAttribute.BoundaryNodeCount; i < nodePoints.Count; i++)
                //{
                //    // 让boundaryNodeIndexList中的元素都是设定好的负值，以便后面的元素匹配
                //    boundaryNodeIndexList.Add(i - 2 * NodeAttribute.BoundaryNodeCount);
                //}

                for (int i = 0; i < NodeAttribute.BoundaryNodeCount; i++)
                {
                    // 让boundaryNodeIndexList中的元素都是设定好的负值，以便后面的元素匹配
                    boundaryNodeIndexList.Add(i - NodeAttribute.BoundaryNodeCount);
                }

                // List<int> list14 = new List<int>();                     // list14 ->
                // List<List<int>> list15 = new List<List<int>>();         // list15 -> volumeAdjacencyBoundary
                // volumeAdjacencyBoundary是每一支表示该支对应的点与哪些BoundaryNode相连
                DataTree<int> volumeAdjacencyBoundary = new DataTree<int>();

                for (int i = 0; i < volumeNodeGraph.BranchCount; i++)
                {
                    volumeAdjacencyBoundary.EnsurePath(i);
                    for (int j = 0; j < volumeNodeGraph.Branch(i).Count; j++)
                    {
                        if (boundaryNodeIndexList.Contains(volumeNodeGraph.Branch(i)[j]))
                        {
                            volumeAdjacencyBoundary.Branch(i).Add(volumeNodeGraph.Branch(i)[j]);
                        }
                    }
                }

                //List<List<double>> list16 = new List<List<double>>();   // list16 -> boundaryAdjacencyGraphPointX
                //List<List<double>> list17 = new List<List<double>>();   // list17 -> boundaryAdjacencyGraphPointY
                DataTree<double> boundaryAdjacencyGraphPointX = new DataTree<double>();
                DataTree<double> boundaryAdjacencyGraphPointY = new DataTree<double>();

                Point3d boundaryPoint;
                for (int i = 0; i < volumeNodeGraph.BranchCount; i++)
                {
                    boundaryAdjacencyGraphPointX.EnsurePath(i);
                    boundaryAdjacencyGraphPointY.EnsurePath(i);

                    for (int j = 0; j < volumeAdjacencyBoundary.Branch(i).Count; j++)
                    {
                        boundaryPoint = nodePoints[volumeAdjacencyBoundary.Branch(i)[j] + NodeAttribute.BoundaryNodeCount + NodeAttribute.VolumeNodeCount];
                        boundaryAdjacencyGraphPointX.Add(boundaryPoint.X, boundaryAdjacencyGraphPointX.Path(i));
                        boundaryAdjacencyGraphPointY.Add(boundaryPoint.Y, boundaryAdjacencyGraphPointY.Path(i));
                    }
                }


                List<double> boundaryAdjacencyGraphPointXSum = new List<double>();      // list20 -> boundaryAdjacencyGraphPointXSum
                List<double> boundaryAdjacencyGraphPointYSum = new List<double>();      // list21 -> boundaryAdjacencyGraphPointYSum
                Matrix MatrixXSum = new Matrix(volumeNodeGraph.BranchCount, 1);
                Matrix MatrixYSum = new Matrix(volumeNodeGraph.BranchCount, 1);

                for (int i = 0; i < volumeNodeGraph.BranchCount; i++)
                {
                    boundaryAdjacencyGraphPointXSum.Add(Sum(boundaryAdjacencyGraphPointX.Branch(i)));
                    MatrixXSum[i, 0] = Sum(boundaryAdjacencyGraphPointX.Branch(i));
                    boundaryAdjacencyGraphPointYSum.Add(Sum(boundaryAdjacencyGraphPointY.Branch(i)));
                    MatrixYSum[i, 0] = Sum(boundaryAdjacencyGraphPointY.Branch(i));
                }

                // volumeAdjacencyVolume是每一支表示该支对应的点与哪些volumeNode相连
                DataTree<int> volumeConnectivityVolume = new DataTree<int>();            // dataTree -> volumeConnectivityVolume

                for (int i = 0; i < nodePoints.Count - NodeAttribute.BoundaryNodeCount; i++)
                {
                    volumeConnectivityVolume.EnsurePath(i);

                    for (int j = 0; j < nodePoints.Count - NodeAttribute.BoundaryNodeCount; j++)
                    {
                        // 从子矩阵中找到1（表示有连接，0表示没有连接）
                        if (SubAdjacencyMatrix(GraphToAdjacencyMatrix(graph), volumeNodeIndexList).Branch(i)[j] == 1)
                        {
                            volumeConnectivityVolume.Branch(i).Add(j);
                        }
                    }
                }


                List<Point3d> list22 = new List<Point3d>();                             // list22 -> 
                // SubGraph输出
                DA.SetDataTree(2, volumeConnectivityVolume);


                // matrix3 -> harmonicMatrix
                Matrix harmonicMatrix = ToGHMatrix(HarmonicMatrix(SubAdjacencyMatrix(GraphToAdjacencyMatrix(graph), volumeNodeIndexList), volumeNodeBranchCountList));
                harmonicMatrix.Invert(0.0);
                Matrix inverseHarmonicMatrix = harmonicMatrix;                          // matrix4 -> inverseHarmonicMatrix
                Matrix matrix5 = inverseHarmonicMatrix * MatrixXSum;                    // matrix5 ->
                Matrix matrix6 = inverseHarmonicMatrix * MatrixYSum;                    // matrix6 ->
                List<Point3d> newBoundaryPoints = new List<Point3d>();                             // list23 -> newBoundaryPoints

                Point3d point;

                for (int i = 0; i < volumeNodeBranchCountList.Count; i++)
                {
                    // List<Point3d> list24 = newBoundaryPoints;                                      // list24 ->
                    point = new Point3d(matrix5[i, 0], matrix6[i, 0], 0.0);
                    newBoundaryPoints.Add(point);
                }

                for (int i = 0; i < volumeNodeIndexList.Count; i++)
                {
                    nodePoints[volumeNodeIndexList[i]] = newBoundaryPoints[i];
                }

                List<Point3d> newNodePointsRelocated = UtilityFunctions.Relocate(nodePoints, Plane.WorldXY, worldXY);
                DA.SetDataList("New Graph Vertices", newNodePointsRelocated);
                DA.SetDataList("ConvexFaceBorders", GMeshFaceBoundaries(newNodePointsRelocated, GraphEdgeList(graph, newNodePointsRelocated), worldXY));


                for (int i = 0; i < nodePoints.Count; i++)
                {
                    list22.Add(nodePoints[i]);
                }
                DA.SetDataList("SubVertices", list22);


                // subAttributes先不考虑


                //if (DA.GetDataList("NodeAttributes", list2))
                //{
                //    List<NodeAttribute> list26 = ;
                //    DA.SetDataList("SubAttributes", list26);
                //}
            }

        }

        /// <summary>
        /// double列表求和
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double Sum(List<double> x)
        {
            //int num = 0;
            //int num2 = x.Count - 1;
            //int num3 = num;
            double sum = 0;
            for (int i = 0; i < x.Count; i++)
            {
                sum += x[i];
            }
            return sum;
        }

        //public List<NodeAttribute> SubAttributes(List<NodeAttribute> attributes, List<int> subIndices)
        //{
        //    List<NodeAttribute> list = new List<NodeAttribute>();





        //    return list;
        //}

        /// <summary>
        /// 从一个矩阵中，取出与indices有关的那几行和几列数据，形成一个子矩阵
        /// </summary>
        /// <param name="adjacencyMatrixOfGraph"></param>
        /// <param name="indicesForSubMatrix"></param>
        /// <returns></returns>
        public DataTree<int> SubAdjacencyMatrix(DataTree<int> adjacencyMatrixOfGraph, List<int> indicesForSubMatrix)
        {
            // indicesRowData中存储MatrixOfGraph中，有indicesForSubMatrix中序号的，那几行0,1数据
            DataTree<int> indicesRowData = new DataTree<int>();                           // list -> indicesRowData
            // indicesRowAndColumeData中存储MatrixOfGraph中，有indicesForSubMatrix中序号的那几行和有indicesForSubMatrix中序号的那几列0,1数据
            DataTree<int> indicesRowAndColumeData = new DataTree<int>();                // list2 -> indicesRowAndColumeData

            // 取MatrixOfGraph中，有indicesForSubMatrix中序号的那几行
            for (int i = 0; i < indicesForSubMatrix.Count; i++)
            {
                indicesRowData.EnsurePath(i);
                indicesRowData.Branch(i).AddRange(adjacencyMatrixOfGraph.Branch(indicesForSubMatrix[i]));
            }

            // 取上面取到的那几行中，有indicesForSubMatrix中序号的那几列
            for (int i = 0; i < indicesForSubMatrix.Count; i++)
            {
                indicesRowAndColumeData.EnsurePath(i);
                for (int j = 0; j < indicesForSubMatrix.Count; j++)
                {
                    indicesRowAndColumeData.Branch(i).Add(indicesRowData.Branch(i)[indicesForSubMatrix[j]]);
                }
            }
            return indicesRowAndColumeData;
        }

        /// <summary>
        /// 将一个graph转化为有0,1的表示关系的邻接矩阵 Adjacency Matrix
        /// </summary>
        /// <param name="graph"></param>
        /// <returns></returns>
        public DataTree<int> GraphToAdjacencyMatrix(DataTree<int> graph)
        {
            // List<List<int>> list = new List<List<int>>();

            DataTree<int> dataTree = new DataTree<int>();

            //int num = graph.Count - 1;
            //int num2 = 0;
            //int num3 = num;
            //int num4 = num2;
            for (int i = 0; i < graph.BranchCount; i++)
            {
                dataTree.EnsurePath(i);
                //list.Add(new List<int>());

                //int num7 = 0;
                //int num8 = num;
                //int num9 = num7;
                for (int j = 0; j < graph.BranchCount; j++)
                {
                    // connectivity连接关系时
                    if (j < NodeAttribute.VolumeNodeCount)
                    {
                        if (graph.Branch(i).Contains(j))
                        {
                            //list[i].Add(1);
                            dataTree.Branch(i).Add(1);
                        }
                        else
                        {
                            //list[i].Add(0);
                            dataTree.Branch(i).Add(0);
                        }
                    }
                    // adjacency连接关系时
                    else
                    {
                        if (graph.Branch(i).Contains(j - NodeAttribute.VolumeNodeCount - NodeAttribute.BoundaryNodeCount))
                        {
                            //list[i].Add(1);
                            dataTree.Branch(i).Add(1);
                        }
                        else
                        {
                            //list[i].Add(0);
                            dataTree.Branch(i).Add(0);
                        }
                    }
                }
            }

            return dataTree;
        }
        
        /// <summary>
        /// 将DataTree<int>结构的矩阵转化为Rhino.Geometry.Matrix结构
        /// </summary>
        /// <param name="dataTree"></param>
        /// <returns></returns>
        public Matrix ToGHMatrix(DataTree<int> dataTree)
        {
            Matrix matrix = new Matrix(dataTree.Branch(0).Count, dataTree.BranchCount);

            if (dataTree == null)
            {
                matrix = null;
            }
            else
            {
                for (int i = 0; i < dataTree.BranchCount; i++)
                {
                    if (dataTree.Branch(i).Count != dataTree.Branch(0).Count)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Matrix must be rectagular! Your input graph is not valid");
                        matrix = null;
                    }

                    for (int j = 0; j < dataTree.Branch(0).Count; j++)
                    {
                        for (int k = 0; k < dataTree.BranchCount; k++)
                        {
                            matrix[j, k] = (double)dataTree.Branch(k)[j];
                        }
                    }
                }
            }
            return matrix;
        }

        /// <summary>
        /// 调和矩阵L Harmonic Matrix，又称拉普拉斯矩阵 Laplacian Matrix。是图的矩阵表示
        /// 在图论和计算机科学中，邻接矩阵A（英语：adjacency matrix）是一种方阵，用来表示有限图。它的每个元素代表各点之间是否有边相连。
        /// 在数学领域图论中，度数矩阵D是一个对角矩阵 ，其中包含的信息为的每一个顶点的度数，也就是说，每个顶点相邻的边数[1] 它可以和邻接矩阵一起使用以构造图的拉普拉斯算子矩阵。
        /// L = D - A
        /// </summary>
        /// <param name="adjacencyMatrix">邻接矩阵 Adjacency Matrix，对角线都是0，其他部分有数据0或1</param>
        /// <param name="degreesMatrix">度数矩阵 Degrees Matrix，因为是只有对角线有数据，其他部分都是零，所以可以用list存储</param>
        /// <returns></returns>
        public DataTree<int> HarmonicMatrix(DataTree<int> adjacencyMatrix, List<int> degreesMatrix)
        {
            // 调和矩阵L Harmonic Matrix，又称拉普拉斯矩阵 Laplacian Matrix
            DataTree<int> harmonicMatrix = new DataTree<int>();

            for (int i = 0; i < adjacencyMatrix.BranchCount; i++)
            {
                harmonicMatrix.EnsurePath(i);
                for (int j = 0; j < adjacencyMatrix.BranchCount; j++)
                {
                    // L = D - A
                    // 如果是对角线上的位置
                    if (j == i)
                    {
                        harmonicMatrix.Branch(i).Add(degreesMatrix[j] - 0);
                    }
                    // 其他位置
                    else
                    {
                        harmonicMatrix.Branch(i).Add(0 - adjacencyMatrix.Branch(i)[j]);
                    }
                }
            }

            return harmonicMatrix;
        }

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

        public List<Line> GraphEdgeList(DataTree<int> graph, List<Point3d> graphVertices)
        {
            List<Line> list = new List<Line>();

            for (int i = 0; i < graph.BranchCount; i++)
            {
                for (int j = 0; j < graph.Branch(i).Count; j++)
                {
                    if (graph.Branch(i)[j] >= 0)
                    {
                        Line item = new Line(graphVertices[i], graphVertices[graph.Branch(i)[j]]);
                        list.Add(item);
                    }
                    else
                    {
                        Line item = new Line(graphVertices[i], graphVertices[graph.Branch(i)[j] + NodeAttribute.BoundaryNodeCount + NodeAttribute.VolumeNodeCount]);
                        list.Add(item);
                    }
                    
                    
                    
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