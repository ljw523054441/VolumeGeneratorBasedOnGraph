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
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddPlaneParameter("BasePlane", "BasePlane", "A plane on which you want to get the Tutte convex drawing", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddIntegerParameter("Graph", "Graph", "A graph describing which node is connected to which other nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("NodePoints", "Vertices", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            pManager.AddGenericParameter("NodeAttributes", "Attributes", "图结构中节点所包含的属性", GH_ParamAccess.list);
            pManager[4].Optional = true;
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

            pManager.AddLineParameter("graphEdge", "", "", GH_ParamAccess.list);
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

                for (int i = 0; i < volumeNodeCount; i++)
                {
                    volumeNodeIndexList.Add(i);
                }

                //// List<int> list9 = new List<int>();                                      // list9 ->
                //// 每个节点（包括Volume和Boundary）的度的列表
                //List<int> graphBranchCountList = new List<int>();                       // list10 -> graphBranchCountList
                //for (int i = 0; i < graph.BranchCount; i++)
                //{
                //    graphBranchCountList.Add(graph.Branch(i).Count);
                //}

                // 每个volume点，其连接的数量，即该节点的度
                List<int> volumeNodeDegreeList = new List<int>();                     // list11 -> volumeNodeBranchCountList
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    volumeNodeDegreeList.Add(graph.Branch(i).Count);
                }

                // List<List<int>> volumeNodeGraph = new List<List<int>>();         // list12 -> volumeNodeGraph
                // 每个volume点与其他点的连接关系，包括其他volume点和boundary点
                DataTree<int> volumeNodeGraph = new DataTree<int>();
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    volumeNodeGraph.EnsurePath(i);
                    volumeNodeGraph.Branch(i).AddRange(graph.Branch(volumeNodeIndexList[i]));
                }

                //List<int> boundaryNodeIndexList = new List<int>();                     // list13 -> boundaryNodeIndexList

                //for (int i = 0; i < boundaryNodeCount; i++)
                //{
                //    // 让boundaryNodeIndexList中的元素都是 volumeNodeCount + i,以便后面的元素匹配。N:9(0)  E:10(1)  W:11(2)  S:12(3)  以此类推....
                //    boundaryNodeIndexList.Add(volumeNodeCount + i);
                //}

                // List<int> list14 = new List<int>();                     // list14 ->
                // List<List<int>> list15 = new List<List<int>>();         // list15 -> volumeAdjacencyBoundary



                //List<List<double>> list16 = new List<List<double>>();   // list16 -> volumeAdjacencyBoundaryPointX
                //List<List<double>> list17 = new List<List<double>>();   // list17 -> volumeAdjacencyBoundaryPointY

                //矩阵 P(outer)
                Matrix P_outer = new Matrix(boundaryNodeCount, 2);
                for (int i = 0; i < boundaryNodeCount; i++)
                {
                    P_outer[i, 0] = nodePoints[volumeNodeCount + i].X;
                    P_outer[i, 1] = nodePoints[volumeNodeCount + i].Y;
                }

                


                // 整个大的矩阵，包含四个部分  inner行-inner列  inner行-outer列：Q  outer行-inner列：Q(上标T)  outer行-outer列
                DataTree<int> wholeAdjacencyDataTree = GraphToAdjacencyMatrix_NxN(graph, globalParameter);
                Matrix wholeAdjacencyMatrix = ToGHMatrix(wholeAdjacencyDataTree);



                // 矩阵Q：inner_OuterAdjacencyMatrix
                DataTree<int> volumeBoundaryGraph = new DataTree<int>();
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    volumeBoundaryGraph.EnsurePath(i);
                    for (int j = 0; j < volumeNodeGraph.Branch(i).Count; j++)
                    {
                        // 选出outer的部分
                        if (volumeNodeGraph.Branch(i)[j] >= volumeNodeCount)
                        {
                            volumeBoundaryGraph.Branch(i).Add(volumeNodeGraph.Branch(i)[j]);
                        }
                    }
                }
                DataTree<int> inner_OuterAdjacencyDataTree = GraphToAdjacencyMatrix_NxM(volumeBoundaryGraph, globalParameter);
                // 这里直接得到的矩阵是行列反掉的，需要转置
                Matrix inner_OuterAdjacencyMatrix = ToGHMatrix(inner_OuterAdjacencyDataTree);
                inner_OuterAdjacencyMatrix.Transpose();



                DataTree<int> inner_InnerAdjacencyDataTree = SubAdjacencyMatrix(wholeAdjacencyDataTree, volumeNodeIndexList);
                Matrix inner_InnerAdjacencyMatrix = ToGHMatrix(inner_InnerAdjacencyDataTree);
                // 矩阵L(下标1)：inner_InnerLaplacianMatrix
                DataTree<int> inner_InnerLaplacianDataTree = LaplacianMatrix(inner_InnerAdjacencyDataTree, volumeNodeDegreeList);
                Matrix inner_InnerLaplacianMatrix = ToGHMatrix(inner_InnerLaplacianDataTree);       // matrix3 -> volumeAndVolumeLaplacianMatrix



                inner_InnerLaplacianMatrix.Invert(0.0);
                // 拉普拉斯矩阵的逆 L(下标1)(上标-1)：inverse_Inner_InnerLaplacianMatrix
                Matrix inverse_Inner_InnerLaplacianMatrix = inner_InnerLaplacianMatrix;                              // matrix4 -> inverseLaplacianMatrix



                // 求解矩阵 P(inner)
                Matrix P_inner = new Matrix(volumeNodeCount, 2);
                P_inner = inverse_Inner_InnerLaplacianMatrix * inner_OuterAdjacencyMatrix * P_outer;
                // 原插件这里没乘负号
                // P_inner.Scale(-1.0);


                // 从P_inner矩阵中生成新的inner点的坐标
                List<Point3d> newInnerPoints = new List<Point3d>();                             // list23 -> newBoundaryPoints
                Point3d innerPoint;
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    innerPoint = new Point3d(P_inner[i, 0], P_inner[i, 1], 0.0);
                    newInnerPoints.Add(innerPoint);
                }
                DA.SetDataList("SubVertices", newInnerPoints);

                // 将新的inner点与原来的outer点合并成一个新的列表
                List<Point3d> newNodePoints = new List<Point3d>();
                newNodePoints.AddRange(newInnerPoints);
                List<Point3d> sortedBoundaryNodePoints = new List<Point3d>();
                for (int i = 0; i < boundaryNodeCount; i++)
                {
                    sortedBoundaryNodePoints.Add(nodePoints[volumeNodeCount + i]);
                }
                newNodePoints.AddRange(sortedBoundaryNodePoints);
                
                // 将新的inner+outer列表重新定位在基于worldXY变量的另一个地方
                List<Point3d> newNodePointsRelocated = UtilityFunctions.Relocate(newNodePoints, Plane.WorldXY, worldXY);
                DA.SetDataList("New Graph Vertices", newNodePointsRelocated);
                
                // Graph的所有Edge，在新的位置上画Edge
                List<Line> graphEdge = GraphEdgeList(graph, newNodePointsRelocated, globalParameter);
                DA.SetDataList("graphEdge", graphEdge);

                // Graph形成的每个convex，在新的位置上画Convex
                List<Point3d> sortedBoundaryNodePointsRelocated = new List<Point3d>();
                for (int i = volumeNodeCount; i < newNodePointsRelocated.Count; i++)
                {
                    sortedBoundaryNodePointsRelocated.Add(newNodePointsRelocated[i]);
                }
                List<Polyline> allSplitFaceBoundarys = GMeshFaceBoundaries(sortedBoundaryNodePointsRelocated, graphEdge, worldXY);
                DA.SetDataList("ConvexFaceBorders", allSplitFaceBoundarys);

                // subGraph的每一支表示一个volumeNode与哪些其他的volumeNode相连
                DataTree<int> subGraph = new DataTree<int>();            // dataTree -> volumeConnectivityVolume
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    subGraph.EnsurePath(i);

                    for (int j = 0; j < volumeNodeCount; j++)
                    {
                        if (inner_InnerAdjacencyDataTree.Branch(i)[j] == 1)
                        {
                            subGraph.Branch(i).Add(j);
                        }
                    }
                }
                // SubGraph输出
                DA.SetDataTree(2, subGraph);


                // SubAttribute输出
                DA.GetDataList("NodeAttributes", nodeAttributes);
                List<NodeAttribute> innerNodeAttributes = new List<NodeAttribute>();
                for (int i = 0; i < volumeNodeCount; i++)
                {
                    innerNodeAttributes.Add(nodeAttributes[i]);
                }
                DA.SetDataList("SubAttributes", innerNodeAttributes);





                //Matrix matrix5 = inverse_Inner_InnerLaplacianMatrix * MatrixXSum;                    // matrix5 ->
                //Matrix matrix6 = inverse_Inner_InnerLaplacianMatrix * MatrixYSum;                    // matrix6 ->


                //Point3d point;

                //for (int i = 0; i < volumeNodeDegreeList.Count; i++)
                //{
                //    // List<Point3d> list24 = newBoundaryPoints;                                      // list24 ->
                //    point = new Point3d(matrix5[i, 0], matrix6[i, 0], 0.0);
                //    newInnerPoints.Add(point);
                //}

                //for (int i = 0; i < volumeNodeIndexList.Count; i++)
                //{
                //    nodePoints[volumeNodeIndexList[i]] = newInnerPoints[i];
                //}

                //List<Point3d> newNodePointsRelocated = UtilityFunctions.Relocate(nodePoints, Plane.WorldXY, worldXY);
                //DA.SetDataList("New Graph Vertices", newNodePointsRelocated);





                // DA.SetDataList("ConvexFaceBorders", GMeshFaceBoundaries(newNodePointsRelocated, GraphEdgeList(graph, newNodePointsRelocated, globalParameter), worldXY));

                //List<Point3d> subVertices = new List<Point3d>();                             // list22 -> 
                //for (int i = 0; i < nodePoints.Count; i++)
                //{
                //    subVertices.Add(nodePoints[i]);
                //}



                //List<Point3d> subVerticesRelocated = UtilityFunctions.Relocate(nodePoints, Plane.WorldXY, worldXY);
                //DA.SetDataList("SubVertices", newNodePointsRelocated);


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

        public List<NodeAttribute> SubAttributes(List<NodeAttribute> attributes, List<int> subIndices)
        {
            List<NodeAttribute> list = new List<NodeAttribute>();





            return list;
        }

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
        /// <param name="graphForNxN"></param>
        /// <returns></returns>
        public DataTree<int> GraphToAdjacencyMatrix_NxN(DataTree<int> graphForNxN, GlobalParameter globalParameter)
        {
            DataTree<int> dataTreeForNxN = new DataTree<int>();

            for (int i = 0; i < graphForNxN.BranchCount; i++)
            {
                dataTreeForNxN.EnsurePath(i);
                for (int j = 0; j < graphForNxN.BranchCount; j++)
                {
                    if (graphForNxN.Branch(i).Contains(j))
                    {
                        //list[i].Add(1);
                        dataTreeForNxN.Branch(i).Add(1);
                    }
                    else
                    {
                        //list[i].Add(0);
                        dataTreeForNxN.Branch(i).Add(0);
                    }

                    //// connectivity连接关系时
                    //if (j < globalParameter.VolumeNodeCount)
                    //{
                    //    if (graphForNxN.Branch(i).Contains(j))
                    //    {
                    //        //list[i].Add(1);
                    //        dataTreeForNxN.Branch(i).Add(1);
                    //    }
                    //    else
                    //    {
                    //        //list[i].Add(0);
                    //        dataTreeForNxN.Branch(i).Add(0);
                    //    }
                    //}
                    //// adjacency连接关系时
                    //else
                    //{
                    //    if (graphForNxN.Branch(i).Contains(j))
                    //    {
                    //        //list[i].Add(1);
                    //        dataTreeForNxN.Branch(i).Add(1);
                    //    }
                    //    else
                    //    {
                    //        //list[i].Add(0);
                    //        dataTreeForNxN.Branch(i).Add(0);
                    //    }
                    //}
                }
            }
            return dataTreeForNxN;
        }

        public DataTree<int> GraphToAdjacencyMatrix_NxM(DataTree<int> graphForNxM, GlobalParameter globalParameter)
        {
            DataTree<int> dataTreeForNxM = new DataTree<int>();

            for (int i = 0; i < graphForNxM.BranchCount; i++)
            {
                dataTreeForNxM.EnsurePath(i);
                for (int j = 0; j < globalParameter.BoundaryNodeCount; j++)
                {
                    if (graphForNxM.Branch(i).Contains(j + globalParameter.VolumeNodeCount))
                    {
                        dataTreeForNxM.Branch(i).Add(1);
                    }
                    else
                    {
                        dataTreeForNxM.Branch(i).Add(0);
                    }
                }
            }
            return dataTreeForNxM;
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
        public DataTree<int> LaplacianMatrix(DataTree<int> adjacencyMatrix, List<int> degreesMatrix)
        {
            // 调和矩阵L Harmonic Matrix，又称拉普拉斯矩阵 Laplacian Matrix
            DataTree<int> laplacianMatrix = new DataTree<int>();

            for (int i = 0; i < adjacencyMatrix.BranchCount; i++)
            {
                laplacianMatrix.EnsurePath(i);
                for (int j = 0; j < adjacencyMatrix.BranchCount; j++)
                {
                    // L = D - A
                    // 如果是对角线上的位置
                    if (j == i)
                    {
                        laplacianMatrix.Branch(i).Add(degreesMatrix[j] - 0);
                    }
                    // 其他位置
                    else
                    {
                        laplacianMatrix.Branch(i).Add(0 - adjacencyMatrix.Branch(i)[j]);
                    }
                }
            }

            return laplacianMatrix;
        }

        /// <summary>
        /// 将排好序的点连接成多段线（NurbsCurve类型）
        /// </summary>
        /// <param name="sortedBoundaryPoints"></param>
        /// <param name="basePlane"></param>
        /// <returns></returns>
        public Polyline SortedPointsToPolyline(List<Point3d> sortedBoundaryPoints, Plane basePlane)
        {
            // 把startPoint放到最后作为endPoint，这样polyline才会闭合
            sortedBoundaryPoints.Add(sortedBoundaryPoints[0]);

            Polyline convex = new Polyline(sortedBoundaryPoints);
            return convex;
        }


        //public Curve ConvexBorder(List<Point3d> vertices, Plane basePlane)
        //{
        //    Point3dList point3dList = new Point3dList(vertices);
        //    point3dList.Sort();
        //    Plane plane = basePlane;
        //    // 旋转 PI / 4
        //    plane.Rotate(0.7853981633974483, plane.ZAxis);
        //    Rectangle3d rectangle3d = new Rectangle3d(plane, point3dList.First, point3dList.Last);
        //    return rectangle3d.ToNurbsCurve();
        //}

        public List<Polyline> GMeshFaceBoundaries(List<Point3d> sortedBoundaryPoints, List<Line> edges, Plane basePlane)
        {
            List<Curve> splitLines = new List<Curve>();
            foreach (Line line in edges)
            {
                splitLines.Add(line.ToNurbsCurve());
            }

            Polyline polyline = SortedPointsToPolyline(sortedBoundaryPoints, basePlane);
            Curve curve = polyline.ToNurbsCurve();

            Brep[] PlanarBrepArray = Brep.CreatePlanarBreps(curve, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            Brep planarBreps = PlanarBrepArray[0];
            BrepFace planarBrepFace = planarBreps.Faces[0];
            Brep splittedPlanarBrepFaces = planarBrepFace.Split(splitLines, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            
            List<Polyline> allConvexPolyline = new List<Polyline>();
            foreach (BrepFace brepFace in splittedPlanarBrepFaces.Faces)
            {
                Brep singleSplitFaceBrep = brepFace.DuplicateFace(true);
                Point3d[] singleSplitFaceVertices = singleSplitFaceBrep.DuplicateVertices();
                Curve splitFaceConvexCure = Curve.JoinCurves(singleSplitFaceBrep.Faces[0].DuplicateFace(true).DuplicateEdgeCurves())[0];
                Polyline splitFaceConvexPolyline = null;
                splitFaceConvexCure.TryGetPolyline(out splitFaceConvexPolyline);
                allConvexPolyline.Add(splitFaceConvexPolyline);
            }

            return allConvexPolyline;
        }

        public List<Line> GraphEdgeList(DataTree<int> graph, List<Point3d> graphVertices, GlobalParameter globalParameter)
        {
            List<Line> list = new List<Line>();

            for (int i = 0; i < graph.BranchCount; i++)
            {
                for (int j = 0; j < graph.Branch(i).Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graph.Branch(i)[j]]);
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