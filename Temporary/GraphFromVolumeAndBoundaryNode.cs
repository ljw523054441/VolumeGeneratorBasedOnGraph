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
// using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class GraphFromVolumeAndBoundaryNode : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CSVtoGraph class.
        /// </summary>
        public GraphFromVolumeAndBoundaryNode()
          : base("GraphFromVolumeAndBoundaryNode", "GraphFromVolumeAndBoundaryNode",
              "基于体量节点和边界节点的图结构",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
        }

        private List<string> Texts = new List<string>();
        private List<Point3d> Locations = new List<Point3d>();
        private double Fontsize;
        private List<Polyline> NEWS_Polygons = new List<Polyline>();
        private List<Brep> NEWS_Surfaces = new List<Brep>();
        private List<string> NEWS_Tags = new List<string>();
        private List<Point3d> NEWS_Location = new List<Point3d>();
        private double NEWS_Size;
        private Plane BasePlane;


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);
            
            pManager.AddPointParameter("VolumeNodePoints", "VolumeNodePoints", "能够代表volumeNode节点的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePointsExceptNEWS", "BoundaryNodePointsExceptNEWS", "能够代表除了NEWS点外的其他BoundaryNode节点的抽象点", GH_ParamAccess.list);
            pManager.AddIntegerParameter("VolumeConnectivityTree", "ConnectivityTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "AdjacencyTree", "体量邻近关系树", GH_ParamAccess.tree);
            pManager.AddGenericParameter("VolumeLabelList", "LabelList", "体量标签列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("VolumeAttributeList", "AttributeList", "体量属性列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("BoundaryLabelListExceptNEWS", "BoundaryLabelList", "除了NEWS外的其他BoundaryNode的标签列表", GH_ParamAccess.list);

            // pManager.AddGenericParameter("BoundaryAttributeListExceptNEWS", "BoundaryAttributeList", "除了NEWS外的其他BoundaryNode的属性列表", GH_ParamAccess.list);
            pManager[2].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Graph", "Graph", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);
            // pManager.AddIntegerParameter("BoundaryAdjacencyGraph")
            //pManager.AddPointParameter("VolumeNode", "GraphNode", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            //pManager.AddPointParameter("BoundaryNode", "BoundaryNode", "用抽象点（point）表示的边界（node）", GH_ParamAccess.list);
            //pManager.AddGenericParameter("VolumeNodeAttributes", "VolumeAttributes", "图结构中节点所包含的属性", GH_ParamAccess.list);
            //pManager.AddGenericParameter("BoundaryNodeAttributes", "BoundaryAttributes", "边界节点所包含的属性", GH_ParamAccess.list);
            pManager.AddPointParameter("NodePoints", "GraphNode", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            pManager.AddGenericParameter("NodeAttributes", "NodeAttributes", "图结构中节点所包含的属性", GH_ParamAccess.list);
            // pManager.AddIntegerParameter("DebugIndex", "index", "", GH_ParamAccess.list);


        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Texts.Clear();
            Locations.Clear();
            NEWS_Polygons.Clear();
            NEWS_Surfaces.Clear();
            NEWS_Tags.Clear();
            NEWS_Location.Clear();

            // 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int volumeNodeCount = globalParameter.VolumeNodeCount;
            int boundaryNodeCount = globalParameter.BoundaryNodeCount;


            // Local Varible 初始化
            List<Point3d> volumeNodeList = new List<Point3d>();                                   // list -> volumeNodeList
            List<Point3d> boundaryNodeListExceptNEWS = new List<Point3d>();

            //List<Line> allEdgesList = new List<Line>();                                     // list2 -> allEdgesList  allLineList uniqueLineList
            //List<Line> allAdjacencyList = new List<Line>();
            GH_Structure<GH_Integer> gh_Structure_ConnectivityTree = null;
            GH_Structure<GH_Integer> gh_Structure_AdjacencyTree = null;
            DataTree<int> volumeConnectivityTree = new DataTree<int>();
            DataTree<int> boundaryAdjacencyTree = new DataTree<int>();


            List<string> volumeLabelList = new List<string>();                              // list3 -> volumeLabelList
            // List<double> volumeAreaAttributeList = new List<double>();                      // list4 -> volumeAreaAttributeList
            List<NodeAttribute> volumeNodeAttributeList = new List<NodeAttribute>();

            List<string> boundaryLabelList = new List<string>();
            boundaryLabelList.AddRange(new string[]
                {
                    "N",
                    "E",
                    "W",
                    "S"
                });
            // List<double> boundaryAreaAttributeList = new List<double>();
            List<NodeAttribute> boundaryNodeAttributeList = new List<NodeAttribute>();

            List<string> boundaryLabelListExceptNEWS = new List<string>();
            // List<double> boundaryAreaAttributeListExceptNEWS = new List<double>();
            List<NodeAttribute> boundaryNodeAttributeListExceptNEWS = new List<NodeAttribute>();

            // double totalArea = 0.0;                                                               // num -> totalArea
            // Plane worldXY = Plane.WorldXY;

            List<Polyline> NEWS_Regions = new List<Polyline>();                             // list5 -> NEWS_Regions
            List<Point3d> NEWS_Vertices = new List<Point3d>();                              // collection -> NEWS_Vertices

            DataTree<int> graph = new DataTree<int>();
            List<Point3d> nodePoints = new List<Point3d>();
            List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();



            // 获取NEWS四个方向的Region，同时将四个Region的中心添加到节点列表nodeList中
            if (DA.GetDataList<Point3d>("VolumeNodePoints", volumeNodeList))
            {
                DA.GetDataList<Point3d>("BoundaryNodePointsExceptNEWS", boundaryNodeListExceptNEWS);

                // 获取NEWS四个方向的Region，同时将四个Region的中心添加到节点列表nodeList中
                SketchBoxNEWS(PadSquare(volumeNodeList, Plane.WorldXY), ref NEWS_Regions, ref NEWS_Vertices);

                
                if (DA.GetDataList<string>("VolumeLabelList", volumeLabelList) && volumeLabelList.Count != volumeNodeList.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "number of nodes and labels must be the same!");
                    return;
                }
                if (DA.GetDataList<NodeAttribute>("VolumeAttributeList", volumeNodeAttributeList) && volumeNodeAttributeList.Count != volumeNodeList.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "number of nodes and areas must be the same!");
                    return;
                }
                if (DA.GetDataList<string>("BoundaryLabelListExceptNEWS", boundaryLabelListExceptNEWS) && boundaryLabelListExceptNEWS.Count != boundaryNodeListExceptNEWS.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "number of nodes and labels must be the same!");
                    return;
                }


                // 输出Graph的树形结构
                if (DA.GetDataTree<GH_Integer>("VolumeConnectivityTree", out gh_Structure_ConnectivityTree)
                & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyTree", out gh_Structure_AdjacencyTree))
                {
                    // 自定义的boundary节点按照角度进行排序，插入到NEWS中，得到所有boundary点之间的邻接关系
                    Sphere NEWSSphere = Sphere.FitSphereToPoints(NEWS_Vertices);
                    Point3d centerPoint = NEWSSphere.Center;

                    List<Point3d> boundaryNodeList = new List<Point3d>();
                    boundaryNodeList.AddRange(NEWS_Vertices);
                    boundaryNodeList.AddRange(boundaryNodeListExceptNEWS);

                    List<int> sortedBoundaryNodeIndexList = new List<int>();

                    List<Point3d> SortedBoundaryNodeList = UtilityFunctions.SortPolyPoints(boundaryNodeList, centerPoint);
                    // SortedBoundaryNodeList.Reverse();

                    for (int i = 0; i < SortedBoundaryNodeList.Count; i++)
                    {
                        int boundaryPointsIndex = boundaryNodeList.IndexOf(SortedBoundaryNodeList[i]);
                        sortedBoundaryNodeIndexList.Add(boundaryPointsIndex);
                    }

                    // DA.SetDataList("DebugIndex", sortedBoundaryNodeIndexList);


                    // 构造包含所有node节点的列表，顺序是 inner node + outer node
                    nodePoints.AddRange(volumeNodeList);
                    //nodePoints.AddRange(NEWS_Vertices);
                    //nodePoints.AddRange(boundaryNodeListExceptNEWS);
                    nodePoints.AddRange(SortedBoundaryNodeList);

                    // 输出NodePoints
                    DA.SetDataList("NodePoints", nodePoints);


                    List<Point3d> volumeAndNEWSNodes = volumeNodeList;
                    volumeAndNEWSNodes.AddRange(NEWS_Vertices);
                    // 包含volume节点和NEWS节点的大球
                    Sphere sphere = Sphere.FitSphereToPoints(volumeAndNEWSNodes);

                    // 计算每个volume点的面积占比
                    UtilityFunctions.CalculateAreaProportion(volumeNodeAttributeList);
                    


                    // NEWS节点的属性列表
                    for (int i = 0; i < NEWS_Vertices.Count; i++)
                    {
                        NodeAttribute NEWSNodeAttribute = new NodeAttribute(boundaryLabelList[i], Math.PI * Math.Pow(sphere.Radius, 2) / (double)(4 * volumeAndNEWSNodes.Count));
                        boundaryNodeAttributeList.Add(NEWSNodeAttribute);
                    }

                    // 除NEWS节点以外的其他BoundaryNode的属性列表
                    for (int i = 0; i < boundaryNodeListExceptNEWS.Count; i++)
                    {
                        NodeAttribute boundaryNodeExceptNEWSAttribute = new NodeAttribute(boundaryLabelListExceptNEWS[i], Math.PI * Math.Pow(sphere.Radius, 2) / (double)(4 * volumeAndNEWSNodes.Count));
                        boundaryNodeAttributeList.Add(boundaryNodeExceptNEWSAttribute);
                    }

                    // 用前面生成的对于Boundary点的排序列表来对BoundaryAttribute列表进行对应的排序
                    List<NodeAttribute> SortedBoundaryNodeAttributes = new List<NodeAttribute>();
                    foreach (int index in sortedBoundaryNodeIndexList)
                    {
                        SortedBoundaryNodeAttributes.Add(boundaryNodeAttributeList[index]);
                    }

                    // 构造包含所有node节点的属性列表，顺序是 inner node + outer node
                    nodeAttributes.AddRange(volumeNodeAttributeList);
                    nodeAttributes.AddRange(SortedBoundaryNodeAttributes);
                    // 输出NodeAttributes
                    DA.SetDataList("NodeAttributes", nodeAttributes);



                    // 将ConnectivityTree从GH_Structure<GH_Integer>转化为DataTree<int>
                    UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_ConnectivityTree, ref volumeConnectivityTree);
                    // 将AdjacencyTree从GH_Structure<GH_Integer>转化为DataTree<int>
                    UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_AdjacencyTree, ref boundaryAdjacencyTree);

                    // 将每个VolumeNode对于其他 inner + outer 点的邻接关系，添加到每个VolumeNode分支上
                    for (int i = 0; i < volumeConnectivityTree.BranchCount; i++)
                    {
                        graph.EnsurePath(i);
                        graph.AddRange(volumeConnectivityTree.Branches[i], graph.Path(i));
                        for (int j = 0; j < boundaryAdjacencyTree.Branches[i].Count; j++)
                        {
                            // 注意这里add的item应该是排过序后的NEWS对应的index
                            graph.Add(sortedBoundaryNodeIndexList.IndexOf(boundaryAdjacencyTree.Branch(i)[j]) + volumeNodeCount, graph.Path(i));
                        }
                    }


                    // 对每个boundary点（outer），向树形数据中添加新的path（BoundaryNode分支）
                    for (int i = 0; i < boundaryNodeCount; i++)
                    {
                        graph.EnsurePath(volumeNodeCount + i);
                    }

                    // 每个BoundaryNode对于其他 inner 点的邻接关系，添加到每个BoundaryNode分支上，i相当于是volumeNode的序号了
                    for (int i = 0; i < boundaryAdjacencyTree.BranchCount; i++)
                    {
                        for (int j = 0; j < boundaryAdjacencyTree.Branch(i).Count; j++)
                        {
                            graph.Add(i, graph.Path(volumeNodeCount + (sortedBoundaryNodeIndexList.IndexOf(boundaryAdjacencyTree.Branch(i)[j]))));
                        }
                    }

                    // 每个BoundaryNode对于其他 outer 点的邻接关系（outer点顺时针排序，每个outer点与其左右的两个outer点分别连接），添加到每个BoundaryNode分支上
                    for (int i = 0; i < sortedBoundaryNodeIndexList.Count; i++)
                    {
                        if (i == 0)
                        {
                            graph.Add((sortedBoundaryNodeIndexList.Count - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                            graph.Add((i + 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                            continue;
                        }
                        if (i == sortedBoundaryNodeIndexList.Count - 1)
                        {
                            graph.Add((i - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                            graph.Add((0) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                            continue;
                        }

                        graph.Add((i - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                        graph.Add((i + 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));

                    }

                    DA.SetDataTree(0, graph);
                }
            }
        }

        /// <summary>
        /// 将点的列表按照从正北开始，逆时针排序
        /// </summary>
        /// <param name="vPoints"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        //public List<Point3d> SortPolyPoints(List<Point3d> vPoints, Point3d center)
        //{
        //    List<Point3d> cloneList = new List<Point3d>();
        //    cloneList.AddRange(vPoints);


        //    if (cloneList == null || cloneList.Count == 0)
        //    {
        //        return null;
        //    }

        //    for (int i = 0; i < cloneList.Count - 1; i++)
        //    {
        //        for (int j = 0; j < cloneList.Count - i - 1; j++)
        //        {
        //            bool flag = PointCmp(cloneList[j], cloneList[j + 1], center);
        //            if (flag)
        //            {
        //                Point3d tmp = cloneList[j];
        //                cloneList[j] = cloneList[j + 1];
        //                cloneList[j + 1] = tmp;
        //            }
        //        }
        //    }
        //    return cloneList;
        //}

        ///// <summary>
        ///// 若点a大于点b，即点a在点b的顺时针方向，返回true，否则返回false
        ///// </summary>
        ///// <param name="a"></param>
        ///// <param name="b"></param>
        ///// <param name="center"></param>
        ///// <returns></returns>
        //private bool PointCmp(Point3d a, Point3d b, Point3d center)
        //{
        //    Vector3d vectorOA = new Vector3d(a) - new Vector3d(center);
        //    Vector3d vectorOB = new Vector3d(b) - new Vector3d(center);

        //    // OA,OB分别与0,1,0的夹角的弧度值
        //    double angleOA = Vector3d.VectorAngle(new Vector3d(0, 1, 0), vectorOA);
        //    double angleOB = Vector3d.VectorAngle(new Vector3d(0, 1, 0), vectorOB);
            
        //    // 向量0,1,0和向量OA的叉积
        //    Vector3d vectorZOA = Vector3d.CrossProduct(new Vector3d(0, 1, 0), vectorOA);
        //    if (vectorZOA.Z < 0)
        //    {
        //        angleOA = 2 * Math.PI - angleOA;
        //    }
        //    Vector3d vectorZOB = Vector3d.CrossProduct(new Vector3d(0, 1, 0), vectorOB);
        //    if (vectorZOB.Z < 0)
        //    {
        //        angleOB = 2 * Math.PI - angleOB;
        //    }


        //    if (angleOA < angleOB)
        //    {
        //        return false;
        //    }
        //    return true;
            
        //}

        /// <summary>
        /// 输入skechpad正方形，得到对应的NEWS四个区域及NEWS四个中心点
        /// </summary>
        /// <param name="rectangleRegion"></param>
        /// <param name="NEWS_Regions"></param>
        /// <param name="NEWS_Nodes"></param>
        public void SketchBoxNEWS(Rectangle3d rectangleRegion, ref List<Polyline> NEWS_Regions, ref List<Point3d> NEWS_Nodes)
        {
            Point3d centerPoint = rectangleRegion.Center;
            double area = rectangleRegion.Area;
            List<Point3d> cornerPointList = new List<Point3d>();
            List<Point3d> transformedCornerPointList = new List<Point3d>();
            double num = Math.Sqrt(area) / 4;

            for (int i = 0; i < 4; i++)
            {
                /* Index of corner, valid values are:
                    0 = lower left(min - x, min - y)
                    1 = lower right(max - x, min - y)
                    2 = upper right(max - x, max - y)
                    3 = upper left(min - x, max - y) */
                Point3d cornerPoint = rectangleRegion.Corner(i);
                cornerPointList.Add(cornerPoint);
                Vector3d centerToCorner = new Vector3d(cornerPoint.X - centerPoint.X, cornerPoint.Y - centerPoint.Y, cornerPoint.Z - centerPoint.Z);
                centerToCorner.Unitize();
                // 对角线长度的四分之一
                centerToCorner = centerToCorner * num * Math.Sqrt(2.0);
                Transform matrixforMotion = Transform.Translation(centerToCorner);
                cornerPoint.Transform(matrixforMotion);
                transformedCornerPointList.Add(cornerPoint);
            }

            List<Point3d> NEWSPolylineCenterList = new List<Point3d>();
            List<Polyline> NEWSPolylineList = new List<Polyline>();
            // Npolyline and Npolyline Center
            Polyline Npolyline = new Polyline(new Point3d[]
            {
                cornerPointList[2],
                cornerPointList[3],
                transformedCornerPointList[3],
                transformedCornerPointList[2],
                cornerPointList[2]

            });
            NEWSPolylineCenterList.Add(Npolyline.CenterPoint());
            // Epolyline and Epolyline Center
            Polyline Epolyline = new Polyline(new Point3d[]
            {
                transformedCornerPointList[2],
                cornerPointList[2],
                cornerPointList[1],
                transformedCornerPointList[1],
                transformedCornerPointList[2]
            });
            NEWSPolylineCenterList.Add(Epolyline.CenterPoint());
            // Wpolyline and Wpolyline Center
            Polyline Wpolyline = new Polyline(new Point3d[]
            {
                cornerPointList[3],
                transformedCornerPointList[3],
                transformedCornerPointList[0],
                cornerPointList[0],
                cornerPointList[3]
            });
            NEWSPolylineCenterList.Add(Wpolyline.CenterPoint());
            // Spolyline and Spolyline Center
            Polyline Spolyline = new Polyline(new Point3d[]
            {
                cornerPointList[1],
                cornerPointList[0],
                transformedCornerPointList[0],
                transformedCornerPointList[1],
                cornerPointList[1]
            });
            NEWSPolylineCenterList.Add(Spolyline.CenterPoint());

            NEWSPolylineList.AddRange(new Polyline[]
            {
                Npolyline,
                Epolyline,
                Wpolyline,
                Spolyline
            });
            NEWS_Regions = NEWSPolylineList;
            NEWS_Nodes = NEWSPolylineCenterList;
        }

        /// <summary>
        /// 获得能够包围输入节点的长方形Sketchpad
        /// </summary>
        /// <param name="points"></param>
        /// <param name="basePlane"></param>
        /// <returns></returns>
        public Rectangle3d PadSquare(List<Point3d> points, Plane basePlane)
        {
            // 如果没有节点输入，直接返回默认的长方形
            if (points == null | points.Count <= 1)
            {
                return default(Rectangle3d);
            }

            // 取输入节点的FitSphere，可以理解成包围球
            Sphere fitSphere = Sphere.FitSphereToPoints(points);
            Point3d fitSphereCenter = fitSphere.Center;
            // 将平面原点移动到球心
            basePlane.Origin = fitSphere.Center;
            // 对包围球取BoundingBox包围盒
            Box fitSphereBoundingBox = new Box(fitSphere.BoundingBox);
            // Interval intervalX = new Interval(fitSphereBoundingBox.X.T0 * 2.0, fitSphereBoundingBox.X.T1 * 2.0);
            // Interval intervalY = new Interval(fitSphereBoundingBox.Y.T0 * 2.0, fitSphereBoundingBox.Y.T1 * 2.0);
            // 用包围盒的x范围，y范围来做长方形
            Rectangle3d result = new Rectangle3d(Plane.WorldXY, fitSphereBoundingBox.X, fitSphereBoundingBox.Y);
            // 将包围盒得到的长方形放大两倍
            Transform scaleXYTransform = Transform.Scale(basePlane, 2.0, 2.0, 1.0);
            result.Transform(scaleXYTransform);
            return result;
        }

        /// <summary>
        /// 检测输入的连接（line）中，哪些是表示节点与NEWS的邻接关系的，将表示该邻接关系的线变为节点到NEWS中心的线，并且添加NEWS四个节点之间的连接
        /// </summary>
        /// <param name="edges"></param>
        /// <param name="rectangles"></param>
        public void NodeLinkToNEWS(ref List<Line> edges, List<Polyline> rectangles)
        {
            if (rectangles.Count != 4)
            {
                return;
            }
            if (edges == null)
            {
                return;
            }

            for (int i = 0; i < edges.Count; i++)
            {
                // bool flag = false;
                for (int j = 0; j < rectangles.Count; j++)
                {
                    
                    if (rectangles[j].ToNurbsCurve().Contains(edges[i].To, Plane.WorldXY, RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == PointContainment.Inside)
                    {
                        edges[i] = new Line(edges[i].From, rectangles[j].CenterPoint());
                        // flag = true;
                    }
                    else if (rectangles[j].ToNurbsCurve().Contains(edges[i].From, Plane.WorldXY,RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == PointContainment.Inside)
                    {
                        edges[i] = new Line(edges[i].To, rectangles[j].CenterPoint());
                        // flag = true;
                    }
                    else
                    {
                        continue;
                    }
                }
                // if (!flag)
                // {
                //     edges[i] = edges[i];
                // }
            }

            // 添加NEWS四个节点之间的连接
            edges.Add(new Line(rectangles[0].CenterPoint(), rectangles[1].CenterPoint()));
            edges.Add(new Line(rectangles[1].CenterPoint(), rectangles[3].CenterPoint()));
            edges.Add(new Line(rectangles[3].CenterPoint(), rectangles[2].CenterPoint()));
            edges.Add(new Line(rectangles[2].CenterPoint(), rectangles[0].CenterPoint()));
        }

        /// <summary>
        /// 将体量连接关系树转化为edge列表
        /// </summary>
        /// <param name="nodeList"></param>
        /// <param name="nodeConnectivityTree"></param>
        /// <returns></returns>
        public List<Line> ConnectivityTreeToEdges(List<Point3d> nodeList, DataTree<int> nodeConnectivityTree)
        {
            List<Line> edges = new List<Line>();

            for (int i = 0; i < nodeConnectivityTree.BranchCount; i++)
            {
                for (int j = 0; j < nodeConnectivityTree.Branch(i).Count; j++)
                {
                    Line edge = new Line(nodeList[i], nodeList[nodeConnectivityTree.Branch(i)[j]]);
                    edges.Add(edge);
                }
            }

            return edges;
        }

        /// <summary>
        /// 将体量与边界的邻接关系树转化为edge列表
        /// </summary>
        /// <param name="nodeList"></param>
        /// <param name="NEWS_Nodes"></param>
        /// <param name="boundaryAdjacencyTree"></param>
        /// <returns></returns>
        public List<Line> BoundaryAdjacencyTreeToEdges(List<Point3d> nodeList, List<Point3d> NEWS_Nodes, DataTree<int> boundaryAdjacencyTree)
        {
            List<Line> edges = new List<Line>();
            for (int i = 0; i < boundaryAdjacencyTree.BranchCount; i++)
            {
                for (int j = 0; j < boundaryAdjacencyTree.Branch(i).Count; j++)
                {
                    Line edge = new Line(nodeList[i], NEWS_Nodes[boundaryAdjacencyTree.Branch(i)[j]]);
                }
            }

            return edges;
        }

        /// <summary>
        /// 剔除重复的线
        /// </summary>
        /// <param name="lines">输入要检测有无重复线的线列表。</param>
        public List<Line> CleanLinesUD(List<Line> lines)
        {
            List<Line> uniqueLineList = new List<Line>();
            double tolerance = RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
            foreach (Line line in lines)
            {
                Line LF = line;
                LF.Flip();
                Line TL = line;
                
                if (!uniqueLineList.Exists((Line x) => x.EpsilonEquals(TL,tolerance)) & !uniqueLineList.Exists((Line x) => x.EpsilonEquals(LF,tolerance)))
                {
                    uniqueLineList.Add(line);
                }
            }
            return uniqueLineList;
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
        /// 求double列表中每项相加之和
        /// </summary>
        /// <param name="numbers"></param>
        /// <returns>double</returns>
        public double Total(List<double> numbers)
        {
            double num = 0;
            foreach (double item in numbers)
            {
                num += item;
            }
            return num;
        }

        /// <summary>
        /// 根据浮点数列表返回颜色值列表
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public List<Color> RainbowColors(List<double> x)
        {
            GH_Gradient gh_Gradient = GH_Gradient.Spectrum();
            gh_Gradient.Linear = true;
            List<Color> colorList = new List<Color>();
            for (int i = 0; i < x.Count; i++)
            {
                colorList.Add(gh_Gradient.ColourAt(x[i]));
            }
            return colorList;
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
            get { return new Guid("26d894de-f241-4ef3-9cb6-ec7d2db730cb"); }
        }
    }
}