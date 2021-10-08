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
            pManager.AddPointParameter("VolumeNodePoints", "VolumeNodePoints", "能够代表volumeNode节点的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePointsExceptNEWS", "BoundaryNodePointsExceptNEWS", "能够代表除了NEWS点外的其他BoundaryNode节点的抽象点", GH_ParamAccess.list);
            pManager.AddIntegerParameter("VolumeConnectivityTree", "ConnectivityTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "AdjacencyTree", "体量邻近关系树", GH_ParamAccess.tree);
            pManager.AddGenericParameter("VolumeLabelList", "LabelList", "体量标签列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("VolumeAttributeList", "AttributeList", "体量属性列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("BoundaryLabelListExceptNEWS", "BoundaryLabelList", "除了NEWS外的其他BoundaryNode的标签列表", GH_ParamAccess.list);
            // pManager.AddGenericParameter("BoundaryAttributeListExceptNEWS", "BoundaryAttributeList", "除了NEWS外的其他BoundaryNode的属性列表", GH_ParamAccess.list);
            pManager[1].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
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
            pManager.AddIntegerParameter("index", "index", "", GH_ParamAccess.list);


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

            List<string> boundaryLabelListExceptNEWS = new List<string>();
            // List<double> boundaryAreaAttributeListExceptNEWS = new List<double>();
            List<NodeAttribute> boundaryNodeAttributeListExceptNEWS = new List<NodeAttribute>();

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


            // double totalArea = 0.0;                                                               // num -> totalArea
            // Plane worldXY = Plane.WorldXY;

            List<Polyline> NEWS_Regions = new List<Polyline>();                             // list5 -> NEWS_Regions
            List<Point3d> NEWS_Vertices = new List<Point3d>();                              // collection -> NEWS_Vertices

            DataTree<int> graph = new DataTree<int>();
            List<Point3d> nodePoints = new List<Point3d>();
            List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();

            //List<NodeAttribute> volumeNodeAttributes = new List<NodeAttribute>();
            //List<NodeAttribute> boundaryNodeAttributes = new List<NodeAttribute>();


            // 获取NEWS四个方向的Region，同时将四个Region的中心添加到节点列表nodeList中
            if (DA.GetDataList<Point3d>("VolumeNodePoints", volumeNodeList))//& DA.GetDataList<Point3d>("BoundaryNodePointsExceptNEWS",boundaryNodeListExceptNEWS)
            {
                // 获取NEWS四个方向的Region，同时将四个Region的中心添加到节点列表nodeList中
                SketchBoxNEWS(PadSquare(volumeNodeList, Plane.WorldXY), ref NEWS_Regions, ref NEWS_Vertices);
                // 按NEWS的顺序，将四个点添加到nodeList末尾
                // nodeList.AddRange(NEWS_Vertices);

                nodePoints.AddRange(volumeNodeList);
                nodePoints.AddRange(NEWS_Vertices);
                nodePoints.AddRange(boundaryNodeListExceptNEWS);

                // 输出NodePoints
                DA.SetDataList("NodePoints", nodePoints);
                //DA.SetDataList("BoundaryNode", NEWS_Vertices);

                // NodeAttribute.BoundaryNodeCount = 4 + boundaryNodeListExceptNEWS.Count;


                
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

                //// 获取NodeAttribute中的NodeArea信息
                //for (int i = 0; i < volumeNodeAttributeList.Count; i++)
                //{
                //    volumeAreaAttributeList.Add(volumeNodeAttributeList[i].NodeArea);
                //}

                //// 计算总面积
                //totalArea = Total(volumeAreaAttributeList);
                //// 面积修正
                //List<double> modifiedAreaAttributeList = new List<double>();
                //for (int i = 0; i < volumeAreaAttributeList.Count; i++)
                //{
                //    modifiedAreaAttributeList.Add(volumeAreaAttributeList[i] / totalArea * totalArea);
                //}

                //// 构造NodeAttribute列表
                //volumeNodeAttributeList.Clear();
                //for (int i = 0; i < volumeLabelList.Count; i++)
                //{
                //    NodeAttribute nodeAttribute = new NodeAttribute(volumeLabelList[i], modifiedAreaAttributeList[i]);
                //    volumeNodeAttributes.Add(nodeAttribute);
                //}

                List<Point3d> allNodeList = volumeNodeList;
                allNodeList.AddRange(NEWS_Vertices);
                Sphere sphere = Sphere.FitSphereToPoints(allNodeList);// 包含NEWS节点的大球

                for (int i = 0; i < NEWS_Vertices.Count; i++)
                {
                    NodeAttribute NEWSNodeAttribute = new NodeAttribute(boundaryLabelList[i], Math.PI * Math.Pow(sphere.Radius, 2) / (double)(4 * allNodeList.Count));
                    boundaryNodeAttributeList.Add(NEWSNodeAttribute);
                }

                for (int i = 0; i < boundaryNodeListExceptNEWS.Count; i++)
                {
                    NodeAttribute boundaryNodeExceptNEWSAttribute = new NodeAttribute(boundaryLabelListExceptNEWS[i], Math.PI * Math.Pow(sphere.Radius, 2) / (double)(4 * allNodeList.Count));
                    boundaryNodeAttributeListExceptNEWS.Add(boundaryNodeExceptNEWSAttribute);
                }

                //DA.SetDataList("VolumeNodeAttributes", volumeNodeAttributes);
                //DA.SetDataList("BoundaryNodeAttributes", boundaryNodeAttributes);

                nodeAttributes.AddRange(volumeNodeAttributeList);
                nodeAttributes.AddRange(boundaryNodeAttributeList);
                nodeAttributes.AddRange(boundaryNodeAttributeListExceptNEWS);

                DA.SetDataList("NodeAttributes", nodeAttributes);

                if (DA.GetDataTree<GH_Integer>("VolumeConnectivityTree", out gh_Structure_ConnectivityTree)
                & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyTree", out gh_Structure_AdjacencyTree))
                {
                    // 将ConnectivityTree从GH_Structure<GH_Integer>转化为DataTree<int>
                    UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_ConnectivityTree, ref volumeConnectivityTree);
                    // 将AdjacencyTree从GH_Structure<GH_Integer>转化为DataTree<int>
                    UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_AdjacencyTree, ref boundaryAdjacencyTree);

                    
                    // 将每个volumeNode的邻接表添加到树形数据中volumeNode对应的那一支上
                    for (int i = 0; i < volumeConnectivityTree.BranchCount; i++)
                    {
                        graph.EnsurePath(i);
                        graph.AddRange(volumeConnectivityTree.Branches[i], graph.Path(i));
                        graph.AddRange(boundaryAdjacencyTree.Branches[i], graph.Path(i));
                    }

                    // 对每个boundary点，向树形数据中添加新的path
                    for (int i = 0; i < NEWS_Vertices.Count + boundaryNodeListExceptNEWS.Count; i++)
                    {
                        graph.EnsurePath(volumeConnectivityTree.BranchCount + i);
                    }

                    // 将每个volume点对于boundary点的邻接表，转化成boundary对于volume的，并且添加到上面新增的boundary的对应path中
                    for (int i = 0; i < boundaryAdjacencyTree.BranchCount; i++)
                    {
                        for (int j = 0; j < boundaryAdjacencyTree.Branch(i).Count; j++)
                        {
                            if (boundaryAdjacencyTree.Branch(i)[j] + NodeAttribute.BoundaryNodeCount < 0)
                            {
                                continue;
                            }
                            graph.Add(i, graph.Path(volumeConnectivityTree.BranchCount + (boundaryAdjacencyTree.Branch(i)[j] + NodeAttribute.BoundaryNodeCount)));
                        }
                    }

                    // 自定义的boundary节点按照角度进行排序，插入到NEWS中，得到所有boundary点之间的邻接关系
                    Sphere NEWSSphere = Sphere.FitSphereToPoints(NEWS_Vertices);
                    Point3d centerPoint = NEWSSphere.Center;

                    List<Point3d> boundaryNodeList = new List<Point3d>();
                    boundaryNodeList.AddRange(NEWS_Vertices);
                    boundaryNodeList.AddRange(boundaryNodeListExceptNEWS);

                    List<int> indexList = new List<int>();

                    List<Point3d> SortedBoundaryNodeList = SortPolyPoints(boundaryNodeList, centerPoint);
                    SortedBoundaryNodeList.Reverse();

                    for (int i = 0; i < SortedBoundaryNodeList.Count; i++)
                    {
                        int index = boundaryNodeList.IndexOf(SortedBoundaryNodeList[i]);
                        indexList.Add(index);
                    }

                    for (int i = 0; i < indexList.Count; i++)
                    {
                        if (i == 0)
                        {
                            graph.Add(indexList[indexList.Count - 1] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                            graph.Add(indexList[i + 1] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                            continue;
                        }
                        if (i == indexList.Count - 1)
                        {
                            graph.Add(indexList[i - 1] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                            graph.Add(indexList[0] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                            continue;
                        }

                        graph.Add(indexList[i - 1] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                        graph.Add(indexList[i + 1] - NodeAttribute.BoundaryNodeCount, graph.Path(volumeConnectivityTree.BranchCount + indexList[i]));
                    }

                    //List<int> indexList = new List<int>();
                    //List<Vector3d> vectorList = new List<Vector3d>();
                    //var keyValueList = new List<KeyValuePair<int, Vector3d>>();

                    //for (int i = 0; i < NEWS_Vertices.Count; i++)
                    //{
                    //    Vector3d vector = new Vector3d(NEWS_Vertices[i]) - new Vector3d(centerPoint);
                    //    vector.Unitize();
                    //    vectorList.Add(vector);

                    //    keyValueList.Add(new KeyValuePair<int, Vector3d>(i, vector));
                    //}
                    //for (int i = 0; i < boundaryNodeListExceptNEWS.Count; i++)
                    //{
                    //    Vector3d vector = new Vector3d(boundaryNodeListExceptNEWS[i]) - new Vector3d(centerPoint);
                    //    vector.Unitize();
                    //    vectorList.Add(vector);

                    //    keyValueList.Add(new KeyValuePair<int, Vector3d>(i + NEWS_Vertices.Count, vector));
                    //}

                    //keyValueList.Sort((a, b) => (a.Value.CompareTo(b.Value)));
                    //indexList = keyValueList.Select(a => a.Key).ToList();




                    //vectorList.Sort();

                    
                    //for (int i = 0; i < vectorList.Count; i++)
                    //{
                    //    int index = vectorList.IndexOf(vectorList[i]);
                    //    indexList.Add(index);
                    //}


                    DA.SetDataList("index", indexList);
                    DA.SetDataTree(0, graph);

                    //Dictionary<int, Vector3d> dictionary = new Dictionary<int, Vector3d>();
                    //Vector3d vector;
                    //for (int i = 0; i < NEWS_Vertices.Count; i++)
                    //{
                    //    vector = new Vector3d(NEWS_Vertices[i]) - new Vector3d(centerPoint);
                    //    vector.Unitize();
                    //    dictionary.Add(i, vector);
                    //}


                    //for (int i = 0; i < boundaryNodeListExceptNEWS.Count; i++)
                    //{
                    //    vector = new Vector3d(boundaryNodeListExceptNEWS[i]) - new Vector3d(centerPoint);
                    //    vector.Unitize();
                    //    dictionary.Add(4 + i, vector);
                    //}

                    //List<KeyValuePair<int, Vector3d>> list = new List<KeyValuePair<int, Vector3d>>(dictionary);
                    //list.Sort(
                    //    delegate (KeyValuePair<int, Vector3d> s1, KeyValuePair<int, Vector3d> s2)
                    //    {
                    //        return s2.Value.CompareTo(s1.Value);
                    //    });
                    //dictionary.Clear();
                    //foreach (KeyValuePair<int,Vector3d> ss in list)
                    //{
                    //    dictionary[ss.Key] = ss.Value;
                    //}



                    //DA.SetDataTree(0, graph);
                }
            }


            // Component中属性的设置
            // 暂不知道Fontsize作用
            //Fontsize = Sphere.FitSphereToPoints(nodeList).Radius / 10.0;

            //List<Brep> NEWS_Brep = new List<Brep>();                                            // list6 -> NEWS_Brep
            //for (int i = 0; i < 4; i++)
            //{
            //    Brep item = Brep.CreatePlanarBreps(NEWS_Regions[i].ToNurbsCurve(), RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)[0];
            //    NEWS_Brep.Add(item);
            //}

            //NEWS_Polygons.AddRange(NEWS_Regions);
            //NEWS_Surfaces.AddRange(NEWS_Brep);
            //NEWS_Tags.AddRange(new string[]
            //{
            //    "N",
            //    "E",
            //    "W",
            //    "S"
            //});
            //NEWS_Location.AddRange(NEWS_Vertices);
            //NEWS_Size = Math.Sqrt(PadSquare(nodeList, worldXY).Area) / 16.0;
            //BasePlane = worldXY;

            //if (volumeLabelList.Count == nodeList.Count - 4)
            //{
            //    Texts.AddRange(volumeLabelList);
            //    Texts.AddRange(new string[]
            //    {
            //        "",
            //        "",
            //        "",
            //        ""
            //    });
            //    Locations.AddRange(nodeList);
            //}
            //else
            //{
            //    List<string> debugIndexList = new List<string>();                                    // list7 -> 
            //    for (int i = 0; i <= nodeList.Count - 5; i++)
            //    {
            //        debugIndexList.Add(i.ToString());
            //    }
            //    Texts.AddRange(debugIndexList);
            //    Texts.AddRange(new string[]
            //    {
            //        "",
            //        "",
            //        "",
            //        ""
            //    });
            //    Locations.AddRange(nodeList);
            //}


            // 主体计算部分
            // Convert GH_Interger to Int, and return a new DataTree<int> volumeConnectivityAttributeTree
            //if (DA.GetDataTree<GH_Integer>(1, out gh_Structure_ConnectivityTree) & gh_Structure_ConnectivityTree.DataCount > 1 & DA.GetDataTree<GH_Integer>(2, out gh_Structure_AdjacencyTree) & gh_Structure_AdjacencyTree.DataCount > 1)
            //{
            //    // 将ConnectivityTree从GH_Structure<GH_Integer>转化为DataTree<int>
            //    //foreach (GH_Path gh_Path in gh_Structure_ConnectivityTree.Paths)
            //    //{
            //    //    List<int> branchItems = new List<int>();
            //    //    foreach (GH_Integer element in gh_Structure_ConnectivityTree.get_Branch(gh_Path))
            //    //    {
            //    //        branchItems.Add(Convert.ToInt32(element.Value));
            //    //    }
            //    //    volumeConnectivityTree.AddRange(branchItems, gh_Path);
            //    //}
            //    //// 将AdjacencyTree从GH_Structure<GH_Integer>转化为DataTree<int>
            //    //foreach (GH_Path gh_Path in gh_Structure_AdjacencyTree.Paths)
            //    //{
            //    //    List<int> branchItems = new List<int>();
            //    //    foreach (GH_Integer element in gh_Structure_AdjacencyTree.get_Branch(gh_Path))
            //    //    {
            //    //        branchItems.Add(Convert.ToInt32(element.Value));
            //    //    }
            //    //    // 
            //    //}

            //    ///* 
            //    //需要修改，原来是List<Line>输入，我这里只是DataTree<int> 
            //    //*/
            //    //allEdgesList = ConnectivityTreeToEdges(nodeList, volumeConnectivityTree);





            //    //NodeLinkToNEWS(ref allEdgesList, NEWS_Regions);
            //    //List<Line> uniqueLineList = CleanLinesUD(allEdgesList);
            //    //Point3dList nodePoint3dList = new Point3dList();                                // point3dList放置nodeList

            //    //// DA.SetDataList(1, nodeList);// nodeList变量的值中包含NEWS节点
                

            //    //// 建好类似DataTree的结构
            //    //List<List<int>> eachNodeConnectivity = new List<List<int>>();                              // list8 -> 每个节点与其他哪些节点相连
            //    //List<List<double>> eachNodeConnectivityLength = new List<List<double>>();                  // list9 -> 与list8对应的连接的线段长度
            //    //for (int i = 0; i < nodeList.Count; i++)
            //    //{
            //    //    eachNodeConnectivity.Add(new List<int>());
            //    //    eachNodeConnectivityLength.Add(new List<double>());
            //    //}

            //    //nodePoint3dList.Clear();
            //    //nodePoint3dList.AddRange(nodeList);

            //    //double modeAbsoluteTolerance = RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;

            //    //// 向建好的类似DataTree结构中塞入数据
            //    //// 因为Edge的输入是没有顺序的，这里是通过下面的计算，找到edge对应的起始点及对应关系
            //    //foreach (Line line in uniqueLineList)
            //    //{
            //    //    if (EpsilonContains(nodePoint3dList, line.From, modeAbsoluteTolerance) & EpsilonContains(nodePoint3dList, line.To, modeAbsoluteTolerance))
            //    //    {
            //    //        int fromClosestIndex = nodePoint3dList.ClosestIndex(line.From);         // num2 -> fromClosestIndex
            //    //        int toClosestIndex = nodePoint3dList.ClosestIndex(line.To);             // num3 -> toClosestIndex
            //    //        if (!eachNodeConnectivity[fromClosestIndex].Contains(toClosestIndex) & fromClosestIndex != toClosestIndex)
            //    //        {
            //    //            eachNodeConnectivity[fromClosestIndex].Add(toClosestIndex);
            //    //            eachNodeConnectivityLength[fromClosestIndex].Add(line.Length);
            //    //        }
            //    //        if (!eachNodeConnectivity[toClosestIndex].Contains(fromClosestIndex) & toClosestIndex != fromClosestIndex)
            //    //        {
            //    //            eachNodeConnectivity[toClosestIndex].Add(fromClosestIndex);
            //    //            eachNodeConnectivityLength[toClosestIndex].Add(line.Length);
            //    //        }

            //    //    }
            //    //}

            //    //// 将连接关系变为DataTree形式输出
            //    //DataTree<int> dataTree = new DataTree<int>();                               // dataTree -> 
            //    //for (int i = 0; i < eachNodeConnectivity.Count; i++)
            //    //{
            //    //    // 按照Length的大小，给nodeConnectivity排序
            //    //    Array.Sort(eachNodeConnectivityLength[i].ToArray(), eachNodeConnectivity[i].ToArray(), Comparer<double>.Default);
            //    //}
            //    //for (int i = 0; i < eachNodeConnectivity.Count; i++)
            //    //{
            //    //    // 先尝试将datatree的分支补满（如果有空分支的话）
            //    //    dataTree.EnsurePath(new int[]
            //    //    {
            //    //        i
            //    //    });
            //    //    if (eachNodeConnectivity[i].Count > 0)
            //    //    {
            //    //        // 生成用来表示连接关系的dataTree
            //    //        dataTree.Branch(i).AddRange(eachNodeConnectivity[i]);
            //    //    }
            //    //    else
            //    //    {
            //    //        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "There is a isolated point in the input. 输入的点中有一个孤立的点");
            //    //    }
            //    //}
            //    //DA.SetDataTree(0, dataTree);


            //    List<object> list10 = new List<object>();                                   // list10 -> 
                
            //    // double num4 = sphere.Radius / 10.0;                                         // num4 -> 
            //    List<double> valuesForColor = new List<double>();                           // list11 ->valuesForColor
            //    for (int i = 1; i <= nodeList.Count - 4; i++)
            //    {
            //        valuesForColor.Add((double)i / (double)(nodeList.Count - 4 + 1));
            //    }
                
            //    // 如果标签列表是空的，或者数量跟nodeList匹配不上，就让标签列表变成序号式的，0,1,2,3...
            //    if (volumeLabelList == null | volumeLabelList.Count != nodeList.Count - 4)
            //    {
            //        for (int i = 0; i <= nodeList.Count - 5; i++)
            //        {
            //            volumeLabelList.Add(Convert.ToString(i));
            //        }
            //    }
            //    // 如果面积列表是空的，或者数量跟nodeList匹配不上，就让面积列表全部变成 Pi * R方 * 0.5 * nodeList.Count
            //    if (volumeAreaAttributeList == null | volumeAreaAttributeList.Count != nodeList.Count - 4)
            //    {
            //        for (int i = 0; i <= nodeList.Count - 5; i++)
            //        {
            //            volumeAreaAttributeList.Add(3.141592653589793 * Math.Pow(sphere.Radius, 2.0) / (double)(2 * nodeList.Count));
            //        }
            //    }
                


            //    volumeLabelList.AddRange(new string[]
            //    {
            //        "N",
            //        "E",
            //        "W",
            //        "S"
            //    });
            //    list10.Add(volumeLabelList);
            //    double totalArea1 = Total(volumeAreaAttributeList);                          // num6 -> totalArea1
            //    List<double> list12 = new List<double>();                                    // list12 -> 
            //    for (int i = 0; i < volumeAreaAttributeList.Count; i++)
            //    {
            //        list12.Add(volumeAreaAttributeList[i] / totalArea1 * totalArea);
            //    }
            //    for (int i = 0; i < 4; i++)
            //    {
            //        list12.Add(3.141592653589793 * Math.Pow(sphere.Radius, 2.0) / (double)(4 * nodeList.Count));
            //    }


            //    //List<Color> list13 = RainbowColors(valuesForColor);
            //    //for (int i = 0; i < 4; i++)
            //    //{
            //    //    list13.Add(Color.Black);
            //    //}
            //    //list10.Add(list12);
            //    //list10.Add(list13);
            //    //List<object> list14 = new List<object>();
            //    //List<string> item2 = new List<string>();
            //    //List<double> item3 = new List<double>();
            //    //List<Color> item4 = new List<Color>();
            //    //list14.Add(item2);
            //    //list14.Add(item3);
            //    //list14.Add(item4);
            //    //DA.SetDataList(2, list10);
            //}


        }

        public List<Point3d> SortPolyPoints(List<Point3d> vPoints, Point3d center)
        {
            List<Point3d> cloneList = new List<Point3d>();
            cloneList.AddRange(vPoints);


            if (cloneList == null || cloneList.Count == 0)
            {
                return null;
            }

            for (int i = 0; i < cloneList.Count - 1; i++)
            {
                for (int j = 0; j < cloneList.Count - i - 1; j++)
                {
                    bool flag = PointCmp(cloneList[j], cloneList[j + 1], center);
                    if (flag)
                    {
                        Point3d tmp = cloneList[j];
                        cloneList[j] = cloneList[j + 1];
                        cloneList[j + 1] = tmp;
                    }
                }
            }
            return cloneList;
        }

        /// <summary>
        /// 若点a大于点b，即点a在点b的顺时针方向，返回true，否则返回false
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        private bool PointCmp(Point3d a, Point3d b, Point3d center)
        {
            Vector3d vectorOA = new Vector3d(a) - new Vector3d(center);
            Vector3d vectorOB = new Vector3d(b) - new Vector3d(center);

            if (vectorOA.X >= 0 && vectorOB.X < 0)
            {
                return true;
            }
            else if (vectorOA.X == 0 && vectorOB.X == 0)
            {
                return vectorOA.Y > vectorOB.Y;
            }
            // 向量OA和向量OB的叉积
            

            // double det = (a.X - center.X) * (b.Y - center.Y) - (b.X - center.X) * (a.Y - center.Y);
            double cross = vectorOA.X * vectorOB.Y - vectorOA.Y * vectorOB.X;

            Vector3d vector = Vector3d.CrossProduct(vectorOA, vectorOB);

            if (vector.Z < 0)
            {
                return true;
            }
            if (vector.Z > 0)
            {
                return false;
            }

            // 这里出了问题，E和W，叉积为0，并且两个的d是一样的
            double d1 = Math.Pow((a.X - center.X), 2.0) + Math.Pow((a.Y - center.Y), 2.0);
            double d2 = Math.Pow((b.X - center.X), 2.0) + Math.Pow((b.Y - center.Y), 2.0);
            return d1 > d2;
        }

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