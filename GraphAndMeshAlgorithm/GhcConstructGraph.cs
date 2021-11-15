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
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcConstructGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ConstructGraph class.
        /// </summary>
        public GhcConstructGraph()
          : base("ConstructGraph", "ConstructGraph",
              "构建图结构",
              "VolumeGeneratorBasedOnGraph", "Construct Graph")
        {
            
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // pManager.AddGenericParameter("GlobalParameter", "GP", "全局参数传递", GH_ParamAccess.item);

            // pManager.AddIntegerParameter("VolumeConnectivityTree", "CTree", "体量连接关系树", GH_ParamAccess.tree);
            // pManager.AddIntegerParameter("BoundaryAdjacencyTree", "ATree", "体量邻近关系树", GH_ParamAccess.tree);

            // pManager.AddGenericParameter("VolumeNodeAttributes", "VNAttributes", "体量属性列表", GH_ParamAccess.list);

            pManager.AddGenericParameter("CSVFilePath", "Path", "CSV文件的路径", GH_ParamAccess.item);

            pManager.AddPointParameter("VolumeNodePoints", "VNPoints", "能够代表volumeNode节点的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePoints", "BNPoints", "能够代表BoundaryNode节点的抽象点（point）", GH_ParamAccess.list);

            pManager.AddGenericParameter("BoundaryNodeLabels", "BNLabels", "每个BoundaryNode所对应的标签组成的列表", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            // pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);
            // pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);

            pManager.AddGenericParameter("Graph", "G", "图结构", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //#region 全局参数传递
            //GlobalParameter globalParameter = new GlobalParameter();
            //DA.GetData("GlobalParameter", ref globalParameter);
            //int volumeNodeCount = globalParameter.VolumeNodeCount;
            //int boundaryNodeCount = globalParameter.BoundaryNodeCount;
            //List<Point3d> volumeNodeList = globalParameter.VolumeNodePointLocations.ToList<Point3d>();
            //List<Point3d> boundaryNodeList = globalParameter.BoundaryNodePointLocations.ToList<Point3d>();
            //#endregion

            string csvPath = null;
            List<Point3d> volumeNodes = new List<Point3d>();
            List<Point3d> boundaryNodes = new List<Point3d>();
            List<string> boundaryLabelList = new List<string>();

            int volumeNodeCount = 0;
            int boundaryNodeCount = 0;

            if (DA.GetData("CSVFilePath", ref csvPath)
                && DA.GetDataList<Point3d>("VolumeNodePoints", volumeNodes)
                && DA.GetDataList<Point3d>("BoundaryNodePoints", boundaryNodes)
                && DA.GetDataList<string>("BoundaryNodeLabels", boundaryLabelList))
            {
                volumeNodeCount = volumeNodes.Count;
                boundaryNodeCount = boundaryNodes.Count;

                #region 读取CSV数据
                #region Initialize variables 初始化变量
                List<string> innerNodeLabelList = new List<string>();
                List<NodeAttribute> innerNodeAttributesList = new List<NodeAttribute>();
                DataTree<int> volumeConnectivityTree = new DataTree<int>();
                DataTree<int> boundaryAdjacencyTree = new DataTree<int>();
                #endregion

                #region Read the data of the CSV files as individual lines 按行读取csv文件中的数据
                string[] csvLines = null;
                try
                {
                    
                    csvLines = System.IO.File.ReadAllLines(csvPath);
                    
                }
                catch (System.IO.DirectoryNotFoundException)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "没有找到CSV文件，请检查文件路径");
                    return;
                }
                #endregion

                #region Parse all lines 解析所有行的数据
                for (int i = 1; i < csvLines.Length; i++)
                {
                    #region Split rowData with "," 按逗号分隔每一行的字符串
                    string[] rowData = csvLines[i].Split(',');
                    #endregion

                    #region 读取时可能出现的报错
                    if (rowData[0] == "")
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "节点必须要有标签，请检查CSV文件");
                        return;
                    }
                    if (rowData[1] == "")
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "节点必须要有属性，请检查CSV文件");
                        return;
                    }
                    #endregion

                    #region Add first and second colume data into corresponding list 将得到的第一列和第二列数据分别添加到对应的列表中
                    innerNodeLabelList.Add(rowData[0]);
                    innerNodeAttributesList.Add(new NodeAttribute(rowData[0],
                                                             Convert.ToDouble(rowData[1]),
                                                             Convert.ToInt32(rowData[4]),
                                                             Convert.ToInt32(rowData[5]),
                                                             Convert.ToInt32(rowData[6]),
                                                             Convert.ToDouble(rowData[7]),
                                                             Convert.ToDouble(rowData[8]),
                                                             Convert.ToDouble(rowData[9])
                                                             ));
                    #endregion

                    #region Split ConnectivityArribute with "-" and add it to corresponding datatree with correct path 用-号将第三列数据分隔，并按照对应的路径添加到树形结构中
                    string[] connectivity = rowData[2].Split('-');
                    GH_Path path = new GH_Path(0, i - 1);
                    for (int j = 0; j < connectivity.Length; j++)
                    {
                        volumeConnectivityTree.EnsurePath(path);

                        if (connectivity[j] == "")
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "节点必须有连接关系, 请检查CSV文件。");
                            return;
                        }
                        if (Convert.ToInt32(connectivity[j]) >= csvLines.Length - 1)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, connectivity[j].ToString() + "节点的连接对象的序号不在节点列表中。");
                            return;
                        }
                        else
                        {
                            volumeConnectivityTree.Add(Convert.ToInt32(connectivity[j]), path);
                        }
                    }
                    innerNodeAttributesList[i - 1].ConnectivityTable = volumeConnectivityTree.Branch(path).ToArray();
                    #endregion

                    #region Split AjacencyArribute with "-" and add it to corresponding datatree with correct path 用-号将第四列数据分隔，并按照对应的路径添加到树形结构中
                    string[] adjacency = rowData[3].Split('-');
                    for (int j = 0; j < adjacency.Length; j++)
                    {
                        boundaryAdjacencyTree.EnsurePath(path);

                        // 如果是空，就跳过它继续后面的循环
                        if (adjacency[j] == "")
                        {
                            boundaryAdjacencyTree.AddRange(new List<int>(), path);
                            continue;
                        }
                        else if (Convert.ToInt32(adjacency[j]) >= boundaryLabelList.Count)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查adjacency中序号的输入，是否超过了BoundaryNode数量的上限。");
                            return;
                        }
                        else if (Convert.ToInt32(adjacency[j]) < 0)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查adjacency中序号的输入，是否超过了BoundaryNode数量的下限。");
                            return;
                        }
                        else
                        {
                            boundaryAdjacencyTree.Add(Convert.ToInt32(adjacency[j]), path);
                        }
                    }
                    innerNodeAttributesList[i - 1].AdjacencyTable = boundaryAdjacencyTree.Branch(path).ToArray();
                    #endregion
                }
                #endregion

                #region 排除数量不等的错误
                if (volumeConnectivityTree.BranchCount != volumeNodes.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的体量点的数量，与输入的代表体量节点的点的数量不相等");
                    return;
                }
                if (boundaryLabelList.Count != boundaryNodes.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的边界点的数量，与输入的代表边界节点的点的数量不相等");
                    return;
                }
                #endregion

                #endregion



                #region 构建GraphLoL

                #region 局部变量初始化
                List<List<int>> boundaryNodeConnectivityTable = new List<List<int>>();
                List<List<int>> boundaryNodeAdjacencyTable = new List<List<int>>();

                DataTree<int> graphDT = new DataTree<int>();
                #endregion

                #region 对边界节点进行排序
                Point3d centerPoint = UtilityFunctions.CalCenterPoint(boundaryNodes);
                List<int> sortedBoundaryNodeIndexList = new List<int>();
                List<Point3d> SortedBoundaryNodes = UtilityFunctions.SortPolyPoints(boundaryNodes, centerPoint);

                // 根据point的排序，获取对应的index的排序
                for (int i = 0; i < SortedBoundaryNodes.Count; i++)
                {
                    int boundaryPointsIndex = boundaryNodes.IndexOf(SortedBoundaryNodes[i]);
                    sortedBoundaryNodeIndexList.Add(boundaryPointsIndex);
                }
                #endregion

                #region 将每个VolumeNode对于其他 inner + outer 点的邻接关系，添加到每个VolumeNode分支上
                for (int i = 0; i < volumeConnectivityTree.BranchCount; i++)
                {
                    graphDT.EnsurePath(i);
                    graphDT.AddRange(volumeConnectivityTree.Branches[i], graphDT.Path(i));
                    for (int j = 0; j < boundaryAdjacencyTree.Branches[i].Count; j++)
                    {
                        // 注意这里add的item应该是排过序后的NEWS对应的index
                        graphDT.Add(sortedBoundaryNodeIndexList.IndexOf(boundaryAdjacencyTree.Branch(i)[j]) + volumeNodeCount, graphDT.Path(i));
                    }
                }
                #endregion

                // 对每个boundary点（outer），向树形数据中添加新的path（BoundaryNode分支）
                for (int i = 0; i < boundaryNodeCount; i++)
                {
                    graphDT.EnsurePath(volumeNodeCount + i);
                }

                #region 每个BoundaryNode对于其他 inner 点的邻接关系，添加到每个BoundaryNode分支上，i相当于是volumeNode的序号了
                for (int i = 0; i < boundaryAdjacencyTree.BranchCount; i++)
                {
                    for (int j = 0; j < boundaryAdjacencyTree.Branch(i).Count; j++)
                    {
                        graphDT.Add(i, graphDT.Path(volumeNodeCount + (sortedBoundaryNodeIndexList.IndexOf(boundaryAdjacencyTree.Branch(i)[j]))));
                    }
                }

                for (int i = 0; i < sortedBoundaryNodeIndexList.Count; i++)
                {
                    boundaryNodeConnectivityTable.Add(new List<int>());
                    for (int j = 0; j < boundaryAdjacencyTree.BranchCount; j++)
                    {
                        if (boundaryAdjacencyTree.Branch(j).Contains(i))
                        {
                            boundaryNodeConnectivityTable[i].Add(j);
                        }
                    }
                }
                #endregion

                #region 每个BoundaryNode对于其他 outer 点的邻接关系（outer点顺时针排序，每个outer点与其左右的两个outer点分别连接），添加到每个BoundaryNode分支上
                for (int i = 0; i < sortedBoundaryNodeIndexList.Count; i++)
                {
                    boundaryNodeAdjacencyTable.Add(new List<int>());

                    // 对第一个点的特殊处理
                    if (i == 0)
                    {
                        graphDT.Add((sortedBoundaryNodeIndexList.Count - 1) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));
                        graphDT.Add((i + 1) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));

                        boundaryNodeAdjacencyTable[i].Add((sortedBoundaryNodeIndexList.Count - 1) + volumeNodeCount);
                        boundaryNodeAdjacencyTable[i].Add((i + 1) + volumeNodeCount);
                        continue;
                    }
                    // 对最后一个点的特殊处理
                    if (i == sortedBoundaryNodeIndexList.Count - 1)
                    {
                        graphDT.Add((i - 1) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));
                        graphDT.Add((0) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));

                        boundaryNodeAdjacencyTable[i].Add((i - 1) + volumeNodeCount);
                        boundaryNodeAdjacencyTable[i].Add((0) + volumeNodeCount);
                        continue;
                    }

                    // 对普通点的处理
                    graphDT.Add((i - 1) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));
                    graphDT.Add((i + 1) + volumeNodeCount, graphDT.Path(volumeNodeCount + i));
                    boundaryNodeAdjacencyTable[i].Add((i - 1) + volumeNodeCount);
                    boundaryNodeAdjacencyTable[i].Add((i + 1) + volumeNodeCount);
                }
                #endregion

                List<List<int>> graphLoL = UtilityFunctions.DataTreeToLoL<int>(graphDT);

                #endregion

                #region 构建GraphNode

                #region 局部变量初始化
                List<NodeAttribute> outerNodeAttributeList = new List<NodeAttribute>();

                List<Point3d> nodePoints = new List<Point3d>();
                List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();

                List<GraphNode> graphNodes = new List<GraphNode>();
                #endregion

                #region 构造包含所有nodePoint的列表，顺序是 inner node + outer node
                nodePoints.AddRange(volumeNodes);
                nodePoints.AddRange(SortedBoundaryNodes);
                #endregion

                #region 构造包含所有nodeAttributes的列表，顺序是inner node + outer node
                // 计算每个volume点的面积占比，并写入NodeAttribute中的NodeAreaProportion属性
                UtilityFunctions.CalculateAreaProportion(innerNodeAttributesList);

                // 构造BoundaryNodeAttribute列表
                for (int i = 0; i < boundaryNodeCount; i++)
                {
                    NodeAttribute boundaryNodeAttribute = new NodeAttribute(boundaryLabelList[i]);
                    boundaryNodeAttribute.ConnectivityTable = boundaryNodeConnectivityTable[i].ToArray();
                    boundaryNodeAttribute.AdjacencyTable = boundaryNodeAdjacencyTable[i].ToArray();
                    outerNodeAttributeList.Add(boundaryNodeAttribute);
                }

                // 构造包含所有node节点的属性列表，顺序是 inner node + outer node
                nodeAttributes.AddRange(innerNodeAttributesList);
                nodeAttributes.AddRange(outerNodeAttributeList);
                #endregion

                #region 构造Node类的列表，并输出
                for (int i = 0; i < nodePoints.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        graphNodes.Add(new GraphNode(nodePoints[i], nodeAttributes[i], true));
                    }
                    if (i >= volumeNodeCount && i < nodePoints.Count)
                    {
                        graphNodes.Add(new GraphNode(nodePoints[i], nodeAttributes[i], false));
                    }
                }

                #endregion

                Graph graph = new Graph(graphNodes, graphLoL);
                DA.SetData("Graph", graph);

                #endregion
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
            get { return new Guid("95f16532-4629-4afd-8227-75f14332553e"); }
        }
    }
}