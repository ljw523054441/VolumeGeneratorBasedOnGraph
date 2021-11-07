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

namespace VolumeGeneratorBasedOnGraph
{
    public class ConstructGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ConstructGraph class.
        /// </summary>
        public ConstructGraph()
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
            pManager.AddGenericParameter("GlobalParameter", "GP", "全局参数传递", GH_ParamAccess.item);

            pManager.AddIntegerParameter("VolumeConnectivityTree", "CTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "ATree", "体量邻近关系树", GH_ParamAccess.tree);

            pManager.AddGenericParameter("VolumeNodeAttributes", "VNAttributes", "体量属性列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("BoundaryNodeLabels", "BNLabels", "每个BoundaryNode所对应的标签组成的列表", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Graph", "G", "描述VolumeNode和BoundaryNode的所有连接关系的图结构", GH_ParamAccess.tree);
            pManager.AddGenericParameter("GraphNode", "GNode", "图结构中的节点", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int volumeNodeCount = globalParameter.VolumeNodeCount;
            int boundaryNodeCount = globalParameter.BoundaryNodeCount;
            List<Point3d> volumeNodeList = globalParameter.VolumeNodePointLocations.ToList<Point3d>();
            List<Point3d> boundaryNodeList = globalParameter.BoundaryNodePointLocations.ToList<Point3d>();
            #endregion

            #region 局部变量初始化
            GH_Structure<GH_Integer> gh_Structure_ConnectivityTree = null;
            GH_Structure<GH_Integer> gh_Structure_AdjacencyTree = null;
            DataTree<int> volumeConnectivityTree = new DataTree<int>();
            DataTree<int> boundaryAdjacencyTree = new DataTree<int>();

            List<string> volumeLabelList = new List<string>();
            List<NodeAttribute> volumeNodeAttributeList = new List<NodeAttribute>();
            List<string> boundaryLabelList = new List<string>();
            List<NodeAttribute> boundaryNodeAttributeList = new List<NodeAttribute>();

            DataTree<int> graph = new DataTree<int>();
            List<Point3d> nodePoints = new List<Point3d>();
            List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();

            List<Node> nodes = new List<Node>();
            #endregion

            #region 可能的报错
            if (DA.GetDataList<NodeAttribute>("VolumeNodeAttributes", volumeNodeAttributeList) && volumeNodeAttributeList.Count != volumeNodeCount)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "体量节点的数量跟体量节点属性的数量要一致！");
                return;
            }
            if (DA.GetDataList<string>("BoundaryNodeLabels", boundaryLabelList) && boundaryLabelList.Count != boundaryNodeCount)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "边界节点的数量与边界节点的标签数量要一致！");
                return;
            }
            #endregion

            // 输出Graph的树形结构
            if (DA.GetDataTree<GH_Integer>("VolumeConnectivityTree", out gh_Structure_ConnectivityTree)
            & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyTree", out gh_Structure_AdjacencyTree))
            {
                #region 对边界节点进行排序
                Point3d centerPoint = UtilityFunctions.CalCenterPoint(boundaryNodeList);
                List<int> sortedBoundaryNodeIndexList = new List<int>();
                List<Point3d> SortedBoundaryNodeList = UtilityFunctions.SortPolyPoints(boundaryNodeList, centerPoint);

                // 根据point的排序，获取对应的index的排序
                for (int i = 0; i < SortedBoundaryNodeList.Count; i++)
                {
                    int boundaryPointsIndex = boundaryNodeList.IndexOf(SortedBoundaryNodeList[i]);
                    sortedBoundaryNodeIndexList.Add(boundaryPointsIndex);
                }
                #endregion

                #region 邻接表转化为包括所有node关系的图结构
                #region 类型转换：GH_Structure<GH_Integer>转化为DataTree<int>
                // 将ConnectivityTree从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_ConnectivityTree, ref volumeConnectivityTree);
                // 将AdjacencyTree从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_AdjacencyTree, ref boundaryAdjacencyTree);
                #endregion

                #region 将每个VolumeNode对于其他 inner + outer 点的邻接关系，添加到每个VolumeNode分支上
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
                #endregion

                List<List<int>> boundaryNodeConnectivityTable = new List<List<int>>();
                List<List<int>> boundaryNodeAdjacencyTable = new List<List<int>>();

                // 对每个boundary点（outer），向树形数据中添加新的path（BoundaryNode分支）
                for (int i = 0; i < boundaryNodeCount; i++)
                {
                    graph.EnsurePath(volumeNodeCount + i);
                }

                #region 每个BoundaryNode对于其他 inner 点的邻接关系，添加到每个BoundaryNode分支上，i相当于是volumeNode的序号了
                for (int i = 0; i < boundaryAdjacencyTree.BranchCount; i++)
                {
                    for (int j = 0; j < boundaryAdjacencyTree.Branch(i).Count; j++)
                    {
                        graph.Add(i, graph.Path(volumeNodeCount + (sortedBoundaryNodeIndexList.IndexOf(boundaryAdjacencyTree.Branch(i)[j]))));
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
                        
                        
                        //for (int k = 0; k < boundaryAdjacencyTree.Branch(j).Count; k++)
                        //{
                            
                        //    //if (boundaryAdjacencyTree.Branch(j)[k] == i)
                        //    //{
                        //    //    boundaryNodeConnectivityTable[i].Add(j);
                        //    //}
                        //}
                        
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
                        graph.Add((sortedBoundaryNodeIndexList.Count - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                        graph.Add((i + 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));

                        boundaryNodeAdjacencyTable[i].Add((sortedBoundaryNodeIndexList.Count - 1) + volumeNodeCount);
                        boundaryNodeAdjacencyTable[i].Add((i + 1) + volumeNodeCount);
                        continue;
                    }
                    // 对最后一个点的特殊处理
                    if (i == sortedBoundaryNodeIndexList.Count - 1)
                    {
                        graph.Add((i - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                        graph.Add((0) + volumeNodeCount, graph.Path(volumeNodeCount + i));

                        boundaryNodeAdjacencyTable[i].Add((i - 1) + volumeNodeCount);
                        boundaryNodeAdjacencyTable[i].Add((0) + volumeNodeCount);
                        continue;
                    }

                    // 对普通点的处理
                    graph.Add((i - 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                    graph.Add((i + 1) + volumeNodeCount, graph.Path(volumeNodeCount + i));
                    boundaryNodeAdjacencyTable[i].Add((i - 1) + volumeNodeCount);
                    boundaryNodeAdjacencyTable[i].Add((i + 1) + volumeNodeCount);
                }
                #endregion
                DA.SetDataTree(0, graph);
                #endregion

                #region 构造包含所有nodePoint的列表，顺序是 inner node + outer node
                nodePoints.AddRange(volumeNodeList);
                nodePoints.AddRange(SortedBoundaryNodeList);
                #endregion

                #region 构造包含所有nodeAttributes的列表，顺序是inner node + outer node
                // 计算每个volume点的面积占比，并写入NodeAttribute中的NodeAreaProportion属性
                UtilityFunctions.CalculateAreaProportion(volumeNodeAttributeList);

                // 构造BoundaryNodeAttribute列表
                for (int i = 0; i < boundaryNodeList.Count; i++)
                {
                    NodeAttribute boundaryNodeAttribute = new NodeAttribute(boundaryLabelList[i]);
                    boundaryNodeAttribute.ConnectivityTable = boundaryNodeConnectivityTable[i].ToArray();
                    boundaryNodeAttribute.AdjacencyTable = boundaryNodeAdjacencyTable[i].ToArray();
                    boundaryNodeAttributeList.Add(boundaryNodeAttribute);
                }

                // 构造包含所有node节点的属性列表，顺序是 inner node + outer node
                nodeAttributes.AddRange(volumeNodeAttributeList);
                nodeAttributes.AddRange(boundaryNodeAttributeList);
                #endregion

                #region 构造Node类的列表，并输出
                for (int i = 0; i < nodePoints.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        nodes.Add(new Node(nodePoints[i], nodeAttributes[i], true));
                    }
                    if (i >= volumeNodeCount && i < nodePoints.Count)
                    {
                        nodes.Add(new Node(nodePoints[i], nodeAttributes[i], false));
                    }
                }
                DA.SetDataList("GraphNode", nodes);
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