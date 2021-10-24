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
            // pManager.AddPointParameter("GraphNodePoints", "GNodePoints", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            // pManager.AddGenericParameter("GraphNodeAttributes", "GNodeAttributes", "图结构中节点所包含的属性", GH_ParamAccess.list);
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
            List<Point3d> volumeNodeList = globalParameter.VolumeNodePointLocations.ToList<Point3d>();
            List<Point3d> boundaryNodeList = globalParameter.BoundaryNodePointLocations.ToList<Point3d>();


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

            // 输出Graph的树形结构
            if (DA.GetDataTree<GH_Integer>("VolumeConnectivityTree", out gh_Structure_ConnectivityTree)
            & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyTree", out gh_Structure_AdjacencyTree))
            {
                Point3d centerPoint = UtilityFunctions.CalCenterPoint(boundaryNodeList);

                List<int> sortedBoundaryNodeIndexList = new List<int>();

                // 对边界节点进行排序
                List<Point3d> SortedBoundaryNodeList = UtilityFunctions.SortPolyPoints(boundaryNodeList, centerPoint);

                for (int i = 0; i < SortedBoundaryNodeList.Count; i++)
                {
                    int boundaryPointsIndex = boundaryNodeList.IndexOf(SortedBoundaryNodeList[i]);
                    sortedBoundaryNodeIndexList.Add(boundaryPointsIndex);
                }

                // 构造包含所有node节点的列表，顺序是 inner node + outer node
                nodePoints.AddRange(volumeNodeList);
                nodePoints.AddRange(SortedBoundaryNodeList);
                // 输出NodePoints
                // DA.SetDataList("GraphNodePoints", nodePoints);

                // 计算每个volume点的面积占比，并写入NodeAttribute中的NodeAreaProportion属性
                UtilityFunctions.CalculateAreaProportion(volumeNodeAttributeList);


                // 构造BoundaryNodeAttribute列表
                for (int i = 0; i < boundaryNodeList.Count; i++)
                {
                    NodeAttribute boundaryNodeAttribute = new NodeAttribute(boundaryLabelList[i]);
                    boundaryNodeAttributeList.Add(boundaryNodeAttribute);
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
                // DA.SetDataList("GraphNodeAttributes", nodeAttributes);


                // 构造Node类的列表，并输出
                for (int i = 0; i < nodePoints.Count; i++)
                {
                    nodes.Add(new Node(nodePoints[i], nodeAttributes[i]));
                }

                DA.SetDataList("GraphNode", nodes);



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