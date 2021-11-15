using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.GUI.Gradient;
using Rhino;
using Rhino.Geometry;
using Rhino.Collections;
using System;
using System.Drawing;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class CoinGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CoinGraph class.
        /// </summary>
        public CoinGraph()
          : base("ForceDirectedCoinGraph", "CoinGraph",
              "Finds locations of vertices for a kissing-disk/coin-drawing of a graph. Note that you need to connect the output vertices to a general graph drawing component to actually draw the graph",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddIntegerParameter("ConnectivityGraph", "ConnectivityGraph", "体量连接关系", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyGraph", "AdjacencyGraph", "体量与边界的邻接关系", GH_ParamAccess.tree);
            pManager.AddPointParameter("VolumeNode", "VolumeNode", "用抽象点（point）表示的图结构节点（node）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNode", "BoundaryNode", "用抽象点（point）表示的边界（node）", GH_ParamAccess.list);
            pManager.AddGenericParameter("VolumeNodeAttributes", "VolumeAttributes", "体量属性列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("BoundaryNodeAttributes", "BoundaryAttributes", "边界属性列表", GH_ParamAccess.list);
            pManager.AddNumberParameter("AttractionStrength", "AttractionStrength", "吸引强度", GH_ParamAccess.item, 0.25);
            pManager.AddNumberParameter("RepulsionStrength", "RepulsionStrength", "排斥强度", GH_ParamAccess.item, 1.4);
            pManager.AddNumberParameter("ConvergenceTolerance", "ConvergenceTolerance", "收敛公差 Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 0.001);
            pManager.AddIntegerParameter("MaximumIteration", "MaximumIteration", "最大迭代次数", GH_ParamAccess.item, 10000);
            pManager.AddBooleanParameter("Active", "Active", "是否开启", GH_ParamAccess.item, false);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("NeswVolumeNode", "NewVolumeNode", "A set of new vertices on which a neat coin/kissing disk drawing can be drawn", GH_ParamAccess.list);
            pManager.AddIntegerParameter("IterationCount", "IterationCount", "How many iterations (recusrions) have happend?", GH_ParamAccess.item);
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


            GH_Structure<GH_Integer> gh_Structure_ConnectivityGraph = null;                         // gh_Structure -> gh_Structure_ConnectivityGraph
            GH_Structure<GH_Integer> gh_Structure_AdjacencyGraph = null;                            // 
            DataTree<int> connectivityTree = new DataTree<int>();
            DataTree<int> adjacencyTree = new DataTree<int>();

            List<Point3d> volumeNodeList = new List<Point3d>();                                           // list -> nodeList
            List<Point3d> boundaryNodeList = new List<Point3d>();
            List<GraphNodeAttribute> volumeNodeAttributes = new List<GraphNodeAttribute>();        // list2 -> nodeAttributeList
            List<GraphNodeAttribute> boundaryNodeAttributes = new List<GraphNodeAttribute>();

            double st_Att = 0.0;
            double st_Rep = 0.0;
            double tolerance = 0.01;
            int iteration = 0;                                                        // num -> iteration
            bool activeFlag = false;

            DataTree<Vector3d> attrationForce = new DataTree<Vector3d>();
            DataTree<Vector3d> repulsionForce = new DataTree<Vector3d>();

            if (DA.GetDataTree<GH_Integer>("ConnectivityGraph", out gh_Structure_ConnectivityGraph) 
                & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyGraph", out gh_Structure_AdjacencyGraph) 
                & DA.GetDataList<Point3d>("VolumeNode", volumeNodeList) 
                & DA.GetDataList<GraphNodeAttribute>("VolumeNodeAttributes", volumeNodeAttributes) 
                & DA.GetDataList<Point3d>("BoundaryNode", boundaryNodeList) 
                & DA.GetDataList<GraphNodeAttribute>("BoundaryNodeAttributes", boundaryNodeAttributes))
            {
                // 将ConnectivityTree从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_ConnectivityGraph, ref connectivityTree);
                // 将AdjacencyTree从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_AdjacencyGraph, ref adjacencyTree);

                // 要考虑两个树合并的问题
                connectivityTree.MergeTree(adjacencyTree);

                DA.GetData<double>("AttractionStrength", ref st_Att);
                DA.GetData<double>("RepulsionStrength", ref st_Rep);
                DA.GetData<double>("ConvergenceTolerance", ref tolerance);
                DA.GetData<int>("MaximumIteration", ref iteration);
                DA.GetData<bool>("Active", ref activeFlag);

                // 主体计算部分
                List<Point3d> newVolumeNodeList = new List<Point3d>();                                      // list7 ->
                bool repeatFlag = false;                                                             // flag2 -> repeatFlag

                if (volumeNodeList.Count == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "something went wrong! there are no vertices to operate on!");
                    return;
                }

                List<Point3d> currentIterationNodeList = volumeNodeList;                                                 // list8 -> newNodeList
                int currentIteration = 0;                                                          // num2 -> currentIteration
                do
                {
                    if (activeFlag)
                    {
                        newVolumeNodeList = currentIterationNodeList;
                    }
                    DA.SetData("IterationCount", currentIteration);

                    attrationForce = AttractionForce(connectivityTree, volumeNodeList, boundaryNodeList, volumeNodeAttributes, boundaryNodeAttributes, st_Att, globalParameter);
                    repulsionForce = RepulsionForce(connectivityTree, volumeNodeList, boundaryNodeList, volumeNodeAttributes, boundaryNodeAttributes, st_Rep);

                    Result_Force(
                        volumeNodeList, 
                        attrationForce, 
                        repulsionForce, 
                        tolerance, 
                        ref currentIterationNodeList, 
                        ref repeatFlag);
                    currentIteration++;

                } 
                while (!(!repeatFlag | !activeFlag | currentIteration > iteration));

                if (activeFlag)
                {
                    DA.SetDataList("NewVolumeNode", newVolumeNodeList);
                    return;
                }
                DA.SetDataList("NewVolumeNode", volumeNodeList);
            }
        }

        /// <summary>
        /// 计算作用在每个代表体量volume的点node上的合力
        /// </summary>
        /// <param name="volumeNodeList">每个代表体量的点node的列表</param>
        /// <param name="attractionForce">每个代表体量volume的点node的吸引力列表（包括graph上的也包括NEWS上的）</param>
        /// <param name="repulsionForce">每个代表体量volume的点node的排斥力列表（包括graph上的也包括NEWS上的）</param>
        /// <param name="tolerance">用来判断的公差</param>
        /// <param name="newNodeList">输出的合力作用后的新点列表</param>
        /// <param name="repeat">决定是否重复的布尔值</param>
        /// /// <returns>返回每个点在合力作用后形成的新点的位置</returns>
        public void Result_Force(List<Point3d> volumeNodeList, DataTree<Vector3d> attractionForce, DataTree<Vector3d> repulsionForce, double tolerance, ref List<Point3d> newNodeList, ref bool repeat)
        {
            if (volumeNodeList.Count == 0 | attractionForce.DataCount == 0 | repulsionForce.DataCount == 0 | tolerance == 0.0)
            {
                return;
            }
            // 每个代表volume的点的合力组成的列表
            List<Vector3d> resultantForceList = new List<Vector3d>();
            // 对每个代表volume的点分别求合吸引力，合斥力，最后求合力
            for (int i = 0; i < volumeNodeList.Count; i++)
            {
                Vector3d allAttractionForce = new Vector3d(0,0,0);
                Vector3d allRepulsionForce = new Vector3d(0, 0, 0);
                for (int j = 0; j < attractionForce.Branches[i].Count; j++)
                {
                    allAttractionForce = Vector3d.Add(allAttractionForce, attractionForce.Branches[i][j]);
                }
                for (int j = 0; j < repulsionForce.Branches[i].Count; j++)
                {
                    allRepulsionForce = Vector3d.Add(allRepulsionForce, repulsionForce.Branches[i][j]);
                }

                resultantForceList.Add(Vector3d.Add(allAttractionForce, allRepulsionForce));
            }

            // 将合力作用在每个代表volume的点上，并且返回当前公差来决定是否继续解算
            List<Point3d> movedNodeList= new List<Point3d>();                                         // list2 -> 
            double currentTolerance = 0.0;                                                       // num2 -> currentTolerance
            for (int i = 0; i < volumeNodeList.Count; i++)
            {
                movedNodeList.Add(volumeNodeList[i] + resultantForceList[i]);
                currentTolerance += resultantForceList[i].Length;
            }
            newNodeList = movedNodeList;
            repeat = currentTolerance > tolerance;
        }

        /// <summary>
        /// 计算每个代表体量volume的点node的吸引力列表（包括graph上的也包括NEWS上的）
        /// </summary>
        /// <param name="mergedGraph">每个代表体量的点node的connectivity（与体量点）和adjacent（与边界NEWS）</param>
        /// <param name="volumeNodeList">代表体量volume的点的列表</param>
        /// <param name="boundaryNodeList">代表边界NEWS的点的列表</param>
        /// <param name="volumeNodeAttributes">代表体量volume的点的属性</param>
        /// <param name="boundaryNodeAttributes">代表边界NEWS的点的属性</param>
        /// <param name="st_Att">吸引力强度</param>
        /// <returns>返回每个点（树形数据的一支branch）的volume点对它的吸引力，以及NEWS点对它的吸引力，两者加起来（共有mergeGraph.Branch[i].Count个）</returns>
        public DataTree<Vector3d> AttractionForce(DataTree<int> mergedGraph, List<Point3d> volumeNodeList, List<Point3d> boundaryNodeList, List<GraphNodeAttribute> volumeNodeAttributes, List<GraphNodeAttribute> boundaryNodeAttributes, double st_Att, GlobalParameter globalParameter)
        {
            if (mergedGraph.DataCount != 0 & volumeNodeList.Count != 0 & volumeNodeAttributes.Count != 0 & st_Att != 0.0)
            {
                DataTree<Vector3d> attractionForce = new DataTree<Vector3d>();                      // list2 -> attractionForce

                List<double> volumeNodeCoinRadius = new List<double>();                                       // list -> coinRadius
                List<double> boundaryNodeCoinRadius = new List<double>();

                // 计算每个点由面积带来的半径
                for (int i = 0; i < volumeNodeAttributes.Count; i++)
                {
                    // 半径等于面积除以PI，然后开根号
                    volumeNodeCoinRadius.Add(Math.Sqrt(volumeNodeAttributes[i].NodeArea / Math.PI));
                }
                for (int i = 0; i < boundaryNodeAttributes.Count; i++)
                {
                    boundaryNodeCoinRadius.Add(Math.Sqrt(boundaryNodeAttributes[i].NodeArea / Math.PI));
                }

                // 对于每个代表volume的点（不含有NEWS点）
                for (int i = 0; i < mergedGraph.BranchCount; i++)
                {
                    List<int> branchItems = mergedGraph.Branch(i);
                    // 对于每一支的数据
                    for (int j = 0; j < branchItems.Count; j++)
                    {
                        // 当是adjacent连接关系时
                        if (branchItems[j] < 0)
                        {
                            // 对于adjacent连接关系为空值时（-1 - 4 = -5），其吸引力是0
                            if (branchItems[j] == -1 - globalParameter.BoundaryNodeCount)
                            {
                                attractionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                            }
                            // 对于其他adjacent情况
                            else
                            {
                                // 连线，算距离
                                Line line = new Line(volumeNodeList[i], boundaryNodeList[branchItems[j]+ globalParameter.BoundaryNodeCount]);
                                double distance = line.Length - (volumeNodeCoinRadius[i] + boundaryNodeCoinRadius[branchItems[j] + globalParameter.BoundaryNodeCount]);
                                // 如果相离（>）时，那么在吸引力的树形数据中，这个点的这一支上，添加一个NEWS点对于这个点的Vector3d
                                if (distance > 0)
                                {
                                    double vectorLength = st_Att * distance;
                                    Vector3d direction = line.Direction;
                                    direction.Unitize();
                                    attractionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                                }
                                else
                                {
                                    attractionForce.Add(new Vector3d(0.0,0.0,0.0), mergedGraph.Paths[i]);
                                }
                            }
                        }
                        // 当对于connectivity连接关系时
                        else
                        {
                            // 连线，算距离
                            Line line = new Line(volumeNodeList[i], volumeNodeList[branchItems[j]]);
                            double distance = line.Length - (volumeNodeCoinRadius[i] + volumeNodeCoinRadius[branchItems[j]]);
                            // 如果相离（>）时，那么在吸引力的树形数据中，这个点的这一支上，添加一个volume点对于这个点的Vector3d
                            if (distance > 0)
                            {
                                double vectorLength = st_Att * distance;
                                Vector3d direction = line.Direction;
                                direction.Unitize();
                                attractionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                            }
                            else
                            {
                                attractionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                            }
                        }
                    }
                }
                // 返回每个点（树形数据的一支branch）的volume点对它的吸引力，以及NEWS点对它的吸引力，两者加起来（共有mergeGraph.Branch[i].Count个）
                return attractionForce;
            }
            return null;
        }

        /// <summary>
        /// 计算每个代表体量volume的点node的排斥力列表（包括graph上的也包括NEWS上的）
        /// </summary>
        /// <param name="mergedGraph">每个代表体量的点node的connectivity（与体量点）和adjacent（与边界NEWS）</param>
        /// <param name="volumeNodeList">代表体量volume的点的列表</param>
        /// <param name="boundaryNodeList">代表边界NEWS的点的列表</param>
        /// <param name="volumeNodeAttributes">代表体量volume的点的属性</param>
        /// <param name="boundaryNodeAttributes">代表边界NEWS的点的属性</param>
        /// <param name="st_Rep"></param>
        /// <returns>返回每个点（树形数据的一支branch）的volume点对它的斥力（共有volumeNodeList.Count个），以及NEWS点对它的斥力（共有boundaryNodeList.Count个）</returns>
        public DataTree<Vector3d> RepulsionForce(DataTree<int> mergedGraph, List<Point3d> volumeNodeList, List<Point3d> boundaryNodeList, List<GraphNodeAttribute> volumeNodeAttributes, List<GraphNodeAttribute> boundaryNodeAttributes, double st_Rep)
        {
            if (mergedGraph.DataCount != 0 & volumeNodeList.Count != 0 & volumeNodeAttributes.Count != 0  & st_Rep != 0.0)
            {
                DataTree<Vector3d> repulsionForce = new DataTree<Vector3d>();

                List<double> volumeNodeCoinRadius = new List<double>();                                       // list -> coinRadius
                List<double> boundaryNodeCoinRadius = new List<double>();

                // 计算每个点由面积带来的半径
                for (int i = 0; i < volumeNodeAttributes.Count; i++)
                {
                    // 半径等于面积除以PI，然后开根号
                    volumeNodeCoinRadius.Add(Math.Sqrt(volumeNodeAttributes[i].NodeArea / Math.PI));
                }
                for (int i = 0; i < boundaryNodeAttributes.Count; i++)
                {
                    boundaryNodeCoinRadius.Add(Math.Sqrt(boundaryNodeAttributes[i].NodeArea / Math.PI));
                }

                // 对于每个代表volume的点，其他volume点对它的斥力
                for (int i = 0; i < volumeNodeList.Count; i++)
                {
                    for (int j = 0; j < volumeNodeList.Count; j++)
                    {
                        // 对于除去它自己的其他点
                        if (i == j)
                        {
                            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                            continue;
                        }
                        // 连线，算距离
                        Line line = new Line(volumeNodeList[i], volumeNodeList[j]);
                        double distance = line.Length - (volumeNodeCoinRadius[i] + volumeNodeCoinRadius[j]);
                        // 如果相交（<），或相切（=），那么在排斥力的树形数据中，这个点的这一支上，添加一个volume点对于这个点的Vector3d
                        if (distance <= 0)
                        {
                            double vectorLength = -st_Rep * Math.Pow(distance, 2);
                            Vector3d direction = line.Direction;
                            direction.Unitize();
                            repulsionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                        }
                        else
                        {
                            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                        }
                    }
                }
                // 对于每个代表volume的点，boundary点对它的斥力
                for (int i = 0; i < volumeNodeList.Count; i++)
                {
                    for (int j = 0; j < boundaryNodeList.Count; j++)
                    {
                        // 连线，算距离
                        Line line = new Line(volumeNodeList[i], boundaryNodeList[j]);
                        double distance = line.Length - (volumeNodeCoinRadius[i] + boundaryNodeCoinRadius[j]);
                        // 如果相交（<），或相切（=），那么在排斥力的树形数据中，这个点的这一支上，添加一个NEWS点对于这个点的Vector3d
                        if (distance <= 0)
                        {
                            double vectorLength = -st_Rep * Math.Pow(distance, 2);
                            Vector3d direction = line.Direction;
                            direction.Unitize();
                            repulsionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                        }
                        else
                        {
                            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                        }
                    }
                }

                //// 对于每个代表volume的点（不含有NEWS点）
                //for (int i = 0; i < mergedGraph.BranchCount; i++)
                //{
                //    for (int j = 0; j < mergedGraph.BranchCount; j++)
                //    {
                //        // 对于除去它自己的其他点
                //        if (i == j)
                //        {
                //            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                //        }
                //        // 连线，算距离
                //        Line line = new Line(volumeNodeList[i], volumeNodeList[j]);
                //        double distance = line.Length - (volumeNodeCoinRadius[i] + volumeNodeCoinRadius[j]);
                //        // 如果相交（<），或相切（=），那么在排斥力的树形数据中，这个点的这一支上，添加一个volume点对于这个点的Vector3d
                //        if (distance <= 0)
                //        {
                //            double vectorLength = -st_Rep * Math.Pow(distance, 2);
                //            Vector3d direction = line.Direction;
                //            direction.Unitize();
                //            repulsionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                //        }
                //        else
                //        {
                //            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                //        }
                //    }
                //}
                //// 对于每个代表volume的点
                //for (int i = 0; i < mergedGraph.BranchCount; i++)
                //{
                //    //对于每个NEWS点，来计算NEWS点对于它的排斥力
                //    for (int j = 0; j < boundaryNodeList.Count; j++)
                //    {
                //        // 连线，算距离
                //        Line line = new Line(volumeNodeList[i], boundaryNodeList[j]);
                //        double distance = line.Length - (volumeNodeCoinRadius[i] + boundaryNodeCoinRadius[j]);
                //        // 如果相交（<），或相切（=），那么在排斥力的树形数据中，这个点的这一支上，添加一个NEWS点对于这个点的Vector3d
                //        if (distance <= 0)
                //        {
                //            double vectorLength = -st_Rep * Math.Pow(distance, 2);
                //            Vector3d direction = line.Direction;
                //            direction.Unitize();
                //            repulsionForce.Add(vectorLength * direction, mergedGraph.Paths[i]);
                //        }
                //        else
                //        {
                //            repulsionForce.Add(new Vector3d(0.0, 0.0, 0.0), mergedGraph.Paths[i]);
                //        }
                //    }
                //}


                // 返回每个点（树形数据的一支branch）的volume点对它的斥力（共有volumeNodeList.Count个），以及NEWS点对它的斥力（共有boundaryNodeList.Count个）
                return repulsionForce;
            }
            return null;
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
            get { return new Guid("315a6da3-7af2-4c03-b21a-7effa3dfb3dd"); }
        }
    }
}