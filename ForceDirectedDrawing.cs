using Grasshopper;
using Grasshopper.GUI.Gradient;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;

namespace VolumeGeneratorBasedOnGraph
{
    public class ForceDirectedDrawing : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ForceDirectedDrawing class.
        /// </summary>
        public ForceDirectedDrawing()
          : base("ForceDirectedDrawing2", "ForceDirectedDrawing2",
              "Finds locations of vertices for a kissing-disk/coin-drawing of a graph. Note that you need to connect the output vertices to a general graph drawing component to actually draw the graph",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
            GraphEdges = new List<Line>();

            NodeBubbles = new List<Brep>();

            NodeMats = new List<DisplayMaterial>();

            NodeTexts = new List<TextDot>();
            NodeLocations = new List<Plane>();
            // NodeTextSize = new List<double>();

            ResultPoints = new List<Point3d>();
        }

        private bool ActiveFlag;

        private int Thickness;

        private List<Line> GraphEdges;

        private List<Brep> NodeBubbles;

        private List<DisplayMaterial> NodeMats;

        // private List<string> NodeTexts;
        private List<TextDot> NodeTexts;
        private List<Plane> NodeLocations;
        // private List<double> NodeTextSize;

        private List<Point3d> ResultPoints;



        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddIntegerParameter("Graph", "G", "A graph describing which node is connected to which other nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("Vertices", "GV", "Graph vertices (points representing nodes in space)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Attributes", "Att", "A list of lists containing name, area and color attributes of the disks", GH_ParamAccess.list);
            pManager.AddNumberParameter("AttractionStrength", "St_Atr", "Strength of attraction forces. Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 0.25);
            pManager.AddNumberParameter("RepulsionStrength", "St_Rep", "Strength of repulsion forces. Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 1.4);
            pManager.AddNumberParameter("ConvergenceTolearnce", "Tol", "How much error you would tolerate for considering a drawing converged/relaxed Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 0.001);
            pManager.AddIntegerParameter("MaximumIterations", "It_#", "Maximum number of iterations (recursions) you allow this algorithm to perform", GH_ParamAccess.item, 4000);
            pManager.AddBooleanParameter("Active", "Act", "Do you want this to perform any operation? If so, set it to true or false otherwise.", GH_ParamAccess.item, false);
            pManager[4].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("NewVolumeNode", "NGV", "A set of new vertices on which a neat coin/kissing disk drawing can be drawn", GH_ParamAccess.list);
            pManager.AddIntegerParameter("IterationCount", "Count", "How many iterations (recusrions) have happend?", GH_ParamAccess.item);
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

            GH_Structure<GH_Integer> gh_Structure_Graph = null;
            DataTree<int> graphTree = new DataTree<int>();

            List<Point3d> nodeList = new List<Point3d>();
            List<NodeAttribute> nodeAttributes = new List<NodeAttribute>();

            double st_Att = 0.0;
            double st_Rep = 0.0;
            double tolerance = 0.01;
            int iteration = 0;
            bool activeFlag = false;

            List<List<Vector3d>> attrationForce = new List<List<Vector3d>>();
            List<List<Vector3d>> repulsionForce = new List<List<Vector3d>>();

            if (DA.GetDataTree<GH_Integer>("Graph", out gh_Structure_Graph) 
                & DA.GetDataList<Point3d>("Vertices", nodeList)
                & DA.GetDataList<NodeAttribute>("Attributes", nodeAttributes))
            {
                // 将Graph从GH_Structure<GH_Integer>转化为DataTree<int>
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_Graph, ref graphTree);

                DA.GetData<double>("AttractionStrength", ref st_Att);
                DA.GetData<double>("RepulsionStrength", ref st_Rep);
                DA.GetData<double>("ConvergenceTolearnce", ref tolerance);
                DA.GetData<int>("MaximumIterations", ref iteration);
                DA.GetData<bool>("Active", ref activeFlag);

                ActiveFlag = activeFlag;

                if (ActiveFlag == false)
                {
                    ResultPoints.Clear();
                }

                // 主体计算部分
                List<Point3d> newVolumeNodeList = new List<Point3d>();
                bool repeatFlag = false;

                if (nodeList.Count == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "something went wrong! there are no vertices to operate on!");
                    return;
                }

                List<Point3d> currentIterationNodeList = nodeList;
                int currentIteration = 0;

                do
                {
                    if (activeFlag)
                    {
                        newVolumeNodeList = currentIterationNodeList;
                    }
                    DA.SetData("IterationCount", currentIteration);

                    attrationForce = AttractionForce(graphTree, nodeList, nodeAttributes, st_Att, globalParameter);
                    repulsionForce = RepulsionForce(graphTree, nodeList, nodeAttributes, st_Rep);

                    Result_Force(
                        nodeList, 
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
                    ResultPoints.AddRange(newVolumeNodeList);
                    
                }
                DA.SetDataList("NewVolumeNode", nodeList);
            }

            if (ResultPoints.Count != 0)
            {
                List<Line> edges = new List<Line>();
                for (int i = 0; i < graphTree.BranchCount; i++)
                {
                    for (int j = 0; j < graphTree.Branch(i).Count; j++)
                    {
                        edges.Add(new Line(ResultPoints[i], ResultPoints[graphTree.Branch(i)[j]]));
                    }
                }

                List<Brep> bubbles = new List<Brep>();
                for (int i = 0; i < nodeList.Count; i++)
                {
                    Plane centerPlane = Plane.WorldXY;
                    centerPlane.Origin = ResultPoints[i];
                    Circle circle = new Circle(centerPlane, Math.Sqrt(nodeAttributes[i].NodeArea / Math.PI));
                    bubbles.Add(Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(circle.ToNurbsCurve(), ResultPoints[i])));
                }

                List<DisplayMaterial> nodeMats = new List<DisplayMaterial>();
                List<double> radius = new List<double>();
                List<TextDot> texts = new List<TextDot>();
                List<double> remapRadius = new List<double>();

                for (int i = 0; i < nodeList.Count; i++)
                {
                    radius.Add(Math.Sqrt(nodeAttributes[i].NodeArea / Math.PI));
                    texts.Add(new TextDot(nodeAttributes[i].NodeLabel, ResultPoints[i]));
                }

                radius.Sort();
                radius.Reverse();
                for (int i = 0; i < nodeList.Count; i++)
                {
                    remapRadius.Add(ReMapNumber(Math.Sqrt(nodeAttributes[i].NodeArea / Math.PI), radius[0], radius[radius.Count - 1], 0, 1));
                }

                List<Color> colors = PickColor(remapRadius);
                foreach (Color color in colors)
                {
                    nodeMats.Add(new DisplayMaterial(color));
                }

                GraphEdges.Clear();
                NodeBubbles.Clear();
                NodeMats.Clear();
                NodeTexts.Clear();
                // NodeTextSize.Clear();

                Thickness = 2;
                GraphEdges.AddRange(edges);
                NodeBubbles.AddRange(bubbles);
                NodeMats.AddRange(nodeMats);
                NodeTexts.AddRange(texts);
            }
        }

        public void Result_Force(
            List<Point3d> nodeList,
            List<List<Vector3d>> attractionForce,
            List<List<Vector3d>> repulsionForce,
            double tolerance,
            ref List<Point3d> newNodeList,
            ref bool repeat)
        {
            if (!(nodeList.Count == 0
                | attractionForce.Count == 0
                | repulsionForce.Count == 0
                | tolerance == 0.0))
            {
                List<List<Vector3d>> forceList = new List<List<Vector3d>>();
                // List<Vector3d> list2 = new List<Vector3d>();
                List<Point3d> newPointLocations = new List<Point3d>();

                // 每个点的合力组成的列表
                List<Vector3d> resultantForceList = new List<Vector3d>();

                // 对每个点分别求合吸引力，合排斥力
                for (int i = 0; i < nodeList.Count; i++)
                {
                    //Vector3d allAttractionForce = new Vector3d(0, 0, 0);
                    //Vector3d allRepulsionForce = new Vector3d(0, 0, 0);
                    //for (int j = 0; j < attractionForce.Branches[i].Count; j++)
                    //{
                    //    allAttractionForce = Vector3d.Add(allAttractionForce, attractionForce.Branches[i][j]);
                    //}
                    //for (int j = 0; j < repulsionForce.Branches[i].Count; j++)
                    //{
                    //    allRepulsionForce = Vector3d.Add(allRepulsionForce, repulsionForce.Branches[i][j]);
                    //}

                    //resultantForceList.Add(Vector3d.Add(allAttractionForce, allRepulsionForce));

                    forceList.Add(new List<Vector3d>());
                    for (int j = 0; j < nodeList.Count; j++)
                    {
                        bool flag1 = attractionForce[i][j].IsValid;
                        bool flag2 = repulsionForce[i][j].IsValid;

                        if (!flag1 & !flag2)
                        {
                            Vector3d item = new Vector3d(0, 0, 0);
                            forceList[i].Add(item);
                        }
                        else
                        {
                            forceList[i].Add(Vector3d.Add(attractionForce[i][j], repulsionForce[i][j]));
                        }
                    }
                }

                double currentTolerance = 0.0;
                for (int i = 0; i < nodeList.Count; i++)
                {
                    Vector3d movement = default(Vector3d);
                    for (int j = 0; j < nodeList.Count; j++)
                    {
                        movement = forceList[i][j] + movement;
                    }

                    // list2.Add(vector3d);
                    newPointLocations.Add(nodeList[i] + movement);

                    currentTolerance += movement.Length;
                }

                // 
                List<Point3d> movedNodeList = new List<Point3d>();

                //for (int i = 0; i < nodeList.Count; i++)
                //{
                //    movedNodeList.Add(nodeList[i] + resultantForceList[i]);
                //    currentTolerance += resultantForceList[i].Length;
                //}
                //newNodeList = movedNodeList;
                //repeat = currentTolerance > tolerance;

                newNodeList = newPointLocations;
                repeat = currentTolerance > tolerance;
            }
        }

        /// <summary>
        /// 计算每个node的吸引力列表
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="nodeList"></param>
        /// <param name="nodeAttributes"></param>
        /// <param name="st_Att"></param>
        /// <param name="globalParameter"></param>
        /// <returns></returns>
        public List<List<Vector3d>> AttractionForce(
            DataTree<int> graph,
            List<Point3d> nodeList,
            List<NodeAttribute> nodeAttributes,
            double st_Att,
            GlobalParameter globalParameter)
        {
            if (graph.DataCount != 0
                & nodeList.Count != 0
                & nodeAttributes.Count != 0
                & st_Att != 0.0)
            {
                List<List<Vector3d>> attractionForce = new List<List<Vector3d>>();

                List<double> nodeCoinRadius = new List<double>();

                // 计算每个node由面积占比所带来的的半径
                for (int i = 0; i < nodeAttributes.Count; i++)
                {
                    // 半径等于面积除以PI，然后开根号
                    nodeCoinRadius.Add(Math.Sqrt(nodeAttributes[i].NodeArea / Math.PI));
                }

                // 对于每个点，它所受到的相邻的其他点给它的力就是，这一支上所有的元素
                for (int i = 0; i < graph.BranchCount; i++)
                {
                    attractionForce.Add(new List<Vector3d>());
                    for (int j = 0; j < graph.BranchCount; j++)
                    {
                        if (graph.Branch(i).Contains(j))
                        {
                            // 连线，算距离
                            Line line = new Line(nodeList[i], nodeList[j]);
                            double distance = line.Length - (nodeCoinRadius[i] + nodeCoinRadius[j]);

                            // 如果相离（>）时，那么在吸引力的树形数据中，这个点的这一支上，添加上一个相邻的点对它的吸引力Vecter3d
                            if (distance > 0)
                            {
                                double vectorLength = st_Att * distance;
                                Vector3d direction = line.Direction;
                                direction.Unitize();
                                attractionForce[i].Add(vectorLength * direction);
                            }
                            else
                            {
                                attractionForce[i].Add(new Vector3d(0.0, 0.0, 0.0));
                            }
                        }
                        else
                        {
                            attractionForce[i].Add(new Vector3d(0.0, 0.0, 0.0));
                        }
                    }
                }
                return attractionForce;
            }
            return null;
        }

        /// <summary>
        /// 计算每个node的排斥力列表
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="nodeList"></param>
        /// <param name="nodeAttributes"></param>
        /// <param name="st_Rep"></param>
        /// <returns></returns>
        public List<List<Vector3d>> RepulsionForce(
            DataTree<int> graph,
            List<Point3d> nodeList,
            List<NodeAttribute> nodeAttributes,
            double st_Rep)
        {
            if (graph.DataCount != 0
                & nodeList.Count != 0
                & nodeAttributes.Count != 0
                & st_Rep != 0.0)
            {
                List<List<Vector3d>> repulsionForce = new List<List<Vector3d>>();

                List<double> nodeCoinRadius = new List<double>();

                // 计算每个node由面积占比所带来的的半径
                for (int i = 0; i < graph.BranchCount; i++)
                {
                    nodeCoinRadius.Add(Math.Sqrt(nodeAttributes[i].NodeArea / Math.PI));
                }

                for (int i = 0; i < graph.BranchCount; i++)
                {
                    repulsionForce.Add(new List<Vector3d>());
                    for (int j = 0; j < graph.BranchCount; j++)
                    {

                        // 连线，算距离
                        Line line = new Line(nodeList[i], nodeList[j]);
                        double distance = line.Length - (nodeCoinRadius[i] + nodeCoinRadius[j]);
                        // 如果相交（<），或相切（=），那么在排斥力的树形数据中，这个点的这一支上，添加上一个相邻的点对它的排斥力Vecter3d
                        if (distance <= 0.0)
                        {
                            Vector3d direction = line.Direction;
                            direction.Unitize();
                            repulsionForce[i].Add(-st_Rep * Math.Pow(distance, 2.0) * direction);
                        }
                        else
                        {
                            repulsionForce[i].Add(new Vector3d(0.0, 0.0, 0.0));
                        }
                    }
                }
                return repulsionForce;
            }
            return null;
        }

        public double ReMapNumber(double x, double iMin, double iMax, double oMin, double oMax)
        {
            return (iMax - iMin) / (oMax - oMin) * (x - oMin) + iMin;
        }


        public List<Color> PickColor(List<double> x)
        {
            GH_Gradient gh_Gradient = GH_Gradient.Heat();
            gh_Gradient.Linear = true;
            List<Color> colorList = new List<Color>();
            for (int i = 0; i < x.Count; i++)
            {
                colorList.Add(gh_Gradient.ColourAt(x[i]));
            }
            return colorList;
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);
            
            if (ActiveFlag == true)
            {
                args.Display.DrawLines(GraphEdges, Color.Gray, Thickness);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportMeshes(args);


            if (ActiveFlag == true)
            {
                for (int i = 0; i < NodeBubbles.Count; i++)
                {
                    
                    args.Display.DrawBrepShaded(NodeBubbles[i], NodeMats[i]);

                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(NodeTexts[i], Color.Black, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);

                }
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
            get { return new Guid("eb1599d6-cfc1-4421-a3e4-f764bc13ce2f"); }
        }
    }
}