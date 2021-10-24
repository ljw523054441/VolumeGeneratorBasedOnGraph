using Grasshopper;
using Grasshopper.GUI.Gradient;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

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
            ShowForceDiagram = true;
            
            GraphEdges = new List<Line>();

            NodeBubbles = new List<Brep>();

            NodeMats = new List<DisplayMaterial>();

            NodeTexts = new List<TextDot>();
            NodeLocations = new List<Plane>();
            // NodeTextSize = new List<double>();

            ResultPoints = new List<Point3d>();

            base.Message = "开发中";
        }

        private bool ShowForceDiagram;
        
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

            pManager.AddCurveParameter("BoundingCurve", "BC", "A polyline confining the area of the diagram", GH_ParamAccess.item);

            pManager.AddIntegerParameter("Graph", "G", "A graph describing which node is connected to which other nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("Vertices", "GV", "Graph vertices (points representing nodes in space)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Attributes", "Att", "A list of lists containing name, area and color attributes of the disks", GH_ParamAccess.list);
            pManager.AddNumberParameter("AttractionStrength", "St_Atr", "Strength of attraction forces. Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 0.25);
            pManager.AddNumberParameter("RepulsionStrength", "St_Rep", "Strength of repulsion forces. Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 1.4);
            pManager.AddNumberParameter("ConvergenceTolearnce", "Tol", "How much error you would tolerate for considering a drawing converged/relaxed Do not change this unless you exactly know what you are doing!", GH_ParamAccess.item, 0.001);
            pManager.AddIntegerParameter("MaximumIterations", "It_#", "Maximum number of iterations (recursions) you allow this algorithm to perform", GH_ParamAccess.item, 4000);
            pManager.AddBooleanParameter("Active", "Act", "Do you want this to perform any operation? If so, set it to true or false otherwise.", GH_ParamAccess.item, false);
            pManager[5].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
            pManager[8].Optional = true;
            pManager.AddNumberParameter("Scale", "S", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("NewVolumeNode", "NGV", "A set of new vertices on which a neat coin/kissing disk drawing can be drawn", GH_ParamAccess.list);
            pManager.AddIntegerParameter("IterationCount", "Count", "How many iterations (recusrions) have happend?", GH_ParamAccess.item);
            pManager.AddCurveParameter("VoronoiCells", "VC", "Bubbles in the form of Voronoi cells", GH_ParamAccess.list);
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

            // 场地边界
            Curve boundingCurve = null;
            DA.GetData<Curve>("BoundingCurve", ref boundingCurve);
            if (boundingCurve == null || !boundingCurve.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the bounding curve is not closed!");
                return;
            }
            if (!boundingCurve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the bounding curve is not planar!");
                return;
            }
            // 场地总面积groundArea
            double groundArea = AreaMassProperties.Compute(boundingCurve).Area;


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
                    return;
                }


                List<double> nodeAreas = new List<double>();
                for (int i = 0; i < nodeAttributes.Count; i++)
                {
                    nodeAreas.Add(nodeAttributes[i].NodeAreaProportion * groundArea);
                }

                List<Point3d> newVolumeNodeList = new List<Point3d>();
                if (nodeList.Count == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "something went wrong! there are no vertices to operate on!");
                }
                else
                {
                    // 主体计算部分
                    List<Point3d> currentIterationNodeList = nodeList;
                    int currentIteration = 0;
                    bool repeatFlag = false;

                    do
                    {
                        if (activeFlag)
                        {
                            newVolumeNodeList = currentIterationNodeList;
                        }
                        DA.SetData("IterationCount", currentIteration);

                        attrationForce = AttractionForce(graphTree, currentIterationNodeList, nodeAreas, st_Att);
                        repulsionForce = RepulsionForce(graphTree, currentIterationNodeList, nodeAreas, st_Rep);

                        Result_Force(
                            currentIterationNodeList,
                            attrationForce,
                            repulsionForce,
                            tolerance,
                            ref currentIterationNodeList,
                            ref repeatFlag);
                        currentIteration++;
                    }
                    while (!(!repeatFlag || !activeFlag || currentIteration > iteration));

                    if (activeFlag)
                    {
                        DA.SetDataList("NewVolumeNode", newVolumeNodeList);
                        ResultPoints.AddRange(newVolumeNodeList);

                    }
                    else
                    {
                        DA.SetDataList("NewVolumeNode", nodeList);
                        ResultPoints.AddRange(nodeList);
                    }
                }
            }

            CalDrawParameters(graphTree, ResultPoints, nodeAttributes, groundArea);



            // 计算加权Voronoi
            double scale = 0.0;
            DA.GetData("Scale", ref scale);

            List<double> radius = new List<double>();

            for (int i = 0; i < nodeAttributes.Count; i++)
            {
                radius.Add(Math.Sqrt(nodeAttributes[i].NodeAreaProportion * groundArea / Math.PI) * scale);
            }

            List<Curve> cellCurve = AlphaVoronoi(Plane.WorldXY, ResultPoints, radius, boundingCurve);

            DA.SetDataList("VoronoiCells", cellCurve);
            // DA.SetDataTree(0, adjacencyGraph);



        }

        private void CalDrawParameters(DataTree<int> graphTree, List<Point3d> nodeList, List<NodeAttribute> nodeAttributes, double groundArea)
        {
            // 当运行产生结果时，绘制图像
            if (ResultPoints.Count != 0)
            {
                // edge
                List<Line> edges = new List<Line>();
                for (int i = 0; i < graphTree.BranchCount; i++)
                {
                    for (int j = 0; j < graphTree.Branch(i).Count; j++)
                    {
                        edges.Add(new Line(ResultPoints[i], ResultPoints[graphTree.Branch(i)[j]]));
                    }
                }

                // bubble
                List<Brep> bubbles = new List<Brep>();
                for (int i = 0; i < nodeList.Count; i++)
                {
                    Plane centerPlane = Plane.WorldXY;
                    centerPlane.Origin = ResultPoints[i];
                    Circle circle = new Circle(centerPlane, Math.Sqrt(nodeAttributes[i].NodeAreaProportion * groundArea / Math.PI));
                    bubbles.Add(Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(circle.ToNurbsCurve(), ResultPoints[i])));
                }

                // mat and textdot
                List<DisplayMaterial> nodeMats = new List<DisplayMaterial>();
                List<double> radius = new List<double>();
                List<TextDot> texts = new List<TextDot>();
                List<double> remapRadius = new List<double>();

                for (int i = 0; i < nodeList.Count; i++)
                {
                    radius.Add(nodeAttributes[i].NodeAreaProportion);
                    texts.Add(new TextDot(nodeAttributes[i].NodeLabel, ResultPoints[i]));
                }

                double iMax = radius.Max();
                double iMin = radius.Min();


                for (int i = 0; i < nodeList.Count; i++)
                {
                    remapRadius.Add(ReMapNumber(radius[i], iMin, iMax, 1, 0));
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

        /// <summary>
        /// 求每个node受到的吸引力和排斥力的合力
        /// </summary>
        /// <param name="nodeList"></param>
        /// <param name="attractionForce"></param>
        /// <param name="repulsionForce"></param>
        /// <param name="tolerance"></param>
        /// <param name="newNodeList"></param>
        /// <param name="repeat"></param>
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
        /// <param name="nodeAreas"></param>
        /// <param name="st_Att"></param>
        /// <param name="globalParameter"></param>
        /// <returns></returns>
        public List<List<Vector3d>> AttractionForce(
            DataTree<int> graph,
            List<Point3d> nodeList,
            List<double> nodeAreas,
            double st_Att)
        {
            if (graph.DataCount != 0
                & nodeList.Count != 0
                & nodeAreas.Count != 0
                & st_Att != 0.0)
            {
                List<List<Vector3d>> attractionForce = new List<List<Vector3d>>();

                List<double> nodeCoinRadius = new List<double>();

                // 计算每个node由面积占比所带来的的半径
                for (int i = 0; i < nodeAreas.Count; i++)
                {
                    // 半径等于面积除以PI，然后开根号
                    nodeCoinRadius.Add(Math.Sqrt(nodeAreas[i] / Math.PI));
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
        /// <param name="nodeAreas"></param>
        /// <param name="st_Rep"></param>
        /// <returns></returns>
        public List<List<Vector3d>> RepulsionForce(
            DataTree<int> graph,
            List<Point3d> nodeList,
            List<double> nodeAreas,
            double st_Rep)
        {
            if (graph.DataCount != 0
                & nodeList.Count != 0
                & nodeAreas.Count != 0
                & st_Rep != 0.0)
            {
                List<List<Vector3d>> repulsionForce = new List<List<Vector3d>>();

                List<double> nodeCoinRadius = new List<double>();

                // 计算每个node由面积占比所带来的的半径
                for (int i = 0; i < graph.BranchCount; i++)
                {
                    nodeCoinRadius.Add(Math.Sqrt(nodeAreas[i] / Math.PI));
                }

                for (int i = 0; i < graph.BranchCount; i++)
                {
                    repulsionForce.Add(new List<Vector3d>());
                    for (int j = 0; j < graph.BranchCount; j++)
                    {
                        if (i != j)
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



        /// <summary>
        /// 计算加权Voronoi
        /// </summary>
        /// <param name="basePlane"></param>
        /// <param name="points"></param>
        /// <param name="radius"></param>
        /// <param name="boundingCurve"></param>
        /// <param name="adjacencyGraph"></param>
        /// <returns></returns>
        public List<Curve> AlphaVoronoi(
            Plane basePlane,
            List<Point3d> points,
            List<double> radius,
            Curve boundingCurve
            )
        {
            List<Curve> result;

            if (!basePlane.IsValid || points == null || radius == null)
            {
                result = null;
            }
            else
            {
                // 如果radius中有小于0的数，跳出
                if (radius.FindAll(x => x <= 0.0).Count > 0)
                {
                    result = null;
                }
                else
                {
                    if (points.Count != radius.Count)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the number of attribute and vertices do not match!");
                        result = null;
                    }
                    else
                    {
                        // Vertices.Clear();
                        // Indices.Clear();

                        List<Curve> pointBubbleCurves = new List<Curve>();
                        List<List<Curve>> pointBubbleCurvesLoL = new List<List<Curve>>();
                        // DataTree<int> adjacencyDataTree = new DataTree<int>();
                        List<List<bool>> dominatedFlagaLoL = new List<List<bool>>();
                        List<Curve> cellCurveList = new List<Curve>();
                        List<Line> verticalLineList = new List<Line>();

                        BoundingBox pointsBoundingBox = new BoundingBox(points);

                        Plane worldXY = Plane.WorldXY;
                        for (int i = 0; i < points.Count; i++)
                        {
                            worldXY.Origin = points[i];
                            // List<Curve> list6 = list;
                            Circle bubble = new Circle(worldXY, radius[i]);
                            // Circle circle2 = circle;
                            pointBubbleCurves.Add(bubble.ToNurbsCurve());

                            pointsBoundingBox = BoundingBox.Union(pointBubbleCurves[i].GetBoundingBox(false), pointsBoundingBox);

                            // Indices.Add(i);

                            // adjacencyDataTree.EnsurePath(i);

                            pointBubbleCurvesLoL.Add(new List<Curve>());
                            pointBubbleCurvesLoL[i].Add(pointBubbleCurves[i]);

                            dominatedFlagaLoL.Add(new List<bool>());
                        }

                        // pointsBoundingBox上x最大y最大的点，到，x最小y最小的点的矩形
                        Rectangle3d pointsPlanarBoundingBox = new Rectangle3d(basePlane, pointsBoundingBox.Corner(true, true, true), pointsBoundingBox.Corner(false, false, true));
                        // Rectangle3d rectangle3d2 = rectangle3d;
                        Curve pointsPlanarBoundingBoxCurve = pointsPlanarBoundingBox.ToNurbsCurve();

                        // 
                        if (boundingCurve == null)
                        {
                            boundingCurve = pointsPlanarBoundingBoxCurve;
                        }
                        else
                        {
                            if (boundingCurve != null & !boundingCurve.IsClosed)
                            {
                                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the bounding curve is not closed!");
                                return null;
                            }
                        }


                        // 
                        for (int i = 0; i < points.Count; i++)
                        {
                            for (int j = 0; j < points.Count; j++)
                            {
                                // 两点间距离
                                double distance = points[j].DistanceTo(points[i]);
                                // 如果两点间距离小于两个半径的和，即两圆相交，并且i!=j时
                                if (distance < (radius[j] + radius[i]) & i != j)
                                {
                                    // 
                                    // adjacencyDataTree.Branch(i).Add(j);
                                    bool dominated = false;
                                    bool outside = false;
                                    double distanceBetweenP1andVertical = 0;

                                    // 求短的verticalLines
                                    Line verticalLine = DrawVerticalLineBetweenTwoCells(
                                        points[j],
                                        points[i],
                                        radius[j],
                                        radius[i],
                                        basePlane,
                                        ref distanceBetweenP1andVertical,
                                        boundingCurve,
                                        ref dominated,
                                        ref outside);
                                    verticalLineList.Add(verticalLine);

                                    if (outside)
                                    {
                                        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "The boundary curve is so close to the points that at least one of the cell boundaries is outside of it. try enlarging the boundary or changing the weights such that the power boundary lies inside.");
                                        return null;
                                    }

                                    pointBubbleCurvesLoL[i].Add(
                                        VoronoiRegionPart(
                                        points[i],
                                        points[j],
                                        boundingCurve,
                                        verticalLine,
                                        0.01,
                                        distanceBetweenP1andVertical,
                                        dominated));
                                    dominatedFlagaLoL[i].Add(dominated);
                                }
                            }


                            if (!dominatedFlagaLoL[i].Contains(true) & pointBubbleCurvesLoL[i] != null)
                            {
                                cellCurveList.Add(UnionIntersectCells(pointBubbleCurvesLoL[i], boundingCurve, basePlane, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance));
                            }
                            else
                            {
                                cellCurveList.Add(null);
                            }
                        }


                        // adjacencyGraph = adjacencyDataTree;
                        result = cellCurveList;
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// 绘制两点之间连线的权重垂直线
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="radius1"></param>
        /// <param name="radius2"></param>
        /// <param name="basePlane"></param>
        /// <param name="distanceBetweenP1andVertical"></param>
        /// <param name="boundingCurve"></param>
        /// <param name="dominated"></param>
        /// <param name="outside"></param>
        /// <returns></returns>
        public Line DrawVerticalLineBetweenTwoCells(
            Point3d point1,
            Point3d point2,
            double radius1,
            double radius2,
            Plane basePlane,
            ref double distanceBetweenP1andVertical,
            Curve boundingCurve,
            ref bool dominated,
            ref bool outside)
        {
            //两点之间的间距
            double distance = point1.DistanceTo(point2);
            // 点1的权重和点2的权重
            double x1 = distance / 2.0 + (Math.Pow(radius1, 2.0) - Math.Pow(radius2, 2.0)) / (2.0 * distance);
            double x2 = distance / 2.0 - (Math.Pow(radius1, 2.0) - Math.Pow(radius2, 2.0)) / (2.0 * distance);

            Line line = new Line(point1, point2);
            distanceBetweenP1andVertical = x1;
            // 权重垂直线的端点
            Point3d pointBetweenTwoCell = line.PointAt(distanceBetweenP1andVertical / distance);

            // 如果这个端点在boundingCurve外或者上，那么把outside设为true
            bool flag = boundingCurve.Contains(pointBetweenTwoCell, basePlane, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) != PointContainment.Inside
                || boundingCurve.Contains(pointBetweenTwoCell, basePlane, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == PointContainment.Coincident;
            if (flag)
            {
                outside = true;
            }
            // 如果
            flag = distanceBetweenP1andVertical > Math.Max(radius1, radius2);
            if (flag)
            {
                dominated = true;
            }

            // Point3d point3d2 = pointBetweenTwoCell;
            Vector3d lineDirection = new Vector3d(point2 - point1);
            // 叉积求出垂直方向，在这个方向上画长度为1的直线
            Line verticalLines = new Line(pointBetweenTwoCell, Vector3d.CrossProduct(lineDirection, basePlane.ZAxis), 1.0);
            return verticalLines;
        }

        /// <summary>
        /// 求cell
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="boundingCurve"></param>
        /// <param name="verticalLine"></param>
        /// <param name="tolerance"></param>
        /// <param name="distanceBetweenP1andVertical"></param>
        /// <param name="dominated"></param>
        /// <returns></returns>
        public Curve VoronoiRegionPart(
            Point3d point1,
            Point3d point2,
            Curve boundingCurve,
            Line verticalLine,
            double tolerance,
            double distanceBetweenP1andVertical,
            bool dominated)
        {
            Curve voronoiCellCurve;

            if (dominated)
            {
                voronoiCellCurve = null;
            }
            else
            {
                // 对于刚才的长度为1的verticalCurve进行延伸
                Curve verticalCurve = verticalLine.ToNurbsCurve();
                Curve extendedVerticalCurve = verticalCurve.Extend(CurveEnd.Both, CurveExtensionStyle.Line, new GeometryBase[]
                {
                    boundingCurve
                });

                if (extendedVerticalCurve == null)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "The boundary curve is so close to the points that at least one of the cell boundaries is outside of it.");
                }

                // 找到boundingCurve上的交点，一定有两个
                CurveIntersections curveIntersections = Intersection.CurveCurve(extendedVerticalCurve, boundingCurve, tolerance, 0.0);

                // 找到交点
                List<Point3d> intersectPointOnCurveB = new List<Point3d>();
                List<double> intersectParamOnCurveB = new List<double>();
                for (int i = 0; i < curveIntersections.Count; i++)
                {
                    IntersectionEvent intersectionEvent = curveIntersections[i];
                    intersectPointOnCurveB.Add(intersectionEvent.PointB);
                    intersectParamOnCurveB.Add(intersectionEvent.ParameterB);
                }

                // 在交点处Trim
                Curve halfBoundingCurve01 = boundingCurve.Trim(intersectParamOnCurveB[0], intersectParamOnCurveB[1]);
                Curve halfBoundingCurve10 = boundingCurve.Trim(intersectParamOnCurveB[1], intersectParamOnCurveB[0]);

                // 得到线
                Line line01 = new Line(halfBoundingCurve01.PointAtStart, halfBoundingCurve01.PointAtEnd);
                Curve curve01 = line01.ToNurbsCurve();
                Line line10 = new Line(halfBoundingCurve10.PointAtStart, halfBoundingCurve10.PointAtEnd);
                Curve curve10 = line10.ToNurbsCurve();

                // 合并线，形成闭合polyline
                Curve splittedBoundingCurve1 = Curve.JoinCurves(new Curve[]
                {
                    halfBoundingCurve01,
                    curve01
                })[0];
                Curve splittedBoundingCurve2 = Curve.JoinCurves(new Curve[]
                {
                    halfBoundingCurve10,
                    curve10
                })[0];

                Polyline splittedBoundingPolyline1;
                splittedBoundingCurve1.TryGetPolyline(out splittedBoundingPolyline1);
                Polyline splittedBoundingPolyline2;
                splittedBoundingCurve2.TryGetPolyline(out splittedBoundingPolyline2);

                // 两个子区域各自的中心点
                Point3d center1 = splittedBoundingPolyline1.CenterPoint();
                Point3d center2 = splittedBoundingPolyline2.CenterPoint();
                Vector3d vector = new Vector3d(center2 - center1);

                Curve curve9 = null;

                Line line = new Line(point1, point2);
                // double length = line.Length;

                double dotProduct = vector * line.Direction;
                if (dotProduct > 0.0)
                {
                    curve9 = splittedBoundingCurve1;
                }
                else
                {
                    if (dotProduct < 0.0)
                    {
                        curve9 = splittedBoundingCurve2;
                    }
                }
                voronoiCellCurve = curve9;
            }
            return voronoiCellCurve;
        }

        public Curve UnionIntersectCells(List<Curve> cells, Curve boundingCurve, Plane basePlane, double tolerance)
        {
            Curve curve;
            if (cells.Count == 0)
            {
                curve = null;
            }
            else
            {
                curve = cells[1];
                for (int i = 0; i < cells.Count; i++)
                {
                    if (Curve.PlanarClosedCurveRelationship(curve, cells[i], basePlane, tolerance) == RegionContainment.Disjoint)
                    {
                        break;
                    }
                    curve = Curve.CreateBooleanIntersection(curve, cells[i], tolerance)[0];
                }

                Curve[] array = Curve.CreateBooleanIntersection(curve, cells[0], tolerance);
                if (array.Length == 0)
                {
                    return null;
                }
                curve = array[0];
            }
            return curve;
        }


        public double ReMapNumber(double x, double targetMin, double targetMax, double sourceMin, double sourceMax)
        {
            //return (iMax - iMin) / (oMax - oMin) * (x - oMin) + iMin;

            return (x - targetMin) / (targetMax - targetMin) * (sourceMax - sourceMin) + sourceMin;
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