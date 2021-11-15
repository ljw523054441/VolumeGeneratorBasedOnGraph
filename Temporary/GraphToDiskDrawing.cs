using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Drawing;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class GraphToDiskDrawing : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GraphToDiskDrawing class.
        /// </summary>
        public GraphToDiskDrawing()
          : base("GraphToDiskDrawing", "GraphDiskDrawing",
              "展示图结构",
              "VolumeGeneratorBasedOnGraph", "GraphDrawing")
        {
            InnerEdges = new List<Line>();
            OuterEdges = new List<Line>();
            InnerBubbles = new List<Brep>();
            OuterBubbles = new List<Brep>();
            Mats = new List<DisplayMaterial>();
            Mats2 = new List<DisplayMaterial>();

            Vertices = new PointCloud();

            InnerTexts = new List<string>();
            OuterTexts = new List<string>();
            InnerLocations = new List<Plane>();
            OuterLocations = new List<Plane>();
            InnerTextSize = new List<double>();
            OuterTextSize = new List<double>();
        }

        private int Thickness;

        private List<Line> InnerEdges;
        private List<Line> OuterEdges;

        private List<Brep> InnerBubbles;
        private List<Brep> OuterBubbles;

        private List<DisplayMaterial> Mats;
        private List<DisplayMaterial> Mats2;

        private PointCloud Vertices;

        private List<string> InnerTexts;
        private List<string> OuterTexts;
        private List<Plane> InnerLocations;
        private List<Plane> OuterLocations;
        private List<double> InnerTextSize;
        private List<double> OuterTextSize;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddPlaneParameter("Plane", "P", "A plane for placing the disk-graph drawing", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddIntegerParameter("Graph", "G", "A graph describing which node is connected to which other nodes", GH_ParamAccess.tree);
            pManager.AddPointParameter("Vertices", "GV", "Graph vertices (points representing nodes in space)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Attributes", "Att", "A list of lists containing name, area and color attributes of the disks", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Thickness", "Th", "Thickness of graph links in the diagram", GH_ParamAccess.item, 1);
            pManager[1].Optional = true;
            pManager[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
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


            Plane location = default(Plane);
            int thickness = 0;                                                    // num -> thickness

            GH_Structure<GH_Integer> gh_Structue_Graph = new GH_Structure<GH_Integer>();

            List<Point3d> nodePoints = new List<Point3d>();                       // list -> 
            List<GraphNodeAttribute> nodeAttributes = new List<GraphNodeAttribute>();          // list2 ->
            List<string> nodeLabels = new List<string>();                    // list3 -> nodeLabels
            List<double> nodeAreas = new List<double>();                    // list4 -> nodeAreas
            // List<Color> nodeMatColor = new List<Color>();                      // list5 -> nodeMatColor

            List<Point3d> innerNodePoints = new List<Point3d>();                       // list -> 
            List<GraphNodeAttribute> innerNodeAttributes = new List<GraphNodeAttribute>();          // list2 ->
            List<string> innerNodeLabels = new List<string>();                    // list3 -> nodeLabels
            List<double> innerNodeAreas = new List<double>();                    // list4 -> nodeAreas
            // List<Color> innernodeMatColor = new List<Color>();                      // list5 -> nodeMatColor

            List<Point3d> outerNodePoints = new List<Point3d>();                       // list -> 
            List<GraphNodeAttribute> outerNodeAttributes = new List<GraphNodeAttribute>();          // list2 ->
            List<string> outerNodeLabels = new List<string>();                    // list3 -> nodeLabels
            List<double> outerNodeAreas = new List<double>();                    // list4 -> nodeAreas
            // List<Color> outerNodeMatColor = new List<Color>();                      // list5 -> nodeMatColor

            if (DA.GetDataTree<GH_Integer>("Graph", out gh_Structue_Graph) 
                & DA.GetDataList<Point3d>("Vertices", nodePoints) 
                & DA.GetDataList<GraphNodeAttribute>("Attributes", nodeAttributes))
            {
                DA.GetData<Plane>("Plane", ref location);
                DA.GetData<int>("Thickness", ref thickness);

                foreach (GraphNodeAttribute attribute in nodeAttributes)
                {
                    nodeLabels.Add(attribute.NodeLabel);
                    nodeAreas.Add(attribute.NodeArea);
                }

                // List<List<int>> list6 = new List<List<int>>();              // list6 -> graphDataTree;
                DataTree<int> graphDataTree = new DataTree<int>();
                DataTree<int> innerGraph;
                DataTree<int> outerGraph;

                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structue_Graph, ref graphDataTree);
                UtilityFunctions.GraphToSubGraph(graphDataTree, globalParameter, out innerGraph, out outerGraph);

                List<Point3d> nodePointsRelocated = UtilityFunctions.Relocate(nodePoints, Plane.WorldXY, location);    // list8 -> nodePointsRelocated
                
                for (int i = 0; i < nodePointsRelocated.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        innerNodePoints.Add(nodePointsRelocated[i]);
                        innerNodeAttributes.Add(nodeAttributes[i]);
                        innerNodeLabels.Add(nodeLabels[i]);
                        innerNodeAreas.Add(nodeAreas[i]);
                    }
                    else
                    {
                        outerNodePoints.Add(nodePointsRelocated[i]);
                        outerNodeAttributes.Add(nodeAttributes[i]);
                        outerNodeLabels.Add(nodeLabels[i]);
                        outerNodeAreas.Add(nodeAreas[i]);
                    }
                }

                // List<List<Line>> list9 = new List<List<Line>>();         // list9 -> graphLineDataTree
                // DataTree<Line> graphLineDataTree = new DataTree<Line>();
                List<Line> innerLineList = new List<Line>();                       // list10 -> graphLineList
                List<Line> outerLineList = new List<Line>();
                for (int i = 0; i < innerGraph.BranchCount; i++)
                {
                    for (int j = 0; j < innerGraph.Branch(i).Count; j++)
                    {
                        Line line = new Line(nodePointsRelocated[i], nodePointsRelocated[innerGraph.Branch(i)[j]]);
                        // graphLineDataTree.Branch(i).Add(line);
                        innerLineList.Add(line);
                    }
                }
                for (int i = 0; i < outerGraph.BranchCount; i++)
                {
                    for (int j = 0; j < outerGraph.Branch(i).Count; j++)
                    {
                        Line line = new Line(nodePointsRelocated[i + volumeNodeCount], nodePointsRelocated[outerGraph.Branch(i)[j]]);
                        // graphLineDataTree.Branch(i).Add(line);
                        outerLineList.Add(line);
                    }
                }


                List<double> innerNodeRadiusListForLabelSize = new List<double>();           // list11 -> nodeRadiusListForLabelSize
                List<double> outerNodeRadiusListForLabelSize = new List<double>();
                for (int i = 0; i < nodePointsRelocated.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        innerNodeRadiusListForLabelSize.Add(Math.Sqrt(0.2 * nodeAreas[i]) / Math.PI);
                    }
                    else
                    {
                        outerNodeRadiusListForLabelSize.Add(Math.Sqrt(0.2 * nodeAreas[i]) / Math.PI);
                    }
                }


                List<Plane> planesBasedInnerNodePoints = new List<Plane>();                                 // list12 -> planesBasedNodePoints
                List<Plane> planesBasedOuterNodePoints = new List<Plane>();
                for (int i = 0; i < nodePointsRelocated.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        Plane plane = new Plane(nodePointsRelocated[i], location.XAxis, location.YAxis);
                        planesBasedInnerNodePoints.Add(plane);
                    }
                    else
                    {
                        Plane plane = new Plane(nodePointsRelocated[i], location.XAxis, location.YAxis);
                        planesBasedOuterNodePoints.Add(plane);
                    }
                }


                List<double> innerNodeRadius = new List<double>();                           // list14 -> nodeRadius
                List<Curve> innerNodeCircleList = new List<Curve>();                         // list15 -> nodeCircleList
                List<Brep> innerNodeCircleBrepList = new List<Brep>();                       // list16 -> nodeCircleBrepList
                // List<DisplayMaterial> list17 = new List<DisplayMaterial>();
                List<double> outerNodeRadius = new List<double>();
                List<Curve> outerNodeCircleList = new List<Curve>();
                List<Brep> outerNodeCircleBrepList = new List<Brep>();
                for (int i = 0; i < nodePointsRelocated.Count; i++)
                {
                    if (i < volumeNodeCount)
                    {
                        innerNodeRadius.Add(Math.Sqrt(nodeAreas[i] / Math.PI));

                        Plane worldXY = Plane.WorldXY;
                        worldXY.Origin = nodePointsRelocated[i];

                        Circle circle = new Circle(worldXY, innerNodeRadius[i]);
                        innerNodeCircleList.Add(circle.ToNurbsCurve());
                        Brep brep = Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(innerNodeCircleList[i], nodePointsRelocated[i]));
                        brep.Flip();
                        innerNodeCircleBrepList.Add(brep);
                    }
                    else
                    {
                        outerNodeRadius.Add(Math.Sqrt(nodeAreas[i] / Math.PI));

                        Plane worldXY = Plane.WorldXY;
                        worldXY.Origin = nodePointsRelocated[i];

                        Circle circle = new Circle(worldXY, outerNodeRadius[i - volumeNodeCount]);
                        outerNodeCircleList.Add(circle.ToNurbsCurve());
                        Brep brep = Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(outerNodeCircleList[i - volumeNodeCount], nodePointsRelocated[i]));
                        brep.Flip();
                        outerNodeCircleBrepList.Add(brep);
                    }
                }



                Thickness = 0;

                InnerEdges.Clear();
                OuterEdges.Clear();
                InnerBubbles.Clear();
                OuterBubbles.Clear();

                InnerLocations.Clear();
                OuterLocations.Clear();
                InnerTextSize.Clear();
                OuterTextSize.Clear();
                InnerTexts.Clear();
                OuterTexts.Clear();


                Thickness = thickness;

                InnerEdges.AddRange(innerLineList);
                OuterEdges.AddRange(outerLineList);
                InnerBubbles.AddRange(innerNodeCircleBrepList);
                OuterBubbles.AddRange(outerNodeCircleBrepList);

                InnerLocations.AddRange(planesBasedInnerNodePoints);
                OuterLocations.AddRange(planesBasedOuterNodePoints);
                InnerTextSize.AddRange(innerNodeRadiusListForLabelSize);
                OuterTextSize.AddRange(outerNodeRadiusListForLabelSize);
                InnerTexts.AddRange(innerNodeLabels);
                OuterTexts.AddRange(outerNodeLabels);
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            args.Display.DrawLines(InnerEdges, Color.Black, Thickness);
            args.Display.DrawLines(OuterEdges, Color.Gray, Thickness);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            if (InnerBubbles.Count != 0)
            {
                for (int i = 0; i < InnerBubbles.Count; i++)
                {
                    // DisplayMaterial displayMaterial = Mats[i];
                    args.Display.DrawBrepShaded(InnerBubbles[i], new DisplayMaterial(Color.Red));
                    args.Display.Draw3dText(InnerTexts[i], Color.White, InnerLocations[i], InnerTextSize[i], "Arial",
                        true, false,
                        Rhino.DocObjects.TextHorizontalAlignment.Center,
                        Rhino.DocObjects.TextVerticalAlignment.Middle);
                }
            }
            if (OuterBubbles.Count != 0)
            {
                for (int i = 0; i < OuterBubbles.Count; i++)
                {
                    args.Display.DrawBrepWires(OuterBubbles[i], Color.Blue);
                    args.Display.Draw3dText(OuterTexts[i], Color.White, OuterLocations[i], OuterTextSize[i], "Arial",
                        false, false,
                        Rhino.DocObjects.TextHorizontalAlignment.Center,
                        Rhino.DocObjects.TextVerticalAlignment.Middle);
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
            get { return new Guid("5c1f85b2-b769-4e54-9729-5ae7c42cc666"); }
        }
    }
}