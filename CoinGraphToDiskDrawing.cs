using Grasshopper;
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
    public class CoinGraphToDiskDrawing : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CoinGraphToDiskDrawing class.
        /// </summary>
        public CoinGraphToDiskDrawing()
          : base("CoinGraphToDiskDrawing", "DiskDrawing",
              "绘制硬币图",
              "VolumeGeneratorBasedOnGraph", "GraphDrawing")
        {
            ConnectivityEdges = new List<Line>();
            AdjacencyEdges = new List<Line>();

            VolumeNodeBubbles = new List<Brep>();
            BoundaryNodeBubbles = new List<Brep>();

            VolumeNodeMat = new List<DisplayMaterial>();
            BoundaryNodeMat = new List<DisplayMaterial>();

            // Vertices = new PointCloud();

            VolumeNodeTexts = new List<string>();
            VolumeNodeLocations = new List<Plane>();
            VolumeNodeTextSize = new List<double>();

            BoundaryNodeTexts = new List<string>();
            BoundaryNodeLocations = new List<Plane>();
            BoundaryNodeTextSize = new List<double>();
        }

        private int Thickness;

        private List<Line> ConnectivityEdges;
        private List<Line> AdjacencyEdges;

        private List<Brep> VolumeNodeBubbles;
        private List<Brep> BoundaryNodeBubbles;

        private List<DisplayMaterial> VolumeNodeMat;
        private List<DisplayMaterial> BoundaryNodeMat;

        // private PointCloud Vertices;

        private List<string> VolumeNodeTexts;
        private List<Plane> VolumeNodeLocations;
        private List<double> VolumeNodeTextSize;
        private List<string> BoundaryNodeTexts;
        private List<Plane> BoundaryNodeLocations;
        private List<double> BoundaryNodeTextSize;

        public override BoundingBox ClippingBox
        {
            get
            {
                return BoundingBox.Empty;
            }

        }



        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPlaneParameter("Plane", "P", "A plane for placing the disk-graph drawing", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddIntegerParameter("ConnectivityArributeTree", "ConnectivityTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "BoundaryAdjacencyTree", "体量与边界的邻接关系树", GH_ParamAccess.tree);
            pManager.AddPointParameter("VolumeNodePoints", "VolumeNode", "能够代表图（graph）中节点（node）的抽象点（point）", GH_ParamAccess.list);
            pManager.AddPointParameter("BoundaryNodePoints", "BoundaryNode", "用抽象点（point）表示的边界（node）", GH_ParamAccess.list);
            pManager.AddGenericParameter("VolumeNodeAttributes", "VolumeNodeAttributes", "A list of lists containing name, area and color attributes of the disks", GH_ParamAccess.list);
            pManager.AddGenericParameter("BoundaryNodeAttributes", "BoundaryAttributes", "边界节点所包含的属性", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Thickness", "Th", "Thickness of graph links in the diagram", GH_ParamAccess.item, 1);
            pManager[0].Optional = true;
            pManager[7].Optional = true;
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
            Plane location = default(Plane);

            GH_Structure<GH_Integer> gh_Structure_ConnectivityGraph = null;                     // gh_Structure -> gh_Structure_ConnectivityGraph
            GH_Structure<GH_Integer> gh_Structure_AdjacencyGraph = null;
            DataTree<int> connectivityTree = new DataTree<int>();
            DataTree<int> adjacencyTree = new DataTree<int>();

            List<Point3d> volumeNodeList = new List<Point3d>();                               // list -> volumeNodeList
            List<Point3d> boundaryNodeList = new List<Point3d>();
            List<NodeAttribute> volunmeNodeAttributes = new List<NodeAttribute>();            // list2 -> volunmeNodeAttributes
            List<NodeAttribute> boundaryNodeAttributes = new List<NodeAttribute>();

            List<string> volumeLabelList = new List<string>();                              // list3 -> volumeLabelList
            List<string> boundaryLabelList = new List<string>();
            List<double> volumeAreaList = new List<double>();                               // list4 -> volumeAreaList
            List<double> boundaryAreaList = new List<double>();

            int thickness = 0;                                                            // num -> thickness
            List<Color> brepColor = new List<Color>();                                          // list5 -> brepColor

            if (DA.GetDataTree<GH_Integer>("ConnectivityArributeTree", out gh_Structure_ConnectivityGraph)
                & DA.GetDataTree<GH_Integer>("BoundaryAdjacencyTree", out gh_Structure_AdjacencyGraph)
                & DA.GetDataList<Point3d>("VolumeNodePoints", volumeNodeList)
                & DA.GetDataList<Point3d>("BoundaryNodePoints", boundaryNodeList)
                & DA.GetDataList("VolumeNodeAttributes", volunmeNodeAttributes)
                & DA.GetDataList("BoundaryNodeAttributes", boundaryNodeAttributes))
            {
                DA.GetData<Plane>("Plane", ref location);
                DA.GetData<int>("Thickness", ref thickness);

                for (int i = 0; i < volunmeNodeAttributes.Count; i++)
                {
                    volumeLabelList.Add(volunmeNodeAttributes[i].NodeLabel);
                    volumeAreaList.Add(volunmeNodeAttributes[i].NodeArea);
                }
                for (int i = 0; i < boundaryNodeAttributes.Count; i++)
                {
                    boundaryLabelList.Add(boundaryNodeAttributes[i].NodeLabel);
                    boundaryAreaList.Add(boundaryNodeAttributes[i].NodeArea);
                }

                // List<List<int>> list6 = new List<List<int>>();                      // list6 -> (DataTree) connectivityTree
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_ConnectivityGraph, ref connectivityTree);
                UtilityFunctions.GH_StructureToDataTree_Int(gh_Structure_AdjacencyGraph, ref adjacencyTree);

                List<Point3d> volumeNodeListRelocated = new List<Point3d>();                              // list8 -> volumeNodeListRelocated
                List<Point3d> boundaryNodeListRelocated = new List<Point3d>();
                volumeNodeListRelocated = UtilityFunctions.Relocate(volumeNodeList, Plane.WorldXY, location);
                boundaryNodeListRelocated = UtilityFunctions.Relocate(boundaryNodeList, Plane.WorldXY, location);


                // List<List<Line>> listOfListForLine = new List<List<Line>>();                    // list9 -> listOfListForLine

                List<Line> connectivityEdgeList = new List<Line>();                               // list10 -> graphEdgeList
                List<Line> adjacencyEdgeList = new List<Line>();

                // 对Connectivity中的Edge进行连线
                connectivityEdgeList = UtilityFunctions.ConnectivityGraphEdgeToLine(connectivityTree, volumeNodeListRelocated);
                // 对adjacency中的Edge进行连线
                adjacencyEdgeList = UtilityFunctions.AdjacencyGraphEdgeToLine(adjacencyTree, volumeNodeListRelocated, boundaryNodeListRelocated);


                List<double> volumeNodeRadiusListForLabelSize = new List<double>();                                       // list11 -> nodeRadiusList
                for (int i = 0; i < connectivityTree.BranchCount; i++)
                {
                    for (int j = 0; j < connectivityTree.Branch(i).Count; j++)
                    {
                        volumeNodeRadiusListForLabelSize.Add(Math.Sqrt(0.2 * volumeAreaList[connectivityTree.Branch(i)[j]]) / Math.PI);
                    }
                }
                List<double> boundaryNodeRadiusListForLabelSize = new List<double>();
                for (int i = 0; i < adjacencyTree.BranchCount; i++)
                {
                    for (int j = 0; j < adjacencyTree.Branch(j).Count; j++)
                    {
                        boundaryNodeRadiusListForLabelSize.Add(Math.Sqrt(0.2 * boundaryAreaList[adjacencyTree.Branch(j)[j] + NodeAttribute.BoundaryNodeCount] / Math.PI));
                    }
                }


                List<Plane> planeBasedVolumeNodeList = new List<Plane>();                                 // list12 -> planeBasedVolumeNodeList
                foreach (Point3d point3d in volumeNodeListRelocated)
                {
                    Plane plane = new Plane(point3d, location.XAxis, location.YAxis);
                    planeBasedVolumeNodeList.Add(plane);
                }
                List<Plane> planeBasedBoundaryNodeList = new List<Plane>();
                foreach (Point3d point3d in boundaryNodeListRelocated)
                {
                    Plane plane = new Plane(point3d, location.XAxis, location.YAxis);
                    planeBasedBoundaryNodeList.Add(plane);
                }


                List<double> volumeNodeRadius = new List<double>();                                      // list14 -> nodeRadiusList2
                List<Curve> volumeNodeCircleList = new List<Curve>();                                             // list15 -> circleList
                List<Brep> volumeNodeCircleBrepList = new List<Brep>();                                           // list16 -> circleBrepList
                // List<DisplayMaterial> volumeNodeCircleMatList = new List<DisplayMaterial>();                      // list17 -> circleMatList

                for (int i = 0; i < volumeLabelList.Count; i++)
                {
                    volumeNodeRadius.Add(Math.Sqrt(volumeAreaList[i] / Math.PI));

                    Plane worldXY = Plane.WorldXY;
                    worldXY.Origin = volumeNodeListRelocated[i];

                    Circle circle = new Circle(worldXY, volumeNodeRadius[i]);
                    volumeNodeCircleList.Add(circle.ToNurbsCurve());
                    Brep brep = Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(volumeNodeCircleList[i], volumeNodeListRelocated[i]));
                    brep.Flip();
                    volumeNodeCircleBrepList.Add(brep);
                    // circleMatList.Add(new DisplayMaterial(brepColor[i]));
                }


                List<double> boundaryNodeRadius = new List<double>();
                List<Curve> boundaryNodeCircleList = new List<Curve>();
                List<Brep> boundaryNodeCircleBrepList = new List<Brep>();

                for (int i = 0; i < boundaryLabelList.Count; i++)
                {
                    boundaryNodeRadius.Add(Math.Sqrt(boundaryAreaList[i] / Math.PI));

                    Plane worldXY = Plane.WorldXY;
                    worldXY.Origin = boundaryNodeListRelocated[i];

                    Circle circle = new Circle(worldXY, boundaryNodeRadius[i]);
                    boundaryNodeCircleList.Add(circle.ToNurbsCurve());
                    Brep brep = Brep.CreateFromSurface(Surface.CreateExtrusionToPoint(boundaryNodeCircleList[i], boundaryNodeListRelocated[i]));
                    brep.Flip();
                    boundaryNodeCircleBrepList.Add(brep);
                }

                Thickness = 0;

                ConnectivityEdges.Clear();
                AdjacencyEdges.Clear();

                VolumeNodeBubbles.Clear();
                BoundaryNodeBubbles.Clear();

                VolumeNodeMat.Clear();
                BoundaryNodeMat.Clear();

                VolumeNodeLocations.Clear();
                VolumeNodeTextSize.Clear();
                VolumeNodeTexts.Clear();
                BoundaryNodeLocations.Clear();
                BoundaryNodeTextSize.Clear();
                BoundaryNodeTexts.Clear();


                Thickness = thickness;
                ConnectivityEdges.AddRange(connectivityEdgeList);
                AdjacencyEdges.AddRange(adjacencyEdgeList);

                VolumeNodeBubbles.AddRange(volumeNodeCircleBrepList);
                BoundaryNodeBubbles.AddRange(boundaryNodeCircleBrepList);

                // Mats.AddRange(circleMatList);

                VolumeNodeLocations.AddRange(planeBasedVolumeNodeList);
                VolumeNodeTextSize.AddRange(volumeNodeRadiusListForLabelSize);
                VolumeNodeTexts.AddRange(volumeLabelList);
                BoundaryNodeLocations.AddRange(planeBasedBoundaryNodeList);
                BoundaryNodeTextSize.AddRange(boundaryNodeRadiusListForLabelSize);
                BoundaryNodeTexts.AddRange(boundaryLabelList);
            }
        }

        //public List<Point3d> Relocate(List<Point3d> p, Plane location1, Plane location2)
        //{
        //    Point3d point3d = new Point3d(0.0, 0.0, 0.0);
        //    List<Point3d> list = new List<Point3d>();
        //    foreach (Point3d point3d2 in p)
        //    {
        //        point3d += point3d2;
        //    }
        //    point3d /= (double)p.Count;
        //    location1.Origin = point3d;
        //    foreach (Point3d item in p)
        //    {
        //        item.Transform(Transform.PlaneToPlane(location1, location2));
        //        list.Add(item);
        //    }
        //    return list;
        //}

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            args.Display.DrawLines(ConnectivityEdges, Color.Black, Thickness);
            args.Display.DrawLines(AdjacencyEdges, Color.Gray, Thickness);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            if (VolumeNodeBubbles.Count != 0)
            {
                for (int i = 0; i < VolumeNodeBubbles.Count; i++)
                {
                    // DisplayMaterial displayMaterial = Mats[i];
                    args.Display.DrawBrepShaded(VolumeNodeBubbles[i], new DisplayMaterial(Color.Red));
                    args.Display.Draw3dText(VolumeNodeTexts[i], Color.White, VolumeNodeLocations[i], VolumeNodeTextSize[i], "Arial",
                        true, false,
                        Rhino.DocObjects.TextHorizontalAlignment.Center,
                        Rhino.DocObjects.TextVerticalAlignment.Middle);
                }
            }
            if (BoundaryNodeBubbles.Count != 0)
            {
                for (int i = 0; i < BoundaryNodeBubbles.Count; i++)
                {
                    args.Display.DrawBrepWires(BoundaryNodeBubbles[i], Color.Blue);
                    args.Display.Draw3dText(BoundaryNodeTexts[i], Color.White, BoundaryNodeLocations[i], BoundaryNodeTextSize[i], "Arial",
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
            get { return new Guid("162c046f-7649-4e30-945f-3dc0765808e2"); }
        }
    }
}