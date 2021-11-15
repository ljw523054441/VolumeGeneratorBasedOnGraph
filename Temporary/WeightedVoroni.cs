using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class WeightedVoroni : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the WeightedVoroni class.
        /// </summary>
        public WeightedVoroni()
          : base("WeightedVoroni", "WeightedVoroni",
              "Draws a Voroni diagram for points with various radii on a base plane; this diagram will be sort of a bubble diagram",
              "VolumeGeneratorBasedOnGraph", "Graph")
        {
            Vertices = new List<Point3d>();
            Indices = new List<int>();
            Height = 1.0;
            base.Message = "开发中";
        }

        private List<Point3d> Vertices;
        private List<int> Indices;
        private double Height;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("BoundingCurve", "BC", "A polyline confining the area of the diagram", GH_ParamAccess.item);

            pManager.AddPlaneParameter("BasePlane", "BP", "The plane on which the diagram is to be drawn", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddPointParameter("Points", "P", "The centre points of bubbles/Voronoi cells", GH_ParamAccess.list);
            pManager.AddGenericParameter("Attributes", "R", "The list of radii of cells correspoding to the list of points", GH_ParamAccess.list);
            pManager.AddNumberParameter("Scale", "S", "", GH_ParamAccess.item);
            // pManager.AddCurveParameter("BoundingCurve", "BC", "A polyline confining the area of the diagram", GH_ParamAccess.item);
            // pManager[4].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("AdjacencyGraph", "AG", "A graph (adjacency lists stored in a datatree) describing how the cells are adjacent to each other", GH_ParamAccess.tree);
            pManager.AddCurveParameter("VoronoiCells", "VC", "Bubbles in the form of Voronoi cells", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 场地边界
            Curve boundaryCurve = null;
            DA.GetData<Curve>("BoundingCurve", ref boundaryCurve);
            if (boundaryCurve == null || !boundaryCurve.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the bounding curve is not closed!");
                return;
            }
            if (!boundaryCurve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "the bounding curve is not planar!");
                return;
            }
            // 场地总面积groundArea
            double groundArea = AreaMassProperties.Compute(boundaryCurve).Area;



            List<Point3d> pointList = new List<Point3d>();
            List<GraphNodeAttribute> nodeAttributes = new List<GraphNodeAttribute>();
            Plane basePlace = Plane.WorldXY;
            Curve boundingCurve = null;
            DataTree<int> adjacencyGraph = null;

            List<double> radius = new List<double>();

            double scale = 0.0;

            if (DA.GetDataList<Point3d>("Points", pointList) 
                & DA.GetDataList<GraphNodeAttribute>("Attributes", nodeAttributes))
            {
                DA.GetData<Plane>("BasePlane", ref basePlace);
                DA.GetData<Curve>("BoundingCurve", ref boundingCurve);
                DA.GetData("Scale", ref scale);

                for (int i = 0; i < nodeAttributes.Count; i++)
                {
                    // radius.Add(Math.Sqrt(nodeAttributes[i].NodeAreaProportion * groundArea / Math.PI) * scale);
                    radius.Add(Math.Sqrt(nodeAttributes[i].NodeAreaProportion * groundArea / Math.PI));
                }

                List<Curve> cellCurve = AlphaVoronoi(basePlace, pointList, radius, boundingCurve, ref adjacencyGraph, scale);

                DA.SetDataList("VoronoiCells", cellCurve);
                DA.SetDataTree(0, adjacencyGraph);

            }
        }

        /// <summary>
        /// 生成加权Voronoi
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
            Curve boundingCurve, 
            ref DataTree<int> adjacencyGraph,
            double scale)
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
                        Vertices.Clear();
                        Indices.Clear();

                        List<Curve> pointBubbleCurves = new List<Curve>();
                        List<List<Curve>> pointBubbleCurvesLoL = new List<List<Curve>>();
                        DataTree<int> adjacencyDataTree = new DataTree<int>();
                        List<List<bool>> dominatedFlagaLoL = new List<List<bool>>();
                        List<Curve> cellCurveList = new List<Curve>();
                        List<Line> verticalLineList = new List<Line>();

                        BoundingBox pointsBoundingBox = new BoundingBox(points);

                        Plane worldXY = Plane.WorldXY;
                        for (int i = 0; i < points.Count; i++)
                        {
                            worldXY.Origin = points[i];
                            // List<Curve> list6 = list;
                            Circle bubble = new Circle(worldXY, radius[i] * scale);
                            // Circle circle2 = circle;
                            pointBubbleCurves.Add(bubble.ToNurbsCurve());

                            pointsBoundingBox = BoundingBox.Union(pointBubbleCurves[i].GetBoundingBox(false), pointsBoundingBox);

                            Indices.Add(i);

                            adjacencyDataTree.EnsurePath(i);

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
                                    adjacencyDataTree.Branch(i).Add(j);
                                    bool dominated = false;
                                    bool outside = false;
                                    double distanceBetweenP1andVertical = 0;

                                    // 求短的verticalLines
                                    Line verticalLine = DrawVerticalLinesBetweenTwoCells(
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
                                        Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, 
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


                        adjacencyGraph = adjacencyDataTree;
                        result = cellCurveList;
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// 绘制两个cell之间的voronoi垂直线
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
        public Line DrawVerticalLinesBetweenTwoCells(
            Point3d point1, 
            Point3d point2 ,
            double radius1,
            double radius2, 
            Plane basePlane, 
            ref double distanceBetweenP1andVertical, 
            Curve boundingCurve, 
            ref bool dominated, 
            ref bool outside)
        {
            double distance = point1.DistanceTo(point2);
            double x1 = distance / 2.0 + (Math.Pow(radius1, 2.0) - Math.Pow(radius2, 2.0)) / (2.0 * distance);
            double x2 = distance / 2.0 - (Math.Pow(radius1, 2.0) - Math.Pow(radius2, 2.0)) / (2.0 * distance);

            Line line = new Line(point1, point2);
            distanceBetweenP1andVertical = x1;
            Point3d pointBetweenTwoCell = line.PointAt(distanceBetweenP1andVertical / distance);

            bool flag = boundingCurve.Contains(pointBetweenTwoCell, basePlane, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) != PointContainment.Inside
                || boundingCurve.Contains(pointBetweenTwoCell, basePlane, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == PointContainment.Coincident;
            if (flag)
            {
                outside = true;
            }

            flag = distanceBetweenP1andVertical > Math.Max(radius1, radius2);
            if (flag)
            {
                dominated = true;
            }

            // Point3d point3d2 = pointBetweenTwoCell;
            Vector3d lineDirection = new Vector3d(point2-point1);
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
            get { return new Guid("18dd9201-de43-4ae4-88c4-be8a5db150ec"); }
        }
    }
}