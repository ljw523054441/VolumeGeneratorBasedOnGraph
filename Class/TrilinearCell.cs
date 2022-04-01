using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class TrilinearCell : Cell
    {
        public Line MiddleBaseLine { get; set; } = Line.Unset;
        public Line WestBaseLine1 { get; set; } = Line.Unset;
        public Line EastBaseLine1 { get; set; } = Line.Unset;

        public Line HSouthBaseLine { get; set; } = Line.Unset;
        public Line VSouthBaseLine { get; set; } = Line.Unset;
        public Line HNorthBaseLine { get; set; } = Line.Unset;
        public Line VNorthBaseLine { get; set; } = Line.Unset;

        public List<Interval> Middle_Interval { get; set; }
        public List<Interval> West1_Interval { get; set; }
        public List<Interval> East1_Interval { get; set; }
        public List<Interval> HSouth_Interval { get; set; }
        public List<Interval> HNouth_Interval { get; set; }
        public List<Interval> VSouth_Interval { get; set; }
        public List<Interval> VNorth_Interval { get; set; }

        public double CenterT0Value { get; set; }
        public double CenterT1Value { get; set; }

        public int BoundaryBaseLineCount { get; set; }
        public int CenterBaseLineCount { get; set; }


        public Curve[] crvForScaled { get; set; }

        public double[] MinScaleFactors { get; set; }

        public enum TrilinearShapeType
        {
            Unset = 0,
            EShape = 1,
            EVarientShape = 2,
            ZShape = 3,
            EightVariantShape = 4,
            EightShape = 5,
            SixVariantShape = 6,
            SixShape = 7,

            Scaled = 8
        }

        public TrilinearShapeType ShapeType { get; set; }

        public TrilinearCell(TrilinearCell source)
        {
            #region baseLine
            this.SouthBaseLine = source.SouthBaseLine;
            this.NorthBaseLine = source.NorthBaseLine;
            this.WestBaseLine = source.WestBaseLine;
            this.EastBaseLine = source.EastBaseLine;

            this.MiddleBaseLine = source.MiddleBaseLine;
            this.WestBaseLine1 = source.WestBaseLine1;
            this.EastBaseLine1 = source.EastBaseLine1;
            this.HSouthBaseLine = source.HSouthBaseLine;
            this.VSouthBaseLine = source.HSouthBaseLine;
            this.HNorthBaseLine = source.HNorthBaseLine;
            this.VNorthBaseLine = source.VNorthBaseLine;
            #endregion

            #region interval
            if (source.South_Interval == null)
            {
                this.South_Interval = null;
            }
            else
            {
                this.South_Interval = new List<Interval>();
                this.South_Interval.AddRange(source.South_Interval);
            }
            if (source.North_Interval == null)
            {
                this.North_Interval = null;
            }
            else
            {
                this.North_Interval = new List<Interval>();
                this.North_Interval.AddRange(source.North_Interval);
            }
            if (source.West_Interval == null)
            {
                this.West_Interval = null;
            }
            else
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.AddRange(source.West_Interval);
            }
            if (source.East_Interval == null)
            {
                this.East_Interval = null;
            }
            else
            {
                this.East_Interval = new List<Interval>();
                this.East_Interval.AddRange(source.East_Interval);
            }

            if (source.Middle_Interval == null)
            {
                this.Middle_Interval = null;
            }
            else
            {
                this.Middle_Interval = new List<Interval>();
                this.Middle_Interval.AddRange(source.Middle_Interval);
            }
            if (source.West1_Interval == null)
            {
                this.West1_Interval = null;
            }
            else
            {
                this.West1_Interval = new List<Interval>();
                this.West1_Interval.AddRange(source.West1_Interval);
            }
            if (source.East1_Interval == null)
            {
                this.East1_Interval = null;
            }
            else
            {
                this.East1_Interval = new List<Interval>();
                this.East1_Interval.AddRange(source.East1_Interval);
            }
            if (source.HSouth_Interval == null)
            {
                this.HSouth_Interval = null;
            }
            else
            {
                this.HSouth_Interval = new List<Interval>();
                this.HSouth_Interval.AddRange(source.HSouth_Interval);
            }
            if (source.VSouth_Interval == null)
            {
                this.VSouth_Interval = null;
            }
            else
            {
                this.VSouth_Interval = new List<Interval>();
                this.VSouth_Interval.AddRange(source.VSouth_Interval);
            }
            if (source.HNouth_Interval == null)
            {
                this.HNouth_Interval = null;
            }
            else
            {
                this.HNouth_Interval = new List<Interval>();
                this.HNouth_Interval.AddRange(source.HNouth_Interval);
            }
            if (source.VNorth_Interval == null)
            {
                this.VNorth_Interval = null;
            }
            else
            {
                this.VNorth_Interval = new List<Interval>();
                this.VNorth_Interval.AddRange(source.VNorth_Interval);
            }
            #endregion

            this.CenterT0Value = source.CenterT0Value;
            this.CenterT1Value = source.CenterT1Value;
            this.BoundaryBaseLineCount = source.BoundaryBaseLineCount;
            this.CenterBaseLineCount = source.CenterBaseLineCount;
            this.ShapeType = source.ShapeType;

            if (source.CellBoundary == null)
            {
                this.CellBoundary = null;
            }
            else
            {
                this.CellBoundary = source.CellBoundary.DuplicateCurve();
            }
            

            if (source.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[source.CellBreps_Multistorey.Length];
                for (int i = 0; i < source.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = source.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (source.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = source.crvForHighRise.DuplicateCurve();
            }

            if (source.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[source.CellBreps_Highrise.Length];
                for (int i = 0; i < source.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = source.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = source.IsLengthConstraint;

            this.PublicBaseLineIndex = source.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = source.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(source.TurningPoint);
        }

        public TrilinearCell(Curve southLine, 
                             Curve middleLine, 
                             Curve northLine,
                             Curve westLine, 
                             Curve westLine1, 
                             Curve eastLine,
                             Curve eastLine1,

                             Curve southBoundary,
                             Curve northBoundary,
                             Curve westBoundary,
                             Curve eastBoundary,
                             int turningPointIndex)
        {
            #region BaseLine 和 interval
            int count = 0;
            if (southLine != null)
            {
                this.SouthBaseLine = new Line(southLine.PointAtStart, southLine.PointAtEnd);
                count++;

                this.South_Interval = new List<Interval>();
                this.South_Interval.Add(new Interval(0, 1));
            }
            if (northLine != null)
            {
                this.NorthBaseLine = new Line(northLine.PointAtStart, northLine.PointAtEnd);
                count++;

                this.North_Interval = new List<Interval>();
                this.North_Interval.Add(new Interval(0, 1));
            }
            if (middleLine != null)
            {
                this.MiddleBaseLine = new Line(middleLine.PointAtStart, middleLine.PointAtEnd);
                count++;

                this.Middle_Interval = new List<Interval>();
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            if (westLine != null)
            {
                this.WestBaseLine = new Line(westLine.PointAtStart, westLine.PointAtEnd);
                count++;

                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
            }
            if (eastLine != null)
            {
                this.EastBaseLine = new Line(eastLine.PointAtStart, eastLine.PointAtEnd);
                count++;

                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));
            }
            if (westLine1 != null)
            {
                this.WestBaseLine1 = new Line(westLine1.PointAtStart, westLine1.PointAtEnd);
                count++;

                this.West1_Interval = new List<Interval>();
                this.West1_Interval.Add(new Interval(0, 1));
            }
            if (eastLine1 != null)
            {
                this.EastBaseLine1 = new Line(eastLine1.PointAtStart, eastLine1.PointAtEnd);
                count++;

                this.East1_Interval = new List<Interval>();
                this.East1_Interval.Add(new Interval(0, 1));
            }
            #endregion

            this.ShapeType = TrilinearShapeType.EightShape;
            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            List<Curve> extendedBoundaryList = new List<Curve>();
            extendedBoundaryList.Add(southBoundary);
            extendedBoundaryList.Add(eastBoundary);
            //extendedBoundaryList.Add(east1Boundary);
            extendedBoundaryList.Add(northBoundary);
            //extendedBoundaryList.Add(west1Boundary);
            extendedBoundaryList.Add(westBoundary);

            List<Curve> boundaryList = new List<Curve>();
            for (int k = 0; k < extendedBoundaryList.Count; k++)
            {
                CurveIntersections intersectionEvents0 = Intersection.CurveCurve(extendedBoundaryList[k], extendedBoundaryList[((k - 1) + extendedBoundaryList.Count) % extendedBoundaryList.Count], Cell.Tolerance, Cell.Tolerance);
                CurveIntersections intersectionEvents1 = Intersection.CurveCurve(extendedBoundaryList[k], extendedBoundaryList[(k + 1) % extendedBoundaryList.Count], Cell.Tolerance, Cell.Tolerance);

                Point3d p0;
                Point3d p1;
                if (intersectionEvents0[0].IsOverlap)
                {
                    Point3d a0 = intersectionEvents0[0].PointA;
                    Point3d a1 = intersectionEvents0[0].PointA2;
                    p0 = (a0 + a1) / 2;
                }
                else
                {
                    p0 = intersectionEvents0[0].PointA;
                }

                if (intersectionEvents1[0].IsOverlap)
                {
                    Point3d a0 = intersectionEvents1[0].PointA;
                    Point3d a1 = intersectionEvents1[0].PointA2;
                    p1 = (a0 + a1) / 2;
                }
                else
                {
                    p1 = intersectionEvents1[0].PointA;
                }

                boundaryList.Add(new Line(p0, p1).ToNurbsCurve());
            }

            List<Point3d> pts = new List<Point3d>();
            for (int k = 0; k < boundaryList.Count; k++)
            {
                pts.Add(boundaryList[k].PointAtStart);
            }
            pts.Add(pts[0]);

            //Curve[] joinedCrv = Curve.JoinCurves(boundaryList);
            //this.CellBoundary = joinedCrv[0];
            this.CellBoundary = new Polyline(pts).ToNurbsCurve();
            #endregion

            if (turningPointIndex == 0)
            {
                this.TurningPoint = pts[0];
            }
            else if (turningPointIndex == 1)
            {
                this.TurningPoint = pts[1];
            }
            else if (turningPointIndex == 2)
            {
                this.TurningPoint = pts[2];
            }
            else
            {
                this.TurningPoint = pts[3];
            }
        }

        /// <summary>
        /// 生成EShape或ZShape
        /// </summary>
        private TrilinearCell(TrilinearCell initialCell, bool isEShape, int directionCode)
        {
            #region Interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval.AddRange(initialCell.South_Interval);
            }

            this.Middle_Interval = new List<Interval>();
            if (initialCell.Middle_Interval == null)
            {
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.Middle_Interval.AddRange(initialCell.Middle_Interval);
            }

            this.North_Interval = new List<Interval>();
            if (initialCell.North_Interval == null)
            {
                this.North_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.North_Interval.AddRange(initialCell.North_Interval);
            }


            if (isEShape)
            {
                if (directionCode < 2)
                {
                    // 西侧
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.West1_Interval = new List<Interval>();
                    this.West1_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    // 东侧
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));
                    this.East1_Interval = new List<Interval>();
                    this.East1_Interval.Add(new Interval(0, 1));
                }
            }
            else
            {
                if (directionCode < 2)
                {
                    // 西南
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.East1_Interval = new List<Interval>();
                    this.East1_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    // 东北
                    this.West1_Interval = new List<Interval>();
                    this.West1_Interval.Add(new Interval(0, 1));
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));
                }
            }
            #endregion

            #region BaseLine
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;
            #endregion

            if (isEShape)
            {
                this.ShapeType = TrilinearShapeType.EShape;
            }
            else
            {
                this.ShapeType = TrilinearShapeType.ZShape;
            }
            
            this.BoundaryBaseLineCount = 5;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            #endregion

            if (initialCell.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[initialCell.CellBreps_Multistorey.Length];
                for (int i = 0; i < initialCell.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = initialCell.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (initialCell.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = initialCell.crvForHighRise.DuplicateCurve();
            }

            if (initialCell.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[initialCell.CellBreps_Highrise.Length];
                for (int i = 0; i < initialCell.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = initialCell.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = initialCell.IsLengthConstraint;

            this.PublicBaseLineIndex = initialCell.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = initialCell.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(initialCell.TurningPoint);
            this.PrevBaseLineIndexRelatedToCutPoint = initialCell.PrevBaseLineIndexRelatedToCutPoint;
            this.NextBaseLineIndexRelatedToCutPoint = initialCell.NextBaseLineIndexRelatedToCutPoint;
        }

        /// <summary>
        /// 生成EVariantShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="t0"></param>
        /// <param name="t1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, int directionCode, double t0, double t1)
        {
            #region Interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval.AddRange(initialCell.South_Interval);
            }

            this.Middle_Interval = new List<Interval>();
            if (initialCell.Middle_Interval == null)
            {
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.Middle_Interval.AddRange(initialCell.Middle_Interval);
            }

            this.North_Interval = new List<Interval>();
            if (initialCell.North_Interval == null)
            {
                this.North_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.North_Interval.AddRange(initialCell.North_Interval);
            }

            this.VSouth_Interval = new List<Interval>();
            this.VSouth_Interval.Add(new Interval(0, 1));
            this.VNorth_Interval = new List<Interval>();
            this.VNorth_Interval.Add(new Interval(0, 1));
            #endregion

            #region BaseLine
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;
            #region VBaseLine
            if (directionCode < 2)
            {
                // 两根中线不共线
                Point3d vStart0 = this.SouthBaseLine.PointAt(t0);
                Point3d vEnd0 = this.MiddleBaseLine.PointAt(t0);
                Point3d vStart1 = this.MiddleBaseLine.PointAt(t1);
                Point3d vEnd1 = this.NorthBaseLine.PointAt(t1);

                this.VSouthBaseLine = new Line(vStart0, vEnd0);
                this.VNorthBaseLine = new Line(vStart1, vEnd1);
            }
            else
            {
                // 两根中线共线
                Point3d vStart0 = this.SouthBaseLine.PointAt(t0);
                Point3d vEnd0 = this.MiddleBaseLine.PointAt(t0);
                Point3d vStart1 = this.MiddleBaseLine.PointAt(t0);
                Point3d vEnd1 = this.NorthBaseLine.PointAt(t0);

                this.VSouthBaseLine = new Line(vStart0, vEnd0);
                this.VNorthBaseLine = new Line(vStart1, vEnd1);
            }
            #endregion
            #endregion

            this.ShapeType = TrilinearShapeType.EVarientShape;
            this.BoundaryBaseLineCount = 3;
            this.CenterBaseLineCount = 2;

            #region CellBoundary
            this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            #endregion

            if (initialCell.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[initialCell.CellBreps_Multistorey.Length];
                for (int i = 0; i < initialCell.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = initialCell.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (initialCell.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = initialCell.crvForHighRise.DuplicateCurve();
            }

            if (initialCell.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[initialCell.CellBreps_Highrise.Length];
                for (int i = 0; i < initialCell.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = initialCell.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = initialCell.IsLengthConstraint;

            this.PublicBaseLineIndex = initialCell.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = initialCell.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(initialCell.TurningPoint);
            this.PrevBaseLineIndexRelatedToCutPoint = initialCell.PrevBaseLineIndexRelatedToCutPoint;
            this.NextBaseLineIndexRelatedToCutPoint = initialCell.NextBaseLineIndexRelatedToCutPoint;
        }

        /// <summary>
        /// 生成EightShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, TrilinearShapeType typeForEight)
        {
            #region Interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval.AddRange(initialCell.South_Interval);
            }
            this.Middle_Interval = new List<Interval>();
            if (initialCell.Middle_Interval == null)
            {
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.Middle_Interval.AddRange(initialCell.Middle_Interval);
            }
            this.North_Interval = new List<Interval>();
            if (initialCell.North_Interval == null)
            {
                this.North_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.North_Interval.AddRange(initialCell.North_Interval);
            }


            this.West_Interval = new List<Interval>();
            this.West_Interval.Add(new Interval(0, 1));
            this.West1_Interval = new List<Interval>();
            this.West1_Interval.Add(new Interval(0, 1));
            this.East_Interval = new List<Interval>();
            this.East_Interval.Add(new Interval(0, 1));
            this.East1_Interval = new List<Interval>();
            this.East1_Interval.Add(new Interval(0, 1));
            #endregion

            #region baseLine
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;
            #endregion

            this.ShapeType = TrilinearShapeType.EightShape;
            this.BoundaryBaseLineCount = 6;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            #endregion

            if (initialCell.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[initialCell.CellBreps_Multistorey.Length];
                for (int i = 0; i < initialCell.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = initialCell.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (initialCell.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = initialCell.crvForHighRise.DuplicateCurve();
            }

            if (initialCell.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[initialCell.CellBreps_Highrise.Length];
                for (int i = 0; i < initialCell.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = initialCell.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = initialCell.IsLengthConstraint;

            this.PublicBaseLineIndex = initialCell.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = initialCell.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(initialCell.TurningPoint);
            this.PrevBaseLineIndexRelatedToCutPoint = initialCell.PrevBaseLineIndexRelatedToCutPoint;
            this.NextBaseLineIndexRelatedToCutPoint = initialCell.NextBaseLineIndexRelatedToCutPoint;
        }

        /// <summary>
        /// 生成EightVariantShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, double t)
        {
            #region Interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval = new List<Interval>();
                this.South_Interval.AddRange(initialCell.South_Interval);
            }

            this.Middle_Interval = new List<Interval>();
            if (initialCell.Middle_Interval == null)
            {
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.Middle_Interval.AddRange(initialCell.Middle_Interval);
            }

            this.North_Interval = new List<Interval>();
            if (initialCell.North_Interval == null)
            {
                this.North_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.North_Interval.AddRange(initialCell.North_Interval);
            }

            this.West_Interval = new List<Interval>();
            this.West_Interval.Add(new Interval(0, 1));
            this.East_Interval = new List<Interval>();
            this.East_Interval.Add(new Interval(0, 1));
            this.West1_Interval = new List<Interval>();
            this.West1_Interval.Add(new Interval(0, 1));
            this.East1_Interval = new List<Interval>();
            this.East1_Interval.Add(new Interval(0, 1));
            #endregion

            #region BaseLine
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;

            Point3d pointOnSouth0 = this.SouthBaseLine.PointAt(t);
            Point3d pointOnMiddel0 = this.MiddleBaseLine.PointAt(t);
            Point3d pointOnNorth0 = this.NorthBaseLine.PointAt(t);
            Point3d pointOnSouth1 = this.SouthBaseLine.PointAt(1 - t);
            Point3d pointOnMiddel1 = this.MiddleBaseLine.PointAt(1 - t);
            Point3d pointOnNorth1 = this.NorthBaseLine.PointAt(1 - t);
            this.WestBaseLine = new Line(pointOnSouth0, pointOnMiddel0);
            this.WestBaseLine1 = new Line(pointOnMiddel0, pointOnNorth0);
            this.EastBaseLine = new Line(pointOnSouth1, pointOnMiddel1);
            this.EastBaseLine1 = new Line(pointOnMiddel1, pointOnNorth1);

            #endregion

            this.ShapeType = TrilinearShapeType.EightVariantShape;
            this.BoundaryBaseLineCount = 7;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            #endregion

            if (initialCell.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[initialCell.CellBreps_Multistorey.Length];
                for (int i = 0; i < initialCell.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = initialCell.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (initialCell.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = initialCell.crvForHighRise.DuplicateCurve();
            }

            if (initialCell.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[initialCell.CellBreps_Highrise.Length];
                for (int i = 0; i < initialCell.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = initialCell.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = initialCell.IsLengthConstraint;

            this.PublicBaseLineIndex = initialCell.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = initialCell.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(initialCell.TurningPoint);
            this.PrevBaseLineIndexRelatedToCutPoint = initialCell.PrevBaseLineIndexRelatedToCutPoint;
            this.NextBaseLineIndexRelatedToCutPoint = initialCell.NextBaseLineIndexRelatedToCutPoint;
        }

        /// <summary>
        /// 生成SixShape或SixVariantShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="t"></param>
        private TrilinearCell(TrilinearCell initialCell, int directionCode, double t)
        {
            #region interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval = new List<Interval>();
                this.South_Interval.AddRange(initialCell.South_Interval);
            }

            this.Middle_Interval = new List<Interval>();
            if (initialCell.Middle_Interval == null)
            {
                this.Middle_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.Middle_Interval.AddRange(initialCell.Middle_Interval);
            }

            this.North_Interval = new List<Interval>();
            if (initialCell.North_Interval == null)
            {
                this.North_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.North_Interval.AddRange(initialCell.North_Interval);
            }

            if (directionCode == 0)
            {
                //this.West_Interval = new List<Interval>();
                //this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));
                this.West1_Interval = new List<Interval>();
                this.West1_Interval.Add(new Interval(0, 1));
                this.East1_Interval = new List<Interval>();
                this.East1_Interval.Add(new Interval(0, 1));
            }
            else if (directionCode == 1)
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                //this.East_Interval = new List<Interval>();
                //this.East_Interval.Add(new Interval(0, 1));
                this.West1_Interval = new List<Interval>();
                this.West1_Interval.Add(new Interval(0, 1));
                this.East1_Interval = new List<Interval>();
                this.East1_Interval.Add(new Interval(0, 1));
            }
            else if (directionCode == 2)
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));
                //this.West1_Interval = new List<Interval>();
                //this.West1_Interval.Add(new Interval(0, 1));
                this.East1_Interval = new List<Interval>();
                this.East1_Interval.Add(new Interval(0, 1));
            }
            else if (directionCode == 3)
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));
                this.West1_Interval = new List<Interval>();
                this.West1_Interval.Add(new Interval(0, 1));
                //this.East1_Interval = new List<Interval>();
                //this.East1_Interval.Add(new Interval(0, 1));
            }
            else if (directionCode != 4)
            {
                //this.West_Interval = new List<Interval>();
                //this.West_Interval.Add(new Interval(0, 1));
                //this.East_Interval = new List<Interval>();
                //this.East_Interval.Add(new Interval(0, 1));
                this.West1_Interval = new List<Interval>();
                this.West1_Interval.Add(new Interval(0, 1));
                this.East1_Interval = new List<Interval>();
                this.East1_Interval.Add(new Interval(0, 1));

                this.VSouth_Interval = new List<Interval>();
                this.VSouth_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));
                //this.West1_Interval = new List<Interval>();
                //this.West1_Interval.Add(new Interval(0, 1));
                //this.East1_Interval = new List<Interval>();
                //this.East1_Interval.Add(new Interval(0, 1));

                this.VNorth_Interval = new List<Interval>();
                this.VNorth_Interval.Add(new Interval(0, 1));
            }
            #endregion

            #region baseline
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;

            if (directionCode == 4)
            {
                Point3d vStart0 = this.SouthBaseLine.PointAt(t);
                Point3d vEnd0 = this.MiddleBaseLine.PointAt(t);

                this.VSouthBaseLine = new Line(vStart0, vEnd0);
                
            }
            else if (directionCode == 5)
            {
                Point3d vStart1 = this.MiddleBaseLine.PointAt(t);
                Point3d vEnd1 = this.NorthBaseLine.PointAt(t);
                this.VNorthBaseLine = new Line(vStart1, vEnd1);
            }
            #endregion

            if (directionCode < 4)
            {
                this.ShapeType = TrilinearShapeType.SixShape;
                this.BoundaryBaseLineCount = 6;
                this.CenterBaseLineCount = 1;
            }
            else
            {
                this.ShapeType = TrilinearShapeType.SixVariantShape;
                this.BoundaryBaseLineCount = 4;
                this.CenterBaseLineCount = 2;
            }

            #region CellBoundary
            this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            #endregion

            if (initialCell.CellBreps_Multistorey == null)
            {
                this.CellBreps_Multistorey = null;
            }
            else
            {
                this.CellBreps_Multistorey = new Brep[initialCell.CellBreps_Multistorey.Length];
                for (int i = 0; i < initialCell.CellBreps_Multistorey.Length; i++)
                {
                    this.CellBreps_Multistorey[i] = initialCell.CellBreps_Multistorey[i].DuplicateBrep();
                }
            }

            if (initialCell.crvForHighRise == null)
            {
                this.crvForHighRise = null;
            }
            else
            {
                this.crvForHighRise = initialCell.crvForHighRise.DuplicateCurve();
            }

            if (initialCell.CellBreps_Highrise == null)
            {
                this.CellBreps_Highrise = null;
            }
            else
            {
                this.CellBreps_Highrise = new Brep[initialCell.CellBreps_Highrise.Length];
                for (int i = 0; i < initialCell.CellBreps_Highrise.Length; i++)
                {
                    this.CellBreps_Highrise[i] = initialCell.CellBreps_Highrise[i].DuplicateBrep();
                }
            }

            this.IsLengthConstraint = initialCell.IsLengthConstraint;

            this.PublicBaseLineIndex = initialCell.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = initialCell.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(initialCell.TurningPoint);
            this.PrevBaseLineIndexRelatedToCutPoint = initialCell.PrevBaseLineIndexRelatedToCutPoint;
            this.NextBaseLineIndexRelatedToCutPoint = initialCell.NextBaseLineIndexRelatedToCutPoint;
        }

        public TrilinearCell GenerateEShape(int directionCode)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, true, directionCode);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateEVariantShape(int directionCode, double randomT0, double randomT1, double w)
        {
            TrilinearCell newGeneratedCell = null;

            double length0 = this.SouthBaseLine.Length;
            double length1 = this.NorthBaseLine.Length;

            double length = length0 < length1 ? length0 : length1;
            double tMin = w / length;
            double tMax = (length - w) / length;
            double t0 = Math.Abs(tMin + (tMax - tMin) * randomT0);
            double t1 = Math.Abs(tMin + (tMax - tMin) * randomT1);

            newGeneratedCell = new TrilinearCell(this, directionCode, t0, t1);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateZShape(int directionCode)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, false, directionCode);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateEightAndEightVariantShape(double randomT, int directionCode, double lMin, double w, bool flag)
        {
            TrilinearCell newGeneratedCell;
            if (directionCode < 2)
            {
                // 执行GenerateEightShape函数
                newGeneratedCell = new TrilinearCell(this, TrilinearShapeType.EightShape);
            }
            else
            {
                if (flag)
                {
                    // 执行GenerateEightShape函数
                    newGeneratedCell = new TrilinearCell(this, TrilinearShapeType.EightShape);
                }
                else
                {
                    double length0 = this.SouthBaseLine.Length;
                    double length1 = this.NorthBaseLine.Length;

                    double length = length0 < length1 ? length0 : length1;
                    if (length >= lMin + 3 * w)
                    {
                        double deltaW = length - 2 * w - w - lMin;
                        double t0 = w / length;
                        double t1 = (w + 0.5 * deltaW) / length;
                        double t = t0 + (t1 - t0) * randomT;

                        newGeneratedCell = new TrilinearCell(this, t);
                    }
                    else
                    {
                        // 执行GenerateEightShape函数
                        newGeneratedCell = new TrilinearCell(this, TrilinearShapeType.EightShape);
                    }
                }
            }
            return newGeneratedCell;
        }

        public TrilinearCell GenerateSixAndSixVariantShape(double randomT,int directionCode, double lMin, double w)
        {
            TrilinearCell newGeneratedCell;

            if (directionCode < 4)
            {
                // 生成SixShape
                newGeneratedCell = new TrilinearCell(this, directionCode, randomT);
            }
            else
            {
                double length0 = this.SouthBaseLine.Length;
                double length1 = this.NorthBaseLine.Length;

                double length = length0 < length1 ? length0 : length1;
                if (length >= lMin + 3 * w)
                {
                    // 生成SixVariantShape
                    double deltaW = length - 2 * w - w - lMin;
                    double t0 = w / length;
                    double t1 = (w + 0.5 * deltaW) / length;
                    double t = t0 + (t1 - t0) * randomT;

                    newGeneratedCell = new TrilinearCell(this, directionCode, t);
                }
                else
                {
                    // 生成SixShape
                    int newDirectionCode = directionCode - 4;
                    newGeneratedCell = new TrilinearCell(this, newDirectionCode, randomT);
                }
            }

            return newGeneratedCell;
        }

        public TrilinearCell GenerateOneRemovedEShapeOrZShapeOrEVariantShape(int directionCode)
        {
            if (this.ShapeType != TrilinearShapeType.ZShape
                || this.ShapeType != TrilinearShapeType.EShape
                || this.ShapeType != TrilinearShapeType.EVarientShape)
            {
                return this;
            }
            else
            {
                if (this.ShapeType == TrilinearShapeType.ZShape)
                {
                    if (this.West_Interval == null)
                    {
                        if (directionCode < 2)
                        {
                            this.West1_Interval = new List<Interval>();
                        }
                        else
                        {
                            this.East_Interval = new List<Interval>();
                        }
                    }
                    else if (this.East_Interval == null)
                    {
                        if (directionCode < 2)
                        {
                            this.East1_Interval = new List<Interval>();
                        }
                        else
                        {
                            this.West_Interval = new List<Interval>();
                        }
                    }
                }
                else if (this.ShapeType == TrilinearShapeType.EShape)
                {
                    if (this.West_Interval == null)
                    {
                        if (directionCode < 2)
                        {
                            this.East_Interval = new List<Interval>();
                        }
                        else
                        {
                            this.East1_Interval = new List<Interval>();
                        }
                    }
                    else if (this.East_Interval == null)
                    {
                        if (directionCode < 2)
                        {
                            this.West_Interval = new List<Interval>();
                        }
                        else
                        {
                            this.West1_Interval = new List<Interval>();
                        }
                    }
                }
                else
                {
                    if (directionCode < 2)
                    {
                        this.HSouth_Interval = new List<Interval>();
                    }
                    else
                    {
                        this.HNouth_Interval = new List<Interval>();
                    }
                }

                return this;
            }
        }

        public TrilinearCell GenerateTwoRemovedEShapeOrZShapeOrEVariantShape()
        {
            if (this.ShapeType != TrilinearShapeType.ZShape
                && this.ShapeType != TrilinearShapeType.EShape
                && this.ShapeType != TrilinearShapeType.EVarientShape)
            {
                return this;
            }
            else
            {
                if (this.ShapeType == TrilinearShapeType.ZShape)
                {
                    // 去掉两条竖
                    if (this.West_Interval == null)
                    {
                        this.West1_Interval = new List<Interval>();
                        this.East_Interval = new List<Interval>();
                    }
                    else if (this.East_Interval == null)
                    {
                        this.East1_Interval = new List<Interval>();
                        this.West_Interval = new List<Interval>();
                    }
                }
                else if (this.ShapeType == TrilinearShapeType.EShape)
                {
                    // 去掉两条竖
                    if (this.West_Interval == null)
                    {
                        this.East_Interval = new List<Interval>();
                        this.East1_Interval = new List<Interval>();
                    }
                    else if (this.East_Interval == null)
                    {
                        this.West_Interval = new List<Interval>();
                        this.West1_Interval = new List<Interval>();
                    }
                }
                else
                {
                    //  去掉中间两条竖
                    this.HSouth_Interval = new List<Interval>();
                    this.HNouth_Interval = new List<Interval>();
                }

                return this;
            }
        }

        public TrilinearCell GenerateRemovedH(int removeCount, int directionCode)
        {
            if (this.ShapeType != TrilinearShapeType.ZShape
                && this.ShapeType != TrilinearShapeType.EShape
                && this.ShapeType != TrilinearShapeType.EVarientShape)
            {
                return this;
            }
            else
            {
                if (removeCount == 3)
                {
                    this.South_Interval = new List<Interval>();
                    this.Middle_Interval = new List<Interval>();
                    this.North_Interval = new List<Interval>();
                }
                else if (removeCount == 2)
                {
                    if (directionCode == 0)
                    {
                        this.South_Interval = new List<Interval>();
                        this.Middle_Interval = new List<Interval>();
                    }
                    else if (directionCode == 1)
                    {
                        this.Middle_Interval = new List<Interval>();
                        this.North_Interval = new List<Interval>();
                    }
                    else
                    {
                        this.South_Interval = new List<Interval>();
                        this.North_Interval = new List<Interval>();
                    }
                }
                else
                {
                    if (directionCode == 0)
                    {
                        this.South_Interval = new List<Interval>();
                    }
                    else if (directionCode == 1)
                    {
                        this.Middle_Interval = new List<Interval>();
                    }
                    else
                    {
                        this.North_Interval = new List<Interval>();
                    }
                }
            }

            return this;
        }

        public TrilinearCell GenerateScaledShape(double lForHighRise, double scaleFactor,double lMin, double w)
        {
            
            //Curve crv0 = this.CellBoundary.DuplicateCurve();
            //crv0.Transform(transform);

            Polyline poly;
            this.CellBoundary.TryGetPolyline(out poly);
            List<Point3d> pts = poly.ToList();
            pts.RemoveAt(pts.Count - 1);
            int basePointIndex = pts.IndexOf(this.TurningPoint);

            Point3d pt0 = Point3d.Unset;
            Point3d pt1 = Point3d.Unset;
            List<Point3d> ptList0 = new List<Point3d>();
            List<Point3d> ptList1 = new List<Point3d>();
            if (basePointIndex == 0 || basePointIndex == 2)
            {
                pt0 = (pts[((basePointIndex - 1) + pts.Count) % pts.Count] + this.TurningPoint) / 2;
                pt1 = (pts[((basePointIndex - 2) + pts.Count) % pts.Count] + pts[(basePointIndex + 1) % pts.Count]) / 2;

                ptList0.Add(this.TurningPoint);
                ptList0.Add(pts[(basePointIndex + 1) % pts.Count]);
                ptList0.Add(pt1);
                ptList0.Add(pt0);
                ptList0.Add(this.TurningPoint);

                ptList1.Add(pt0);
                ptList1.Add(pt1);
                ptList1.Add(pts[((basePointIndex - 2) + pts.Count) % pts.Count]);
                ptList1.Add(pts[((basePointIndex - 1) + pts.Count) % pts.Count]);
                ptList1.Add(pt0);

                //Curve crv0 = new Polyline(ptList0).ToNurbsCurve();
                ptList0.RemoveAt(ptList0.Count - 1);
                //Curve crv1 = new Polyline(ptList1).ToNurbsCurve();
                ptList1.RemoveAt(ptList1.Count - 1);

                #region 构造两个Brep，并且得到surf
                Curve brep0_railCrv0;
                Curve brep0_railCrv1;
                List<Curve> brep0_shapes = new List<Curve>();

                brep0_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 3) % ptList0.Count] }, 1);
                brep0_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 1) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                Curve brep0_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                Curve brep0_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                brep0_shapes.Add(brep0_shape0);
                brep0_shapes.Add(brep0_shape1);

                Curve brep1_railCrv0;
                Curve brep1_railCrv1;
                List<Curve> brep1_shapes = new List<Curve>();
                brep1_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 3) % ptList1.Count] }, 1);
                brep1_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 1) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                Curve brep1_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                Curve brep1_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                brep1_shapes.Add(brep1_shape0);
                brep1_shapes.Add(brep1_shape1);

                Brep brp0 = Brep.CreateFromSweep(brep0_railCrv0, brep0_railCrv1, brep0_shapes, false, Tolerance)[0];
                Surface surf0 = brp0.Surfaces[0];
                surf0.SetDomain(0, new Interval(0, 1));
                surf0.SetDomain(1, new Interval(0, 1));

                Brep brp1 = Brep.CreateFromSweep(brep1_railCrv0, brep1_railCrv1, brep1_shapes, false, Tolerance)[0];
                Surface surf1 = brp1.Surfaces[0];
                surf1.SetDomain(0, new Interval(0, 1));
                surf1.SetDomain(1, new Interval(0, 1));
                #endregion

                #region 按照scale的值，计算此时XY方向上的占比
                // 各自的占比就是scale
                double scale_X0 = scaleFactor;
                double scale_Y0 = scaleFactor;

                double scale_X1 = scaleFactor;
                double scale_Y1 = scaleFactor;
                #endregion

                #region 求XY方向缩放比的下限
                #region 0
                Line lineX0 = new Line(ptList0[0], ptList0[1]);
                Line lineY0 = new Line(ptList0[0], ptList0[3]);
                double lengthX0 = lineX0.Length;
                double lengthY0 = lineY0.Length;

                double scaleFactor_X0_Min;
                double scaleFactor_Y0_Min;
                if (lengthX0 > w)
                {
                    scaleFactor_X0_Min = w / lengthX0;
                }
                else
                {
                    scaleFactor_X0_Min = 1;
                }
                if (lengthY0 > w)
                {
                    scaleFactor_Y0_Min = w / lengthY0;
                }
                else
                {
                    scaleFactor_Y0_Min = 1;
                }

                double scaleFactor_X0;
                double scaleFactor_Y0;
                if (lengthX0 > lForHighRise)
                {
                    scaleFactor_X0 = lForHighRise / lengthX0;
                }
                else
                {
                    scaleFactor_X0 = 1;
                }
                if (lengthY0 > lForHighRise)
                {
                    scaleFactor_Y0 = lForHighRise / lengthY0;
                }
                else
                {
                    scaleFactor_Y0 = 1;
                }
                #endregion

                #region 1
                Line lineX1 = new Line(ptList1[2], ptList1[3]);
                Line lineY1 = new Line(ptList1[2], ptList1[1]);
                double lengthX1 = lineX1.Length;
                double lengthY1 = lineY1.Length;

                double scaleFactor_X1_Min;
                double scaleFactor_Y1_Min;
                if (lengthX1 > w)
                {
                    scaleFactor_X1_Min = w / lengthX1;
                }
                else
                {
                    scaleFactor_X1_Min = 1;
                }
                if (lengthY1 > w)
                {
                    scaleFactor_Y1_Min = w / lengthY1;
                }
                else
                {
                    scaleFactor_Y1_Min = 1;
                }

                double scaleFactor_X1;
                double scaleFactor_Y1;
                if (lengthX1 > lForHighRise)
                {
                    scaleFactor_X1 = lForHighRise / lengthX1;
                }
                else
                {
                    scaleFactor_X1 = 1;
                }
                if (lengthY1 > lForHighRise)
                {
                    scaleFactor_Y1 = lForHighRise / lengthY1;
                }
                else
                {
                    scaleFactor_Y1 = 1;
                }
                #endregion
                #endregion

                #region 比较
                #region 0
                if (scale_X0 < scaleFactor_X0)
                {
                    scale_X0 = scaleFactor_X0;
                }
                if (scale_Y0 < scaleFactor_Y0)
                {
                    scale_Y0 = scaleFactor_Y0;
                }
                if (scale_X0 < scaleFactor_X0_Min)
                {
                    scale_X0 = scaleFactor_X0_Min;
                }
                if (scale_Y0 < scaleFactor_Y0_Min)
                {
                    scale_Y0 = scaleFactor_Y0_Min;
                }
                #endregion

                #region 1
                if (scale_X1 < scaleFactor_X1)
                {
                    scale_X1 = scaleFactor_X1;
                }
                if (scale_Y1 < scaleFactor_Y1)
                {
                    scale_Y1 = scaleFactor_Y1;
                }
                if (scale_X1 < scaleFactor_X1_Min)
                {
                    scale_X1 = scaleFactor_X1_Min;
                }
                if (scale_Y1 < scaleFactor_Y1_Min)
                {
                    scale_Y1 = scaleFactor_Y1_Min;
                }
                #endregion
                #endregion

                #region 取子面
                #region 0
                UVInterval uvinterval0;
                if (basePointIndex == 0)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_Y0), new Interval(0, scale_X0));
                }
                else if (basePointIndex == 1)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_X0), new Interval((1 - scale_Y0), 1));
                }
                else if (basePointIndex == 2)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_Y0), new Interval(0, scale_X0));
                }
                else
                {
                    uvinterval0 = new UVInterval(new Interval((1 - scale_X0), 1), new Interval(0, scale_Y0));
                }
                Surface trimedSurf0 = UtilityFunctions.IsoTrim(surf0, uvinterval0);
                Curve[] crvs0 = trimedSurf0.ToBrep().DuplicateEdgeCurves();
                Curve crv0 = Curve.JoinCurves(crvs0)[0];
                #endregion

                #region 1
                UVInterval uvinterval1;
                if (basePointIndex == 0)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_Y1), new Interval(0, scale_X1));
                }
                else if (basePointIndex == 1)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_X1), new Interval((1 - scale_Y1), 1));
                }
                else if (basePointIndex == 2)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_Y1), new Interval(0, scale_X1));
                }
                else
                {
                    uvinterval1 = new UVInterval(new Interval((1 - scale_X1), 1), new Interval(0, scale_Y1));
                }
                Surface trimedSurf1 = UtilityFunctions.IsoTrim(surf1, uvinterval1);
                Curve[] crvs1 = trimedSurf1.ToBrep().DuplicateEdgeCurves();
                Curve crv1 = Curve.JoinCurves(crvs1)[0];
                #endregion
                #endregion

                this.crvForScaled = new Curve[2] { crv0, crv1 };
            }
            else if (basePointIndex == 1 || basePointIndex == 3)
            {
                pt1 = (pts[(basePointIndex + 1) % pts.Count] + this.TurningPoint) / 2;
                pt0 = (pts[(basePointIndex + 2) % pts.Count] + pts[(basePointIndex + 3) % pts.Count]) / 2;

                ptList0.Add(pts[(basePointIndex + 3) % pts.Count]);
                ptList0.Add(this.TurningPoint);
                ptList0.Add(pt1);
                ptList0.Add(pt0);
                ptList0.Add(pts[(basePointIndex + 3) % pts.Count]);

                ptList1.Add(pt0);
                ptList1.Add(pt1);
                ptList1.Add(pts[(basePointIndex + 1) % pts.Count]);
                ptList1.Add(pts[(basePointIndex + 2) % pts.Count]);
                ptList1.Add(pt0);

                //Curve crv0 = new Polyline(ptList0).ToNurbsCurve();
                ptList0.RemoveAt(ptList0.Count - 1);
                //Curve crv1 = new Polyline(ptList1).ToNurbsCurve();
                ptList1.RemoveAt(ptList1.Count - 1);

                #region 构造两个Brep，并且得到surf
                Curve brep0_railCrv0;
                Curve brep0_railCrv1;
                List<Curve> brep0_shapes = new List<Curve>();

                brep0_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                brep0_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                Curve brep0_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[basePointIndex] }, 1);
                Curve brep0_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 2) % ptList0.Count], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                brep0_shapes.Add(brep0_shape0);
                brep0_shapes.Add(brep0_shape1);

                Curve brep1_railCrv0;
                Curve brep1_railCrv1;
                List<Curve> brep1_shapes = new List<Curve>();

                brep1_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                brep1_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                Curve brep1_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[basePointIndex] }, 1);
                Curve brep1_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 2) % ptList1.Count], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                brep1_shapes.Add(brep1_shape0);
                brep1_shapes.Add(brep1_shape1);

                Brep brp0 = Brep.CreateFromSweep(brep0_railCrv0, brep0_railCrv1, brep0_shapes, false, Tolerance)[0];
                Surface surf0 = brp0.Surfaces[0];
                surf0.SetDomain(0, new Interval(0, 1));
                surf0.SetDomain(1, new Interval(0, 1));

                Brep brp1 = Brep.CreateFromSweep(brep1_railCrv0, brep1_railCrv1, brep1_shapes, false, Tolerance)[0];
                Surface surf1 = brp1.Surfaces[0];
                surf1.SetDomain(0, new Interval(0, 1));
                surf1.SetDomain(1, new Interval(0, 1));
                #endregion

                #region 按照scale的值，计算此时XY方向上的占比
                // 各自的占比就是scale
                double scale_X0 = scaleFactor;
                double scale_Y0 = scaleFactor;

                double scale_X1 = scaleFactor;
                double scale_Y1 = scaleFactor;
                #endregion

                #region 求XY方向缩放比的下限
                #region 0
                Line lineX0 = new Line(ptList0[1], ptList0[2]);
                Line lineY0 = new Line(ptList0[1], ptList0[0]);
                double lengthX0 = lineX0.Length;
                double lengthY0 = lineY0.Length;

                double scaleFactor_X0_Min;
                double scaleFactor_Y0_Min;
                if (lengthX0 > w)
                {
                    scaleFactor_X0_Min = w / lengthX0;
                }
                else
                {
                    scaleFactor_X0_Min = 1;
                }
                if (lengthY0 > w)
                {
                    scaleFactor_Y0_Min = w / lengthY0;
                }
                else
                {
                    scaleFactor_Y0_Min = 1;
                }

                double scaleFactor_X0;
                double scaleFactor_Y0;
                if (lengthX0 > lForHighRise)
                {
                    scaleFactor_X0 = lForHighRise / lengthX0;
                }
                else
                {
                    scaleFactor_X0 = 1;
                }
                if (lengthY0 > lForHighRise)
                {
                    scaleFactor_Y0 = lForHighRise / lengthY0;
                }
                else
                {
                    scaleFactor_Y0 = 1;
                }
                #endregion

                #region 1
                Line lineX1 = new Line(ptList1[3], ptList1[0]);
                Line lineY1 = new Line(ptList1[3], ptList1[2]);
                double lengthX1 = lineX1.Length;
                double lengthY1 = lineY1.Length;

                double scaleFactor_X1_Min;
                double scaleFactor_Y1_Min;
                if (lengthX1 > w)
                {
                    scaleFactor_X1_Min = w / lengthX1;
                }
                else
                {
                    scaleFactor_X1_Min = 1;
                }
                if (lengthY1 > w)
                {
                    scaleFactor_Y1_Min = w / lengthY1;
                }
                else
                {
                    scaleFactor_Y1_Min = 1;
                }

                double scaleFactor_X1;
                double scaleFactor_Y1;
                if (lengthX1 > lForHighRise)
                {
                    scaleFactor_X1 = lForHighRise / lengthX1;
                }
                else
                {
                    scaleFactor_X1 = 1;
                }
                if (lengthY1 > lForHighRise)
                {
                    scaleFactor_Y1 = lForHighRise / lengthY1;
                }
                else
                {
                    scaleFactor_Y1 = 1;
                }
                #endregion
                #endregion

                #region 比较
                #region 0
                if (scale_X0 < scaleFactor_X0)
                {
                    scale_X0 = scaleFactor_X0;
                }
                if (scale_Y0 < scaleFactor_Y0)
                {
                    scale_Y0 = scaleFactor_Y0;
                }
                if (scale_X0 < scaleFactor_X0_Min)
                {
                    scale_X0 = scaleFactor_X0_Min;
                }
                if (scale_Y0 < scaleFactor_Y0_Min)
                {
                    scale_Y0 = scaleFactor_Y0_Min;
                }
                #endregion

                #region 1
                if (scale_X1 < scaleFactor_X1)
                {
                    scale_X1 = scaleFactor_X1;
                }
                if (scale_Y1 < scaleFactor_Y1)
                {
                    scale_Y1 = scaleFactor_Y1;
                }
                if (scale_X1 < scaleFactor_X1_Min)
                {
                    scale_X1 = scaleFactor_X1_Min;
                }
                if (scale_Y1 < scaleFactor_Y1_Min)
                {
                    scale_Y1 = scaleFactor_Y1_Min;
                }
                #endregion
                #endregion

                #region 取子面
                #region 0
                UVInterval uvinterval0;
                if (basePointIndex == 0)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_Y0), new Interval(0, scale_X0));
                }
                else if (basePointIndex == 1)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_X0), new Interval((1 - scale_Y0), 1));
                }
                else if (basePointIndex == 2)
                {
                    uvinterval0 = new UVInterval(new Interval(0, scale_Y0), new Interval(0, scale_X0));
                }
                else
                {
                    uvinterval0 = new UVInterval(new Interval((1 - scale_X0), 1), new Interval(0, scale_Y0));
                }
                Surface trimedSurf0 = UtilityFunctions.IsoTrim(surf0, uvinterval0);
                Curve[] crvs0 = trimedSurf0.ToBrep().DuplicateEdgeCurves();
                Curve crv0 = Curve.JoinCurves(crvs0)[0];
                #endregion

                #region 1
                UVInterval uvinterval1;
                if (basePointIndex == 0)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_Y1), new Interval(0, scale_X1));
                }
                else if (basePointIndex == 1)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_X1), new Interval((1 - scale_Y1), 1));
                }
                else if (basePointIndex == 2)
                {
                    uvinterval1 = new UVInterval(new Interval(0, scale_Y1), new Interval(0, scale_X1));
                }
                else
                {
                    uvinterval1 = new UVInterval(new Interval((1 - scale_X1), 1), new Interval(0, scale_Y1));
                }
                Surface trimedSurf1 = UtilityFunctions.IsoTrim(surf1, uvinterval1);
                Curve[] crvs1 = trimedSurf1.ToBrep().DuplicateEdgeCurves();
                Curve crv1 = Curve.JoinCurves(crvs1)[0];
                #endregion
                #endregion

                this.crvForScaled = new Curve[2] { crv0, crv1 };
            }

            #region 清空所有的interval
            this.South_Interval = null;
            this.North_Interval = null;
            this.East_Interval = null;
            this.East1_Interval = null;
            this.West_Interval = null;
            this.West1_Interval = null;
            this.Middle_Interval = null;
            this.HSouth_Interval = null;
            this.HNouth_Interval = null;
            #endregion
            this.ShapeType = TrilinearShapeType.Scaled;

            return this;
        }

        public TrilinearCell GenerateHighRise(double lForHighRise, double lMin, double w, int randomCode0, Random m_random)
        {
            Polyline poly;
            this.CellBoundary.TryGetPolyline(out poly);
            List<Point3d> pts = poly.ToList();
            pts.RemoveAt(pts.Count - 1);

            if (this.ShapeType == TrilinearShapeType.Scaled)
            {
                // 找turningPoint
                int basePointIndex = pts.IndexOf(this.TurningPoint);

                Point3d pt0 = Point3d.Unset;
                Point3d pt1 = Point3d.Unset;
                List<Point3d> ptList0 = new List<Point3d>();
                List<Point3d> ptList1 = new List<Point3d>();
                Point3d basePoint = Point3d.Unset;
                if (basePointIndex == 0 || basePointIndex == 2)
                {
                    pt0 = (pts[((basePointIndex - 1) + pts.Count) % pts.Count] + this.TurningPoint) / 2;
                    pt1 = (pts[((basePointIndex - 2) + pts.Count) % pts.Count] + pts[(basePointIndex + 1) % pts.Count]) / 2;

                    ptList0.Add(this.TurningPoint);
                    ptList0.Add(pts[(basePointIndex + 1) % pts.Count]);
                    ptList0.Add(pt1);
                    ptList0.Add(pt0);
                    ptList0.Add(this.TurningPoint);

                    ptList1.Add(pt0);
                    ptList1.Add(pt1);
                    ptList1.Add(pts[((basePointIndex - 2) + pts.Count) % pts.Count]);
                    ptList1.Add(pts[((basePointIndex - 1) + pts.Count) % pts.Count]);
                    ptList1.Add(pt0);

                    if (randomCode0 < 2)
                    {
                        basePoint = ptList0[0];
                        ptList0.RemoveAt(ptList0.Count - 1);

                        #region 构造第0个Brep，并且得到surf
                        Curve brep0_railCrv0;
                        Curve brep0_railCrv1;
                        List<Curve> brep0_shapes = new List<Curve>();

                        brep0_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 3) % ptList0.Count] }, 1);
                        brep0_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 1) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                        Curve brep0_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                        Curve brep0_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                        brep0_shapes.Add(brep0_shape0);
                        brep0_shapes.Add(brep0_shape1);

                        Brep brp0 = Brep.CreateFromSweep(brep0_railCrv0, brep0_railCrv1, brep0_shapes, false, Tolerance)[0];
                        Surface surf0 = brp0.Surfaces[0];
                        surf0.SetDomain(0, new Interval(0, 1));
                        surf0.SetDomain(1, new Interval(0, 1));
                        #endregion

                        #region 求XY方向缩放比的下限
                        #region 0
                        Line lineX0 = new Line(ptList0[0], ptList0[1]);
                        Line lineY0 = new Line(ptList0[0], ptList0[3]);
                        double lengthX0 = lineX0.Length;
                        double lengthY0 = lineY0.Length;

                        double scaleFactor_X0_Min;
                        double scaleFactor_Y0_Min;
                        if (lengthX0 > w)
                        {
                            scaleFactor_X0_Min = w / lengthX0;
                        }
                        else
                        {
                            scaleFactor_X0_Min = 1;
                        }
                        if (lengthY0 > w)
                        {
                            scaleFactor_Y0_Min = w / lengthY0;
                        }
                        else
                        {
                            scaleFactor_Y0_Min = 1;
                        }

                        double scaleFactor_X0;
                        double scaleFactor_Y0;
                        if (lengthX0 > lForHighRise)
                        {
                            scaleFactor_X0 = lForHighRise / lengthX0;
                        }
                        else
                        {
                            scaleFactor_X0 = 1;
                        }
                        if (lengthY0 > lForHighRise)
                        {
                            scaleFactor_Y0 = lForHighRise / lengthY0;
                        }
                        else
                        {
                            scaleFactor_Y0 = 1;
                        }
                        #endregion
                        #endregion

                        #region XY方向缩放比与下限比较
                        #region 0
                        if (scaleFactor_X0 < scaleFactor_X0_Min)
                        {
                            scaleFactor_X0 = scaleFactor_X0_Min;
                        }
                        if (scaleFactor_Y0 < scaleFactor_Y0_Min)
                        {
                            scaleFactor_Y0 = scaleFactor_Y0_Min;
                        }
                        #endregion
                        #endregion

                        #region 取子面
                        #region 0
                        UVInterval uvinterval0;
                        if (basePointIndex == 0)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                        }
                        else if (basePointIndex == 1)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_X0), new Interval((1 - scaleFactor_Y0), 1));
                        }
                        else if (basePointIndex == 2)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                        }
                        else
                        {
                            uvinterval0 = new UVInterval(new Interval((1 - scaleFactor_X0), 1), new Interval(0, scaleFactor_Y0));
                        }
                        Surface trimedSurf0 = UtilityFunctions.IsoTrim(surf0, uvinterval0);
                        Curve[] crvs0 = trimedSurf0.ToBrep().DuplicateEdgeCurves();
                        Curve crv0 = Curve.JoinCurves(crvs0)[0];
                        #endregion
                        #endregion

                        this.crvForHighRise = crv0.DuplicateCurve();
                    }
                    else
                    {
                        basePoint = ptList1[2];
                        ptList1.RemoveAt(ptList1.Count - 1);

                        #region 构造第1个Brep，并且得到surf
                        Curve brep1_railCrv0;
                        Curve brep1_railCrv1;
                        List<Curve> brep1_shapes = new List<Curve>();

                        brep1_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 3) % ptList1.Count] }, 1);
                        brep1_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 1) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                        Curve brep1_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                        Curve brep1_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                        brep1_shapes.Add(brep1_shape0);
                        brep1_shapes.Add(brep1_shape1);

                        Brep brp1 = Brep.CreateFromSweep(brep1_railCrv0, brep1_railCrv1, brep1_shapes, false, Tolerance)[0];
                        Surface surf1 = brp1.Surfaces[0];
                        surf1.SetDomain(0, new Interval(0, 1));
                        surf1.SetDomain(1, new Interval(0, 1));
                        #endregion

                        #region 求XY方向缩放比的下限
                        #region 1
                        Line lineX1 = new Line(ptList1[2], ptList1[3]);
                        Line lineY1 = new Line(ptList1[2], ptList1[1]);
                        double lengthX1 = lineX1.Length;
                        double lengthY1 = lineY1.Length;

                        double scaleFactor_X1_Min;
                        double scaleFactor_Y1_Min;
                        if (lengthX1 > w)
                        {
                            scaleFactor_X1_Min = w / lengthX1;
                        }
                        else
                        {
                            scaleFactor_X1_Min = 1;
                        }
                        if (lengthY1 > w)
                        {
                            scaleFactor_Y1_Min = w / lengthY1;
                        }
                        else
                        {
                            scaleFactor_Y1_Min = 1;
                        }

                        double scaleFactor_X1;
                        double scaleFactor_Y1;
                        if (lengthX1 > lForHighRise)
                        {
                            scaleFactor_X1 = lForHighRise / lengthX1;
                        }
                        else
                        {
                            scaleFactor_X1 = 1;
                        }
                        if (lengthY1 > lForHighRise)
                        {
                            scaleFactor_Y1 = lForHighRise / lengthY1;
                        }
                        else
                        {
                            scaleFactor_Y1 = 1;
                        }
                        #endregion
                        #endregion

                        #region XY方向缩放比与下限比较
                        #region 1
                        if (scaleFactor_X1 < scaleFactor_X1_Min)
                        {
                            scaleFactor_X1 = scaleFactor_X1_Min;
                        }
                        if (scaleFactor_Y1 < scaleFactor_Y1_Min)
                        {
                            scaleFactor_Y1 = scaleFactor_Y1_Min;
                        }
                        #endregion
                        #endregion

                        #region 取子面
                        #region 1
                        UVInterval uvinterval1;
                        if (basePointIndex == 0)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                        }
                        else if (basePointIndex == 1)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_X1), new Interval((1 - scaleFactor_Y1), 1));
                        }
                        else if (basePointIndex == 2)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                        }
                        else
                        {
                            uvinterval1 = new UVInterval(new Interval((1 - scaleFactor_X1), 1), new Interval(0, scaleFactor_Y1));
                        }
                        Surface trimedSurf1 = UtilityFunctions.IsoTrim(surf1, uvinterval1);
                        Curve[] crvs1 = trimedSurf1.ToBrep().DuplicateEdgeCurves();
                        Curve crv1 = Curve.JoinCurves(crvs1)[0];
                        #endregion
                        #endregion

                        this.crvForHighRise = crv1.DuplicateCurve();
                    }
                }
                else if (basePointIndex == 1 || basePointIndex == 3)
                {
                    pt1 = (pts[(basePointIndex + 1) % pts.Count] + this.TurningPoint) / 2;
                    pt0 = (pts[(basePointIndex + 2) % pts.Count] + pts[(basePointIndex + 3) % pts.Count]) / 2;

                    ptList0.Add(pts[(basePointIndex + 3) % pts.Count]);
                    ptList0.Add(this.TurningPoint);
                    ptList0.Add(pt1);
                    ptList0.Add(pt0);
                    ptList0.Add(pts[(basePointIndex + 3) % pts.Count]);

                    ptList1.Add(pt0);
                    ptList1.Add(pt1);
                    ptList1.Add(pts[(basePointIndex + 1) % pts.Count]);
                    ptList1.Add(pts[(basePointIndex + 2) % pts.Count]);
                    ptList1.Add(pt0);

                    if (randomCode0 < 2)
                    {
                        basePoint = ptList0[1];
                        ptList0.RemoveAt(ptList0.Count - 1);

                        #region 构造第0个Brep，并且得到surf
                        Curve brep0_railCrv0;
                        Curve brep0_railCrv1;
                        List<Curve> brep0_shapes = new List<Curve>();

                        brep0_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[(basePointIndex + 2) % ptList0.Count] }, 1);
                        brep0_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[basePointIndex], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                        Curve brep0_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 3) % ptList0.Count], ptList0[basePointIndex] }, 1);
                        Curve brep0_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[(basePointIndex + 2) % ptList0.Count], ptList0[(basePointIndex + 1) % ptList0.Count] }, 1);
                        brep0_shapes.Add(brep0_shape0);
                        brep0_shapes.Add(brep0_shape1);

                        Brep brp0 = Brep.CreateFromSweep(brep0_railCrv0, brep0_railCrv1, brep0_shapes, false, Tolerance)[0];
                        Surface surf0 = brp0.Surfaces[0];
                        surf0.SetDomain(0, new Interval(0, 1));
                        surf0.SetDomain(1, new Interval(0, 1));
                        #endregion

                        #region 求XY方向缩放比的下限
                        #region 0
                        Line lineX0 = new Line(ptList0[0], ptList0[1]);
                        Line lineY0 = new Line(ptList0[0], ptList0[3]);
                        double lengthX0 = lineX0.Length;
                        double lengthY0 = lineY0.Length;

                        double scaleFactor_X0_Min;
                        double scaleFactor_Y0_Min;
                        if (lengthX0 > w)
                        {
                            scaleFactor_X0_Min = w / lengthX0;
                        }
                        else
                        {
                            scaleFactor_X0_Min = 1;
                        }
                        if (lengthY0 > w)
                        {
                            scaleFactor_Y0_Min = w / lengthY0;
                        }
                        else
                        {
                            scaleFactor_Y0_Min = 1;
                        }

                        double scaleFactor_X0;
                        double scaleFactor_Y0;
                        if (lengthX0 > lForHighRise)
                        {
                            scaleFactor_X0 = lForHighRise / lengthX0;
                        }
                        else
                        {
                            scaleFactor_X0 = 1;
                        }
                        if (lengthY0 > lForHighRise)
                        {
                            scaleFactor_Y0 = lForHighRise / lengthY0;
                        }
                        else
                        {
                            scaleFactor_Y0 = 1;
                        }
                        #endregion
                        #endregion

                        #region XY方向缩放比与下限比较
                        #region 0
                        if (scaleFactor_X0 < scaleFactor_X0_Min)
                        {
                            scaleFactor_X0 = scaleFactor_X0_Min;
                        }
                        if (scaleFactor_Y0 < scaleFactor_Y0_Min)
                        {
                            scaleFactor_Y0 = scaleFactor_Y0_Min;
                        }
                        #endregion
                        #endregion

                        #region 取子面
                        #region 0
                        UVInterval uvinterval0;
                        if (basePointIndex == 0)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                        }
                        else if (basePointIndex == 1)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_X0), new Interval((1 - scaleFactor_Y0), 1));
                        }
                        else if (basePointIndex == 2)
                        {
                            uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                        }
                        else
                        {
                            uvinterval0 = new UVInterval(new Interval((1 - scaleFactor_X0), 1), new Interval(0, scaleFactor_Y0));
                        }
                        Surface trimedSurf0 = UtilityFunctions.IsoTrim(surf0, uvinterval0);
                        Curve[] crvs0 = trimedSurf0.ToBrep().DuplicateEdgeCurves();
                        Curve crv0 = Curve.JoinCurves(crvs0)[0];
                        #endregion
                        #endregion

                        this.crvForHighRise = crv0.DuplicateCurve();
                    }
                    else
                    {
                        basePoint = ptList1[3];
                        ptList1.RemoveAt(ptList1.Count - 1);

                        #region 构造第1个Brep，并且得到surf
                        Curve brep1_railCrv0;
                        Curve brep1_railCrv1;
                        List<Curve> brep1_shapes = new List<Curve>();

                        brep1_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[(basePointIndex + 2) % ptList1.Count] }, 1);
                        brep1_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[basePointIndex], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                        Curve brep1_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 3) % ptList1.Count], ptList1[basePointIndex] }, 1);
                        Curve brep1_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[(basePointIndex + 2) % ptList1.Count], ptList1[(basePointIndex + 1) % ptList1.Count] }, 1);
                        brep1_shapes.Add(brep1_shape0);
                        brep1_shapes.Add(brep1_shape1);

                        Brep brp1 = Brep.CreateFromSweep(brep1_railCrv0, brep1_railCrv1, brep1_shapes, false, Tolerance)[0];
                        Surface surf1 = brp1.Surfaces[0];
                        surf1.SetDomain(0, new Interval(0, 1));
                        surf1.SetDomain(1, new Interval(0, 1));
                        #endregion

                        #region 求XY方向缩放比的下限
                        #region 1
                        Line lineX1 = new Line(ptList1[2], ptList1[3]);
                        Line lineY1 = new Line(ptList1[2], ptList1[1]);
                        double lengthX1 = lineX1.Length;
                        double lengthY1 = lineY1.Length;

                        double scaleFactor_X1_Min;
                        double scaleFactor_Y1_Min;
                        if (lengthX1 > w)
                        {
                            scaleFactor_X1_Min = w / lengthX1;
                        }
                        else
                        {
                            scaleFactor_X1_Min = 1;
                        }
                        if (lengthY1 > w)
                        {
                            scaleFactor_Y1_Min = w / lengthY1;
                        }
                        else
                        {
                            scaleFactor_Y1_Min = 1;
                        }

                        double scaleFactor_X1;
                        double scaleFactor_Y1;
                        if (lengthX1 > lForHighRise)
                        {
                            scaleFactor_X1 = lForHighRise / lengthX1;
                        }
                        else
                        {
                            scaleFactor_X1 = 1;
                        }
                        if (lengthY1 > lForHighRise)
                        {
                            scaleFactor_Y1 = lForHighRise / lengthY1;
                        }
                        else
                        {
                            scaleFactor_Y1 = 1;
                        }
                        #endregion
                        #endregion

                        #region XY方向缩放比与下限比较
                        #region 1
                        if (scaleFactor_X1 < scaleFactor_X1_Min)
                        {
                            scaleFactor_X1 = scaleFactor_X1_Min;
                        }
                        if (scaleFactor_Y1 < scaleFactor_Y1_Min)
                        {
                            scaleFactor_Y1 = scaleFactor_Y1_Min;
                        }
                        #endregion
                        #endregion

                        #region 取子面
                        #region 1
                        UVInterval uvinterval1;
                        if (basePointIndex == 0)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                        }
                        else if (basePointIndex == 1)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_X1), new Interval((1 - scaleFactor_Y1), 1));
                        }
                        else if (basePointIndex == 2)
                        {
                            uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                        }
                        else
                        {
                            uvinterval1 = new UVInterval(new Interval((1 - scaleFactor_X1), 1), new Interval(0, scaleFactor_Y1));
                        }
                        Surface trimedSurf1 = UtilityFunctions.IsoTrim(surf1, uvinterval1);
                        Curve[] crvs1 = trimedSurf1.ToBrep().DuplicateEdgeCurves();
                        Curve crv1 = Curve.JoinCurves(crvs1)[0];
                        #endregion
                        #endregion

                        this.crvForHighRise = crv1.DuplicateCurve();
                    }
                }
            }
            else
            {
                // 随机
                int turningPointIndex = pts.IndexOf(this.TurningPoint);

                Point3d pt0 = Point3d.Unset;
                Point3d pt1 = Point3d.Unset;
                List<Point3d> ptList0 = new List<Point3d>();
                List<Point3d> ptList1 = new List<Point3d>();
                if (turningPointIndex == 0 || turningPointIndex == 2)
                {
                    pt0 = (pts[((turningPointIndex - 1) + pts.Count) % pts.Count] + this.TurningPoint) / 2;
                    pt1 = (pts[((turningPointIndex - 2) + pts.Count) % pts.Count] + pts[(turningPointIndex + 1) % pts.Count]) / 2;

                    ptList0.Add(this.TurningPoint);
                    ptList0.Add(pts[(turningPointIndex + 1) % pts.Count]);
                    ptList0.Add(pt1);
                    ptList0.Add(pt0);
                    ptList0.Add(this.TurningPoint);

                    ptList1.Add(pt0);
                    ptList1.Add(pt1);
                    ptList1.Add(pts[((turningPointIndex - 2) + pts.Count) % pts.Count]);
                    ptList1.Add(pts[((turningPointIndex - 1) + pts.Count) % pts.Count]);
                    ptList1.Add(pt0);
                }
                else if (turningPointIndex == 1 || turningPointIndex == 3)
                {
                    pt1 = (pts[(turningPointIndex + 1) % pts.Count] + this.TurningPoint) / 2;
                    pt0 = (pts[(turningPointIndex + 2) % pts.Count] + pts[(turningPointIndex + 3) % pts.Count]) / 2;

                    ptList0.Add(pts[(turningPointIndex + 3) % pts.Count]);
                    ptList0.Add(this.TurningPoint);
                    ptList0.Add(pt1);
                    ptList0.Add(pt0);
                    ptList0.Add(pts[(turningPointIndex + 3) % pts.Count]);

                    ptList1.Add(pt0);
                    ptList1.Add(pt1);
                    ptList1.Add(pts[(turningPointIndex + 1) % pts.Count]);
                    ptList1.Add(pts[(turningPointIndex + 2) % pts.Count]);
                    ptList1.Add(pt0);
                }

                Point3d basePoint = Point3d.Unset;
                if (randomCode0 < 2)
                {
                    if (turningPointIndex == 0 || turningPointIndex == 1)
                    {
                        // 下 ptList0
                        
                        // 找到可以放置高层的位置
                        //List<Point3d> pointsForCanPlaceHighrise_Two = new List<Point3d>();
                        //List<Point3d> pointsForCanPlaceHighrise_One = new List<Point3d>();

                        List<int> indexForCanPlaceHighrise_Two = new List<int>();
                        List<int> indexForCanPlaceHighrise_One = new List<int>();

                        #region west
                        bool west_0 = false;
                        bool west_1 = false;

                        if (this.West_Interval != null)
                        {
                            for (int i = 0; i < this.West_Interval.Count; i++)
                            {
                                if (this.West_Interval[i].IncludesParameter(0))
                                {
                                    west_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.West_Interval.Count; i++)
                            {
                                if (this.West_Interval[i].IncludesParameter(1))
                                {
                                    west_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region south
                        bool south_0 = false;
                        bool south_1 = false;

                        if (this.South_Interval != null)
                        {
                            for (int i = 0; i < this.South_Interval.Count; i++)
                            {
                                if (this.South_Interval[i].IncludesParameter(0))
                                {
                                    south_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.South_Interval.Count; i++)
                            {
                                if (this.South_Interval[i].IncludesParameter(1))
                                {
                                    south_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region east
                        bool east_0 = false;
                        bool east_1 = false;

                        if (this.East_Interval != null)
                        {
                            for (int i = 0; i < this.East_Interval.Count; i++)
                            {
                                if (this.East_Interval[i].IncludesParameter(0))
                                {
                                    east_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.East_Interval.Count; i++)
                            {
                                if (this.East_Interval[i].IncludesParameter(1))
                                {
                                    east_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region Mid
                        bool mid_0 = false;
                        bool mid_1 = false;

                        if (this.Middle_Interval != null)
                        {
                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(0))
                                {
                                    mid_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(1))
                                {
                                    mid_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region 西南角点
                        if (west_0 && south_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(0);
                        }
                        else if (west_0 || south_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(0);
                        }
                        #endregion

                        #region 东南角点
                        if (south_1 && east_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(1);
                        }
                        else if (south_1 || east_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(1);
                        }
                        #endregion

                        #region 东北角点
                        if (east_1 && mid_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(2);
                        }
                        else if (east_1 || mid_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(2);
                        }
                        #endregion

                        #region 西北角点
                        if (mid_0 && west_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(3);
                        }
                        else if (mid_0 || west_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(3);
                        }
                        #endregion

                        #region 随机取点
                        if (indexForCanPlaceHighrise_Two.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_Two.Count);
                            basePoint = ptList0[indexForCanPlaceHighrise_Two[randomIndex]];
                        }
                        else if (indexForCanPlaceHighrise_One.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_One.Count);
                            basePoint = ptList0[indexForCanPlaceHighrise_One[randomIndex]];
                        }
                        else
                        {
                            // 应该不存在此种情况
                            int randomIndex = m_random.Next(4);
                            basePoint = ptList0[randomIndex];
                        }
                        #endregion
                    }
                    else
                    {
                        // turningPointIndex == 2 || turningPointIndex == 3
                        // 上 ptList0

                        // 找到可以放置高层的位置
                        //List<Point3d> pointsForCanPlaceHighrise_Two = new List<Point3d>();
                        //List<Point3d> pointsForCanPlaceHighrise_One = new List<Point3d>();

                        List<int> indexForCanPlaceHighrise_Two = new List<int>();
                        List<int> indexForCanPlaceHighrise_One = new List<int>();

                        #region west
                        bool west_0 = false;
                        bool west_1 = false;

                        if (this.West1_Interval != null)
                        {
                            for (int i = 0; i < this.West1_Interval.Count; i++)
                            {
                                if (this.West1_Interval[i].IncludesParameter(0))
                                {
                                    west_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.West1_Interval.Count; i++)
                            {
                                if (this.West1_Interval[i].IncludesParameter(1))
                                {
                                    west_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region Mid
                        bool mid_0 = false;
                        bool mid_1 = false;

                        if (this.Middle_Interval != null)
                        {
                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(0))
                                {
                                    mid_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(1))
                                {
                                    mid_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region east
                        bool east_0 = false;
                        bool east_1 = false;

                        if (this.East1_Interval != null)
                        {
                            for (int i = 0; i < this.East1_Interval.Count; i++)
                            {
                                if (this.East1_Interval[i].IncludesParameter(0))
                                {
                                    east_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.East1_Interval.Count; i++)
                            {
                                if (this.East1_Interval[i].IncludesParameter(1))
                                {
                                    east_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region north
                        bool north_0 = false;
                        bool north_1 = false;

                        if (this.North_Interval != null)
                        {
                            for (int i = 0; i < this.North_Interval.Count; i++)
                            {
                                if (this.North_Interval[i].IncludesParameter(0))
                                {
                                    north_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.North_Interval.Count; i++)
                            {
                                if (this.North_Interval[i].IncludesParameter(1))
                                {
                                    north_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region 西南角点
                        if (west_0 && mid_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(0);
                        }
                        else if (west_0 || mid_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(0);
                        }
                        #endregion

                        #region 东南角点
                        if (mid_1 && east_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(1);
                        }
                        else if (mid_1 || east_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(1);
                        }
                        #endregion

                        #region 东北角点
                        if (east_1 && north_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(2);
                        }
                        else if (east_1 || north_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(2);
                        }
                        #endregion

                        #region 西北角点
                        if (north_0 && west_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(3);
                        }
                        else if (north_0 || west_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(3);
                        }
                        #endregion

                        #region 随机取点
                        if (indexForCanPlaceHighrise_Two.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_Two.Count);
                            basePoint = ptList0[indexForCanPlaceHighrise_Two[randomIndex]];
                        }
                        else if (indexForCanPlaceHighrise_One.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_One.Count);
                            basePoint = ptList0[indexForCanPlaceHighrise_One[randomIndex]];
                        }
                        else
                        {
                            // 应该不存在此种情况
                            int randomIndex = m_random.Next(4);
                            basePoint = ptList0[randomIndex];
                        }
                        #endregion
                    }

                    //if (randomCode1 == 0)
                    //{
                    //    basePoint = ptList0[0];
                    //}
                    //else if (randomCode1 == 1)
                    //{
                    //    basePoint = ptList0[1];
                    //}
                    //else if (randomCode1 == 2)
                    //{
                    //    basePoint = ptList0[2];
                    //}
                    //else
                    //{
                    //    basePoint = ptList0[3];
                    //}

                    ptList0.RemoveAt(ptList0.Count - 1);
                    int basePointIndex = ptList0.IndexOf(basePoint);

                    #region 构造第0个Brep，并且得到surf
                    Curve brep0_railCrv0;
                    Curve brep0_railCrv1;
                    List<Curve> brep0_shapes = new List<Curve>();

                    brep0_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[0], ptList0[3] }, 1);
                    brep0_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[1], ptList0[2] }, 1);
                    Curve brep0_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[0], ptList0[1] }, 1);
                    Curve brep0_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList0[3], ptList0[2] }, 1);
                    brep0_shapes.Add(brep0_shape0);
                    brep0_shapes.Add(brep0_shape1);

                    Brep brp0 = Brep.CreateFromSweep(brep0_railCrv0, brep0_railCrv1, brep0_shapes, false, Tolerance)[0];
                    Surface surf0 = brp0.Surfaces[0];
                    surf0.SetDomain(0, new Interval(0, 1));
                    surf0.SetDomain(1, new Interval(0, 1));
                    #endregion

                    #region 求XY方向缩放比的下限
                    #region 0
                    Line lineX0 = new Line(basePoint, ptList0[(basePointIndex+1)%ptList0.Count]);
                    Line lineY0 = new Line(basePoint, ptList0[((basePointIndex - 1) + ptList0.Count) % ptList0.Count]);
                    double lengthX0 = lineX0.Length;
                    double lengthY0 = lineY0.Length;

                    double scaleFactor_X0_Min;
                    double scaleFactor_Y0_Min;
                    if (lengthX0 > w)
                    {
                        scaleFactor_X0_Min = w / lengthX0;
                    }
                    else
                    {
                        scaleFactor_X0_Min = 1;
                    }
                    if (lengthY0 > w)
                    {
                        scaleFactor_Y0_Min = w / lengthY0;
                    }
                    else
                    {
                        scaleFactor_Y0_Min = 1;
                    }

                    double scaleFactor_X0;
                    double scaleFactor_Y0;
                    if (lengthX0 > lForHighRise)
                    {
                        scaleFactor_X0 = lForHighRise / lengthX0;
                    }
                    else
                    {
                        scaleFactor_X0 = 1;
                    }
                    if (lengthY0 > lForHighRise)
                    {
                        scaleFactor_Y0 = lForHighRise / lengthY0;
                    }
                    else
                    {
                        scaleFactor_Y0 = 1;
                    }
                    #endregion
                    #endregion

                    #region XY方向缩放比与下限比较
                    #region 0
                    if (scaleFactor_X0 < scaleFactor_X0_Min)
                    {
                        scaleFactor_X0 = scaleFactor_X0_Min;
                    }
                    if (scaleFactor_Y0 < scaleFactor_Y0_Min)
                    {
                        scaleFactor_Y0 = scaleFactor_Y0_Min;
                    }
                    #endregion
                    #endregion

                    #region 取子面
                    #region 0
                    UVInterval uvinterval0;
                    if (basePointIndex == 0)
                    {
                        uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                    }
                    else if (basePointIndex == 1)
                    {
                        uvinterval0 = new UVInterval(new Interval(0, scaleFactor_X0), new Interval((1 - scaleFactor_Y0), 1));
                    }
                    else if (basePointIndex == 2)
                    {
                        uvinterval0 = new UVInterval(new Interval(0, scaleFactor_Y0), new Interval(0, scaleFactor_X0));
                    }
                    else
                    {
                        uvinterval0 = new UVInterval(new Interval((1 - scaleFactor_X0), 1), new Interval(0, scaleFactor_Y0));
                    }
                    Surface trimedSurf0 = UtilityFunctions.IsoTrim(surf0, uvinterval0);
                    Curve[] crvs0 = trimedSurf0.ToBrep().DuplicateEdgeCurves();
                    Curve crv0 = Curve.JoinCurves(crvs0)[0];
                    #endregion
                    #endregion

                    this.crvForHighRise = crv0.DuplicateCurve();
                }
                else
                {
                    if (turningPointIndex == 0 || turningPointIndex == 1)
                    {
                        // 上 ptList1

                        // 找到可以放置高层的位置
                        //List<Point3d> pointsForCanPlaceHighrise_Two = new List<Point3d>();
                        //List<Point3d> pointsForCanPlaceHighrise_One = new List<Point3d>();

                        List<int> indexForCanPlaceHighrise_Two = new List<int>();
                        List<int> indexForCanPlaceHighrise_One = new List<int>();

                        #region west
                        bool west_0 = false;
                        bool west_1 = false;

                        if (this.West1_Interval != null)
                        {
                            for (int i = 0; i < this.West1_Interval.Count; i++)
                            {
                                if (this.West1_Interval[i].IncludesParameter(0))
                                {
                                    west_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.West1_Interval.Count; i++)
                            {
                                if (this.West1_Interval[i].IncludesParameter(1))
                                {
                                    west_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region Mid
                        bool mid_0 = false;
                        bool mid_1 = false;

                        if (this.Middle_Interval != null)
                        {
                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(0))
                                {
                                    mid_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(1))
                                {
                                    mid_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region east
                        bool east_0 = false;
                        bool east_1 = false;

                        if (this.East1_Interval != null)
                        {
                            for (int i = 0; i < this.East1_Interval.Count; i++)
                            {
                                if (this.East1_Interval[i].IncludesParameter(0))
                                {
                                    east_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.East1_Interval.Count; i++)
                            {
                                if (this.East1_Interval[i].IncludesParameter(1))
                                {
                                    east_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region north
                        bool north_0 = false;
                        bool north_1 = false;

                        if (this.North_Interval != null)
                        {
                            for (int i = 0; i < this.North_Interval.Count; i++)
                            {
                                if (this.North_Interval[i].IncludesParameter(0))
                                {
                                    north_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.North_Interval.Count; i++)
                            {
                                if (this.North_Interval[i].IncludesParameter(1))
                                {
                                    north_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region 西南角点
                        if (west_0 && mid_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(0);
                        }
                        else if (west_0 || mid_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(0);
                        }
                        #endregion

                        #region 东南角点
                        if (mid_1 && east_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(1);
                        }
                        else if (mid_1 || east_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(1);
                        }
                        #endregion

                        #region 东北角点
                        if (east_1 && north_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(2);
                        }
                        else if (east_1 || north_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(2);
                        }
                        #endregion

                        #region 西北角点
                        if (north_0 && west_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(3);
                        }
                        else if (north_0 || west_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(3);
                        }
                        #endregion

                        #region 随机取点
                        if (indexForCanPlaceHighrise_Two.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_Two.Count);
                            basePoint = ptList1[indexForCanPlaceHighrise_Two[randomIndex]];
                        }
                        else if (indexForCanPlaceHighrise_One.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_One.Count);
                            basePoint = ptList1[indexForCanPlaceHighrise_One[randomIndex]];
                        }
                        else
                        {
                            // 应该不存在此种情况
                            int randomIndex = m_random.Next(4);
                            basePoint = ptList1[randomIndex];
                        }
                        #endregion
                    }
                    else
                    {
                        // turningPointIndex == 2 || turningPointIndex == 3
                        // 下 ptList1

                        // 找到可以放置高层的位置
                        //List<Point3d> pointsForCanPlaceHighrise_Two = new List<Point3d>();
                        //List<Point3d> pointsForCanPlaceHighrise_One = new List<Point3d>();

                        List<int> indexForCanPlaceHighrise_Two = new List<int>();
                        List<int> indexForCanPlaceHighrise_One = new List<int>();

                        #region west
                        bool west_0 = false;
                        bool west_1 = false;

                        if (this.West_Interval != null)
                        {
                            for (int i = 0; i < this.West_Interval.Count; i++)
                            {
                                if (this.West_Interval[i].IncludesParameter(0))
                                {
                                    west_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.West_Interval.Count; i++)
                            {
                                if (this.West_Interval[i].IncludesParameter(1))
                                {
                                    west_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region south
                        bool south_0 = false;
                        bool south_1 = false;

                        if (this.South_Interval != null)
                        {
                            for (int i = 0; i < this.South_Interval.Count; i++)
                            {
                                if (this.South_Interval[i].IncludesParameter(0))
                                {
                                    south_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.South_Interval.Count; i++)
                            {
                                if (this.South_Interval[i].IncludesParameter(1))
                                {
                                    south_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region east
                        bool east_0 = false;
                        bool east_1 = false;

                        if (this.East_Interval != null)
                        {
                            for (int i = 0; i < this.East_Interval.Count; i++)
                            {
                                if (this.East_Interval[i].IncludesParameter(0))
                                {
                                    east_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.East_Interval.Count; i++)
                            {
                                if (this.East_Interval[i].IncludesParameter(1))
                                {
                                    east_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region Mid
                        bool mid_0 = false;
                        bool mid_1 = false;

                        if (this.Middle_Interval != null)
                        {
                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(0))
                                {
                                    mid_0 = true;
                                    break;
                                }
                            }

                            for (int i = 0; i < this.Middle_Interval.Count; i++)
                            {
                                if (this.Middle_Interval[i].IncludesParameter(1))
                                {
                                    mid_1 = true;
                                    break;
                                }
                            }
                        }
                        #endregion

                        #region 西南角点
                        if (west_0 && south_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(0);
                        }
                        else if (west_0 || south_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.SouthBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(0);
                        }
                        #endregion

                        #region 东南角点
                        if (south_1 && east_0)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_Two.Add(1);
                        }
                        else if (south_1 || east_0)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.EastBaseLine.From);
                            indexForCanPlaceHighrise_One.Add(1);
                        }
                        #endregion

                        #region 东北角点
                        if (east_1 && mid_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(2);
                        }
                        else if (east_1 || mid_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.MiddleBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(2);
                        }
                        #endregion

                        #region 西北角点
                        if (mid_0 && west_1)
                        {
                            //pointsForCanPlaceHighrise_Two.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_Two.Add(3);
                        }
                        else if (mid_0 || west_1)
                        {
                            //pointsForCanPlaceHighrise_One.Add(this.WestBaseLine.To);
                            indexForCanPlaceHighrise_One.Add(3);
                        }
                        #endregion

                        #region 随机取点
                        if (indexForCanPlaceHighrise_Two.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_Two.Count);
                            basePoint = ptList1[indexForCanPlaceHighrise_Two[randomIndex]];
                        }
                        else if (indexForCanPlaceHighrise_Two.Count != 0)
                        {
                            int randomIndex = m_random.Next(indexForCanPlaceHighrise_Two.Count);
                            basePoint = ptList1[indexForCanPlaceHighrise_Two[randomIndex]];
                        }
                        else
                        {
                            // 应该不存在此种情况
                            int randomIndex = m_random.Next(4);
                            basePoint = ptList1[randomIndex];
                        }
                        #endregion
                    }

                    //if (randomCode1 == 0)
                    //{
                    //    basePoint = ptList1[0];
                    //}
                    //else if (randomCode1 == 1)
                    //{
                    //    basePoint = ptList1[1];
                    //}
                    //else if (randomCode1 == 2)
                    //{
                    //    basePoint = ptList1[2];
                    //}
                    //else
                    //{
                    //    basePoint = ptList1[3];
                    //}

                    ptList1.RemoveAt(ptList1.Count - 1);
                    int basePointIndex = ptList1.IndexOf(basePoint);

                    #region 构造第1个Brep，并且得到surf
                    Curve brep1_railCrv0;
                    Curve brep1_railCrv1;
                    List<Curve> brep1_shapes = new List<Curve>();

                    brep1_railCrv0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[0], ptList1[3] }, 1);
                    brep1_railCrv1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[1], ptList1[2] }, 1);
                    Curve brep1_shape0 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[0], ptList1[1] }, 1);
                    Curve brep1_shape1 = Curve.CreateInterpolatedCurve(new Point3d[] { ptList1[3], ptList1[2] }, 1);
                    brep1_shapes.Add(brep1_shape0);
                    brep1_shapes.Add(brep1_shape1);

                    Brep brp1 = Brep.CreateFromSweep(brep1_railCrv0, brep1_railCrv1, brep1_shapes, false, Tolerance)[0];
                    Surface surf1 = brp1.Surfaces[0];
                    surf1.SetDomain(0, new Interval(0, 1));
                    surf1.SetDomain(1, new Interval(0, 1));
                    #endregion

                    #region 求XY方向缩放比的下限
                    #region 1
                    Line lineX1 = new Line(basePoint, ptList1[(basePointIndex + 1) % ptList1.Count]);
                    Line lineY1 = new Line(basePoint, ptList1[((basePointIndex - 1) + ptList1.Count) % ptList1.Count]);
                    double lengthX1 = lineX1.Length;
                    double lengthY1 = lineY1.Length;

                    double scaleFactor_X1_Min;
                    double scaleFactor_Y1_Min;
                    if (lengthX1 > w)
                    {
                        scaleFactor_X1_Min = w / lengthX1;
                    }
                    else
                    {
                        scaleFactor_X1_Min = 1;
                    }
                    if (lengthY1 > w)
                    {
                        scaleFactor_Y1_Min = w / lengthY1;
                    }
                    else
                    {
                        scaleFactor_Y1_Min = 1;
                    }

                    double scaleFactor_X1;
                    double scaleFactor_Y1;
                    if (lengthX1 > lForHighRise)
                    {
                        scaleFactor_X1 = lForHighRise / lengthX1;
                    }
                    else
                    {
                        scaleFactor_X1 = 1;
                    }
                    if (lengthY1 > lForHighRise)
                    {
                        scaleFactor_Y1 = lForHighRise / lengthY1;
                    }
                    else
                    {
                        scaleFactor_Y1 = 1;
                    }
                    #endregion
                    #endregion

                    #region XY方向缩放比与下限比较
                    #region 1
                    if (scaleFactor_X1 < scaleFactor_X1_Min)
                    {
                        scaleFactor_X1 = scaleFactor_X1_Min;
                    }
                    if (scaleFactor_Y1 < scaleFactor_Y1_Min)
                    {
                        scaleFactor_Y1 = scaleFactor_Y1_Min;
                    }
                    #endregion
                    #endregion

                    #region 取子面
                    #region 1
                    UVInterval uvinterval1;
                    if (basePointIndex == 0)
                    {
                        uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                    }
                    else if (basePointIndex == 1)
                    {
                        uvinterval1 = new UVInterval(new Interval(0, scaleFactor_X1), new Interval((1 - scaleFactor_Y1), 1));
                    }
                    else if (basePointIndex == 2)
                    {
                        uvinterval1 = new UVInterval(new Interval(0, scaleFactor_Y1), new Interval(0, scaleFactor_X1));
                    }
                    else
                    {
                        uvinterval1 = new UVInterval(new Interval((1 - scaleFactor_X1), 1), new Interval(0, scaleFactor_Y1));
                    }
                    Surface trimedSurf1 = UtilityFunctions.IsoTrim(surf1, uvinterval1);
                    Curve[] crvs1 = trimedSurf1.ToBrep().DuplicateEdgeCurves();
                    Curve crv1 = Curve.JoinCurves(crvs1)[0];
                    #endregion
                    #endregion

                    this.crvForHighRise = crv1.DuplicateCurve();
                }
            }

            return this;
        }

        public override List<Curve> GetAllHorizontal()
        {
            List<Curve> horizontalLines = new List<Curve>();
            if (this.SouthBaseLine != Line.Unset && this.South_Interval != null)
            {
                for (int i = 0; i < this.South_Interval.Count; i++)
                {
                    Point3d p0 = this.SouthBaseLine.PointAt(this.South_Interval[i].T0);
                    Point3d p1 = this.SouthBaseLine.PointAt(this.South_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.MiddleBaseLine != Line.Unset && this.Middle_Interval != null)
            {
                for (int i = 0; i < this.Middle_Interval.Count; i++)
                {
                    Point3d p0 = this.MiddleBaseLine.PointAt(this.Middle_Interval[i].T0);
                    Point3d p1 = this.MiddleBaseLine.PointAt(this.Middle_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.NorthBaseLine != Line.Unset && this.North_Interval != null)
            {
                for (int i = 0; i < this.North_Interval.Count; i++)
                {
                    Point3d p0 = this.NorthBaseLine.PointAt(this.North_Interval[i].T0);
                    Point3d p1 = this.NorthBaseLine.PointAt(this.North_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.HSouthBaseLine != Line.Unset && this.HSouth_Interval != null)
            {
                for (int i = 0; i < this.HSouth_Interval.Count; i++)
                {
                    Point3d p0 = this.HSouthBaseLine.PointAt(this.HSouth_Interval[i].T0);
                    Point3d p1 = this.HSouthBaseLine.PointAt(this.HSouth_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.HNorthBaseLine != Line.Unset && this.HNouth_Interval != null)
            {
                for (int i = 0; i < this.HNouth_Interval.Count; i++)
                {
                    Point3d p0 = this.HNorthBaseLine.PointAt(this.HNouth_Interval[i].T0);
                    Point3d p1 = this.HNorthBaseLine.PointAt(this.HNouth_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            return horizontalLines;
        }

        public override List<Curve> GetAllVertical()
        {
            List<Curve> verticalLines = new List<Curve>();
            if (this.WestBaseLine != Line.Unset && this.West_Interval != null)
            {
                for (int i = 0; i < this.West_Interval.Count; i++)
                {
                    Point3d p0 = this.WestBaseLine.PointAt(this.West_Interval[i].T0);
                    Point3d p1 = this.WestBaseLine.PointAt(this.West_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.WestBaseLine1 != Line.Unset && this.West1_Interval != null)
            {
                for (int i = 0; i < this.West1_Interval.Count; i++)
                {
                    Point3d p0 = this.WestBaseLine1.PointAt(this.West1_Interval[i].T0);
                    Point3d p1 = this.WestBaseLine1.PointAt(this.West1_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.EastBaseLine != Line.Unset && this.East_Interval != null)
            {
                for (int i = 0; i < this.East_Interval.Count; i++)
                {
                    Point3d p0 = this.EastBaseLine.PointAt(this.East_Interval[i].T0);
                    Point3d p1 = this.EastBaseLine.PointAt(this.East_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.EastBaseLine1 != Line.Unset && this.East1_Interval != null)
            {
                for (int i = 0; i < this.East1_Interval.Count; i++)
                {
                    Point3d p0 = this.EastBaseLine1.PointAt(this.East1_Interval[i].T0);
                    Point3d p1 = this.EastBaseLine1.PointAt(this.East1_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.VSouthBaseLine != Line.Unset && this.VSouth_Interval != null)
            {
                for (int i = 0; i < this.VSouth_Interval.Count; i++)
                {
                    Point3d p0 = this.VSouthBaseLine.PointAt(this.VSouth_Interval[i].T0);
                    Point3d p1 = this.VSouthBaseLine.PointAt(this.VSouth_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.VNorthBaseLine != Line.Unset && this.VNorth_Interval != null)
            {
                for (int i = 0; i < this.VNorth_Interval.Count; i++)
                {
                    Point3d p0 = this.VNorthBaseLine.PointAt(this.VNorth_Interval[i].T0);
                    Point3d p1 = this.VNorthBaseLine.PointAt(this.VNorth_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            return verticalLines;
        }

        public void LengthConstraint(double lMin, double lMax, double d, double w, Random random)
        {
            if (this.SouthBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.SouthBaseLine.PointAt(0), this.SouthBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd;
                        if (length - lMax > dw)
                        {
                            newEnd = line.PointAtLength(lMax);
                        }
                        else
                        {
                            newEnd = line.PointAtLength(length - dw);
                        }
                        
                        double tEnd = line.ClosestParameter(newEnd);
                        this.South_Interval = new List<Interval>();
                        if (random.Next(0, 2) == 0)
                        {
                            this.South_Interval.Add(new Interval(0, tEnd));
                        }
                        else
                        {
                            this.South_Interval.Add(new Interval(1 - tEnd, 1));
                        }
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.South_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.South_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                    else
                    {
                        // 按照 lMaxCount 和 lMinCount 中的小的数目走
                        int count = lMinCount > lMaxCount ? lMaxCount : lMinCount;
                        bool flag = lMinCount > lMaxCount ? true : false;
                        if (flag && lMaxCount == 1)
                        {
                            count = lMinCount;
                            flag = false;
                        }
                        List<double> factors = GetFactors(flag, lMin, lMax, dw, count);

                        this.South_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.South_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                }
            }
            if (this.NorthBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.NorthBaseLine.PointAt(0), this.NorthBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd;
                        if (length - lMax > dw)
                        {
                            newEnd = line.PointAtLength(lMax);
                        }
                        else
                        {
                            newEnd = line.PointAtLength(length - dw);
                        }

                        double tEnd = line.ClosestParameter(newEnd);
                        this.North_Interval = new List<Interval>();
                        if (random.Next(0, 2) == 0)
                        {
                            this.North_Interval.Add(new Interval(0, tEnd));
                        }
                        else
                        {
                            this.North_Interval.Add(new Interval(1 - tEnd, 1));
                        }
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.North_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.North_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                    else
                    {
                        // 按照 lMaxCount 和 lMinCount 中的小的数目走
                        int count = lMinCount > lMaxCount ? lMaxCount : lMinCount;
                        bool flag = lMinCount > lMaxCount ? true : false;
                        if (flag && lMaxCount == 1)
                        {
                            count = lMinCount;
                            flag = false;
                        }
                        List<double> factors = GetFactors(flag, lMin, lMax, dw, count);

                        this.North_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.North_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                }
            }
            if (this.MiddleBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.MiddleBaseLine.PointAt(0), this.MiddleBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd;
                        if (length - lMax > dw)
                        {
                            newEnd = line.PointAtLength(lMax);
                        }
                        else
                        {
                            newEnd = line.PointAtLength(length - dw);
                        }

                        double tEnd = line.ClosestParameter(newEnd);
                        this.Middle_Interval = new List<Interval>();
                        if (random.Next(0, 2) == 0)
                        {
                            this.Middle_Interval.Add(new Interval(0, tEnd));
                        }
                        else
                        {
                            this.Middle_Interval.Add(new Interval(1 - tEnd, 1));
                        }
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.Middle_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.Middle_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                    else
                    {
                        // 按照 lMaxCount 和 lMinCount 中的小的数目走
                        int count = lMinCount > lMaxCount ? lMaxCount : lMinCount;
                        bool flag = lMinCount > lMaxCount ? true : false;
                        if (flag && lMaxCount == 1)
                        {
                            count = lMinCount;
                            flag = false;
                        }
                        List<double> factors = GetFactors(flag, lMin, lMax, dw, count);

                        this.Middle_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.Middle_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                }
            }
            if (this.HSouthBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.HSouthBaseLine.PointAt(0), this.HSouthBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd;
                        if (length - lMax > dw)
                        {
                            newEnd = line.PointAtLength(lMax);
                        }
                        else
                        {
                            newEnd = line.PointAtLength(length - dw);
                        }

                        double tEnd = line.ClosestParameter(newEnd);
                        this.HSouth_Interval = new List<Interval>();
                        if (random.Next(0, 2) == 0)
                        {
                            this.HSouth_Interval.Add(new Interval(0, tEnd));
                        }
                        else
                        {
                            this.HSouth_Interval.Add(new Interval(1 - tEnd, 1));
                        }
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.HSouth_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.HSouth_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                    else
                    {
                        // 按照 lMaxCount 和 lMinCount 中的小的数目走
                        int count = lMinCount > lMaxCount ? lMaxCount : lMinCount;
                        bool flag = lMinCount > lMaxCount ? true : false;
                        if (flag && lMaxCount == 1)
                        {
                            count = lMinCount;
                            flag = false;
                        }
                        List<double> factors = GetFactors(flag, lMin, lMax, dw, count);

                        this.HSouth_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.HSouth_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                }
            }
            if (this.HNorthBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.HNorthBaseLine.PointAt(0), this.HNorthBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd;
                        if (length - lMax > dw)
                        {
                            newEnd = line.PointAtLength(lMax);
                        }
                        else
                        {
                            newEnd = line.PointAtLength(length - dw);
                        }

                        double tEnd = line.ClosestParameter(newEnd);
                        this.HNouth_Interval = new List<Interval>();
                        this.HNouth_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.HNouth_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.HNouth_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                    else
                    {
                        // 按照 lMaxCount 和 lMinCount 中的小的数目走
                        int count = lMinCount > lMaxCount ? lMaxCount : lMinCount;
                        bool flag = lMinCount > lMaxCount ? true : false;
                        if (flag && lMaxCount == 1)
                        {
                            count = lMinCount;
                            flag = false;
                        }
                        List<double> factors = GetFactors(flag, lMin, lMax, dw, count);

                        this.HNouth_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.HNouth_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
                            iter += 2;
                        }
                    }
                }
            }
        }

        private static List<double> GetFactors(bool flag, double lMin, double lMax, double dw, int count)
        {
            double sum = 0;
            if (flag)
            {
                sum = (lMax + dw) * count - dw;
            }
            else
            {
                sum = (lMin + dw) * count - dw;
            }
            List<double> eachLength = new List<double>();
            //eachLength.Add(0);
            for (int i = 0; i < 2 * count - 1; i++)
            {
                if (i % 2 == 0)
                {
                    if (flag)
                    {
                        eachLength.Add(lMax);
                    }
                    else
                    {
                        eachLength.Add(lMin);
                    }
                }
                else
                {
                    eachLength.Add(dw);
                }
            }
            List<double> cumulatives = new List<double>();
            cumulatives.Add(eachLength[0]);
            for (int i = 1; i < eachLength.Count; i++)
            {
                cumulatives.Add(cumulatives[i - 1] + eachLength[i]);
            }
            //cumulatives.Add(0);
            //cumulatives.Add(lMin);
            //cumulatives.Add(lMin + dw);
            //cumulatives.Add((lMin + dw) * 2 - dw);
            List<double> factors = new List<double>();
            factors.Add(0);
            for (int i = 0; i < 2 * count - 1; i++)
            {
                double factor = cumulatives[i] / sum;
                factors.Add(factor);
            }
            return factors;
        }

        public override double GetAllHCurveLength(out int hCount)
        {
            double hLength = 0;
            List<Curve> hCurves = this.GetAllHorizontal();
            hCount = hCurves.Count;
            for (int i = 0; i < hCurves.Count; i++)
            {
                hLength += hCurves[i].GetLength();
            }
            return hLength;
        }

        public override double GetAllVCurveLength(out int vCount)
        {
            double vLength = 0;
            List<Curve> vCurves = this.GetAllVertical();
            vCount = vCurves.Count;
            for (int i = 0; i < vCurves.Count; i++)
            {
                vLength += vCurves[i].GetLength();
            }
            return vLength;
        }
    }
}
