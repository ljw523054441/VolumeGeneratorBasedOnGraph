using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class BilinearCell : Cell
    {
        public Line HCenterBaseLine { get; set; } = Line.Unset;
        public Line VCenterBaseLine { get; set; } = Line.Unset;

        //public List<Line> HLines { get; set; }
        //public List<Line> VLines { get; set; }

        public List<Interval> HCenter_Interval { get; set; }
        public List<Interval> VCenter_Interval { get; set; }

        public double CenterTValue { get; set; }

        public enum MainDirection
        {
            Unset = 0,
            SouthNorth = 1,
            WestEast = 2
        }

        public MainDirection MainVolumeDirection { get; set; }

        public int BoundaryBaseLineCount { get; set; }
        public int CenterBaseLineCount { get; set; }

        List<List<double>> KeyPointTLoL { get; set; }

        List<List<bool>> KeyPointTypeLoL { get; set; }

        public enum BilinearShapeType
        {
            Unset = 0,
            CShape = 1,
            IShape = 2,
            RecShape = 3,
            RecVariantShape = 4,
            SingleRegion = 5,
            SinglePolyline = 6
        }

        public BilinearShapeType ShapeType { get; set; }

        public BilinearCell(Curve southLine, 
                            Curve northLine, 
                            Curve westLine, 
                            Curve eastLine, 

                            Curve southBoundary,
                            Curve northBoundary,
                            Curve westBoundary,
                            Curve eastBoundary,

                            //double w, 
                            int i, 
                            int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region baseLine
            int count = 0;
            if (southLine != null)
            {
                this.SouthBaseLine = new Line(southLine.PointAtStart, southLine.PointAtEnd);
                count++;
            }
            if (northLine != null)
            {
                this.NorthBaseLine = new Line(northLine.PointAtStart, northLine.PointAtEnd);
                count++;
            }
            if (westLine != null)
            {
                this.WestBaseLine = new Line(westLine.PointAtStart, westLine.PointAtEnd);
                count++;
            }
            if (eastLine != null)
            {
                this.EastBaseLine = new Line(eastLine.PointAtStart, eastLine.PointAtEnd);
                count++;
            }
            #endregion

            if (count == 4)
            {
                this.ShapeType = BilinearShapeType.RecShape;
            }
            else if (count == 3)
            {
                this.ShapeType = BilinearShapeType.CShape;
            }

            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            List<Curve> extendedBoundaryList = new List<Curve>();
            extendedBoundaryList.Add(southBoundary);
            extendedBoundaryList.Add(eastBoundary);
            extendedBoundaryList.Add(northBoundary);
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

            Curve[] joinedCrv = Curve.JoinCurves(boundaryList);
            this.CellBoundary = joinedCrv[0];

            #endregion
        }

        /// <summary>
        /// SingleRegion 以及 SinglePolyline
        /// </summary>
        /// <param name="southLine"></param>
        /// <param name="northLine"></param>
        /// <param name="westLine"></param>
        /// <param name="eastLine"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public BilinearCell(Curve boundary, Curve southLine,Curve northLine, Curve westLine, Curve eastLine,int i,int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region baseLine 和 interval
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
            #endregion

            if (count != 4)
            {
                this.ShapeType = BilinearShapeType.SinglePolyline;

                this.BoundaryBaseLineCount = 2;
                this.CenterBaseLineCount = 0;

                #region CellBoundary
                this.CellBoundary = boundary.DuplicateCurve();
                #endregion
            }
            else
            {
                this.ShapeType = BilinearShapeType.SingleRegion;

                this.BoundaryBaseLineCount = 4;
                this.CenterBaseLineCount = 0;

                #region CellBoundary
                this.CellBoundary = boundary.DuplicateCurve();
                #endregion
            }
        }


        public BilinearCell(Curve centerH, 
                            Curve westLine, 
                            Curve eastLine, 

                            Curve southBoundary,
                            Curve northBoundary,
                            Curve westBoundary,
                            Curve eastBoundary,

                            int i, 
                            int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region baseLine
            int centerCount = 0;
            int count = 0;
            if (centerH != null)
            {
                this.HCenterBaseLine = new Line(centerH.PointAtStart, centerH.PointAtEnd);
                centerCount++;
            }
            if (westLine != null)
            {
                this.WestBaseLine = new Line(westLine.PointAtStart, westLine.PointAtEnd);
                count++;
            }
            if (eastLine != null)
            {
                this.EastBaseLine = new Line(eastLine.PointAtStart, eastLine.PointAtEnd);
                count++;
            }
            #endregion

            this.ShapeType = BilinearShapeType.IShape;
            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = count;

            #region CellBoundary
            List<Curve> extendedBoundaryList = new List<Curve>();
            extendedBoundaryList.Add(southBoundary);
            extendedBoundaryList.Add(eastBoundary);
            extendedBoundaryList.Add(northBoundary);
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

            Curve[] joinedCrv = Curve.JoinCurves(boundaryList);
            this.CellBoundary = joinedCrv[0];

            #endregion
        }

        /// <summary>
        /// 生成IShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="isSouthNorth"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, bool isSouthNorth, double t, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region Interval
            if (isSouthNorth)
            {
                this.South_Interval = new List<Interval>();
                if (initialCell.South_Interval == null)
                {
                    this.South_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    this.South_Interval.AddRange(initialCell.South_Interval);
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

                this.VCenter_Interval = new List<Interval>();
                this.VCenter_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));

                this.HCenter_Interval = new List<Interval>();
                if (initialCell.HCenter_Interval == null)
                {
                    this.HCenter_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    this.HCenter_Interval.AddRange(initialCell.HCenter_Interval);
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
            if (isSouthNorth)
            {
                // 南北体量为主
                Line south = this.SouthBaseLine;
                Line north = this.NorthBaseLine;
                Point3d pointOnSouth = south.PointAt(t);
                Point3d pointOnNorth = north.PointAt(t);
                this.VCenterBaseLine = new Line(pointOnSouth, pointOnNorth);

                this.HCenterBaseLine = initialCell.HCenterBaseLine;
            }
            else
            {
                // 东西体量为主
                Line west = this.WestBaseLine;
                Line east = this.EastBaseLine;
                Point3d pointOnWest = west.PointAt(t);
                Point3d pointOnEast = east.PointAt(t);
                this.HCenterBaseLine = new Line(pointOnWest, pointOnEast);

                this.VCenterBaseLine = initialCell.VCenterBaseLine;
            }
            #endregion

            if (isSouthNorth)
            {
                this.MainVolumeDirection = MainDirection.SouthNorth;
            }
            else
            {
                this.MainVolumeDirection = MainDirection.WestEast;
            }

            this.CenterTValue = t;
            this.ShapeType = BilinearShapeType.IShape;

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();

            #region CellBoundary
            if (initialCell.CellBoundary != null)
            {
                this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            }
            
            #endregion
        }

        /// <summary>
        /// 生成CShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, int directionCode, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region Interval
            switch (directionCode)
            {
                case 0:
                    this.South_Interval = new List<Interval>();
                    if (initialCell.South_Interval == null)
                    {
                        this.South_Interval.Add(new Interval(0, 1));
                    }
                    else
                    {
                        this.South_Interval.AddRange(initialCell.South_Interval);
                    }

                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 0;
                    break;
                case 1:
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

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 1;
                    break;
                case 2:
                    this.South_Interval = new List<Interval>();
                    if (initialCell.South_Interval == null)
                    {
                        this.South_Interval.Add(new Interval(0, 1));
                    }
                    else
                    {
                        this.South_Interval.AddRange(initialCell.South_Interval);
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

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 0;
                    break;
                case 3:
                    this.South_Interval = new List<Interval>();
                    if (initialCell.South_Interval == null)
                    {
                        this.South_Interval.Add(new Interval(0, 1));
                    }
                    else
                    {
                        this.South_Interval.AddRange(initialCell.South_Interval);
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
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 1;
                    break;
            }
            #endregion

            #region BaseLine
            // 父类属性
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            // 子类属性
            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;
            #endregion

            this.ShapeType = BilinearShapeType.CShape;

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();

            #region CellBoundary
            if (initialCell.CellBoundary != null)
            {
                this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            }
            
            #endregion
        }

        /// <summary>
        /// 生成RecShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

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
            #endregion

            #region baseLine
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;

            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;
            #endregion

            this.MainVolumeDirection = MainDirection.Unset;
            this.CenterTValue = 0;

            this.ShapeType = BilinearShapeType.RecShape;

            this.BoundaryBaseLineCount = 4;
            this.CenterBaseLineCount = 0;

            #region CellBoundary
            if (initialCell.CellBoundary != null)
            {
                this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            }
            
            #endregion
        }

        /// <summary>
        /// 生成RecVariantShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell ,double t, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region interval
            this.South_Interval = new List<Interval>();
            if (initialCell.South_Interval == null)
            {
                this.South_Interval.Add(new Interval(0, 1));
            }
            else
            {
                this.South_Interval.AddRange(initialCell.South_Interval);
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
            #endregion

            #region baseLine
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;

            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;

            Point3d pointOnSouth0 = this.SouthBaseLine.PointAt(t);
            Point3d pointOnSouth1 = this.SouthBaseLine.PointAt(1 - t);
            Point3d pointOnNorth0 = this.NorthBaseLine.PointAt(t);
            Point3d pointOnNorth1 = this.NorthBaseLine.PointAt(1 - t);
            this.WestBaseLine = new Line(pointOnSouth0, pointOnNorth0);
            this.EastBaseLine = new Line(pointOnSouth1, pointOnNorth1);
            #endregion

            this.MainVolumeDirection = MainDirection.SouthNorth;
            this.CenterTValue = t;
            this.ShapeType = BilinearShapeType.RecVariantShape;

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();

            #region CellBoundary
            if (initialCell.CellBoundary != null)
            {
                this.CellBoundary = initialCell.CellBoundary.DuplicateCurve();
            }

            #endregion
        }

        public BilinearCell GenerateIAndIVariantShape(int randomT,int directionCode, double w)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                double length0 = this.WestBaseLine.Length;
                double length1 = this.EastBaseLine.Length;
                double length = length0 < length1 ? length0 : length1;
                double tMin = w / length;
                double tMax = (length - w) / length;
                double t = Math.Abs(tMin + (tMax - tMin) * randomT);

                newGeneratedCell = new BilinearCell(this, false, t, this.CurrentIndex[0], this.CurrentIndex[1]);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCAndCVariantShape(randomT, directionCode, w);
            }
            else
            {
                if (directionCode >= 2)
                {
                    // 南北体量为主
                    double length0 = this.SouthBaseLine.Length;
                    double length1 = this.NorthBaseLine.Length;
                    double length = length0 < length1 ? length0 : length1;
                    double tMin = w / length;
                    double tMax = (length - w) / length;
                    double t = Math.Abs(tMin + (tMax - tMin) * randomT);

                    newGeneratedCell = new BilinearCell(this, true, t, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    // 东西体量为主
                    double length0 = this.WestBaseLine.Length;
                    double length1 = this.EastBaseLine.Length;
                    double length = length0 < length1 ? length0 : length1;
                    double tMin = w / length;
                    double tMax = (length - w) / length;
                    double t = Math.Abs(tMin + (tMax - tMin) * randomT);

                    newGeneratedCell = new BilinearCell(this, false, t, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateCAndCVariantShape(int randomT,int directionCode, double w)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIAndIVariantShape(randomT, directionCode, w);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                if (this.SouthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 0, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.NorthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 1, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.WestBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 2, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    newGeneratedCell = new BilinearCell(this, 3, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            else
            {
                switch (directionCode)
                {
                    case 0:
                        newGeneratedCell = new BilinearCell(this, 0, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 1:
                        newGeneratedCell = new BilinearCell(this, 1, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 2:
                        newGeneratedCell = new BilinearCell(this, 2, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 3:
                        newGeneratedCell = new BilinearCell(this, 3, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    default:
                        newGeneratedCell = new BilinearCell(this, 0, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateRecAndRecVariantShape(int randomT, int directionCode, double lMin, double w, bool flag)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIAndIVariantShape(randomT, directionCode, w);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCAndCVariantShape(randomT, directionCode, w);
            }
            else
            {
                if (directionCode < 2)
                {
                    // 执行GenerateRecShape函数
                    newGeneratedCell = new BilinearCell(this, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    if (flag)
                    {
                        // 执行GenerateRecShape函数
                        newGeneratedCell = new BilinearCell(this, this.CurrentIndex[0], this.CurrentIndex[1]);
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

                            newGeneratedCell = new BilinearCell(this, t, this.CurrentIndex[0], this.CurrentIndex[1]);
                        }
                        else
                        {
                            // 执行GenerateRecShape函数
                            newGeneratedCell = new BilinearCell(this, this.CurrentIndex[0], this.CurrentIndex[1]);
                        }
                    }
                }
            }
            return newGeneratedCell;
        }

        private int GetBoundaryBaseLinesCount()
        {
            int count = 0;
            if (this.SouthBaseLine != Line.Unset)
            {
                count++;
            }
            if (this.NorthBaseLine != Line.Unset)
            {
                count++;
            }
            if (this.WestBaseLine != Line.Unset)
            {
                count++;
            }
            if (this.EastBaseLine != Line.Unset)
            {
                count++;
            }
            return count;
        }

        private int GetCenterBaseLineCount()
        {
            int count = 0;
            if (this.HCenterBaseLine != Line.Unset)
            {
                count++;
            }
            if (this.VCenterBaseLine != Line.Unset)
            {
                count++;
            }
            return count;
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
            if (this.NorthBaseLine != Line.Unset && this.North_Interval != null)
            {
                for (int i = 0; i < this.North_Interval.Count; i++)
                {
                    Point3d p0 = this.NorthBaseLine.PointAt(this.North_Interval[i].T0);
                    Point3d p1 = this.NorthBaseLine.PointAt(this.North_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.HCenterBaseLine != Line.Unset && this.HCenter_Interval != null)
            {
                for (int i = 0; i < this.HCenter_Interval.Count; i++)
                {
                    Point3d p0 = this.HCenterBaseLine.PointAt(this.HCenter_Interval[i].T0);
                    Point3d p1 = this.HCenterBaseLine.PointAt(this.HCenter_Interval[i].T1);
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
            if (this.EastBaseLine != Line.Unset && this.East_Interval != null)
            {
                for (int i = 0; i < this.East_Interval.Count; i++)
                {
                    Point3d p0 = this.EastBaseLine.PointAt(this.East_Interval[i].T0);
                    Point3d p1 = this.EastBaseLine.PointAt(this.East_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.VCenterBaseLine != Line.Unset && this.VCenter_Interval != null)
            {
                for (int i = 0; i < this.VCenter_Interval.Count; i++)
                {
                    Point3d p0 = this.VCenterBaseLine.PointAt(this.VCenter_Interval[i].T0);
                    Point3d p1 = this.VCenterBaseLine.PointAt(this.VCenter_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            return verticalLines;
        }

        public void LengthConstraint(double lMin, double lMax, double d, double w)
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
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.South_Interval = new List<Interval>();
                        this.South_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.South_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.South_Interval.Add(new Interval(factors[i], factors[i + 1]));
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
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.North_Interval = new List<Interval>();
                        this.North_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.North_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.North_Interval.Add(new Interval(factors[i], factors[i + 1]));
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
            if (this.HCenterBaseLine != Line.Unset)
            {
                double dw = d + w;
                Line line = new Line(this.HCenterBaseLine.PointAt(0), this.HCenterBaseLine.PointAt(1));
                double length = line.Length;

                if (length > lMax)
                {
                    int lMaxCount = (int)Math.Floor((length - lMax) / (lMax + dw)) + 1;
                    int lMinCount = (int)Math.Floor((length - lMin) / (lMin + dw)) + 1;

                    if (lMinCount == 1 && lMaxCount == 1)
                    {
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.HCenter_Interval = new List<Interval>();
                        this.HCenter_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.HCenter_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.HCenter_Interval.Add(new Interval(factors[i], factors[i + 1]));
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

                        this.HCenter_Interval = new List<Interval>();
                        int iter = 0;
                        while (iter < factors.Count)
                        {
                            this.HCenter_Interval.Add(new Interval(factors[iter], factors[iter + 1]));
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
            factors.Add(0.0);
            for (int i = 0; i < 2 * count - 1; i++)
            {
                double factor = cumulatives[i] / sum;
                factors.Add(factor);
            }
            return factors;
        }
    }
}
