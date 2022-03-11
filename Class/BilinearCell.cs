﻿using Grasshopper.Kernel;
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

        //Polyline polylineForSingleRegion { get; set; }

        // south:0;east:1;north:2;west:3
        //List<int> WhichBoundaryToOffset { get; set; }

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

        public BilinearCell(BilinearCell source)
        {
            #region baseLine
            this.SouthBaseLine = source.SouthBaseLine;
            this.NorthBaseLine = source.NorthBaseLine;
            this.WestBaseLine = source.WestBaseLine;
            this.EastBaseLine = source.EastBaseLine;

            this.HCenterBaseLine = source.HCenterBaseLine;
            this.VCenterBaseLine = source.VCenterBaseLine;
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

            if (source.HCenter_Interval == null)
            {
                this.HCenter_Interval = null;
            }
            else
            {
                this.HCenter_Interval = new List<Interval>();
                this.HCenter_Interval.AddRange(source.HCenter_Interval);
            }
            if (source.VCenter_Interval == null)
            {
                this.VCenter_Interval = null;
            }
            else
            {
                this.VCenter_Interval = new List<Interval>();
                this.VCenter_Interval.AddRange(source.VCenter_Interval);
            }
            #endregion

            this.CenterTValue = source.CenterTValue;
            this.MainVolumeDirection = source.MainVolumeDirection;

            this.ShapeType = source.ShapeType;
            this.BoundaryBaseLineCount = source.BoundaryBaseLineCount;
            this.CenterBaseLineCount = source.CenterBaseLineCount;

            if (source.CellBoundary == null)
            {
                this.CellBoundary = null;
            }
            else
            {
                this.CellBoundary = source.CellBoundary.DuplicateCurve();
            }

            if (source.CellBreps == null)
            {
                this.CellBreps = null;
            }
            else
            {
                this.CellBreps = new Brep[source.CellBreps.Length];
                for (int i = 0; i < source.CellBreps.Length; i++)
                {
                    this.CellBreps[i] = source.CellBreps[i].DuplicateBrep();
                }
            }

            this.PublicBaseLineIndex = source.PublicBaseLineIndex;
            this.AnotherBaseLineIndexRelatedToCutPoint = source.AnotherBaseLineIndexRelatedToCutPoint;
            this.TurningPoint = new Point3d(source.TurningPoint);
        }

        public BilinearCell(Curve southLine, 
                            Curve northLine, 
                            Curve westLine, 
                            Curve eastLine, 

                            Curve southBoundary,
                            Curve northBoundary,
                            Curve westBoundary,
                            Curve eastBoundary,
                            
                            int turningPointIndex)
        {
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

            Polyline poly;
            this.CellBoundary.TryGetPolyline(out poly);
            if (turningPointIndex == 0)
            {
                this.TurningPoint = poly[0];
            }
            else if (turningPointIndex == 1)
            {
                this.TurningPoint = poly[1];
            }
            else if (turningPointIndex == 2)
            {
                this.TurningPoint = poly[2];
            }
            else
            {
                this.TurningPoint = poly[3];
            }

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
        public BilinearCell(Curve boundary, 
                            Curve southLine,
                            Curve northLine, 
                            Curve westLine, 
                            Curve eastLine,
                            int publicBaseLineIndex,
                            int anotherBaseLineIndexRelatedToCutPoint,
                            Point3d turningPoint)
        {
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

                this.TurningPoint = new Point3d(turningPoint);
            }
            else
            {
                this.ShapeType = BilinearShapeType.SingleRegion;

                this.BoundaryBaseLineCount = 4;
                this.CenterBaseLineCount = 0;

                #region CellBoundary
                this.CellBoundary = boundary.DuplicateCurve();
                #endregion

                this.PublicBaseLineIndex = publicBaseLineIndex;
                this.AnotherBaseLineIndexRelatedToCutPoint = anotherBaseLineIndexRelatedToCutPoint;

                this.TurningPoint = new Point3d(turningPoint);
            }
        }


        public BilinearCell(Curve centerH, 
                            Curve westLine, 
                            Curve eastLine, 

                            Curve southBoundary,
                            Curve northBoundary,
                            Curve westBoundary,
                            Curve eastBoundary,

                            int turningPointIndex)
        {
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

            Polyline poly;
            this.CellBoundary.TryGetPolyline(out poly);
            if (turningPointIndex == 0)
            {
                this.TurningPoint = poly[0];
            }
            else if (turningPointIndex == 1)
            {
                this.TurningPoint = poly[1];
            }
            else if (turningPointIndex == 2)
            {
                this.TurningPoint = poly[2];
            }
            else
            {
                this.TurningPoint = poly[3];
            }
        }

        /// <summary>
        /// 生成IShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="isSouthNorth"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, bool isSouthNorth, double t)
        {
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
        private BilinearCell(BilinearCell initialCell, int directionCode)
        {
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
        private BilinearCell(BilinearCell initialCell, BilinearShapeType typeForRec)
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
        private BilinearCell(BilinearCell initialCell ,double t)
        {
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

        public BilinearCell GenerateIAndIVariantShape(double randomT,int directionCode, double w)
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

                newGeneratedCell = new BilinearCell(this, false, t);
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

                    newGeneratedCell = new BilinearCell(this, true, t);
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

                    newGeneratedCell = new BilinearCell(this, false, t);
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateCAndCVariantShape(double randomT,int directionCode, double w)
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
                    newGeneratedCell = new BilinearCell(this, 0);
                }
                else if (this.NorthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 1);
                }
                else if (this.WestBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 2);
                }
                else
                {
                    newGeneratedCell = new BilinearCell(this, 3);
                }
            }
            else
            {
                switch (directionCode)
                {
                    case 0:
                        newGeneratedCell = new BilinearCell(this, 0);
                        break;
                    case 1:
                        newGeneratedCell = new BilinearCell(this, 1);
                        break;
                    case 2:
                        newGeneratedCell = new BilinearCell(this, 2);
                        break;
                    case 3:
                        newGeneratedCell = new BilinearCell(this, 3);
                        break;
                    default:
                        newGeneratedCell = new BilinearCell(this, 0);
                        break;
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateRecAndRecVariantShape(double randomT, int directionCode, double lMin, double w, bool flag)
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
                    newGeneratedCell = new BilinearCell(this, BilinearShapeType.RecShape);
                }
                else
                {
                    if (flag)
                    {
                        // 执行GenerateRecShape函数
                        newGeneratedCell = new BilinearCell(this, BilinearShapeType.RecShape);
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

                            newGeneratedCell = new BilinearCell(this, t);
                        }
                        else
                        {
                            // 执行GenerateRecShape函数
                            newGeneratedCell = new BilinearCell(this, BilinearShapeType.RecShape);
                        }
                    }
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateRemovedCShapeOrIShapeOrSinglePolyline()
        {
            if (this.ShapeType != BilinearShapeType.SinglePolyline 
                || this.ShapeType != BilinearShapeType.IShape
                || this.ShapeType != BilinearShapeType.CShape)
            {
                return this;
            }
            else
            {
                if (this.ShapeType == BilinearShapeType.SinglePolyline)
                {
                    // 去掉非连接的那一条
                    if (AnotherBaseLineIndexRelatedToCutPoint == 0)
                    {
                        this.South_Interval = new List<Interval>();
                    }
                    else if (AnotherBaseLineIndexRelatedToCutPoint == 1)
                    {
                        this.East_Interval = new List<Interval>();
                    }
                    else if (AnotherBaseLineIndexRelatedToCutPoint == 2)
                    {
                        this.North_Interval = new List<Interval>();
                    }
                    else
                    {
                        this.West_Interval = new List<Interval>();
                    }
                }
                else if (this.ShapeType == BilinearShapeType.CShape)
                {
                    // 去掉中间的那一条
                    if (this.South_Interval.Count == 0)
                    {
                        this.North_Interval = new List<Interval>();
                    }
                    else if (this.East_Interval.Count == 0)
                    {
                        this.West_Interval = new List<Interval>();
                    }
                    else if (this.North_Interval.Count == 0)
                    {
                        this.South_Interval = new List<Interval>();
                    }
                    else
                    {
                        this.East_Interval = new List<Interval>();
                    }
                }
                else
                {
                    // 去掉中间的那一条
                    if (this.HCenter_Interval.Count == 1)
                    {
                        this.HCenter_Interval = new List<Interval>();
                    }
                    if (this.VCenter_Interval.Count == 1)
                    {
                        this.VCenter_Interval = new List<Interval>();
                    }
                }
                return this;
            }
        }

        public BilinearCell GenerateOffsetedSingleRegion(double scaleFactor)
        {
            if (this.ShapeType != BilinearShapeType.SingleRegion)
            {
                return this;
            }
            else
            {

                //if (this.PublicBaseLineIndex == -1 && this.AnotherBaseLineIndexRelatedToCutPoint == -1)
                //{
                //    this.PublicBaseLineIndex = randomIndex;
                //    this.AnotherBaseLineIndexRelatedToCutPoint = (randomIndex + 1) % 4;
                //}

                //int index0ForBLPublic = this.PublicBaseLineIndex;
                //int index1ForBLPublic = this.AnotherBaseLineIndexRelatedToCutPoint;
                //int index0ForBLNeedToOffset = (this.PublicBaseLineIndex + 2) % 4;
                //int index1ForBLNeedToOffset = (this.AnotherBaseLineIndexRelatedToCutPoint + 2) % 4;

                //#region 对于index0ForBLPublic
                //Line line0Public = Line.Unset;
                //if (index0ForBLPublic == 0)
                //{
                //    line0Public = this.SouthBaseLine;
                //}
                //else if (index0ForBLPublic == 1)
                //{
                //    line0Public = this.EastBaseLine;
                //}
                //else if (index0ForBLPublic == 2)
                //{
                //    line0Public = this.NorthBaseLine;
                //}
                //else
                //{
                //    line0Public = this.WestBaseLine;
                //}
                //#endregion

                //#region 对于index1ForBLPublic
                //Line line1Public = Line.Unset;
                //if (index1ForBLPublic == 0)
                //{
                //    line1Public = this.SouthBaseLine;
                //}
                //else if (index1ForBLPublic == 1)
                //{
                //    line1Public = this.EastBaseLine;
                //}
                //else if (index1ForBLPublic == 2)
                //{
                //    line1Public = this.NorthBaseLine;
                //}
                //else
                //{
                //    line1Public = this.WestBaseLine;
                //}
                //#endregion

                //#region 对于index0ForBLNeedToOffset
                //Line line0 = Line.Unset;
                //if (index0ForBLNeedToOffset == 0)
                //{
                //    line0 = this.SouthBaseLine;
                //}
                //else if (index0ForBLNeedToOffset == 1)
                //{
                //    line0 = this.EastBaseLine;
                //}
                //else if (index0ForBLNeedToOffset == 2)
                //{
                //    line0 = this.NorthBaseLine;
                //}
                //else
                //{
                //    line0 = this.WestBaseLine;
                //}
                //#endregion

                //#region 对于index1ForBLNeedToOffset
                //Line line1 = Line.Unset;
                //if (index1ForBLNeedToOffset == 0)
                //{
                //    line1 = this.SouthBaseLine;
                //}
                //else if (index1ForBLNeedToOffset == 1)
                //{
                //    line1 = this.EastBaseLine;
                //}
                //else if (index1ForBLNeedToOffset == 2)
                //{
                //    line1 = this.NorthBaseLine;
                //}
                //else
                //{
                //    line1 = this.WestBaseLine;
                //}
                //#endregion

                //#region offsetedIndex0
                //Vector3d vec = line0.Direction;
                //Vector3d verticalVec = new Vector3d(-vec.Y, vec.X, vec.Z);
                //verticalVec.Unitize();
                //verticalVec *= offset0;
                //Line offsetdLine0 = new Line(line0.From + verticalVec, line0.To + verticalVec);
                //#endregion
                //#region offsetedIndex1
                //vec = line1.Direction;
                //verticalVec = new Vector3d(-vec.Y, vec.X, vec.Z);
                //verticalVec.Unitize();
                //verticalVec *= offset1;
                //#endregion

                Transform transform = Transform.Scale(this.TurningPoint, scaleFactor);
                Curve crv = this.CellBoundary.DuplicateCurve();
                crv.Transform(transform);
                Curve[] segments = crv.DuplicateSegments();
                this.SouthBaseLine = new Line(segments[0].PointAtStart, segments[0].PointAtEnd);
                this.EastBaseLine = new Line(segments[1].PointAtStart, segments[1].PointAtEnd);
                this.NorthBaseLine = new Line(segments[2].PointAtEnd, segments[2].PointAtStart);
                this.WestBaseLine = new Line(segments[3].PointAtEnd, segments[3].PointAtStart);

                return this;
            }
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
