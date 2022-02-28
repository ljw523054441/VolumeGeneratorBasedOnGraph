using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
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

        public List<Interval> HCenter_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> VCenter_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };

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
            // 奇数排的形式种类
            CShape = 1,
            IShape = 2,
            //TShape = 3,
            RecShape = 3
            // 偶数排的形式种类
            //NullShape = 4,
            //LineShape = 5,
            //ParallelShape = 6,
            // 单个体量的形式种类
            //SingleRegion = 7
        }

        public BilinearShapeType ShapeType { get; set; }

        public BilinearCell(Curve southLine, Curve northLine, Curve westLine, Curve eastLine, int i, int j)
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
        }

        public BilinearCell(Curve centerH, Curve westLine, Curve eastLine, int i, int j)
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
        }
    
        //private BilinearCell(BilinearCell initialCell, int i ,int j)
        //{
        //    this.SouthBaseLine = initialCell.SouthBaseLine;
        //    this.NorthBaseLine = initialCell.NorthBaseLine;
        //    this.WestBaseLine = initialCell.WestBaseLine;
        //    this.EastBaseLine = initialCell.EastBaseLine;
        //    this.HCenterBaseLine = initialCell.HCenterBaseLine;
        //    this.VCenterBaseLine = initialCell.VCenterBaseLine;

        //    this.CellShapeType = ShapeType.NullShape;
        //}

        /// <summary>
        /// 生成IShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="isSouthNorth"></param>
        /// <param name="centerLine"></param>
        /// <param name="line0"></param>
        /// <param name="line1"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, bool isSouthNorth, Line centerLine, Line line0, Line line1, double t, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
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

            #region Lines
            if (isSouthNorth)
            {
                //// 父类属性
                //this.SouthLines = new List<Line>();
                //this.SouthLines.Add(new Line(line0.From, line0.To));
                //this.NorthLines = new List<Line>();
                //this.NorthLines.Add(new Line(line1.From, line1.To));
                //// 子类属性
                //this.VLines = new List<Line>();
                //this.VLines.Add(new Line(centerLine.From, centerLine.To));

                this.South_Interval = new List<Interval>();
                this.South_Interval.Add(new Interval(0, 1));
                this.North_Interval = new List<Interval>();
                this.North_Interval.Add(new Interval(0, 1));

                this.VCenter_Interval = new List<Interval>();
                this.VCenter_Interval.Add(new Interval(0, 1));
            }
            else
            {
                //// 父类属性
                //this.WestLines = new List<Line>();
                //this.WestLines.Add(new Line(line0.From, line0.To));
                //this.EastLines = new List<Line>();
                //this.EastLines.Add(new Line(line1.From, line1.To));
                //// 子类属性
                //this.HLines = new List<Line>();
                //this.HLines.Add(new Line(centerLine.From, centerLine.To));

                this.West_Interval = new List<Interval>();
                this.West_Interval.Add(new Interval(0, 1));
                this.East_Interval = new List<Interval>();
                this.East_Interval.Add(new Interval(0, 1));

                this.HCenter_Interval = new List<Interval>();
                this.HCenter_Interval.Add(new Interval(0, 1));
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
        }

        /// <summary>
        /// 生成CShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="columnLine"></param>
        /// <param name="line0"></param>
        /// <param name="line1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, int directionCode, Line columnLine, Line line0, Line line1, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
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

            #region Lines
            switch (directionCode)
            {
                case 0:
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(columnLine.From, columnLine.To));
                    //this.WestLines = new List<Line>();
                    //this.WestLines.Add(new Line(line0.From, line0.To));
                    //this.EastLines = new List<Line>();
                    //this.EastLines.Add(new Line(line1.From, line1.To));

                    this.South_Interval = new List<Interval>();
                    this.South_Interval.Add(new Interval(0, 1));
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 0;
                    break;
                case 1:
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(columnLine.From, columnLine.To));
                    //this.WestLines = new List<Line>();
                    //this.WestLines.Add(new Line(line0.From, line0.To));
                    //this.EastLines = new List<Line>();
                    //this.EastLines.Add(new Line(line1.From, line1.To));

                    this.North_Interval = new List<Interval>();
                    this.North_Interval.Add(new Interval(0, 1));
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval());
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval());

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 1;
                    break;
                case 2:
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(line0.From, line0.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(line1.From, line1.To));
                    //this.WestLines = new List<Line>();
                    //this.WestLines.Add(new Line(columnLine.From, columnLine.To));

                    this.South_Interval = new List<Interval>();
                    this.South_Interval.Add(new Interval(0, 1));
                    this.North_Interval = new List<Interval>();
                    this.North_Interval.Add(new Interval(0, 1));
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 0;
                    break;
                case 3:
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(line0.From, line0.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(line1.From, line1.To));
                    //this.EastLines = new List<Line>();
                    //this.EastLines.Add(new Line(columnLine.From, columnLine.To));

                    this.South_Interval = new List<Interval>();
                    this.South_Interval.Add(new Interval(0, 1));
                    this.North_Interval = new List<Interval>();
                    this.North_Interval.Add(new Interval(0, 1));
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 1;
                    break;
            }
            #endregion

            this.ShapeType = BilinearShapeType.CShape;

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();
        }

        /// <summary>
        /// 生成RecShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="south"></param>
        /// <param name="north"></param>
        /// <param name="west"></param>
        /// <param name="east"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private BilinearCell(BilinearCell initialCell, Line south, Line north, Line west, Line east, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region baseLine
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;

            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;
            #endregion

            #region Lines
            //this.SouthLines = new List<Line>();
            //this.SouthLines.Add(new Line(south.From, south.To));
            //this.NorthLines = new List<Line>();
            //this.NorthLines.Add(new Line(north.From, north.To));
            //this.WestLines = new List<Line>();
            //this.WestLines.Add(new Line(west.From, west.To));
            //this.EastLines = new List<Line>();
            //this.EastLines.Add(new Line(east.From, east.To));

            this.South_Interval = new List<Interval>();
            this.South_Interval.Add(new Interval(0, 1));
            this.North_Interval = new List<Interval>();
            this.North_Interval.Add(new Interval(0, 1));
            this.West_Interval = new List<Interval>();
            this.West_Interval.Add(new Interval(0, 1));
            this.East_Interval = new List<Interval>();
            this.East_Interval.Add(new Interval(0, 1));
            #endregion

            this.MainVolumeDirection = MainDirection.Unset;
            this.CenterTValue = 0;

            this.ShapeType = BilinearShapeType.RecShape;

            this.BoundaryBaseLineCount = 4;
            this.CenterBaseLineCount = 0;
        }

        public BilinearCell GenerateIShape(int randomT,int directionCode)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                Line west = this.WestBaseLine;
                Line east = this.EastBaseLine;

                newGeneratedCell = new BilinearCell(this, false, this.HCenterBaseLine, west, east, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCShape(randomT, directionCode);
            }
            else
            {
                if (directionCode >= 2)
                {
                    // 南北体量为主
                    Line south = this.SouthBaseLine;
                    Line north = this.NorthBaseLine;
                    Point3d pointOnSouth = south.PointAt(0.2 + randomT * 0.2);
                    Point3d pointOnNorth = north.PointAt(0.2 + randomT * 0.2);
                    Line newVCenter = new Line(pointOnSouth, pointOnNorth);

                    newGeneratedCell = new BilinearCell(this, true, newVCenter, south, north, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    // 东西体量为主
                    Line west = this.WestBaseLine;
                    Line east = this.EastBaseLine;
                    Point3d pointOnWest = west.PointAt(0.2 + randomT * 0.2);
                    Point3d pointOnEast = east.PointAt(0.2 + randomT * 0.2);
                    Line newHCenter = new Line(pointOnWest, pointOnEast);

                    newGeneratedCell = new BilinearCell(this, false, newHCenter, west, east, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateCShape(int randomT,int directionCode)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIShape(randomT, directionCode);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                if (this.SouthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.NorthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 1, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.WestBaseLine != Line.Unset)
                {
                    newGeneratedCell = new BilinearCell(this, 2, this.WestBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    newGeneratedCell = new BilinearCell(this, 3, this.EastBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            else
            {
                switch (directionCode)
                {
                    case 0:
                        newGeneratedCell = new BilinearCell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 1:
                        newGeneratedCell = new BilinearCell(this, 1, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 2:
                        newGeneratedCell = new BilinearCell(this, 2, this.WestBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 3:
                        newGeneratedCell = new BilinearCell(this, 3, this.EastBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    default:
                        newGeneratedCell = new BilinearCell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                }
            }
            return newGeneratedCell;
        }

        public BilinearCell GenerateRecShape(int randomT, int directionCode)
        {
            BilinearCell newGeneratedCell;
            if (this.ShapeType == BilinearShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIShape(randomT, directionCode);
            }
            else if (this.ShapeType == BilinearShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCShape(randomT, directionCode);
            }
            else
            {
                newGeneratedCell = new BilinearCell(this, this.SouthBaseLine, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
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
            if (this.SouthBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.South_Interval.Count; i++)
                {
                    Point3d p0 = this.SouthBaseLine.PointAt(this.South_Interval[i].T0);
                    Point3d p1 = this.SouthBaseLine.PointAt(this.South_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.NorthBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.North_Interval.Count; i++)
                {
                    Point3d p0 = this.NorthBaseLine.PointAt(this.North_Interval[i].T0);
                    Point3d p1 = this.NorthBaseLine.PointAt(this.North_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.HCenterBaseLine != Line.Unset)
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
            if (this.WestBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.West_Interval.Count; i++)
                {
                    Point3d p0 = this.WestBaseLine.PointAt(this.West_Interval[i].T0);
                    Point3d p1 = this.WestBaseLine.PointAt(this.West_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.EastBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.East_Interval.Count; i++)
                {
                    Point3d p0 = this.EastBaseLine.PointAt(this.East_Interval[i].T0);
                    Point3d p1 = this.EastBaseLine.PointAt(this.East_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.VCenterBaseLine != Line.Unset)
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
    }
}
