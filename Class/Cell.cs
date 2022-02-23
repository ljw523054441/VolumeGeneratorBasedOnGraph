using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    class Cell
    {
        public enum MainDirection
        {
            Unset = 0,
            SouthNorth = 1,
            WestEast = 2
        }

        public enum ShapeType
        {
            Unset = 0,
            // 奇数排的形式种类
            CShape = 1,
            IShape = 2,
            //TShape = 3,
            RecShape = 3,
            // 偶数排的形式种类
            NullShape = 4,
            LineShape = 5,
            ParallelShape = 6,
            // 单个体量的形式种类
            SingleRegion = 7
        }

        public MainDirection MainVolumeDirection { get; set; }

        public ShapeType CellShapeType { get; set; }
        
        public Line SouthBaseLine { get; set; }
        public Line NorthBaseLine { get; set; }
        public Line EastBaseLine { get; set; }
        public Line WestBaseLine { get; set; }
        public Line HCenterBaseLine { get; set; }
        public Line VCenterBaseLine { get; set; }

        public List<Line> SouthLines { get; set; }
        public List<Line> NorthLines { get; set; }
        public List<Line> EastLines { get; set; }
        public List<Line> WestLines { get; set; }
        public List<Line> HLines { get; set; }
        public List<Line> VLines { get; set; }

        public int BoundaryBaseLineCount { get; set; }
        public int CenterBaseLineCount { get; set; }

        public int[] CurrentIndex { get; set; }
        public int[] PrevSouthIndex { get; set; }
        public int[] PrevWestIndex { get; set; }
        public int[] NexEastIndex { get; set; }
        public int[] NexNorthIndex { get; set; }

        public double CenterTValue { get; set; }


        /// <summary>
        /// 生成SingleRegion
        /// </summary>
        /// <param name="polyline"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Curve polyline, int i, int j)
        {
            // todo
        }

        public Cell(Curve southLine, Curve northLine, Curve westLine, Curve eastLine, int i, int j)
        {
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

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            if (count == 4 )
            {
                this.CellShapeType = ShapeType.RecShape;
            }
            else if (count == 3)
            {
                this.CellShapeType = ShapeType.CShape;
            }

            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = 0;
        }

        public Cell(Curve centerH, Curve westLine, Curve eastLine, int i, int j)
        {
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

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.CellShapeType = ShapeType.IShape;
            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = count;
        }


        /// <summary>
        /// 生成NullShape
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, int i, int j)
        {
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;

            this.CellShapeType = ShapeType.NullShape;
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
        }

        /// <summary>
        /// 生成IShape
        /// </summary>
        /// <param name="isSouthNorth"></param>
        /// <param name="centerLine"></param>
        /// <param name="line0"></param>
        /// <param name="line1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, bool isSouthNorth, Line centerLine, Line line0, Line line1, double t, int i, int j)
        {
            if (isSouthNorth)
            {
                this.SouthBaseLine = initialCell.SouthBaseLine;
                this.NorthBaseLine = initialCell.NorthBaseLine;
                this.WestBaseLine = initialCell.WestBaseLine;
                this.EastBaseLine = initialCell.EastBaseLine;
                this.HCenterBaseLine = initialCell.HCenterBaseLine;
                this.VCenterBaseLine = initialCell.VCenterBaseLine;
                this.MainVolumeDirection = MainDirection.SouthNorth;
                this.CenterTValue = t;

                this.SouthLines = new List<Line>();
                this.SouthLines.Add(new Line(line0.From, line0.To));
                this.NorthLines = new List<Line>();
                this.NorthLines.Add(new Line(line1.From, line1.To));
                this.VLines = new List<Line>();
                this.VLines.Add(new Line(centerLine.From, centerLine.To));
            }
            else
            {
                this.SouthBaseLine = initialCell.SouthBaseLine;
                this.NorthBaseLine = initialCell.NorthBaseLine;
                this.WestBaseLine = initialCell.WestBaseLine;
                this.EastBaseLine = initialCell.EastBaseLine;
                this.HCenterBaseLine = initialCell.HCenterBaseLine;
                this.VCenterBaseLine = initialCell.VCenterBaseLine;

                this.WestLines = new List<Line>();
                this.WestLines.Add(new Line(line0.From, line0.To));
                this.EastLines = new List<Line>();
                this.EastLines.Add(new Line(line1.From, line1.To));
                this.HLines = new List<Line>();
                this.HLines.Add(new Line(centerLine.From, centerLine.To));

                this.MainVolumeDirection = MainDirection.WestEast;
                this.CenterTValue = t;
            }
            this.CellShapeType = ShapeType.IShape;

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();
        }

        /// <summary>
        /// 生成CShape
        /// </summary>
        /// <param name="directionCode">0:S; 1:N; 2:W; 3:E;</param>
        /// <param name="columnLine"></param>
        /// <param name="line0"></param>
        /// <param name="line1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, int directionCode, Line columnLine, Line line0, Line line1,int i, int j)
        {
            switch (directionCode)
            {
                case 0:
                    this.SouthBaseLine = initialCell.SouthBaseLine;
                    this.NorthBaseLine = initialCell.NorthBaseLine;
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.HCenterBaseLine = initialCell.HCenterBaseLine;
                    this.VCenterBaseLine = initialCell.VCenterBaseLine;

                    this.SouthLines = new List<Line>();
                    this.SouthLines.Add(new Line(columnLine.From, columnLine.To));
                    this.WestLines = new List<Line>();
                    this.WestLines.Add(new Line(line0.From, line0.To));
                    this.EastLines = new List<Line>();
                    this.EastLines.Add(new Line(line1.From, line1.To));

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 0;
                    break;
                case 1:
                    this.SouthBaseLine = initialCell.SouthBaseLine;
                    this.NorthBaseLine = initialCell.NorthBaseLine;
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.HCenterBaseLine = initialCell.HCenterBaseLine;
                    this.VCenterBaseLine = initialCell.VCenterBaseLine;

                    this.NorthLines = new List<Line>();
                    this.NorthLines.Add(new Line(columnLine.From, columnLine.To));
                    this.WestLines = new List<Line>();
                    this.WestLines.Add(new Line(line0.From, line0.To));
                    this.EastLines = new List<Line>();
                    this.EastLines.Add(new Line(line1.From, line1.To));

                    this.MainVolumeDirection = MainDirection.WestEast;
                    this.CenterTValue = 1;
                    break;
                case 2:
                    this.SouthBaseLine = initialCell.SouthBaseLine;
                    this.NorthBaseLine = initialCell.NorthBaseLine;
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.HCenterBaseLine = initialCell.HCenterBaseLine;
                    this.VCenterBaseLine = initialCell.VCenterBaseLine;

                    this.SouthLines = new List<Line>();
                    this.SouthLines.Add(new Line(line0.From, line0.To));
                    this.NorthLines = new List<Line>();
                    this.NorthLines.Add(new Line(line1.From, line1.To));
                    this.WestLines = new List<Line>();
                    this.WestLines.Add(new Line(columnLine.From, columnLine.To));

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 0;
                    break;
                case 3:
                    this.SouthBaseLine = initialCell.SouthBaseLine;
                    this.NorthBaseLine = initialCell.NorthBaseLine;
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.HCenterBaseLine = initialCell.HCenterBaseLine;
                    this.VCenterBaseLine = initialCell.VCenterBaseLine;

                    this.SouthLines = new List<Line>();
                    this.SouthLines.Add(new Line(line0.From, line0.To));
                    this.NorthLines = new List<Line>();
                    this.NorthLines.Add(new Line(line1.From, line1.To));
                    this.EastLines = new List<Line>();
                    this.EastLines.Add(new Line(columnLine.From, columnLine.To));

                    this.MainVolumeDirection = MainDirection.SouthNorth;
                    this.CenterTValue = 1;
                    break;
            }
            this.CellShapeType = ShapeType.CShape;

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();
        }

        /// <summary>
        /// 生成RecShape
        /// </summary>
        /// <param name="south"></param>
        /// <param name="north"></param>
        /// <param name="west"></param>
        /// <param name="east"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, Line south, Line north, Line west, Line east, int i, int j)
        {
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;

            this.SouthLines = new List<Line>();
            this.SouthLines.Add(new Line(south.From, south.To));
            this.NorthLines = new List<Line>();
            this.NorthLines.Add(new Line(north.From, north.To));
            this.WestLines = new List<Line>();
            this.WestLines.Add(new Line(west.From, west.To));
            this.EastLines = new List<Line>();
            this.EastLines.Add(new Line(east.From, east.To));

            this.MainVolumeDirection = MainDirection.Unset;
            this.CenterTValue = 0;
            this.CellShapeType = ShapeType.RecShape;

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.BoundaryBaseLineCount = 4;
            this.CenterBaseLineCount = 0;
        }

        /// <summary>
        /// 生成ParallelShape
        /// </summary>
        /// <param name="line0"></param>
        /// <param name="line1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, Line line0, Line line1, int i, int j)
        {
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;

            this.WestLines = new List<Line>();
            this.WestLines.Add(new Line(line0.From, line0.To));
            this.EastLines = new List<Line>();
            this.EastLines.Add(new Line(line1.From, line1.To));

            this.MainVolumeDirection = MainDirection.WestEast;
            this.CenterTValue = 0;

            this.CellShapeType = ShapeType.ParallelShape;
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();
        }

        /// <summary>
        /// 生成LineShape
        /// </summary>
        /// <param name="line"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        public Cell(Cell initialCell, Line line, double t, int i, int j)
        {
            this.SouthBaseLine = initialCell.SouthBaseLine;
            this.NorthBaseLine = initialCell.NorthBaseLine;
            this.WestBaseLine = initialCell.WestBaseLine;
            this.EastBaseLine = initialCell.EastBaseLine;
            this.HCenterBaseLine = initialCell.HCenterBaseLine;
            this.VCenterBaseLine = initialCell.VCenterBaseLine;

            if (t == 0)
            {
                this.WestLines = new List<Line>();
                this.WestLines.Add(new Line(line.From, line.To));
                this.CenterTValue = 0;
            }
            else if (t == 1)
            {
                this.EastLines = new List<Line>();
                this.EastLines.Add(new Line(line.From, line.To));
                this.CenterTValue = 1;
            }
            else
            {
                this.VLines = new List<Line>();
                this.VLines.Add(new Line(line.From, line.To));
                this.CenterTValue = t;
            }
            this.MainVolumeDirection = MainDirection.WestEast;
            this.CellShapeType = ShapeType.LineShape;
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };

            this.BoundaryBaseLineCount = GetBoundaryBaseLinesCount();
            this.CenterBaseLineCount = GetCenterBaseLineCount();
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

        public List<Curve> GetAllLines()
        {
            List<Curve> allLines = new List<Curve>();
            if (this.SouthBaseLine != Line.Unset)
            {
                allLines.Add(this.SouthBaseLine.ToNurbsCurve());
            }
            if (this.NorthBaseLine != Line.Unset)
            {
                allLines.Add(this.NorthBaseLine.ToNurbsCurve());
            }
            if (this.WestBaseLine != Line.Unset)
            {
                allLines.Add(this.WestBaseLine.ToNurbsCurve());
            }
            if (this.EastBaseLine != Line.Unset)
            {
                allLines.Add(this.EastBaseLine.ToNurbsCurve());
            }
            if (this.HCenterBaseLine != Line.Unset)
            {
                allLines.Add(this.HCenterBaseLine.ToNurbsCurve());
            }
            if (this.VCenterBaseLine != Line.Unset)
            {
                allLines.Add(this.VCenterBaseLine.ToNurbsCurve());
            }

            return allLines;
        }

        public Cell GenerateIShape(int randomT, int directionCode)
        {
            Cell newGeneratedCell;
            if (this.CellShapeType == ShapeType.IShape)
            {
                Line west = this.WestBaseLine;
                Line east = this.EastBaseLine;

                newGeneratedCell = new Cell(this, false, this.HCenterBaseLine, west, east, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
            }
            else if (this.CellShapeType == ShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCShape(randomT, directionCode);
            }
            else
            {
                // 对于偶数排的cell，偶数排才会生成IShape
                if (directionCode >= 2)
                {
                    // 南北体量为主
                    Line south = this.SouthBaseLine;
                    Line north = this.NorthBaseLine;
                    Point3d pointOnSouth = south.PointAt(0.2 + randomT * 0.2);
                    Point3d pointOnNorth = north.PointAt(0.2 + randomT * 0.2);
                    Line newVCenter = new Line(pointOnSouth, pointOnNorth);

                    newGeneratedCell = new Cell(this, true, newVCenter, south, north, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    // 东西体量为主
                    Line west = this.WestBaseLine;
                    Line east = this.EastBaseLine;
                    Point3d pointOnWest = west.PointAt(0.2 + randomT * 0.2);
                    Point3d pointOnEast = east.PointAt(0.2 + randomT * 0.2);
                    Line newHCenter = new Line(pointOnWest, pointOnEast);

                    newGeneratedCell = new Cell(this, false, newHCenter, west, east, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            return newGeneratedCell;
        }

        public Cell GenerateCShape(int randomT, int directionCode)
        {
            Cell newGeneratedCell;
            if (this.CellShapeType == ShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIShape(randomT, directionCode);
            }
            else if (this.CellShapeType == ShapeType.CShape)
            {
                if (this.SouthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new Cell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.NorthBaseLine != Line.Unset)
                {
                    newGeneratedCell = new Cell(this, 1, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else if (this.WestBaseLine != Line.Unset)
                {
                    newGeneratedCell = new Cell(this, 2, this.WestBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
                else
                {
                    newGeneratedCell = new Cell(this, 3, this.EastBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                }
            }
            else
            {
                // 对于偶数排的cell，偶数排才会生成CShape
                switch (directionCode)
                {
                    case 0:
                        newGeneratedCell = new Cell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 1:
                        newGeneratedCell = new Cell(this, 1, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 2:
                        newGeneratedCell = new Cell(this, 2, this.WestBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    case 3:
                        newGeneratedCell = new Cell(this, 3, this.EastBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                    default:
                        newGeneratedCell = new Cell(this, 0, this.SouthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                        break;
                }
            }
            return newGeneratedCell;
        }

        public Cell GenerateRecShape(int randomT, int directionCode)
        {
            Cell newGeneratedCell;
            if (this.CellShapeType == ShapeType.IShape)
            {
                // 执行GenerateIShape函数
                newGeneratedCell = GenerateIShape(randomT, directionCode);
            }
            else if (this.CellShapeType == ShapeType.CShape)
            {
                // 执行GenerateCShape函数
                newGeneratedCell = GenerateCShape(randomT, directionCode);
            }
            else
            {
                newGeneratedCell = new Cell(this, this.SouthBaseLine, this.NorthBaseLine, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
            }
            return newGeneratedCell;
        }

        public Cell GenerateLineShapeAndVariants(int randomT, int directionCode)
        {
            Cell newGeneratedCell;
            // 对于奇数排的cell
            switch (directionCode)
            {
                case 0:
                    // 生成NullShape
                    newGeneratedCell = new Cell(this, this.CurrentIndex[0], this.CurrentIndex[1]);
                    break;
                case 1:
                    // 生成左侧或右侧的LineShape
                    Line line;
                    double t;
                    if (randomT < 2)
                    {
                        line = new Line(this.WestBaseLine.From, this.WestBaseLine.To);
                        t = 0;
                    }
                    else
                    {
                        line = new Line(this.EastBaseLine.From, this.EastBaseLine.To);
                        t = 1;
                    }
                    newGeneratedCell = new Cell(this, line, t, this.CurrentIndex[0], this.CurrentIndex[1]);
                    break;
                case 2:
                    // 生成中间的LineShape
                    Line south = this.SouthBaseLine;
                    Line north = this.NorthBaseLine;
                    Point3d pointOnSouth = south.PointAt(0.2 + randomT * 0.2);
                    Point3d pointOnNorth = north.PointAt(0.2 + randomT * 0.2);
                    Line newVCenter = new Line(pointOnSouth, pointOnNorth);
                    newGeneratedCell = new Cell(this, newVCenter, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
                    break;
                case 3:
                    // 生成ParallelShape
                    newGeneratedCell = new Cell(this, this.WestBaseLine, this.EastBaseLine, this.CurrentIndex[0], this.CurrentIndex[1]);
                    break;
                default:
                    // 生成NullShape
                    newGeneratedCell = new Cell(this, this.CurrentIndex[0], this.CurrentIndex[1]);
                    break;
            }

            return newGeneratedCell;
        }

        public bool ConstraintLength(int random, double lMax, double lMin, double d)
        {
            bool flag0 = false;
            if (this.SouthLines != null)
            {
                if (this.SouthLines[0].Length > lMax)
                {
                    double length = this.SouthLines[0].Length;
                    int volumeLMinCount = (int)Math.Floor((length - lMin) % (lMin + d)) + 1;
                    int volumeLMaxCount = (int)Math.Floor((length - lMax) % (lMax + d)) + 1;

                    if (volumeLMinCount == 1 && volumeLMaxCount == 1)
                    {

                        if (this.CellShapeType == ShapeType.CShape)
                        {
                            // CShape时
                            if (this.NorthLines == null)
                            {
                                if (random == 0)
                                {
                                    //Point3d start = this.SouthLines[0].From;
                                    double t = (length - lMax) / length;
                                    Line line = this.SouthLines[0];
                                    Curve crv = line.ToNurbsCurve();
                                    crv = crv.Trim(t, 1.0);
                                    this.SouthLines.Clear();
                                    this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                                }
                                else
                                {
                                    double t = (length -(length - lMax)) / length;
                                    Line line = this.SouthLines[0];
                                    Curve crv = line.ToNurbsCurve();
                                    crv = crv.Trim(0.0, t);
                                    this.SouthLines.Clear();
                                    this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                                }
                            }
                            else if (this.WestLines == null)
                            {
                                double t = (length - lMax) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(t, 1.0);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                            else if (this.EastLines == null)
                            {
                                double t = (length - (length - lMax)) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(0.0, t);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                        }
                        else if (this.CellShapeType == ShapeType.IShape)
                        {
                            // IShape时
                            if (random == 0)
                            {
                                //Point3d start = this.SouthLines[0].From;
                                double t = (length - lMax) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(t, 1.0);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                            else
                            {
                                double t = (length - (length - lMax)) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(0.0, t);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                        }
                        else
                        {
                            // RecShape时
                            if (random == 0)
                            {
                                //Point3d start = this.SouthLines[0].From;
                                double t = (length - lMax) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(t, 1.0);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                            else
                            {
                                double t = (length - (length - lMax)) / length;
                                Line line = this.SouthLines[0];
                                Curve crv = line.ToNurbsCurve();
                                crv = crv.Trim(0.0, t);
                                this.SouthLines.Clear();
                                this.SouthLines.Add(new Line(crv.PointAtStart, crv.PointAtEnd));
                            }
                        }
                    }
                    else
                    {
                        // 按照 volumeLMinCount 和 volumeLMaxCount 中较大的那个做切割
                        int cutCount = volumeLMinCount > volumeLMaxCount ? volumeLMinCount : volumeLMaxCount;
                        // 特别注意要处理IShape时的切割
                        if (this.CellShapeType == ShapeType.CShape)
                        {
                            // CShape时
                            if (this.NorthLines == null)
                            {

                            }
                            else if (this.WestLines == null)
                            {

                            }
                            else if (this.EastLines == null)
                            {

                            }
                        }
                        else if (this.CellShapeType == ShapeType.IShape)
                        {

                        }
                        else 
                        {
                            // RecShape时
                        }
                    }

                    flag0 = true;
                }
            }

            bool flag1 = false;
            if (this.NorthLines != null)
            {
                
                
                flag1 = true;
            }

            if (flag0 == true || flag1 == true)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public void GetAllLines(out List<Curve> hCurves, out List<Curve> vCurves)
        {
            hCurves = new List<Curve>();
            vCurves = new List<Curve>();
            if (this.SouthLines != null)
            {
                for (int i = 0; i < this.SouthLines.Count; i++)
                {
                    hCurves.Add(this.SouthLines[i].ToNurbsCurve());
                }
            }
            if (this.NorthLines != null)
            {
                for (int i = 0; i < this.NorthLines.Count; i++)
                {
                    hCurves.Add(this.NorthLines[i].ToNurbsCurve());
                }
            }
            if (this.HLines != null)
            {
                for (int i = 0; i < this.HLines.Count; i++)
                {
                    hCurves.Add(this.HLines[i].ToNurbsCurve());
                }
            }
            if (this.WestLines != null)
            {
                for (int i = 0; i < this.WestLines.Count; i++)
                {
                    vCurves.Add(this.WestLines[i].ToNurbsCurve());
                }
            }
            if (this.EastLines != null)
            {
                for (int i = 0; i < this.EastLines.Count; i++)
                {
                    vCurves.Add(this.EastLines[i].ToNurbsCurve());
                }
            }
            if (this.VLines != null)
            {
                for (int i = 0; i < this.VLines.Count; i++)
                {
                    vCurves.Add(this.VLines[i].ToNurbsCurve());
                }
            }

        }
    }
}
