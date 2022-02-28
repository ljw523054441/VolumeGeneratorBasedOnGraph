using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
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

        //public List<Line> MiddleLines { get; set; }
        //public List<Line> WestLines1 { get; set; }
        //public List<Line> EastLines1 { get; set; }

        //public List<Line> HLinesSouth { get; set; }
        //public List<Line> HLinesNorth { get; set; }

        //public List<Line> VLinesSouth { get; set; }
        //public List<Line> VLinesNorth { get; set; }

        public List<Interval> Middle_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> West1_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> East1_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> HSouth_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> HNouth_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> VSouth_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };
        public List<Interval> VNorth_Interval { get; set; } = new List<Interval>() { new Interval(0, 1) };

        public double CenterT0Value { get; set; }
        public double CenterT1Value { get; set; }

        public int BoundaryBaseLineCount { get; set; }
        public int CenterBaseLineCount { get; set; }

        public enum TrilinearShapeType
        {
            Unset = 0,
            EShape = 1,
            ZShape = 2,
            WangShape = 3,
            IIShape = 4,
            EightShape = 5
        }

        public TrilinearShapeType ShapeType { get; set; }

        List<List<double>> KeyPointTLoL { get; set; }

        List<List<bool>> KeyPointTypeLoL { get; set; }

        public TrilinearCell(Curve southLine, Curve middleLine, Curve northLine,
                             Curve westLine, Curve westLine1, Curve eastLine,
                             Curve eastLine1, int i, int j)
        {
            #region index
            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
            #endregion

            #region BaseLine
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
            if (middleLine != null)
            {
                this.MiddleBaseLine = new Line(middleLine.PointAtStart, middleLine.PointAtEnd);
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
            if (westLine1 != null)
            {
                this.WestBaseLine1 = new Line(westLine1.PointAtStart, westLine1.PointAtEnd);
                count++;
            }
            if (eastLine1 != null)
            {
                this.EastBaseLine1 = new Line(eastLine1.PointAtStart, eastLine1.PointAtEnd);
                count++;
            }
            #endregion

            this.ShapeType = TrilinearShapeType.EightShape;
            this.BoundaryBaseLineCount = count;
            this.CenterBaseLineCount = 0;
        }

        /// <summary>
        /// 生成EShape或ZShape
        /// </summary>
        private TrilinearCell(TrilinearCell initialCell, bool isEShape, int directionCode, int i, int j)
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
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;
            #endregion

            #region Lines
            this.South_Interval = new List<Interval>();
            this.South_Interval.Add(new Interval(0, 1));
            this.Middle_Interval = new List<Interval>();
            this.Middle_Interval.Add(new Interval(0, 1));
            this.North_Interval = new List<Interval>();
            this.North_Interval.Add(new Interval(0, 1));
            if (isEShape)
            {
                if (directionCode < 2)
                {
                    // 西侧
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(initialCell.SouthBaseLine.From, initialCell.SouthBaseLine.To));
                    //this.MiddleLines = new List<Line>();
                    //this.MiddleLines.Add(new Line(initialCell.MiddleBaseLine.From, initialCell.MiddleBaseLine.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(initialCell.NorthBaseLine.From, initialCell.NorthBaseLine.To));
                    //this.WestLines = new List<Line>();
                    //this.WestLines.Add(new Line(initialCell.WestBaseLine.From, initialCell.WestBaseLine.To));
                    //this.WestLines1 = new List<Line>();
                    //this.WestLines1.Add(new Line(initialCell.WestBaseLine1.From, initialCell.WestBaseLine1.To));

                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.West1_Interval = new List<Interval>();
                    this.West1_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    // 东侧
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(initialCell.SouthBaseLine.From, initialCell.SouthBaseLine.To));
                    //this.MiddleLines = new List<Line>();
                    //this.MiddleLines.Add(new Line(initialCell.MiddleBaseLine.From, initialCell.MiddleBaseLine.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(initialCell.NorthBaseLine.From, initialCell.NorthBaseLine.To));
                    //this.EastLines = new List<Line>();
                    //this.EastLines.Add(new Line(initialCell.EastBaseLine.From, initialCell.EastBaseLine.To));
                    //this.EastLines1 = new List<Line>();
                    //this.EastLines1.Add(new Line(initialCell.EastBaseLine1.From, initialCell.EastBaseLine1.To));

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
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(initialCell.SouthBaseLine.From, initialCell.SouthBaseLine.To));
                    //this.MiddleLines = new List<Line>();
                    //this.MiddleLines.Add(new Line(initialCell.MiddleBaseLine.From, initialCell.MiddleBaseLine.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(initialCell.NorthBaseLine.From, initialCell.NorthBaseLine.To));
                    //this.WestLines = new List<Line>();
                    //this.WestLines.Add(new Line(initialCell.WestBaseLine.From, initialCell.WestBaseLine.To));
                    //this.EastLines1 = new List<Line>();
                    //this.EastLines1.Add(new Line(initialCell.EastBaseLine1.From, initialCell.EastBaseLine1.To));
                    this.West_Interval = new List<Interval>();
                    this.West_Interval.Add(new Interval(0, 1));
                    this.East1_Interval = new List<Interval>();
                    this.East1_Interval.Add(new Interval(0, 1));
                }
                else
                {
                    // 东北
                    //this.SouthLines = new List<Line>();
                    //this.SouthLines.Add(new Line(initialCell.SouthBaseLine.From, initialCell.SouthBaseLine.To));
                    //this.MiddleLines = new List<Line>();
                    //this.MiddleLines.Add(new Line(initialCell.MiddleBaseLine.From, initialCell.MiddleBaseLine.To));
                    //this.NorthLines = new List<Line>();
                    //this.NorthLines.Add(new Line(initialCell.NorthBaseLine.From, initialCell.NorthBaseLine.To));
                    //this.WestLines1 = new List<Line>();
                    //this.WestLines1.Add(new Line(initialCell.WestBaseLine1.From, initialCell.WestBaseLine1.To));
                    //this.EastLines = new List<Line>();
                    //this.EastLines.Add(new Line(initialCell.EastBaseLine.From, initialCell.EastBaseLine.To));
                    this.West1_Interval = new List<Interval>();
                    this.West1_Interval.Add(new Interval(0, 1));
                    this.East_Interval = new List<Interval>();
                    this.East_Interval.Add(new Interval(0, 1));
                }
            }
            #endregion

            this.ShapeType = TrilinearShapeType.EShape;
            this.BoundaryBaseLineCount = 5;
            this.CenterBaseLineCount = 0;
        }

        /// <summary>
        /// 生成WangShape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="t0"></param>
        /// <param name="t1"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, int directionCode, double t0, double t1, int i, int j)
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
            this.MiddleBaseLine = initialCell.MiddleBaseLine;
            this.WestBaseLine1 = initialCell.WestBaseLine1;
            this.EastBaseLine1 = initialCell.EastBaseLine1;
            #endregion

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

            #region Interval
            this.South_Interval = new List<Interval>();
            this.South_Interval.Add(new Interval(0, 1));
            this.Middle_Interval = new List<Interval>();
            this.Middle_Interval.Add(new Interval(0, 1));
            this.North_Interval = new List<Interval>();
            this.North_Interval.Add(new Interval(0, 1));

            this.VSouth_Interval = new List<Interval>();
            this.VSouth_Interval.Add(new Interval(0, 1));
            this.VNorth_Interval = new List<Interval>();
            this.VNorth_Interval.Add(new Interval(0, 1));
            #endregion

            this.ShapeType = TrilinearShapeType.WangShape;
            this.BoundaryBaseLineCount = 3;
            this.CenterBaseLineCount = 2;
        }

        /// <summary>
        /// 生成IIshape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, int directionCode, double t,int i, int j)
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
            
            // 子类属性
            this.MiddleBaseLine = initialCell.MiddleBaseLine;

            /* t要小于0.3 */
            if (t > 0.3)
            {
                t = 0.3;
            }
            switch (directionCode)
            {
                case 0:
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.WestBaseLine1 = initialCell.WestBaseLine1;
                    this.EastBaseLine1 = initialCell.EastBaseLine1;
                    break;
                case 1:
                    Point3d vStart00 = this.SouthBaseLine.PointAt(t);
                    Point3d vEnd01 = this.MiddleBaseLine.PointAt(t);
                    Point3d vEnd02 = this.NorthBaseLine.PointAt(t);
                    this.WestBaseLine = new Line(vStart00, vEnd01);
                    this.WestBaseLine1 = new Line(vEnd01, vEnd02);
                    this.EastBaseLine = initialCell.EastBaseLine;
                    this.EastBaseLine1 = initialCell.EastBaseLine1;
                    break;
                case 2:
                    Point3d vStart10 = this.SouthBaseLine.PointAt(1 - t);
                    Point3d vEnd11 = this.MiddleBaseLine.PointAt(1 - t);
                    Point3d vEnd12 = this.NorthBaseLine.PointAt(1 - t);
                    this.WestBaseLine = initialCell.WestBaseLine;
                    this.WestBaseLine1 = initialCell.WestBaseLine1;
                    this.EastBaseLine = new Line(vStart10, vEnd11);
                    this.EastBaseLine = new Line(vEnd11, vEnd12);
                    break;
                case 3:
                    vStart00 = this.SouthBaseLine.PointAt(t);
                    vEnd01 = this.MiddleBaseLine.PointAt(t);
                    vEnd02 = this.NorthBaseLine.PointAt(t);
                    vStart10 = this.SouthBaseLine.PointAt(1 - t);
                    vEnd11 = this.MiddleBaseLine.PointAt(1 - t);
                    vEnd12 = this.NorthBaseLine.PointAt(1 - t);
                    this.WestBaseLine = new Line(vStart00, vEnd01);
                    this.WestBaseLine1 = new Line(vEnd01, vEnd02);
                    this.EastBaseLine = new Line(vStart10, vEnd11);
                    this.EastBaseLine = new Line(vEnd11, vEnd12);
                    break;
            }
            #endregion

            #region Lines
            this.South_Interval = new List<Interval>();
            this.South_Interval.Add(new Interval(0, 1));
            this.Middle_Interval = new List<Interval>();
            this.Middle_Interval.Add(new Interval(0, 1));
            this.North_Interval = new List<Interval>();
            this.North_Interval.Add(new Interval(0, 1));

            this.West_Interval = new List<Interval>();
            this.West_Interval.Add(new Interval(0, 1));
            this.East_Interval = new List<Interval>();
            this.East_Interval.Add(new Interval(0, 1));
            this.West1_Interval = new List<Interval>();
            this.West1_Interval.Add(new Interval(0, 1));
            this.East1_Interval = new List<Interval>();
            this.East1_Interval.Add(new Interval(0, 1));
            #endregion
        }

        public TrilinearCell GenerateEShape(int directionCode, int i, int j)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, true, directionCode, i, j);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateZShape(int directionCode, int i, int j)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, false, directionCode, i, j);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateWangShape(int directionCode, double t0, double t1, int i, int j)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, directionCode, t0, t1, i, j);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateIIShape(int directionCode, double t, int i, int j)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, directionCode, t, i, j);
            return newGeneratedCell;
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
            if (this.MiddleBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.Middle_Interval.Count; i++)
                {
                    Point3d p0 = this.MiddleBaseLine.PointAt(this.Middle_Interval[i].T0);
                    Point3d p1 = this.MiddleBaseLine.PointAt(this.Middle_Interval[i].T1);
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

            if (this.HSouthBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.HSouth_Interval.Count; i++)
                {
                    Point3d p0 = this.HSouthBaseLine.PointAt(this.HSouth_Interval[i].T0);
                    Point3d p1 = this.HSouthBaseLine.PointAt(this.HSouth_Interval[i].T1);
                    horizontalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.HNorthBaseLine != Line.Unset)
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
            if (this.WestBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.West_Interval.Count; i++)
                {
                    Point3d p0 = this.WestBaseLine.PointAt(this.West_Interval[i].T0);
                    Point3d p1 = this.WestBaseLine.PointAt(this.West_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.WestBaseLine1 != Line.Unset)
            {
                for (int i = 0; i < this.West1_Interval.Count; i++)
                {
                    Point3d p0 = this.WestBaseLine1.PointAt(this.West1_Interval[i].T0);
                    Point3d p1 = this.WestBaseLine1.PointAt(this.West1_Interval[i].T1);
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
            if (this.EastBaseLine1 != Line.Unset)
            {
                for (int i = 0; i < this.East1_Interval.Count; i++)
                {
                    Point3d p0 = this.EastBaseLine1.PointAt(this.East1_Interval[i].T0);
                    Point3d p1 = this.EastBaseLine1.PointAt(this.East1_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }

            if (this.VSouthBaseLine != Line.Unset)
            {
                for (int i = 0; i < this.VSouth_Interval.Count; i++)
                {
                    Point3d p0 = this.VSouthBaseLine.PointAt(this.VSouth_Interval[i].T0);
                    Point3d p1 = this.VSouthBaseLine.PointAt(this.VSouth_Interval[i].T1);
                    verticalLines.Add(new Line(p0, p1).ToNurbsCurve());
                }
            }
            if (this.VNorthBaseLine != Line.Unset)
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
    }
}
