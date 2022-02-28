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

        public enum TrilinearShapeType
        {
            Unset = 0,
            EShape = 1,
            EVarientShape = 2,
            ZShape = 3,
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

            this.ShapeType = TrilinearShapeType.EShape;
            this.BoundaryBaseLineCount = 5;
            this.CenterBaseLineCount = 0;
        }

        // todo:要新增对于面宽方向上分解过的体量线 interval 值的规避
        /// <summary>
        /// 生成EVariantShape
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
        }

        // todo:要新增对于面宽方向上分解过的体量线 interval 值的规避
        /// <summary>
        /// 生成IIshape
        /// </summary>
        /// <param name="initialCell"></param>
        /// <param name="directionCode"></param>
        /// <param name="t"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private TrilinearCell(TrilinearCell initialCell, double t, int i, int j)
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

            /* t要小于0.3 */
            //if (t > 0.1)
            //{
            //    t = t * 0.1;
            //}
            double remapT = Cell.MinT + (Cell.MaxT - Cell.MinT) * t;
            Point3d vStart00 = this.SouthBaseLine.PointAt(remapT);
            Point3d vEnd01 = this.MiddleBaseLine.PointAt(remapT);
            Point3d vEnd02 = this.NorthBaseLine.PointAt(remapT);
            Point3d vStart10 = this.SouthBaseLine.PointAt(1 - remapT);
            Point3d vEnd11 = this.MiddleBaseLine.PointAt(1 - remapT);
            Point3d vEnd12 = this.NorthBaseLine.PointAt(1 - remapT);
            this.WestBaseLine = new Line(vStart00, vEnd01);
            this.WestBaseLine1 = new Line(vEnd01, vEnd02);
            this.EastBaseLine = new Line(vStart10, vEnd11);
            this.EastBaseLine1 = new Line(vEnd11, vEnd12);

            #endregion

            this.ShapeType = TrilinearShapeType.IIShape;
            this.BoundaryBaseLineCount = 7;
            this.CenterBaseLineCount = 0;
        }

        public TrilinearCell GenerateEShape(int directionCode)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, true, directionCode, this.CurrentIndex[0], this.CurrentIndex[1]);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateEVariantShape(int directionCode, int randomT0, int randomT1)
        {
            TrilinearCell newGeneratedCell = null;
            newGeneratedCell = new TrilinearCell(this, directionCode, 0.2 + randomT0 * 0.2, 0.2 + randomT1 * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);

            // todo：基于水平分割进行Vertical位置的修改
            //if (this.South_Interval.Count != 0 || this.Middle_Interval.Count != 0 || this.North_Interval.Count != 0)
            //{
                  // todo
            //}
            //else
            //{
            //    newGeneratedCell = new TrilinearCell(this, directionCode, 0.2 + randomT0 * 0.2, 0.2 + randomT1 * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
            //}


            return newGeneratedCell;
        }

        public TrilinearCell GenerateZShape(int directionCode)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, false, directionCode, this.CurrentIndex[0], this.CurrentIndex[1]);
            return newGeneratedCell;
        }

        public TrilinearCell GenerateIIShape(int directionCode, int randomT)
        {
            TrilinearCell newGeneratedCell = new TrilinearCell(this, 0.2 + randomT * 0.2, this.CurrentIndex[0], this.CurrentIndex[1]);
            return newGeneratedCell;
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
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.Middle_Interval = new List<Interval>();
                        this.Middle_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.Middle_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.Middle_Interval.Add(new Interval(factors[i], factors[i + 1]));
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
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.HSouth_Interval = new List<Interval>();
                        this.HSouth_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.HSouth_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.HSouth_Interval.Add(new Interval(factors[i], factors[i + 1]));
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
                        Point3d newEnd = line.PointAtLength(lMax);
                        double tEnd = line.ClosestParameter(newEnd);
                        this.HNouth_Interval = new List<Interval>();
                        this.HNouth_Interval.Add(new Interval(0, tEnd));
                    }
                    else if (lMinCount == 2 && lMaxCount == 1)
                    {
                        List<double> factors = GetFactors(false, lMin, lMax, dw, lMinCount);

                        this.HNouth_Interval = new List<Interval>();
                        for (int i = 0; i < lMinCount; i++)
                        {
                            this.HNouth_Interval.Add(new Interval(factors[i], factors[i + 1]));
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
    }
}
