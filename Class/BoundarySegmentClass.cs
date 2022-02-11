using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class BoundarySegment
    {
        //public Line Line { get; set; }

        public PolylineCurve PolylineCurve { get; set; }

        public string Label { get; set; }

        public Point3d From { get; set; }
        public Point3d To { get; set; }

        public List<Point3d> Points { get; set; }

        public List<Line> Lines { get; set; }

        public List<double> TurningTs { get; set; }

        public List<int> PointOnWhichSegments { get; set; }

        public int[] HIndex { get; set; }
        public int FIndex { get; set; }
        public int AdjacentFIndex { get; set; }

        //public int AdjacentNodeIndex { get; set; }

        public bool IsZeroCurvatureAtStart { get; set; }

        public bool IsZeroCurvatureAtEnd { get; set; }

        public enum Location
        {
            Unset,
            N,
            W,
            S,
            E,
            Inner,
            adjacent
        }

        public Location LocationValue { get; set; } = Location.Unset;

        public double setback { get; set; }

        public BoundarySegment()
        {

        }

        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="segment"></param>
        /// <param name="label"></param>
        public BoundarySegment(Line segment, string label)
        {
            this.Points = new List<Point3d>();
            this.Points.Add(segment.From);
            this.Points.Add(segment.To);
            this.PolylineCurve = new PolylineCurve(this.Points);
            this.PolylineCurve.Domain = new Interval(0, 1);
            this.Label = label;

            this.From = Points[0];
            this.To = Points.Last();

            this.Lines = new List<Line>();
            this.Lines.Add(new Line(segment.From, segment.To));

            this.TurningTs = new List<double>();

            this.PointOnWhichSegments = new List<int>();

            this.HIndex = new int[1];
            this.FIndex = -1;
            this.AdjacentFIndex = -1;
        }

        public BoundarySegment(List<Point3d> points, int[] hIndex, int fIndex, int adjacentFIndex)
        {
            this.Points = new List<Point3d>();
            this.Points.AddRange(points);
            this.PolylineCurve = new PolylineCurve(this.Points);
            this.PolylineCurve.Domain = new Interval(0, 1);
            this.Label = null;

            this.From = this.Points[0];
            this.To = this.Points.Last();

            this.Lines = new List<Line>();
            Point3d start = this.Points[0];
            for (int i = 1; i < this.Points.Count; i++)
            {
                this.Lines.Add(new Line(start, this.Points[i]));
                start = this.Points[i];
            }

            this.TurningTs = new List<double>();
            for (int i = 0; i < this.Points.Count; i++)
            {
                double t;
                this.PolylineCurve.ClosestPoint(this.Points[i], out t);
                this.TurningTs.Add(t);
            }

            this.PointOnWhichSegments = new List<int>();

            this.HIndex = new int[hIndex.Length];
            for (int i = 0; i < hIndex.Length; i++)
            {
                this.HIndex[i] = hIndex[i];
            }
            this.FIndex = fIndex;
            this.AdjacentFIndex = adjacentFIndex;
        }

        public BoundarySegment(List<Line> segments, string label)
        {
            if (segments.Count == 0)
            {
                throw new Exception("必须有至少一段线来构成BoundarySegment");
            }
            if (segments.Count == 1)
            {
                this.Points = new List<Point3d>();
                this.Points.Add(segments[0].From);
                this.Points.Add(segments[0].To);
                this.PolylineCurve = new PolylineCurve(this.Points);
                this.PolylineCurve.Domain = new Interval(0, 1);
                this.Label = label;

                this.From = Points[0];
                this.To = Points.Last();

                this.Lines = new List<Line>();
                this.Lines.Add(new Line(segments[0].From, segments[0].To));

                this.TurningTs = new List<double>();
                this.PointOnWhichSegments = new List<int>();

                this.HIndex = new int[1];
                this.FIndex = -1;
                this.AdjacentFIndex = -1;
            }
            
            this.Points = new List<Point3d>();
            this.Points.Add(segments[0].From);
            for (int i = 0; i < segments.Count; i++)
            {
                this.Points.Add(segments[i].To);
            }
            this.PolylineCurve = new PolylineCurve(this.Points);
            this.PolylineCurve.Domain = new Interval(0, 1);
            this.Label = label;

            this.From = Points[0];
            this.To = Points.Last();

            this.Lines = new List<Line>();
            for (int i = 0; i < segments.Count; i++)
            {
                this.Lines.Add(new Line(segments[i].From, segments[i].To));
            }

            this.TurningTs = new List<double>();
            for (int i = 0; i < this.Points.Count; i++)
            {
                double t;
                this.PolylineCurve.ClosestPoint(this.Points[i], out t);
                this.TurningTs.Add(t);
            }

            this.PointOnWhichSegments = new List<int>();

            this.HIndex = new int[segments.Count];
            this.FIndex = -1;
            this.AdjacentFIndex = -1;
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="source"></param>
        public BoundarySegment(BoundarySegment source)
        {
            this.Points = new List<Point3d>();
            for (int i = 0; i < source.Points.Count; i++)
            {
                this.Points.Add(new Point3d(source.Points[i]));
            }
            this.PolylineCurve = new PolylineCurve(source.PolylineCurve);
            this.PolylineCurve.Domain = new Interval(0, 1);
            this.Label = source.Label;

            this.From = new Point3d(source.From);
            this.To = new Point3d(source.To);

            this.Lines = new List<Line>();
            for (int i = 0; i < source.Lines.Count; i++)
            {
                this.Lines.Add(new Line(source.Lines[i].From, source.Lines[i].To));
            }

            
            if (source.TurningTs == null)
            {
                this.TurningTs = null;
            }
            else
            {
                this.TurningTs = new List<double>();
                this.TurningTs.AddRange(source.TurningTs);
            }

            this.PointOnWhichSegments = new List<int>();
            this.PointOnWhichSegments.AddRange(source.PointOnWhichSegments);

            this.HIndex = new int[source.HIndex.Length];
            for (int i = 0; i < source.HIndex.Length; i++)
            {
                this.HIndex[i] = source.HIndex[i];
            }
            this.FIndex = source.FIndex;
            this.AdjacentFIndex = source.AdjacentFIndex;

            this.IsZeroCurvatureAtStart = source.IsZeroCurvatureAtStart;
            this.IsZeroCurvatureAtEnd = source.IsZeroCurvatureAtEnd;
            this.LocationValue = source.LocationValue;
            this.setback = source.setback;
        }

        public virtual void Reverse()
        {
            this.Points.Reverse();
            this.PolylineCurve.Reverse();

            Point3d tmp = this.From;
            this.From = this.To;
            this.To = tmp;

            List<Line> newLines = new List<Line>();
            for (int i = 0; i < this.Lines.Count; i++)
            {
                Line newLine = new Line(this.Lines[i].To, this.Lines[i].From);
                newLines.Add(newLine);
            }
            newLines.Reverse();
            this.Lines = new List<Line>();
            this.Lines.AddRange(newLines);
        }

        /// <summary>
        /// 重写ToString()方法
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            if (this.Label == null)
            {
                return "BS:" + this.LocationValue.ToString();
            }
            
            return "BS:" + this.Label;
        }

        ///// <summary>
        ///// 调用ICloneable接口实现深拷贝
        ///// </summary>
        ///// <returns></returns>
        //public object Clone()
        //{
        //    BoundarySegment boundarySegment = (BoundarySegment)MemberwiseClone();
        //    boundarySegment.From = (Point3d)MemberwiseClone();
        //    boundarySegment.To = (Point3d)MemberwiseClone();
        //    boundarySegment.Line = new Line(boundarySegment.From, boundarySegment.To);
        //    return boundarySegment;
        //}
    }

    //public class BoundarySegmentGoo : GH_GeometricGoo<BoundarySegment>
    //{
    //    public 
    //}
}