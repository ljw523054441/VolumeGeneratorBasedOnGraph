using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class BoundarySegment
    {
        public Line Line { get; set; }

        public string Label { get; set; }

        public Point3d From { get; set; }
        public Point3d To { get; set; }

        public List<int> IncludedDVertice { get; set; }

        public int HIndex { get; set; }

        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="segment"></param>
        /// <param name="label"></param>
        public BoundarySegment(Line segment, string label)
        {
            this.Line = segment;
            this.Label = label;

            this.From = segment.From;
            this.To = segment.To;

            this.IncludedDVertice = null;
            this.HIndex = -1;
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="source"></param>
        public BoundarySegment(BoundarySegment source)
        {
            
            this.Label = source.Label;
            this.From = new Point3d(source.From);
            this.To = new Point3d(source.To);
            this.Line = new Line(this.From, this.To);

            if (source.IncludedDVertice == null)
            {
                this.IncludedDVertice = null;
            }
            else
            {
                this.IncludedDVertice = new List<int>();
                this.IncludedDVertice.AddRange(source.IncludedDVertice);
            }

            this.HIndex = source.HIndex;
        }

        public void Reverse()
        {
            Point3d tmp = this.From;
            this.From = this.To;
            this.To = tmp;

            this.Line = new Line(this.From, this.To);

            this.IncludedDVertice.Reverse();
        }

        /// <summary>
        /// 重写ToString()方法
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
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