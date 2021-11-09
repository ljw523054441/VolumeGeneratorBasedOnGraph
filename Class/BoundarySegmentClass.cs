using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class BoundarySegment : ICloneable
    {
        public Line Line { get; set; }

        public string Label { get; set; }

        public Point3d From { get; set; }
        public Point3d To { get; set; }

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
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="boundarySegment"></param>
        public BoundarySegment(BoundarySegment boundarySegment)
        {
            
            this.Label = boundarySegment.Label;
            this.From = new Point3d(boundarySegment.From);
            this.To = new Point3d(boundarySegment.To);
            this.Line = new Line(this.From, this.To);
        }

        public void Reverse()
        {
            Point3d tmp = this.From;
            this.From = this.To;
            this.To = tmp;

            this.Line = new Line(this.From, this.To);
        }

        /// <summary>
        /// 调用ICloneable接口实现深拷贝
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            BoundarySegment boundarySegment = (BoundarySegment)MemberwiseClone();
            boundarySegment.From = (Point3d)MemberwiseClone();
            boundarySegment.To = (Point3d)MemberwiseClone();
            boundarySegment.Line = new Line(boundarySegment.From, boundarySegment.To);
            return boundarySegment;
        }
    }
}