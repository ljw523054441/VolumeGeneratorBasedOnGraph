using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class BoundarySegment
    {
        public Line Line;

        public string Label;

        public Point3d From;
        public Point3d To;

        public BoundarySegment(Line segment, string label)
        {
            this.Line = segment;
            this.Label = label;

            this.From = segment.From;
            this.To = segment.To;
        }

        public void Reverse()
        {
            Point3d tmp = this.From;
            this.From = this.To;
            this.To = tmp;

            this.Line = new Line(this.From, this.To);
        }
    }
}