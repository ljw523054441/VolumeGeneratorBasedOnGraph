using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    class FaceEdgeSegment : BoundarySegment
    {
        //public List<int> IncludedDVertice { get; set; }

        //public int HIndex { get; set; }

        //public FaceEdgeSegment(Line segment,string label, List<int> includedDVertice,int hIndex)
        //{
        //    #region 父类属性
        //    this.Line = segment;
        //    this.Label = label;

        //    this.From = segment.From;
        //    this.To = segment.To;
        //    #endregion

        //    this.IncludedDVertice = new List<int>();
        //    this.IncludedDVertice.AddRange(includedDVertice);

        //    this.HIndex = hIndex;
        //}

        //public FaceEdgeSegment(FaceEdgeSegment source)
        //{
        //    #region 父类属性
        //    this.Label = source.Label;
        //    this.From = new Point3d(source.From);
        //    this.To = new Point3d(source.To);
        //    this.Line = new Line(this.From, this.To);
        //    #endregion

        //    this.IncludedDVertice = new List<int>();
        //    this.IncludedDVertice.AddRange(source.IncludedDVertice);

        //    this.HIndex = source.HIndex;
        //}

        //public override void Reverse()
        //{
        //    Point3d tmp = this.From;
        //    this.From = this.To;
        //    this.To = tmp;

        //    this.Line = new Line(this.From, this.To);

        //    this.IncludedDVertice.Reverse();
        //}
    }
}
