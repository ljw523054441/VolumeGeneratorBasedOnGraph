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
        Line SouthBaseLine { get; set; }
        Line NorthBaseLine { get; set; }
        Line EastBaseLine { get; set; }
        Line WestBaseLine { get; set; }
        Line HCenterBaseLine { get; set; }
        Line VCenterBaseLine { get; set; }

        int[] CurrentIndex { get; set; }
        int[] PrevSouthIndex { get; set; }
        int[] PrevWestIndex { get; set; }
        int[] NexEastIndex { get; set; }
        int[] NexNorthIndex { get; set; }

        public Cell(Curve southLine, Curve northLine, Curve westLine, Curve eastLine, int i, int j)
        {
            if (southLine != null)
            {
                this.SouthBaseLine = new Line(southLine.PointAtStart, southLine.PointAtEnd);
            }
            if (northLine != null)
            {
                this.NorthBaseLine = new Line(northLine.PointAtStart, northLine.PointAtEnd);
            }
            if (westLine != null)
            {
                this.WestBaseLine = new Line(westLine.PointAtStart, westLine.PointAtEnd);
            }
            if (eastLine != null)
            {
                this.EastBaseLine = new Line(eastLine.PointAtStart, eastLine.PointAtEnd);
            }

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
        }

        public Cell(Curve centerH, Curve westLine, Curve eastLine, int i, int j)
        {
            if (centerH != null)
            {
                this.HCenterBaseLine = new Line(centerH.PointAtStart, centerH.PointAtEnd);
            }
            if (westLine != null)
            {
                this.WestBaseLine = new Line(westLine.PointAtStart, westLine.PointAtEnd);
            }
            if (eastLine != null)
            {
                this.EastBaseLine = new Line(eastLine.PointAtStart, eastLine.PointAtEnd);
            }

            this.CurrentIndex = new int[2] { i, j };
            this.PrevSouthIndex = new int[2] { i - 1, j };
            this.NexNorthIndex = new int[2] { i + 1, j };
            this.PrevWestIndex = new int[2] { i, j - 1 };
            this.NexEastIndex = new int[2] { i, j + 1 };
        }
    }
}
