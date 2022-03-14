using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class Cell
    {
        public static double Tolerance { get; set; }

        public Line SouthBaseLine { get; set; } = Line.Unset;
        public Line NorthBaseLine { get; set; } = Line.Unset;
        public Line EastBaseLine { get; set; } = Line.Unset;
        public Line WestBaseLine { get; set; } = Line.Unset;

        //public List<Line> SouthLines { get; set; }
        //public List<Line> NorthLines { get; set; }
        //public List<Line> EastLines { get; set; }
        //public List<Line> WestLines { get; set; }

        public List<Interval> South_Interval { get; set; }
        public List<Interval> North_Interval { get; set; }
        public List<Interval> East_Interval { get; set; }
        public List<Interval> West_Interval { get; set; }

        //public int[] CurrentIndex { get; set; }
        //public int[] PrevSouthIndex { get; set; }
        //public int[] PrevWestIndex { get; set; }
        //public int[] NexEastIndex { get; set; }
        //public int[] NexNorthIndex { get; set; }

        // 注意 singlePolyline 在使用 PublicBaseLineIndex 需要 +2，才是所希望得到的数字
        public int PublicBaseLineIndex { get; set; } = -1;
        public int AnotherBaseLineIndexRelatedToCutPoint { get; set; } = -1;

        public int PrevBaseLineIndexRelatedToCutPoint { get; set; } = -1;
        public int NextBaseLineIndexRelatedToCutPoint { get; set; } = -1;

        public Point3d TurningPoint { get; set; } = Point3d.Unset;

        public Curve CellBoundary { get; set; }

        public Brep[] CellBreps { get; set; }

        /// <summary>
        /// 构造空的Cell对象
        /// </summary>
        public Cell()
        {

        }

        public bool IsBilinearCell()
        {
            if (this is BilinearCell)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public virtual List<Curve> GetAllHorizontal()
        {
            if (this.IsBilinearCell())
            {
                BilinearCell b = this as BilinearCell;
                return b.GetAllHorizontal();
            }
            else
            {
                TrilinearCell t = this as TrilinearCell;
                return t.GetAllHorizontal();
            }
        }

        public virtual List<Curve> GetAllVertical()
        {
            if (this.IsBilinearCell())
            {
                BilinearCell b = this as BilinearCell;
                return b.GetAllVertical();
            }
            else
            {
                TrilinearCell t = this as TrilinearCell;
                return t.GetAllVertical();
            }
        }

        public virtual double GetAllHCurveLength(out int hCount)
        {
            double hLength = 0;
            List<Curve> hCurves = this.GetAllHorizontal();
            hCount = hCurves.Count;
            for (int i = 0; i < hCurves.Count; i++)
            {
                hLength += hCurves[i].GetLength();
            }
            return hLength;
        }

        public virtual double GetAllVCurveLength(out int vCount)
        {
            double vLength = 0;
            List<Curve> vCurves = this.GetAllVertical();
            vCount = vCurves.Count;
            for (int i = 0; i < vCurves.Count; i++)
            {
                vLength += vCurves[i].GetLength();
            }
            return vLength;
        }
    }
}
