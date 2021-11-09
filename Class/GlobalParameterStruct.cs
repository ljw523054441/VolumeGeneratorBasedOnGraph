using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public struct GlobalParameter
    {
        public int BoundaryNodeCount;
        public int VolumeNodeCount;

        public Point3d[] VolumeNodePointLocations { get; set; }
        public Point3d[] BoundaryNodePointLocations { get; set; }

        public GlobalParameter(int volumeNodeCount, int boundaryNodeCount)
        {
            this.VolumeNodeCount = volumeNodeCount;
            this.BoundaryNodeCount = boundaryNodeCount;

            VolumeNodePointLocations = null;
            BoundaryNodePointLocations = null;
        }
    }
}