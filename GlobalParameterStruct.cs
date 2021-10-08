using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public struct GlobalParameter
    {
        public int BoundaryNodeCount;
        public int VolumeNodeCount;

        public GlobalParameter(int volumeNodeCount,int boundaryNodeCount)
        {
            VolumeNodeCount = volumeNodeCount;
            BoundaryNodeCount = boundaryNodeCount;
        }
    }
}