using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public  struct NodeAttribute
    {
        // 包括NEWS点的boundaryNodeCount
        public static int BoundaryNodeCount = 4;
        
        public string NodeLabel;
        public double NodeArea;

        public NodeAttribute(string nodeLabel, double nodeArea)
        {
            NodeLabel = nodeLabel;
            NodeArea = nodeArea;
        }
    }
}