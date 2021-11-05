using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class Node
    {
        public Point3d NodeVertex { get; set; }

        public NodeAttribute NodeAttribute { get; set; }

        public bool IsInner { get; }

        public Node(Point3d nodeVertex, NodeAttribute nodeAttribute, bool isInner)
        {
            this.NodeVertex = nodeVertex;
            this.NodeAttribute = nodeAttribute;
            this.IsInner = isInner;
        }
    }
}
