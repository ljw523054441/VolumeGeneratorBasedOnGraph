using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    class Node
    {
        public Point3d NodeVertex { get; set; }

        public NodeAttribute NodeAttribute { get; set; }

        public Node(Point3d nodeVertex, NodeAttribute nodeAttribute)
        {
            this.NodeVertex = nodeVertex;
            this.NodeAttribute = nodeAttribute;
        }
    }
}
