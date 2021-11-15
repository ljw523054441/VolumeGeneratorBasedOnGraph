using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class Node
    {
        public Point3d NodeVertex { get; set; }

        public Node()
        {

        }

        public Node(Node node)
        {
            this.NodeVertex = new Point3d(node.NodeVertex);
        }

        public Node(Point3d nodeVertex)
        {
            this.NodeVertex = nodeVertex;
        }
    }
}
