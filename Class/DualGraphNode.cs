using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    class DualGraphNode : Node, ICloneable
    {
        public DualGraphNodeAttribute NodeAttribute { get; set; }

        public DualGraphNode()
        {

        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="node"></param>
        public DualGraphNode(DualGraphNode node)
        {
            this.NodeVertex = new Point3d(node.NodeVertex);
            this.NodeAttribute = new DualGraphNodeAttribute(node.NodeAttribute);
        }

        public DualGraphNode(Point3d nodeVertex, DualGraphNodeAttribute nodeAttribute)
        {
            this.NodeVertex = nodeVertex;
            this.NodeAttribute = nodeAttribute;
        }

        public object Clone()
        {
            DualGraphNode node = (DualGraphNode)MemberwiseClone();
            node.NodeAttribute = (DualGraphNodeAttribute)NodeAttribute.Clone();
            return node;
        }
    }
}
