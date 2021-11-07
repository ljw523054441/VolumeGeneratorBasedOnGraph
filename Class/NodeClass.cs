using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class Node : ICloneable
    {
        public Point3d NodeVertex { get; set; }

        public NodeAttribute NodeAttribute { get; set; }

        public bool IsInner { get; }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="node"></param>
        public Node(Node node)
        {
            this.NodeVertex = new Point3d(node.NodeVertex);
            this.NodeAttribute = new NodeAttribute(node.NodeAttribute);

            this.IsInner = node.IsInner;
        }

        public Node(Point3d nodeVertex, NodeAttribute nodeAttribute, bool isInner)
        {
            this.NodeVertex = nodeVertex;
            this.NodeAttribute = nodeAttribute;
            this.IsInner = isInner;
        }

        /// <summary>
        /// 调用ICloneable接口实现深拷贝
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            Node node = (Node)MemberwiseClone();
            node.NodeAttribute = (NodeAttribute)NodeAttribute.Clone();
            return node;
        }
    }
}
