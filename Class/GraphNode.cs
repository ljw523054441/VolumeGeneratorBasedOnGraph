using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class GraphNode : Node, ICloneable
    {
        // public Point3d NodeVertex { get; set; }

        public GraphNodeAttribute NodeAttribute { get; set; }

        public bool IsInner { get; }

        // public bool IsDecomposed { get; }


        public GraphNode()
        {

        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="node"></param>
        public GraphNode(GraphNode node)
        {
            this.NodeVertex = new Point3d(node.NodeVertex);
            this.NodeAttribute = new GraphNodeAttribute(node.NodeAttribute);

            this.IsInner = node.IsInner;
            // this.IsDecomposed = node.IsDecomposed;
        }

        public GraphNode(Point3d nodeVertex, GraphNodeAttribute nodeAttribute, bool isInner)
        {
            this.NodeVertex = nodeVertex;
            this.NodeAttribute = nodeAttribute;
            this.IsInner = isInner;
            // this.IsDecomposed = false;
        }

        /// <summary>
        /// 调用ICloneable接口实现深拷贝
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            GraphNode node = (GraphNode)MemberwiseClone();
            node.NodeAttribute = (GraphNodeAttribute)NodeAttribute.Clone();
            return node;
        }
    }
}
