using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public  class NodeAttribute : ICloneable
    {
        /// <summary>
        /// 体量的名称标签
        /// </summary>
        public string NodeLabel { get; set; }

        /// <summary>
        /// 体量的输入面积，需要除以全部体量的面积和，再乘以场地面积，来得到占地比例
        /// </summary>
        public double NodeArea { get; set; }

        /// <summary>
        /// 每个体量占场地的面积比例，[0,1]区间
        /// </summary>
        public double NodeAreaProportion { get; set; }

        /// <summary>
        /// 节点的Connectivity表
        /// </summary>
        public int[] ConnectivityTable { get; set; }
        /// <summary>
        /// 节点的Adjacency表
        /// </summary>
        public int[] AdjacencyTable { get; set; }

        /// <summary>
        /// 开间数
        /// </summary>
        public int SpanNumX { get; set; }
        /// <summary>
        /// 进深数
        /// </summary>
        public int SpanNumY { get; set; }

        /// <summary>
        /// 开间方向的柱跨尺寸
        /// </summary>
        public double SpanSizeX { get; set; }
        /// <summary>
        /// 进深方向的柱跨尺寸
        /// </summary>
        public double SpanSizeY { get; set; }

        /// <summary>
        /// 层数
        /// </summary>
        public int FloorNumZ { get; set; }
        /// <summary>
        /// 层高
        /// </summary>
        public double FloorHeightZ { get; set; }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="nodeAttribute"></param>
        public NodeAttribute(NodeAttribute nodeAttribute)
        {
            this.NodeLabel = nodeAttribute.NodeLabel;
            this.NodeArea = nodeAttribute.NodeArea;
            this.NodeAreaProportion = nodeAttribute.NodeAreaProportion;

            this.ConnectivityTable = (int[])nodeAttribute.ConnectivityTable.Clone();
            this.AdjacencyTable = (int[])nodeAttribute.AdjacencyTable.Clone();

            this.SpanNumX = nodeAttribute.SpanNumX;
            this.SpanNumY = nodeAttribute.SpanNumY;
            this.SpanSizeX = nodeAttribute.SpanSizeX;
            this.SpanSizeY = nodeAttribute.SpanSizeY;
            this.FloorNumZ = nodeAttribute.FloorNumZ;
            this.FloorHeightZ = nodeAttribute.FloorHeightZ;
        }

        /// <summary>
        /// 构造边界node
        /// </summary>
        /// <param name="nodeLabel"></param>
        public NodeAttribute(string nodeLabel)
        {
            this.NodeLabel = nodeLabel;

            this.NodeArea = 0;
            this.NodeAreaProportion = 0;
            //this.ConnectivityTable = new List<int>();
            //this.AdjacencyTable = new List<int>();
        }
        
        public NodeAttribute(string nodeLabel, double nodeArea)
        {
            this.NodeLabel = nodeLabel;
            this.NodeArea = nodeArea;

            this.NodeAreaProportion = 0;
        }

        /// <summary>
        /// 构造体量node
        /// </summary>
        /// <param name="nodeLabel"></param>
        /// <param name="nodeArea"></param>
        /// <param name="spanNumX"></param>
        /// <param name="spanNumY"></param>
        /// <param name="floorNumZ"></param>
        /// <param name="spanSizeX"></param>
        /// <param name="spanSizeY"></param>
        /// <param name="floorHeightZ"></param>
        public NodeAttribute(string nodeLabel, 
                             double nodeArea, 
                             int spanNumX, 
                             int spanNumY, 
                             int floorNumZ, 
                             double spanSizeX, 
                             double spanSizeY,
                             double floorHeightZ)
        {
            this.NodeLabel = nodeLabel;
            this.NodeArea = nodeArea;

            this.NodeAreaProportion = 0;

            //this.ConnectivityTable = n
            //this.AdjacencyTable = new List<int>();

            this.SpanNumX = spanNumX;
            this.SpanNumY = spanNumY;
            this.FloorNumZ = floorNumZ;
            this.SpanSizeX = spanSizeX;
            this.SpanSizeY = spanSizeY;
            this.FloorHeightZ = floorHeightZ;
        }

        /// <summary>
        /// 调用ICloneable接口实现深拷贝
        /// </summary>
        /// <returns></returns>
        public object Clone()
        {
            NodeAttribute nodeAttribute = (NodeAttribute)MemberwiseClone();
            nodeAttribute.ConnectivityTable = (int[])ConnectivityTable.Clone();
            nodeAttribute.AdjacencyTable = (int[])AdjacencyTable.Clone();
            return nodeAttribute;
        }


    }
}