using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public  class NodeAttribute
    {
        // 包括NEWS点的boundaryNodeCount
        //public static int BoundaryNodeCount = 4;

        //public static int VolumeNodeCount;
        
        /// <summary>
        /// 体量的名称标签
        /// </summary>
        public string NodeLabel;

        /// <summary>
        /// 体量的输入面积，需要除以全部体量的面积和，再乘以场地面积，来得到占地比例
        /// </summary>
        public double NodeArea;

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

            //this.ConnectivityTable = n
            //this.AdjacencyTable = new List<int>();
        }

        //public void SetNodeAreaProportion(double x)
        //{
        //    this.NodeAreaProportion = x;
        //}


    }
}