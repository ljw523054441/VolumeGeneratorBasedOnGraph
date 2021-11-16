using Grasshopper;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
// using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class DualGraphWithHM : GraphWithHM
    {
        public List<List<int>> DVertexBelongsToWhichInnerNode { get; }

        public List<List<int>> DVertexBelongsToWhichVolume { get; }

        /// <summary>
        /// 构造空的DualGraphWithHM对象
        /// </summary>
        public DualGraphWithHM() : base()
        {
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="dualGraphWithHM"></param>
        public DualGraphWithHM(DualGraphWithHM source)
        {
            // 深拷贝List<List<int>>
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < source.GraphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(source.GraphTables[i]);
            }
            // 深拷贝List<Node>
            this.GraphNodes = new List<GraphNode>();
            for (int i = 0; i < source.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new GraphNode(source.GraphNodes[i]));
            }

            this.InnerNodeCount = source.InnerNodeCount;
            this.OuterNodeCount = source.OuterNodeCount;
            this.InnerNodeIndexList = source.InnerNodeIndexList;
            this.OuterNodeIndexList = source.OuterNodeIndexList;

            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(source.PlanktonMesh);
            this.TreeNodeLabel = (string)source.TreeNodeLabel.Clone();

            // 深拷贝List<List<int>>
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            for (int i = 0; i < source.DVertexBelongsToWhichInnerNode.Count; i++)
            {
                this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
                this.DVertexBelongsToWhichInnerNode[i].AddRange(source.DVertexBelongsToWhichInnerNode[i]);
            }
        }

        /// <summary>
        /// 用于构造DualGraphWithHM
        /// </summary>
        /// <param name="planktonMesh"></param>
        /// <param name="graphNodes"></param>
        /// <param name="graphTables"></param>
        /// <param name="dVertexBelongsToWhichInnerNode"></param>
        public DualGraphWithHM(PlanktonMesh planktonMesh,
                               List<GraphNode> graphNodes,
                               List<List<int>> graphTables,
                               string label,
                               List<List<int>> dVertexBelongsToWhichInnerNode) 
        {
            // 深拷贝List<List<int>>
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }
            // 深拷贝List<Node>
            this.GraphNodes = new List<GraphNode>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                this.GraphNodes.Add(new GraphNode(graphNodes[i]));
            }

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                if (graphNodes[i].IsInner)
                {
                    this.InnerNodeCount++;
                    this.InnerNodeIndexList.Add(i);
                }
                else
                {
                    this.OuterNodeCount++;
                    this.OuterNodeIndexList.Add(i);
                }
            }

            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(planktonMesh);
            this.TreeNodeLabel = label;

            // 深拷贝List<List<int>>
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            for (int i = 0; i < dVertexBelongsToWhichInnerNode.Count; i++)
            {
                this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
                this.DVertexBelongsToWhichInnerNode[i].AddRange(dVertexBelongsToWhichInnerNode[i]);
            }

            // 
        }
    }
}
