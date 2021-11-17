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
    public class DualGraphWithHM : Graph
    {
        public PlanktonMesh PlanktonMesh { get; set; }

        public List<int> DualNodeCorrespondingInnerNodeIndex { get; }

        public List<List<int>> FaceIndexsAroundOuterNodes { get; }

        public List<List<int>> DVertexBelongsToWhichInnerNode { get; }

        public Dictionary<int,List<int>> DVertexBelongsToWhichVolume { get; }



        /// <summary>
        /// 构造空的DualGraphWithHM对象
        /// </summary>
        public DualGraphWithHM() : base()
        {
            this.PlanktonMesh = new PlanktonMesh();
            this.DualNodeCorrespondingInnerNodeIndex = new List<int>();
            this.FaceIndexsAroundOuterNodes = new List<List<int>>();
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
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

            // 深拷贝DualNodeCorrespondingInnerNodeIndex
            this.DualNodeCorrespondingInnerNodeIndex = new List<int>();
            this.DualNodeCorrespondingInnerNodeIndex.AddRange(source.DualNodeCorrespondingInnerNodeIndex);

            // 深拷贝FaceIndexsAroundOuterNodes
            this.FaceIndexsAroundOuterNodes = new List<List<int>>();
            for (int i = 0; i < source.FaceIndexsAroundOuterNodes.Count; i++)
            {
                this.FaceIndexsAroundOuterNodes.Add(new List<int>());
                this.FaceIndexsAroundOuterNodes[i].AddRange(source.FaceIndexsAroundOuterNodes[i]);
            }

            // 深拷贝DVertexBelongsToWhichInnerNode
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            for (int i = 0; i < source.DVertexBelongsToWhichInnerNode.Count; i++)
            {
                this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
                this.DVertexBelongsToWhichInnerNode[i].AddRange(source.DVertexBelongsToWhichInnerNode[i]);
            }

            // 深拷贝DVertexBelongsToWhichVolume
            this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
            foreach (KeyValuePair<int,List<int>> pair in source.DVertexBelongsToWhichVolume)
            {
                this.DVertexBelongsToWhichVolume.Add(pair.Key, pair.Value);
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
                               List<int> dualNodeCorrespondingInnerNodeIndex,
                               List<List<int>> faceIndexsAroundOuterNodes,
                               List<List<int>> dVertexBelongsToWhichInnerNode,
                               Dictionary<int, List<int>> dVertexBelongsToWhichVolume) 
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

            // 深拷贝DualNodeCorrespondingInnerNodeIndex
            this.DualNodeCorrespondingInnerNodeIndex = new List<int>();
            this.DualNodeCorrespondingInnerNodeIndex.AddRange(dualNodeCorrespondingInnerNodeIndex);

            // 深拷贝FaceIndexsAroundOuterNodes
            this.FaceIndexsAroundOuterNodes = new List<List<int>>();
            for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
            {
                this.FaceIndexsAroundOuterNodes.Add(new List<int>());
                this.FaceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNodes[i]);
            }

            // 深拷贝DVertexBelongsToWhichInnerNode
            this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            for (int i = 0; i < dVertexBelongsToWhichInnerNode.Count; i++)
            {
                this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
                this.DVertexBelongsToWhichInnerNode[i].AddRange(dVertexBelongsToWhichInnerNode[i]);
            }

            // 深拷贝DVertexBelongsToWhichVolume
            this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
            foreach (KeyValuePair<int, List<int>> pair in dVertexBelongsToWhichVolume)
            {
                this.DVertexBelongsToWhichVolume.Add(pair.Key, pair.Value);
            }
        }
    }
}
