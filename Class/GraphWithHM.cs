using Grasshopper.Kernel;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
//using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class GraphWithHM : Graph
    {
        /// <summary>
        /// 对应的半边网格
        /// </summary>
        public PlanktonMesh PlanktonMesh { get; set; }

        public string TreeNodeLabel { get; set; }

        /// <summary>
        /// 构造空的GraphWithHFMesh对象
        /// </summary>
        public GraphWithHM():base()
        {
            this.PlanktonMesh = new PlanktonMesh();
            this.TreeNodeLabel = null;
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="source"></param>
        public GraphWithHM(GraphWithHM source)
        {
            // 深拷贝List<List<int>>
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < source.GraphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(source.GraphTables[i]);
            }
            // 深拷贝List<Node>
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < source.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(source.GraphNodes[i]));
            }

            this.InnerNodeCount = source.InnerNodeCount;
            this.OuterNodeCount = source.OuterNodeCount;
            this.InnerNodeIndexList = source.InnerNodeIndexList;
            this.OuterNodeIndexList = source.OuterNodeIndexList;


            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(source.PlanktonMesh);

            this.TreeNodeLabel = (string)source.TreeNodeLabel.Clone();
        }

        /// <summary>
        /// 用于在还没有进行Decompose时，构造GraphWithHFMesh对象
        /// </summary>
        public GraphWithHM(PlanktonMesh planktonMesh,
                           List<Node> graphNodes,
                           List<List<int>> graphTables)
        {
            // 深拷贝List<List<int>>
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }
            // 深拷贝List<Node>
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(graphNodes[i]));
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

            this.TreeNodeLabel = "尚未进行过Decompose的GraphWithHFMesh对象";
        }

        /// <summary>
        /// 用于在进行Decompose后，构造GraphWithHFMesh对象
        /// </summary>
        /// <param name="planktonMesh"></param>
        /// <param name="graphTables"></param>
        /// <param name="graphNodes"></param>
        /// <param name="decomposeTheithVertex"></param>
        /// <param name="decomposeTheithPairHFVertexIndex"></param>
        /// <param name="decomposeTheithPairHFResult"></param>
        public GraphWithHM(PlanktonMesh planktonMesh,
                           List<Node> graphNodes,
                           List<List<int>> graphTables, 
                           string treeNodeLabel)
        {
            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(planktonMesh);
            // 深拷贝List<List<int>>
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }
            // 深拷贝List<Node>
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(graphNodes[i]));
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

            this.TreeNodeLabel = treeNodeLabel;
        }

        



    }
}