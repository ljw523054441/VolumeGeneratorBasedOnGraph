using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class Graph
    {
        public List<Node> GraphNodes 
        {
            get
            {
                return this.GraphNodes;
            }
            set
            {
                this.GraphNodes = new List<Node>();
                for (int i = 0; i < value.Count; i++)
                {
                    this.GraphNodes.Add(new Node(value[i]));
                }

                this.InnerNodeCount = 0;
                this.OuterNodeCount = 0;
                this.InnerNodeIndexList = new List<int>();
                this.OuterNodeIndexList = new List<int>();
                for (int i = 0; i < this.GraphNodes.Count; i++)
                {
                    if (this.GraphNodes[i].IsInner)
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
            }
        }
        public List<List<int>> GraphTables 
        {
            get
            {
                return this.GraphTables;
            }
            set
            {
                if (value.Count != this.GraphNodes.Count)
                {
                    throw new Exception("未能设置GraphTables，因为GraphTables.Count与当前Graph中的Node数量不同");
                }
                else
                {
                    this.GraphTables = new List<List<int>>();
                    for (int i = 0; i < value.Count; i++)
                    {
                        this.GraphTables.Add(new List<int>());
                        this.GraphTables[i].AddRange(value[i]);
                    }
                }
            }
        }

        public int InnerNodeCount { get; internal set; }
        public int OuterNodeCount { get; internal set; }

        public List<int> InnerNodeIndexList { get; internal set; }
        public List<int> OuterNodeIndexList { get; internal set; }

        public Graph()
        {
            this.GraphNodes = new List<Node>();
            this.GraphTables = new List<List<int>>();
        }

        public Graph(List<Node> graphNodes, List<List<int>> graphTables)
        {
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(graphNodes[i]));
            }

            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            for (int i = 0; i < this.GraphNodes.Count; i++)
            {
                if (this.GraphNodes[i].IsInner)
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
        }

        public Graph(Graph graph)
        {
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < graph.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(graph.GraphNodes[i]));
            }

            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graph.GraphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graph.GraphTables[i]);
            }

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            for (int i = 0; i < this.GraphNodes.Count; i++)
            {
                if (this.GraphNodes[i].IsInner)
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
        }

    }
}