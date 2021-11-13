using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class Graph
    {
        private List<Node> _graphNodes;
        private List<List<int>> _graphTables;
        private int _innerNodeCount;
        private int _outerNodeCount;
        private List<int> _innerNodeIndexList;
        private List<int> _outerNodeIndexList;
        
        public List<Node> GraphNodes
        {
            get
            {
                return _graphNodes;
            }
            set
            {
                _graphNodes = new List<Node>();
                for (int i = 0; i < value.Count; i++)
                {
                    _graphNodes.Add(new Node(value[i]));
                }

                _innerNodeCount = 0;
                _outerNodeCount = 0;
                _innerNodeIndexList = new List<int>();
                _outerNodeIndexList = new List<int>();
                for (int i = 0; i < _graphNodes.Count; i++)
                {
                    if (_graphNodes[i].IsInner)
                    {
                        _innerNodeCount++;
                        _innerNodeIndexList.Add(i);
                    }
                    else
                    {
                        _outerNodeCount++;
                        _outerNodeIndexList.Add(i);
                    }
                }

            }
        }
        public List<List<int>> GraphTables
        {
            get
            {
                return _graphTables;
            }
            set
            {
                if (value.Count != _graphNodes.Count && _graphNodes.Count != 0)
                {
                    throw new Exception("未能设置GraphTables，因为GraphTables.Count与当前Graph中的Node数量不同");
                }
                if (value.Count != _graphNodes.Count && _graphNodes.Count == 0)
                {
                    _graphTables = new List<List<int>>();
                    for (int i = 0; i < value.Count; i++)
                    {
                        _graphTables.Add(new List<int>());
                        _graphTables[i].AddRange(value[i]);
                    }
                }
                else
                {
                    _graphTables = new List<List<int>>();
                    for (int i = 0; i < value.Count; i++)
                    {
                        _graphTables.Add(new List<int>());
                        _graphTables[i].AddRange(value[i]);
                    }
                }
                
            }
        }

        public int InnerNodeCount { get; set; }
        public int OuterNodeCount { get; set; }

        public List<int> InnerNodeIndexList { get; set; }
        public List<int> OuterNodeIndexList { get; set; }

        /// <summary>
        /// 构造空的Graph对象
        /// </summary>
        public Graph()
        {
            this.GraphNodes = new List<Node>();
            this.GraphTables = new List<List<int>>();
            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
        }

        /// <summary>
        /// 用拷贝构造函数来实现深拷贝
        /// </summary>
        /// <param name="source"></param>
        public Graph(Graph source)
        {
            List<Node> nodes = new List<Node>();
            for (int i = 0; i < source.GraphNodes.Count; i++)
            {
                nodes.Add(new Node(source.GraphNodes[i]));
            }
            this.GraphNodes = nodes;

            List<List<int>> tables = new List<List<int>>();
            for (int i = 0; i < source.GraphTables.Count; i++)
            {
                tables.Add(new List<int>());
                tables[i].AddRange(source.GraphTables[i]);
            }
            this.GraphTables = tables;

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            for (int i = 0; i < nodes.Count; i++)
            {
                if (nodes[i].IsInner)
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

        /// <summary>
        /// 用GraphNodes和GraphTables来组装新的Graph对象
        /// </summary>
        /// <param name="graphNodes"></param>
        /// <param name="graphTables"></param>
        public Graph(List<Node> graphNodes, List<List<int>> graphTables)
        {
            List<Node> nodes = new List<Node>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                nodes.Add(new Node(graphNodes[i]));
            }
            this.GraphNodes = nodes;

            List<List<int>> tables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                tables.Add(new List<int>());
                tables[i].AddRange(graphTables[i]);
            }
            this.GraphTables = tables;

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            for (int i = 0; i < nodes.Count; i++)
            {
                if (nodes[i].IsInner)
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