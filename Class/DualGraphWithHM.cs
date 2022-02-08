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
        public List<GraphNode> UndividedGraphNodes { get; set; }
        public List<List<int>> UndividedGraphTable { get; set; }

        public PlanktonMesh PlanktonMesh { get; set; }
        
        public PlanktonMesh DualPlanktonMesh { get; set; }
        public PlanktonMesh IntegrateDualPlanktonMesh { get; set; }

        public SortedDictionary<int,GraphNode> Volume_VolumeNode { get; set; }
        public SortedDictionary<int,List<int>> Volume_Inner { get; set; }
        public SortedDictionary<int,int> Inner_DFace { get; set; }
        public SortedDictionary<int,List<int>> DFace_DVertice { get; set; }

        public SortedDictionary<int,GraphNode> Inner_InnerNode { get; set; }
        public SortedDictionary<int,GraphNode> DFace_InnerNode { get; set; }
        public SortedDictionary<int,List<int>> DVertice_Inner { get; set; }
        public SortedDictionary<int,List<int>> DVertice_Volume { get; set; }

        public SortedDictionary<int,List<int>> Volume_DVertice { get; set; }

        public SortedDictionary<int,int> Volume_VFace { get; set; }

        public List<List<int>> DFaceIndexsAroundOuterNodes { get; }


        public List<int> FaceCorrespondingGraphNodeIndex { get; set; }


        /// <summary>
        /// 构造空的DualGraphWithHM对象
        /// </summary>
        public DualGraphWithHM() : base()
        {
            this.PlanktonMesh = new PlanktonMesh();
            this.DualPlanktonMesh = new PlanktonMesh();
            this.IntegrateDualPlanktonMesh = new PlanktonMesh();
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="dualGraphWithHM"></param>
        public DualGraphWithHM(DualGraphWithHM source)
        {
            #region 父类属性
            // 深拷贝GraphTables
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < source.GraphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(source.GraphTables[i]);
            }
            // 深拷贝GraphNodes
            this.GraphNodes = new List<GraphNode>();
            for (int i = 0; i < source.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new GraphNode(source.GraphNodes[i]));
            }

            this.InnerNodeCount = source.InnerNodeCount;
            this.OuterNodeCount = source.OuterNodeCount;
            this.InnerNodeIndexList = source.InnerNodeIndexList;
            this.OuterNodeIndexList = source.OuterNodeIndexList;
            #endregion


            // 深拷贝PlanktonMesh
            this.DualPlanktonMesh = new PlanktonMesh(source.DualPlanktonMesh);
            // 深拷贝IntegratePlanktonMesh
            if (source.IntegrateDualPlanktonMesh != null)
            {
                this.IntegrateDualPlanktonMesh = new PlanktonMesh(source.IntegrateDualPlanktonMesh);
            }

            if (source.UndividedGraphNodes != null && source.UndividedGraphTable != null)
            {
                this.UndividedGraphNodes = new List<GraphNode>();
                this.UndividedGraphNodes.AddRange(source.UndividedGraphNodes);
                this.UndividedGraphTable = new List<List<int>>();
                for (int i = 0; i < source.UndividedGraphTable.Count; i++)
                {
                    this.UndividedGraphTable.Add(new List<int>());
                    this.UndividedGraphTable[i].AddRange(source.UndividedGraphTable[i]);
                }
            }
            

            if (source.Volume_VolumeNode != null)
            {
                this.Volume_VolumeNode = new SortedDictionary<int, GraphNode>();
                foreach (var item in source.Volume_VolumeNode)
                {
                    this.Volume_VolumeNode.Add(item.Key, item.Value);
                }
            }

            if (source.Volume_Inner != null)
            {
                this.Volume_Inner = new SortedDictionary<int, List<int>>();
                foreach (var item in source.Volume_Inner)
                {
                    this.Volume_Inner.Add(item.Key, item.Value);
                }
            }

            if (source.Inner_DFace != null)
            {
                this.Inner_DFace = new SortedDictionary<int, int>();
                foreach (var item in source.Inner_DFace)
                {
                    this.Inner_DFace.Add(item.Key, item.Value);
                }
            }

            if (source.DFace_DVertice != null)
            {
                this.DFace_DVertice = new SortedDictionary<int, List<int>>();
                foreach (var item in source.DFace_DVertice)
                {
                    this.DFace_DVertice.Add(item.Key, item.Value);
                }
            }

            if (source.DVertice_Inner != null)
            {
                this.DVertice_Inner = new SortedDictionary<int, List<int>>();
                foreach (var item in source.DVertice_Inner)
                {
                    this.DVertice_Inner.Add(item.Key, item.Value);
                }
            }

            if (source.DVertice_Volume != null)
            {
                this.DVertice_Volume = new SortedDictionary<int, List<int>>();
                foreach (var item in source.DVertice_Volume)
                {
                    this.DVertice_Volume.Add(item.Key, item.Value);
                }
            }

            if (source.Volume_DVertice != null)
            {
                this.Volume_DVertice = new SortedDictionary<int, List<int>>();
                foreach (var item in source.Volume_DVertice)
                {
                    this.Volume_DVertice.Add(item.Key, item.Value);
                }
            }

            if (source.Inner_InnerNode != null)
            {
                this.Inner_InnerNode = new SortedDictionary<int, GraphNode>();
                foreach (var item in source.Inner_InnerNode)
                {
                    this.Inner_InnerNode.Add(item.Key, item.Value);
                }
            }

            if (source.DFace_InnerNode != null)
            {
                this.DFace_InnerNode = new SortedDictionary<int, GraphNode>();
                foreach (var item in source.DFace_InnerNode)
                {
                    this.DFace_InnerNode.Add(item.Key, item.Value);
                }
            }

            if (source.Volume_VFace != null)
            {
                this.Volume_VFace = new SortedDictionary<int, int>();
                foreach (var item in source.Volume_VFace)
                {
                    this.Volume_VFace.Add(item.Key, item.Value);
                }
            }

            // 深拷贝FaceIndexsAroundOuterNodes
            if (source.DFaceIndexsAroundOuterNodes != null)
            {
                this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
                for (int i = 0; i < source.DFaceIndexsAroundOuterNodes.Count; i++)
                {
                    this.DFaceIndexsAroundOuterNodes.Add(new List<int>());
                    this.DFaceIndexsAroundOuterNodes[i].AddRange(source.DFaceIndexsAroundOuterNodes[i]);
                }
            }
        }

        public DualGraphWithHM(PlanktonMesh dualPlanktonMesh,
                               List<GraphNode> graphNodes,
                               List<List<int>> graphTables,
                               List<List<int>> faceIndexsAroundOuterNodes)
        {
            #region 父类属性
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }
            this.GraphNodes = new List<GraphNode>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                this.GraphNodes.Add(new GraphNode(graphNodes[i]));
            }
            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            this.Inner_InnerNode = new SortedDictionary<int, GraphNode>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                if (graphNodes[i].IsInner)
                {
                    this.InnerNodeCount++;
                    this.InnerNodeIndexList.Add(i);
                    this.Inner_InnerNode.Add(i, graphNodes[i]);
                }
                else
                {
                    this.OuterNodeCount++;
                    this.OuterNodeIndexList.Add(i);
                }
            }
            #endregion

            // 深拷贝PlanktonMesh
            this.DualPlanktonMesh = new PlanktonMesh(dualPlanktonMesh);

            this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
            for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
            {
                this.DFaceIndexsAroundOuterNodes.Add(new List<int>());
                this.DFaceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNodes[i]);
            }
        }

        /// <summary>
        /// 用于构造DualGraphWithHM
        /// </summary>
        /// <param name="dualPlanktonMesh"></param>
        /// <param name="decomposedPGraphNodes"></param>
        /// <param name="decomposedPGraphTables"></param>
        /// <param name="volume_volumeNode"></param>
        /// <param name="volume_inner"></param>
        /// <param name="inner_dFace"></param>
        /// <param name="faceIndexsAroundOuterNodes"></param>
        public DualGraphWithHM(PlanktonMesh dualPlanktonMesh,

                               List<GraphNode> decomposedPGraphNodes,
                               List<List<int>> decomposedPGraphTables,

                               SortedDictionary<int,GraphNode> volume_volumeNode,
                               SortedDictionary<int,List<int>> volume_inner,
                               SortedDictionary<int, int> inner_dFace,

                               List<List<int>> faceIndexsAroundOuterNodes
                               
                               ) 
        {
            #region 父类属性
            // 深拷贝GraphTables
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < decomposedPGraphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(decomposedPGraphTables[i]);
            }
            // 深拷贝GraphNodes
            this.GraphNodes = new List<GraphNode>();
            for (int i = 0; i < decomposedPGraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new GraphNode(decomposedPGraphNodes[i]));
            }

            this.InnerNodeCount = 0;
            this.OuterNodeCount = 0;
            this.InnerNodeIndexList = new List<int>();
            this.OuterNodeIndexList = new List<int>();
            this.Inner_InnerNode = new SortedDictionary<int, GraphNode>();
            for (int i = 0; i < decomposedPGraphNodes.Count; i++)
            {
                if (decomposedPGraphNodes[i].IsInner)
                {
                    this.InnerNodeCount++;
                    this.InnerNodeIndexList.Add(i);
                    this.Inner_InnerNode.Add(i, decomposedPGraphNodes[i]);
                }
                else
                {
                    this.OuterNodeCount++;
                    this.OuterNodeIndexList.Add(i);
                }
            }
            #endregion

            // 深拷贝PlanktonMesh
            this.DualPlanktonMesh = new PlanktonMesh(dualPlanktonMesh);


            // 构造Volume_VolumeNode
            this.Volume_VolumeNode = new SortedDictionary<int, GraphNode>();
            foreach (var item in volume_volumeNode)
            {
                this.Volume_VolumeNode.Add(item.Key, item.Value);
            }
            // 构造Volume_Inner
            this.Volume_Inner = new SortedDictionary<int, List<int>>();
            foreach (var item in volume_inner)
            {
                this.Volume_Inner.Add(item.Key, item.Value);
            }
            // 构造Inner_DFace
            this.Inner_DFace = new SortedDictionary<int, int>();
            foreach (var item in inner_dFace)
            {
                this.Inner_DFace.Add(item.Key, item.Value);
            }
            // 构造DFace_DVertice
            this.DFace_DVertice = new SortedDictionary<int, List<int>>();
            for (int i = 0; i < dualPlanktonMesh.Faces.Count; i++)
            {
                this.DFace_DVertice.Add(i, dualPlanktonMesh.Faces.GetFaceVertices(i).ToList());
            }


            this.DFace_InnerNode = new SortedDictionary<int, GraphNode>();
            foreach (var item1 in inner_dFace)
            {
                foreach (var item2 in Inner_InnerNode)
                {
                    if (item1.Key == item2.Key)
                    {
                        this.DFace_InnerNode.Add(item1.Value, item2.Value);
                    }
                }
            }



            SortedDictionary<int, List<int>> dVertice_dFace = new SortedDictionary<int, List<int>>();
            for (int i = 0; i < dualPlanktonMesh.Vertices.Count; i++)
            {
                foreach (var item in this.DFace_DVertice)
                {
                    if (item.Value.Contains(i))
                    {
                        if (dVertice_dFace.ContainsKey(i))
                        {
                            dVertice_dFace[i].Add(item.Key);
                        }
                        else
                        {
                            dVertice_dFace.Add(i, new List<int>());
                            dVertice_dFace[i].Add(item.Key);
                        }
                    }
                }
            }

            this.DVertice_Inner = new SortedDictionary<int, List<int>>();
            foreach (var item1 in dVertice_dFace)
            {
                for (int i = 0; i < item1.Value.Count; i++)
                {
                    foreach (var item2 in inner_dFace)
                    {
                        if (item2.Value == item1.Value[i])
                        {
                            if (this.DVertice_Inner.ContainsKey(item1.Key))
                            {
                                this.DVertice_Inner[item1.Key].Add(item2.Key);
                            }
                            else
                            {
                                this.DVertice_Inner.Add(item1.Key, new List<int>());
                                this.DVertice_Inner[item1.Key].Add(item2.Key);
                            }
                        }
                    }
                }
            }

            this.DVertice_Volume = new SortedDictionary<int, List<int>>();
            foreach (var item1 in this.DVertice_Inner)
            {
                for (int i = 0; i < item1.Value.Count; i++)
                {
                    foreach (var item2 in volume_inner)
                    {
                        if (item2.Value.Contains(item1.Value[i]))
                        {
                            if (this.DVertice_Volume.ContainsKey(item1.Key))
                            {
                                if (!this.DVertice_Volume[item1.Key].Contains(item2.Key))
                                {
                                    this.DVertice_Volume[item1.Key].Add(item2.Key);
                                }
                            }
                            else
                            {
                                this.DVertice_Volume.Add(item1.Key, new List<int>());
                                this.DVertice_Volume[item1.Key].Add(item2.Key);
                            }
                        }
                    }
                }
            }

            this.Volume_DVertice = new SortedDictionary<int, List<int>>();
            foreach (var item in this.Volume_VolumeNode)
            {
                List<int> dVertice = new List<int>();
                foreach (var item2 in this.DVertice_Volume)
                {
                    if (item2.Value.Contains(item.Key))
                    {
                        dVertice.Add(item2.Key);
                    }
                }
                this.Volume_DVertice.Add(item.Key, dVertice);
            }

            SortedDictionary<int, List<int>> volume_dFace = new SortedDictionary<int, List<int>>();
            foreach (var item1 in this.Volume_Inner)
            {
                volume_dFace.Add(item1.Key, new List<int>());
                for (int i = 0; i < item1.Value.Count; i++)
                {
                    volume_dFace[item1.Key].Add(this.Inner_DFace[item1.Value[i]]);
                }
            }



            PlanktonMesh integratePlanktonMesh = new PlanktonMesh(this.DualPlanktonMesh);

            List<int> halfedgeForMerge = new List<int>();
            foreach (var item in volume_dFace)
            {
                List<int> publicHalfedge = this.FindPublicHalfedge(this.DualPlanktonMesh, item.Value);
                halfedgeForMerge.AddRange(publicHalfedge);
            }


            for (int i = 0; i < halfedgeForMerge.Count; i++)
            {
                integratePlanktonMesh.Faces.MergeFaces(halfedgeForMerge[i]);
            }
            integratePlanktonMesh.Compact();
            this.IntegrateDualPlanktonMesh = integratePlanktonMesh;


            this.Volume_VFace = new SortedDictionary<int, int>();
            foreach (var item in this.Volume_DVertice)
            {
                for (int i = 0; i < this.IntegrateDualPlanktonMesh.Faces.Count; i++)
                {
                    int[] set = this.IntegrateDualPlanktonMesh.Faces.GetFaceVertices(i);
                    if (set.All(t => item.Value.Any(b => b == t)))
                    {
                        this.Volume_VFace.Add(item.Key, i);
                    }
                }
            }


            this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
            for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
            {
                this.DFaceIndexsAroundOuterNodes.Add(new List<int>());
                this.DFaceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNodes[i]);
            }
        }

        public List<int> FindPublicHalfedge(PlanktonMesh p, List<int> faceIndexs)
        {
            List<int[]> facePairs = new List<int[]>();
            for (int i = 0; i < faceIndexs.Count; i++)
            {
                for (int j = i + 1; j < faceIndexs.Count; j++)
                {
                    facePairs.Add(new int[2] { faceIndexs[i], faceIndexs[j] });
                }
            }

            List<int> publicHalfedge = new List<int>();
            for (int i = 0; i < facePairs.Count; i++)
            {
                int secondFace = facePairs[i][1];
                int[] secondFaceEdges = p.Faces.GetHalfedges(secondFace);
                foreach (int firstFaceEdge in p.Faces.GetHalfedges(facePairs[i][0]))
                {
                    if (secondFaceEdges.Contains(p.Halfedges.GetPairHalfedge(firstFaceEdge)))
                    {
                        publicHalfedge.Add(firstFaceEdge);
                    }
                }
            }

            return publicHalfedge;
        }

    }
}
