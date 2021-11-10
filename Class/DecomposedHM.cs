﻿using Grasshopper.Kernel;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
//using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class DecomposedHM
    {
        /// <summary>
        /// 对应的图结构的邻接表
        /// </summary>
        public List<List<int>> GraphLoL { get; set; }

        /// <summary>
        /// 对应的图结构中的Node
        /// </summary>
        public List<Node> GraphNodes { get; set; }

        /// <summary>
        /// 对应的半边网格
        /// </summary>
        public PlanktonMesh PlanktonMesh { get; set; }


        public int DecomposeTheithVertex { get; set; }

        public int[,] DecomposeTheithPairHFVertexIndex { get; set; }

        public int DecomposeTheithPairHFResult { get; set; }

        public string TreeNodeName { get; set; }

        // public List<int> VertexIndexHasBeenDecomposed { get; set; }

        public DecomposedHM()
        {
            this.GraphLoL = new List<List<int>>();
            this.GraphNodes = new List<Node>();
            this.PlanktonMesh = new PlanktonMesh();
        }

        public DecomposedHM(PlanktonMesh planktonMesh, 
                            List<List<int>> graphLoL, 
                            List<Node> graphNode,
                            int decomposeTheithVertex,
                            int[,] decomposeTheithPairHFVertexIndex,
                            int decomposeTheithPairHFResult)
        {
            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(planktonMesh);
            // 深拷贝List<List<int>>
            this.GraphLoL = new List<List<int>>();
            for (int i = 0; i < graphLoL.Count; i++)
            {
                this.GraphLoL.Add(new List<int>());
                this.GraphLoL[i].AddRange(graphLoL[i]);
            }
            // 深拷贝Node
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < graphNode.Count; i++)
            {
                this.GraphNodes.Add(new Node(graphNode[i]));
            }

            this.DecomposeTheithVertex = decomposeTheithVertex;
            this.DecomposeTheithPairHFVertexIndex = decomposeTheithPairHFVertexIndex;
            this.DecomposeTheithPairHFResult = decomposeTheithPairHFResult;

            this.TreeNodeName = "分裂了第" + decomposeTheithVertex.ToString() + "个Vertex" +
                                "，所选的两个半边的起点终点是" + 
                                decomposeTheithPairHFVertexIndex[0, 0].ToString() + ";" +
                                decomposeTheithPairHFVertexIndex[0, 1].ToString() + ";" +
                                decomposeTheithPairHFVertexIndex[1, 0].ToString() + ";" +
                                decomposeTheithPairHFVertexIndex[1, 1].ToString() +
                                "，并且是第" + decomposeTheithPairHFResult.ToString() + "分裂可能性";

            //this.VertexIndexHasBeenDecomposed = new List<int>();
            //this.VertexIndexHasBeenDecomposed.Add(decomposeTheithVertex);
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="decomposedHM"></param>
        public DecomposedHM(DecomposedHM decomposedHM)
        {
            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(decomposedHM.PlanktonMesh);
            // 深拷贝List<List<int>>
            this.GraphLoL = new List<List<int>>();
            for (int i = 0; i < decomposedHM.GraphLoL.Count; i++)
            {
                this.GraphLoL.Add(new List<int>());
                this.GraphLoL[i].AddRange(decomposedHM.GraphLoL[i]);
            }
            // 深拷贝Node
            this.GraphNodes = new List<Node>();
            for (int i = 0; i < decomposedHM.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(decomposedHM.GraphNodes[i]));
            }

            this.DecomposeTheithVertex = decomposedHM.DecomposeTheithVertex;
            this.DecomposeTheithPairHFVertexIndex = (int[,])decomposedHM.DecomposeTheithPairHFVertexIndex.Clone();
            this.DecomposeTheithPairHFResult = decomposedHM.DecomposeTheithPairHFResult;
            this.TreeNodeName = (string)decomposedHM.TreeNodeName.Clone();
            //this.VertexIndexHasBeenDecomposed = new List<int>();
            //for (int i = 0; i < decomposedHM.VertexIndexHasBeenDecomposed.Count; i++)
            //{
            //    this.VertexIndexHasBeenDecomposed.Add(decomposedHM.VertexIndexHasBeenDecomposed[i]);
            //}
        }

        //public object Clone()
        //{
        //    PlanktonMesh planktonMesh = (PlanktonMesh)MemberwiseClone();
        //}
    }
}