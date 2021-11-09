using Grasshopper.Kernel;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph
{
    public class DecomposedHalfedgeMesh
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

        public DecomposedHalfedgeMesh(PlanktonMesh planktonMesh, 
                                      List<List<int>> graphLoL, 
                                      List<Node> graphNode)
        {
            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(planktonMesh);
            // 深拷贝List<List<int>>
            for (int i = 0; i < graphLoL.Count; i++)
            {
                this.GraphLoL.Add(new List<int>());
                this.GraphLoL[i].AddRange(graphLoL[i]);
            }
            // 深拷贝Node
            for (int i = 0; i < graphNode.Count; i++)
            {
                this.GraphNodes.Add(new Node(graphNode[i]));
            }
        }

        /// <summary>
        /// 利用拷贝构造函数实现深拷贝
        /// </summary>
        /// <param name="decomposedHM"></param>
        public DecomposedHalfedgeMesh(DecomposedHalfedgeMesh decomposedHM)
        {
            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(decomposedHM.PlanktonMesh);
            // 深拷贝List<List<int>>
            for (int i = 0; i < decomposedHM.GraphLoL.Count; i++)
            {
                this.GraphLoL.Add(new List<int>());
                this.GraphLoL[i].AddRange(decomposedHM.GraphLoL[i]);
            }
            // 深拷贝Node
            for (int i = 0; i < decomposedHM.GraphNodes.Count; i++)
            {
                this.GraphNodes.Add(new Node(decomposedHM.GraphNodes[i]));
            }
        }

        //public object Clone()
        //{
        //    PlanktonMesh planktonMesh = (PlanktonMesh)MemberwiseClone();
        //}
    }
}