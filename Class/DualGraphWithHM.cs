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

        public PlanktonMesh IntegratePlanktonMesh { get; set; }


        public List<int> VolumeNodesIndex { get; set; }

        // public List<int> DualFaceCorrespondingInnerNodeIndexs { get; set; }

        // public List<int> InnerNodeCorrespondingDualFace { get; set; }

        public Matrix DFaceInnerNodeTable { get; set; }

        public List<int> DFaceInnerNodeList { get; set; }
        public List<int> InnerNodeDFaceList { get; set; }

        public Matrix DFaceVolumeNodeTable { get; set; }

        // public List<List<int>> InnerNodeCorrespondingDualFaceVertices { get; set; }

        // public List<List<int>> VolumeNodeCorrespondingDualFaceVertices { get; set; }
        public Matrix DVerticeInnerNodeTable { get; set; }

        public Matrix DVerticeVolumeNodeTable { get; set; }


        public List<GraphNode> DualFaceCorrespondingInnerNodes { get; set; }
        public List<GraphNode> IntegrateFaceCorrespondingVolumeNodes { get; set; }

        // public List<int> IntegrateFaceCorrespondingVolumeNodeIndexs { get; set; }



        public List<List<int>> DFaceIndexsAroundOuterNodes { get; }

        //public List<List<int>> DVertexBelongsToWhichInnerNode { get; }

        //public Dictionary<int,List<int>> DVertexBelongsToWhichVolume { get; }



        /// <summary>
        /// 构造空的DualGraphWithHM对象
        /// </summary>
        public DualGraphWithHM() : base()
        {
            this.PlanktonMesh = new PlanktonMesh();
            this.IntegratePlanktonMesh = new PlanktonMesh();


            // this.DualFaceCorrespondingInnerNodeIndexs = new List<int>();

            // this.InnerNodeCorrespondingDualFace = new List<int>();
            // this.InnerNodeCorrespondingDualFaceVertices = new List<List<int>>();
            // this.VolumeNodeCorrespondingDualFaceVertices = new List<List<int>>();

            this.VolumeNodesIndex = null;
            
            this.DFaceInnerNodeTable = null;
            this.DFaceVolumeNodeTable = null;

            this.DVerticeInnerNodeTable = null;
            this.DVerticeVolumeNodeTable = null;


            this.DualFaceCorrespondingInnerNodes = new List<GraphNode>();
            this.IntegrateFaceCorrespondingVolumeNodes = new List<GraphNode>();
            // this.IntegrateFaceCorrespondingVolumeNodeIndexs = new List<int>();
            

            this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
            //this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            //this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
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
            this.PlanktonMesh = new PlanktonMesh(source.PlanktonMesh);
            // 深拷贝IntegratePlanktonMesh
            this.IntegratePlanktonMesh = new PlanktonMesh(source.IntegratePlanktonMesh);


            //// 深拷贝DualFaceCorrespondingInnerNodeIndexs
            //this.DualFaceCorrespondingInnerNodeIndexs = new List<int>();
            //this.DualFaceCorrespondingInnerNodeIndexs.AddRange(source.DualFaceCorrespondingInnerNodeIndexs);

            //// 深拷贝InnerNodeCorrespondingDualFace
            //this.InnerNodeCorrespondingDualFace = new List<int>();
            //this.InnerNodeCorrespondingDualFace.AddRange(source.InnerNodeCorrespondingDualFace);


            this.VolumeNodesIndex = new List<int>();
            this.VolumeNodesIndex.AddRange(source.VolumeNodesIndex);

            // 深拷贝 Matrix DFaceInnerNodeTable
            this.DFaceInnerNodeTable = new Matrix(source.DFaceInnerNodeTable.RowCount, source.DFaceInnerNodeTable.ColumnCount);
            for (int i = 0; i < source.DFaceInnerNodeTable.RowCount; i++)
            {
                for (int j = 0; j < source.DFaceInnerNodeTable.ColumnCount; j++)
                {
                    this.DFaceInnerNodeTable[i, j] = source.DVerticeInnerNodeTable[i, j];
                }
            }
            // 深拷贝 Matrix DFaceVolumeNodeTable
            this.DFaceVolumeNodeTable = new Matrix(source.DFaceVolumeNodeTable.RowCount, source.DFaceVolumeNodeTable.ColumnCount);
            for (int i = 0; i < source.DFaceVolumeNodeTable.RowCount; i++)
            {
                for (int j = 0; j < source.DFaceVolumeNodeTable.ColumnCount; j++)
                {
                    this.DFaceVolumeNodeTable[i, j] = source.DFaceVolumeNodeTable[i, j];
                }
            }
            // 深拷贝 Matrix DVerticeInnerNodeTable
            this.DVerticeInnerNodeTable = new Matrix(source.DVerticeInnerNodeTable.RowCount, source.DVerticeInnerNodeTable.ColumnCount);
            for (int i = 0; i < source.DVerticeInnerNodeTable.RowCount; i++)
            {
                for (int j = 0; j < source.DVerticeInnerNodeTable.ColumnCount; j++)
                {
                    this.DVerticeInnerNodeTable[i, j] = source.DVerticeInnerNodeTable[i, j];
                }
            }
            // 深拷贝 Matrix DVerticeVolumeNodeTable
            this.DVerticeVolumeNodeTable = new Matrix(source.DVerticeVolumeNodeTable.RowCount, source.DVerticeVolumeNodeTable.ColumnCount);
            for (int i = 0; i < source.DVerticeVolumeNodeTable.RowCount; i++)
            {
                for (int j = 0; j < source.DVerticeVolumeNodeTable.ColumnCount; j++)
                {
                    this.DVerticeVolumeNodeTable[i, j] = source.DVerticeVolumeNodeTable[i, j];
                }
            }


            // 深拷贝 List<Graph> DualFaceCorrespondingInnerNodes
            this.DualFaceCorrespondingInnerNodes = new List<GraphNode>();
            for (int i = 0; i < source.DualFaceCorrespondingInnerNodes.Count; i++)
            {
                this.DualFaceCorrespondingInnerNodes.Add(new GraphNode(source.DualFaceCorrespondingInnerNodes[i]));
            }

            // 深拷贝 List<Graph> IntegrateFaceCorrespondingVolumeNodes
            this.IntegrateFaceCorrespondingVolumeNodes = new List<GraphNode>();
            for (int i = 0; i < source.IntegrateFaceCorrespondingVolumeNodes.Count; i++)
            {
                this.IntegrateFaceCorrespondingVolumeNodes.Add(new GraphNode(source.IntegrateFaceCorrespondingVolumeNodes[i]));
            }
            //// 深拷贝IntegrateFaceCorrespondingVolumeNodeIndexs
            //this.IntegrateFaceCorrespondingVolumeNodeIndexs = new List<int>();
            //this.IntegrateFaceCorrespondingVolumeNodeIndexs.AddRange(source.IntegrateFaceCorrespondingVolumeNodeIndexs);





            // 深拷贝FaceIndexsAroundOuterNodes
            this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
            for (int i = 0; i < source.DFaceIndexsAroundOuterNodes.Count; i++)
            {
                this.DFaceIndexsAroundOuterNodes.Add(new List<int>());
                this.DFaceIndexsAroundOuterNodes[i].AddRange(source.DFaceIndexsAroundOuterNodes[i]);
            }

            //// 深拷贝DVertexBelongsToWhichInnerNode
            //this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            //for (int i = 0; i < source.DVertexBelongsToWhichInnerNode.Count; i++)
            //{
            //    this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
            //    this.DVertexBelongsToWhichInnerNode[i].AddRange(source.DVertexBelongsToWhichInnerNode[i]);
            //}

            //// 深拷贝DVertexBelongsToWhichVolume
            //this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
            //foreach (KeyValuePair<int,List<int>> pair in source.DVertexBelongsToWhichVolume)
            //{
            //    this.DVertexBelongsToWhichVolume.Add(pair.Key, pair.Value);
            //}
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

                               List<int> faceCorrespondingInnerNodeIndexs,

                               List<int> volumeNodesIndex,
                               List<GraphNode> volumeNodes,
                               Dictionary<int,List<int>> volumeContainsWhichInnerNode,


                               List<List<int>> faceIndexsAroundOuterNodes
                               
                               ) 
        {
            #region 父类属性
            // 深拷贝GraphTables
            this.GraphTables = new List<List<int>>();
            for (int i = 0; i < graphTables.Count; i++)
            {
                this.GraphTables.Add(new List<int>());
                this.GraphTables[i].AddRange(graphTables[i]);
            }
            // 深拷贝GraphNodes
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
            #endregion

            this.VolumeNodesIndex = new List<int>();
            this.VolumeNodesIndex.AddRange(volumeNodesIndex);

            // 深拷贝PlanktonMesh
            this.PlanktonMesh = new PlanktonMesh(planktonMesh);

            // 构造DFaceInnerNodeTable，表头是DFace的序号列表和InnerNodeIndexList列表

            this.DFaceInnerNodeList = faceCorrespondingInnerNodeIndexs;
            this.InnerNodeDFaceList = new List<int>();
            foreach (var item in faceCorrespondingInnerNodeIndexs)
            {
                this.InnerNodeDFaceList.Add(faceCorrespondingInnerNodeIndexs.IndexOf(item));
            }

            this.DFaceInnerNodeTable = new Matrix(this.PlanktonMesh.Faces.Count, this.InnerNodeCount);
            for (int i = 0; i < faceCorrespondingInnerNodeIndexs.Count; i++)
            {
                for (int j = 0; j < this.InnerNodeCount; j++)
                {
                    if (faceCorrespondingInnerNodeIndexs[i] == this.InnerNodeIndexList[j])
                    {
                        this.DFaceInnerNodeTable[i, j] = 1;
                    }
                    else
                    {
                        this.DFaceInnerNodeTable[i, j] = 0;
                    }
                }
            }

            /////debug
            List<string> debugDFaceInnerNodeTable = DualGraphWithHM.PrintMatrix(this.DFaceInnerNodeTable);


            // 根据InnerNode所对应的对偶Face序号，求得InnerNode所对应的对偶FaceVertice序号
            List<List<int>> innerNodeCorrespondingDualFaceVertices = new List<List<int>>();
            for (int i = 0; i < this.InnerNodeIndexList.Count; i++)
            {
                innerNodeCorrespondingDualFaceVertices.Add(new List<int>());
                int debugIndex = this.InnerNodeCorrespondingWhichDFace()[i];
                innerNodeCorrespondingDualFaceVertices[i].AddRange(this.PlanktonMesh.Faces.GetFaceVertices(debugIndex));
            }

            this.DVerticeInnerNodeTable = new Matrix(this.PlanktonMesh.Vertices.Count, this.InnerNodeCount);
            for (int i = 0; i < innerNodeCorrespondingDualFaceVertices.Count; i++)
            {
                for (int j = 0; j < this.PlanktonMesh.Vertices.Count; j++)
                {
                    if (innerNodeCorrespondingDualFaceVertices[i].Contains(j))
                    {
                        this.DVerticeInnerNodeTable[j, i] = 1;
                    }
                    else
                    {
                        this.DVerticeInnerNodeTable[j, i] = 0;
                    }
                }
            }

            /////debug
            List<string> debugDVerticeInnerNodeTable = DualGraphWithHM.PrintMatrix(this.DVerticeInnerNodeTable);

            // 根据FaceCorrespondingInnerNodeIndexs，求得FaceCorrespondingInnerNodes
            this.DualFaceCorrespondingInnerNodes = new List<GraphNode>();
            for (int i = 0; i < faceCorrespondingInnerNodeIndexs.Count; i++)
            {
                this.DualFaceCorrespondingInnerNodes.Add(this.GraphNodes[faceCorrespondingInnerNodeIndexs[i]]);
            }



            List<List<int>> volumeNodeCorrespondingInnerNodes = new List<List<int>>();
            for (int i = 0; i < volumeNodesIndex.Count; i++)
            {
                volumeNodeCorrespondingInnerNodes.Add(new List<int>());
                volumeNodeCorrespondingInnerNodes[i].AddRange(volumeContainsWhichInnerNode[volumeNodesIndex[i]]);
            }

            List<List<int>> volumeNodeCorrespondingDFace = new List<List<int>>();
            for (int i = 0; i < volumeNodeCorrespondingInnerNodes.Count; i++)
            {
                volumeNodeCorrespondingDFace.Add(new List<int>());
                for (int j = 0; j < volumeNodeCorrespondingInnerNodes[i].Count; j++)
                {
                    List<int> innerNodeCorrespondingWhichDFace = this.InnerNodeCorrespondingWhichDFace();
                    volumeNodeCorrespondingDFace[i].Add(innerNodeCorrespondingWhichDFace[InnerNodeIndexList.IndexOf(volumeNodeCorrespondingInnerNodes[i][j])]);
                }
            }

            // 构造DFaceVolumeNodeTable，表头是DFace的序号列表和volumeNodesIndex列表
            this.DFaceVolumeNodeTable = new Matrix(this.PlanktonMesh.Faces.Count, volumeNodes.Count);
            for (int i = 0; i < volumeNodesIndex.Count; i++)
            {
                for (int j = 0; j < this.PlanktonMesh.Faces.Count; j++)
                {
                    if (volumeNodeCorrespondingInnerNodes[volumeNodesIndex[i]].Contains(faceCorrespondingInnerNodeIndexs[j]))
                    {
                        this.DFaceVolumeNodeTable[j, volumeNodesIndex[i]] = 1;
                    }
                    else
                    {
                        this.DFaceVolumeNodeTable[j, volumeNodesIndex[i]] = 0;
                    }
                }
            }

            /////debug
            List<string> debugDFaceVolumeNodeTable = DualGraphWithHM.PrintMatrix(this.DFaceVolumeNodeTable);


            List<List<int>> volumeNodeCorrespondingDVertice = new List<List<int>>();
            for (int i = 0; i < volumeNodeCorrespondingInnerNodes.Count; i++)
            {
                volumeNodeCorrespondingDVertice.Add(new List<int>());
                for (int j = 0; j < volumeNodeCorrespondingInnerNodes[i].Count; j++)
                {
                    int innerNodeIndex = volumeNodeCorrespondingInnerNodes[i][j];
                    int faceIndex = faceCorrespondingInnerNodeIndexs.IndexOf(innerNodeIndex);
                    volumeNodeCorrespondingDVertice[i].AddRange(this.PlanktonMesh.Faces.GetFaceVertices(faceIndex));
                }
                volumeNodeCorrespondingDVertice[i] = volumeNodeCorrespondingDVertice[i].Distinct().ToList();
            }

            this.DVerticeVolumeNodeTable = new Matrix(this.PlanktonMesh.Vertices.Count, volumeNodesIndex.Count);
            for (int i = 0; i < this.PlanktonMesh.Vertices.Count; i++)
            {
                for (int j = 0; j < volumeNodeCorrespondingDVertice.Count; j++)
                {
                    if (volumeNodeCorrespondingInnerNodes[j].Contains(i))
                    {
                        this.DVerticeVolumeNodeTable[i, j] = 1;
                    }
                    else
                    {
                        this.DVerticeVolumeNodeTable[i, j] = 0;
                    }
                }
            }

            /////debug
            List<string> debugDVerticeVolumeNodeTable = DualGraphWithHM.PrintMatrix(this.DVerticeVolumeNodeTable);

            //List<List<int>> volumeNodeCorrespondingDualFaceVertices = new List<List<int>>();
            //for (int i = 0; i < volumeNodeCorrespondingInnerNodes.Count; i++)
            //{
            //    volumeNodeCorrespondingDualFaceVertices.Add(new List<int>());
            //    for (int j = 0; j < volumeNodeCorrespondingInnerNodes[i].Count; j++)
            //    {
            //        int debug = volumeNodeCorrespondingInnerNodes[i][j];
            //        List<List<int>> debug2 = volumeNodeCorrespondingDVertice;
            //        List<int> dVertice = debug2[debug];
            //        volumeNodeCorrespondingDualFaceVertices[i].AddRange(dVertice);
            //    }
            //}


            PlanktonMesh integratePlanktonMesh = new PlanktonMesh(this.PlanktonMesh);



            List<int> halfedgeForMerge = new List<int>();
            for (int i = 0; i < volumeNodeCorrespondingDFace.Count; i++)
            {
                List<int> publicHalfedge = this.FindPublicHalfedge(this.PlanktonMesh, volumeNodeCorrespondingDFace[i]);
                halfedgeForMerge.AddRange(publicHalfedge);
            }


            for (int i = 0; i < halfedgeForMerge.Count; i++)
            {
                integratePlanktonMesh.Faces.MergeFaces(halfedgeForMerge[i]);
            }
            integratePlanktonMesh.Compact();
            this.IntegratePlanktonMesh = integratePlanktonMesh;


            this.IntegrateFaceCorrespondingVolumeNodes = new List<GraphNode>();
            this.IntegrateFaceCorrespondingVolumeNodes.AddRange(volumeNodes);

            //List<PlanktonXYZ> planktonVertex = new List<PlanktonXYZ>();
            //for (int i = 0; i < this.PlanktonMesh.Vertices.Count; i++)
            //{
            //    planktonVertex.Add(this.PlanktonMesh.Vertices[i].ToXYZ());
            //}

            //this.IntegratePlanktonMesh = new PlanktonMesh();
            //this.IntegratePlanktonMesh.Vertices.AddVertices(planktonVertex);
            //this.IntegratePlanktonMesh.Faces.AddFaces(volumeNodeCorrespondingDVertice);



            //// 深拷贝FaceIndexsAroundOuterNodes
            //this.DFaceIndexsAroundOuterNodes = new List<List<int>>();
            //for (int i = 0; i < faceIndexsAroundOuterNodes.Count; i++)
            //{
            //    this.DFaceIndexsAroundOuterNodes.Add(new List<int>());
            //    this.DFaceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNodes[i]);
            //}





            //// 深拷贝DVertexBelongsToWhichInnerNode
            //this.DVertexBelongsToWhichInnerNode = new List<List<int>>();
            //for (int i = 0; i < dVertexBelongsToWhichInnerNode.Count; i++)
            //{
            //    this.DVertexBelongsToWhichInnerNode.Add(new List<int>());
            //    this.DVertexBelongsToWhichInnerNode[i].AddRange(dVertexBelongsToWhichInnerNode[i]);
            //}

            //// 深拷贝DVertexBelongsToWhichVolume
            //this.DVertexBelongsToWhichVolume = new Dictionary<int, List<int>>();
            //foreach (KeyValuePair<int, List<int>> pair in dVertexBelongsToWhichVolume)
            //{
            //    this.DVertexBelongsToWhichVolume.Add(pair.Key, pair.Value);
            //}
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

        public List<int> InnerNodeCorrespondingWhichDFace()
        {
            List<int> innerNodeCorrespondingWhichDFace = new List<int>();
            for (int i = 0; i < this.DFaceInnerNodeTable.RowCount; i++)
            {
                for (int j = 0; j < this.DFaceInnerNodeTable.ColumnCount; j++)
                {
                    if (this.DFaceInnerNodeTable[i, j] == 1)
                    {
                        innerNodeCorrespondingWhichDFace.Add(i);
                    }
                }
            }
            return innerNodeCorrespondingWhichDFace;
        }


        public List<List<int>> DVertexBelongsToWhichInnerNode()
        {
            List<List<int>> dVertexBelongsToWhichInnerNode = new List<List<int>>();
            for (int i = 0; i < this.DVerticeInnerNodeTable.RowCount; i++)
            {
                dVertexBelongsToWhichInnerNode.Add(new List<int>());
                for (int j = 0; j < this.DVerticeInnerNodeTable.ColumnCount; j++)
                {
                    if (this.DVerticeInnerNodeTable[i, j] == 1)
                    {
                        dVertexBelongsToWhichInnerNode[i].Add(this.InnerNodeIndexList[j]);
                    }
                }
            }
            return dVertexBelongsToWhichInnerNode;
        }

        public List<List<int>> DVertexBelongsToWhichVolumeNode()
        {
            List<List<int>> dVertexBelongsToWhichVolumeNode = new List<List<int>>();
            for (int i = 0; i < this.DVerticeVolumeNodeTable.RowCount; i++)
            {
                dVertexBelongsToWhichVolumeNode.Add(new List<int>());
                for (int j = 0; j < this.DVerticeVolumeNodeTable.ColumnCount; j++)
                {
                    if (this.DVerticeVolumeNodeTable[i, j] == 1)
                    {
                        dVertexBelongsToWhichVolumeNode[i].Add(this.VolumeNodesIndex[j]);
                    }
                }
            }
            return dVertexBelongsToWhichVolumeNode;
        }

        public List<int> DFaceBelongsToWhichVolumeNode()
        {
            List<int> dFaceBelongsToWhichVolumeNode = new List<int>();
            for (int i = 0; i < this.DFaceVolumeNodeTable.RowCount; i++)
            {
                for (int j = 0; j < this.DFaceVolumeNodeTable.ColumnCount; j++)
                {
                    if (this.DFaceVolumeNodeTable[i,j] == 1)
                    {
                        dFaceBelongsToWhichVolumeNode.Add(this.VolumeNodesIndex[j]);
                    }
                }
            }
            return dFaceBelongsToWhichVolumeNode;
        }


        public static List<string> PrintMatrix(Matrix M)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < M.RowCount; i++)
            {
                string str = string.Format("{0}行:", i);
                for (int j = 0; j < M.ColumnCount; j++)
                {
                    str += string.Format("{0},", M[i, j]);
                }
                output.Add(str);
            }

            return output;
        }
    }
}
