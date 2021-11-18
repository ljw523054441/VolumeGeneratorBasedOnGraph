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

        public List<string> TreeNodeHistory { get; set; }

        public List<int> VolumeNodesIndex { get; set; }

        public List<GraphNode> VolumeNodes { get; set; }

        public Dictionary<int, List<int>> VolumeContainsWhichInnerNode { get; set; }


        public List<GraphNode> OriginGraphNodes { get; set; }
        public List<List<int>> OriginGraphTables { get; set; }


        /// <summary>
        /// 构造空的GraphWithHFMesh对象
        /// </summary>
        public GraphWithHM():base()
        {
            this.PlanktonMesh = new PlanktonMesh();

            this.TreeNodeLabel = null;

            this.TreeNodeHistory = null;
            this.VolumeNodesIndex = null;
            this.VolumeNodes = null;

            this.VolumeContainsWhichInnerNode = null;
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
            // 深拷贝TreeNodeLabel
            this.TreeNodeLabel = (string)source.TreeNodeLabel.Clone();


            // 深拷贝TreeNodeHistory，深拷贝VolumeNodesIndex，深拷贝VolumeNodes
            if (source.TreeNodeHistory != null && source.VolumeNodesIndex != null && source.VolumeNodes != null)
            {
                this.TreeNodeHistory = new List<string>();
                this.VolumeNodesIndex = new List<int>();
                this.VolumeNodes = new List<GraphNode>();
                for (int i = 0; i < source.TreeNodeHistory.Count; i++)
                {
                    this.TreeNodeHistory.Add(source.TreeNodeHistory[i]);
                    this.VolumeNodesIndex.Add(source.VolumeNodesIndex[i]);
                    this.VolumeNodes.Add(new GraphNode(source.VolumeNodes[i]));
                }
            }
            else
            {
                this.TreeNodeHistory = null;
                this.VolumeNodesIndex = null;
                this.VolumeNodes = null;
            }

            // 深拷贝VolumeContainsWhichInnerNode
            this.VolumeContainsWhichInnerNode = new Dictionary<int, List<int>>();
            if (source.VolumeContainsWhichInnerNode != null)
            {
                foreach (KeyValuePair<int, List<int>> pair in source.VolumeContainsWhichInnerNode)
                {
                    this.VolumeContainsWhichInnerNode.Add(pair.Key, pair.Value);
                }
            }
        }

        /// <summary>
        /// 用于在还没有进行Decompose时，构造GraphWithHFMesh对象
        /// </summary>
        public GraphWithHM(PlanktonMesh planktonMesh,
                           List<GraphNode> graphNodes,
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
            // 设置空的TreeNodeLabel
            this.TreeNodeLabel = "尚未进行过Decompose的GraphWithHFMesh对象";
            // 设置空的TreeNodeHistory，设置空的VolumeNodesIndex，设置空的VolumeNode
            this.TreeNodeHistory = null;
            this.VolumeNodesIndex = null;
            this.VolumeNodes = null;
            // 设置空的VolumeContainsWhichInnerNode
            this.VolumeContainsWhichInnerNode = null;
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
                           List<GraphNode> graphNodes,
                           List<List<int>> graphTables, 
                           string treeNodeLabel,
                           Dictionary<int,List<int>> volumeContainsWhichInnerNode)
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
            // 设置新的TreeNodeLabel
            this.TreeNodeLabel = treeNodeLabel;
            // 设置空的TreeNodeHistory，设置空的VolumeNodesIndex，设置空的VolumeNode
            this.TreeNodeHistory = null;
            this.VolumeNodesIndex = null;
            this.VolumeNodes = null;
            // 设置新的VolumeContainsWhichInnerNode
            this.VolumeContainsWhichInnerNode = new Dictionary<int, List<int>>();
            foreach (KeyValuePair<int, List<int>> pair in volumeContainsWhichInnerNode)
            {
                this.VolumeContainsWhichInnerNode.Add(pair.Key, pair.Value);
            }
        }




        public static GraphWithHM ReTutteEmbeding(GraphWithHM source)
        {
            #region 计算相关矩阵的前置条件
            GraphWithHM newGraphWithHM = new GraphWithHM(source);

            // 需要修改的相关属性
            PlanktonMesh newPlanktonMesh = newGraphWithHM.PlanktonMesh;
            List<GraphNode> newGraphNodes = newGraphWithHM.GraphNodes;

            // 不需要修改的相关属性
            List<List<int>> originGraphTables = newGraphWithHM.GraphTables;
            string originTreeNodeLabel = newGraphWithHM.TreeNodeLabel;


            int innerNodeCount = newGraphWithHM.InnerNodeCount;
            int outerNodeCount = newGraphWithHM.OuterNodeCount;
            List<int> outerNodeIndexList = newGraphWithHM.OuterNodeIndexList;
            List<int> innerNodeIndexList = newGraphWithHM.InnerNodeIndexList;

            List<List<int>> triangleGraphTables = new List<List<int>>();
            List<int> innerNodeDegreeList = new List<int>();
            for (int i = 0; i < newPlanktonMesh.Vertices.Count; i++)
            {
                triangleGraphTables.Add(new List<int>());
                triangleGraphTables[i].AddRange(newPlanktonMesh.Vertices.GetVertexNeighbours(i));
                if (newGraphNodes[i].IsInner)
                {
                    innerNodeDegreeList.Add(newPlanktonMesh.Vertices.GetVertexNeighbours(i).Length);
                }
            }
            #endregion

            #region 计算矩阵 P(outer)
            Matrix P_outer = new Matrix(outerNodeCount, 2);
            for (int i = 0; i < newGraphWithHM.GraphNodes.Count; i++)
            {
                if (!newGraphWithHM.GraphNodes[i].IsInner)
                {
                    P_outer[outerNodeIndexList.IndexOf(i), 0] = newGraphWithHM.GraphNodes[i].NodeVertex.X;
                    P_outer[outerNodeIndexList.IndexOf(i), 1] = newGraphWithHM.GraphNodes[i].NodeVertex.Y;
                }
            }
            #endregion
            #region 计算inner，outer相关矩阵
            Matrix inner_innerM;
            // 矩阵Q：inner_OuterAdjacencyMatrix
            Matrix inner_outerM;
            Matrix outer_innerM;
            Matrix outer_outerM;
            Matrix wholeM = GraphLoLToMatrix(triangleGraphTables,
                             innerNodeIndexList,
                             outerNodeIndexList,
                             out inner_innerM,
                             out inner_outerM,
                             out outer_innerM,
                             out outer_outerM);

            /**/
            //List<string> debugWhole = PrintMatrix(wholeM);
            //List<string> debugii = PrintMatrix(inner_innerM);
            //List<string> debugio = PrintMatrix(inner_outerM);
            //List<string> debugoo = PrintMatrix(outer_outerM);

            Matrix inner_InnerLaplacianM = LaplacianMatrix(inner_innerM, innerNodeDegreeList);

            /**/
            //List<string> printLaplacianM = PrintMatrix(inner_InnerLaplacianM);


            bool flag = inner_InnerLaplacianM.Invert(0.0);
            Matrix inverse_Inner_InnerLaplacianM = inner_InnerLaplacianM;

            #endregion

            #region 求解矩阵 P(inner)
            // 求解矩阵 P(inner)
            Matrix P_inner = new Matrix(innerNodeCount, 2);
            P_inner = inverse_Inner_InnerLaplacianM * inner_outerM * P_outer;
            // P_inner.Scale(-1.0);
            #endregion

            #region 新的NodePoint列表
            List<Point3d> newNodePoints = new List<Point3d>();
            for (int i = 0; i < newGraphNodes.Count; i++)
            {
                if (newGraphNodes[i].IsInner)
                {
                    newNodePoints.Add(new Point3d(P_inner[innerNodeIndexList.IndexOf(i), 0], P_inner[innerNodeIndexList.IndexOf(i), 1], 0.0));
                }
                else
                {
                    newNodePoints.Add(new Point3d(P_outer[outerNodeIndexList.IndexOf(i), 0], P_outer[outerNodeIndexList.IndexOf(i), 1], 0.0));
                }
            }

            for (int i = 0; i < newGraphNodes.Count; i++)
            {
                newGraphNodes[i].NodeVertex = newNodePoints[i];
            }
            #endregion

            #region 修改对应的PlanktonMesh

            List<List<int>> faceVertexOrder = new List<List<int>>();
            for (int i = 0; i < newPlanktonMesh.Faces.Count; i++)
            {
                faceVertexOrder.Add(new List<int>());
                faceVertexOrder[i].AddRange(newPlanktonMesh.Faces.GetFaceVertices(i));
            }

            // 构造修改后的PlanktonMesh
            PlanktonMesh embededPlanktonMesh = new PlanktonMesh();
            // Vertex位置信息改变
            for (int i = 0; i < newGraphNodes.Count; i++)
            {
                embededPlanktonMesh.Vertices.Add(newGraphNodes[i].NodeVertex.X, newGraphNodes[i].NodeVertex.Y, newGraphNodes[i].NodeVertex.Z);
            }
            // face信息没变
            embededPlanktonMesh.Faces.AddFaces(faceVertexOrder);

            #endregion

            // 输出新的GraphWithHM
            GraphWithHM embededGraphWithHM = new GraphWithHM(embededPlanktonMesh, newGraphNodes, originGraphTables);

            embededGraphWithHM.TreeNodeLabel = newGraphWithHM.TreeNodeLabel;
            embededGraphWithHM.TreeNodeHistory = newGraphWithHM.TreeNodeHistory;
            embededGraphWithHM.VolumeNodes = newGraphWithHM.VolumeNodes;
            embededGraphWithHM.VolumeNodesIndex = newGraphWithHM.VolumeNodesIndex;
            embededGraphWithHM.VolumeContainsWhichInnerNode = newGraphWithHM.VolumeContainsWhichInnerNode;

            return embededGraphWithHM;
        }

        public static Matrix GraphLoLToMatrix(List<List<int>> graphLoL,
                                     List<int> innerNodeIndexList,
                                     List<int> outerNodeIndexList,
                                     out Matrix inner_innerM,
                                     out Matrix inner_outerM,
                                     out Matrix outer_innerM,
                                     out Matrix outer_outerM)
        {
            Matrix wholeMatrix = new Matrix(graphLoL.Count, graphLoL.Count);

            List<List<int>> sortedGraphLoL = new List<List<int>>();

            List<int> sortedIndex = new List<int>();

            // 先outer再inner

            for (int i = 0; i < outerNodeIndexList.Count; i++)
            {
                sortedGraphLoL.Add(new List<int>(graphLoL[outerNodeIndexList[i]]));
                sortedIndex.Add(outerNodeIndexList[i]);
            }
            for (int i = 0; i < innerNodeIndexList.Count; i++)
            {
                sortedGraphLoL.Add(new List<int>(graphLoL[innerNodeIndexList[i]]));
                sortedIndex.Add(innerNodeIndexList[i]);
            }

            // 如果前面先outer再inner打乱顺序的话，就不能再两层for循环，然后i，j了，列的序号会错

            for (int i = 0; i < sortedGraphLoL.Count; i++)
            {
                for (int j = 0; j < sortedGraphLoL.Count; j++)
                {
                    if (sortedGraphLoL[i].Contains(sortedIndex[j]))
                    {
                        wholeMatrix[i, j] = 1;
                    }
                    else
                    {
                        wholeMatrix[i, j] = 0;
                    }
                }
            }

            outer_outerM = new Matrix(outerNodeIndexList.Count, outerNodeIndexList.Count);
            for (int i = 0; i < outerNodeIndexList.Count; i++)
            {
                for (int j = 0; j < outerNodeIndexList.Count; j++)
                {
                    outer_outerM[i, j] = wholeMatrix[i, j];
                }
            }

            outer_innerM = new Matrix(outerNodeIndexList.Count, innerNodeIndexList.Count);
            for (int i = 0; i < outerNodeIndexList.Count; i++)
            {
                for (int j = outerNodeIndexList.Count; j < graphLoL.Count; j++)
                {
                    outer_innerM[i, j - outerNodeIndexList.Count] = wholeMatrix[i, j];
                }
            }

            inner_outerM = new Matrix(innerNodeIndexList.Count, outerNodeIndexList.Count);
            for (int i = outerNodeIndexList.Count; i < graphLoL.Count; i++)
            {
                for (int j = 0; j < outerNodeIndexList.Count; j++)
                {
                    inner_outerM[i - outerNodeIndexList.Count, j] = wholeMatrix[i, j];
                }
            }

            inner_innerM = new Matrix(innerNodeIndexList.Count, innerNodeIndexList.Count);
            for (int i = outerNodeIndexList.Count; i < graphLoL.Count; i++)
            {
                for (int j = outerNodeIndexList.Count; j < graphLoL.Count; j++)
                {
                    inner_innerM[i - outerNodeIndexList.Count, j - outerNodeIndexList.Count] = wholeMatrix[i, j];
                }
            }

            return wholeMatrix;
        }

        /// <summary>
        /// 调和矩阵L Harmonic Matrix，又称拉普拉斯矩阵 Laplacian Matrix。是图的矩阵表示
        /// 在图论和计算机科学中，邻接矩阵A（英语：adjacency matrix）是一种方阵，用来表示有限图。它的每个元素代表各点之间是否有边相连。
        /// 在数学领域图论中，度数矩阵D是一个对角矩阵 ，其中包含的信息为的每一个顶点的度数，也就是说，每个顶点相邻的边数[1] 它可以和邻接矩阵一起使用以构造图的拉普拉斯算子矩阵。
        /// L = D - A
        /// </summary>
        /// <param name="M">邻接矩阵 Adjacency Matrix，对角线都是0，其他部分有数据0或1</param>
        /// <param name="degrees">度数矩阵 Degrees Matrix，因为是只有对角线有数据，其他部分都是零，所以可以用list存储</param>
        /// <returns></returns>
        public static Matrix LaplacianMatrix(Matrix M, List<int> degrees)
        {
            if (!M.IsSquare)
            {
                return null;
            }
            if (M.RowCount != degrees.Count)
            {
                return null;
            }

            Matrix LaplacianM = new Matrix(M.RowCount, M.ColumnCount);
            for (int i = 0; i < M.RowCount; i++)
            {
                for (int j = 0; j < M.ColumnCount; j++)
                {
                    if (i == j)
                    {
                        LaplacianM[i, j] = degrees[i] - 0;
                    }
                    else
                    {
                        LaplacianM[i, j] = 0 - M[i, j];
                    }
                }
            }

            return LaplacianM;
        }


        //public List<string> PrintMatrix(Matrix M)
        //{
        //    List<string> output = new List<string>();
        //    for (int i = 0; i < M.RowCount; i++)
        //    {
        //        string str = string.Format("{0}行:", i);
        //        for (int j = 0; j < M.ColumnCount; j++)
        //        {
        //            str += string.Format("{0},", M[i, j]);
        //        }
        //        output.Add(str);
        //    }

        //    return output;
        //}


    }
}