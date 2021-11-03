using Grasshopper.Kernel;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Reflection;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDecomposeVertex : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDecomposeVertex class.
        /// </summary>
        public GhcDecomposeVertex()
          : base("DecomposeVertex", "DecomposeVertex",
              "拆解度大于4的顶点",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);
            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("NewTriangleHalfedgeMesh", "NTHMesh", "新生成的三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);

            pManager.AddGenericParameter("DebugFacesHalfedges", "DebugFH", "Debug结果面的半边", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int innerNodeCount = globalParameter.VolumeNodeCount;
            int outerNodeCount = globalParameter.BoundaryNodeCount;

            List<int> innerNodeIndexs = new List<int>();
            for (int i = 0; i < innerNodeCount; i++)
            {
                innerNodeIndexs.Add(i);
            }
            List<int> outerNodeIndexs = new List<int>();
            for (int i = 0; i < outerNodeCount; i++)
            {
                outerNodeIndexs.Add(i + innerNodeCount);
            }

            PlanktonMesh P = new PlanktonMesh();

            PlanktonMesh newPlanktonMesh = new PlanktonMesh();

            if (DA.GetData("TheChosenTriangleHalfedgeMesh", ref P))
            {
                // 利用PlanktonMesh为复制准备的构造函数，进行深拷贝
                PlanktonMesh PDeepCopy = new PlanktonMesh(P);
                

                // 对于每个InnerNode来说，判断它的度是否大于4
                for (int i = 0; i < PDeepCopy.Vertices.Count - outerNodeCount; i++)
                {
                    
                }

                int degree = PDeepCopy.Vertices.GetValence(1);

                // 点1周围的halfedgeindex列表，顺时针
                int[] halfEdgesIndexStartFromVertex = PDeepCopy.Vertices.GetHalfedges(1);

                switch (degree)
                {
                    case 5:
                        // SplitVertex(PDeepCopy, 1, innerNodeIndexs);
                        // P.Vertices.SetVertex(P.Vertices.Count-1, P.Vertices[0].ToXYZ().X, P.Vertices[0].ToXYZ().Y + (float)0.5, P.Vertices[0].ToXYZ().Z);

                        //Point3d startPoint = PDeepCopy.Vertices[1].ToPoint3d();
                        //Point3d endPoint = PDeepCopy.Vertices[PDeepCopy.Halfedges.EndVertex(15)].ToPoint3d();

                        //int newHalfEdgeIndex = PDeepCopy.Vertices.SplitVertex(15, 16);
                        // PDeepCopy.Vertices.SetVertex(PDeepCopy.Halfedges[15].StartVertex, (startPoint + endPoint) / 2);

                        newPlanktonMesh = SplitEdge(PDeepCopy, 15);

                        break;

                    default:
                        break;
                }

                DA.SetData("NewTriangleHalfedgeMesh", newPlanktonMesh);


                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(newPlanktonMesh);
                DA.SetDataList("DebugVerticesOutput", printVertices);

                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(newPlanktonMesh);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);

                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFaces(newPlanktonMesh);
                DA.SetDataList("DebugFacesOutput", printFaces);

                List<string> printFacesHalfedge = new List<string>();
                printFacesHalfedge = UtilityFunctions.PrintFacesHalfedges(newPlanktonMesh);
                DA.SetDataList("DebugFacesHalfedges", printFacesHalfedge);

            }
        }


        public PlanktonMesh SplitEdge(PlanktonMesh P, int halfedgeIndex)
        {
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);
            
            PlanktonVertexList vertices = PDeepCopy.Vertices;
            PlanktonHalfEdgeList halfedges = PDeepCopy.Halfedges;
            PlanktonFaceList faces = PDeepCopy.Faces;

            int adjacentFaceIndex = PDeepCopy.Halfedges[halfedgeIndex].AdjacentFace;
            int pairAdjacentFaceIndex = PDeepCopy.Halfedges[PDeepCopy.Halfedges.GetPairHalfedge(halfedgeIndex)].AdjacentFace;

            // 包围AdjacentFace的边
            List<int> hiAroundAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[adjacentFaceIndex].FirstHalfedge).ToList<int>();
            // 包围PairAdjacentFace的边
            List<int> hiAroundPairAdjacentFace = PDeepCopy.Halfedges.GetFaceCirculator(PDeepCopy.Faces[pairAdjacentFaceIndex].FirstHalfedge).ToList<int>();

            // 包围AdjacentFace的顶点
            List<int> viAroundAdjacentFace = PDeepCopy.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();
            // 包围PairAdjcentFace的顶点
            List<int> viAroundPairAdjacentFace = PDeepCopy.Faces.GetFaceVertices(pairAdjacentFaceIndex).ToList<int>();

            int startVertexIndex = PDeepCopy.Halfedges[halfedgeIndex].StartVertex;
            int endVertexIndex = PDeepCopy.Halfedges.EndVertex(halfedgeIndex);

            Point3d startVertex = PDeepCopy.Vertices[startVertexIndex].ToPoint3d();
            Point3d endVertex = PDeepCopy.Vertices[endVertexIndex].ToPoint3d();
            Point3d midVertex = (startVertex + endVertex) / 2;

            // 向P中增加顶点
            PDeepCopy.Vertices.Add(midVertex);
            int midVertexIndex = PDeepCopy.Vertices.Count - 1;

            // 构造添加顶点后的新Face
            for (int i = 0; i < viAroundAdjacentFace.Count; i++)
            {
                if (viAroundAdjacentFace[i] == startVertexIndex)
                {
                    viAroundAdjacentFace.Insert(i + 1, midVertexIndex);
                }
            }
            for (int i = 0; i < viAroundPairAdjacentFace.Count; i++)
            {
                if (viAroundPairAdjacentFace[i] == endVertexIndex)
                {
                    viAroundPairAdjacentFace.Insert(i + 1, midVertexIndex);
                }
            }

            // 建立需要进行置换的面的index列表和VertexIndex列表，方便转移P的Face属性
            List<int> needChangeFaceIndexs = new List<int>();
            needChangeFaceIndexs.Add(adjacentFaceIndex);
            needChangeFaceIndexs.Add(pairAdjacentFaceIndex);
            List<List<int>> needChangeViAround = new List<List<int>>();
            needChangeViAround.Add(viAroundAdjacentFace);
            needChangeViAround.Add(viAroundPairAdjacentFace);

            // 转移P的Vertex属性
            List<PlanktonXYZ> pPlanktonVertex = new List<PlanktonXYZ>();
            for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
            {
                pPlanktonVertex.Add(PDeepCopy.Vertices[i].ToXYZ());
            }

            // 转移P的Face属性，同时替换掉两个应该删除的面
            List<List<int>> pFaceVertexOrder = new List<List<int>>();

            for (int i = 0; i < PDeepCopy.Faces.Count; i++)
            {
                if (needChangeFaceIndexs.Contains(i))
                {
                    pFaceVertexOrder.Add(new List<int>());
                    pFaceVertexOrder[i].AddRange(needChangeViAround[needChangeFaceIndexs.IndexOf(i)]);
                }
                else
                {
                    pFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = PDeepCopy.Faces.GetFaceVertices(i);
                    pFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }


            // 用转移的P的Vertex属性和P的Face属性来构造新的PlanktonMeshnewP
            PlanktonMesh newP = new PlanktonMesh();

            newP.Vertices.AddVertices(pPlanktonVertex);
            newP.Faces.AddFaces(pFaceVertexOrder);

            return newP;
        }




        //public PlanktonMesh SplitVertex(PlanktonMesh P, int vertexIndexToSplit, List<int> innerNodeIndexs)
        //{
            
            
        //    int[] halfEdgesIndexStartFromVertex = P.Vertices.GetHalfedges(vertexIndexToSplit);


        //    /* 把for循环去掉直接测试看一下结果
        //     * i变为0
        //     */
        //    for (int i = 0; i < halfEdgesIndexStartFromVertex.Length; i++)
        //    {
        //        // int pairHalfEdgeIndex = P.Halfedges.GetPairHalfedge(halfEdgesIndexStartFromVertex[i]);

                
        //    }

        //    // 如果半边i的终点是innerNode，那么将这个半边split开始加点
        //    int endVertexIndex = P.Halfedges.EndVertex(halfEdgesIndexStartFromVertex[1]);
        //    if (innerNodeIndexs.Contains(endVertexIndex))
        //    {
        //        Point3d startPoint = P.Vertices[vertexIndexToSplit].ToPoint3d();
        //        Point3d endPoint = P.Vertices[P.Halfedges.EndVertex(halfEdgesIndexStartFromVertex[1])].ToPoint3d();
        //        int index = P.Halfedges.SplitEdge(halfEdgesIndexStartFromVertex[1]);
        //        // P.Vertices.SetVertex(P.Halfedges.EndVertex(index), (startPoint + endPoint) / 2);
        //        P.Vertices.SetVertex(P.Halfedges[index].StartVertex, (startPoint + endPoint) / 2);
        //    }

        //    PlanktonMesh newPlanktonMesh = new PlanktonMesh(P);

        //    return newPlanktonMesh;
        //}

        

        //public int SplitEdge(PlanktonMesh P, int index)
        //{
        //    int pair = P.Halfedges.GetPairHalfedge(index);

        //    // Create a copy of the existing vertex (user can move it afterwards if needs be)
        //    // 将中点作为新添加的顶点
        //    int start_vertex = P.Halfedges[index].StartVertex;
        //    int end_vertex = P.Halfedges[pair].StartVertex;
        //    int new_vertex_index = P.Vertices.Add((P.Vertices[start_vertex].ToPoint3d() + P.Vertices[end_vertex].ToPoint3d()) / 2);

        //    // Add a new halfedge pair
        //    int new_halfedge1 = P.Halfedges.AddPair(new_vertex_index, P.Halfedges.EndVertex(index), P.Halfedges[index].AdjacentFace);
        //    int new_halfedge2 = P.Halfedges.GetPairHalfedge(new_halfedge1);
        //    P.Halfedges[new_halfedge2].AdjacentFace = P.Halfedges[pair].AdjacentFace;
        //}

        //public int AddPair(PlanktonMesh P, int start, int end, int face)
        //{
        //    // he->next = he->pair
        //    int i = this.Count;
        //    P.Halfedges.Add(new PlanktonHalfedge(start, face, i + 1));
        //    P.Halfedges.Add(new PlanktonHalfedge(end, -1, i));
        //    return i;
        //}

        //public int SplitVertex(PlanktonMesh P, int first, int second)
        //{
        //    PlanktonHalfEdgeList hs = P.Halfedges;

        //    // 输入的两个半边是否都从同一个顶点开始，如果不是就返回-1
        //    int v_old = hs[first].StartVertex;
        //    if (v_old != hs[second].StartVertex)
        //    {
        //        return -1;
        //    }

        //    // 创建一个v_old的副本，并添加到顶点列表中
        //    int v_new = P.Vertices.Add(RhinoSupport.ToPoint3d(P.Vertices[v_old].ToXYZ()));

        //    // 从second到first
        //    // 
        //    bool reset_v_old = false;

        //    List<int> vertexCirculator = hs.GetVertexCirculator(second).ToList();

        //    // 与second的起始顶点相关的可枚举的halfedge索引。围绕顶点顺时针排列。返回的可枚举项将从指定的半边（即second）开始。
        //    foreach (int h in vertexCirculator)
        //    {
        //        if (h == first)
        //        {
        //            break;
        //        }

        //        hs[h].StartVertex = v_new;
        //        // 如果新顶点还没有传出的半边，并且他现在是裸露的，那么
        //        if (P.Vertices[v_new].OutgoingHalfedge == -1
        //            && hs[h].AdjacentFace == -1)
        //        {
        //            P.Vertices[v_new].OutgoingHalfedge = h;
        //        }
        //        // 
        //        if (h == P.Vertices[v_old].OutgoingHalfedge)
        //        {
        //            reset_v_old = true;
        //        }
        //    }

        //    // 如果没有裸露的半边，就用second
        //    if (P.Vertices[v_new].OutgoingHalfedge == -1)
        //    {
        //        P.Vertices[v_new].OutgoingHalfedge = second;
        //    }



        //}


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("8505b065-3ef0-4f90-bcb8-9cbf18874dca"); }
        }
    }
}