using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.GUI.Gradient;
using Rhino;
using Rhino.Geometry;
using Rhino.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class UtilityFunctions
    {

        #region 数据类型转化

        /// <summary>
        /// 将GH_Structure<GH_Integer>类型的树形数据，转化为DataTree<int>类型的树形数据
        /// </summary>
        /// <param name="gh_Structure_GraphTree">GH_Structure<GH_Integer>类型的树形数据</param>
        /// <param name="dataTree">DataTree<int>类型的树形数据</param>
        internal static void GH_StructureToDataTree_Int(GH_Structure<GH_Integer> gh_Structure_GraphTree, ref DataTree<int> dataTree)
        {
            foreach (GH_Path gh_Path in gh_Structure_GraphTree.Paths)
            {
                List<int> branchitems = new List<int>();
                foreach (GH_Integer element in gh_Structure_GraphTree.get_Branch(gh_Path))
                {
                    branchitems.Add(element.Value);
                }
                dataTree.AddRange(branchitems, gh_Path);
            }
        }

        internal static GH_Structure<GH_Integer> DataTreeToGH_Structure_Int(DataTree<int> dataTree)
        {
            GH_Structure<GH_Integer> gh_Structure_GraphTree = new GH_Structure<GH_Integer>();
            foreach (GH_Path ghPath in dataTree.Paths)
            {
                List<GH_Integer> branchitems = new List<GH_Integer>();
                foreach (int element in dataTree.Branch(ghPath))
                {
                    GH_Integer ghInteger = new GH_Integer(element);
                    branchitems.Add(ghInteger);
                }
                gh_Structure_GraphTree.AppendRange(branchitems);
            }
            return gh_Structure_GraphTree;
        }

        internal static DataTree<double> GH_StructureToDataTree_Double(GH_Structure<GH_Number> gh_Structure)
        {
            DataTree<double> dataTree = new DataTree<double>();
            foreach (GH_Path gh_Path in gh_Structure.Paths)
            {
                List<double> branchitems = new List<double>();
                foreach (GH_Number element in gh_Structure.get_Branch(gh_Path))
                {
                    branchitems.Add(element.Value);
                }
                dataTree.AddRange(branchitems, gh_Path);
            }

            return dataTree;
        }

        internal static List<List<double>> GH_StructureToLoL_Double(GH_Structure<GH_Number> gh_Structure)
        {
            List<List<double>> lol = new List<List<double>>();
            foreach (GH_Path gh_Path in gh_Structure.Paths)
            {
                lol.Add(new List<double>());
                List<double> branchitems = new List<double>();
                foreach (GH_Number element in gh_Structure.get_Branch(gh_Path))
                {
                    branchitems.Add(element.Value);
                }
                lol.Last().AddRange(branchitems);
            }
            return lol;
        }

        internal static List<List<int>> GH_StructureToLoL_Int(GH_Structure<GH_Integer> gh_Structure)
        {
            List<List<int>> lol = new List<List<int>>();
            foreach (GH_Path gh_Path in gh_Structure.Paths)
            {
                lol.Add(new List<int>());
                List<int> branchItems = new List<int>();
                foreach (GH_Integer element in gh_Structure.get_Branch(gh_Path))
                {
                    branchItems.Add(element.Value);
                }
                lol.Last().AddRange(branchItems);
            }
            return lol;
        }

        /// <summary>
        /// 由LoL转化成DataTree的泛型方法
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="LoL"></param>
        /// <returns></returns>
        internal static DataTree<T> LoLToDataTree<T>(List<List<T>> LoL)
        {
            DataTree<T> dataTree = new DataTree<T>();

            for (int i = 0; i < LoL.Count; i++)
            {
                dataTree.EnsurePath(i);
                dataTree.Branch(i).AddRange(LoL[i]);
            }

            return dataTree;
        }

        /// <summary>
        /// 由DataTree转化为LoL的泛型方法
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="dataTree"></param>
        /// <returns></returns>
        internal static List<List<T>> DataTreeToLoL<T>(DataTree<T> dataTree)
        {
            List<List<T>> LoL = new List<List<T>>();

            for (int i = 0; i < dataTree.BranchCount; i++)
            {
                LoL.Add(new List<T>());
                LoL[i].AddRange(dataTree.Branch(i));
            }

            return LoL;
        }

        internal static List<T> Shift<T>(List<T> originList, int shift, bool isLeft)
        {
            if (isLeft)
            {
                // 深拷贝
                List<T> newList = new List<T>();
                newList.AddRange(originList);

                int iter = 0;
                while (iter < shift)
                {
                    T item = newList[0];
                    newList.RemoveAt(0);
                    newList.Add(item);
                    iter++;
                }

                return newList;
            }
            else
            {
                // 深拷贝
                List<T> newList = new List<T>();
                newList.AddRange(originList);

                int iter = 0;
                while (iter < shift)
                {
                    T item = newList[newList.Count - 1];
                    newList.Insert(0, item);
                    newList.RemoveAt(newList.Count - 1);
                    iter++;
                }

                return newList;
            }
        }

        /// <summary>
        /// 二维数组转置函数
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array"></param>
        /// <returns></returns>
        internal static T[,] Transpose<T>(T[,] array)
        {
            int x = array.GetUpperBound(0);// 一维
            int y = array.GetUpperBound(1);// 二维
            T[,] newArray = new T[y + 1, x + 1];// 构造转置二维数组
            for (int i = 0; i <= x; i++)
            {
                for (int j = 0; j <= y; j++)
                {
                    newArray[j, i] = array[i, j];
                }
            }

            return newArray;
        }

        /// <summary>
        /// 将二维列表(List)转换成二维数组，二维数组转置，然后将二维数组转换成列表
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="original"></param>
        /// <returns></returns>
        internal static List<List<T>> Transpose<T>(List<List<T>> original)
        {
            List<T>[] array = original.ToArray();
            List<List<T>> lists = new List<List<T>>();
            if (array.Length == 0)
            {
                throw new IndexOutOfRangeException("Index Out Of Range");
            }
            int x = array[0].Count;
            int y = original.Count;

            // 将列表转换成数组
            T[,] twoArray = new T[y, x];
            for (int i = 0; i < y; i++)
            {
                int j = 0;
                foreach (T item in array[i])
                {
                    twoArray[i, j] = item;
                    j++;
                }
            }

            T[,] newTwoArray = new T[x, y];
            newTwoArray = Transpose<T>(twoArray);// 转置

            // 二维数组转换成二维List集合
            for (int i = 0; i < x; i++)
            {
                List<T> list = new List<T>();
                for (int j = 0; j < y; j++)
                {
                    list.Add(newTwoArray[i, j]);
                }
                lists.Add(list);
            }

            return lists;
        }

        internal static DataTree<T> FlipMatrix<T>(DataTree<T> source)
        {
            DataTree<T> target = new DataTree<T>();
            int rowCount = source.BranchCount;
            int columnCount = 0;
            for (int i = 0; i < source.BranchCount; i++)
            {
                int num = source.Branch(i).Count;
                if (num > columnCount)
                {
                    columnCount = num;
                }
            }
            for (int i = 0; i < columnCount; i++)
            {
                target.EnsurePath(i);
                for (int j = 0; j < rowCount; j++)
                {
                    if (source.ItemExists(new GH_Path(j),i))
                    {
                        target.Branch(i).Add(source.Branch(j)[i]);
                    }
                    else
                    {
                        target.Branch(i).Add(default(T));
                    }
                }
            }

            return target;
        }

        #endregion

        #region 图形绘制

        /// <summary>
        /// 将Connectivity图结构中的Edge转化为可以显示的Line线段
        /// </summary>
        /// <param name="connectivityGraphTree">connectivity图结构</param>
        /// <param name="volumeNodeList">顶点列表</param>
        /// <returns>可以显示的Line线段</returns>
        internal static List<Line> ConnectivityGraphEdgeToLine(DataTree<int> connectivityGraphTree, List<Point3d> volumeNodeList)
        {
            List<Line> lineList = new List<Line>();

            for (int i = 0; i < connectivityGraphTree.BranchCount; i++)
            {
                for (int j = 0; j < connectivityGraphTree.Branch(i).Count; j++)
                {
                    Line line = new Line(volumeNodeList[i], volumeNodeList[connectivityGraphTree.Branch(i)[j]]);
                    lineList.Add(line);
                }
            }

            return lineList;
        }

        /// <summary>
        /// 绘制表示图结构关系的Line
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="graphVertices"></param>
        /// <returns></returns>
        public static List<Line> GraphEdgeLine(List<List<int>> graph, List<GraphNode> graphNodes)
        {
            List<Point3d> graphVertices = new List<Point3d>();
            for (int i = 0; i < graphNodes.Count; i++)
            {
                graphVertices.Add(graphNodes[i].NodeVertex);
            }

            List<Line> list = new List<Line>();
            for (int i = 0; i < graph.Count; i++)
            {
                for (int j = 0; j < graph[i].Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graph[i][j]]);
                    list.Add(item);
                }
            }
            return list;
        }

        /// <summary>
        /// 用来画虚线，详见Triangulate电池
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="pattern"></param>
        /// <returns></returns>
        internal static IEnumerable<Curve> ApplyDashPattern(Curve curve, double[] pattern)
        {
            if (pattern == null || pattern.Length == 0)
            {
                return new Curve[] { curve };
            }

            double curveLength = curve.GetLength();
            List<Curve> dashes = new List<Curve>();

            double offset0 = 0.0;
            int index = 0;
            while (true)
            {
                double dashLength = pattern[index++];
                if (index >= pattern.Length)
                    index = 0;

                // Compute the offset of the current dash from the curve start.
                double offset1 = offset0 + dashLength;
                if (offset1 > curveLength)
                    offset1 = curveLength;

                // Solve the curve parameters at the current dash start and end.
                double t0, t1;
                curve.LengthParameter(offset0, out t0);
                curve.LengthParameter(offset1, out t1);

                Curve dash = curve.Trim(t0, t1);
                if (dash != null)
                    dashes.Add(dash);

                // Get the current gap length.
                double gapLength = pattern[index++];
                if (index >= pattern.Length)
                    index = 0;

                // Set the start of the next dash to be the end of the current
                // dash + the length of the adjacent gap.
                offset0 = offset1 + gapLength;

                // Abort when we've reached the end of the curve.
                if (offset0 >= curveLength)
                    break;
            }

            return dashes;
        }

        #endregion

        #region Point相关计算

        /// <summary>
        /// 将一个列表的点，由一个坐标系整体平移到另一个坐标系
        /// </summary>
        /// <param name="p"></param>
        /// <param name="location1"></param>
        /// <param name="location2"></param>
        /// <returns></returns>
        internal static List<Point3d> Relocate(List<Point3d> p, Plane location1, Plane location2)
        {
            //Point3d point3d = new Point3d(0.0, 0.0, 0.0);
            List<Point3d> list = new List<Point3d>();
            //foreach (Point3d point3d2 in p)
            //{
            //    point3d += point3d2;
            //}
            //point3d /= (double)p.Count;
            //location1.Origin = point3d;
            foreach (Point3d item in p)
            {
                item.Transform(Transform.PlaneToPlane(location1, location2));
                list.Add(item);
            }
            return list;
        }

        /// <summary>
        /// 计算一堆点的中点
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        internal static Point3d CalCenterPoint(List<Point3d> p)
        {
            Point3d centerPoint = new Point3d();

            for (int i = 0; i < p.Count; i++)
            {
                centerPoint += p[i];
            }

            centerPoint /= p.Count;

            return centerPoint;
        }

        /// <summary>
        /// 将点的列表按照从正北开始，逆时针排序
        /// </summary>
        /// <param name="vPoints"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        internal static List<Point3d> SortPolyPoints(List<Point3d> vPoints, Point3d center)
        {
            List<Point3d> cloneList = new List<Point3d>();
            cloneList.AddRange(vPoints);


            if (cloneList == null || cloneList.Count == 0)
            {
                return null;
            }

            for (int i = 0; i < cloneList.Count - 1; i++)
            {
                for (int j = 0; j < cloneList.Count - i - 1; j++)
                {
                    bool flag = PointCompare(cloneList[j], cloneList[j + 1], center, cloneList[0]);
                    if (flag)
                    {
                        Point3d tmp = cloneList[j];
                        cloneList[j] = cloneList[j + 1];
                        cloneList[j + 1] = tmp;
                    }
                }
            }
            return cloneList;
        }

        /// <summary>
        /// 若点a大于点b，即点a在点b的顺时针方向，返回true，否则返回false
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        public static bool PointCompare(Point3d a, Point3d b, Point3d center, Point3d firstInputPoint)
        {
            Vector3d vectorOA = new Vector3d(a) - new Vector3d(center);
            Vector3d vectorOB = new Vector3d(b) - new Vector3d(center);
            Vector3d vectorOC = new Vector3d(firstInputPoint) - new Vector3d(center);

            // OA,OB分别与第一个输入的outerNodePoint的夹角的弧度值
            double angleOA = Vector3d.VectorAngle(vectorOC, vectorOA);
            double angleOB = Vector3d.VectorAngle(vectorOC, vectorOB);

            // 向量OC和向量OA的叉积
            Vector3d vectorZOA = Vector3d.CrossProduct(vectorOC, vectorOA);
            if (vectorZOA.Z < 0)
            {
                angleOA = 2 * Math.PI - angleOA;
            }
            Vector3d vectorZOB = Vector3d.CrossProduct(vectorOC, vectorOB);
            if (vectorZOB.Z < 0)
            {
                angleOB = 2 * Math.PI - angleOB;
            }


            if (angleOA < angleOB)
            {
                return false;
            }
            return true;

        }

        #endregion

        #region 图结构相关操作


        #endregion

        /// <summary>
        /// 计算体量占地面积的比例
        /// </summary>
        /// <param name="nodeAttributes"></param>
        internal static void CalculateAreaProportion(List<GraphNodeAttribute> nodeAttributes)
        {
            List<double> nodeAreaProportions = new List<double>();
            for (int i = 0; i < nodeAttributes.Count; i++)
            {
                nodeAreaProportions.Add(nodeAttributes[i].NodeArea);
            }

            for (int i = 0; i < nodeAttributes.Count; i++)
            {
                nodeAttributes[i].NodeAreaProportion =  nodeAttributes[i].NodeArea / nodeAreaProportions.Sum();
            }
        }

        #region 半边结构DebugPrint
        /// <summary>
        /// HalfedgeMesh中每个顶点的属性
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        internal static List<string> PrintVertices(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                string str = "Vertices[" + i.ToString() + "]=";
                str += mesh.Vertices[i].X.ToString() +
                    "," + mesh.Vertices[i].Y.ToString() +
                    "," + mesh.Vertices[i].Z.ToString();
                output.Add(str);
            }
            return output;
        }

        internal static List<string> PrintVertexConnection(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                string str = "Vertices[" + i.ToString() + "]=";
                for (int j = 0; j < mesh.Vertices.GetVertexNeighbours(i).Length; j++)
                {
                    str += mesh.Vertices.GetVertexNeighbours(i)[j].ToString() + ",";
                }
                output.Add(str);
            }
            return output;
        } 

        /// <summary>
        /// HalfedgeMesh中每条半边的属性
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        internal static List<string> PrintHalfedges(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            output.Add("Format: StartVertex,AdjacentFace,NextHalfedge,PrevHalfedge");
            for (int i = 0; i < mesh.Halfedges.Count; i++)
            {
                string str = "Halfedges[" + i.ToString() + "]=";
                str += mesh.Halfedges[i].StartVertex.ToString() + "," +
                     mesh.Halfedges[i].AdjacentFace.ToString() + "," +
                      mesh.Halfedges[i].NextHalfedge.ToString() + "," +
                       mesh.Halfedges[i].PrevHalfedge.ToString();
                output.Add(str);
            }
            return output;
        }

        internal static List<string> PrintHalfedgeStartAndEnd(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < mesh.Halfedges.Count; i++)
            {
                string str = "Halfedges[" + i.ToString() + "]=";
                str += mesh.Halfedges[i].StartVertex.ToString() + "," +
                       mesh.Halfedges.EndVertex(i).ToString();
                output.Add(str);
            }
            return output;
        }

        /// <summary>
        /// HalfedgeMesh的面由哪些顶点组成
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        internal static List<string> PrintFacesVertices(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                string str = "Faces[" + i.ToString() + "]=";
                int[] findex = mesh.Faces.GetFaceVertices(i);
                for (int j = 0; j < findex.Length; j++)
                {
                    if (j > 0) str += ",";
                    str += findex[j].ToString();
                }
                output.Add(str);
            }
            return output;
        }

        /// <summary>
        /// HalfedgeMesh的面由哪些半边组成
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        internal static List<string> PrintFacesHalfedges(PlanktonMesh mesh)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                string str = "Faces[" + i.ToString() + "]=";
                int[] fhindex = mesh.Faces.GetHalfedges(i);
                for (int j = 0; j < fhindex.Length; j++)
                {
                    if (j > 0) str += ",";
                    str += fhindex[j].ToString();
                }
                output.Add(str);
            }
            return output;
        }
        #endregion

        #region 半边结构相关操作

        /// <summary>
        /// 获取半边数据结构中，每个面所邻接的面的序号
        /// </summary>
        /// <param name="D"></param>
        /// <returns></returns>
        internal static List<List<int>> GetAdjacencyFaceIndexs(PlanktonMesh D)
        {
            List<List<int>> faceAdjacency = new List<List<int>>();
            for (int i = 0; i < D.Faces.Count; i++)
            {
                faceAdjacency.Add(new List<int>());

                HashSet<int> adjacentFace = new HashSet<int>();
                int[] halfedges = D.Faces.GetHalfedges(i);
                // length = halfedges.Length;
                for (int j = 0; j < halfedges.Length; j++)
                {
                    // 找D.Halfedges[halfedges[j]]的对边，然后再.AdjacentFace
                    int index = D.Halfedges[D.Halfedges.GetPairHalfedge(halfedges[j])].AdjacentFace;
                    if (index < 0)
                    {
                        continue;
                    }
                    else
                    {
                        adjacentFace.Add(index);
                    }
                }

                List<int> list = adjacentFace.ToList();
                faceAdjacency[i].AddRange(list);
            }
            return faceAdjacency;
        }

        /// <summary>
        /// 拉普拉斯网格平滑
        /// </summary>
        /// <param name="P"></param>
        /// <param name="W"></param>
        /// <param name="Strength"></param>
        /// <returns></returns>
        internal static Vector3d[] LaplacianSmooth(PlanktonMesh P, int W, double Strength)
        {
            int vertCount = P.Vertices.Count;
            Vector3d[] smooth = new Vector3d[vertCount];

            for (int i = 0; i < vertCount; i++)
            {
                if (P.Vertices[i].IsUnused == false
                    && P.Vertices.IsBoundary(i) == false)
                {
                    int[] neighbours = P.Vertices.GetVertexNeighbours(i);
                    Point3d vertex = P.Vertices[i].ToPoint3d();
                    Point3d centroid = new Point3d();
                    if (W == 0)
                    {
                        for (int j = 0; j < neighbours.Length; j++)
                        {
                            centroid = centroid + P.Vertices[neighbours[j]].ToPoint3d();
                        }
                        smooth[i] = ((centroid * (1.0 / P.Vertices.GetValence(i))) - vertex) * Strength;
                    }
                    if (W == 1)
                    {
                        //get the radial vectors of the 1-ring
                        //get the vectors around the 1-ring
                        //get the cotangent weights for each edge

                        int valence = neighbours.Length;

                        Point3d[] neighbourPts = new Point3d[valence];
                        Vector3d[] radial = new Vector3d[valence];
                        Vector3d[] around = new Vector3d[valence];
                        double[] cotWeight = new double[valence];
                        double weightSum = 0;

                        for (int j = 0; j < valence; j++)
                        {
                            neighbourPts[j] = P.Vertices[neighbours[j]].ToPoint3d();
                            radial[j] = neighbourPts[j] - vertex;
                        }

                        for (int j = 0; j < valence; j++)
                        {
                            around[j] = neighbourPts[(j + 1) % valence] - neighbourPts[j];
                        }

                        for (int j = 0; j < neighbours.Length; j++)
                        {
                            //get the cotangent weights
                            int previous = (j + valence - 1) % valence;
                            Vector3d cross1 = Vector3d.CrossProduct(radial[previous], around[previous]);
                            double cross1Length = cross1.Length;
                            double dot1 = radial[previous] * around[previous];

                            int next = (j + 1) % valence;
                            Vector3d cross2 = Vector3d.CrossProduct(radial[next], around[j]);
                            double cross2Length = cross2.Length;
                            double dot2 = radial[next] * around[j];

                            cotWeight[j] = Math.Abs(dot1 / cross1Length) + Math.Abs(dot2 / cross2Length);
                            weightSum += cotWeight[j];
                        }

                        double invWeightSum = 1.0 / weightSum;
                        Vector3d thisSmooth = new Vector3d();
                        for (int j = 0; j < neighbours.Length; j++)
                        {
                            thisSmooth = thisSmooth + radial[j] * cotWeight[j];
                        }

                        smooth[i] = thisSmooth * invWeightSum * Strength;
                    }
                }
            }
            return smooth;
        }

        #endregion

        #region Segment操作

        /// <summary>
        /// 将给定超想的segment按照对应朝向的规则进行排序
        /// </summary>
        /// <param name="sameOrientationSegments"></param>
        /// <param name="orientation"></param>
        /// <returns></returns>
        internal static List<Line> SortSameOrientationSegments(List<Line> sameOrientationSegments, string orientation)
        {
            List<Line> cloneList = new List<Line>();
            cloneList.AddRange(sameOrientationSegments);

            if (cloneList == null || cloneList.Count == 0)
            {
                return null;
            }

            for (int i = 0; i < cloneList.Count - 1; i++)
            {
                for (int j = 0; j < cloneList.Count - i - 1; j++)
                {
                    bool flag = SegmentCenterCompare(cloneList[j], cloneList[j + 1], orientation);
                    if (flag)
                    {
                        Line tml = cloneList[j];
                        cloneList[j] = cloneList[j + 1];
                        cloneList[j + 1] = tml;
                    }
                }
            }
            return cloneList;
        }

        /// <summary>
        /// NEWS每个方向都是逆时针排序;
        /// 对于N方向，如果i在i+1左边，那么要交换，把i+1放在排序好的列表前面;
        /// 对于W方向，如果i在i+1下面，那么要交换，把i+1放在排序好的列表前面;
        /// 对于S方向，如果i在i+1左边，那么不用交换，把i放在排序好的列表前面;
        /// 对于E方向，如果i在i+1下面，那么不用交换，把i放在排序好的列表前面;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="orientation"></param>
        /// <returns></returns>
        private static bool SegmentCenterCompare(Line a, Line b, string orientation)
        {
            Point3d aCenter = a.PointAt(0.5);
            Point3d bCenter = b.PointAt(0.5);

            bool flag = false;

            switch (orientation)
            {
                // NEWS每个方向都是逆时针排序
                case "N":
                    // 对于N方向，如果i在i+1左边，那么要交换，把i+1放在排序好的列表前面
                    if (aCenter.X < bCenter.X)
                    {
                        flag = true;
                    }
                    else
                    {
                        flag = false;
                    }
                    break;
                case "W":
                    // 对于W方向，如果i在i+1下面，那么要交换，把i+1放在排序好的列表前面
                    if (aCenter.Y < bCenter.Y)
                    {
                        flag = true;
                    }
                    else
                    {
                        flag = false;
                    }
                    break;
                case "S":
                    // 对于S方向，如果i在i+1左边，那么不用交换，把i放在排序好的列表前面
                    if (aCenter.X < bCenter.X)
                    {
                        flag = false;
                    }
                    else
                    {
                        flag = true;
                    }
                    break;
                case "E":
                    // 对于E方向，如果i在i+1下面，那么不用交换，把i放在排序好的列表前面
                    if (aCenter.Y < bCenter.Y)
                    {
                        flag = false;
                    }
                    else
                    {
                        flag = true;
                    }
                    break;

                default:
                    break;
            }

            return flag;
        }

        #endregion

        #region 图结构DebugPrint

        public static List<string> PrintGraphNode(Graph graph)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < graph.GraphNodes.Count; i++)
            {
                string str = string.Format("Node[{0}]:{1}--{2}", i, graph.GraphNodes[i].NodeAttribute.NodeLabel, graph.GraphNodes[i].IsInner);
                output.Add(str);
            }
            return output;
        }

        public static List<string> PrintGraphLoL(Graph graph)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < graph.GraphTables.Count; i++)
            {
                string str = string.Format("Node[{0}]:", i);
                for (int j = 0; j < graph.GraphTables[i].Count; j++)
                {
                    if (j == graph.GraphTables[i].Count - 1)
                    {
                        str += string.Format("{0}.", graph.GraphTables[i][j]);
                    }
                    else
                    {
                        str += string.Format("{0},", graph.GraphTables[i][j]);
                    }
                    
                }
                output.Add(str);
            }
            return output;
        }

        public static List<string> PrintGraphVertex(Graph graph)
        {
            List<string> output = new List<string>();
            for (int i = 0; i < graph.GraphNodes.Count; i++)
            {
                string str = string.Format("Node[{0}]:", i);
                str += graph.GraphNodes[i].NodeVertex.ToString();
                output.Add(str);
            }
            return output;
        }

        #endregion

    }
}