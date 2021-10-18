using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.GUI.Gradient;
using Rhino;
using Rhino.Geometry;
using Rhino.Collections;
using System;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph
{
    public class UtilityFunctions
    {

        /// <summary>
        /// 将GH_Structure<GH_Integer>类型的树形数据，转化为DataTree<int>类型的树形数据
        /// </summary>
        /// <param name="gh_Structure_GraphTree">GH_Structure<GH_Integer>类型的树形数据</param>
        /// <param name="dataTree">DataTree<int>类型的树形数据</param>
        public static void GH_StructureToDataTree_Int(GH_Structure<GH_Integer> gh_Structure_GraphTree, ref DataTree<int> dataTree)
        {
            foreach (GH_Path gh_Path in gh_Structure_GraphTree.Paths)
            {
                List<int> branchitems = new List<int>();
                foreach (GH_Integer element in gh_Structure_GraphTree.get_Branch(gh_Path))
                {
                    branchitems.Add(Convert.ToInt32(element.Value));
                }
                dataTree.AddRange(branchitems, gh_Path);
            }
        }

        public static void GraphToSubGraph(DataTree<int> graph, GlobalParameter globalParameter, out DataTree<int> connectivityGraph, out DataTree<int> adjacencyGraph)
        {
            connectivityGraph = new DataTree<int>();
            adjacencyGraph = new DataTree<int>();

            for (int i = 0; i < graph.BranchCount; i++)
            {
                for (int j = 0; j < graph.Branch(i).Count; j++)
                {
                    if (i < globalParameter.VolumeNodeCount)
                    {
                        connectivityGraph.EnsurePath(i);
                        if (graph.Branch(i)[j] < globalParameter.VolumeNodeCount)
                        {
                            connectivityGraph.Branch(i).Add(graph.Branch(i)[j]);
                        }
                    }
                    else
                    {
                        adjacencyGraph.EnsurePath(i - globalParameter.VolumeNodeCount);
                        adjacencyGraph.Branch(i - globalParameter.VolumeNodeCount).Add(graph.Branch(i)[j]);
                    }
                }
            }
        }

        /// <summary>
        /// 将Connectivity图结构中的Edge转化为可以显示的Line线段
        /// </summary>
        /// <param name="connectivityGraphTree">connectivity图结构</param>
        /// <param name="volumeNodeList">顶点列表</param>
        /// <returns>可以显示的Line线段</returns>
        public static List<Line> ConnectivityGraphEdgeToLine(DataTree<int> connectivityGraphTree, List<Point3d> volumeNodeList)
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
        /// 将Adjacency图结构中的Edge转化为可以显示的Line线段
        /// </summary>
        /// <param name="adjacencyGraphTree">adjacent图结构</param>
        /// <param name="volumeNodeList">volume顶点列表</param>
        /// <param name="boundaryNodeList">boundary顶点列表</param>
        /// <returns>可以显示的Line线段</returns>
        public static List<Line> AdjacencyGraphEdgeToLine(DataTree<int> adjacencyGraphTree, List<Point3d> volumeNodeList, List<Point3d> boundaryNodeList, GlobalParameter globalParameter)
        {
            List<Line> lineList = new List<Line>();

            for (int i = 0; i < adjacencyGraphTree.BranchCount; i++)
            {
                for (int j = 0; j < adjacencyGraphTree.Branch(i).Count; j++)
                {
                    if (adjacencyGraphTree.Branch(i)[j] == -1 - globalParameter.BoundaryNodeCount)
                    {
                        continue;
                    }
                    Line line = new Line(volumeNodeList[i], boundaryNodeList[adjacencyGraphTree.Branch(i)[j] + globalParameter.BoundaryNodeCount]);
                    lineList.Add(line);
                }
            }

            return lineList;
        }





        /// <summary>
        /// 将一个列表的点，由一个坐标系整体平移到另一个坐标系
        /// </summary>
        /// <param name="p"></param>
        /// <param name="location1"></param>
        /// <param name="location2"></param>
        /// <returns></returns>
        public static List<Point3d> Relocate(List<Point3d> p, Plane location1, Plane location2)
        {
            Point3d point3d = new Point3d(0.0, 0.0, 0.0);
            List<Point3d> list = new List<Point3d>();
            foreach (Point3d point3d2 in p)
            {
                point3d += point3d2;
            }
            point3d /= (double)p.Count;
            location1.Origin = point3d;
            foreach (Point3d item in p)
            {
                item.Transform(Transform.PlaneToPlane(location1, location2));
                list.Add(item);
            }
            return list;
        }

        /// <summary>
        /// 将点的列表按照从正北开始，逆时针排序
        /// </summary>
        /// <param name="vPoints"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        public static List<Point3d> SortPolyPoints(List<Point3d> vPoints, Point3d center)
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
                    bool flag = PointCompare(cloneList[j], cloneList[j + 1], center);
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
        private static bool PointCompare(Point3d a, Point3d b, Point3d center)
        {
            Vector3d vectorOA = new Vector3d(a) - new Vector3d(center);
            Vector3d vectorOB = new Vector3d(b) - new Vector3d(center);

            // OA,OB分别与0,1,0的夹角的弧度值
            double angleOA = Vector3d.VectorAngle(new Vector3d(0, 1, 0), vectorOA);
            double angleOB = Vector3d.VectorAngle(new Vector3d(0, 1, 0), vectorOB);

            // 向量0,1,0和向量OA的叉积
            Vector3d vectorZOA = Vector3d.CrossProduct(new Vector3d(0, 1, 0), vectorOA);
            if (vectorZOA.Z < 0)
            {
                angleOA = 2 * Math.PI - angleOA;
            }
            Vector3d vectorZOB = Vector3d.CrossProduct(new Vector3d(0, 1, 0), vectorOB);
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

        /// <summary>
        /// 计算体量占地面积的比例
        /// </summary>
        /// <param name="nodeAttributes"></param>
        public static void CalculateAreaProportion(List<NodeAttribute> nodeAttributes)
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



        /// <summary>
        /// 用来画虚线，详见Triangulate电池
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="pattern"></param>
        /// <returns></returns>
        public static IEnumerable<Curve> ApplyDashPattern(Curve curve, double[] pattern)
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





    }
}