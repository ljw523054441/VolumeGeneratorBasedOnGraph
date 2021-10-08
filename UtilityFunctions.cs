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

        // public static void 

        //public static void DataTreeMerge_Int(ref DataTree<int> mainDataTree, DataTree<int> secondaryDataTree)
        //{
        //    if (mainDataTree.BranchCount == secondaryDataTree.BranchCount)
        //    {
        //        foreach (GH_Path gh_Path in mainDataTree.Paths)
        //        {
        //            mainDataTree.AddRange(secondaryDataTree.Branch(gh_Path), gh_Path);
        //        }
        //    }
        //}

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
    }
}