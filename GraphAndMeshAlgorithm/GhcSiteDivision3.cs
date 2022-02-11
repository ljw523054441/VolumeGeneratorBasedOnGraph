using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Attributes;
using Plankton;
using PlanktonGh;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;
using Grasshopper.GUI.Canvas;
using Grasshopper.GUI;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcSiteDivision3 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcSiteDivision3 class.
        /// </summary>
        public GhcSiteDivision3()
          : base("GhcSiteDivision3", "SiteDivision2",
              "进行场地划分",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            //DualVertices = new List<Point3d>();
            DualPolylines = new List<List<Point3d>>();
            DualVertexTextDots = new List<TextDot>();

            //GraphNodePoints = new List<Point3d>();
            GraphEdges = new List<Line>();
            InnerNodeTextDots = new List<TextDot>();
            OuterNodeTextDots = new List<TextDot>();

            MeshForVisualize = null;
            DottedCurve = new List<Curve>();
        }

        private int Thickness;

        //private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        //private List<Point3d> GraphNodePoints;
        private List<Line> GraphEdges;
        private List<TextDot> InnerNodeTextDots;
        private List<TextDot> OuterNodeTextDots;

        private Mesh MeshForVisualize;
        private List<Curve> DottedCurve;

        public enum ShowMode { ShowTopology, NotShowTopology };/* 定义一个enum类型 */
        public ShowMode CompWorkMode { get; set; } = ShowMode.ShowTopology;/* 使用这个enum类型来定义一个代表电池工作状态的变量 */

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);

            pManager.AddGenericParameter("VerticesIndexForEachBS", "VIFBS", "每个BS上点对应的对偶图中的index", GH_ParamAccess.list);

            pManager.AddGenericParameter("IndexOnEachBS", "IOBS", "每个BS上的VolumeJunctionIndex", GH_ParamAccess.list);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectVolumeJunctionsIndex", "N_VJI", "需要构成分界线的两个volumeJunction的Index", GH_ParamAccess.list);
            pManager.AddGenericParameter("NeedToConnectBSIndex", "N_BSI", "这两个volumeJunction所对应的bsIndex", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectPairCorrespondingFaceIndex", "N_CFI", "需要构成分界线的两个volumeJunction所对对应的FaceIndex", GH_ParamAccess.list);

            pManager.AddIntegerParameter("BoundaryTCounts", "BoundaryTCounts", "在BS上的点所需要的t值", GH_ParamAccess.list);
            pManager.AddIntegerParameter("InnerTMaxCounts", "InnerTMaxCount", "由BS上的点生成的转折所需要的t值的最大数量", GH_ParamAccess.list);

            pManager.AddIntegerParameter("TotalLayerCount", "TLC", "对于D进行拆解的总层数", GH_ParamAccess.item);

            pManager.AddNumberParameter("BoundaryTDT", "BoundaryTDT", "每个边界分裂点对应的t值", GH_ParamAccess.tree);
            pManager.AddNumberParameter("InnerTMaxDT", "InnerTMaxDT", "", GH_ParamAccess.tree);
            //pManager[9].Optional = true;
            //pManager[10].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("tXDT", "tXDT", "", GH_ParamAccess.tree);
            pManager.AddNumberParameter("tYDT", "tYDT", "", GH_ParamAccess.tree);

            pManager.AddGenericParameter("DualGraphWithHM", "DualGraphWithHM", "", GH_ParamAccess.item);

            pManager.AddGenericParameter("DualPlanktonMesh", "DualPlanktonMesh", "", GH_ParamAccess.item);
            
            //pManager.AddCurveParameter("DeBug Polyline", "DP", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Thickness = 2;
            
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            
            List<BoundarySegment> sortedBoundarySegments = new List<BoundarySegment>();

            List<string> volumeJunctionsTexts = new List<string>();

            //List<List<List<int>>> allLayerVerticesIndexForEachBSLoL = new List<List<List<int>>>();
            //List<List<List<int>>> allLayerVolumeJunctionsIndexOnEachBSLoL = new List<List<List<int>>>();
            List<DataTree<int>> allLayerVerticesIndexForEachBS = new List<DataTree<int>>();
            List<DataTree<int>> allLayerVolumeJunctionsIndexOnEachBS = new List<DataTree<int>>();

            List<List<int[]>> allLayerVertexIndexPairs = new List<List<int[]>>();
            List<List<int[]>> allLayerBsIndexPairs = new List<List<int[]>>();

            List<List<int>> allLayerPairCorrespondingFaceIndexLoL_OnD = new List<List<int>>();

            List<int> allLayerBoundaryTCounts = new List<int>();
            List<int> allLayerInnerTMaxCounts = new List<int>();

            GH_Structure<GH_Number> allLayerBoundaryT_GHS = new GH_Structure<GH_Number>();
            GH_Structure<GH_Number> allLayerInnerTMax_GHS = new GH_Structure<GH_Number>();
            //DataTree<double> allLayerBoundaryTDT = new DataTree<double>();
            //DataTree<double> allLayerInnerTMaxDT = new DataTree<double>();

            List<List<double>> allLayerBoundaryTLoL = new List<List<double>>();
            List<List<double>> allLayerInnerTMaxLoL = new List<List<double>>();

            int totalLayerCount = 0;

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetDataList("SortedBoundarySegments", sortedBoundarySegments)
                //&& DA.GetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions)
                && DA.GetDataList("VerticesIndexForEachBS", allLayerVerticesIndexForEachBS)
                && DA.GetDataList("IndexOnEachBS", allLayerVolumeJunctionsIndexOnEachBS)

                && DA.GetDataList("VolumeJunctionsTexts", volumeJunctionsTexts)

                && DA.GetDataList("NeedToConnectVolumeJunctionsIndex", allLayerVertexIndexPairs)
                && DA.GetDataList("NeedToConnectBSIndex", allLayerBsIndexPairs)

                && DA.GetDataList("NeedToConnectPairCorrespondingFaceIndex", allLayerPairCorrespondingFaceIndexLoL_OnD)

                && DA.GetDataList("BoundaryTCounts", allLayerBoundaryTCounts)
                && DA.GetData("TotalLayerCount", ref totalLayerCount))
            {
                DA.GetDataTree<GH_Number>("BoundaryTDT", out allLayerBoundaryT_GHS);
                DA.GetDataTree<GH_Number>("InnerTMaxDT", out allLayerInnerTMax_GHS);
                allLayerBoundaryTLoL = UtilityFunctions.GH_StructureToLoL_Double(allLayerBoundaryT_GHS);
                allLayerInnerTMaxLoL = UtilityFunctions.GH_StructureToLoL_Double(allLayerInnerTMax_GHS);

                List<BoundarySegment> sortedBoundarySegmentsDeepCopy = new List<BoundarySegment>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    sortedBoundarySegmentsDeepCopy.Add(new BoundarySegment(sortedBoundarySegments[i]));
                }

                #region 计算boundingBox
                List<Point3d> cornerPoints = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegmentsDeepCopy.Count; i++)
                {
                    cornerPoints.Add(sortedBoundarySegmentsDeepCopy[i].From);
                }

                BoundingBox boundingBox = new BoundingBox(cornerPoints);
                #endregion

                #region 最外层
                int currentLayer = 0;
                #region 构造DataTree<int> verticesIndexForEachBSDT
                DataTree<int> verticesIndexForEachBSDT = allLayerVerticesIndexForEachBS[0];
                #endregion
                #region 构造indexOnEachBS
                DataTree<int> volumeJunctionsIndexOnEachBS = allLayerVolumeJunctionsIndexOnEachBS[0];
                #endregion
                #region 构造DataTree<double> tOnEachBS
                // 每个包含volumeJunction的bs为了确定volumeJunction位置所需要的t的datatree
                if (volumeJunctionsIndexOnEachBS.DataCount > allLayerBoundaryTLoL[0].Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("boundaryTList[{0}]的数量不足，请增加boundaryT的数量", currentLayer));
                }

                DataTree<double> boundaryTOnEachBS = new DataTree<double>();
                int index = 0;
                for (int i = 0; i < volumeJunctionsIndexOnEachBS.BranchCount; i++)
                {
                    boundaryTOnEachBS.EnsurePath(i);
                    for (int j = 0; j < volumeJunctionsIndexOnEachBS.Branch(i).Count; j++)
                    {
                        boundaryTOnEachBS.Branch(i).Add(allLayerBoundaryTLoL[0][index]);
                        index++;
                    }
                }


                for (int i = 0; i < boundaryTOnEachBS.BranchCount; i++)
                {
                    List<double> ts = boundaryTOnEachBS.Branch(i);
                    for (int j = 0; j < ts.Count - 1; j++)
                    {
                        if (ts[j] > ts[j + 1])
                        {
                            // 后一个t不能大于前一个t
                            return;
                        }
                    }
                }
                #endregion

                #region DataTree<Point3d> pointOnEachBS
                DataTree<Point3d> pointOnEachBS = new DataTree<Point3d>();
                for (int i = 0; i < boundaryTOnEachBS.BranchCount; i++)
                {
                    pointOnEachBS.EnsurePath(i);
                    for (int j = 0; j < boundaryTOnEachBS.Branch(i).Count; j++)
                    {
                        double t = boundaryTOnEachBS.Branch(i)[j];
                        if (sortedBoundarySegmentsDeepCopy[i].TurningTs.Count != 0)
                        {
                            int lineIndex = 0;
                            for (int k = 0; k < sortedBoundarySegmentsDeepCopy[i].TurningTs.Count - 1; k++)
                            {
                                if (t > sortedBoundarySegmentsDeepCopy[i].TurningTs[k] && t < sortedBoundarySegmentsDeepCopy[i].TurningTs[k + 1])
                                {
                                    // 注意这里一定要添加PointOnWhichSegments信息，即这个点是在BS的哪一段上
                                    sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments.Add(k);
                                    break;
                                }
                                else if (t > sortedBoundarySegmentsDeepCopy[i].TurningTs.Last() && t < 1)
                                {
                                    // 大于sortedInnerBSDeepCopy[i].TurningTs.Last() 并且 小于 1
                                    sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments.Add(sortedBoundarySegmentsDeepCopy[i].TurningTs.Count - 1);
                                    lineIndex = sortedBoundarySegmentsDeepCopy[i].TurningTs.Count - 1;
                                    break;
                                }
                                else
                                {
                                    continue;
                                }
                            }

                            pointOnEachBS.Branch(i).Add(sortedBoundarySegmentsDeepCopy[i].PolylineCurve.PointAt(t));
                        }
                        else
                        {
                            sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments.Add(0);

                            pointOnEachBS.Branch(i).Add(sortedBoundarySegmentsDeepCopy[i].PolylineCurve.PointAt(t));
                        }
                    }
                }

                #endregion

                #region 计算tXDT和tYDT的数量
                //if (tXCounts[0] != tXList[0].Count)
                //{
                //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "tX的数量与预期的不符，请增加或减少tX的数量");
                //}
                //if (tYCounts[0] != InnerTList[0].Count)
                //{
                //    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "tY的数量与预期的不符，请增加或减少tY的数量");
                //}

                //DataTree<double> tXDT = new DataTree<double>();
                //DataTree<double> tYDT = new DataTree<double>();

                //int txCount = 0;
                //int tyCount = 0;
                //for (int i = 0; i < allLayerVertexIndexPairs[currentLayer].Count; i++)
                //{
                //    int interval = CalFirstBSInterval(allLayerBsIndexPairs[currentLayer][i], sortedBoundarySegments.Count);
                //    int[] tCount = CalTCount(interval,
                //                             sortedBoundarySegments[allLayerBsIndexPairs[currentLayer][i][0]],
                //                             sortedBoundarySegments[allLayerBsIndexPairs[currentLayer][i][1]]);


                //    tXDT.EnsurePath(i);
                //    tXDT.Branch(i).AddRange(tXList[0].Skip(txCount).Take(tCount[0]));
                //    tYDT.EnsurePath(i);
                //    tYDT.Branch(i).AddRange(InnerTList[0].Skip(tyCount).Take(tCount[1]));

                //    txCount += tCount[0];
                //    tyCount += tCount[1];
                //}
                //DA.SetDataTree(0, tXDT);
                //DA.SetDataTree(1, tYDT);
                #endregion

                #region 原来的tXDT和tYDT合并为innerTForEachBS
                int currentInnerTCount = 0;
                DataTree<double> innerTOnEachBSForCurrentLayer = new DataTree<double>();

                //int currentBranchTCount = 0;
                for (int i = 0; i < allLayerVertexIndexPairs[currentLayer].Count; i++)
                {
                    int interval = CalFirstBSInterval(allLayerBsIndexPairs[currentLayer][i], sortedBoundarySegmentsDeepCopy.Count);
                    int innerTCount = CalTCount(interval);
                                                //sortedBoundarySegmentsDeepCopy[allLayerBsIndexPairs[currentLayer][i][0]],
                                                //sortedBoundarySegmentsDeepCopy[allLayerBsIndexPairs[currentLayer][i][1]],
                                                //0,
                                                //0);

                    if (currentInnerTCount + innerTCount > allLayerInnerTMaxLoL[currentLayer].Count)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("InnerTMaxDT.Branch({0})的数量不足，请增加innerT的数量", currentLayer));
                    }

                    innerTOnEachBSForCurrentLayer.EnsurePath(i);
                    innerTOnEachBSForCurrentLayer.Branch(i).AddRange(allLayerInnerTMaxLoL[currentLayer].Skip(currentInnerTCount).Take(innerTCount));
                    //currentBranchTCount += innerTCount;
                    currentInnerTCount += innerTCount;
                }
                //// 增加目前currentInnerTCount的值，以备下次使用
                //currentInnerTCount += currentBranchTCount;
                #endregion

                List<int[]> onP_vertexIndexPairs;
                int iterCount;
                bool needToCalNextLayer;
                #region 构造边缘一圈的PlanktonXYZ列表
                #region 生成bs上的一圈点的列表
                List<Point3d> additionalPoint3ds = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegmentsDeepCopy.Count; i++)
                {
                    if (pointOnEachBS.Branch(i).Count == 0)
                    {
                        // 如果当前这个bs没有Point在它上面，就直接把每段的From点加上即可
                        for (int j = 0; j < sortedBoundarySegmentsDeepCopy[i].Lines.Count; j++)
                        {
                            additionalPoint3ds.Add(sortedBoundarySegmentsDeepCopy[i].Lines[j].From);
                        }
                    }
                    else
                    {
                        // 如果当前这个bs有Point在上面
                        for (int j = 0; j < sortedBoundarySegmentsDeepCopy[i].Lines.Count; j++)
                        {
                            // 首先先加上每个Line的起点
                            additionalPoint3ds.Add(sortedBoundarySegmentsDeepCopy[i].Lines[j].From);

                            if (sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments.Contains(j))
                            {
                                // 如果这段Line上有Point就添加Point
                                // 对于同一个segment上有几个点就循环几次
                                List<int> list = new List<int>();
                                for (int k = 0; k < sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments.Count; k++)
                                {
                                    if (sortedBoundarySegmentsDeepCopy[i].PointOnWhichSegments[k] == j)
                                    {
                                        list.Add(k);
                                    }
                                }

                                for (int k = 0; k < list.Count; k++)
                                {
                                    additionalPoint3ds.Add(pointOnEachBS.Branch(i)[list[k]]);
                                }
                            }
                            else
                            {
                                continue;
                            }
                        }
                    }
                }
                #endregion
                #region 构造PlanktonXYZ列表
                List<PlanktonXYZ> additionalPlanktonXYZs = new List<PlanktonXYZ>();
                for (int i = 0; i < additionalPoint3ds.Count; i++)
                {
                    additionalPlanktonXYZs.Add(new PlanktonXYZ((float)additionalPoint3ds[i].X, (float)additionalPoint3ds[i].Y, (float)additionalPoint3ds[i].Z));
                }
                #endregion
                #region 得到onD_vertexOrderFromBS
                List<int> onD_vertexOrderFromBS = new List<int>();
                for (int i = 0; i < verticesIndexForEachBSDT.BranchCount; i++)
                {
                    List<int> verticesIndex = new List<int>();
                    verticesIndex.AddRange(verticesIndexForEachBSDT.Branch(i));
                    verticesIndex.RemoveAt(verticesIndex.Count - 1);
                    onD_vertexOrderFromBS.AddRange(verticesIndex);
                }
                #endregion
                #region 得到onP_vertexOrderFromBS，即faceVertexOrder
                List<int> onP_vertexOrderFromBS = new List<int>();
                // originP目前没有Face，即第一次构造
                for (int i = 0; i < additionalPlanktonXYZs.Count; i++)
                {
                    onP_vertexOrderFromBS.Add(i);
                }

                List<List<int>> faceVertexOrder = new List<List<int>>();
                faceVertexOrder.Add(new List<int>());
                faceVertexOrder[0].AddRange(onP_vertexOrderFromBS);
                #endregion
                #endregion
                #region 得到最外圈的最初始的originP
                PlanktonMesh originP = new PlanktonMesh();
                originP.Vertices.AddVertices(additionalPlanktonXYZs);
                originP.Faces.AddFaces(faceVertexOrder);
                #endregion

                #region 对originP进行划分
                
                #region 得到当前这层的onP_vertexIndexPairs
                List<int[]> onD_vertexIndexPairs = allLayerVertexIndexPairs[currentLayer];
                onP_vertexIndexPairs = new List<int[]>();
                for (int i = 0; i < onD_vertexIndexPairs.Count; i++)
                {
                    int[] onP_pair = new int[2] { 0, 0 };
                    for (int j = 0; j < onD_vertexIndexPairs[i].Length; j++)
                    {
                        onP_pair[j] = onD_vertexOrderFromBS.IndexOf(onD_vertexIndexPairs[i][j]);
                    }
                    onP_vertexIndexPairs.Add(onP_pair);
                }
                #endregion

                #region 控制循环次数
                iterCount = 0;
                needToCalNextLayer = false;
                if (currentLayer == totalLayerCount - 1)
                {
                    // 注意这里 onP_vertexIndexPairs.Count - 1，是为了不计算最后一个pair
                    iterCount = onP_vertexIndexPairs.Count - 1;
                    needToCalNextLayer = false;
                }
                else
                {
                    iterCount = onP_vertexIndexPairs.Count;
                    needToCalNextLayer = true;
                }
                #endregion

                PlanktonMesh P = new PlanktonMesh(originP);
                int currentPFaceCount = P.Faces.Count;
                for (int i = 0; i < iterCount; i++)
                {
                    int interval = CalFirstBSInterval(allLayerBsIndexPairs[currentLayer][i], sortedBoundarySegmentsDeepCopy.Count);

                    bool isLastPair = false;
                    int thisLayerFirstFaceIndex_OnP = currentPFaceCount - 1;
                    if (needToCalNextLayer && i == iterCount - 1)
                    {
                        isLastPair = true;
                    }

                    P = MakeStepLine(P,
                                     interval,
                                     onP_vertexIndexPairs[i][0],
                                     onP_vertexIndexPairs[i][1],
                                     sortedBoundarySegmentsDeepCopy[allLayerBsIndexPairs[currentLayer][i][0]],
                                     sortedBoundarySegmentsDeepCopy[allLayerBsIndexPairs[currentLayer][i][1]],
                                     0,
                                     0,
                                     -1,
                                     boundingBox,
                                     innerTOnEachBSForCurrentLayer,
                                     i,
                                     isLastPair,
                                     thisLayerFirstFaceIndex_OnP);

                    List<string> debugString2 = UtilityFunctions.PrintFacesVertices(P);
                    List<string> debugString3 = UtilityFunctions.PrintHalfedges(P);
                    List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(P);
                }
                #endregion
                List<string> debugString5 = UtilityFunctions.PrintFacesVertices(P);
                List<string> debugString6 = UtilityFunctions.PrintHalfedges(P);
                List<string> debugString7 = UtilityFunctions.PrintHalfedgeStartAndEnd(P);

                #endregion

                while (needToCalNextLayer)
                {
                    currentLayer++;

                    

                    // 首先针对现有的PlanktonMesh中，上一层的划分中，每个vertexIndexPair所生成的面，所产生的新的内部的bs
                    // onP_vertexIndexPairsCorrespondingFaceIndex与onD_vertexIndexPairsCorrespondingFaceIndex存在映射关系
                    List<int> onP_vertexIndexPairsCorrespondingFaceIndex = new List<int>();
                    for (int i = 0; i < onP_vertexIndexPairs.Count; i++)
                    {
                        for (int j = 0; j < P.Faces.Count; j++)
                        {
                            int[] faceVertice = P.Faces.GetFaceVertices(j);
                            int[] intersectResult = faceVertice.Intersect(onP_vertexIndexPairs[i]).ToArray();
                            if (intersectResult.Length == onP_vertexIndexPairs[i].Length)
                            {
                                onP_vertexIndexPairsCorrespondingFaceIndex.Add(j);
                            }
                        }
                    }

                    #region 找到每个边缘面和中心面的公共边的列表
                    List<List<int>> onP_publicHsList = new List<List<int>>();
                    for (int i = 0; i < onP_vertexIndexPairsCorrespondingFaceIndex.Count; i++)
                    {
                        List<int> publicHs = FindPublicHalfEdges(P, P.Faces.Count - 1, onP_vertexIndexPairsCorrespondingFaceIndex[i]);
                        onP_publicHsList.Add(publicHs);
                    }
                    #endregion

                    #region 每个公共边对应的onP_vertexOrderFromBS
                    List<List<int>> onP_vertexOrderForEachBS = new List<List<int>>();
                    for (int i = 0; i < onP_publicHsList.Count; i++)
                    {
                        onP_vertexOrderForEachBS.Add(new List<int>());

                        List<int> vertexOrder = new List<int>();
                        for (int j = 0; j < onP_publicHsList[i].Count; j++)
                        {
                            int end = P.Halfedges.EndVertex(onP_publicHsList[i][j]);
                            vertexOrder.Add(end);
                        }
                        int start = P.Halfedges[onP_publicHsList[i].Last()].StartVertex;
                        vertexOrder.Add(start);

                        onP_vertexOrderForEachBS[i].AddRange(vertexOrder);
                    }
                    #endregion

                    #region 由公共边的列表，得到innerBS的列表
                    List<BoundarySegment> sortedInnerBS = new List<BoundarySegment>();
                    for (int i = 0; i < onP_vertexIndexPairsCorrespondingFaceIndex.Count; i++)
                    {
                        List<int> publicHs = onP_publicHsList[i];

                        // 对innerBS中的Lines大于1的转折BS进行翻转
                        if (publicHs.Count > 1)
                        {
                            // publicHs.Reverse();
                            // 对innerBS中的Lines大于1的转折BS进行翻转
                            for (int j = 0; j < publicHs.Count; j++)
                            {
                                if (j % 2 == 0 && j != publicHs.Count - 1)
                                {
                                    // 对于偶数序号的publicH，找到对应的start序号
                                    int vertexToFlip = P.Halfedges[publicHs[j]].StartVertex;
                                    Point3d pointToFlip = P.Vertices[vertexToFlip].ToPoint3d();

                                    int prevVertexToFlip = P.Halfedges[P.Halfedges[publicHs[j]].PrevHalfedge].StartVertex;
                                    int nextVertexToFlip = P.Halfedges.EndVertex(publicHs[j]);
                                    Point3d prevPointToFlip = P.Vertices[prevVertexToFlip].ToPoint3d();
                                    Point3d nextPointToFlip = P.Vertices[nextVertexToFlip].ToPoint3d();

                                    Point3d newPointToFlip;
                                    if (pointToFlip.X == nextPointToFlip.X)
                                    {
                                        newPointToFlip = new Point3d(prevPointToFlip.X, nextPointToFlip.Y, prevPointToFlip.Z);
                                    }
                                    else
                                    {
                                        newPointToFlip = new Point3d(nextPointToFlip.X, prevPointToFlip.Y, prevPointToFlip.Z);
                                    }

                                    P.Vertices.SetVertex(vertexToFlip, newPointToFlip);
                                }
                            }
                        }


                        List<Line> lines = new List<Line>();
                        for (int j = 0; j < publicHs.Count; j++)
                        {
                            int start = P.Halfedges.EndVertex(publicHs[j]);
                            int end = P.Halfedges[publicHs[j]].StartVertex;
                            Point3d startPoint = P.Vertices[start].ToPoint3d();
                            Point3d endPoint = P.Vertices[end].ToPoint3d();
                            Line line = new Line(startPoint, endPoint);
                            lines.Add(line);
                        }
                        BoundarySegment bs = new BoundarySegment(lines, string.Format("PFace:{0}", onP_vertexIndexPairsCorrespondingFaceIndex[i]));
                        sortedInnerBS.Add(bs);
                    }
                    #endregion

                    if (allLayerVerticesIndexForEachBS[currentLayer].DataCount == 0)
                    {
                        break;
                    }


                    P = GenerateInnerPart(allLayerVerticesIndexForEachBS[currentLayer],
                                          allLayerVolumeJunctionsIndexOnEachBS[currentLayer],
                                          allLayerBoundaryTLoL[currentLayer],
                                          allLayerInnerTMaxLoL[currentLayer],
                                          //tYCounts[currentLayer],
                                          //tXList[currentLayer],

                                          allLayerVertexIndexPairs[currentLayer],
                                          allLayerBsIndexPairs[currentLayer],

                                          onP_vertexOrderForEachBS,
                                          onP_vertexIndexPairsCorrespondingFaceIndex,

                                          P,
                                          currentLayer,
                                          totalLayerCount,
                                          sortedInnerBS,
                                          ref needToCalNextLayer);
                }

                DA.SetData("DualPlanktonMesh", P);

                List<GraphNode> pGraphNodes = dualGraphWithHM.GraphNodes;
                List<List<int>> pGraphTables = dualGraphWithHM.GraphTables;
                int innerNodeCount = dualGraphWithHM.InnerNodeCount;
                int outerNodeCount = dualGraphWithHM.OuterNodeCount;

                #region 更新原来图的相关信息
                #region 添加对偶图所有面的中心点作为innerGraphNodePoints的坐标位置
                Polyline[] polylines = P.ToPolylines();
                List<Point3d> innerGraphNodePoints = new List<Point3d>();
                for (int i = 0; i < polylines.Length; i++)
                {
                    List<Point3d> points = polylines[i].ToList();
                    Point3d gravityPoint = GetCenterOfGravityPoint(points);
                    innerGraphNodePoints.Add(gravityPoint);
                    //GraphNodePoints.Add(gravityPoint);
                }
                #endregion

                #region outerGraphNodePoints的坐标位置
                // 表示NEWS的四个点
                List<Point3d> outerGraphNodePoints = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    Vector3d verticalVector = Vector3d.CrossProduct(sortedBoundarySegments[i].Lines[0].Direction, Vector3d.ZAxis);
                    verticalVector.Unitize();

                    Point3d point0 = (sortedBoundarySegments[i].From + sortedBoundarySegments[i].To) / 2;
                    point0 += verticalVector * 3 * sortedBoundarySegments[i].Lines[0].Length / 4;
                    outerGraphNodePoints.Add(point0);
                }
                outerGraphNodePoints.Insert(0, outerGraphNodePoints.Last());
                outerGraphNodePoints.RemoveAt(outerGraphNodePoints.Count - 1);
                #endregion

                #region 把LOL形式的allLayerPairCorrespondingFaceIndexLoL转换为list形式
                List<int> allLayerPairCorrespondingFaceIndex_OnD = new List<int>();
                for (int i = 0; i < allLayerPairCorrespondingFaceIndexLoL_OnD.Count; i++)
                {
                    allLayerPairCorrespondingFaceIndex_OnD.AddRange(allLayerPairCorrespondingFaceIndexLoL_OnD[i]);
                }
                #endregion
                if (allLayerVerticesIndexForEachBS.Last().DataCount == 0)
                {
                    allLayerPairCorrespondingFaceIndex_OnD.Add(P.Faces.Count - 1);
                }

                #region newNodePoints
                List<Point3d> newNodePoints = new List<Point3d>();
                for (int i = 0; i < pGraphNodes.Count; i++)
                {
                    if (pGraphNodes[i].IsInner)
                    {
                        newNodePoints.Add(innerGraphNodePoints[allLayerPairCorrespondingFaceIndex_OnD.IndexOf(i)]);
                    }
                    else
                    {
                        newNodePoints.Add(outerGraphNodePoints[i - innerNodeCount]);
                    }
                }
                #endregion
                for (int i = 0; i < pGraphNodes.Count; i++)
                {
                    pGraphNodes[i].NodeVertex = newNodePoints[i];
                }

                Graph newGraph = new Graph(pGraphNodes, pGraphTables);
                DataTree<int> newGraphDT = UtilityFunctions.LoLToDataTree<int>(newGraph.GraphTables);
                List<Line> graphEdge = GraphEdgeList(newGraphDT, newNodePoints);
                //List<Polyline> graphPolyline = GMeshFaceBoundaries(outerGraphNodePoints, graphEdge, Plane.WorldXY);

                //List<List<int>> planktonMeshFaceVertexOrder = new List<List<int>>();
                //for (int i = 0; i < graphPolyline.Count; i++)
                //{
                //    planktonMeshFaceVertexOrder.Add(new List<int>());
                //    for (int j = 0; j < graphPolyline[i].Count; j++)
                //    {
                //        int vertexIndex = newNodePoints.IndexOf(graphPolyline[i][j]);
                //        planktonMeshFaceVertexOrder[i].Add(vertexIndex);
                //    }
                //}
                //PlanktonMesh planktonMeshForGraph = new PlanktonMesh();
                //for (int i = 0; i < newNodePoints.Count; i++)
                //{
                //    planktonMeshForGraph.Vertices.Add(newNodePoints[i].X, newNodePoints[i].Y, newNodePoints[i].Z);
                //}
                //planktonMeshForGraph.Faces.AddFaces(planktonMeshFaceVertexOrder);

                // 更新DualGraphWithHM对象中的PlanktonMesh属性
                PlanktonMesh originPlanktonMesh = new PlanktonMesh(dualGraphWithHM.PlanktonMesh);
                PlanktonMesh newPlanktonMesh = PlanktonGh.RhinoSupport.ReplaceVertices(originPlanktonMesh, newNodePoints);

                DualGraphWithHM newDualGraphWithHM = new DualGraphWithHM(P, newPlanktonMesh, pGraphNodes, pGraphTables, dualGraphWithHM.DFaceIndexsAroundOuterNodes);
                List<int> DFace_PVertice = allLayerPairCorrespondingFaceIndex_OnD;
                newDualGraphWithHM.DFIndex_PVIndex = new List<int>();
                newDualGraphWithHM.DFIndex_PVIndex.AddRange(DFace_PVertice);

                //GraphWithHM newGraphWithHM = new GraphWithHM(planktonMeshForGraph, newGraph.GraphNodes, newGraph.GraphTables);
                DA.SetData("DualGraphWithHM", newDualGraphWithHM);
                #endregion


                #region 可视化部分
                //PlanktonMesh meshForVisualize = new PlanktonMesh();
                //meshForVisualize = P;
                #region clear
                // 对偶图
                //DualVertices.Clear();
                DualPolylines.Clear();
                DualVertexTextDots.Clear();
                // 原来的图
                //GraphNodePoints.Clear();
                GraphEdges.Clear();
                InnerNodeTextDots.Clear();

                MeshForVisualize = null;
                DottedCurve.Clear();
                #endregion

                #region 对偶图
                //DualVertices = PlanktonGh.RhinoSupport.GetPositions(MeshForVisualize).ToList();
                List<Point3d> dualVertices = PlanktonGh.RhinoSupport.GetPositions(P).ToList();

                Polyline[] dualPolylines = PlanktonGh.RhinoSupport.ToPolylines(P);
                for (int i = 0; i < dualPolylines.Length; i++)
                {
                    DualPolylines.Add(new List<Point3d>());
                    DualPolylines[i].AddRange(dualPolylines[i]);
                }


                for (int i = 0; i < P.Vertices.Count; i++)
                {
                    List<int> facesAroundVertex = P.Vertices.GetVertexFaces(i).ToList();
                    while (facesAroundVertex.Contains(-1))
                    {
                        facesAroundVertex.Remove(-1);
                    }

                    List<int> VertexCorrespondingFaceIndice_OnP = new List<int>();
                    for (int j = 0; j < facesAroundVertex.Count; j++)
                    {
                        int correspondingDFaceIndex = allLayerPairCorrespondingFaceIndex_OnD[facesAroundVertex[j]];
                        VertexCorrespondingFaceIndice_OnP.Add(correspondingDFaceIndex);
                    }

                    string arg = string.Join(";", VertexCorrespondingFaceIndice_OnP);
                    
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), dualVertices[i]);
                    DualVertexTextDots.Add(textDot);
                }



                #endregion

                #region 原来的图

                #region 实线和虚线
                // 实线
                // 绘制Graph的Edge
                GraphEdges = graphEdge;
                // 绘制Graph形成的Convex
                //List<Polyline> graphConvex = GMeshFaceBoundaries(outerGraphNodePoints, graphEdge, Plane.WorldXY);

                MeshForVisualize = newPlanktonMesh.ToRhinoMesh();

                // 虚线
                List<Line> meshEdges = new List<Line>();
                for (int i = 0; i < MeshForVisualize.TopologyEdges.Count; i++)
                {
                    meshEdges.Add(MeshForVisualize.TopologyEdges.EdgeLine(i));
                }
                double[] pattern = { 1.0 };
                for (int i = 0; i < meshEdges.Count; i++)
                {
                    IEnumerable<Curve> segments = UtilityFunctions.ApplyDashPattern(meshEdges[i].ToNurbsCurve(), pattern);
                    DottedCurve.AddRange(segments);
                }


                #endregion

                //GraphNodePoints = newNodePoints;


                #region TextDot
                // outer
                //List<string> outerNodeLabels = new List<string>();
                //for (int i = 0; i < pGraphNodes.Count; i++)
                //{
                //    if (!pGraphNodes[i].IsInner)
                //    {
                //        outerNodeLabels.Add(pGraphNodes[i].NodeAttribute.NodeLabel);
                //    }
                //}
                //for (int i = 0; i < sortedBoundarySegments.Count; i++)
                //{
                //    TextDot boundarySegmentTextDot = new TextDot(string.Format("{0} | {1}", i + innerNodeCount, outerNodeLabels[i]), outerGraphNodePoints[i]);
                //    OuterNodeTextDots.Add(boundarySegmentTextDot);
                //}

                for (int i = 0; i < newGraph.GraphNodes.Count; i++)
                {
                    if (!newGraph.GraphNodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, newGraph.GraphNodes[i].NodeAttribute.NodeLabel), newGraph.GraphNodes[i].NodeVertex);
                        OuterNodeTextDots.Add(textDot);
                    }
                }

                // inner
                for (int i = 0; i < newGraph.GraphNodes.Count; i++)
                {
                    if (newGraph.GraphNodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, newGraph.GraphNodes[i].NodeAttribute.NodeLabel), newGraph.GraphNodes[i].NodeVertex);
                        InnerNodeTextDots.Add(textDot);
                    }
                }



                //List<GraphNode> graphNodes = new List<GraphNode>();
                //for (int i = 0; i < MeshForVisualize.Faces.Count; i++)
                //{
                //    // FaceCorrespondingGraphNodeIndex是OnD
                //    int nodeIndex = allLayerPairCorrespondingFaceIndex_OnD[dualGraphWithHM.FaceCorrespondingGraphNodeIndex[i]];
                //    graphNodes.Add(dualGraphWithHM.GraphNodes[nodeIndex]);
                //}

                //for (int i = 0; i < MeshForVisualize.Faces.Count; i++)
                //{
                //    InnerNodeTextDots.Add(new TextDot(string.Format("{0} | {1}", allLayerPairCorrespondingFaceIndex_OnD[i], graphNodes[i].NodeAttribute.NodeLabel), GraphNodePoints[i]));
                //}
                #endregion

                #endregion
                #endregion
            }
        }

        private PlanktonMesh GenerateInnerPart(DataTree<int> verticesIndexForEachInnerBSDT,
                                               DataTree<int> volumeJunctionsIndexOnEachInnerBS,
                                               List<double> boundaryTList,
                                               List<double> innerTListCurrentLayer,
                                               List<int[]> onD_innerVertexIndexPairs,
                                               List<int[]> currentBSIndexPairs,
                                               List<List<int>> onP_vertexOrderForEachBS,
                                               List<int> onP_correspondingFaceIndex,
                                               PlanktonMesh P,
                                               int currentLayer,
                                               int totalLayerCount,
                                               List<BoundarySegment> sortedInnerBS,
                                               ref bool needToCalNextLayer)
        {
            // 深拷贝P
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);

            // 深拷贝sortedInnerBS
            List<BoundarySegment> sortedInnerBSDeepCopy = new List<BoundarySegment>();
            for (int i = 0; i < sortedInnerBS.Count; i++)
            {
                sortedInnerBSDeepCopy.Add(new BoundarySegment(sortedInnerBS[i]));
            }

            #region BoundingBox
            List<Point3d> cornerPoints = new List<Point3d>();
            for (int i = 0; i < sortedInnerBSDeepCopy.Count; i++)
            {
                cornerPoints.AddRange(sortedInnerBSDeepCopy[i].Points);
            }

            BoundingBox currentBoundingBox = new BoundingBox(cornerPoints);
            #endregion

            #region 构造DataTree<double> tOnEachBS
            // 每个包含volumeJunction的bs为了确定volumeJunction位置所需要的t的datatree
            if (volumeJunctionsIndexOnEachInnerBS.DataCount > boundaryTList.Count)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("boundaryTList[{0}]的数量不足，请增加boundaryT的数量", currentLayer));
            }

            DataTree<double> boundaryTOnEachInnerBS = new DataTree<double>();
            int indexInner = 0;
            for (int i = 0; i < volumeJunctionsIndexOnEachInnerBS.BranchCount; i++)
            {
                boundaryTOnEachInnerBS.EnsurePath(i);
                for (int j = 0; j < volumeJunctionsIndexOnEachInnerBS.Branch(i).Count; j++)
                {
                    boundaryTOnEachInnerBS.Branch(i).Add(boundaryTList[indexInner]);
                    indexInner++;
                }
            }

            for (int i = 0; i < boundaryTOnEachInnerBS.BranchCount; i++)
            {
                List<double> ts = boundaryTOnEachInnerBS.Branch(i);
                for (int j = 0; j < ts.Count - 1; j++)
                {
                    if (ts[j] > ts[j + 1])
                    {
                        // 后一个t不能大于前一个t
                        return null;
                    }
                }
            }
            #endregion

            #region DataTree<Point3d> pointOnEachBS
            int currentLastVertexIndex = P.Vertices.Count;
            List<int> onP_newInnerVolumeJunctionVertexIndex = new List<int>();

            DataTree<Point3d> pointOnEachBS = new DataTree<Point3d>();
            for (int i = 0; i < boundaryTOnEachInnerBS.BranchCount; i++)
            {
                pointOnEachBS.EnsurePath(i);
                for (int j = 0; j < boundaryTOnEachInnerBS.Branch(i).Count; j++)
                {
                    double t = boundaryTOnEachInnerBS.Branch(i)[j];
                    if (sortedInnerBSDeepCopy[i].TurningTs.Count != 0)
                    {
                        int lineIndex = 0;
                        for (int k = 0; k < sortedInnerBSDeepCopy[i].TurningTs.Count - 1; k++)
                        {
                            if (t > sortedInnerBSDeepCopy[i].TurningTs[k] && t < sortedInnerBSDeepCopy[i].TurningTs[k + 1])
                            {
                                // 注意这里一定要添加PointOnWhichSegments信息，即这个点是在BS的哪一段上
                                sortedInnerBSDeepCopy[i].PointOnWhichSegments.Add(k);
                                lineIndex = k;
                                break;
                            }
                            else if (t > sortedInnerBSDeepCopy[i].TurningTs.Last() && t < 1)
                            {
                                // 大于sortedInnerBSDeepCopy[i].TurningTs.Last() 并且 小于 1
                                sortedInnerBSDeepCopy[i].PointOnWhichSegments.Add(sortedInnerBSDeepCopy[i].TurningTs.Count - 1);
                                lineIndex = sortedInnerBSDeepCopy[i].TurningTs.Count - 1;
                                break;
                            }
                            else
                            {
                                continue;
                            }
                        }

                        pointOnEachBS.Branch(i).Add(sortedInnerBSDeepCopy[i].PolylineCurve.PointAt(t));

                        onP_newInnerVolumeJunctionVertexIndex.Add(currentLastVertexIndex);
                        currentLastVertexIndex++;
                    }
                    else
                    {
                        sortedInnerBSDeepCopy[i].PointOnWhichSegments.Add(0);

                        pointOnEachBS.Branch(i).Add(sortedInnerBSDeepCopy[i].PolylineCurve.PointAt(t));

                        onP_newInnerVolumeJunctionVertexIndex.Add(currentLastVertexIndex);
                        currentLastVertexIndex++;
                    }
                }
            }
            #endregion

            List<int> shifted_onP_correspondingFaceIndex = new List<int>();
            int shift = 0;
            for (int i = 0; i < pointOnEachBS.BranchCount; i++)
            {
                if (pointOnEachBS.Branch(i).Count != 0)
                {
                    shift = i;
                    break;
                }
            }
            shifted_onP_correspondingFaceIndex = LeftShift<int>(onP_correspondingFaceIndex, shift);


            #region 改造一下sortedInnerBSDeepCopy
            List<Line> sortedInnerLines = new List<Line>();
            // Line是累加，segment是单独的
            List<int> pointOnWhichLineIndexList = new List<int>();
            List<int> pointOnWhichSegmentIndexList = new List<int>();
            int prevLineCount = 0;
            for (int i = 0; i < sortedInnerBSDeepCopy.Count; i++)
            {
                // 如果bs上有volumeJunction点，就记录line和segment的序号（累加或单独的）
                if (sortedInnerBSDeepCopy[i].PointOnWhichSegments.Count != 0)
                {
                    for (int j = 0; j < sortedInnerBSDeepCopy[i].PointOnWhichSegments.Count; j++)
                    {
                        int lineIndex = prevLineCount + sortedInnerBSDeepCopy[i].PointOnWhichSegments[j];
                        pointOnWhichLineIndexList.Add(lineIndex);

                        pointOnWhichSegmentIndexList.Add(sortedInnerBSDeepCopy[i].PointOnWhichSegments[j]);
                    }
                }

                // 对每个bs，把bs中的每个line都取出来
                for (int j = 0; j < sortedInnerBSDeepCopy[i].Lines.Count; j++)
                {
                    Line line = sortedInnerBSDeepCopy[i].Lines[j];
                    sortedInnerLines.Add(line);
                }
                prevLineCount += sortedInnerBSDeepCopy[i].Lines.Count;
            }

            List<int[]> pointOnWhichLineIndexPairs = new List<int[]>();
            List<int[]> pointOnWhichSegmentIndexPairs = new List<int[]>();
            for (int i = 0; i < pointOnWhichLineIndexList.Count; i++)
            {
                int[] pair = new int[2] { pointOnWhichLineIndexList[i], pointOnWhichLineIndexList[(i + 1) % pointOnWhichLineIndexList.Count] };
                pointOnWhichLineIndexPairs.Add(pair);

                int[] pair2 = new int[2] { pointOnWhichSegmentIndexList[i], pointOnWhichSegmentIndexList[(i + 1) % pointOnWhichSegmentIndexList.Count] };
                pointOnWhichSegmentIndexPairs.Add(pair2);
            }

            // 找到转折点turningPoint
            List<int> turningIndex = new List<int>();
            for (int i = 0; i < sortedInnerLines.Count - 1; i++)
            {
                Vector3d vector1 = new Vector3d(sortedInnerLines[i].To - sortedInnerLines[i].From);
                Vector3d vector2 = new Vector3d(sortedInnerLines[i + 1].To - sortedInnerLines[i + 1].From);

                // 当前直线和下一条直线不平行时，当前的index为turningIndex
                if (vector1.X * vector2.Y - vector1.Y * vector2.X != 0)
                {
                    turningIndex.Add(i);
                }
            }
            #endregion

            #region 得到onD_vertexOrderFromBS
            List<int> onD_vertexOrderFromInnerBS = new List<int>();
            for (int i = 0; i < verticesIndexForEachInnerBSDT.BranchCount; i++)
            {
                List<int> verticesIndex = new List<int>();
                verticesIndex.AddRange(verticesIndexForEachInnerBSDT.Branch(i));
                verticesIndex.RemoveAt(verticesIndex.Count - 1);
                onD_vertexOrderFromInnerBS.AddRange(verticesIndex);
            }
            #endregion

            #region 得到onP_vertexOrderFromInnerBS
            int currentInnerVolumeJunctionVertexIndex = 0;
            List<int> onP_vertexOrderFromInnerBS = new List<int>();
            for (int i = 0; i < sortedInnerBSDeepCopy.Count; i++)
            {
                if (sortedInnerBSDeepCopy[i].PointOnWhichSegments.Count > 0)
                {
                    // 如果该innerBS上有volumeJunctionPoint
                    for (int j = 0; j < sortedInnerBSDeepCopy[i].PointOnWhichSegments.Count; j++)
                    {
                        int insertIndex = sortedInnerBSDeepCopy[i].PointOnWhichSegments[j] + 1;
                        int insertItem = onP_newInnerVolumeJunctionVertexIndex[currentInnerVolumeJunctionVertexIndex];

                        onP_vertexOrderForEachBS[i].Insert(insertIndex, insertItem);
                        currentInnerVolumeJunctionVertexIndex++;
                    }
                    List<int> list = new List<int>();
                    list.AddRange(onP_vertexOrderForEachBS[i]);
                    list.RemoveAt(list.Count - 1);
                    onP_vertexOrderFromInnerBS.AddRange(list);
                }
                else
                {
                    // 如果该innerBS上没有volumeJunctionPoint，就直接加全部的vertex序号就好了
                    List<int> list = new List<int>();
                    list.AddRange(onP_vertexOrderForEachBS[i]);
                    list.RemoveAt(list.Count - 1);
                    onP_vertexOrderFromInnerBS.AddRange(list);
                }
            }
            #endregion

            #region 对当前的P进行划分

            #region 向当前P的InnerBS那一圈添加volumeJunction点
            for (int i = 0; i < pointOnEachBS.BranchCount; i++)
            {
                for (int j = 0; j < pointOnEachBS.Branch(i).Count; j++)
                {
                    P.Vertices.Add(pointOnEachBS.Branch(i)[j]);
                }
            }
            List<PlanktonXYZ> planktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < P.Vertices.Count; i++)
            {
                planktonXYZ.Add(P.Vertices[i].ToXYZ());
            }

            List<List<int>> faceVertexOrder = new List<List<int>>();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                faceVertexOrder.Add(new List<int>());
                if (onP_correspondingFaceIndex.Contains(i))
                {
                    List<int> vi_old, vi_new;
                    int vertex1 = onP_vertexOrderForEachBS[i][0];
                    int vertex2 = onP_vertexOrderForEachBS[i].Last();
                    GetInitialViLists(vertex1, vertex2, P, i, out vi_old, out vi_new);

                    List<int> all = new List<int>();
                    vi_old.RemoveAt(0);
                    vi_old.RemoveAt(vi_old.Count - 1);

                    List<int> add = new List<int>();
                    add.AddRange(onP_vertexOrderForEachBS[i]);
                    add.Reverse();

                    all.AddRange(vi_old);
                    all.AddRange(add);

                    faceVertexOrder[i].AddRange(all);
                }
                else if (i == P.Faces.Count - 1)
                {
                    // 中心最新的Face
                    List<int> all = new List<int>();
                    all.AddRange(onP_vertexOrderFromInnerBS);
                    faceVertexOrder[i].AddRange(all);
                }
                else
                {
                    List<int> all = P.Faces.GetFaceVertices(i).ToList();
                    faceVertexOrder[i].AddRange(all);
                }
            }

            #region 构造新的P
            PDeepCopy = new PlanktonMesh();
            PDeepCopy.Vertices.AddVertices(planktonXYZ);
            PDeepCopy.Faces.AddFaces(faceVertexOrder);
            #endregion
            #endregion

            List<string> debugString8 = UtilityFunctions.PrintFacesVertices(PDeepCopy);
            List<string> debugString9 = UtilityFunctions.PrintHalfedges(PDeepCopy);
            List<string> debugString10 = UtilityFunctions.PrintHalfedgeStartAndEnd(PDeepCopy);

            #region 得到当前这层的onP_vertexIndexPairs
            //List<int[]> onD_innerVertexIndexPairs = currentVertexIndexPairs;
            List<int[]> onP_innerVertexIndexPairs = new List<int[]>();
            for (int i = 0; i < onP_newInnerVolumeJunctionVertexIndex.Count; i++)
            {
                int[] pair = new int[2] { onP_newInnerVolumeJunctionVertexIndex[i], onP_newInnerVolumeJunctionVertexIndex[(i + 1) % onP_newInnerVolumeJunctionVertexIndex.Count] };
                onP_innerVertexIndexPairs.Add(pair);
            }
            #endregion

            #region 控制循环次数
            int iterCountInner = 0;
            needToCalNextLayer = false;
            if (currentLayer == totalLayerCount - 1)
            {
                // 注意这里 onP_vertexIndexPairs.Count - 1，是为了不计算最后一个pair
                iterCountInner = onP_innerVertexIndexPairs.Count - 1;
                needToCalNextLayer = false;
            }
            else
            {
                iterCountInner = onP_innerVertexIndexPairs.Count;
                needToCalNextLayer = true;
            }
            #endregion

            #region 原来的tXDT和tYDT合并为innerTForEachBS
            // currentInnerTCount是对当前DataTree<double> innerTOnEachBSForCurrentLayer的总记录
            // 用它来作为起始点，在List<double> innerTListCurrentLayer上截取innerTCount长度
            int currentInnerTCount = 0;
            DataTree<double> innerTOnEachBSForCurrentLayer = new DataTree<double>();

            //int currentBranchTCount = 0;
            for (int i = 0; i < iterCountInner; i++)
            {
                int interval = CalInnerBSInterval(pointOnWhichLineIndexPairs[i], turningIndex);
                int innerTCount = CalTCount(interval);

                if (currentInnerTCount + innerTCount > innerTListCurrentLayer.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("innerTList[{0}]的数量不足，请增加innerT的数量", currentLayer));
                }

                innerTOnEachBSForCurrentLayer.EnsurePath(i);
                innerTOnEachBSForCurrentLayer.Branch(i).AddRange(innerTListCurrentLayer.Skip(currentInnerTCount).Take(innerTCount));
                //currentBranchTCount += currentBranchTCount;
                currentInnerTCount += innerTCount;
            }
            //// 增加目前currentInnerTCount的值，以备下次使用
            //currentInnerTCount += currentBranchTCount;
            #endregion

            int currentPFaceCount = PDeepCopy.Faces.Count;
            for (int i = 0; i < iterCountInner; i++)
            {
                int interval = CalInnerBSInterval(pointOnWhichLineIndexPairs[i], turningIndex);

                bool isLastPair = false;
                int thisLayerFirstFaceIndex_OnP = currentPFaceCount - 1;
                if (needToCalNextLayer && i == iterCountInner - 1)
                {
                    isLastPair = true;
                }

                PDeepCopy = MakeStepLine(PDeepCopy,
                                         interval,
                                         onP_innerVertexIndexPairs[i][0],
                                         onP_innerVertexIndexPairs[i][1],

                                         sortedInnerBSDeepCopy[currentBSIndexPairs[i][0]],
                                         sortedInnerBSDeepCopy[currentBSIndexPairs[i][1]],

                                         pointOnWhichSegmentIndexPairs[i][0],
                                         pointOnWhichSegmentIndexPairs[i][1],

                                         shifted_onP_correspondingFaceIndex[i],
                                         currentBoundingBox,
                                         innerTOnEachBSForCurrentLayer,
                                         i,
                                         isLastPair,
                                         thisLayerFirstFaceIndex_OnP);

                List<string> debugString2 = UtilityFunctions.PrintFacesVertices(PDeepCopy);
                List<string> debugString3 = UtilityFunctions.PrintHalfedges(PDeepCopy);
                List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(PDeepCopy);
            }
            #endregion
            List<string> debugString5 = UtilityFunctions.PrintFacesVertices(PDeepCopy);
            List<string> debugString6 = UtilityFunctions.PrintHalfedges(PDeepCopy);
            List<string> debugString7 = UtilityFunctions.PrintHalfedgeStartAndEnd(PDeepCopy);

            return PDeepCopy;
        }

        private List<int> FindPublicHalfEdges(PlanktonMesh P, int face1,int face2)
        {
            List<int> face1Hs = P.Faces.GetHalfedges(face1).ToList();
            List<int> face1PairHs = new List<int>();
            for (int i = 0; i < face1Hs.Count; i++)
            {
                face1PairHs.Add(P.Halfedges.GetPairHalfedge(face1Hs[i]));
            }

            List<int> face2Hs = P.Faces.GetHalfedges(face2).ToList();
            List<int> publicHs = face1PairHs.Intersect(face2Hs).ToList();

            List<int> publicHsIndexInFace1Hs = new List<int>();
            for (int i = 0; i < publicHs.Count; i++)
            {
                publicHsIndexInFace1Hs.Add(face1PairHs.IndexOf(publicHs[i]));
            }

            // 冒泡排序
            for (int i = 0; i < publicHs.Count; i++)
            {
                for (int j = 0; j < publicHs.Count - i - 1; j++)
                {
                    if (publicHsIndexInFace1Hs[j] > publicHsIndexInFace1Hs[j+1])
                    {
                        int temp = publicHs[j];
                        publicHs[j] = publicHs[j + 1];
                        publicHs[j + 1] = temp;
                    }
                }
            }

            // 生成的publicH是属于边缘面的
            return publicHs;
        }


        /// <summary>
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="bsIndexPair"></param>
        /// <param name="bsCount"></param>
        /// <returns></returns>
        private int CalFirstBSInterval(int[] bsIndexPair, int bsCount)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            int index0 = bsIndexPair[0];
            int index1 = bsIndexPair[1];

            if (index0 < index1)
            {
                counterClockwiseCount = index1 - index0;
                clockwiseCount = bsCount - index1 + index0;
            }
            else
            {
                counterClockwiseCount = bsCount - index0 + index1;
                clockwiseCount = index0 - index1;
            }

            if (counterClockwiseCount < clockwiseCount)
            {
                return counterClockwiseCount;
            }
            else
            {
                return clockwiseCount;
            }
        }

        /// <summary>
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="pointOnWhichLineIndexPairs"></param>
        /// <param name="turningIndex"></param>
        /// <returns></returns>
        private int CalInnerBSInterval(int[] pointOnWhichLineIndexPairs, List<int> turningIndex)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            //int index0;
            //int index1;

            List<int[]> turningRange = new List<int[]>();
            for (int i = 0; i < turningIndex.Count - 1; i++)
            {
                int[] pair = new int[2] { turningIndex[i], turningIndex[(i + 1) % turningIndex.Count] };
                turningRange.Add(pair);
            }

            // 需要注意等于的位置
            List<int> indexs = new List<int>();
            for (int i = 0; i < pointOnWhichLineIndexPairs.Length; i++)
            {
                // 一共有四个区间
                // 小于turningIndex[0]是第一个区间，对应四边形的第0条边
                if (pointOnWhichLineIndexPairs[i] <= turningIndex[0])
                {
                    indexs.Add(0);
                }
                // 大于turningIndex.Last()是第四个区间，对应四边形的第3条边
                else if (pointOnWhichLineIndexPairs[i] > turningIndex.Last())
                {
                    indexs.Add(turningRange.Count + 1);
                }
                // 四边形的第1,2条边
                else
                {
                    for (int j = 0; j < turningRange.Count; j++)
                    {
                        if (pointOnWhichLineIndexPairs[i] > turningRange[j][0] && pointOnWhichLineIndexPairs[i] <= turningRange[j][1] )
                        {
                            indexs.Add(j + 1);
                        }
                    }
                }
            }

            if (indexs[0] < indexs[1])
            {
                counterClockwiseCount = indexs[1] - indexs[0];
                clockwiseCount = (turningRange.Count + 2) - indexs[1] + indexs[0];
            }
            else
            {
                counterClockwiseCount = (turningRange.Count + 2) - indexs[0] + indexs[1];
                clockwiseCount = indexs[0] - indexs[1];
            }

            if (counterClockwiseCount < clockwiseCount)
            {
                return counterClockwiseCount;
            }
            else
            {
                return clockwiseCount;
            }
        }

        private int CalTCount(int interval)
        {
            int tCount;
            if (interval == 0)
            {
                tCount = 1;
            }
            else if (interval == 1)
            {
                tCount = 0;
            }
            else
            {
                tCount = 1;
            }

            return tCount;
        }

        private PlanktonMesh MakeStepLine(PlanktonMesh P,
                                           int interval,
                                           int vertex1Index,
                                           int vertex2Index,
                                           //Point3d point1,
                                           //Point3d point2,
                                           BoundarySegment bsForPoint1,
                                           BoundarySegment bsForPoint2,
                                           int point1OnWhichSegment,
                                           int point2OnWhichSegment,
                                           //List<Point3d[]> linesPassingThroughPoint1,

                                           int outerFaceIndex,

                                           BoundingBox boundingBox,
                                           DataTree<double> innerTOnEachBS,
                                           //DataTree<double> tYDT,
                                           int branchIndex,
                                           bool isLastPair,
                                           int thisLayerFirstFaceIndex_OnP)
        {
            PlanktonMesh PDeepCopy = new PlanktonMesh(P);

            List<string> debugString2 = UtilityFunctions.PrintFacesVertices(PDeepCopy);
            List<string> debugString3 = UtilityFunctions.PrintHalfedges(PDeepCopy);
            List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(PDeepCopy);

            #region 基础计算
            #region 求point1和point2所对应的BS的法线（逆时针90度）
            // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
            Vector3d vectorForBS1 = new Vector3d(bsForPoint1.Lines[point1OnWhichSegment].To - bsForPoint1.Lines[point1OnWhichSegment].From);
            Vector3d vectorForBS2 = new Vector3d(bsForPoint2.Lines[point2OnWhichSegment].To - bsForPoint2.Lines[point2OnWhichSegment].From);
            Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.Lines[point1OnWhichSegment].From, bsForPoint1.Lines[point1OnWhichSegment].To);
            Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.Lines[point2OnWhichSegment].From, bsForPoint2.Lines[point2OnWhichSegment].To);

            // 求在x轴正方向或y轴正方向的投影
            Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
            Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
            // 投影后向量的normal，为正X正Y方向
            Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
            Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
            #endregion
            #endregion

            #region 输出point1是否被绘制过
            bool isPoint1Drawn = false;
            int h = -1;
            int pairFaceIndex = -1;
            List<int> hStartAtPoint1 = PDeepCopy.Vertices.GetHalfedges(vertex1Index).ToList();
            for (int i = 0; i < hStartAtPoint1.Count; i++)
            {
                int face1 = PDeepCopy.Halfedges[hStartAtPoint1[i]].AdjacentFace;
                int face2 = PDeepCopy.Halfedges[PDeepCopy.Halfedges.GetPairHalfedge(hStartAtPoint1[i])].AdjacentFace;
                if (face1 != outerFaceIndex && face2 != outerFaceIndex)
                {
                    isPoint1Drawn = true;
                    h = hStartAtPoint1[i];
                    pairFaceIndex = face1;
                    break;
                }
            }

            #endregion

            Point3d point1 = PDeepCopy.Vertices[vertex1Index].ToPoint3d();
            Point3d point2 = PDeepCopy.Vertices[vertex2Index].ToPoint3d();

            List<Point3d> turningPoints = new List<Point3d>();

            //int branchIndex = 0;
            int itemIndex = 0;

            
            List<PlanktonXYZ> newPlanktonXYZ;
            List<List<int>> newFaceVertexOrder;

            if (interval == 0)
            {
                // interval == 0
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = PDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = PDeepCopy.Vertices[endPoint1Index].ToPoint3d();

                    #region 计算turningPoint
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = Math.Abs(endPoint1.Y - point1.Y);

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addY;
                        if (nVectorForBS1.Y > 0)
                        {
                            addY = y * t;
                        }
                        else
                        {
                            addY = -y * t;
                        }

                        turningPoints.Add(new Point3d(point1.X, point1.Y + addY, point1.Z));
                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                    }
                    else
                    {
                        // 法向量为x轴
                        double x = Math.Abs(endPoint1.X - point1.X);

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addX;
                        if (nVectorForBS1.X > 0)
                        {
                            addX = x * t;
                        }
                        else
                        {
                            addX = -x * t;
                        }

                        turningPoints.Add(new Point3d(point1.X + addX, point1.Y, point1.Z));
                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                    }
                    #endregion
                    //branchIndex++;

                    PDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;
                    PDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    if (isLastPair)
                    {
                        viNewFace.Remove(vertex2Index);

                        viNewFace.Remove(vertex1Index);

                        viFace.Add(turningPoint1Index);
                        viFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint1Index);

                        List<int> viPairFace = PDeepCopy.Faces.GetFaceVertices(pairFaceIndex).ToList();
                        int index = viPairFace.IndexOf(vertex1Index);
                        viPairFace.Insert((index + 1) % viPairFace.Count, turningPoint0Index);

                        List<int> viFirstFace = PDeepCopy.Faces.GetFaceVertices(thisLayerFirstFaceIndex_OnP).ToList();
                        int index2 = viFirstFace.IndexOf(vertex2Index);
                        viFirstFace.Insert(index2, turningPoint1Index);

                        #region 生成新的PlanktonMesh
                        GenerateNewPlanktonXYZAndNewFaceVertexOrder4(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, pairFaceIndex, thisLayerFirstFaceIndex_OnP, viFace, viNewFace, viPairFace, viFirstFace);
                        #endregion
                    }
                    else
                    {
                        viNewFace.Remove(vertex1Index);

                        viFace.Add(turningPoint1Index);
                        viFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint1Index);

                        List<int> viPairFace = PDeepCopy.Faces.GetFaceVertices(pairFaceIndex).ToList();
                        int index = viPairFace.IndexOf(vertex1Index);
                        viPairFace.Insert((index + 1) % viPairFace.Count, turningPoint0Index);
                        

                        #region 生成新的PlanktonMesh
                        GenerateNewPlanktonXYZAndNewFaceVertexOrder2(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, pairFaceIndex, viFace, viNewFace, viPairFace);
                        #endregion
                    }
                    #endregion
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = boundingBox.Max.Y - boundingBox.Min.Y;

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addY;
                        if (nVectorForBS1.Y > 0)
                        {
                            addY = y * t;
                        }
                        else
                        {
                            addY = -y * t;
                        }

                        turningPoints.Add(new Point3d(point1.X, point1.Y + addY, point1.Z));
                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                    }
                    else
                    {
                        // 法向量为x轴
                        double x = boundingBox.Max.X - boundingBox.Min.X;

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addX;
                        if (nVectorForBS1.X > 0)
                        {
                            addX = x * t;
                        }
                        else
                        {
                            addX = -x * t;
                        }

                        turningPoints.Add(new Point3d(point1.X + addX, point1.Y, point1.Z));
                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                    }

                    //branchIndex++;

                    PDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;
                    PDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint0Index);
                    viFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }
            else if (interval == 1)
            {
                // interval == 1
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = PDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = PDeepCopy.Vertices[endPoint1Index].ToPoint3d();

                    Point3d intersectPoint = GetLineIntersection(point1, nVectorForBS1, point2, nVectorForBS2);

                    #region 判断交点是否在经过point1的边的内部
                    Vector3d vector1 = new Vector3d(endPoint1 - point1);
                    Vector3d vector2 = new Vector3d(intersectPoint - point1);
                    double t = vector2.Length / vector1.Length;
                    bool isOutOfRange = false;
                    if (t > 1)
                    {
                        isOutOfRange = true;
                    }
                    else
                    {
                        isOutOfRange = false;
                    }
                    #endregion

                    PDeepCopy.Vertices.Add(intersectPoint);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viNewFace.Remove(vertex1Index);
                    if (isOutOfRange)
                    {
                        // 对于OutOfRange时的特殊处理
                        if (isLastPair)
                        {
                            viNewFace.Remove(vertex2Index);

                            viFace.Add(turningPoint0Index);
                            viFace.Add(viNewFace.Last());
                            viNewFace.Add(turningPoint0Index);

                            List<int> viFirstFace = PDeepCopy.Faces.GetFaceVertices(thisLayerFirstFaceIndex_OnP).ToList();
                            int index = viFirstFace.IndexOf(vertex2Index);
                            viFirstFace.Insert(index, turningPoint0Index);

                            #region 生成新的PlanktonMesh
                            GenerateNewPlanktonXYZAndNewFaceVertexOrder3(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, thisLayerFirstFaceIndex_OnP, viFace, viNewFace, viFirstFace);
                            #endregion
                        }
                        else
                        {
                            viFace.Add(turningPoint0Index);
                            viFace.Add(viNewFace.Last());
                            viNewFace.Add(turningPoint0Index);

                            #region 生成新的PlanktonMesh
                            GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                            #endregion
                        }
                    }
                    else
                    {
                        viFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint0Index);

                        List<int> viPairFace = PDeepCopy.Faces.GetFaceVertices(pairFaceIndex).ToList();
                        int index = viPairFace.IndexOf(vertex1Index);
                        viPairFace.Insert((index + 1) % viPairFace.Count, turningPoint0Index);

                        #region 生成新的PlanktonMesh
                        GenerateNewPlanktonXYZAndNewFaceVertexOrder2(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, pairFaceIndex, viFace, viNewFace, viPairFace);
                        #endregion
                    }
                    #endregion
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    if (nVectorForBS1.X == 0)
                    {
                        // 法线是y轴方向上的
                        Point3d turningPoint = new Point3d(point1.X, point2.Y, point1.Z);

                        turningPoints.Add(turningPoint);
                    }
                    else
                    {
                        // 法线是x轴方向上的
                        Point3d turningPoint = new Point3d(point2.X, point1.Y, point1.Z);

                        turningPoints.Add(turningPoint);
                    }

                    PDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }
            else
            {
                // interval == 2
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = PDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = PDeepCopy.Vertices[endPoint1Index].ToPoint3d();

                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = Math.Abs(endPoint1.Y - point1.Y);

                        double t2 = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addY;
                        if (nVectorForBS1.Y > 0)
                        {
                            addY = y * t2;
                        }
                        else
                        {
                            addY = -y * t2;
                        }

                        turningPoints.Add(new Point3d(point1.X, endPoint1.Y + addY, point1.Z));
                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                    }
                    else
                    {
                        // 法向量为x轴
                        double x = Math.Abs(endPoint1.X - point1.X);

                        double t2 = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addX;
                        if (nVectorForBS1.X > 0)
                        {
                            addX = x * t2;
                        }
                        else
                        {
                            addX = -x * t2;
                        }

                        turningPoints.Add(new Point3d(endPoint1.X + addX, point1.Y, point1.Z));
                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                    }

                    PDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;
                    PDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viNewFace.Remove(vertex1Index);

                    int currentEnd = viNewFace.Last();

                    //viFace.Add(viNewFace[0]);
                    viFace.Add(turningPoint1Index);
                    viFace.Add(turningPoint0Index);
                    viFace.Add(currentEnd);

                    if (isLastPair)
                    {
                        viNewFace.Remove(vertex2Index);

                        viNewFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint1Index);

                        List<int> viFirstFace = PDeepCopy.Faces.GetFaceVertices(thisLayerFirstFaceIndex_OnP).ToList();
                        int index = viFirstFace.IndexOf(vertex2Index);
                        viFirstFace.Insert(index, turningPoint1Index);

                        #region 生成新的PlanktonMesh
                        GenerateNewPlanktonXYZAndNewFaceVertexOrder3(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, thisLayerFirstFaceIndex_OnP, viFace, viNewFace, viFirstFace);
                        #endregion
                    }
                    else
                    {
                        viNewFace.Add(turningPoint0Index);
                        viNewFace.Add(turningPoint1Index);
                        #endregion

                        #region 生成新的PlanktonMesh
                        GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                        #endregion
                    }
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = boundingBox.Max.Y - boundingBox.Min.Y;

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addY;
                        if (nVectorForBS1.Y > 0)
                        {
                            addY = y * t;
                        }
                        else
                        {
                            addY = -y * t;
                        }

                        turningPoints.Add(new Point3d(point1.X, point1.Y + addY, point1.Z));
                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                    }
                    else
                    {
                        // 法向量为x轴
                        double x = boundingBox.Max.X - boundingBox.Min.X;

                        double t = innerTOnEachBS.Branch(branchIndex)[itemIndex];
                        double addX;
                        if (nVectorForBS1.X > 0)
                        {
                            addX = x * t;
                        }
                        else
                        {
                            addX = -x * t;
                        }

                        turningPoints.Add(new Point3d(point1.X + addX, point1.Y, point1.Z));
                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                    }

                    //branchIndex++;

                    PDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = PDeepCopy.Vertices.Count - 1;
                    PDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = PDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    //commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, PDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, PDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint1Index);
                    viFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint1Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }

            P = new PlanktonMesh();
            P.Vertices.AddVertices(newPlanktonXYZ);
            P.Faces.AddFaces(newFaceVertexOrder);

            return P;
        }

        private int GetViFaceIndex(int interval, int vertex1Index, int vertex2Index, PlanktonMesh originPDeepCopy)
        {
            List<int> faceList1 = originPDeepCopy.Vertices.GetVertexFaces(vertex1Index).ToList();
            List<int> faceList2 = originPDeepCopy.Vertices.GetVertexFaces(vertex2Index).ToList();
            if (faceList1.Contains(-1))
            {
                faceList1.Remove(-1);
            }
            if (faceList2.Contains(-1))
            {
                faceList2.Remove(-1);
            }

            List<int> intersect = faceList1.Intersect(faceList2).ToList();
            if (intersect.Count == 0 && intersect.Count > 1)
            {
                throw new Exception(string.Format("interval = {0}，point1没有被绘制过时，intersect.Count为{1}", interval, intersect.Count));
            }
            int viFaceIndex = faceList1.Intersect(faceList2).ToList()[0];
            return viFaceIndex;
        }

        private void GetInitialViLists(int vertex1Index, int vertex2Index, PlanktonMesh P, int adjacentFaceIndex, out List<int> viFace, out List<int> viNewFace)
        {
            List<int> viAroundOriginFace = new List<int>();
            viAroundOriginFace = P.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();

            int vertex1InViAroundOriginFaceIndex = viAroundOriginFace.IndexOf(vertex1Index);
            int vertex2InViAroundOriginFaceIndex = viAroundOriginFace.IndexOf(vertex2Index);

            viNewFace = new List<int>();
            int iter1 = vertex1InViAroundOriginFaceIndex;
            while (viAroundOriginFace[(iter1 + viAroundOriginFace.Count) % viAroundOriginFace.Count] != vertex2Index)
            {
                iter1--;
                viNewFace.Add(viAroundOriginFace[(iter1 + viAroundOriginFace.Count) % viAroundOriginFace.Count]);
            }

            viFace = new List<int>();
            iter1 = vertex1InViAroundOriginFaceIndex;
            while (viAroundOriginFace[(iter1 + viAroundOriginFace.Count) % viAroundOriginFace.Count] != vertex2Index)
            {
                iter1++;
                viFace.Add(viAroundOriginFace[(iter1 + viAroundOriginFace.Count) % viAroundOriginFace.Count]);
            }

            viNewFace.Insert(0, vertex1Index);
            viFace.Insert(0, vertex1Index);
            viNewFace.Reverse();
        }

        private void GenerateNewPlanktonXYZAndNewFaceVertexOrder1(PlanktonMesh P, 
                                                                  out List<PlanktonXYZ> newPlanktonXYZ, 
                                                                  out List<List<int>> newFaceVertexOrder, 
                                                                  int viFaceIndex, 
                                                                  List<int> viFace, 
                                                                  List<int> viNewFace)
        {
            // 转移原来originPDeepCopy的属性
            newPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < P.Vertices.Count; i++)
            {
                newPlanktonXYZ.Add(P.Vertices[i].ToXYZ());
            }
            //
            newFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                if (i == viFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFace);
                }
                else
                {
                    newFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = P.Faces.GetFaceVertices(i);
                    newFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            newFaceVertexOrder.Add(viNewFace);
        }

        private void GenerateNewPlanktonXYZAndNewFaceVertexOrder2(PlanktonMesh P, 
                                                                  out List<PlanktonXYZ> newPlanktonXYZ, 
                                                                  out List<List<int>> newFaceVertexOrder, 
                                                                  int viFaceIndex,
                                                                  int viPairFaceIndex, 
                                                                  List<int> viFace, 
                                                                  List<int> viNewFace, 
                                                                  List<int> viPairFace)
        {
            // 转移原来originPDeepCopy的属性
            newPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < P.Vertices.Count; i++)
            {
                newPlanktonXYZ.Add(P.Vertices[i].ToXYZ());
            }
            //
            newFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                if (i == viFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFace);
                }
                else if (i == viPairFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viPairFace);
                }
                else
                {
                    newFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = P.Faces.GetFaceVertices(i);
                    newFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            newFaceVertexOrder.Add(viNewFace);
        }

        private void GenerateNewPlanktonXYZAndNewFaceVertexOrder3(PlanktonMesh P,
                                                                  out List<PlanktonXYZ> newPlanktonXYZ,
                                                                  out List<List<int>> newFaceVertexOrder,
                                                                  int viFaceIndex,
                                                                  int viFirstFaceIndex,
                                                                  List<int> viFace,
                                                                  List<int> viNewFace,
                                                                  List<int> viFirstFace)
        {
            // 转移原来originPDeepCopy的属性
            newPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < P.Vertices.Count; i++)
            {
                newPlanktonXYZ.Add(P.Vertices[i].ToXYZ());
            }
            //
            newFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                if (i == viFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFace);
                }
                else if (i == viFirstFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFirstFace);
                }
                else
                {
                    newFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = P.Faces.GetFaceVertices(i);
                    newFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            newFaceVertexOrder.Add(viNewFace);
        }

        private void GenerateNewPlanktonXYZAndNewFaceVertexOrder4(PlanktonMesh P,
                                                                  out List<PlanktonXYZ> newPlanktonXYZ,
                                                                  out List<List<int>> newFaceVertexOrder,
                                                                  int viFaceIndex,
                                                                  int viPairFaceIndex,
                                                                  int viFirstFaceIndex,
                                                                  List<int> viFace,
                                                                  List<int> viNewFace,
                                                                  List<int> viPairFace,
                                                                  List<int> viFirstFace)
        {
            // 转移原来originPDeepCopy的属性
            newPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < P.Vertices.Count; i++)
            {
                newPlanktonXYZ.Add(P.Vertices[i].ToXYZ());
            }
            //
            newFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                if (i == viFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFace);
                }
                else if (i == viPairFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viPairFace);
                }
                else if (i == viFirstFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFirstFace);
                }
                else
                {
                    newFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = P.Faces.GetFaceVertices(i);
                    newFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            newFaceVertexOrder.Add(viNewFace);
        }

        private List<T> LeftShift<T>(List<T> originList, int shift)
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

        private int WhichQuadrant(Point3d point1, Point3d point2)
        {
            double dx = point2.X - point1.X;
            double dy = point2.Y - point1.Y;

            if (dx == 0)
            {
                // y轴
                if (dy > 0)
                {
                    // y轴正方向，并入第一象限的情况计算
                    return 1;
                }
                else
                {
                    // y轴反方向，并入第三象限的情况计算
                    return 3;
                }
            }
            if (dy == 0)
            {
                // x轴
                if (dx > 0)
                {
                    // x轴正方向，并入第四象限计算
                    return 4;
                }
                else
                {
                    // x轴反方向，并入第二象限计算
                    return 2;
                }
            }

            if (dy / dx > 0 && dx > 0)
            {
                return 1;
            }
            else if (dy / dx > 0 && dx < 0)
            {
                return 3;
            }
            else if (dy / dx < 0 && dx < 0)
            {
                return 2;
            }
            else
            {
                return 4;
            }
        }

        private Vector3d CalProjectVector(Point3d point1, Point3d point2)
        {
            Vector3d projectVector;

            // 要考虑原来的向量是属于哪个象限的问题
            int quadrant = WhichQuadrant(point1, point2);

            if (point2.X - point1.X != 0)
            {
                double k = (point2.Y - point1.Y) / (point2.X - point1.X);
                if (k > 0 && k <= 1)
                {
                    // 斜率小于45度时
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k > 1)
                {
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                }
                else if (k == 0)
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k >= -1 && k < 0)
                {
                    // 斜率小于45度时
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                }
            }
            else
            {
                // 即y轴正方向
                if (quadrant == 1)
                {
                    projectVector = new Vector3d(0, 1, 0);
                }
                else//quadrant == 3
                {
                    projectVector = new Vector3d(0, -1, 0);
                }
            }

            return projectVector;
        }

        /// <summary>
        /// 计算向量逆时针旋转九十度的单位法向量
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        private Vector3d Normal(Vector3d vector)
        {
            double L = Math.Sqrt(vector * vector);
            return new Vector3d(-vector.Y / L, vector.X / L, vector.Z);
        }

        /// <summary>
        /// 找到两个射线的交点
        /// </summary>
        /// <param name="p"></param>
        /// <param name="v"></param>
        /// <param name="q"></param>
        /// <param name="w"></param>
        /// <returns></returns>
        private Point3d GetLineIntersection(Point3d p, Vector3d v, Point3d q, Vector3d w)
        {
            Vector3d u = p - q;
            if (CrossProduct(v, w) == 0)
            {
                return Point3d.Unset;
            }
            else
            {
                double t = CrossProduct(w, u) / CrossProduct(v, w);
                return p + v * t;
            }
        }

        /// <summary>
        /// 不考虑Z值的叉积
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        private double CrossProduct(Vector3d a, Vector3d b)
        {
            return a.X * b.Y - a.Y * b.X;
        }

        private Point3d GetCenterOfGravityPoint(List<Point3d> points)
        {
            // 多边形面积
            double area = 0;
            // 重心的x,y
            double gx = 0;
            double gy = 0;
            for (int i = 1; i <= points.Count; i++)
            {
                double ix = points[(i % points.Count)].X;
                double iy = points[(i % points.Count)].Y;
                double nextx = points[i - 1].X;
                double nexty = points[i - 1].Y;
                double temp = (ix * nexty - iy * nextx) / 2;
                area += temp;
                gx += temp * (ix + nextx) / 3;
                gy += temp * (iy + nexty) / 3;
            }

            gx = gx / area;
            gy = gy / area;

            return new Point3d(gx, gy, points[0].Z);
        }


        public List<Line> GraphEdgeList(DataTree<int> graph, List<Point3d> graphVertices)
        {
            List<Line> list = new List<Line>();
            for (int i = 0; i < graph.BranchCount; i++)
            {
                for (int j = 0; j < graph.Branch(i).Count; j++)
                {
                    Line item = new Line(graphVertices[i], graphVertices[graph.Branch(i)[j]]);
                    list.Add(item);
                }
            }

            return list;
        }

        public List<Polyline> GMeshFaceBoundaries(List<Point3d> sortedBoundaryPoints, List<Line> edges, Plane basePlane)
        {
            List<Curve> splitLines = new List<Curve>();
            foreach (Line line in edges)
            {
                splitLines.Add(line.ToNurbsCurve());
            }

            Polyline polyline = SortedPointsToPolyline(sortedBoundaryPoints, basePlane);
            Curve curve = polyline.ToNurbsCurve();

            Brep[] PlanarBrepArray = Brep.CreatePlanarBreps(curve, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            Brep planarBreps = PlanarBrepArray[0];
            BrepFace planarBrepFace = planarBreps.Faces[0];
            Brep splittedPlanarBrepFaces = planarBrepFace.Split(splitLines, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

            List<Polyline> allConvexPolyline = new List<Polyline>();
            foreach (BrepFace brepFace in splittedPlanarBrepFaces.Faces)
            {
                Brep singleSplitFaceBrep = brepFace.DuplicateFace(true);
                Point3d[] singleSplitFaceVertices = singleSplitFaceBrep.DuplicateVertices();
                Curve splitFaceConvexCure = Curve.JoinCurves(singleSplitFaceBrep.Faces[0].DuplicateFace(true).DuplicateEdgeCurves())[0];
                Polyline splitFaceConvexPolyline = null;
                splitFaceConvexCure.TryGetPolyline(out splitFaceConvexPolyline);
                allConvexPolyline.Add(splitFaceConvexPolyline);
            }

            return allConvexPolyline;
        }

        /// <summary>
        /// 将排好序的点连接成多段线（NurbsCurve类型）
        /// </summary>
        /// <param name="sortedBoundaryPoints"></param>
        /// <param name="basePlane"></param>
        /// <returns></returns>
        public Polyline SortedPointsToPolyline(List<Point3d> sortedBoundaryPoints, Plane basePlane)
        {
            // 把startPoint放到最后作为endPoint，这样polyline才会闭合
            sortedBoundaryPoints.Add(sortedBoundaryPoints[0]);

            Polyline convex = new Polyline(sortedBoundaryPoints);
            return convex;
        }


        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);

            for (int i = 0; i < DualPolylines.Count; i++)
            {
                args.Display.DrawPolyline(DualPolylines[i], Color.DarkOrange, Thickness);
            }

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            if (this.CompWorkMode == ShowMode.ShowTopology)
            {
                for (int i = 0; i < DottedCurve.Count; i++)
                {
                    args.Display.DrawCurve(DottedCurve[i], Color.DarkGreen, Thickness);
                }

                args.Display.DrawLines(GraphEdges, Color.BlueViolet, Thickness);

                for (int i = 0; i < InnerNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(InnerNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
                for (int i = 0; i < OuterNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(OuterNodeTextDots[i], Color.Gray, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);

            for (int i = 0; i < DualPolylines.Count; i++)
            {
                args.Display.DrawPolyline(DualPolylines[i], Color.DarkOrange, Thickness);
            }

            for (int i = 0; i < DualVertexTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(DualVertexTextDots[i], Color.DarkOrange, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            if (this.CompWorkMode == ShowMode.ShowTopology)
            {
                args.Display.DrawMeshShaded(MeshForVisualize, new Rhino.Display.DisplayMaterial(Color.White, 0));

                for (int i = 0; i < DottedCurve.Count; i++)
                {
                    args.Display.DrawCurve(DottedCurve[i], Color.DarkGreen, Thickness);
                }
                
                args.Display.DrawLines(GraphEdges, Color.BlueViolet, Thickness);

                for (int i = 0; i < InnerNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(InnerNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
                for (int i = 0; i < OuterNodeTextDots.Count; i++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(OuterNodeTextDots[i], Color.Gray, Color.White, Color.White);
                    args.Display.EnableDepthTesting(true);
                }
            }
        }

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
            get { return new Guid("1dbd63ea-90a0-4055-8d77-c38df937343e"); }
        }

        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new ShowAttribute(this);
        }

        public class ShowAttribute : GH_ComponentAttributes
        {
            public ShowAttribute(GhcSiteDivision3 component) : base(component) { }

            protected override void Layout()
            {
                base.Layout();
                /* 先执行base.Layout()，可以按GH电池默认方式计算电池的出/入口需要的高度，我们在下面基于这个高度进行更改 */
                Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 20.0f);
            }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);/* 执行基本的电池渲染 */

                /* 额外的电池渲染，仅在“Objects”这个渲染轨道绘制 */
                if (channel == GH_CanvasChannel.Objects)
                {
                    RectangleF buttonRect = /* 按钮的位置 */ new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);

                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect.Inflate(-2.0f, -2.0f);

                    using (GH_Capsule capsule = GH_Capsule.CreateCapsule(buttonRect, GH_Palette.Black))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(((GhcSiteDivision3)Owner).CompWorkMode.ToString(),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.White,
                                        buttonRect,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });
                }
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                RectangleF buttonRect = /* -重新计算按钮的区域大小- */ new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);

                if (e.Button == MouseButtons.Left && buttonRect.Contains(e.CanvasLocation))
                {
                    GhcSiteDivision3 comp = (GhcSiteDivision3)Owner; /* 通过Owner属性来获得电池本身 */

                    /* 依照电池当前工作状态来改变电池 */
                    if (comp.CompWorkMode == GhcSiteDivision3.ShowMode.NotShowTopology)
                        comp.CompWorkMode = GhcSiteDivision3.ShowMode.ShowTopology;
                    else
                        comp.CompWorkMode = GhcSiteDivision3.ShowMode.NotShowTopology;

                    /* 改变完电池后，重启计算 */
                    comp.ExpireSolution(true);

                    /* 结束鼠标事件处理，通知GH已经处理完毕 */
                    return GH_ObjectResponse.Handled;
                }


                return GH_ObjectResponse.Ignore;
            }
        }
    }

    
}