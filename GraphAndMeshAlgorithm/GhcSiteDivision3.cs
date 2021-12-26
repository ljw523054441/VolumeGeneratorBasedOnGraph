using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

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
            DualVertices = new List<Point3d>();
            DualPolylines = new List<List<Point3d>>();
            DualVertexTextDots = new List<TextDot>();
        }

        private int Thickness;

        private List<Point3d> DualVertices;
        private List<List<Point3d>> DualPolylines;
        private List<TextDot> DualVertexTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);

            // pManager.AddGenericParameter("SubDualGraphWithHM", "SubDGHM", "对偶图中心的子图", GH_ParamAccess.list);

            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);
            // pManager.AddIntegerParameter("BSIndexContainVolumeJunctions", "BSI", "包含边界分裂点的BoundarySegment序号", GH_ParamAccess.list);

            pManager.AddIntegerParameter("VerticesIndexForEachBS", "VIFBS", "每个BS上点对应的对偶图中的index", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("indexOnEachBS", "IOBS", "每个BS上的VolumeJunctionIndex", GH_ParamAccess.tree);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);
            // pManager.AddGenericParameter("PairVolumeJunctionsIndex", "PJI", "作为分界线的一对分裂点的序号", GH_ParamAccess.list);
            // pManager.AddBooleanParameter("NeedToConnect", "NTC", "的这一对分裂点是否作为分界线", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectVolumeJunctionsIndex", "N_VJI", "需要构成分界线的两个volumeJunction的Index", GH_ParamAccess.list);
            pManager.AddGenericParameter("VerticeIndexOnNeedToConnectDividingLine", "VI_OnDL", "在分界线上的点的index列表", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectVolumeJunctionsIndexCorrespondingFaceIndex", "N_VJCFI", "需要构成分界线的两个volumeJunction所属于的Face的Index", GH_ParamAccess.list);

            pManager.AddGenericParameter("NeedToConnectBSIndex", "N_BSI", "这两个volumeJunction所对应的bsIndex", GH_ParamAccess.list);

            pManager.AddIntegerParameter("tXCount", "tXCount", "tX的数量", GH_ParamAccess.item);
            pManager.AddIntegerParameter("tYCount", "tYCount", "tY的数量", GH_ParamAccess.item);

            pManager.AddNumberParameter("tList", "ts", "每个边界分裂点对应的t值", GH_ParamAccess.list);

            pManager.AddNumberParameter("tXList", "tXs", "", GH_ParamAccess.list);
            pManager.AddNumberParameter("tYList", "tYs", "", GH_ParamAccess.list);
            pManager[9].Optional = true;
            pManager[10].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("tXDT", "tXDT", "", GH_ParamAccess.tree);
            pManager.AddNumberParameter("tYDT", "tYDT", "", GH_ParamAccess.tree);
            
            pManager.AddCurveParameter("DeBug Polyline", "DP", "", GH_ParamAccess.item);
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
            //List<int> bSIndexContainVolumeJunctions = new List<int>();

            GH_Structure<GH_Integer> gh_structure1 = new GH_Structure<GH_Integer>();

            GH_Structure<GH_Integer> gh_structure2 = new GH_Structure<GH_Integer>();

            List<string> volumeJunctionsTexts = new List<string>();
            //List<int[]> pairVolumeJunctionsIndex = new List<int[]>();


            //List<bool> needToConnect = new List<bool>();
            //List<int[]> vertexIndexPairs = new List<int[]>();
            //List<int[]> vertexIndexOnDividingLine = new List<int[]>();

            List<List<int[]>> vertexIndexPairs = new List<List<int[]>>();
            List<List<int[]>> vertexIndexsOnDividingLine = new List<List<int[]>>();

            List<List<int>> vertexIndexPairsCorrespondingFaceIndex = new List<List<int>>();

            List<int[]> bsIndexPairs = new List<int[]>();

            int tXCount = 0;
            int tYCount = 0;

            List<double> tList = new List<double>();
            List<double> tXList = new List<double>();
            List<double> tYList = new List<double>();




            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetDataList("SortedBoundarySegments", sortedBoundarySegments)
                //&& DA.GetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions)
                && DA.GetDataTree(2, out gh_structure1)
                && DA.GetDataTree(3, out gh_structure2)

                && DA.GetDataList("VolumeJunctionsTexts", volumeJunctionsTexts)

                // && DA.GetDataList("PairVolumeJunctionsIndex", pairVolumeJunctionsIndex)
                //&& DA.GetDataList("NeedToConnect", needToConnect)
                && DA.GetDataList("NeedToConnectVolumeJunctionsIndex", vertexIndexPairs)
                && DA.GetDataList("VerticeIndexOnNeedToConnectDividingLine", vertexIndexsOnDividingLine)

                && DA.GetDataList("NeedToConnectVolumeJunctionsIndexCorrespondingFaceIndex", vertexIndexPairsCorrespondingFaceIndex)

                && DA.GetDataList("NeedToConnectBSIndex", bsIndexPairs)

                && DA.GetData("tXCount", ref tXCount)
                && DA.GetData("tYCount", ref tYCount)

                && DA.GetDataList("tList", tList))
            {
                DA.GetDataList("tXList", tXList);
                DA.GetDataList("tYList", tYList);

                #region 计算boundingBox
                List<Point3d> cornerPoints = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    cornerPoints.Add(sortedBoundarySegments[i].From);
                }

                BoundingBox boundingBox = new BoundingBox(cornerPoints);
                #endregion

                List<DataTree<int>> allLayerIndexOnEachBS = new List<DataTree<int>>();

                DataTree<int> verticesIndexForEachBSDT = new DataTree<int>();
                UtilityFunctions.GH_StructureToDataTree_Int(gh_structure1, ref verticesIndexForEachBSDT);

                #region 构造每一个volumeJunction点

                #region DataTree<int> indexOnEachBS
                DataTree<int> indexOnEachBS = new DataTree<int>();
                UtilityFunctions.GH_StructureToDataTree_Int(gh_structure2, ref indexOnEachBS);
                #endregion







                #region DataTree<double> tOnEachBS
                // 每个包含volumeJunction的bs为了确定volumeJunction位置所需要的t的datatree
                if (indexOnEachBS.DataCount > tList.Count)
                {
                    throw new Exception("tList的数量不足");
                }

                DataTree<double> tOnEachBS = new DataTree<double>();
                int index = 0;
                for (int i = 0; i < indexOnEachBS.BranchCount; i++)
                {
                    tOnEachBS.EnsurePath(i);
                    for (int j = 0; j < indexOnEachBS.Branch(i).Count; j++)
                    {
                        tOnEachBS.Branch(i).Add(tList[index]);
                        index++;
                    }
                }


                for (int i = 0; i < tOnEachBS.BranchCount; i++)
                {
                    List<double> ts = tOnEachBS.Branch(i);
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
                for (int i = 0; i < tOnEachBS.BranchCount; i++)
                {
                    pointOnEachBS.EnsurePath(i);
                    for (int j = 0; j < tOnEachBS.Branch(i).Count; j++)
                    {
                        pointOnEachBS.Branch(i).Add(sortedBoundarySegments[i].Line.PointAt(tOnEachBS.Branch(i)[j]));
                    }
                }

                #endregion

                #endregion

                #region 计算tXDT和tYDT的数量
                if (tXCount != tXList.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "tX的数量与预期的不符，请增加或减少tX的数量");
                }
                if (tYCount != tYList.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "tY的数量与预期的不符，请增加或减少tY的数量");
                }

                DataTree<double> tXDT = new DataTree<double>();
                DataTree<double> tYDT = new DataTree<double>();

                int txCount = 0;
                int tyCount = 0;
                for (int i = 0; i < vertexIndexPairs.Count; i++)
                {
                    int interval = CalInterval(bsIndexPairs[i], 4);
                    int[] tCount = CalTCount(interval,
                                             sortedBoundarySegments[bsIndexPairs[i][0]],
                                             sortedBoundarySegments[bsIndexPairs[i][1]]);


                    tXDT.EnsurePath(i);
                    tXDT.Branch(i).AddRange(tXList.Skip(txCount).Take(tCount[0]));
                    tYDT.EnsurePath(i);
                    tYDT.Branch(i).AddRange(tYList.Skip(tyCount).Take(tCount[1]));

                    txCount += tCount[0];
                    tyCount += tCount[1];
                }
                DA.SetDataTree(0, tXDT);
                DA.SetDataTree(1, tYDT);
                #endregion

                #region 构造最开始的初始originP
                List<PlanktonXYZ> planktonXYZs = new List<PlanktonXYZ>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    Point3d cornerPoint = sortedBoundarySegments[i].From;
                    PlanktonXYZ cornerXYZ = new PlanktonXYZ((float)cornerPoint.X, (float)cornerPoint.Y, (float)cornerPoint.Z);
                    planktonXYZs.Add(cornerXYZ);

                    if (indexOnEachBS.Branch(i).Count > 0)
                    {
                        for (int j = 0; j < indexOnEachBS.Branch(i).Count; j++)
                        {
                            Point3d point = pointOnEachBS.Branch(i)[j];
                            PlanktonXYZ planktonXYZ = new PlanktonXYZ((float)point.X, (float)point.Y, (float)point.Z);
                            planktonXYZs.Add(planktonXYZ);
                        }
                    }
                }

                #region 获取dualPlanktonMesh边界上的点的index，即初始最外圈的一圈点的列表order
                List<int> order = new List<int>();

                for (int i = 0; i < verticesIndexForEachBSDT.BranchCount; i++)
                {
                    List<int> verticesIndex = verticesIndexForEachBSDT.Branch(i);
                    verticesIndex.RemoveAt(verticesIndex.Count - 1);
                    order.AddRange(verticesIndex);
                }
                #endregion

                // order 与 newVertexIndexOrder之间存在映射关系
                List<int> newVertexIndexOrder = new List<int>();
                for (int i = 0; i < planktonXYZs.Count; i++)
                {
                    newVertexIndexOrder.Add(i);
                }

                List<List<int>> faceVertexOrder = new List<List<int>>();
                faceVertexOrder.Add(new List<int>());
                faceVertexOrder[0].AddRange(newVertexIndexOrder);

                PlanktonMesh divisionP = new PlanktonMesh();
                divisionP.Vertices.AddVertices(planktonXYZs);
                divisionP.Faces.AddFaces(faceVertexOrder);
                #endregion
                PlanktonMesh P = new PlanktonMesh(divisionP);


                #region 重新找到对应的P中的vertexIndex
                for (int i = 0; i < vertexIndexPairs.Count; i++)
                {
                    // 对每一层的vertexIndexPair
                    List<int[]> vertexIndexPairsInDivisionP = new List<int[]>();
                    for (int j = 0; j < vertexIndexPairs[i].Count; j++)
                    {
                        int[] vertexIndexPairInDivisionP = new int[2] { 0, 0 };
                        for (int k = 0; k < vertexIndexPairs[i][j].Length; k++)
                        {
                            vertexIndexPairInDivisionP[k] = order.IndexOf(vertexIndexPairs[i][j][k]);
                        }
                        vertexIndexPairsInDivisionP.Add(vertexIndexPairInDivisionP);
                    }

                    #region 控制循环次数
                    int iterCount = 0;
                    bool needToCalNextLayer = false;
                    if (i == vertexIndexPairs.Count - 1)
                    {
                        // 注意这里 vertexIndexPairsInDivisionP.Count - 1，是为了不计算最后一个pair
                        iterCount = vertexIndexPairsInDivisionP.Count - 1;
                        needToCalNextLayer = false;
                    }
                    else
                    {
                        iterCount = vertexIndexPairsInDivisionP.Count;
                        needToCalNextLayer = true;
                    }
                    #endregion

                    List<string> debugString1 = UtilityFunctions.PrintFacesVertices(P);
                    #region 划分
                    for (int j = 0; j < vertexIndexPairsInDivisionP.Count; j++)
                    {
                        int interval = CalInterval(bsIndexPairs[j], 4);
                        P = MakeStepLine(P,
                                         interval,
                                         vertexIndexPairsInDivisionP[j][0],
                                         vertexIndexPairsInDivisionP[j][1],
                                         sortedBoundarySegments[bsIndexPairs[j][0]],
                                         sortedBoundarySegments[bsIndexPairs[j][1]],
                                         boundingBox,
                                         tXDT,
                                         tYDT,
                                         j);

                        List<string> debugString2 = UtilityFunctions.PrintFacesVertices(P);
                        List<string> debugString3 = UtilityFunctions.PrintHalfedges(P);
                        List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(P);
                    }

                    #region 重新计算order列表

                    #region 重新找到对应的P中的vertexIndex

                    #endregion
                    if (needToCalNextLayer)
                    {
                        
                        for (int j = 0; j < vertexIndexsOnDividingLine[i].Count; j++)
                        {
                            int[] vertexIndexOnDivingLine = vertexIndexsOnDividingLine[i][j];
                            for (int k = 0; k < vertexIndexsOnDividingLine[i][j].Length; k++)
                            {
                                vertexIndexOnDivingLine[k] = order.IndexOf(vertexIndexsOnDividingLine[i][j][k]);
                            }
                        }
                    }
                    
                    
                    #endregion
                    #endregion
                }


                // vertexIndexPairs中是每个vertex点在Dual图中的Index，我们需要把它们转换成在divisionP中的Index
                //List<int[]> vertexIndexPairsInDivisionP = new List<int[]>();
                //for (int i = 0; i < vertexIndexPairs.Count; i++)
                //{
                //    int[] vertexIndexPairInDivisionP = new int[2] { 0, 0 };
                //    for (int j = 0; j < vertexIndexPairs[i].Length; j++)
                //    {
                //        vertexIndexPairInDivisionP[j] = order.IndexOf(vertexIndexPairs[i][j]);
                //    }
                //    vertexIndexPairsInDivisionP.Add(vertexIndexPairInDivisionP);
                //}
                #endregion

                #region 划分
                // 注意这里 vertexIndexPairsInDivisionP.Count - 1，是为了不计算最后一个pair
                //for (int i = 0; i < vertexIndexPairsInDivisionP.Count; i++)
                //{
                //    int interval = CalInterval(bsIndexPairs[i], 4);
                //    P = MakeStepLine(P,
                //                     interval,
                //                     vertexIndexPairsInDivisionP[i][0],
                //                     vertexIndexPairsInDivisionP[i][1],
                //                     sortedBoundarySegments[bsIndexPairs[i][0]],
                //                     sortedBoundarySegments[bsIndexPairs[i][1]],
                //                     boundingBox,
                //                     tXDT,
                //                     tYDT,
                //                     i);

                //    List<string> debugString2 = UtilityFunctions.PrintFacesVertices(P);
                //    List<string> debugString3 = UtilityFunctions.PrintHalfedges(P);
                //    List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(P);
                //}
                #endregion

                #region 可视化部分
                PlanktonMesh MeshForVisualize = new PlanktonMesh();

                MeshForVisualize = P;
                #region clear
                // 对偶图
                DualVertices.Clear();
                DualPolylines.Clear();
                DualVertexTextDots.Clear();
                #endregion

                #region 对偶图
                DualVertices = PlanktonGh.RhinoSupport.GetPositions(MeshForVisualize).ToList();

                Polyline[] dualPolylines = PlanktonGh.RhinoSupport.ToPolylines(MeshForVisualize);
                for (int i = 0; i < dualPolylines.Length; i++)
                {
                    DualPolylines.Add(new List<Point3d>());
                    DualPolylines[i].AddRange(dualPolylines[i]);
                }

                for (int i = 0; i < MeshForVisualize.Vertices.Count; i++)
                {
                    string arg;
                    if (i<order.Count)
                    {
                        arg = order[i].ToString();
                    }
                    else
                    {
                        arg = string.Join<int>(";", MeshForVisualize.Vertices.GetVertexFaces(i));
                    }
                    
                    TextDot textDot = new TextDot(string.Format("{0} | {1}", i, arg), DualVertices[i]);
                    DualVertexTextDots.Add(textDot);
                }
                #endregion
                #endregion
            }
        }

        private PlanktonMesh GenerateOriginP(PlanktonMesh originP,
                                             DataTree<int> indexOnEachBS,
                                             List<double> tList,
                                             List<BoundarySegment> sortedBoundarySegments)
        {
            #region DataTree<double> tOnEachBS
            // 每个包含volumeJunction的bs为了确定volumeJunction位置所需要的t的datatree
            if (indexOnEachBS.DataCount > tList.Count)
            {
                throw new Exception("tList的数量不足");
            }

            DataTree<double> tOnEachBS = new DataTree<double>();
            int index = 0;
            for (int i = 0; i < indexOnEachBS.BranchCount; i++)
            {
                tOnEachBS.EnsurePath(i);
                for (int j = 0; j < indexOnEachBS.Branch(i).Count; j++)
                {
                    tOnEachBS.Branch(i).Add(tList[index]);
                    index++;
                }
            }


            for (int i = 0; i < tOnEachBS.BranchCount; i++)
            {
                List<double> ts = tOnEachBS.Branch(i);
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
            DataTree<Point3d> pointOnEachBS = new DataTree<Point3d>();
            for (int i = 0; i < tOnEachBS.BranchCount; i++)
            {
                pointOnEachBS.EnsurePath(i);
                for (int j = 0; j < tOnEachBS.Branch(i).Count; j++)
                {
                    double t = tOnEachBS.Branch(i)[j];
                    if (sortedBoundarySegments[i].TurningTs != null)
                    {
                        int iter = 0;
                        for (int k = 0; k < sortedBoundarySegments[i].TurningTs.Count - 1; k++)
                        {
                            if (t > sortedBoundarySegments[i].TurningTs[k] && t < sortedBoundarySegments[i].TurningTs[k + 1])
                            {
                                iter = k;
                                break;
                            }
                            else
                            {
                                continue;
                            }
                        }

                        // 注意这里一定要添加PointOnWhichSegments信息，即这个点是在BS的哪一段上
                        sortedBoundarySegments[i].PointOnWhichSegments = iter;

                        pointOnEachBS.Branch(i).Add(sortedBoundarySegments[i].Polyline.PointAt(t));
                    }
                    else
                    {
                        pointOnEachBS.Branch(i).Add(sortedBoundarySegments[i].Lines[0].PointAt(t));
                    }
                }
            }
            #endregion

            #region 修正OriginP

            #region 拿到P中已有的PlanktonXYZ和vertexIndex
            List<PlanktonXYZ> originPPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < originP.Vertices.Count; i++)
            {
                originPPlanktonXYZ.Add(originP.Vertices[i].ToXYZ());
            }

            List<List<int>> originPFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < originP.Faces.Count; i++)
            {
                originPFaceVertexOrder.Add(new List<int>());
                originPFaceVertexOrder[i].AddRange(originP.Faces.GetFaceVertices(i));
            }
            #endregion

            #region 需要新增的PlanktonXYZ和vertexIndex
            List<Point3d> additionalPoint3ds = new List<Point3d>();
            for (int i = 0; i < sortedBoundarySegments.Count; i++)
            {
                for (int j = 0; j < sortedBoundarySegments[i].Lines.Count; j++)
                {
                    additionalPoint3ds.Add(sortedBoundarySegments[i].Lines[j].From);
                    if (sortedBoundarySegments)
                    {

                    }
                }
            }



            List<PlanktonXYZ> additionalPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < length; i++)
            {

            }
            #endregion

            #endregion
        }

        private PlanktonMesh Divide(PlanktonMesh P,
                                    PlanktonMesh Dual,
                                    List<int> onP_vertexIndexListFromFirstVertexToDraw,
                                    List<int[]> onD_vertexIndexPairs,
                                    List<int[]> onD_vertexOnDividingLine,

                                    List<int> onD_vertexIndexPairsCorrespondingFaceIndex,

                                    List<int[]> outer_bsIndexPairs,
                                    List<BoundarySegment> outer_bs,

                                    bool needToCalNextLayer,
                                    BoundingBox currentBoundingBox,
                                    DataTree<double> tXDT,
                                    DataTree<double> tYDT)
        {
            // 注意P的PlanktonXYZ的顺序，是从FirstVertexToDraw开始的
            List<int[]> onP_vertexIndexPairs = new List<int[]>();
            for (int i = 0; i < onD_vertexIndexPairs.Count; i++)
            {
                int[] onP_vertexIndexPair = onD_vertexIndexPairs[i];
                for (int j = 0; j < onD_vertexIndexPairs[i].Length; j++)
                {
                    onP_vertexIndexPair[j] = onP_vertexIndexListFromFirstVertexToDraw.IndexOf(onP_vertexIndexPairs[i][j]);
                }
                onP_vertexIndexPairs.Add(onP_vertexIndexPair);
            }

            // 控制循环次数
            int drawCount = 0;
            if (needToCalNextLayer)
            {
                drawCount = onP_vertexIndexPairs.Count;
            }
            else
            {
                drawCount = onP_vertexIndexPairs.Count - 1;
            }

            #region 划分
            List<int> commonVertexIndexList = new List<int>();
            for (int i = 0; i < drawCount; i++)
            {
                int interval = CalInterval(outer_bsIndexPairs[i], outer_bs.Count);
                int commonVertexIndex;
                P = MakeStepLine(P,
                                 interval,
                                 onP_vertexIndexPairs[i][0],
                                 onP_vertexIndexPairs[i][1],
                                 outer_bs[outer_bsIndexPairs[i][0]],
                                 outer_bs[outer_bsIndexPairs[i][1]],
                                 currentBoundingBox,
                                 tXDT,
                                 tYDT,
                                 i,
                                 out commonVertexIndex);

                commonVertexIndexList.Add(commonVertexIndex);
            }


            if (needToCalNextLayer)
            {
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

                // 每个边缘面和中心面的公共边的列表
                List<List<int>> onP_publicHsList = new List<List<int>>();
                for (int i = 0; i < onP_vertexIndexPairsCorrespondingFaceIndex.Count; i++)
                {
                    List<int> publicHs = FindPublicHalfEdges(P, P.Faces.Count - 1, onP_vertexIndexPairsCorrespondingFaceIndex[i]);
                    onP_publicHsList.Add(publicHs);
                }

                // 获得每个边缘面的公共边所对应的一对start和end点的Index（这个Index是新绘制的P中的Index，不是原来dual图中的Index，需要进行映射）
                List<int[]> onP_publicHStartAndEndVertexIndexs = new List<int[]>();
                for (int i = 0; i < onP_publicHsList.Count; i++)
                {
                    if (onP_publicHsList[i].Count > 1)
                    {
                        int start = P.Halfedges[onP_publicHsList[i][0]].StartVertex;
                        int end = P.Halfedges.EndVertex(onP_publicHsList[i].Last());
                        int[] pair = new int[2] { start, end };
                        onP_publicHStartAndEndVertexIndexs.Add(pair);
                    }
                    else
                    {
                        int start = P.Halfedges[onP_publicHsList[i][0]].StartVertex;
                        int end = P.Halfedges.EndVertex(onP_publicHsList[i][0]);
                        int[] pair = new int[2] { start, end };
                        onP_publicHStartAndEndVertexIndexs.Add(pair);
                    }
                }

                // 将上面得到的一对Start和End点的Index进行映射
                List<int> onP_vertexOnInnerBoundary = new List<int>();
                List<int> onD_vertexOnInnerBoundary = new List<int>();
                for (int i = 0; i < onP_publicHStartAndEndVertexIndexs.Count; i++)
                {
                    // 添加每个的EndVertex
                    onP_vertexOnInnerBoundary.Add(onP_publicHStartAndEndVertexIndexs[i].Last());
                }
                for (int i = 0; i < onD_vertexOnDividingLine.Count; i++)
                {
                    // 添加每个的EndVertex
                    onD_vertexOnInnerBoundary.Add(onD_vertexOnDividingLine[i].Last());
                }

                


                List<List<BoundarySegment>> newBSs = new List<List<BoundarySegment>>();

            }
            #endregion
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
        private int CalInterval(int[] bsIndexPair, int bsCount)
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

        private int[] CalTCount(int interval,
                                BoundarySegment bsForPoint1,
                                BoundarySegment bsForPoint2)
        {
            #region 基础计算
            #region 求point1和point2所对应的BS的法线（逆时针90度）
            // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
            Vector3d vectorForBS1 = new Vector3d(bsForPoint1.To - bsForPoint1.From);
            Vector3d vectorForBS2 = new Vector3d(bsForPoint2.To - bsForPoint2.From);
            Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.From, bsForPoint1.To);
            Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.From, bsForPoint2.To);

            // 求在x轴正方向或y轴正方向的投影
            Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
            Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
            // 投影后向量的normal，为正X正Y方向
            Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
            Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
            #endregion
            #endregion

            int[] tCount;
            if (interval == 0)
            {
                // interval == 0
                if (nVectorForBS1.X == 0)
                {
                    // 法向量为y轴
                    tCount = new int[2] { 0, 1 };
                }
                else
                {
                    // 法向量为x轴
                    tCount = new int[2] { 1, 0 };
                }
            }
            else if (interval == 1)
            {
                // interval == 1
                tCount = new int[2] { 0, 0 };
            }
            else
            {
                // interval == 2
                if (nVectorForBS1.X == 0)
                {
                    // 法向量为y轴
                    tCount = new int[2] { 0, 1 };
                }
                else
                {
                    // 法向量为x轴
                    tCount = new int[2] { 1, 0 };
                }
            }

            return tCount;
        }

        private PlanktonMesh MakeStepLine(PlanktonMesh originP,
                                           int interval,
                                           int vertex1Index,
                                           int vertex2Index,
                                           BoundarySegment bsForPoint1,
                                           BoundarySegment bsForPoint2,
                                           //List<Point3d[]> linesPassingThroughPoint1,
                                           BoundingBox boundingBox,
                                           DataTree<double> tXDT,
                                           DataTree<double> tYDT,
                                           int branchIndex,
                                           out int commonVertexIndex)
                                           //out List<Line> lineList,
                                           //out PlanktonMesh P)
        {
            PlanktonMesh originPDeepCopy = new PlanktonMesh(originP);

            List<string> debugString2 = UtilityFunctions.PrintFacesVertices(originPDeepCopy);
            List<string> debugString3 = UtilityFunctions.PrintHalfedges(originPDeepCopy);
            List<string> debugString4 = UtilityFunctions.PrintHalfedgeStartAndEnd(originPDeepCopy);

            #region 基础计算
            #region 求point1和point2所对应的BS的法线（逆时针90度）
            // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
            Vector3d vectorForBS1 = new Vector3d(bsForPoint1.To - bsForPoint1.From);
            Vector3d vectorForBS2 = new Vector3d(bsForPoint2.To - bsForPoint2.From);
            Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.From, bsForPoint1.To);
            Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.From, bsForPoint2.To);

            // 求在x轴正方向或y轴正方向的投影
            Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
            Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
            // 投影后向量的normal，为正X正Y方向
            Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
            Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
            #endregion

            //#region 判断两个法向量的交点，是否都在两个向量的正方向上
            //bool isIntersectOnPositiveDirection = true;
            //// 即由交点到向量原点所构成的向量是否与原来的法向量方向相同
            //// 1.求交点
            //Point3d intersectPoint = GetLineIntersection(point1, nVectorForBS1, point2, nVectorForBS2);
            //if (intersectPoint == Point3d.Unset)
            //{
            //    throw new Exception("interval为奇数时划线，point1和point2的法向量平行");
            //}
            //// 2.构造新向量
            //Vector3d newVector1 = new Vector3d(point1 - intersectPoint);
            //Vector3d newVector2 = new Vector3d(point2 - intersectPoint);
            //// 3.判断新向量与原向量是否同向
            //if (newVector1 * nVectorForBS1 > 0 && newVector2 * nVectorForBS2 > 0)
            //{
            //    isIntersectOnPositiveDirection = true;
            //}
            //else
            //{
            //    isIntersectOnPositiveDirection = false;
            //}
            //#endregion

            //#region 判断由nVectorForBS1到nVectorForBS2是顺时针还是逆时针
            //bool isCounterClockwise = true;
            //Vector3d crossProduct = Vector3d.CrossProduct(nVectorForBS1, nVectorForBS2);
            //if (crossProduct.Z < 0)
            //{
            //    // 逆时针
            //    isCounterClockwise = true;
            //}
            //else
            //{
            //    // 顺时针
            //    isCounterClockwise = false;
            //}
            //#endregion
            #endregion

            #region 输出point1是否被绘制过
            bool isPoint1Drawn = false;
            int h = -1;
            List<int> hStartAtPoint1 = originPDeepCopy.Vertices.GetHalfedges(vertex1Index).ToList();
            for (int i = 0; i < hStartAtPoint1.Count; i++)
            {
                int face1 = originPDeepCopy.Halfedges[hStartAtPoint1[i]].AdjacentFace;
                int face2 = originPDeepCopy.Halfedges[originPDeepCopy.Halfedges.GetPairHalfedge(hStartAtPoint1[i])].AdjacentFace;
                if (face1 != -1 && face2 != -1)
                {
                    isPoint1Drawn = true;
                    h = hStartAtPoint1[i];
                    break;
                }
            }
            #endregion

            Point3d point1 = originPDeepCopy.Vertices[vertex1Index].ToPoint3d();
            Point3d point2 = originPDeepCopy.Vertices[vertex2Index].ToPoint3d();

            List<Point3d> turningPoints = new List<Point3d>();

            //int branchIndex = 0;
            int itemIndex = 0;

            PlanktonMesh P = new PlanktonMesh();
            List<PlanktonXYZ> newPlanktonXYZ;
            List<List<int>> newFaceVertexOrder;

            if (interval == 0)
            {
                // interval == 0
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = originPDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = originP.Vertices[endPoint1Index].ToPoint3d();

                    #region 计算turningPoint
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = Math.Abs(endPoint1.Y - point1.Y);

                        double t = tYDT.Branch(branchIndex)[itemIndex];
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

                        double t = tXDT.Branch(branchIndex)[itemIndex];
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

                    originPDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;
                    originPDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Remove(vertex1Index);
                    viFace.Add(turningPoint0Index);
                    viFace.Add(turningPoint1Index);

                    viNewFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint1Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = boundingBox.Max.Y - boundingBox.Min.Y;

                        double t = tYDT.Branch(branchIndex)[itemIndex];
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

                        double t = tXDT.Branch(branchIndex)[itemIndex];
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

                    originPDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;
                    originPDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint0Index);
                    viFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }
            else if (interval == 1)
            {
                // interval == 1
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = originPDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = originPDeepCopy.Vertices[endPoint1Index].ToPoint3d();

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

                    originPDeepCopy.Vertices.Add(intersectPoint);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viNewFace.Remove(vertex1Index);
                    if (isOutOfRange)
                    {
                        // 对于OutOfRange时的特殊处理
                        viFace.Add(turningPoint0Index);
                        viFace.Add(viNewFace.Last());
                    }
                    viNewFace.Add(turningPoint0Index);
                    viFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
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

                    originPDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }
            else
            {
                // interval == 2
                if (isPoint1Drawn)
                {
                    // 如果point1已经被绘制过
                    int endPoint1Index = originPDeepCopy.Halfedges.EndVertex(h);
                    Point3d endPoint1 = originP.Vertices[endPoint1Index].ToPoint3d();

                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = Math.Abs(endPoint1.Y - point1.Y);

                        double t2 = tYDT.Branch(branchIndex)[itemIndex];
                        double addY;
                        if (nVectorForBS1.Y > 0)
                        {
                            addY = y * t2;
                        }
                        else
                        {
                            addY = -y * t2;
                        }

                        turningPoints.Add(new Point3d(point1.X, point1.Y + addY, point1.Z));
                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                    }
                    else
                    {
                        // 法向量为x轴
                        double x = Math.Abs(endPoint1.X - point1.X);

                        double t2 = tXDT.Branch(branchIndex)[itemIndex];
                        double addX;
                        if (nVectorForBS1.X > 0)
                        {
                            addX = x * t2;
                        }
                        else
                        {
                            addX = -x * t2;
                        }

                        turningPoints.Add(new Point3d(point1.X + addX, point1.Y, point1.Z));
                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                    }

                    //branchIndex++;

                    //#region 判断交点是否在经过point1的边的内部
                    //Vector3d vector1 = new Vector3d(endPoint1 - point1);
                    //Vector3d vector2 = new Vector3d(turningPoints[0] - point1);
                    //double t = vector2.Length / vector1.Length;
                    //bool isOutOfRange = false;
                    //if (t > 1)
                    //{
                    //    isOutOfRange = true;
                    //}
                    //else
                    //{
                    //    isOutOfRange = false;
                    //}
                    //#endregion

                    originPDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;
                    originPDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = turningPoint0Index;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    //if (isOutOfRange)
                    //{
                    //    viNewFace.Remove(vertex1Index);

                    //    viFace.Add(viNewFace[0]);
                    //    viFace.Add(turningPoint0Index);
                    //    viFace.Add(turningPoint1Index);

                    //    viNewFace.Add(turningPoint1Index);
                    //    viNewFace.Add(turningPoint0Index);
                    //}
                    //else
                    //{
                    //    viFace.Add(turningPoint0Index);
                    //    viFace.Add(turningPoint1Index);

                    //    viNewFace.Remove(vertex1Index);
                    //    viNewFace.Insert(0, turningPoint0Index);
                    //    viNewFace.Add(turningPoint1Index);
                    //}
                    viNewFace.Remove(vertex1Index);
                    viNewFace.RemoveAt(0);

                    //viFace.Add(viNewFace[0]);
                    viFace.Add(turningPoint1Index);
                    viFace.Add(turningPoint0Index);

                    viNewFace.Add(turningPoint0Index);
                    viNewFace.Add(turningPoint1Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    if (nVectorForBS1.X == 0)
                    {
                        // 法向量为y轴
                        double y = boundingBox.Max.Y - boundingBox.Min.Y;

                        double t = tYDT.Branch(branchIndex)[itemIndex];
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

                        double t = tXDT.Branch(branchIndex)[itemIndex];
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

                    originPDeepCopy.Vertices.Add(turningPoints[0]);
                    int turningPoint0Index = originPDeepCopy.Vertices.Count - 1;
                    originPDeepCopy.Vertices.Add(turningPoints[1]);
                    int turningPoint1Index = originPDeepCopy.Vertices.Count - 1;

                    // 输出commonVertexIndex
                    commonVertexIndex = -1;

                    #region 找到viAroundOriginFace的两个切片
                    int viFaceIndex = GetViFaceIndex(interval, vertex1Index, vertex2Index, originPDeepCopy);

                    List<int> viFace, viNewFace;
                    GetInitialViLists(vertex1Index, vertex2Index, originPDeepCopy, viFaceIndex, out viFace, out viNewFace);
                    #endregion

                    #region 处理生成最后的viFace和viNewFace
                    viFace.Add(turningPoint0Index);
                    viFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint1Index);
                    viNewFace.Add(turningPoint0Index);
                    #endregion

                    #region 生成新的PlanktonMesh
                    GenerateNewPlanktonXYZAndNewFaceVertexOrder(originPDeepCopy, out newPlanktonXYZ, out newFaceVertexOrder, viFaceIndex, viFace, viNewFace);
                    #endregion
                }
            }

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

        private void GetInitialViLists(int vertex1Index, int vertex2Index, PlanktonMesh originPDeepCopy, int adjacentFaceIndex, out List<int> viFace, out List<int> viNewFace)
        {
            List<int> viAroundOriginFace = new List<int>();
            viAroundOriginFace = originPDeepCopy.Faces.GetFaceVertices(adjacentFaceIndex).ToList<int>();

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

        private void GenerateNewPlanktonXYZAndNewFaceVertexOrder(PlanktonMesh originPDeepCopy, out List<PlanktonXYZ> newPlanktonXYZ, out List<List<int>> newFaceVertexOrder, int viFaceIndex, List<int> viFace, List<int> viNewFace)
        {
            // 转移原来originPDeepCopy的属性
            newPlanktonXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < originPDeepCopy.Vertices.Count; i++)
            {
                newPlanktonXYZ.Add(originPDeepCopy.Vertices[i].ToXYZ());
            }
            //
            newFaceVertexOrder = new List<List<int>>();
            for (int i = 0; i < originPDeepCopy.Faces.Count; i++)
            {
                if (i == viFaceIndex)
                {
                    newFaceVertexOrder.Add(new List<int>());
                    newFaceVertexOrder[i].AddRange(viFace);
                }
                else
                {
                    newFaceVertexOrder.Add(new List<int>());
                    int[] faceVertexOrder = originPDeepCopy.Faces.GetFaceVertices(i);
                    newFaceVertexOrder[i].AddRange(faceVertexOrder);
                }
            }
            newFaceVertexOrder.Add(viNewFace);
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
    }
}