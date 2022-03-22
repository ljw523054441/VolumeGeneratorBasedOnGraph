using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Attributes;
using Plankton;
using PlanktonGh;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;
using Grasshopper.GUI.Canvas;

namespace VolumeGeneratorBasedOnGraph.VolumeAlgorithm
{
    public class GhcVolumeBaseOnTopo : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcVolumeBaseOnTopo class.
        /// </summary>
        public GhcVolumeBaseOnTopo()
          : base("VolumeBaseOnTopo", "VolumeBaseOnTopo",
              "Description",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            InnerResultPolyline = new List<Polyline>();
            InnerNodeTextDots = new List<TextDot>();

            // Temp
            SCrvForShowList = new List<Curve>();
            NCrvForShowList = new List<Curve>();
            ECrvForShowList = new List<Curve>();
            WCrvForShowList = new List<Curve>();
            GapCenterLineList = new List<Curve>();

            SBCrvForShowList = new List<Curve>();
            NBCrvForShowList = new List<Curve>();
            EBCrvForShowList = new List<Curve>();
            WBCrvForShowList = new List<Curve>();

            AllBlockCellBoundaryDT = new DataTree<Curve>();
            //AllBlockCellBoundaryDT2 = new DataTree<Curve>();
        }

        private Random m_random;

        private int Thickness;
        private List<Polyline> InnerResultPolyline;

        private List<TextDot> InnerNodeTextDots;

        // Temp
        private List<Curve> SCrvForShowList;
        private List<Curve> NCrvForShowList;
        private List<Curve> ECrvForShowList;
        private List<Curve> WCrvForShowList;
        private List<Curve> GapCenterLineList;

        private List<Curve> SBCrvForShowList;
        private List<Curve> NBCrvForShowList;
        private List<Curve> EBCrvForShowList;
        private List<Curve> WBCrvForShowList;

        private DataTree<Curve> AllBlockCellBoundaryDT;
        //private DataTree<Curve> AllBlockCellBoundaryDT2;

        private HorizontalVolumeBranchType HBranchType;
        private VerticalGapBranchType VBranchType;
        private enum HorizontalVolumeBranchType
        {
            OnlyOneFull = 0,
            OnlyOneSouth = 1,
            OnlyOneCenter = 2,
            OnlyOneNorth = 3,
            TwoNorthShorten = 4,
            Two = 5,
            AboveThree = 6
        }

        private enum VerticalGapBranchType
        {
            ZeroGapNeedShorten = 0,
            OneGap = 1,
            ZeroGap = 2,
            UncertainGap = 3,
            ZeroGapNorthCull = 4,
            AboveOneGap = 5
        }


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DualGraphWithHM", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("newDualPlanktonMesh", "newDualPlanktonMesh", "", GH_ParamAccess.item);


            pManager.AddIntegerParameter("WhichPairSetbackNeedToChange", "", "", GH_ParamAccess.tree);
            pManager.AddNumberParameter("SetbackValueNeedToChange", "", "", GH_ParamAccess.tree);
            pManager[2].Optional = true;
            pManager[3].Optional = true;

            pManager.AddNumberParameter("SetbackList", "SetbackList", "", GH_ParamAccess.list);

            pManager.AddNumberParameter("w", "", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("W", "", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("lMin", "", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("lMax", "", "", GH_ParamAccess.item);
            pManager.AddNumberParameter("d", "", "", GH_ParamAccess.item);

            pManager.AddNumberParameter("FloorHeight", "", "", GH_ParamAccess.item);

            pManager.AddNumberParameter("Density", "", "", GH_ParamAccess.list);
            pManager.AddNumberParameter("AreaRatio", "", "", GH_ParamAccess.list);

            pManager.AddIntegerParameter("Seed", "", "", GH_ParamAccess.item);
            pManager.AddBooleanParameter("IsJitter", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM with BS","", "", GH_ParamAccess.item);
            
            pManager.AddGeometryParameter("SetBack Polyline", "", "", GH_ParamAccess.list);

            //pManager.AddGenericParameter("BoundarySegments", "", "", GH_ParamAccess.tree);

            // pManager.AddCurveParameter("CutRegion Debug", "", "", GH_ParamAccess.tree);

            // pManager.AddCurveParameter("Public Side", "", "", GH_ParamAccess.tree);

            // pManager.AddBrepParameter("CutsBrepDT", "", "", GH_ParamAccess.tree);
            // pManager.AddBrepParameter("ResultBreps", "", "", GH_ParamAccess.list);
            pManager.AddBrepParameter("Building", "", "", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Contours", "", "", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Holes", "", "", GH_ParamAccess.tree);

            pManager.AddCurveParameter("HCurves", "", "", GH_ParamAccess.tree);
            pManager.AddCurveParameter("VCurves", "", "", GH_ParamAccess.tree);

            //pManager.AddCurveParameter("HCB", "", "", GH_ParamAccess.tree);
            //pManager.AddCurveParameter("VCB", "", "", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            PlanktonMesh dual = new PlanktonMesh();

            PlanktonMesh graph = new PlanktonMesh();
            List<GraphNode> graphNodes = new List<GraphNode>();
            List<List<int>> graphTable = new List<List<int>>();
            List<int> dFIndex_PVIndex = new List<int>();

            GH_Structure<GH_Integer> pairNeedToChange = new GH_Structure<GH_Integer>();
            GH_Structure<GH_Number> setbackValueNeedToChange = new GH_Structure<GH_Number>();

            List<List<int>> pairNeedToChangeLoL = new List<List<int>>();
            List<List<double>> setbackValueNeedToChangeLoL = new List<List<double>>();

            List<double> setbackList = new List<double>();

            double w = 0;
            double W = 0;
            double lMin = 0;
            double lMax = 0;
            double d = 0;

            double floorHeight = 0;

            List<double> densitys = new List<double>();
            List<double> areaRatio = new List<double>();

            int seed = 0;

            bool isJitter = false;

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetData<PlanktonMesh>("newDualPlanktonMesh", ref dual)
                && DA.GetDataList("SetbackList", setbackList)
                && DA.GetData("w", ref w)
                && DA.GetData("W", ref W)
                && DA.GetData("lMin", ref lMin)
                && DA.GetData("lMax", ref lMax)
                && DA.GetData("d", ref d)
                && DA.GetData("FloorHeight", ref floorHeight)
                && DA.GetDataList("Density", densitys)
                && DA.GetDataList("AreaRatio", areaRatio)
                && DA.GetData<int>("Seed", ref seed)
                && DA.GetData<bool>("IsJitter", ref isJitter))
            {
                DA.GetDataTree<GH_Integer>("WhichPairSetbackNeedToChange", out pairNeedToChange);
                DA.GetDataTree<GH_Number>("SetbackValueNeedToChange", out setbackValueNeedToChange);

                // 固定种子随机数生成器
                m_random = new Random(seed);

                Thickness = 2;

                Cell.Tolerance = GH_Component.DocumentTolerance();

                dualGraphWithHM.DualPlanktonMesh = new PlanktonMesh(dual);

                graph = dualGraphWithHM.PlanktonMesh;
                graphNodes = dualGraphWithHM.GraphNodes;
                graphTable = dualGraphWithHM.GraphTables;
                dFIndex_PVIndex = dualGraphWithHM.DFIndex_PVIndex;

                pairNeedToChangeLoL = UtilityFunctions.GH_StructureToLoL_Int(pairNeedToChange);
                setbackValueNeedToChangeLoL = UtilityFunctions.GH_StructureToLoL_Double(setbackValueNeedToChange);
                #region 检查输入异常
                if (pairNeedToChangeLoL.Count != setbackValueNeedToChangeLoL.Count)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查序号对与值对的数量");
                    return;
                }
                for (int i = 0; i < pairNeedToChangeLoL.Count; i++)
                {
                    if (pairNeedToChangeLoL[i].Count != 2
                        || setbackValueNeedToChangeLoL[i].Count != 2)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查序号对与值对的数量");
                        return;
                    }
                }

                for (int i = 0; i < pairNeedToChangeLoL.Count; i++)
                {
                    if (!graphTable[pairNeedToChangeLoL[i][0]].Contains(pairNeedToChangeLoL[i][1]))
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "输入的序号对中包含未定义的拓扑连接");
                        return;
                    }
                }
                #endregion


                // 各边退界
                double setbackW = setbackList[0];
                double setbackS = setbackList[1];
                double setbackE = setbackList[2];
                double setbackN = setbackList[3];
                double setback = setbackList[4];

                List<List<double>> setbackInner = new List<List<double>>();

                #region 得到所有的innerNode的序号和outerNode的序号
                List<int> innerNodeIndexs = dualGraphWithHM.InnerNodeIndexList;
                List<int> outerNodeIndexs = dualGraphWithHM.OuterNodeIndexList;
                int innerNodeCount = dualGraphWithHM.InnerNodeCount;
                int outerNodeCount = dualGraphWithHM.OuterNodeCount;
                #endregion

                List<Polyline> facePolylines = new List<Polyline>();
                for (int i = 0; i < dual.Faces.Count; i++)
                {
                    int[] faceVerticesIndex = dual.Faces.GetFaceVertices(i);
                    List<Point3d> faceVertices = new List<Point3d>();
                    for (int j = 0; j < faceVerticesIndex.Length; j++)
                    {
                        faceVertices.Add(dual.Vertices[faceVerticesIndex[j]].ToPoint3d());
                    }
                    Polyline pl = new Polyline(faceVertices);
                    facePolylines.Add(pl);
                }

                #region 将HIndex分类
                List<List<List<int>>> allFaceHIndex = new List<List<List<int>>>();
                List<List<int>> allFaceFIndex = new List<List<int>>();
                List<List<int>> allFaceAdjacentFIndex = new List<List<int>>();
                for (int i = 0; i < dual.Faces.Count; i++)
                {
                    allFaceHIndex.Add(new List<List<int>>());
                    allFaceFIndex.Add(new List<int>());
                    allFaceAdjacentFIndex.Add(new List<int>());

                    int[] currentFaceHalfedges = dual.Faces.GetHalfedges(i);
                    int lastAdjacentFIndex = -2;
                    for (int j = 0; j < currentFaceHalfedges.Length; j++)
                    {
                        int pairH = dual.Halfedges.GetPairHalfedge(currentFaceHalfedges[j]);
                        int adjacentF = dual.Halfedges[pairH].AdjacentFace;
                        if (adjacentF != lastAdjacentFIndex)
                        {
                            allFaceFIndex[i].Add(i);

                            allFaceHIndex[i].Add(new List<int>());
                            allFaceHIndex[i].Last().Add(currentFaceHalfedges[j]);

                            allFaceAdjacentFIndex[i].Add(adjacentF);
                        }
                        else
                        {
                            if (adjacentF == -1)
                            {
                                allFaceFIndex[i].Add(i);

                                allFaceHIndex[i].Add(new List<int>());
                                allFaceHIndex[i].Last().Add(currentFaceHalfedges[j]);

                                allFaceAdjacentFIndex[i].Add(adjacentF);
                            }
                            else
                            {
                                allFaceHIndex[i].Last().Add(currentFaceHalfedges[j]);
                            }
                        }

                        lastAdjacentFIndex = adjacentF;
                    }
                }
                #endregion

                #region 按照分类情况，构造BS
                List<List<BoundarySegment>> allFaceBS = new List<List<BoundarySegment>>();
                for (int i = 0; i < dual.Faces.Count; i++)
                {
                    allFaceBS.Add(new List<BoundarySegment>());
                    for (int j = 0; j < allFaceHIndex[i].Count; j++)
                    {
                        List<int> currVIndex = new List<int>();
                        currVIndex.Add(dual.Halfedges[allFaceHIndex[i][j][0]].StartVertex);
                        for (int k = 0; k < allFaceHIndex[i][j].Count; k++)
                        {
                            currVIndex.Add(dual.Halfedges.EndVertex(allFaceHIndex[i][j][k]));
                        }

                        List<Point3d> currPoints = new List<Point3d>();
                        for (int k = 0; k < currVIndex.Count; k++)
                        {
                            currPoints.Add(dual.Vertices[currVIndex[k]].ToPoint3d());
                        }

                        BoundarySegment bs = new BoundarySegment(currPoints, allFaceHIndex[i][j].ToArray(), allFaceFIndex[i][j], allFaceAdjacentFIndex[i][j]);

                        // 判断首尾两个点位置处的曲率，即是否在首尾处转折
                        int FH = allFaceHIndex[i][j].First();
                        int LH = allFaceHIndex[i][j].Last();
                        int PFH = dual.Halfedges[FH].PrevHalfedge;
                        int NLH = dual.Halfedges[LH].NextHalfedge;

                        int PFHS = dual.Halfedges[PFH].StartVertex;
                        int NLHE = dual.Halfedges.EndVertex(NLH);

                        Point3d pt0 = dual.Vertices[PFHS].ToPoint3d();
                        Point3d pt1 = currPoints.First();
                        Point3d pt2 = currPoints[1];

                        Point3d pt3 = currPoints[currPoints.Count - 2];
                        Point3d pt4 = currPoints.Last();
                        Point3d pt5 = dual.Vertices[NLHE].ToPoint3d();
                        if (IsThreePointsCollinear(pt0, pt1, pt2))
                        {
                            bs.IsZeroCurvatureAtStart = true;
                        }
                        else
                        {
                            bs.IsZeroCurvatureAtStart = false;
                        }
                        if (IsThreePointsCollinear(pt3, pt4, pt5))
                        {
                            bs.IsZeroCurvatureAtEnd = true;
                        }
                        else
                        {
                            bs.IsZeroCurvatureAtEnd = false;
                        }
                        allFaceBS[i].Add(bs);
                    }
                }
                #endregion

                #region 明确每个Inner点的连接关系，即DFace的所在方位
                // 注意graphTable只包含使用者输入的连接关系，不包含使用者未定义的拓扑关系
                List<List<int>> inner_Outer = new List<List<int>>();
                for (int i = 0; i < graph.Vertices.Count; i++)
                {
                    if (i < innerNodeCount)
                    {
                        inner_Outer.Add(new List<int>());
                        int[] currentConnection = graph.Vertices.GetVertexNeighbours(i);
                        for (int j = 0; j < currentConnection.Length; j++)
                        {
                            if (currentConnection[j] >= innerNodeCount)
                            {
                                inner_Outer[i].Add(currentConnection[j]);
                            }
                        }
                    }
                }

                for (int i = 0; i < inner_Outer.Count; i++)
                {
                    // 按照outer的序号进行升序排序
                    inner_Outer[i].Sort();
                    // 在同时包含首尾时，调整首尾两个的位置，即346变为634，以保证后面与BS的匹配
                    if (inner_Outer[i].Contains(innerNodeCount) && inner_Outer[i].Contains(innerNodeCount + outerNodeCount - 1))
                    {
                        int item = inner_Outer[i].Last();
                        inner_Outer[i].Insert(0, item);
                        inner_Outer[i].RemoveAt(inner_Outer[i].Count - 1);
                    }
                }

                #endregion

                #region 根据DFace的所在方位，确定每个BS要退线的距离，并为BS设置在哪个方位的标签
                for (int i = 0; i < innerNodeCount; i++)
                {
                    int dFIndex = dFIndex_PVIndex.IndexOf(i);

                    List<BoundarySegment> outerBS = new List<BoundarySegment>();
                    int count = 0;
                    for (int j = 0; j < allFaceBS[dFIndex].Count; j++)
                    {
                        if (allFaceBS[dFIndex][j].AdjacentFIndex != -1)
                        {
                            // 不是outerBS时
                            int[] array = new int[2] { i, dFIndex_PVIndex[allFaceBS[dFIndex][j].AdjacentFIndex] };


                            bool flag = false;
                            for (int k = 0; k < pairNeedToChangeLoL.Count; k++)
                            {
                                if (pairNeedToChangeLoL[k].Except(array).ToList().Count == 0)
                                {
                                    int index0 = pairNeedToChangeLoL[k].IndexOf(array[0]);
                                    //int index1 = pairNeedToChangeLoL[k].IndexOf(array[1]);

                                    double finalSetback = setback + setbackValueNeedToChangeLoL[k][index0];
                                    if (finalSetback < 0)
                                    {
                                        finalSetback = 0;
                                    }
                                    //allFaceBSSetBack[dFIndex].Add(finalSetback);
                                    allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.Inner;
                                    allFaceBS[dFIndex][j].setback = finalSetback;
                                    flag = true;
                                    break;
                                }
                            }
                            if (!flag)
                            {
                                //allFaceBSSetBack[dFIndex].Add(setback);
                                allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.Inner;
                                allFaceBS[dFIndex][j].setback = setback;
                            }
                        }
                        else
                        {
                            // 是outerBS时
                            // 需要测试一下这里是否需要排序
                            if (inner_Outer[i][count] - innerNodeCount == 0)
                            {
                                //allFaceBSSetBack[dFIndex].Add(setbackW);
                                allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.W;
                                allFaceBS[dFIndex][j].setback = setbackW;
                                outerBS.Add(allFaceBS[dFIndex][j]);
                            }
                            else if (inner_Outer[i][count] - innerNodeCount == 1)
                            {
                                //allFaceBSSetBack[dFIndex].Add(setbackS);
                                allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.S;
                                allFaceBS[dFIndex][j].setback = setbackS;
                                outerBS.Add(allFaceBS[dFIndex][j]);
                            }
                            else if (inner_Outer[i][count] - innerNodeCount == 2)
                            {
                                //allFaceBSSetBack[dFIndex].Add(setbackE);
                                allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.E;
                                allFaceBS[dFIndex][j].setback = setbackE;
                                outerBS.Add(allFaceBS[dFIndex][j]);
                            }
                            else
                            {
                                //allFaceBSSetBack[dFIndex].Add(setbackN);
                                allFaceBS[dFIndex][j].LocationValue = BoundarySegment.Location.N;
                                allFaceBS[dFIndex][j].setback = setbackN;
                                outerBS.Add(allFaceBS[dFIndex][j]);
                            }

                            count++;
                        }
                    }
                }
                #endregion

                #region 退线
                List<Polyline> innerResultPolylines;
                List<List<BoundarySegment>> newAllFaceBS = OffsetBS(allFaceBS,
                                                                    facePolylines,
                                                                    out innerResultPolylines);
                #endregion
                DA.SetDataList("SetBack Polyline", innerResultPolylines);

                dualGraphWithHM.OffsetedBoundarySegments = new List<List<BoundarySegment>>();
                for (int i = 0; i < newAllFaceBS.Count; i++)
                {
                    dualGraphWithHM.OffsetedBoundarySegments.Add(new List<BoundarySegment>());
                    for (int j = 0; j < newAllFaceBS[i].Count; j++)
                    {
                        dualGraphWithHM.OffsetedBoundarySegments[i].Add(new BoundarySegment(newAllFaceBS[i][j]));
                    }
                }
                DA.SetData("DualGraphWithHM with BS", dualGraphWithHM);


                // 找到整个场地的中心点
                List<Point3d> outerPoints = new List<Point3d>();
                for (int i = 0; i < graph.Halfedges.Count; i++)
                {
                    if (graph.Halfedges[i].AdjacentFace == -1)
                    {
                        outerPoints.Add(graph.Vertices[graph.Halfedges[i].StartVertex].ToPoint3d());
                    }
                }
                Point3d center = new Point3d(0, 0, 0);
                for (int i = 0; i < outerPoints.Count; i++)
                {
                    center += outerPoints[i];
                }
                center /= outerPoints.Count;

                SCrvForShowList.Clear();
                NCrvForShowList.Clear();
                ECrvForShowList.Clear();
                WCrvForShowList.Clear();
                GapCenterLineList.Clear();

                SBCrvForShowList.Clear();
                NBCrvForShowList.Clear();
                EBCrvForShowList.Clear();
                WBCrvForShowList.Clear();

                AllBlockCellBoundaryDT.ClearData();

                DataTree<double> allBlockUnOffsetedArea = new DataTree<double>();

                DataTree<Curve> allBlockHCurvesDT = new DataTree<Curve>();
                DataTree<Curve> allBlockVCurvesDT = new DataTree<Curve>();
                DataTree<Polyline> allBlockContourDT = new DataTree<Polyline>();
                DataTree<Polyline> allBlockHoleDT = new DataTree<Polyline>();

                DataTree<Curve> allBlockCellBoundaryDT = new DataTree<Curve>();

                DataTree<Cell> allBlockBasicCellDT = new DataTree<Cell>();

                DataTree<Polyline> allBlockPolylineDT = new DataTree<Polyline>();

                List<bool> isBlockSplitedList = new List<bool>();
                for (int i = 0; i < newAllFaceBS.Count; i++)
                {
                    #region 向path中添加第一层大Block层级的序号
                    GH_Path ghPathForUnOffsetedArea = new GH_Path(i);

                    //GH_Path ghPathForBrep = new GH_Path(i);
                    GH_Path ghPathForHCurve = new GH_Path(i);
                    GH_Path ghPathForVCurve = new GH_Path(i);
                    GH_Path ghPathForContour = new GH_Path(i);
                    GH_Path ghPathForHole = new GH_Path(i);

                    GH_Path ghPathForCellBoundary = new GH_Path(i);

                    GH_Path ghPathForCell = new GH_Path(i);

                    GH_Path ghPathForBlockPolyline = new GH_Path(i);
                    #endregion

                    bool isHorizontalLayout;
                    Polyline sortedBoundaryPolyline;
                    List<Line> lineForEachEdge;
                    List<BoundarySegment> sortedBS = SortBoundarySegment(newAllFaceBS[i], innerResultPolylines[i], out isHorizontalLayout, out lineForEachEdge, out sortedBoundaryPolyline);

                    List<bool> isGenerateables;
                    // todo:根据BS的属性（两个Face相连的BS的setback值），来生成 List<bool> isGenerateables

                    if (sortedBoundaryPolyline.Count == 5)
                    {
                        isBlockSplitedList.Add(false);

                        #region 向path中添加第二层小Block层级的序号
                        ghPathForUnOffsetedArea = ghPathForUnOffsetedArea.AppendElement(0);

                        ghPathForHCurve = ghPathForHCurve.AppendElement(0);
                        ghPathForVCurve = ghPathForVCurve.AppendElement(0);
                        ghPathForContour = ghPathForContour.AppendElement(0);
                        ghPathForHole = ghPathForHole.AppendElement(0);

                        ghPathForCellBoundary = ghPathForCellBoundary.AppendElement(0);

                        ghPathForCell = ghPathForCell.AppendElement(0);

                        ghPathForBlockPolyline = ghPathForBlockPolyline.AppendElement(0);
                        #endregion

                        #region isGenerateables控制
                        isGenerateables = new List<bool>() { true, true, true, true };
                        #endregion

                        // 设置此时的turningPoint
                        int maxIndex = 0;
                        double maxDis = sortedBoundaryPolyline[0].DistanceTo(center);
                        for (int j = 0; j < sortedBoundaryPolyline.ToArray().Length; j++)
                        {
                            double newMaxDis = sortedBoundaryPolyline[j].DistanceTo(center);
                            if (newMaxDis > maxDis)
                            {
                                maxIndex = j;
                                maxDis = newMaxDis;
                            }
                        }
                        

                        List<List<List<Curve>>> hCurveLoL;
                        List<List<List<Curve>>> vCurveLoL;
                        //List<Polyline> contourList;
                        List<List<List<Polyline>>> holeLoL;
                        List<List<List<Polyline>>> contourLoL;

                        // Temp
                        Curve southBaseLineForShow;
                        Curve northBaseLineForShow;
                        Curve eastBaseLineForShow;
                        Curve westBaseLineForShow;

                        Curve southBoundaryLineForShow;
                        Curve northBoundaryLineForShow;
                        Curve eastBoundaryLineForShow;
                        Curve westBoundaryLineForShow;

                        List<Curve> gap;

                        List<List<Curve>> cellBoundaryLoL;

                        List<List<Cell>> cellLoL = GenerateVolumeBaseLine(w,
                                                                                W,
                                                                                lMin,
                                                                                lMax,
                                                                                d,
                                                                                isJitter,

                                                                                false,
                                                                                sortedBoundaryPolyline[maxIndex],

                                                                                out hCurveLoL,
                                                                                out vCurveLoL,
                                                                                out contourLoL,
                                                                                out holeLoL,

                                                                                out southBaseLineForShow,
                                                                                out northBaseLineForShow,
                                                                                out eastBaseLineForShow,
                                                                                out westBaseLineForShow,

                                                                                out southBoundaryLineForShow,
                                                                                out northBoundaryLineForShow,
                                                                                out eastBoundaryLineForShow,
                                                                                out westBoundaryLineForShow,

                                                                                out gap,

                                                                                out cellBoundaryLoL,

                                                                                sortedBoundaryPolyline,
                                                                                lineForEachEdge,
                                                                                isGenerateables);

                        // Temp
                        SCrvForShowList.Add(southBaseLineForShow);
                        NCrvForShowList.Add(northBaseLineForShow);
                        ECrvForShowList.Add(eastBaseLineForShow);
                        WCrvForShowList.Add(westBaseLineForShow);

                        SBCrvForShowList.Add(southBoundaryLineForShow);
                        NBCrvForShowList.Add(northBoundaryLineForShow);
                        EBCrvForShowList.Add(eastBoundaryLineForShow);
                        WBCrvForShowList.Add(westBoundaryLineForShow);

                        GapCenterLineList.AddRange(gap);

                        #region output dataTree
                        allBlockUnOffsetedArea.EnsurePath(ghPathForUnOffsetedArea);
                        List<Point3d> pts = facePolylines[i].ToList();
                        pts.Add(pts[0]);
                        allBlockUnOffsetedArea.Branch(ghPathForUnOffsetedArea).Add(AreaMassProperties.Compute(new Polyline(pts).ToNurbsCurve()).Area);

                        allBlockPolylineDT.EnsurePath(ghPathForBlockPolyline);
                        allBlockPolylineDT.Branch(ghPathForBlockPolyline).Add(sortedBoundaryPolyline);

                        for (int j = 0; j < cellLoL.Count; j++)
                        {
                            ghPathForCell = ghPathForCell.AppendElement(j);
                            allBlockBasicCellDT.EnsurePath(ghPathForCell);
                            allBlockBasicCellDT.Branch(ghPathForCell).AddRange(cellLoL[j]);
                            ghPathForCell = ghPathForCell.CullElement();
                        }

                        for (int j = 0; j < hCurveLoL.Count; j++)
                        {
                            ghPathForHCurve = ghPathForHCurve.AppendElement(j);
                            for (int k = 0; k < hCurveLoL[j].Count; k++)
                            {
                                ghPathForHCurve = ghPathForHCurve.AppendElement(k);

                                allBlockHCurvesDT.EnsurePath(ghPathForHCurve);
                                allBlockHCurvesDT.Branch(ghPathForHCurve).AddRange(hCurveLoL[j][k]);

                                ghPathForHCurve = ghPathForHCurve.CullElement();
                            }
                            ghPathForHCurve = ghPathForHCurve.CullElement();
                        }

                        for (int j = 0; j < vCurveLoL.Count; j++)
                        {
                            ghPathForVCurve = ghPathForVCurve.AppendElement(j);
                            for (int k = 0; k < vCurveLoL[j].Count; k++)
                            {
                                ghPathForVCurve = ghPathForVCurve.AppendElement(k);

                                allBlockVCurvesDT.EnsurePath(ghPathForVCurve);
                                allBlockVCurvesDT.Branch(ghPathForVCurve).AddRange(vCurveLoL[j][k]);

                                ghPathForVCurve = ghPathForVCurve.CullElement();
                            }
                            ghPathForVCurve = ghPathForVCurve.CullElement();
                        }

                        for (int j = 0; j < contourLoL.Count; j++)
                        {
                            ghPathForContour = ghPathForContour.AppendElement(j);
                            for (int k = 0; k < contourLoL[j].Count; k++)
                            {
                                ghPathForContour = ghPathForContour.AppendElement(k);

                                allBlockContourDT.EnsurePath(ghPathForContour);
                                allBlockContourDT.Branch(ghPathForContour).AddRange(contourLoL[j][k]);

                                ghPathForContour = ghPathForContour.CullElement();
                            }
                            ghPathForContour = ghPathForContour.CullElement();
                        }

                        for (int j = 0; j < holeLoL.Count; j++)
                        {
                            ghPathForHole = ghPathForHole.AppendElement(j);
                            for (int k = 0; k < holeLoL[j].Count; k++)
                            {
                                ghPathForHole = ghPathForHole.AppendElement(k);

                                allBlockHoleDT.EnsurePath(ghPathForHole);
                                allBlockHoleDT.Branch(ghPathForHole).AddRange(holeLoL[j][k]);

                                ghPathForHole = ghPathForHole.CullElement();
                            }
                            ghPathForHole = ghPathForHole.CullElement();
                        }

                        for (int j = 0; j < cellBoundaryLoL.Count; j++)
                        {
                            ghPathForCellBoundary = ghPathForCellBoundary.AppendElement(j);
                            allBlockCellBoundaryDT.EnsurePath(ghPathForCellBoundary);
                            allBlockCellBoundaryDT.Branch(ghPathForCellBoundary).AddRange(cellBoundaryLoL[j]);
                            ghPathForCellBoundary = ghPathForCellBoundary.CullElement();
                        }
                        #endregion

                        #region 清除当前第二层级的序号
                        ghPathForUnOffsetedArea = ghPathForUnOffsetedArea.CullElement();

                        ghPathForHCurve = ghPathForHCurve.CullElement();
                        ghPathForVCurve = ghPathForVCurve.CullElement();
                        ghPathForContour = ghPathForContour.CullElement();
                        ghPathForHole = ghPathForHole.CullElement();

                        ghPathForCellBoundary = ghPathForCellBoundary.CullElement();

                        ghPathForCell = ghPathForCell.CullElement();

                        ghPathForBlockPolyline = ghPathForBlockPolyline.CullElement();
                        #endregion
                    }
                    else
                    {
                        isBlockSplitedList.Add(true);

                        List<double> splitedAreas;

                        List<List<Line>> newLineforEachEdgeLoL;
                        List<List<bool>> isGenerateablesLoL;
                        List<bool> isSmallerList;
                        Point3d turningPt = Point3d.Unset;
                        List<Polyline> splitedPolylines = SplitBlockIntoQuadBlock(sortedBoundaryPolyline, facePolylines[i], out splitedAreas, out newLineforEachEdgeLoL, out isGenerateablesLoL, out isSmallerList, out turningPt);



                        for (int j = 0; j < splitedPolylines.Count; j++)
                        {
                            #region 向path中添加第二层小Block层级的序号
                            ghPathForUnOffsetedArea = ghPathForUnOffsetedArea.AppendElement(j);

                            ghPathForHCurve = ghPathForHCurve.AppendElement(j);
                            ghPathForVCurve = ghPathForVCurve.AppendElement(j);
                            ghPathForContour = ghPathForContour.AppendElement(j);
                            ghPathForHole = ghPathForHole.AppendElement(j);

                            ghPathForCellBoundary = ghPathForCellBoundary.AppendElement(j);

                            ghPathForCell = ghPathForCell.AppendElement(j);

                            ghPathForBlockPolyline = ghPathForBlockPolyline.AppendElement(j);
                            #endregion

                            if (splitedPolylines[j].Count == 5)
                            {
                                List<List<List<Curve>>> hCurveLoL;
                                List<List<List<Curve>>> vCurveLoL;
                                List<List<List<Polyline>>> holeLoL;
                                List<List<List<Polyline>>> contourLoL;

                                // Temp
                                Curve southBaseLineForShow;
                                Curve northBaseLineForShow;
                                Curve eastBaseLineForShow;
                                Curve westBaseLineForShow;

                                Curve southBoundaryLineForShow;
                                Curve northBoundaryLineForShow;
                                Curve eastBoundaryLineForShow;
                                Curve westBoundaryLineForShow;

                                List<Curve> gap;

                                List<List<Curve>> cellBoundaryLoL;

                                List<List<Cell>> cellLoL = GenerateVolumeBaseLine(w,
                                                                             W,
                                                                             lMin,
                                                                             lMax,
                                                                             d,
                                                                             isJitter,

                                                                             isSmallerList[j],
                                                                             turningPt,

                                                                             out hCurveLoL,
                                                                             out vCurveLoL,
                                                                             out contourLoL,
                                                                             out holeLoL,

                                                                             out southBaseLineForShow,
                                                                             out northBaseLineForShow,
                                                                             out eastBaseLineForShow,
                                                                             out westBaseLineForShow,

                                                                             out southBoundaryLineForShow,
                                                                             out northBoundaryLineForShow,
                                                                             out eastBoundaryLineForShow,
                                                                             out westBoundaryLineForShow,

                                                                             out gap,

                                                                             out cellBoundaryLoL,

                                                                             splitedPolylines[j],
                                                                             newLineforEachEdgeLoL[j],
                                                                             isGenerateablesLoL[j]);

                                // Temp
                                SCrvForShowList.Add(southBaseLineForShow);
                                NCrvForShowList.Add(northBaseLineForShow);
                                ECrvForShowList.Add(eastBaseLineForShow);
                                WCrvForShowList.Add(westBaseLineForShow);

                                SBCrvForShowList.Add(southBoundaryLineForShow);
                                NBCrvForShowList.Add(northBoundaryLineForShow);
                                EBCrvForShowList.Add(eastBoundaryLineForShow);
                                WBCrvForShowList.Add(westBoundaryLineForShow);

                                GapCenterLineList.AddRange(gap);

                                #region output dataTree
                                allBlockUnOffsetedArea.EnsurePath(ghPathForUnOffsetedArea);
                                List<Point3d> pts = facePolylines[i].ToList();
                                pts.Add(pts[0]);
                                allBlockUnOffsetedArea.Branch(ghPathForUnOffsetedArea).Add(splitedAreas[j]);

                                allBlockPolylineDT.EnsurePath(ghPathForBlockPolyline);
                                allBlockPolylineDT.Branch(ghPathForBlockPolyline).Add(splitedPolylines[j]);

                                for (int k = 0; k < cellLoL.Count; k++)
                                {
                                    ghPathForCell = ghPathForCell.AppendElement(k);
                                    allBlockBasicCellDT.EnsurePath(ghPathForCell);
                                    allBlockBasicCellDT.Branch(ghPathForCell).AddRange(cellLoL[k]);
                                    ghPathForCell = ghPathForCell.CullElement();
                                }

                                for (int k = 0; k < hCurveLoL.Count; k++)
                                {
                                    ghPathForHCurve = ghPathForHCurve.AppendElement(k);
                                    for (int l = 0; l < hCurveLoL[k].Count; l++)
                                    {
                                        ghPathForHCurve = ghPathForHCurve.AppendElement(l);

                                        allBlockHCurvesDT.EnsurePath(ghPathForHCurve);
                                        allBlockHCurvesDT.Branch(ghPathForHCurve).AddRange(hCurveLoL[k][l]);

                                        ghPathForHCurve = ghPathForHCurve.CullElement();
                                    }
                                    ghPathForHCurve = ghPathForHCurve.CullElement();
                                }

                                for (int k = 0; k < vCurveLoL.Count; k++)
                                {
                                    ghPathForVCurve = ghPathForVCurve.AppendElement(k);
                                    for (int l = 0; l < vCurveLoL[k].Count; l++)
                                    {
                                        ghPathForVCurve = ghPathForVCurve.AppendElement(l);

                                        allBlockVCurvesDT.EnsurePath(ghPathForVCurve);
                                        allBlockVCurvesDT.Branch(ghPathForVCurve).AddRange(vCurveLoL[k][l]);

                                        ghPathForVCurve = ghPathForVCurve.CullElement();
                                    }
                                    ghPathForVCurve = ghPathForVCurve.CullElement();
                                }

                                for (int k = 0; k < contourLoL.Count; k++)
                                {
                                    ghPathForContour = ghPathForContour.AppendElement(k);
                                    for (int l = 0; l < contourLoL[k].Count; l++)
                                    {
                                        ghPathForContour = ghPathForContour.AppendElement(l);

                                        allBlockContourDT.EnsurePath(ghPathForContour);
                                        allBlockContourDT.Branch(ghPathForContour).AddRange(contourLoL[k][l]);

                                        ghPathForContour = ghPathForContour.CullElement();
                                    }
                                    ghPathForContour = ghPathForContour.CullElement();
                                }

                                for (int k = 0; k < holeLoL.Count; k++)
                                {
                                    ghPathForHole = ghPathForHole.AppendElement(k);
                                    for (int l = 0; l < holeLoL[k].Count; l++)
                                    {
                                        ghPathForHole = ghPathForHole.AppendElement(l);

                                        allBlockHoleDT.EnsurePath(ghPathForHole);
                                        allBlockHoleDT.Branch(ghPathForHole).AddRange(holeLoL[k][l]);

                                        ghPathForHole = ghPathForHole.CullElement();
                                    }
                                    ghPathForHole = ghPathForHole.CullElement();
                                }

                                for (int k = 0; k < cellBoundaryLoL.Count; k++)
                                {
                                    ghPathForCellBoundary = ghPathForCellBoundary.AppendElement(k);
                                    allBlockCellBoundaryDT.EnsurePath(ghPathForCellBoundary);
                                    allBlockCellBoundaryDT.Branch(ghPathForCellBoundary).AddRange(cellBoundaryLoL[k]);
                                    ghPathForCellBoundary = ghPathForCellBoundary.CullElement();
                                }
                                #endregion
                            }
                            else
                            {
                                #region output dataTree
                                allBlockUnOffsetedArea.EnsurePath(ghPathForUnOffsetedArea);
                                allBlockHCurvesDT.EnsurePath(ghPathForHCurve);
                                allBlockVCurvesDT.EnsurePath(ghPathForVCurve);
                                allBlockContourDT.EnsurePath(ghPathForContour);
                                allBlockHoleDT.EnsurePath(ghPathForHole);

                                allBlockCellBoundaryDT.EnsurePath(ghPathForCellBoundary);

                                allBlockBasicCellDT.EnsurePath(ghPathForCell);

                                allBlockPolylineDT.EnsurePath(ghPathForBlockPolyline);
                                #endregion
                            }

                            #region 清除当前第二层级的序号
                            ghPathForUnOffsetedArea = ghPathForUnOffsetedArea.CullElement();

                            ghPathForHCurve = ghPathForHCurve.CullElement();
                            ghPathForVCurve = ghPathForVCurve.CullElement();
                            ghPathForContour = ghPathForContour.CullElement();
                            ghPathForHole = ghPathForHole.CullElement();

                            ghPathForCellBoundary = ghPathForCellBoundary.CullElement();

                            ghPathForCell = ghPathForCell.CullElement();

                            ghPathForBlockPolyline = ghPathForBlockPolyline.CullElement();
                            #endregion
                        }

                    }
                }


                

                AllBlockCellBoundaryDT = allBlockCellBoundaryDT;

                //int additionalGenerations = allBlockBasicCellDT.Path(0).Length - 2;
                //DataTree<Cell> shiftedAllBlockCellDT = DataTreeShiftPath<Cell>(allBlockBasicCellDT, additionalGenerations);

                #region 生成每个Cell所对应的Brep[]，注意此时的cell中的brep是所有baseline构成的，口字或8字或I字或C字等，没有经过削减的
                for (int i = 0; i < allBlockBasicCellDT.BranchCount; i++)
                {
                    #region 生成Brep[]，并将 brep[] 添加到对应的Cell对象中
                    //List<List<List<Curve>>> curveForJoinLoL = new List<List<List<Curve>>>();
                    //List<Curve> allCurveForJoin = new List<Curve>();
                    for (int j = 0; j < allBlockBasicCellDT.Branch(i).Count; j++)
                    {
                        List<Curve> allCurves = new List<Curve>();
                        if (allBlockBasicCellDT.Branch(i)[j] is BilinearCell)
                        {
                            BilinearCell bCell = allBlockBasicCellDT.Branch(i)[j] as BilinearCell;
                            if (bCell.ShapeType == BilinearCell.BilinearShapeType.SingleRegion)
                            {
                                Brep[] brepArray = BoundarySurfaces(bCell.CellBoundary);
                                //if (brepArray != null)
                                //{
                                //    allBlockGroundBrepDT.Branch(i).AddRange(brepArray);
                                //}
                                allBlockBasicCellDT.Branch(i)[j].CellBreps = brepArray;
                            }
                            else
                            {
                                List<Curve> hCurves = allBlockBasicCellDT.Branch(i)[j].GetAllHorizontal();
                                List<Curve> vCurves = allBlockBasicCellDT.Branch(i)[j].GetAllVertical();
                                allCurves.AddRange(hCurves);
                                allCurves.AddRange(vCurves);

                                Curve[] joined = Curve.JoinCurves(allCurves);
                                //Curve[] trimedJoined = new Curve[joined.Length];
                                List<Curve> trimedJoined = new List<Curve>();
                                if (joined != null)
                                {
                                    for (int k = 0; k < joined.Length; k++)
                                    {
                                        //Vector3d vec = joined[k].TangentAtStart;
                                        //Vector3d depthVec =  new Vector3d(-vec.Y, vec.X, vec.Z);

                                        if (joined[k] != null)
                                        {
                                            if (!joined[k].IsClosed)
                                            {
                                                Curve crv = TrimBothEnd(allBlockBasicCellDT.Branch(i)[j].CellBoundary, joined[k], 0, w);
                                                trimedJoined.Add(crv);
                                            }
                                            else
                                            {
                                                trimedJoined.Add(joined[k]);
                                            }
                                        }
                                    }

                                    List<Polyline> outContours;
                                    List<Polyline> outHoles;
                                    GenerateBuilding(trimedJoined.ToArray(), w, out outContours, out outHoles);
                                    //GenerateBuilding(joined, w, out outContours, out outHoles);

                                    List<Curve> all = new List<Curve>();
                                    for (int k = 0; k < outContours.Count; k++)
                                    {
                                        all.Add(outContours[k].ToPolylineCurve());
                                    }
                                    for (int k = 0; k < outHoles.Count; k++)
                                    {
                                        all.Add(outHoles[k].ToPolylineCurve());
                                    }

                                    Brep[] brepArray = BoundarySurfaces(all);
                                    //if (brepArray != null)
                                    //{
                                    //    allBlockGroundBrepDT.Branch(i).AddRange(brepArray);
                                    //}
                                    allBlockBasicCellDT.Branch(i)[j].CellBreps = brepArray;
                                }
                            }
                        }
                        else
                        {
                            List<Curve> hCurves = allBlockBasicCellDT.Branch(i)[j].GetAllHorizontal();
                            List<Curve> vCurves = allBlockBasicCellDT.Branch(i)[j].GetAllVertical();
                            allCurves.AddRange(hCurves);
                            allCurves.AddRange(vCurves);

                            Curve[] joined = Curve.JoinCurves(allCurves);
                            //Curve[] trimedJoined = new Curve[joined.Length];
                            List<Curve> trimedJoined = new List<Curve>();
                            if (joined != null)
                            {
                                for (int k = 0; k < joined.Length; k++)
                                {
                                    //Vector3d vec = joined[k].TangentAtStart;
                                    //Vector3d depthVec =  new Vector3d(-vec.Y, vec.X, vec.Z);

                                    if (joined[k] != null)
                                    {
                                        if (!joined[k].IsClosed)
                                        {
                                            Curve crv = TrimBothEnd(allBlockBasicCellDT.Branch(i)[j].CellBoundary, joined[k], 0, w);
                                            trimedJoined.Add(crv);
                                        }
                                        else
                                        {
                                            trimedJoined.Add(joined[k]);
                                        }
                                    }
                                }

                                List<Polyline> outContours;
                                List<Polyline> outHoles;
                                GenerateBuilding(trimedJoined.ToArray(), w, out outContours, out outHoles);
                                //GenerateBuilding(joined, w, out outContours, out outHoles);

                                List<Curve> all = new List<Curve>();
                                for (int k = 0; k < outContours.Count; k++)
                                {
                                    all.Add(outContours[k].ToPolylineCurve());
                                }
                                for (int k = 0; k < outHoles.Count; k++)
                                {
                                    all.Add(outHoles[k].ToPolylineCurve());
                                }

                                Brep[] brepArray = BoundarySurfaces(all);
                                //if (brepArray != null)
                                //{
                                //    allBlockGroundBrepDT.Branch(i).AddRange(brepArray);
                                //}
                                allBlockBasicCellDT.Branch(i)[j].CellBreps = brepArray;
                            }
                        }




                    }

                    #endregion
                }
                #endregion

                #region 对allBlockBasicCellDT的每一个Branch，计算HAverageLength，VAverageLength，AverageLength，MaxFirstFloorArea
                DataTree<double> allBlockHAverageLength = new DataTree<double>();
                DataTree<double> allBlockVAverageLength = new DataTree<double>();
                DataTree<double> allBlockAverageLength = new DataTree<double>();
                DataTree<double> allBlockMaxFirstFloorArea = new DataTree<double>();
                DataTree<double> allBlockCellBoundaryAreaSum = new DataTree<double>();
                foreach (GH_Path path in allBlockBasicCellDT.Paths)
                {
                    allBlockHAverageLength.EnsurePath(new GH_Path(path[0],path[1]));
                    allBlockVAverageLength.EnsurePath(new GH_Path(path[0], path[1]));
                    allBlockAverageLength.EnsurePath(new GH_Path(path[0], path[1]));
                    allBlockMaxFirstFloorArea.EnsurePath(new GH_Path(path[0], path[1]));
                    allBlockCellBoundaryAreaSum.EnsurePath(new GH_Path(path[0], path[1]));

                    double smallBlockFirstFloorAreaSum = 0.0;
                    double smallBlockHLengthSum = 0;
                    double smallBlockVLengthSum = 0;
                    //double smallBlockAllLengthSum = 0;
                    double smallBlockCellBoundarySum = 0;
                    int hCountSum = 0;
                    int vCountSum = 0;
                    for (int i = 0; i < allBlockBasicCellDT.Branch(path).Count; i++)
                    {

                        for (int j = 0; j < allBlockBasicCellDT.Branch(path)[i].CellBreps.Length; j++)
                        {
                            smallBlockFirstFloorAreaSum += allBlockBasicCellDT.Branch(path)[i].CellBreps[j].GetArea();
                        }

                        int hCount;
                        smallBlockHLengthSum += allBlockBasicCellDT.Branch(path)[i].GetAllHCurveLength(out hCount);
                        int vCount;
                        smallBlockVLengthSum += allBlockBasicCellDT.Branch(path)[i].GetAllVCurveLength(out vCount);
                        hCountSum += hCount;
                        vCountSum += vCount;

                        smallBlockCellBoundarySum += AreaMassProperties.Compute(allBlockBasicCellDT.Branch(path)[i].CellBoundary).Area;
                    }

                    allBlockMaxFirstFloorArea.Branch(new GH_Path(path[0], path[1])).Add(smallBlockFirstFloorAreaSum);
                    allBlockHAverageLength.Branch(new GH_Path(path[0], path[1])).Add(smallBlockHLengthSum / hCountSum);
                    allBlockVAverageLength.Branch(new GH_Path(path[0], path[1])).Add(smallBlockVLengthSum / vCountSum);
                    allBlockAverageLength.Branch(new GH_Path(path[0], path[1])).Add((smallBlockHLengthSum + smallBlockVLengthSum) / (hCountSum + vCountSum));
                    allBlockCellBoundaryAreaSum.Branch(new GH_Path(path[0], path[1])).Add(smallBlockCellBoundarySum);
                }

                #endregion

                foreach (GH_Path path in allBlockMaxFirstFloorArea.Paths)
                {
                    double h = 0;
                    for (int i = 0; i < allBlockHAverageLength.Branch(path).Count; i++)
                    {
                        h += allBlockHAverageLength.Branch(path)[i];
                    }
                    h /= allBlockHAverageLength.Branch(path).Count;
                    allBlockHAverageLength.Branch(path).Clear();
                    allBlockHAverageLength.Branch(path).Add(h);

                    double v = 0;
                    for (int i = 0; i < allBlockVAverageLength.Branch(path).Count; i++)
                    {
                        v += allBlockVAverageLength.Branch(path)[i];
                    }
                    v /= allBlockVAverageLength.Branch(path).Count;
                    allBlockVAverageLength.Branch(path).Clear();
                    allBlockVAverageLength.Branch(path).Add(v);

                    double average = 0;
                    for (int i = 0; i < allBlockAverageLength.Branch(path).Count; i++)
                    {
                        average += allBlockAverageLength.Branch(path)[i];
                    }
                    average /= allBlockAverageLength.Branch(path).Count;
                    allBlockAverageLength.Branch(path).Clear();
                    allBlockAverageLength.Branch(path).Add(average);

                    double maxA = 0;
                    for (int i = 0; i < allBlockMaxFirstFloorArea.Branch(path).Count; i++)
                    {
                        maxA += allBlockMaxFirstFloorArea.Branch(path)[i];
                    }
                    // maxA /= allBlockMaxFirstFloorArea.Branch(path).Count;
                    allBlockMaxFirstFloorArea.Branch(path).Clear();
                    allBlockMaxFirstFloorArea.Branch(path).Add(maxA);

                    double sum = 0;
                    for (int i = 0; i < allBlockCellBoundaryAreaSum.Branch(path).Count; i++)
                    {
                        sum += allBlockCellBoundaryAreaSum.Branch(path)[i];
                    }
                    allBlockCellBoundaryAreaSum.Branch(path).Clear();
                    allBlockCellBoundaryAreaSum.Branch(path).Add(sum);
                }

                #region 进行建筑密度约束
                DataTree<double> allBlockArea = new DataTree<double>();
                //DataTree<double> allBlockExpectedFirstFloorArea = new DataTree<double>();
                //DataTree<double> allBlockExpectedBuildingArea = new DataTree<double>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    allBlockArea.EnsurePath(path);
                    //allBlockExpectedFirstFloorArea.EnsurePath(path);
                    //allBlockExpectedBuildingArea.EnsurePath(path);

                    AreaMassProperties blockMass = AreaMassProperties.Compute(facePolylines[path[0]].ToNurbsCurve());
                    allBlockArea.Branch(path).Add(blockMass.Area);
                    //allBlockExpectedFirstFloorArea.Branch(path).Add(densitys[path[0]] * allBlockUnOffsetedArea.Branch(path)[0]);
                    //allBlockExpectedBuildingArea.Branch(path).Add(areaRatio[path[0]] * allBlockUnOffsetedArea.Branch(path)[0]);
                }


                // todo 基于退线后的用地面积的比例，来修改期望建筑面积的分配
                List<double> expectedBuildingArea = new List<double>();
                for (int i = 0; i < areaRatio.Count; i++)
                {
                    expectedBuildingArea.Add(areaRatio[i] * AreaMassProperties.Compute(facePolylines[i].ToNurbsCurve()).Area);
                }
                List<double> expectedFirstFloorArea = new List<double>();
                for (int i = 0; i < areaRatio.Count; i++)
                {
                    expectedFirstFloorArea.Add(densitys[i] * AreaMassProperties.Compute(facePolylines[i].ToNurbsCurve()).Area);
                }

                DataTree<double> allBlockSmallPolylineArea = new DataTree<double>();
                DataTree<double> allBlockPolylineArea = new DataTree<double>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    allBlockSmallPolylineArea.EnsurePath(path);
                    AreaMassProperties blockMass = AreaMassProperties.Compute(allBlockPolylineDT.Branch(path)[0].ToNurbsCurve());
                    allBlockSmallPolylineArea.Branch(path).Add(blockMass.Area);

                    allBlockPolylineArea.EnsurePath(path);
                    AreaMassProperties blockMass2 = AreaMassProperties.Compute(innerResultPolylines[path[0]].ToNurbsCurve());
                    allBlockPolylineArea.Branch(path).Add(blockMass2.Area);
                }

                DataTree<double> allBlockExpectedFirstFloorArea = new DataTree<double>();
                DataTree<double> allBlockExpectedBuildingArea = new DataTree<double>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    allBlockExpectedFirstFloorArea.EnsurePath(path);
                    allBlockExpectedBuildingArea.EnsurePath(path);

                    allBlockExpectedFirstFloorArea.Branch(path).Add(expectedFirstFloorArea[path[0]] * allBlockSmallPolylineArea.Branch(path)[0] / allBlockPolylineArea.Branch(path)[0]);
                    allBlockExpectedBuildingArea.Branch(path).Add(expectedBuildingArea[path[0]] * allBlockSmallPolylineArea.Branch(path)[0] / allBlockPolylineArea.Branch(path)[0]);
                }


                #region 计算每个小block所能生成的最大建筑密度与输入给定的期望建筑密度之间的差值
                DataTree<double> allBlockDelta = new DataTree<double>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    allBlockDelta.EnsurePath(path);
                    if (allBlockMaxFirstFloorArea.Branch(path)[0] < allBlockExpectedFirstFloorArea.Branch(path)[0])
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处所给定的建筑密度，大于所能生成的最大建筑密度！", path.ToString()));
                        // 此时要全部生成 口字 或 8字 
                        allBlockDelta.Branch(path).Add(0);
                    }
                    else
                    {
                        // 此时要修剪
                        allBlockDelta.Branch(path).Add(allBlockMaxFirstFloorArea.Branch(path)[0] - allBlockExpectedFirstFloorArea.Branch(path)[0]);
                    }
                }
                #endregion

                #region 找到属于每个大Block的所有分支 List<List<GH_Path>> pathsBelongsToSameBlock
                DataTree<GH_Path> pathsBelongsToSameBlock = new DataTree<GH_Path>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    pathsBelongsToSameBlock.EnsurePath(path);
                    foreach (var path1 in allBlockBasicCellDT.Paths)
                    {
                        if (path[0] == path1[0] && path[1] == path1[1])
                        {
                            pathsBelongsToSameBlock.Branch(path).Add(path1);
                        }
                    }
                }
                #endregion

                int additionalGenerations = -(allBlockBasicCellDT.Path(0).Length - 2);
                //DataTree<double> allBlockHAverageLength = DataTreeShiftPath<double>(allBlockHAverageLength, additionalGenerations);
                DataTree<Cell> shiftedAllBlockBasicCellDT = DataTreeShiftPath<Cell>(allBlockBasicCellDT, additionalGenerations);

                // 深拷贝CellDT
                DataTree<Cell> allBlockCellDT = new DataTree<Cell>();
                foreach (var path in allBlockBasicCellDT.Paths)
                {
                    allBlockCellDT.EnsurePath(path);
                    foreach (var cell in allBlockBasicCellDT.Branch(path))
                    {
                        if (cell is BilinearCell)
                        {
                            BilinearCell bCell = cell as BilinearCell;
                            allBlockCellDT.Branch(path).Add(new BilinearCell(bCell));
                        }
                        else
                        {
                            TrilinearCell tCell = cell as TrilinearCell;
                            allBlockCellDT.Branch(path).Add(new TrilinearCell(tCell));
                        }
                    }
                }

                #region 整理每个小Block，需要的vReduceCount，offset
                // 整理每个小Block，需要的vReduceCount，offset，每个小Block各自为战，对应shiftedAllBlockExpectedFirstFloorArea即可
                //DataTree<int> vReduceCountDT = new DataTree<int>();
                //DataTree<int> hReduceCountDT = new DataTree<int>();
                //DataTree<double> scaleFactorDT = new DataTree<double>();
                foreach (var path in allBlockPolylineDT.Paths)
                {
                    //vReduceCountDT.EnsurePath(path);
                    //hReduceCountDT.EnsurePath(path);
                    //scaleFactorDT.EnsurePath(path);

                    /* 因为是各自为战，所以不用管当前这个大Block中包不包含小block */

                    //if (pathsBelongsToSameBlock[path[0]].Count == 1)
                    //{
                    //    // 当前这个大Block中不存在小Block
                    //}
                    //else
                    //{
                    //    // 当前这个大Block中存在小Block
                    //}

                    /* 只需先判断当前这个小Block中Cell的情况 */
                    int singleRegionCount = 0;
                    int singlePolylineCount = 0;
                    int IShapeCount = 0;
                    int CShapeCount = 0;
                    int RecShapeCount = 0;
                    int EightShapeCount = 0;
                    //int elseCount = 0;
                    for (int i = 0; i < shiftedAllBlockBasicCellDT.Branch(path).Count; i++)
                    {
                        if (shiftedAllBlockBasicCellDT.Branch(path)[i] is BilinearCell)
                        {
                            // 这个Cell是否是BilinearCell
                            BilinearCell bCell = shiftedAllBlockBasicCellDT.Branch(path)[i] as BilinearCell;
                            if (bCell.ShapeType == BilinearCell.BilinearShapeType.SingleRegion)
                            {
                                singleRegionCount++;
                            }
                            else if (bCell.ShapeType == BilinearCell.BilinearShapeType.IShape)
                            {
                                IShapeCount++;
                            }
                            else if (bCell.ShapeType == BilinearCell.BilinearShapeType.CShape)
                            {
                                CShapeCount++;
                            }
                            else if (bCell.ShapeType == BilinearCell.BilinearShapeType.SinglePolyline)
                            {
                                singlePolylineCount++;
                            }
                            else if (bCell.ShapeType == BilinearCell.BilinearShapeType.RecShape)
                            {
                                RecShapeCount++;
                            }
                            else
                            {
                                RecShapeCount++;
                            }
                        }
                        else
                        {
                            EightShapeCount++;
                        }
                    }

                    if (singleRegionCount == shiftedAllBlockBasicCellDT.Branch(path).Count)
                    {
                        // 全部都是 singleRegion 时
                        // 此时不计算VReduceCount，计算scaleFactor
                        if (singleRegionCount != 1)
                        {
                            // 应该不可能
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("Path为 {0} 的小Block内,存在多个singleRegion的情况，该部分的代码未完善", path.ToString()));
                        }
                        else
                        {
                            //// 不计算VReduceCount
                            //vReduceCountDT.Branch(path).Add(-1);
                            //// 不计算hReduceCount
                            //hReduceCountDT.Branch(path).Add(-1);
                            // 计算scaleFactor
                            double scaleFactor = (allBlockMaxFirstFloorArea.Branch(path)[0] - allBlockDelta.Branch(path)[0]) / allBlockMaxFirstFloorArea.Branch(path)[0];
                            if (scaleFactor > 1)
                            {
                                scaleFactor = 1;
                            }
                            allBlockCellDT = DoDensityReduce_Scale(allBlockCellDT, path, scaleFactor, lMin, w);
                        }
                    }
                    else if (singlePolylineCount == shiftedAllBlockBasicCellDT.Branch(path).Count)
                    {
                        // 此时只有一个 而且它是 singlePolyline
                        if (singlePolylineCount != 1)
                        {
                            // 应该不可能
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("Path为 {0} 的小Block内,存在多个singlePolyline的情况，该部分的代码未完善", path.ToString()));
                        }
                        else
                        {
                            BilinearCell bCell = allBlockCellDT.Branch(new GH_Path(path[0], path[1], 0))[0] as BilinearCell;

                            if (bCell.AnotherBaseLineIndexRelatedToCutPoint == 0 || bCell.AnotherBaseLineIndexRelatedToCutPoint == 2)
                            {
                                // 此时先减v，再减h
                                int vCount = (int)(Math.Ceiling(allBlockDelta.Branch(path)[0] / allBlockVAverageLength.Branch(path)[0]));
                                if (vCount > singlePolylineCount)
                                {
                                    vCount = singlePolylineCount;
                                    double area = allBlockDelta.Branch(path)[0] - vCount * allBlockVAverageLength.Branch(path)[0];
                                    int hCount = 0;
                                    if (area <= 0)
                                    {
                                        hCount = 0;
                                    }
                                    else
                                    {
                                        hCount = (int)Math.Ceiling(area / allBlockHAverageLength.Branch(path)[0]);
                                    }
                                    if (hCount > singlePolylineCount)
                                    {
                                        double scale = allBlockExpectedFirstFloorArea.Branch(path)[0] / allBlockCellBoundaryAreaSum.Branch(path)[0];
                                        allBlockCellDT = DoDensityReduce_Scale(allBlockCellDT, path, scale, lMin, w);
                                    }
                                    else
                                    {
                                        // v直接切除，并且将原来的h变短
                                        hCount = 0;
                                        double sH = area / allBlockHAverageLength.Branch(path)[0];
                                        double sV = 0;
                                        if (bCell.PrevBaseLineIndexRelatedToCutPoint == 3 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                        {
                                            // 从v的0端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 0, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint == 1 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                        {
                                            // 从v的1端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 2, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 1)
                                        {
                                            // 从v的1端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 0, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 3)
                                        {
                                            // 从v的0端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 2, true);
                                        }
                                    }
                                }
                                else
                                {
                                    // 将原来的v变短
                                    vCount = 0;
                                    int hCount = 0;
                                    double sH = 0;
                                    double sV = allBlockDelta.Branch(path)[0] / allBlockVAverageLength.Branch(path)[0];
                                    if (bCell.PrevBaseLineIndexRelatedToCutPoint == 1 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                    {
                                        // 从v的0端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 1, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint == 3 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                    {
                                        // 从v的1端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 3, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 3)
                                    {
                                        // 从v的1端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 3, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 1)
                                    {
                                        // 从v的0端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 1, true);
                                    }
                                }
                            }
                            else
                            {
                                // bCell.AnotherBaseLineIndexRelatedToCutPoint == 1 || bCell.AnotherBaseLineIndexRelatedToCutPoint == 3
                                // 此时先减h，再减v
                                // 即，先判断减 S 还是减 N ，再减West
                                int hCount = (int)(Math.Ceiling(allBlockDelta.Branch(path)[0] / allBlockHAverageLength.Branch(path)[0]));
                                if (hCount > singlePolylineCount)
                                {
                                    hCount = singlePolylineCount;
                                    double area = allBlockDelta.Branch(path)[0] - hCount * allBlockHAverageLength.Branch(path)[0];
                                    int vCount = (int)Math.Ceiling(area / allBlockVAverageLength.Branch(path)[0]);
                                    if (vCount > singlePolylineCount)
                                    {
                                        double scale = allBlockExpectedFirstFloorArea.Branch(path)[0] / allBlockCellBoundaryAreaSum.Branch(path)[0];
                                        allBlockCellDT = DoDensityReduce_Scale(allBlockCellDT, path, scale, lMin, w);
                                    }
                                    else
                                    {
                                        // h直接切除，并且将原来的v变短
                                        vCount = 0;
                                        double sV = area / allBlockVAverageLength.Branch(path)[0];
                                        double sH = 0;
                                        if (bCell.PrevBaseLineIndexRelatedToCutPoint == 0 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                        {
                                            // 从v的0端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 1, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint == 2 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                        {
                                            // 从v的1端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 3, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 2)
                                        {
                                            // 从v的1端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 1, true);
                                        }
                                        else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 0)
                                        {
                                            // 从v的0端开始减
                                            allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 3, true);
                                        }
                                    }
                                }
                                else
                                {
                                    // 将原来的h变短
                                    hCount = 0;
                                    int vCount = 0;
                                    double sV = 0;
                                    double sH = allBlockDelta.Branch(path)[0] / allBlockHAverageLength.Branch(path)[0];
                                    if (bCell.PrevBaseLineIndexRelatedToCutPoint == 2 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                    {
                                        // 从h的0端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 2, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint == 0 && bCell.NextBaseLineIndexRelatedToCutPoint == -1)
                                    {
                                        // 从v的1端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 0, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 0)
                                    {
                                        // 从v的1端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 0, true);
                                    }
                                    else if (bCell.PrevBaseLineIndexRelatedToCutPoint != 1 && bCell.NextBaseLineIndexRelatedToCutPoint == 2)
                                    {
                                        // 从v的0端开始减
                                        allBlockCellDT = DoDensityReduce_SP(allBlockCellDT, path, vCount, hCount, sV, sH, 2, true);
                                    }
                                }
                            }

                        }
                    }
                    else if (singleRegionCount + singlePolylineCount + IShapeCount + CShapeCount == shiftedAllBlockBasicCellDT.Branch(path).Count
                             && singleRegionCount != 0)
                    {
                        // 应该不可能
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("Path为 {0} 的小Block内,存在singleRegion与singlePolyline(或者可能是IShape,CShape)共存的情况，该部分的代码未完善", path.ToString()));
                    }
                    else if (IShapeCount + CShapeCount == shiftedAllBlockBasicCellDT.Branch(path).Count)
                    {
                        // 全部是 singlePolyline 或 IShape 或 CShape 时
                        // 此时先，判断当前形状减去 vReduceCount 个V后，是否仍比期望的首层面积 shiftedAllBlockExpectedFirstFloorArea 大
                        // 如果是的话，就模仿 SingleRegion 进行 offset，不过这个offset操作，与 SingleRegion 的 offset 操作是反的

                        int vCountOnce = (int)(Math.Ceiling(allBlockDelta.Branch(path)[0] / allBlockVAverageLength.Branch(path)[0]));
                        // 注意此时存在的IShape一定是有两条V的，CShape也是一定是有两条V的
                        if (vCountOnce > IShapeCount * 2 + CShapeCount * 2)
                        {
                            // 不用vCountTwice
                            vCountOnce = IShapeCount * 2 + CShapeCount * 2;
                            double area = allBlockDelta.Branch(path)[0] - vCountOnce * allBlockVAverageLength.Branch(path)[0];
                            int hCount = 0;
                            if (area <= 0)
                            {
                                hCount = 0;
                            }
                            else
                            {
                                hCount = (int)Math.Ceiling(area / allBlockHAverageLength.Branch(path)[0]);
                            }
                            if (hCount > IShapeCount * 0 + CShapeCount * 0)
                            {
                                double scale = allBlockExpectedFirstFloorArea.Branch(path)[0] / allBlockCellBoundaryAreaSum.Branch(path)[0];
                                allBlockCellDT = DoDensityReduce_Scale(allBlockCellDT, path, scale, lMin, w);
                            }
                            else
                            {
                                allBlockCellDT = DoDensityReduce_IOrC(allBlockCellDT, path, vCountOnce, hCount, IShapeCount, CShapeCount);
                            }
                        }
                        else
                        {
                            allBlockCellDT = DoDensityReduce_IOrC(allBlockCellDT, path, vCountOnce, 0, IShapeCount, CShapeCount);
                        }
                    }
                    else if (IShapeCount + CShapeCount + RecShapeCount + EightShapeCount == shiftedAllBlockBasicCellDT.Branch(path).Count
                             && IShapeCount + CShapeCount != 0)
                    {
                        // 应该不可能
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("Path为 {0} 的小Block内，存在else与singlePolyline(或者可能是IShape,CShape)共存的情况，该部分的代码未完善", path.ToString()));
                    }
                    else
                    {
                        // 全部都是 elseCount 时
                        // 此时计算VReduceCount，不计算OffsetValue

                        // 计算VReduceCount
                        int vCountOnce = (int)(Math.Ceiling(allBlockDelta.Branch(path)[0] / allBlockVAverageLength.Branch(path)[0]));
                        if (vCountOnce > RecShapeCount * 1 + 2 * EightShapeCount)
                        {
                            int vCountTwice = vCountOnce - (RecShapeCount * 1 + 2 * EightShapeCount);
                            vCountOnce = RecShapeCount * 1 + 2 * EightShapeCount;

                            if (vCountTwice > RecShapeCount * 1 + 2 * EightShapeCount)
                            {
                                vCountTwice = RecShapeCount * 1 + 2 * EightShapeCount;

                                double area = allBlockDelta.Branch(path)[0] - (vCountOnce + vCountTwice) * allBlockVAverageLength.Branch(path)[0];
                                int hCount = 0;
                                if (area <= 0)
                                {
                                    hCount = 0;
                                }
                                else
                                {
                                    hCount = (int)Math.Ceiling(area / allBlockHAverageLength.Branch(path)[0]);
                                }
                                if (hCount > RecShapeCount * 1 + 2 * EightShapeCount)
                                {
                                    double scale = allBlockExpectedFirstFloorArea.Branch(path)[0] / allBlockCellBoundaryAreaSum.Branch(path)[0];
                                    allBlockCellDT = DoDensityReduce_Scale(allBlockCellDT, path, scale, lMin, w);
                                }
                                else
                                {
                                    allBlockCellDT = DoDensityReduce_RecOrEight(allBlockCellDT, path, vCountOnce, vCountTwice, hCount, RecShapeCount, EightShapeCount, lMin, w);
                                }
                            }
                            else
                            {
                                allBlockCellDT = DoDensityReduce_RecOrEight(allBlockCellDT, path, vCountOnce, vCountTwice, 0, RecShapeCount, EightShapeCount, lMin, w);
                            }
                        }
                        else
                        {
                            allBlockCellDT = DoDensityReduce_RecOrEight(allBlockCellDT, path, vCountOnce, 0, 0, RecShapeCount, EightShapeCount, lMin, w);
                        }
                    }
                }
                #endregion
                #endregion

                #region 满足容积率
                #region 计算当前建筑首层面积
                DataTree<Brep> allBlockGroundBrepDT = new DataTree<Brep>();
                #region 对于每个branch上的Cell，按照CellBoundary是否相交，来决定是否将这两个Cell的所有H,V线进行join计算
                for (int i = 0; i < allBlockCellDT.BranchCount; i++)
                {
                    allBlockGroundBrepDT.EnsurePath(allBlockCellDT.Paths[i]);

                    List<List<Curve>> hCurveLoL = new List<List<Curve>>();
                    List<List<Curve>> vCurveLoL = new List<List<Curve>>();

                    // 正常进行裁切和生成brep
                    #region 生成building的底面
                    for (int j = 0; j < allBlockCellDT.Branch(allBlockCellDT.Paths[i]).Count; j++)
                    {
                        List<Curve> allCurves = new List<Curve>();
                        if (allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j] is BilinearCell)
                        {
                            BilinearCell bCell = allBlockCellDT.Branch(i)[j] as BilinearCell;
                            if (bCell.ShapeType == BilinearCell.BilinearShapeType.SingleRegion)
                            {
                                Brep[] brepArray = BoundarySurfaces(bCell.CellBoundary);
                                if (brepArray != null)
                                {
                                    allBlockGroundBrepDT.Branch(allBlockCellDT.Paths[i]).AddRange(brepArray);
                                }
                            }
                            else if (bCell.ShapeType == BilinearCell.BilinearShapeType.Scaled)
                            {
                                Brep[] brepArray = BoundarySurfaces(bCell.FinalRegion);
                                if (brepArray != null)
                                {
                                    allBlockGroundBrepDT.Branch(allBlockCellDT.Paths[i]).AddRange(brepArray);
                                }
                            }
                            else
                            {
                                List<Curve> hCurves = allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].GetAllHorizontal();
                                List<Curve> vCurves = allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].GetAllVertical();
                                allCurves.AddRange(hCurves);
                                allCurves.AddRange(vCurves);

                                Curve[] joined = Curve.JoinCurves(allCurves);
                                //Curve[] trimedJoined = new Curve[joined.Length];
                                List<Curve> trimedJoined = new List<Curve>();
                                if (joined != null)
                                {
                                    for (int k = 0; k < joined.Length; k++)
                                    {
                                        //Vector3d vec = joined[k].TangentAtStart;
                                        //Vector3d depthVec =  new Vector3d(-vec.Y, vec.X, vec.Z);

                                        if (joined[k] != null)
                                        {
                                            if (!joined[k].IsClosed)
                                            {
                                                Curve crv = TrimBothEnd(allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].CellBoundary, joined[k], 0, w);
                                                trimedJoined.Add(crv);
                                            }
                                            else
                                            {
                                                trimedJoined.Add(joined[k]);
                                            }
                                        }
                                    }

                                    List<Polyline> outContours;
                                    List<Polyline> outHoles;
                                    GenerateBuilding(trimedJoined.ToArray(), w, out outContours, out outHoles);
                                    //GenerateBuilding(joined, w, out outContours, out outHoles);

                                    List<Curve> all = new List<Curve>();
                                    for (int k = 0; k < outContours.Count; k++)
                                    {
                                        all.Add(outContours[k].ToPolylineCurve());
                                    }
                                    for (int k = 0; k < outHoles.Count; k++)
                                    {
                                        all.Add(outHoles[k].ToPolylineCurve());
                                    }

                                    Brep[] brepArray = BoundarySurfaces(all);
                                    if (brepArray != null)
                                    {
                                        allBlockGroundBrepDT.Branch(allBlockCellDT.Paths[i]).AddRange(brepArray);
                                    }
                                }
                            }
                        }
                        else
                        {
                            TrilinearCell tCell = allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j] as TrilinearCell;

                            if (tCell.ShapeType == TrilinearCell.TrilinearShapeType.Scaled)
                            {
                                for (int k = 0; k < tCell.FinalRegions.Length; k++)
                                {
                                    Brep[] brepArray = BoundarySurfaces(tCell.FinalRegions[k]);
                                    if (brepArray != null)
                                    {
                                        allBlockGroundBrepDT.Branch(allBlockCellDT.Paths[i]).AddRange(brepArray);
                                    }
                                }

                                //Brep[] brepArray = BoundarySurfaces(tCell.FinalRegions.ToList());
                                //if (brepArray != null)
                                //{
                                //    allBlockGroundBrepDT.Branch(i).AddRange(brepArray);
                                //}
                            }
                            else
                            {
                                List<Curve> hCurves = allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].GetAllHorizontal();
                                List<Curve> vCurves = allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].GetAllVertical();
                                allCurves.AddRange(hCurves);
                                allCurves.AddRange(vCurves);

                                Curve[] joined = Curve.JoinCurves(allCurves);
                                //Curve[] trimedJoined = new Curve[joined.Length];
                                List<Curve> trimedJoined = new List<Curve>();
                                if (joined != null)
                                {
                                    for (int k = 0; k < joined.Length; k++)
                                    {
                                        //Vector3d vec = joined[k].TangentAtStart;
                                        //Vector3d depthVec =  new Vector3d(-vec.Y, vec.X, vec.Z);

                                        if (joined[k] != null)
                                        {
                                            if (!joined[k].IsClosed)
                                            {
                                                Curve crv = TrimBothEnd(allBlockCellDT.Branch(allBlockCellDT.Paths[i])[j].CellBoundary, joined[k], 0, w);
                                                trimedJoined.Add(crv);
                                            }
                                            else
                                            {
                                                trimedJoined.Add(joined[k]);
                                            }
                                        }
                                    }

                                    List<Polyline> outContours;
                                    List<Polyline> outHoles;
                                    GenerateBuilding(trimedJoined.ToArray(), w, out outContours, out outHoles);
                                    //GenerateBuilding(joined, w, out outContours, out outHoles);

                                    List<Curve> all = new List<Curve>();
                                    for (int k = 0; k < outContours.Count; k++)
                                    {
                                        all.Add(outContours[k].ToPolylineCurve());
                                    }
                                    for (int k = 0; k < outHoles.Count; k++)
                                    {
                                        all.Add(outHoles[k].ToPolylineCurve());
                                    }

                                    Brep[] brepArray = BoundarySurfaces(all);
                                    if (brepArray != null)
                                    {
                                        allBlockGroundBrepDT.Branch(allBlockCellDT.Paths[i]).AddRange(brepArray);
                                    }
                                }
                            }
                        }
                    }

                    #endregion
                }
                #endregion

                DataTree<double> allBlockCurrentFirstFloorArea = new DataTree<double>();
                foreach (var path in allBlockCellDT.Paths)
                {
                    allBlockCurrentFirstFloorArea.EnsurePath(new GH_Path(path[0], path[1]));

                    double areaSum = 0;
                    for (int i = 0; i < allBlockCellDT.Branch(path).Count; i++)
                    {
                        allBlockCellDT.Branch(path)[i].CellBreps = allBlockGroundBrepDT.Branch(path).ToArray();

                        for (int j = 0; j < allBlockCellDT.Branch(path)[i].CellBreps.Length; j++)
                        {
                            areaSum += allBlockCellDT.Branch(path)[i].CellBreps[j].GetArea();
                        }
                    }
                    allBlockCurrentFirstFloorArea.Branch(new GH_Path(path[0], path[1])).Add(areaSum);
                }

                foreach (var path in allBlockCurrentFirstFloorArea.Paths)
                {
                    double areaSum = 0;
                    for (int i = 0; i < allBlockCurrentFirstFloorArea.Branch(path).Count; i++)
                    {
                        areaSum += allBlockCurrentFirstFloorArea.Branch(path)[i];
                    }

                    allBlockCurrentFirstFloorArea.Branch(path).Clear();
                    allBlockCurrentFirstFloorArea.Branch(path).Add(areaSum);
                }
                #endregion

                // 每个小block的整体层数的下限
                DataTree<int> allBlockLowestLimitLayerNum = new DataTree<int>();
                DataTree<double> allBlockDeltaBuildingArea = new DataTree<double>();
                foreach (GH_Path path in allBlockPolylineDT.Paths)
                {
                    allBlockLowestLimitLayerNum.EnsurePath(path);
                    allBlockDeltaBuildingArea.EnsurePath(path);
                    int layerNum = (int)Math.Ceiling(allBlockExpectedBuildingArea.Branch(path)[0] / allBlockCurrentFirstFloorArea.Branch(path)[0]);
                    double area = allBlockExpectedBuildingArea.Branch(path)[0] - layerNum * allBlockCurrentFirstFloorArea.Branch(path)[0];
                    allBlockLowestLimitLayerNum.Branch(path).Add(layerNum);
                    allBlockDeltaBuildingArea.Branch(path).Add(area);
                }
                #endregion


                DataTree<Brep> allBlockBrepDT = new DataTree<Brep>();
                foreach (var path in allBlockCellDT.Paths)
                {
                    for (int i = 0; i < allBlockLowestLimitLayerNum.Branch(new GH_Path(path[0],path[1]))[0]; i++)
                    {
                        GH_Path pathForallBlockBrepDT = new GH_Path(path[0], path[1], path[2], i);

                        allBlockBrepDT.EnsurePath(pathForallBlockBrepDT);

                        List<Brep> resultBrps = new List<Brep>();
                        Vector3d dir = new Vector3d(0, 0, floorHeight);
                        for (int j = 0; j < allBlockCellDT.Branch(path).Count; j++)
                        {
                            //Brep[] extrudeBreps = new Brep[allBlockCellDT.Branch(path)[i].CellBreps.Length];
                            for (int k = 0; k < allBlockCellDT.Branch(path)[j].CellBreps.Length; k++)
                            {
                                Brep brp = allBlockCellDT.Branch(path)[j].CellBreps[k];
                                resultBrps.Add(ExtrudeBrep(brp, dir));
                            }
                        }

                        Vector3d translateDir = new Vector3d(0, 0, i * floorHeight);
                        //List<Brep> translatedResultBrps = new List<Brep>();
                        for (int j = 0; j < resultBrps.Count; j++)
                        {
                            resultBrps[j].Translate(translateDir);
                        }

                        allBlockBrepDT.Branch(pathForallBlockBrepDT).AddRange(resultBrps);
                    }
                    
                }




                //DataTree<Polyline> contourDT = UtilityFunctions.LoLToDataTree<Polyline>(contourLoL);
                //DataTree<Polyline> holeDT = UtilityFunctions.LoLToDataTree<Polyline>(holeLoL);
                DA.SetDataTree(3, allBlockCellBoundaryDT);
                DA.SetDataTree(4, allBlockHoleDT);

                //DataTree<Brep> brepDT = UtilityFunctions.LoLToDataTree<Brep>(brepLoL);
                DA.SetDataTree(2, allBlockBrepDT);

                //DataTree<Curve> allHCurvesDT = UtilityFunctions.LoLToDataTree<Curve>(allHCurvesLoL);
                //DataTree<Curve> allVCurvesDT = UtilityFunctions.LoLToDataTree<Curve>(allVCurvesLoL);
                DA.SetDataTree(5, allBlockHCurvesDT);
                DA.SetDataTree(6, allBlockVCurvesDT);

                //DA.SetDataTree(7, allBlockCellBoundaryDT2);
                //DA.SetDataTree(8, allBlockCellBoundaryDT);
                #region 可视化
                InnerResultPolyline.Clear();
                InnerNodeTextDots.Clear();
                InnerResultPolyline = innerResultPolylines;

                for (int i = 0; i < graphNodes.Count; i++)
                {
                    if (graphNodes[i].IsInner)
                    {
                        TextDot textDot = new TextDot(string.Format("{0} | {1}", i, graphNodes[i].NodeAttribute.NodeLabel), innerResultPolylines[dFIndex_PVIndex.IndexOf(i)].CenterPoint());
                        InnerNodeTextDots.Add(textDot);
                    }
                }

                #endregion
            }
        }

        private Brep ExtrudeBrep(Brep brp, Vector3d dir)
        {
            Brep result;
            if (dir.IsZero)
            {
                result = brp;
            }
            else
            {
                List<Curve> list = new List<Curve>();
                int num = brp.Edges.Count;
                for (int i = 0; i < num; i++)
                {
                    if (brp.Edges[i].TrimCount == 1)
                    {
                        Curve item = brp.Edges[i].DuplicateCurve();
                        list.Add(item);
                    }
                }
                if (list.Count == 0)
                {
                    result = brp;
                }
                else
                {
                    List<Brep> list2 = new List<Brep>();
                    list2.Add(brp);
                    int num2 = list.Count;
                    for (int j = 0; j < num2; j++)
                    {
                        Surface surface = Surface.CreateExtrusion(list[j], dir);
                        if (surface != null)
                        {
                            Brep brep = surface.ToBrep();
                            brep.Faces.SplitKinkyFaces();
                            list2.Add(brep);
                        }
                    }
                    Brep brep2 = brp.DuplicateBrep();
                    brep2.Translate(dir);
                    list2.Add(brep2);
                    Brep[] array = Brep.JoinBreps(list2, 0.0001);
                    if (array == null)
                    {
                        result = brp;
                    }
                    else
                    {
                        result = array[0];
                    }
                }
            }
            return result;
        }



        /// <summary>
        /// -1 是从path的后开始减，1 是从path的前开始减
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="sourceDT"></param>
        /// <param name="additionalGenerations"></param>
        /// <returns></returns>
        private DataTree<T> DataTreeShiftPath<T>(DataTree<T> sourceDT, int additionalGenerations)
        {
            int num = additionalGenerations;

            DataTree<T> outDataTree = new DataTree<T>();

            int num2 = sourceDT.Paths.Count - 1;
            int i = 0;
            while (i <= num2)
            {
                GH_Path gh_Path = sourceDT.Path(i);
                List<T> data = sourceDT.Branch(i);
                if (num == 0)
                {
                    goto IL_A9;
                }
                int num3 = Math.Abs(num);
                for (int j = 0; j < num3; j++)
                {
                    if (num < 0)
                    {
                        gh_Path = gh_Path.CullElement();
                    }
                    else
                    {
                        gh_Path = gh_Path.CullFirstElement();
                    }
                }
                if (gh_Path.Length != 0)
                {
                    goto IL_A9;
                }
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, string.Format("Path {0} was shifted out of existence", outDataTree.Path(i)));

            IL_A9:
                outDataTree.AddRange(data, gh_Path);
                goto IL_B3;
            IL_B3:
                i++;
                continue;
            }

            return outDataTree;
        }



        /// <summary>
        /// 此时不计算VReduceCount，不计算HReduceCount，计算ScaleFactor
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="currentSmallBlockPathList"></param>
        /// <param name="scaleFactor"></param>
        /// <returns></returns>
        private DataTree<Cell> DoDensityReduce_Scale(DataTree<Cell> cellDT, GH_Path currentSmallBlockPathList, double scaleFactor, double lMin, double w)
        {
            List<GH_Path> currentSmallBlockContainsCellPath = new List<GH_Path>();
            foreach (var cellPath in cellDT.Paths)
            {
                // 找到所有属于当前这个小Block的所有Cell
                if (cellPath[0] == currentSmallBlockPathList[0] && cellPath[1] == currentSmallBlockPathList[1])
                {
                    currentSmallBlockContainsCellPath.Add(cellPath);
                }
            }

            // 注意此时只会有一个Cell，所以这个For循环也只会循环一次
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                GenerateDensityConstrainedCell(cellDT, currentSmallBlockContainsCellPath[i], scaleFactor, lMin, w);
            }

            return cellDT;
        }

        private DataTree<Cell> DoDensityReduce_SP(DataTree<Cell> cellDT,GH_Path currentSmallBlockPathList, int vReduceCount, int hReduceCount, double sV, double sH, int indexToCut,bool isFromZero)
        {
            #region 得到当前小block的子树
            List<GH_Path> currentSmallBlockContainsCellPath = new List<GH_Path>();
            foreach (var cellPath in cellDT.Paths)
            {
                // 找到所有属于当前这个小Block的所有Cell
                if (cellPath[0] == currentSmallBlockPathList[0] && cellPath[1] == currentSmallBlockPathList[1])
                {
                    currentSmallBlockContainsCellPath.Add(cellPath);
                }
            }

            DataTree<Cell> subTree = new DataTree<Cell>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                subTree.EnsurePath(currentSmallBlockContainsCellPath[i]);
                subTree.Branch(currentSmallBlockContainsCellPath[i]).AddRange(cellDT.Branch(currentSmallBlockContainsCellPath[i]));
            }
            #endregion

            List<Tuple<GH_Path, int>> indexNotChangeListForSinglePolyline_vReduce = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForSinglePolyline_hReduce = new List<Tuple<GH_Path, int>>();

            DataTree<int> singlePolylineIndexDT = new DataTree<int>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                singlePolylineIndexDT.EnsurePath(currentSmallBlockContainsCellPath[i]);
                for (int j = 0; j < subTree.Branch(currentSmallBlockContainsCellPath[i]).Count; j++)
                {
                    BilinearCell cell = subTree.Branch(currentSmallBlockContainsCellPath[i])[j] as BilinearCell;
                    if (cell.ShapeType == BilinearCell.BilinearShapeType.SinglePolyline)
                    {
                        singlePolylineIndexDT.Branch(currentSmallBlockContainsCellPath[i]).Add(j);

                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(currentSmallBlockContainsCellPath[i], j);
                        indexNotChangeListForSinglePolyline_vReduce.Add(path_index);
                        indexNotChangeListForSinglePolyline_hReduce.Add(path_index);
                    }
                }
                
            }

            BilinearCell bCell = cellDT.Branch(singlePolylineIndexDT.Paths[0])[0] as BilinearCell;
            if (hReduceCount == 0 && vReduceCount == 0)
            {
                if (sH == 0)
                {
                    // 将原来的v变短
                    if (indexToCut == 1)
                    {
                        // 减E
                        if (isFromZero)
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 3)
                    {
                        // 减W
                        if (isFromZero)
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                }
                else if (sV == 0)
                {
                    // 将原来的h变短
                    if (indexToCut == 0)
                    {
                        // 减S
                        if (isFromZero)
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 2)
                    {
                        // 减N
                        if (isFromZero)
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }

                }
            }
            else if (hReduceCount == 1 && vReduceCount == 0)
            {
                if (bCell.South_Interval.Count != 0)
                {
                    bCell.South_Interval.Clear();

                    // 将原来的v变短
                    if (indexToCut == 1)
                    {
                        // 减E
                        if (isFromZero)
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 3)
                    {
                        // 减W
                        if (isFromZero)
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                }
                else if (bCell.North_Interval.Count != 0)
                {
                    bCell.North_Interval.Clear();

                    // 将原来的v变短
                    if (indexToCut == 1)
                    {
                        // 减E
                        if (isFromZero)
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.East_Interval.Clear();
                            bCell.East_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 3)
                    {
                        // 减W
                        if (isFromZero)
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.West_Interval.Clear();
                            bCell.West_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                }
            }
            else if (vReduceCount == 1 && hReduceCount == 0)
            {
                if (bCell.West_Interval.Count != 0)
                {
                    bCell.West_Interval.Clear();

                    // 将原来的h变短
                    if (indexToCut == 0)
                    {
                        // 减S
                        if (isFromZero)
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 2)
                    {
                        // 减N
                        if (isFromZero)
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                }
                else if (bCell.East_Interval.Count != 0)
                {
                    bCell.East_Interval.Clear();

                    // 将原来的h变短
                    if (indexToCut == 0)
                    {
                        // 减S
                        if (isFromZero)
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.South_Interval.Clear();
                            bCell.South_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                    else if (indexToCut == 2)
                    {
                        // 减N
                        if (isFromZero)
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(sV, 1));
                        }
                        else
                        {
                            bCell.North_Interval.Clear();
                            bCell.North_Interval.Add(new Interval(0, 1 - sV));
                        }
                    }
                }
            }

            cellDT.Branch(singlePolylineIndexDT.Paths[0])[0] = bCell;

            return cellDT;
        }

        /// <summary>
        /// 全部是 singlePolyline 或 IShape 或 CShape 时
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="currentSmallBlockPathList"></param>
        /// <param name="vReduceCount"></param>
        /// <param name="scaleFactor"></param>
        /// <returns></returns>
        private DataTree<Cell> DoDensityReduce_IOrC(DataTree<Cell> cellDT, GH_Path currentSmallBlockPathList,
                                                    int vReduceCount, int hReduceCount, 
                                                    int IShapeCount, int CShapeCount)
        {
            #region 得到当前小block的子树
            List<GH_Path> currentSmallBlockContainsCellPath = new List<GH_Path>();
            foreach (var cellPath in cellDT.Paths)
            {
                // 找到所有属于当前这个小Block的所有Cell
                if (cellPath[0] == currentSmallBlockPathList[0] && cellPath[1] == currentSmallBlockPathList[1])
                {
                    currentSmallBlockContainsCellPath.Add(cellPath);
                }
            }

            DataTree<Cell> subTree = new DataTree<Cell>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                subTree.EnsurePath(currentSmallBlockContainsCellPath[i]);
                subTree.Branch(currentSmallBlockContainsCellPath[i]).AddRange(cellDT.Branch(currentSmallBlockContainsCellPath[i]));
            }
            #endregion

            //List<Tuple<GH_Path, int>> indexForReduce = new List<Tuple<GH_Path, int>>();
            //List<Tuple<GH_Path, int>> indexForScale = new List<Tuple<GH_Path, int>>();

            List<Tuple<GH_Path, int>> indexNotChangeListForIShape_vReduceOnce = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForCShape_vReduceOnce = new List<Tuple<GH_Path, int>>();

            List<Tuple<GH_Path, int>> indexNotChangeListForIShape_hReduce = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForCShape_hReduce = new List<Tuple<GH_Path, int>>();

            DataTree<int> IShapeIndexDT = new DataTree<int>();
            DataTree<int> CShapeIndexDT = new DataTree<int>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                IShapeIndexDT.EnsurePath(currentSmallBlockContainsCellPath[i]);
                CShapeIndexDT.EnsurePath(currentSmallBlockContainsCellPath[i]);

                for (int j = 0; j < subTree.Branch(currentSmallBlockContainsCellPath[i]).Count; j++)
                {
                    BilinearCell bCell = subTree.Branch(currentSmallBlockContainsCellPath[i])[j] as BilinearCell;
                    if (bCell.ShapeType == BilinearCell.BilinearShapeType.IShape)
                    {
                        IShapeIndexDT.Branch(currentSmallBlockContainsCellPath[i]).Add(j);

                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(currentSmallBlockContainsCellPath[i], j);
                        indexNotChangeListForIShape_vReduceOnce.Add(path_index);
                        indexNotChangeListForIShape_hReduce.Add(path_index);
                    }
                    else if (bCell.ShapeType == BilinearCell.BilinearShapeType.CShape)
                    {
                        CShapeIndexDT.Branch(currentSmallBlockContainsCellPath[i]).Add(j);

                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(currentSmallBlockContainsCellPath[i], j);
                        indexNotChangeListForCShape_vReduceOnce.Add(path_index);
                        indexNotChangeListForCShape_hReduce.Add(path_index);
                    }
                }
            }

            #region 进行vReduce
            #region 求解四元一次不定方程，得到finalCountPairs
            // 求解四元一次不定方程问题：由于参数的特殊性，可以选择扩展二元不定方程的解法
            // x:IShape减2条    y:CShape减1条    z:CShape减2条    w:IShape减1条
            // countPair:[a,b,c,d]
            List<int[]> countPairs_vReduceOnce = new List<int[]>();
            int cd = 0;
            while (vReduceCount - cd >= 0)
            {
                List<int[]> abPairs_vReduce = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, vReduceCount - cd);
                for (int i = 0; i < abPairs_vReduce.Count; i++)
                {
                    List<int[]> countPairs_cd_vReduce = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, cd);
                    for (int j = 0; j < countPairs_cd_vReduce.Count; j++)
                    {
                        int[] abcdPair_vReduce = new int[4] { abPairs_vReduce[i][0], abPairs_vReduce[i][1], countPairs_cd_vReduce[j][0], countPairs_cd_vReduce[j][1] };
                        countPairs_vReduceOnce.Add(abcdPair_vReduce);
                    }
                }
                cd++;
            }

            List<int[]> distinctPairs_vReduceOnce = countPairs_vReduceOnce.Distinct().ToList();

            List<int[]> finalCountPairs_vReduceOnce = new List<int[]>();
            for (int i = 0; i < distinctPairs_vReduceOnce.Count; i++)
            {
                if (distinctPairs_vReduceOnce[i][0] + distinctPairs_vReduceOnce[i][3] <= IShapeCount
                    && distinctPairs_vReduceOnce[i][1] + distinctPairs_vReduceOnce[i][2] <= CShapeCount)
                {
                    bool flag0 = false;
                    bool flag1 = false;
                    bool flag2 = false;
                    bool flag3 = false;
                    if (distinctPairs_vReduceOnce[i][0] >= 0)
                    {
                        flag0 = true;
                    }
                    if (distinctPairs_vReduceOnce[i][1] >= 0)
                    {
                        flag1 = true;
                    }
                    if (distinctPairs_vReduceOnce[i][2] >= 0)
                    {
                        flag2 = true;
                    }
                    if (distinctPairs_vReduceOnce[i][3] >= 0)
                    {
                        flag3 = true;
                    }

                    if (flag0 && flag1 && flag2 && flag3)
                    {
                        finalCountPairs_vReduceOnce.Add(distinctPairs_vReduceOnce[i]);
                    }
                }
            }
            #endregion

            if (finalCountPairs_vReduceOnce.Count != 0)
            {
                int randomCountPairsIndex = m_random.Next(0, finalCountPairs_vReduceOnce.Count);

                int[] finalPair = finalCountPairs_vReduceOnce[randomCountPairsIndex];

                List<Tuple<GH_Path, int>> indexToChangeListForIShapeCutTwo_vReduceOnce = new List<Tuple<GH_Path, int>>();
                List<Tuple<GH_Path, int>> indexToChangeListForCShapeCutOne_vReduceOnce = new List<Tuple<GH_Path, int>>();
                List<Tuple<GH_Path, int>> indexToChangeListForCShapeCutTwo_vReduceOnce = new List<Tuple<GH_Path, int>>();
                List<Tuple<GH_Path, int>> indexToChangeListForIShapeCutOne_vReduceOnce = new List<Tuple<GH_Path, int>>();

                // 对于IShapeCutTwo来说
                int temp = 0;
                int length = finalPair[0];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, IShapeIndexDT.Paths.Count);
                    GH_Path randomPath = IShapeIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, IShapeIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, IShapeIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForIShapeCutTwo_vReduceOnce.Contains(path_index)
                        && !indexToChangeListForIShapeCutOne_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForIShapeCutTwo_vReduceOnce.Add(path_index);
                        indexNotChangeListForIShape_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }
                // 对于CShapeCutOne来说
                temp = 0;
                length = finalPair[1];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, CShapeIndexDT.Paths.Count);
                    GH_Path randomPath = CShapeIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, CShapeIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, CShapeIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForCShapeCutOne_vReduceOnce.Contains(path_index)
                        && !indexToChangeListForCShapeCutTwo_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForCShapeCutOne_vReduceOnce.Add(path_index);
                        indexNotChangeListForCShape_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }
                // 对于CShapeCutTwo来说
                temp = 0;
                length = finalPair[2];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, CShapeIndexDT.Paths.Count);
                    GH_Path randomPath = CShapeIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, CShapeIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, CShapeIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForCShapeCutTwo_vReduceOnce.Contains(path_index)
                        && !indexToChangeListForCShapeCutOne_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForCShapeCutTwo_vReduceOnce.Add(path_index);
                        indexNotChangeListForCShape_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }
                // 对于IShapeCutOne来说
                temp = 0;
                length = finalPair[3];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, IShapeIndexDT.Paths.Count);
                    GH_Path randomPath = IShapeIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, IShapeIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, IShapeIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForIShapeCutOne_vReduceOnce.Contains(path_index)
                        && !indexToChangeListForIShapeCutTwo_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForIShapeCutOne_vReduceOnce.Add(path_index);
                        indexNotChangeListForIShape_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }

                GenerateDensityConstrainedCell_vReduceOnce(cellDT,
                                                           indexToChangeListForIShapeCutTwo_vReduceOnce,
                                                           indexToChangeListForCShapeCutOne_vReduceOnce,
                                                           indexToChangeListForCShapeCutTwo_vReduceOnce,
                                                           indexToChangeListForIShapeCutOne_vReduceOnce,
                                                           indexNotChangeListForIShape_vReduceOnce, indexNotChangeListForCShape_vReduceOnce);
            }
            else
            {
                // 应该不会出现，待测试
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处出现了无法分配的vReduceCountOnce：{1}，此时的IShapeCount：{2}，CShapeCount：{3}", currentSmallBlockPathList.ToString(), vReduceCount, IShapeCount, CShapeCount));
                return null;
            }


            #endregion

            if (hReduceCount != 0)
            {
                #region 进行hReduce
                #region 求解二元一次不定方程，得到finalCountPairs
                List<int[]> countPairs_hReduce = UtilityFunctions.GetAllPositiveIntegerSolution(1, 1, hReduceCount);

                List<int[]> distinctPairs_hReduce = countPairs_hReduce.Distinct().ToList();

                List<int[]> finalCountPairs_hReduce = new List<int[]>();
                for (int i = 0; i < distinctPairs_hReduce.Count; i++)
                {
                    if (distinctPairs_hReduce[i][0] <= IShapeCount
                        && distinctPairs_hReduce[i][1] <= CShapeCount)
                    {
                        bool flag0 = false;
                        bool flag1 = false;
                        if (distinctPairs_hReduce[i][0] >= 0)
                        {
                            flag0 = true;
                        }
                        if (distinctPairs_hReduce[i][1] >= 0)
                        {
                            flag1 = true;
                        }

                        if (flag0 && flag1)
                        {
                            finalCountPairs_hReduce.Add(distinctPairs_hReduce[i]);
                        }
                    }
                }
                #endregion

                if (finalCountPairs_hReduce.Count != 0)
                {
                    int randomCountPairsIndex = m_random.Next(0, finalCountPairs_hReduce.Count);
                    int[] finalPair = finalCountPairs_hReduce[randomCountPairsIndex];

                    List<Tuple<GH_Path, int>> indexToChangeListForIShape_hReduce = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForCShape_hReduce = new List<Tuple<GH_Path, int>>();

                    // 对于IShape来说
                    int temp = 0;
                    int length = finalPair[0];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, IShapeIndexDT.Paths.Count);
                        GH_Path randomPath = IShapeIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, IShapeIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, IShapeIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForIShape_hReduce.Contains(path_index))
                        {
                            indexToChangeListForIShape_hReduce.Add(path_index);
                            indexNotChangeListForIShape_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于CShape来说
                    temp = 0;
                    length = finalPair[1];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, CShapeIndexDT.Paths.Count);
                        GH_Path randomPath = CShapeIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, CShapeIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, CShapeIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForCShape_hReduce.Contains(path_index))
                        {
                            indexToChangeListForCShape_hReduce.Add(path_index);
                            indexNotChangeListForCShape_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                }
                else
                {
                    // 应该不会出现，待测试
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处出现了无法分配的hReduceCount：{1}，此时的IShapeCount：{2}，CShapeCount：{3}", currentSmallBlockPathList.ToString(), hReduceCount, IShapeCount, CShapeCount));
                    return null;
                }
                #endregion
            }

            return cellDT;
        }

        /// <summary>
        /// 全部为RecShape和EightShape时
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="currentSmallBlockPathList"></param>
        /// <param name="vReduceCountOnce"></param>
        /// <param name="vReduceCountTwice"></param>
        /// <param name="hReduceCount"></param>
        /// <param name="RecShapeCount"></param>
        /// <param name="EightShapeCount"></param>
        /// <param name="lMin"></param>
        /// <param name="w"></param>
        /// <returns></returns>
        private DataTree<Cell> DoDensityReduce_RecOrEight(DataTree<Cell> cellDT, GH_Path currentSmallBlockPathList, 
                                                          int vReduceCountOnce, int vReduceCountTwice, int hReduceCount, 
                                                          int RecShapeCount, int EightShapeCount,
                                                          double lMin, double w)
        {
            #region 得到当前小block的子树
            List<GH_Path> currentSmallBlockContainsCellPath = new List<GH_Path>();
            foreach (var cellPath in cellDT.Paths)
            {
                // 找到所有属于当前这个小Block的所有Cell
                if (cellPath[0] == currentSmallBlockPathList[0] && cellPath[1] == currentSmallBlockPathList[1])
                {
                    currentSmallBlockContainsCellPath.Add(cellPath);
                }
            }

            DataTree<Cell> subTree = new DataTree<Cell>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                subTree.EnsurePath(currentSmallBlockContainsCellPath[i]);
                subTree.Branch(currentSmallBlockContainsCellPath[i]).AddRange(cellDT.Branch(currentSmallBlockContainsCellPath[i]));
            }
            #endregion

            List<Tuple<GH_Path, int>> indexNotChangeListForBiCell_vReduceOnce = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForTriCell_vReduceOnce = new List<Tuple<GH_Path, int>>();

            List<Tuple<GH_Path, int>> indexNotChangeListForBiCell_vReduceTwice = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForTriCell_vReduceTwice = new List<Tuple<GH_Path, int>>();

            List<Tuple<GH_Path, int>> indexNotChangeListForBiCell_hReduce = new List<Tuple<GH_Path, int>>();
            List<Tuple<GH_Path, int>> indexNotChangeListForTriCell_hReduce = new List<Tuple<GH_Path, int>>();

            DataTree<int> biCellIndexDT = new DataTree<int>();
            DataTree<int> triCellIndexDT = new DataTree<int>();
            for (int i = 0; i < currentSmallBlockContainsCellPath.Count; i++)
            {
                biCellIndexDT.EnsurePath(currentSmallBlockContainsCellPath[i]);
                triCellIndexDT.EnsurePath(currentSmallBlockContainsCellPath[i]);

                for (int j = 0; j < subTree.Branch(currentSmallBlockContainsCellPath[i]).Count; j++)
                {
                    if (subTree.Branch(currentSmallBlockContainsCellPath[i])[j] is BilinearCell)
                    {
                        biCellIndexDT.Branch(currentSmallBlockContainsCellPath[i]).Add(j);
                        //biCellCount++;

                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(currentSmallBlockContainsCellPath[i], j);
                        indexNotChangeListForBiCell_vReduceOnce.Add(path_index);
                        indexNotChangeListForBiCell_vReduceTwice.Add(path_index);
                        indexNotChangeListForBiCell_hReduce.Add(path_index);
                    }
                    else
                    {
                        triCellIndexDT.Branch(currentSmallBlockContainsCellPath[i]).Add(j);
                        //triCellCount++;

                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(currentSmallBlockContainsCellPath[i], j);
                        indexNotChangeListForTriCell_vReduceOnce.Add(path_index);
                        indexNotChangeListForTriCell_vReduceTwice.Add(path_index);
                        indexNotChangeListForTriCell_hReduce.Add(path_index);
                    }
                }
            }

            #region 进行第一轮vReduce
            #region 求解三元一次不定方程，得到finalCountPairs
            // 求解三元一次不定方程问题：由于参数的特殊性，可以选择扩展二元不定方程的解法
            // x:TriCell减2条    y:BiCell减1条    z:TriCell减1条
            // countPair:[a,b,c]
            List<int[]> countPairs_vReduceOnce = new List<int[]>();
            int v = 0;
            while (vReduceCountOnce - v >= 0)
            {
                List<int[]> abPairs = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, vReduceCountOnce - v);
                //List<int[]> abcPairs = new List<int[]>();
                for (int i = 0; i < abPairs.Count; i++)
                {
                    int[] abcPair = new int[3] { abPairs[i][0], abPairs[i][1], v };
                    countPairs_vReduceOnce.Add(abcPair);
                }
                v++;
            }

            List<int[]> distinctPairs_vReduceOnce = countPairs_vReduceOnce.Distinct().ToList();

            //List<int[]> countPairs = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, vReduceCount);
            List<int[]> finalCountPairs_vReduceOnce = new List<int[]>();
            for (int i = 0; i < distinctPairs_vReduceOnce.Count; i++)
            {
                if (distinctPairs_vReduceOnce[i][0]+ distinctPairs_vReduceOnce[i][2] <= EightShapeCount
                    && distinctPairs_vReduceOnce[i][1] <= RecShapeCount)
                {
                    bool flag0 = false;
                    bool flag1 = false;
                    bool flag2 = false;
                    if (distinctPairs_vReduceOnce[i][0] >= 0)
                    {
                        flag0 = true;
                    }
                    if (distinctPairs_vReduceOnce[i][1] >= 0)
                    {
                        flag1 = true;
                    }
                    if (distinctPairs_vReduceOnce[i][2] >= 0)
                    {
                        flag2 = true;
                    }

                    if (flag0 && flag1 && flag2)
                    {
                        finalCountPairs_vReduceOnce.Add(distinctPairs_vReduceOnce[i]);
                    }
                }
            }
            #endregion

            // 随机选择一组finalCountPairs 对Basic Cell进行削减生成
            //List<Cell> densityConstrainedCells = new List<Cell>();
            if (finalCountPairs_vReduceOnce.Count != 0)
            {
                int randomCountPairsIndex = m_random.Next(0, finalCountPairs_vReduceOnce.Count);

                int[] finalPair = finalCountPairs_vReduceOnce[randomCountPairsIndex];

                List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo_vReduceOnce = new List<Tuple<GH_Path, int>>();
                List<Tuple<GH_Path, int>> indexToChangeListForBiCell_vReduceOnce = new List<Tuple<GH_Path, int>>();
                List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne_vReduceOnce = new List<Tuple<GH_Path, int>>();

                // 对于triCellCutTwo来说
                int temp = 0;
                int length = finalPair[0];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                    GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForTriCellCutTwo_vReduceOnce.Contains(path_index)
                        && !indexToChangeListForTriCellCutOne_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForTriCellCutTwo_vReduceOnce.Add(path_index);
                        indexNotChangeListForTriCell_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }
                // 对于biCell来说
                temp = 0;
                length = finalPair[1];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, biCellIndexDT.Paths.Count);
                    GH_Path randomPath = biCellIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, biCellIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, biCellIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForBiCell_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForBiCell_vReduceOnce.Add(path_index);
                        indexNotChangeListForBiCell_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }
                // 对于TriCellCutOne来说
                temp = 0;
                length = finalPair[2];
                while (temp < length)
                {
                    int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                    GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                    int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                    Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                    if (!indexToChangeListForTriCellCutOne_vReduceOnce.Contains(path_index) 
                        && !indexToChangeListForTriCellCutTwo_vReduceOnce.Contains(path_index))
                    {
                        indexToChangeListForTriCellCutOne_vReduceOnce.Add(path_index);
                        indexNotChangeListForTriCell_vReduceOnce.Remove(path_index);
                        temp++;
                    }
                }

                 GenerateDensityConstrainedCell_vReduceOnce(cellDT,
                                               indexToChangeListForBiCell_vReduceOnce,
                                               indexToChangeListForTriCellCutTwo_vReduceOnce,
                                               indexToChangeListForTriCellCutOne_vReduceOnce,
                                               indexNotChangeListForBiCell_vReduceOnce,
                                               indexNotChangeListForTriCell_vReduceOnce,
                                               lMin, w);
            }
            else
            {
                // 应该不会出现，待测试
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处出现了无法分配的vReduceCountOnce：{1}，此时的TriCellCount：{2}，BiCellCount：{3}", currentSmallBlockPathList.ToString(), vReduceCountOnce, EightShapeCount, RecShapeCount));
                return null;
            }

            #endregion

            if (vReduceCountTwice != 0)
            {
                #region 进行第二轮vReduce
                #region 求解三元一次不定方程，得到finalCountPairs
                // 求解三元一次不定方程问题：由于参数的特殊性，可以选择扩展二元不定方程的解法
                // x:TriCell减2条    y:BiCell减1条    z:TriCell减1条
                // countPair:[a,b,c]
                List<int[]> countPairs_vReduceTwice = new List<int[]>();
                v = 0;
                while (vReduceCountTwice - v >= 0)
                {
                    List<int[]> abPairs_vReduceTwice = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, vReduceCountTwice - v);
                    for (int i = 0; i < abPairs_vReduceTwice.Count; i++)
                    {
                        int[] abcPair_vReduceTwice = new int[3] { abPairs_vReduceTwice[i][0], abPairs_vReduceTwice[i][1], v };
                        countPairs_vReduceTwice.Add(abcPair_vReduceTwice);
                    }
                    v++;
                }

                List<int[]> distinctPairs_vReduceTwice = countPairs_vReduceTwice.Distinct().ToList();

                List<int[]> finalCountPairs_vReduceTwice = new List<int[]>();
                for (int i = 0; i < distinctPairs_vReduceTwice.Count; i++)
                {
                    if (distinctPairs_vReduceTwice[i][0] + distinctPairs_vReduceTwice[i][2] <= EightShapeCount
                        && distinctPairs_vReduceTwice[i][1] <= RecShapeCount)
                    {
                        bool flag0 = false;
                        bool flag1 = false;
                        bool flag2 = false;
                        if (distinctPairs_vReduceTwice[i][0] >= 0)
                        {
                            flag0 = true;
                        }
                        if (distinctPairs_vReduceTwice[i][1] >= 0)
                        {
                            flag1 = true;
                        }
                        if (distinctPairs_vReduceTwice[i][2] >= 0)
                        {
                            flag2 = true;
                        }

                        if (flag0 && flag1 && flag2)
                        {
                            finalCountPairs_vReduceTwice.Add(distinctPairs_vReduceTwice[i]);
                        }
                    }
                }
                #endregion

                // 随机选择一组finalCountPairs_vReduceTwice 对Basic Cell进行削减生成
                if (finalCountPairs_vReduceTwice.Count != 0)
                {
                    int randomCountPairsIndex = m_random.Next(0, finalCountPairs_vReduceTwice.Count);

                    int[] finalPair = finalCountPairs_vReduceTwice[randomCountPairsIndex];

                    List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo_vReduceTwice = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForBiCell_vReduceTwice = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne_vReduceTwice = new List<Tuple<GH_Path, int>>();

                    // 对于triCellCutTwo来说
                    int temp = 0;
                    int length = finalPair[0];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                        GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForTriCellCutTwo_vReduceTwice.Contains(path_index))
                        {
                            indexToChangeListForTriCellCutTwo_vReduceTwice.Add(path_index);
                            indexNotChangeListForTriCell_vReduceTwice.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于biCell来说
                    temp = 0;
                    length = finalPair[1];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, biCellIndexDT.Paths.Count);
                        GH_Path randomPath = biCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, biCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, biCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForBiCell_vReduceTwice.Contains(path_index))
                        {
                            indexToChangeListForBiCell_vReduceTwice.Add(path_index);
                            indexNotChangeListForBiCell_vReduceTwice.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于TriCellCutOne来说
                    temp = 0;
                    length = finalPair[2];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                        GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForTriCellCutOne_vReduceTwice.Contains(path_index) && !indexToChangeListForTriCellCutTwo_vReduceTwice.Contains(path_index))
                        {
                            indexToChangeListForTriCellCutOne_vReduceTwice.Add(path_index);
                            indexNotChangeListForTriCell_vReduceTwice.Remove(path_index);
                            temp++;
                        }
                    }

                    GenerateDensityConstrainedCell_vReduceTwice(cellDT,
                                                                indexToChangeListForBiCell_vReduceTwice,
                                                                indexToChangeListForTriCellCutTwo_vReduceTwice,
                                                                indexToChangeListForTriCellCutOne_vReduceTwice,
                                                                indexNotChangeListForBiCell_vReduceTwice, indexNotChangeListForTriCell_vReduceTwice);
                }
                else
                {
                    // 应该不会出现，待测试
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处出现了无法分配的vReduceCountTwice：{1}，此时的TriCellCount：{2}，BiCellCount：{3}", currentSmallBlockPathList.ToString(), vReduceCountTwice, EightShapeCount, RecShapeCount));
                    return null;
                }
                #endregion
            }

            if (hReduceCount != 0)
            {
                #region 进行hReduce
                #region 求解五元一次不定方程，得到finalCountPairs
                // 求解五元一次不定方程问题：由于参数的特殊性，可以选择扩展二元不定方程的解法
                // x:TriCell减3条    y:BiCell减2条    z:TriCell减2条    w:BiCell减1条    v:TriCell减1条
                // countPair:[a,b,c,d,e]
                List<int[]> countPairs_hReduce = new List<int[]>();
                int cde = 0;
                while (hReduceCount - cde >= 0)
                {
                    List<int[]> abPairs_hReduce = UtilityFunctions.GetAllPositiveIntegerSolution(3, 2, hReduceCount - cde);
                    for (int i = 0; i < abPairs_hReduce.Count; i++)
                    {
                        List<int[]> countPairs_cde_hReduce = new List<int[]>();
                        int e = 0;
                        while (cde - e >= 0)
                        {
                            List<int[]> cdPairs_hReduce = UtilityFunctions.GetAllPositiveIntegerSolution(2, 1, cde - e);
                            for (int j = 0; j < cdPairs_hReduce.Count; j++)
                            {
                                int[] cdePair_hReduce = new int[3] { cdPairs_hReduce[j][0], cdPairs_hReduce[j][1], e };
                                countPairs_cde_hReduce.Add(cdePair_hReduce);
                            }
                            e++;
                        }

                        for (int j = 0; j < countPairs_cde_hReduce.Count; j++)
                        {
                            int[] abcdePair_hReduce = new int[5] { abPairs_hReduce[i][0], abPairs_hReduce[i][1], countPairs_cde_hReduce[j][0], countPairs_cde_hReduce[j][1], countPairs_cde_hReduce[j][2] };
                            countPairs_hReduce.Add(abcdePair_hReduce);
                        }
                    }
                    cde++;
                }

                List<int[]> distinctPairs_hReduce = countPairs_hReduce.Distinct().ToList();

                List<int[]> finalCountPairs_hReduce = new List<int[]>();
                for (int i = 0; i < distinctPairs_hReduce.Count; i++)
                {
                    if (distinctPairs_hReduce[i][0] + distinctPairs_hReduce[i][2] + distinctPairs_hReduce[i][4] <= EightShapeCount
                        && distinctPairs_hReduce[i][1] + distinctPairs_hReduce[i][3] <= RecShapeCount)
                    {
                        bool flag0 = false;
                        bool flag1 = false;
                        bool flag2 = false;
                        bool flag3 = false;
                        bool flag4 = false;
                        if (distinctPairs_hReduce[i][0] >= 0)
                        {
                            flag0 = true;
                        }
                        if (distinctPairs_hReduce[i][1] >= 0)
                        {
                            flag1 = true;
                        }
                        if (distinctPairs_hReduce[i][2] >= 0)
                        {
                            flag2 = true;
                        }
                        if (distinctPairs_hReduce[i][3] >= 0)
                        {
                            flag3 = true;
                        }
                        if (distinctPairs_hReduce[i][4] >= 0)
                        {
                            flag4 = true;
                        }

                        if (flag0 && flag1 && flag2 && flag3 && flag4)
                        {
                            finalCountPairs_hReduce.Add(distinctPairs_hReduce[i]);
                        }
                    }
                }

                #endregion

                if (finalCountPairs_hReduce.Count != 0)
                {
                    int randomCountPairsIndex = m_random.Next(0, finalCountPairs_hReduce.Count);

                    int[] finalPair = finalCountPairs_hReduce[randomCountPairsIndex];

                    List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutThree_hReduce = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForBiCellCutTwo_hReduce = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo_hReduce = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForBiCellCutOne_hReduce = new List<Tuple<GH_Path, int>>();
                    List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne_hReduce = new List<Tuple<GH_Path, int>>();

                    // 对于triCellCutThree来说
                    int temp = 0;
                    int length = finalPair[0];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                        GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForTriCellCutThree_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutTwo_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutOne_hReduce.Contains(path_index))
                        {
                            indexToChangeListForTriCellCutThree_hReduce.Add(path_index);
                            indexNotChangeListForTriCell_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于biCellCutTwo来说
                    temp = 0;
                    length = finalPair[1];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, biCellIndexDT.Paths.Count);
                        GH_Path randomPath = biCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, biCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, biCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForBiCellCutTwo_hReduce.Contains(path_index)
                            && !indexToChangeListForBiCellCutOne_hReduce.Contains(path_index))
                        {
                            indexToChangeListForBiCellCutTwo_hReduce.Add(path_index);
                            indexNotChangeListForBiCell_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于triCellCutTwo来说
                    temp = 0;
                    length = finalPair[2];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                        GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForTriCellCutTwo_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutThree_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutOne_hReduce.Contains(path_index))
                        {
                            indexToChangeListForTriCellCutTwo_hReduce.Add(path_index);
                            indexNotChangeListForTriCell_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于biCellCutOne来说
                    temp = 0;
                    length = finalPair[3];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, biCellIndexDT.Paths.Count);
                        GH_Path randomPath = biCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, biCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, biCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForBiCellCutOne_hReduce.Contains(path_index)
                            && !indexToChangeListForBiCellCutTwo_hReduce.Contains(path_index))
                        {
                            indexToChangeListForBiCellCutOne_hReduce.Add(path_index);
                            indexNotChangeListForBiCell_hReduce.Remove(path_index);
                            temp++;
                        }
                    }
                    // 对于triCellCutOne来说
                    temp = 0;
                    length = finalPair[4];
                    while (temp < length)
                    {
                        int randomPathIndex = m_random.Next(0, triCellIndexDT.Paths.Count);
                        GH_Path randomPath = triCellIndexDT.Paths[randomPathIndex];
                        int randomIndex = m_random.Next(0, triCellIndexDT.Branch(randomPath).Count);
                        Tuple<GH_Path, int> path_index = new Tuple<GH_Path, int>(randomPath, triCellIndexDT.Branch(randomPath)[randomIndex]);
                        if (!indexToChangeListForTriCellCutOne_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutThree_hReduce.Contains(path_index)
                            && !indexToChangeListForTriCellCutTwo_hReduce.Contains(path_index))
                        {
                            indexToChangeListForTriCellCutOne_hReduce.Add(path_index);
                            indexNotChangeListForTriCell_hReduce.Remove(path_index);
                            temp++;
                        }
                    }

                    GenerateDensityConstrainedCell_hReduce(cellDT,
                                                           indexToChangeListForBiCellCutTwo_hReduce,
                                                           indexToChangeListForBiCellCutOne_hReduce,
                                                           indexToChangeListForTriCellCutThree_hReduce,
                                                           indexToChangeListForTriCellCutTwo_hReduce,
                                                           indexToChangeListForTriCellCutOne_hReduce,
                                                           indexNotChangeListForBiCell_hReduce, indexNotChangeListForTriCell_hReduce);
                }
                else
                {
                    // 应该不会出现，待测试
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, string.Format("Path:{0} 处出现了无法分配的hReduceCount：{1}，此时的TriCellCount：{2}，BiCellCount：{3}", currentSmallBlockPathList.ToString(), hReduceCount, EightShapeCount, RecShapeCount));
                    return null;
                }
                #endregion
            }

            return cellDT;
        }



        /// <summary>
        /// 此时不计算VReduceCount，不计算HReduceCount，计算ScaleFactor
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="path"></param>
        /// <param name="scaleFactor"></param>
        private void GenerateDensityConstrainedCell(DataTree<Cell> cellDT, GH_Path path, double scaleFactor, double lMin, double w)
        {
            //BilinearCell bCell = cellDT.Branch(path)[0] as BilinearCell;
            //cellDT.Branch(path)[0] = bCell.GenerateOffsetedSingleRegion(scaleFactor);

            for (int i = 0; i < cellDT.Branch(path).Count; i++)
            {
                if (cellDT.Branch(path)[i] is BilinearCell)
                {
                    BilinearCell bCell = cellDT.Branch(path)[i] as BilinearCell;
                    cellDT.Branch(path)[i] = bCell.GenerateScaledShape(scaleFactor, lMin, w);
                }
                else
                {
                    TrilinearCell tCell = cellDT.Branch(path)[i] as TrilinearCell;
                    cellDT.Branch(path)[i] = tCell.GenerateScaledShape(scaleFactor, lMin, w);
                }
            }
        }

        ///// <summary>
        ///// 全部是 singlePolyline 或 IShape 或 CShape 时
        ///// </summary>
        ///// <param name="cellDT"></param>
        ///// <param name="indexForReduce"></param>
        ///// <param name="indexForScale"></param>
        ///// <param name="scaleFactor"></param>
        //private void GenerateDensityConstrainedCell(DataTree<Cell> cellDT, List<Tuple<GH_Path, int>> indexForReduce, List<Tuple<GH_Path, int>> indexForScale,double scaleFactor, double lMin, double w)
        //{
        //    #region 对于需要进行reduce操作的
        //    foreach (var path_index in indexForReduce)
        //    {
        //        GH_Path path = path_index.Item1;
        //        int index = path_index.Item2;

        //        BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
        //        int 
        //        cellDT.Branch(path)[index] = bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline();
        //    }
        //    #endregion

        //    #region 对于需要进行scale操作的
        //    foreach (var path_index in indexForScale)
        //    {
        //        GH_Path path = path_index.Item1;
        //        int index = path_index.Item2;
        //        BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
        //        cellDT.Branch(path)[index] = bCell.GenerateScaledShape(scaleFactor, lMin, w);
        //    }
        //    #endregion
        //}


        /// <summary>
        /// 全部是 IShape 或 CShape 时
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="indexToChangeListForIShapeCutTwo_vReduceOnce"></param>
        /// <param name="indexToChangeListForCShapeCutOne_vReduceOnce"></param>
        /// <param name="indexToChangeListForCShapeCutTwo_vReduceOnce"></param>
        /// <param name="indexToChangeListForIShapeCutOne_vReduceOnce"></param>
        /// <param name="indexNotChangeListForIShape_vReduceOnce"></param>
        /// <param name="indexNotChangeListForCShape_vReduceOnce"></param>
        private void GenerateDensityConstrainedCell_vReduceOnce(DataTree<Cell> cellDT,
                                                                List<Tuple<GH_Path, int>> indexToChangeListForIShapeCutTwo_vReduceOnce,
                                                                List<Tuple<GH_Path, int>> indexToChangeListForCShapeCutOne_vReduceOnce,
                                                                List<Tuple<GH_Path, int>> indexToChangeListForCShapeCutTwo_vReduceOnce,
                                                                List<Tuple<GH_Path, int>> indexToChangeListForIShapeCutOne_vReduceOnce,
                                                                List<Tuple<GH_Path, int>> indexNotChangeListForIShape_vReduceOnce, List<Tuple<GH_Path, int>> indexNotChangeListForCShape_vReduceOnce)
        {
            #region 对于需要修改的
            foreach (var path_index in indexToChangeListForIShapeCutTwo_vReduceOnce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                cellDT.Branch(path)[index] = bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline(2, 0);
            }

            foreach (var path_index in indexToChangeListForCShapeCutOne_vReduceOnce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                int directionCode0 = m_random.Next(0, 4);
                cellDT.Branch(path)[index] = bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline(1, directionCode0);
            }

            foreach (var path_index in indexToChangeListForCShapeCutTwo_vReduceOnce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                cellDT.Branch(path)[index] = bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline(2, 0);
            }

            foreach (var path_index in indexToChangeListForIShapeCutOne_vReduceOnce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                int directionCode0 = m_random.Next(0, 4);
                cellDT.Branch(path)[index] =  bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline(1, directionCode0);
            }
            #endregion

            #region 其他不需要修改的
            // 不动
            #endregion
        }

        /// <summary>
        /// 全部都是 elseCount 时，此时计算VReduceCount，不计算OffsetValue
        /// </summary>
        /// <param name="cellDT"></param>
        /// <param name="indexToChangeListForBiCell"></param>
        /// <param name="indexToChangeListForTriCellCutTwo"></param>
        /// <param name="indexNotChangeListForBiCell"></param>
        /// <param name="indexNotChangeListForTriCell"></param>
        /// <param name="lMin"></param>
        /// <param name="w"></param>
        private void GenerateDensityConstrainedCell_vReduceOnce(DataTree<Cell> cellDT, 
                                                                List<Tuple<GH_Path, int>> indexToChangeListForBiCell, 
                                                                List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo,
                                                                List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne,
                                                                List<Tuple<GH_Path, int>> indexNotChangeListForBiCell, List<Tuple<GH_Path, int>> indexNotChangeListForTriCell,
                                                                double lMin, double w)
        {
            #region 对于需要修改的
            foreach (var path_index in indexToChangeListForBiCell)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                // 随机生成工,C
                int random0 = m_random.Next(0, 2);// 取值范围是0, 1
                double randomT0 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                int directionCode0 = m_random.Next(2, 4);
                // 如果在面宽方向上被切割过，就不走IShape分支
                if ((bCell.South_Interval != null && bCell.South_Interval.Count > 1)
                    || (bCell.North_Interval != null && bCell.North_Interval.Count > 1))
                {
                    random0 = m_random.Next(1, 2);
                }

                if (random0 == 0)
                {
                    // 此时只会生成 this.MainVolumeDirection = MainDirection.SouthNorth 的
                    cellDT.Branch(path)[index] = bCell.GenerateIAndIVariantShape(randomT0, directionCode0, w);
                }
                else
                {
                    // 此时只会生成 this.MainVolumeDirection = MainDirection.SouthNorth 的
                    cellDT.Branch(path)[index] = bCell.GenerateCAndCVariantShape(randomT0, directionCode0, w);
                }
            }

            foreach (var path_index in indexToChangeListForTriCellCutTwo)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;

                // 随机生成E,2,王,
                int random1 = m_random.Next(0, 3);// 取值范围是0,1,2
                double randomT0 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                double randomT1 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                int directionCode1 = m_random.Next(0, 4);
                // 如果在面宽方向上被切割过，就不走EVariantShape分支
                if ((tCell.South_Interval != null && tCell.South_Interval.Count > 1)
                    || (tCell.Middle_Interval != null && tCell.Middle_Interval.Count > 1)
                    || (tCell.North_Interval != null && tCell.North_Interval.Count > 1))
                {
                    random1 = m_random.Next(1, 3);
                }

                if (random1 == 0)
                {
                    cellDT.Branch(path)[index] = tCell.GenerateEVariantShape(directionCode1, randomT0, randomT1, w);
                }
                else if (random1 == 1)
                {
                    cellDT.Branch(path)[index] = tCell.GenerateEShape(directionCode1);
                }
                else
                {
                    cellDT.Branch(path)[index] = tCell.GenerateZShape(directionCode1);
                }
            }

            foreach (var path_index in indexToChangeListForTriCellCutOne)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                // 随机生成6字形及其变体
                double randomT0 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                int directionCode1 = m_random.Next(0, 6);// 取值范围是0,1,2,3,4,5

                cellDT.Branch(path)[index] = tCell.GenerateSixAndSixVariantShape(randomT0, directionCode1, lMin, w);
            }
            #endregion

            #region 其他不需要修改的Rec或Eight
            foreach (var path_index in indexNotChangeListForBiCell)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;
                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;

                // 随机生成Rec或Rec变体
                // 如果在面宽方向上被切割过，就不走RecVariantShape分支
                double randomT0 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                int directionCode0 = m_random.Next(0, 4);
                bool flag = false;
                if ((bCell.South_Interval != null && bCell.South_Interval.Count > 1)
                        || (bCell.North_Interval != null && bCell.North_Interval.Count > 1))
                {
                    flag = true;
                }

                cellDT.Branch(path)[index] = bCell.GenerateRecAndRecVariantShape(randomT0, directionCode0, lMin, w, flag);
            }

            foreach (var path_index in indexNotChangeListForTriCell)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;
                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;

                // 随机生成Eight或Eight变体
                // 如果在面宽方向上被切割过，就不走EightVariantShape分支
                double randomT1 = m_random.Next(0, 3) / 2; // 取值范围是0, 0.5, 1
                int directionCode1 = m_random.Next(0, 4);
                bool flag = false;
                if ((tCell.South_Interval != null && tCell.South_Interval.Count > 1)
                    || (tCell.Middle_Interval != null && tCell.Middle_Interval.Count > 1)
                    || (tCell.North_Interval != null && tCell.North_Interval.Count > 1))
                {
                    flag = true;
                }

                cellDT.Branch(path)[index] = tCell.GenerateEightAndEightVariantShape(randomT1, directionCode1, lMin, w, flag);
            }
            #endregion
        }

        private void GenerateDensityConstrainedCell_vReduceTwice(DataTree<Cell> cellDT,
                                                                 List<Tuple<GH_Path, int>> indexToChangeListForBiCell_vReduceTwice,
                                                                 List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo_vReduceTwice,
                                                                 List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne_vReduceTwice,
                                                                 List<Tuple<GH_Path, int>> indexNotChangeListForBiCell_vReduceTwice, List<Tuple<GH_Path, int>> indexNotChangeListForTriCell_vReduceTwice)
        {
            #region 对于需要修改的
            foreach (var path_index in indexToChangeListForBiCell_vReduceTwice)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                // 去掉最后一竖
                int directionCode0 = m_random.Next(0, 4);
                cellDT.Branch(path)[index] = bCell.GenerateRemovedCShapeOrIShapeOrSinglePolyline(1, directionCode0);
            }

            foreach (var path_index in indexToChangeListForTriCellCutTwo_vReduceTwice)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                // 去掉最后两竖
                cellDT.Branch(path)[index] = tCell.GenerateTwoRemovedEShapeOrZShapeOrEVariantShape();
            }

            foreach (var path_index in indexToChangeListForTriCellCutOne_vReduceTwice)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                // 随机去掉一竖
                int directionCode1 = m_random.Next(0, 4);// 取值范围是0,1,2,3
                cellDT.Branch(path)[index] = tCell.GenerateOneRemovedEShapeOrZShapeOrEVariantShape(directionCode1);
            }
            #endregion

            #region 其他不需要修改的
            // 不动
            #endregion
        }

        private void GenerateDensityConstrainedCell_hReduce(DataTree<Cell> cellDT,
                                                            
                                                            List<Tuple<GH_Path, int>> indexToChangeListForBiCellCutTwo_hReduce,
                                                            List<Tuple<GH_Path, int>> indexToChangeListForBiCellCutOne_hReduce,
                                                            List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutThree_hReduce,
                                                            List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutTwo_hReduce,
                                                            List<Tuple<GH_Path, int>> indexToChangeListForTriCellCutOne_hReduce,
                                                            List<Tuple<GH_Path, int>> indexNotChangeListForBiCell_hReduce, List<Tuple<GH_Path, int>> indexNotChangeListForTriCell_hReduce)
        {
            #region 对于需要修改的
            foreach (var path_index in indexToChangeListForBiCellCutTwo_hReduce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                // 去掉最后两横
                // int directionCode1 = m_random.Next(0, 4);// 取值范围是0,1,2,3
                cellDT.Branch(path)[index] = bCell.GenerateRemovedH(2, 0);
            }

            foreach (var path_index in indexToChangeListForBiCellCutOne_hReduce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                BilinearCell bCell = cellDT.Branch(path)[index] as BilinearCell;
                // 随机去掉一横
                int directionCode1 = m_random.Next(0, 3);// 取值范围是0,1,2
                cellDT.Branch(path)[index] = bCell.GenerateRemovedH(1, directionCode1);
            }

            foreach (var path_index in indexToChangeListForTriCellCutThree_hReduce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                cellDT.Branch(path)[index] = tCell.GenerateRemovedH(3, 0);
            }
            foreach (var path_index in indexToChangeListForTriCellCutTwo_hReduce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                int directionCode1 = m_random.Next(0, 3);// 取值范围是0,1,2
                cellDT.Branch(path)[index] = tCell.GenerateRemovedH(2, directionCode1);
            }
            foreach (var path_index in indexToChangeListForTriCellCutOne_hReduce)
            {
                GH_Path path = path_index.Item1;
                int index = path_index.Item2;

                TrilinearCell tCell = cellDT.Branch(path)[index] as TrilinearCell;
                int directionCode1 = m_random.Next(0, 3);// 取值范围是0,1,2
                cellDT.Branch(path)[index] = tCell.GenerateRemovedH(1, directionCode1);
            }
            #endregion

            #region 其他不需要修改的
            // 不动
            #endregion
        }
        private List<Polyline> SplitBlockIntoQuadBlock(Polyline offsetedPolyline,
                                                       
                                                       Polyline unOffsetedPolyline,
                                                       out List<double> sortedBrepAreas,

                                                       out List<List<Line>> lineForEachEdgeLoL,
                                                       out List<List<bool>> isGenerateablesLoL,
                                                       out List<bool> isSmallerList,
                                                       out Point3d publicPoint)
        {
            
            List<List<Point3d>> sortedPtsLoL = new List<List<Point3d>>();
            lineForEachEdgeLoL = new List<List<Line>>();
            isGenerateablesLoL = new List<List<bool>>();

            #region 生成分割后的两个sortedPts
            BoundingBox box = new BoundingBox(offsetedPolyline);
            double xMax = box.Max.X - box.Min.X;
            double yMax = box.Max.Y - box.Min.Y;
            double max = xMax > yMax ? xMax : yMax;

            List<Point3d> pts = new List<Point3d>();
            pts.AddRange(offsetedPolyline);
            pts.RemoveAt(pts.Count - 1);

            int turningIndex = 0;
            Line prevLine = Line.Unset;
            Line nextLine = Line.Unset;
            for (int i = 0; i < pts.Count; i++)
            {
                Vector3d vec0 = new Vector3d(pts[i] - pts[((i - 1) + pts.Count) % pts.Count]);
                Vector3d vec1 = new Vector3d(pts[(i + 1) % pts.Count] - pts[i]);

                if (Vector3d.CrossProduct(vec0, vec1).Z < 0)
                {
                    turningIndex = i;
                    prevLine = new Line(pts[((i - 1) + pts.Count) % pts.Count], pts[i]);
                    nextLine = new Line(pts[i], pts[(i + 1) % pts.Count]);
                    break;
                }
            }

            List<Brep> splitedBreps = new List<Brep>();

            isSmallerList = new List<bool>();
            publicPoint = Point3d.Unset;
            Point3d prevTurningPt = Point3d.Unset;
            Point3d turningPt = pts[turningIndex];
            if (prevLine != Line.Unset && nextLine != Line.Unset)
            {
                double length0 = prevLine.Length;
                double length1 = nextLine.Length;
                bool flag = length0 < length1 ? true : false;
                CurveIntersections intersections;
                //Point3d addPt = Point3d.Unset;
                //Point3d turningPt = pts[turningIndex];
                bool flag1 = false;
                bool flag2 = false;
                if (flag)
                {
                    prevLine.Extend(0, max);
                    prevTurningPt = prevLine.From;
                    intersections = Intersection.CurveCurve(prevLine.ToNurbsCurve(), offsetedPolyline.ToNurbsCurve(), 0.5 * GH_Component.DocumentTolerance(), 0.5 * GH_Component.DocumentTolerance());
                    flag2 = true;
                }
                else
                {
                    nextLine.Extend(max, 0);
                    prevTurningPt = nextLine.To;
                    intersections = Intersection.CurveCurve(nextLine.ToNurbsCurve(), offsetedPolyline.ToNurbsCurve(), 0.5 * GH_Component.DocumentTolerance(), 0.5 * GH_Component.DocumentTolerance());
                    flag2 = false;
                }

                for (int i = 0; i < intersections.Count; i++)
                {
                    if (intersections[i].IsPoint)
                    {
                        publicPoint = intersections[i].PointB;
                        flag1 = true;
                        break;
                    }
                }

                #region 得到 splitedBreps
                if (unOffsetedPolyline.Count == 6)
                {
                    #region 得到splitLine
                    Line innerLine;
                    if (flag)
                    {
                        prevLine.Extend(max, 0);
                        innerLine = prevLine;
                    }
                    else
                    {
                        nextLine.Extend(0, max);
                        innerLine = nextLine;
                    }

                    // 找出unOffsetedPolyline 的顶点中，在innerLine上的点
                    Point3d unOffsetedTurningPoint = Point3d.Unset;
                    int unOffsetedTurningIndex = 0;
                    for (int i = 0; i < unOffsetedPolyline.Count; i++)
                    {
                        if (Math.Abs(innerLine.DistanceTo(unOffsetedPolyline[i], false) - 0) < 0.1 * GH_Component.DocumentTolerance())
                        {
                            unOffsetedTurningPoint = unOffsetedPolyline[i];
                            unOffsetedTurningIndex = i;
                        }
                    }
                    // 从这个点做垂线来切割UnOffsetedPolyline
                    Line line = new Line(unOffsetedPolyline[unOffsetedTurningIndex], unOffsetedPolyline[(unOffsetedTurningIndex + 1) % unOffsetedPolyline.Count]);
                    Vector3d direction = line.Direction;
                    Vector3d vertical = new Vector3d(-direction.Y, direction.X, direction.Z);
                    vertical.Unitize();
                    vertical *= 2 * max;
                    Line splitLine = new Line(unOffsetedTurningPoint, unOffsetedTurningPoint + vertical);
                    #endregion

                    #region 得到 splitedBreps
                    //List<Point3d> unOffsetedPts = unOffsetedPolyline.ToList();
                    //List<Point3d> newOffsetedPts = new List<Point3d>();
                    //newOffsetedPts.AddRange(unOffsetedPts);
                    //unOffsetedPts.Add(unOffsetedPts[0]);
                    //Polyline newUnOffsetedPolyline = new Polyline(unOffsetedPts);
                    //intersections = Intersection.CurveCurve(splitLine.ToNurbsCurve(), newUnOffsetedPolyline.ToNurbsCurve(), 0.5 * GH_Component.DocumentTolerance(), 0.5 * GH_Component.DocumentTolerance());

                    //bool flag3 = false;
                    //Point3d unOffsetedPublicPoint = Point3d.Unset;
                    //for (int i = 0; i < intersections.Count; i++)
                    //{
                    //    if (intersections[i].IsPoint)
                    //    {
                    //        unOffsetedPublicPoint = intersections[i].PointB;
                    //        if (intersections[i].PointB != unOffsetedTurningPoint)
                    //        {
                    //            unOffsetedPublicPoint = intersections[i].PointB;
                    //            flag3 = true;
                    //            break;
                    //        }
                    //    }
                    //}

                    //if (flag3)
                    //{
                    //    List<Point3d> unOffsetedPts1 = new List<Point3d>();
                    //    List<Point3d> unOffsetedPts2 = new List<Point3d>();

                    //    newOffsetedPts.Insert(((unOffsetedTurningIndex - 2) + newOffsetedPts.Count) % newOffsetedPts.Count, unOffsetedPublicPoint);

                    //    int newUnOffsetedTurningIndex = newOffsetedPts.IndexOf(unOffsetedTurningPoint);

                    //    int currIndex = ()
                    //}

                    List<Point3d> closedUnOffsetedPts = unOffsetedPolyline.ToList();
                    List<Point3d> unClosedOffsetedPts = new List<Point3d>();
                    unClosedOffsetedPts.AddRange(closedUnOffsetedPts);
                    closedUnOffsetedPts.Add(closedUnOffsetedPts[0]);
                    Polyline closedUnOffsetedPolyline = new Polyline(closedUnOffsetedPts);

                    List<Curve> crvs = new List<Curve>();
                    crvs.Add(splitLine.ToNurbsCurve());
                    Brep brep = BoundarySurfaces(closedUnOffsetedPolyline.ToNurbsCurve())[0];
                    brep = brep.Faces[0].Split(crvs, 0.1 * GH_Component.DocumentTolerance());
                    
                    for (int i = 0; i < brep.Faces.Count; i++)
                    {
                        splitedBreps.Add(brep.Faces[i].DuplicateFace(false));
                    }
                    #endregion

                }
                else
                {
                    // unOffsetedPolyline.Count == 7
                    #region 得到splitLine
                    Line innerLine;
                    Vector3d vec;
                    if (flag)
                    {
                        prevLine.Extend(max, 0);
                        innerLine = prevLine;
                        vec = prevLine.Direction;
                        vec.Unitize();
                    }
                    else
                    {
                        nextLine.Extend(0, max);
                        innerLine = nextLine;
                        vec = nextLine.Direction;
                        vec.Unitize();
                    }

                    // 首先找出unOffsetedPolyline的每一段小边
                    List<Line> allEdges = new List<Line>();
                    for (int j = 0; j < unOffsetedPolyline.Count; j++)
                    {
                        allEdges.Add(new Line(unOffsetedPolyline[j], unOffsetedPolyline[(j + 1) % unOffsetedPolyline.Count]));
                    }
                    // 从小边中找出所有与innerLine平行的小边
                    List<Line> parallelEdges = new List<Line>();
                    for (int j = 0; j < allEdges.Count; j++)
                    {
                        Vector3d cur = allEdges[j].Direction;
                        cur.Unitize();
                        if (Math.Abs(Math.Abs(vec * cur) - 1) < 0.1 * GH_Component.DocumentTolerance())
                        {
                            parallelEdges.Add(allEdges[j]);
                        }
                    }
                    // 从所有平行的小边里找出距离最近的那一条，作为切割UnOffsetedPolyline的起始边
                    Line splitLine = parallelEdges[0];
                    double min = innerLine.DistanceTo(parallelEdges[0].From, false);
                    for (int j = 0; j < parallelEdges.Count; j++)
                    {
                        double distance = innerLine.DistanceTo(parallelEdges[j].From, false);
                        if (distance < min)
                        {
                            min = distance;
                            splitLine = parallelEdges[j];
                        }
                    }

                    splitLine.Extend(max, max);
                    #endregion

                    #region 得到 splitedBreps
                    List<Point3d> closedUnOffsetedPts = unOffsetedPolyline.ToList();
                    List<Point3d> unClosedOffsetedPts = new List<Point3d>();
                    unClosedOffsetedPts.AddRange(closedUnOffsetedPts);
                    closedUnOffsetedPts.Add(closedUnOffsetedPts[0]);
                    Polyline closedUnOffsetedPolyline = new Polyline(closedUnOffsetedPts);

                    List<Curve> crvs = new List<Curve>();
                    crvs.Add(splitLine.ToNurbsCurve());
                    Brep brep = BoundarySurfaces(closedUnOffsetedPolyline.ToNurbsCurve())[0];
                    brep = brep.Faces[0].Split(crvs, 0.1 * GH_Component.DocumentTolerance());

                    for (int i = 0; i < brep.Faces.Count; i++)
                    {
                        splitedBreps.Add(brep.Faces[i].DuplicateFace(false));
                    }
                    #endregion

                    //#region 得到 out List<Polyline> splietedUnOffsetedPolylines
                    //List<Point3d> unOffsetedPts = unOffsetedPolyline.ToList();
                    //List<Point3d> newOffsetedPts = new List<Point3d>();
                    //newOffsetedPts.AddRange(unOffsetedPts);
                    //unOffsetedPts.Add(unOffsetedPts[0]);
                    //Polyline newUnOffsetedPolyline = new Polyline(unOffsetedPts);

                    //int unOffsetedTurningIndex = 0;
                    //Point3d unOffsetedTurningPt = Point3d.Unset;
                    //bool flag3 = false;
                    //bool flag4 = false;
                    //if (flag)
                    //{
                    //    unOffsetedTurningIndex = newOffsetedPts.IndexOf(splitLine.From);

                    //    splitLine.Extend(0, max);
                    //    unOffsetedTurningPt = splitLine.From;
                    //    intersections = Intersection.CurveCurve(splitLine.ToNurbsCurve(), newUnOffsetedPolyline.ToNurbsCurve(), 0.5 * GH_Component.DocumentTolerance(), 0.5 * GH_Component.DocumentTolerance());
                    //    flag4 = true;
                    //}
                    //else
                    //{
                    //    unOffsetedTurningIndex = newOffsetedPts.IndexOf(splitLine.To);

                    //    splitLine.Extend(max, 0);
                    //    unOffsetedTurningPt = splitLine.To;
                    //    intersections = Intersection.CurveCurve(splitLine.ToNurbsCurve(), newUnOffsetedPolyline.ToNurbsCurve(), 0.5 * GH_Component.DocumentTolerance(), 0.5 * GH_Component.DocumentTolerance());
                    //    flag4 = false;
                    //}

                    //Point3d unOffsetedPublicPoint = Point3d.Unset;
                    //for (int i = 0; i < intersections.Count; i++)
                    //{
                    //    if (intersections[i].IsPoint)
                    //    {
                    //        unOffsetedPublicPoint = intersections[i].PointB;
                    //        flag3 = true;
                    //        break;
                    //    }
                    //}

                    //if (flag3)
                    //{
                    //    List<Point3d> unOffsetedPts1 = new List<Point3d>();
                    //    List<Point3d> unOffsetedPts2 = new List<Point3d>();

                    //    if (flag4)
                    //    {
                    //        newOffsetedPts.Insert(((unOffsetedTurningIndex - 1 - 2) + newOffsetedPts.Count) % newOffsetedPts.Count, unOffsetedPublicPoint);

                    //        int newUnOffsetedTurningIndex = newOffsetedPts.IndexOf(unOffsetedTurningPt);

                    //        int currIndex = ((newUnOffsetedTurningIndex - 1 - 3) + newOffsetedPts.Count) % newOffsetedPts.Count;
                    //        while (currIndex != newUnOffsetedTurningIndex)
                    //        {
                    //            unOffsetedPts1.Add(newOffsetedPts[currIndex]);
                    //            currIndex = (currIndex + 1) % newOffsetedPts.Count;
                    //        }

                    //        currIndex = newUnOffsetedTurningIndex;
                    //        while (currIndex != (newUnOffsetedTurningIndex + 1 + 3) % newOffsetedPts.Count)
                    //        {
                    //            unOffsetedPts2.Add(newOffsetedPts[currIndex]);
                    //            currIndex = (currIndex + 1) % newOffsetedPts.Count;
                    //        }
                    //    }
                    //    else
                    //    {
                    //        newOffsetedPts.Insert((unOffsetedTurningIndex - 2) + newOffsetedPts.Count, unOffsetedPublicPoint);

                    //        int newUnOffsetedTurningIndex = newOffsetedPts.IndexOf(unOffsetedTurningPt);

                    //        int currIndex = (newUnOffsetedTurningIndex + 1 + 3) % newOffsetedPts.Count;
                    //        while (currIndex != newUnOffsetedTurningIndex+1)
                    //        {
                    //            unOffsetedPts1.Add(newOffsetedPts[currIndex]);
                    //            currIndex = (currIndex + 1) % newOffsetedPts.Count;
                    //        }

                    //        currIndex = newUnOffsetedTurningIndex + 1;
                    //        while (currIndex != (newUnOffsetedTurningIndex+2+3)%newOffsetedPts.Count)
                    //        {
                    //            unOffsetedPts2.Add(newOffsetedPts[currIndex]);
                    //            currIndex = (currIndex + 1) % newOffsetedPts.Count;
                    //        }
                    //    }

                    //    //List<List<Point3d>> newOffsetedPtsLoL = new List<List<Point3d>>();
                    //    unOffsetedPts1.Add(unOffsetedPts1[0]);
                    //    unOffsetedPts2.Add(unOffsetedPts2[0]);
                    //    //newOffsetedPtsLoL.Add(unOffsetedPts1);
                    //    //newOffsetedPtsLoL.Add(unOffsetedPts2);

                    //    splietedUnOffsetedPolylines = new List<Polyline>();
                    //    splietedUnOffsetedPolylines.Add(new Polyline(unOffsetedPts1));
                    //    splietedUnOffsetedPolylines.Add(new Polyline(unOffsetedPts2));
                    //}
                    //else
                    //{
                    //    splietedUnOffsetedPolylines = null;
                    //}
                    //#endregion
                }
                #endregion

                if (flag1)
                {
                    List<Point3d> pts1 = new List<Point3d>();
                    List<Point3d> pts2 = new List<Point3d>();

                    if (flag2)
                    {
                        pts.Insert(((turningIndex - 1 - 2) + pts.Count) % pts.Count, publicPoint);

                        int newTurningIndex = pts.IndexOf(turningPt);

                        int currIndex = ((newTurningIndex - 1 - 3) + pts.Count) % pts.Count;
                        while (currIndex != newTurningIndex)
                        {
                            pts1.Add(pts[currIndex]);
                            currIndex = (currIndex + 1) % pts.Count;
                        }
                        isSmallerList.Add(false);

                        currIndex = newTurningIndex;
                        while (currIndex != (newTurningIndex + 1 + 3) % pts.Count)
                        {
                            pts2.Add(pts[currIndex]);
                            currIndex = (currIndex + 1) % pts.Count;
                        }
                        isSmallerList.Add(true);
                    }
                    else
                    {
                        pts.Insert((turningIndex - 2) + pts.Count, publicPoint);

                        int newTurningIndex = pts.IndexOf(turningPt);

                        int currIndex = (newTurningIndex + 1 + 3) % pts.Count;
                        while (currIndex != newTurningIndex + 1)
                        {
                            pts1.Add(pts[currIndex]);
                            currIndex = (currIndex + 1) % pts.Count;
                        }
                        isSmallerList.Add(true);

                        currIndex = newTurningIndex + 1;
                        while (currIndex != (newTurningIndex + 2 + 3) % pts.Count)
                        {
                            pts2.Add(pts[currIndex]);
                            currIndex = (currIndex + 1) % pts.Count;
                        }
                        isSmallerList.Add(false);
                    }

                    List<List<Point3d>> ptsLoL = new List<List<Point3d>>();
                    ptsLoL.Add(pts1);
                    ptsLoL.Add(pts2);
                    //polylines = new List<Polyline>();
                    for (int i = 0; i < ptsLoL.Count; i++)
                    {
                        // sort
                        BoundingBox box1 = new BoundingBox(ptsLoL[i]);
                        double minY = box1.Min.Y;
                        List<int> lowestPointIndexs = new List<int>();
                        List<Point3d> lowestPoints = new List<Point3d>();
                        for (int j = 0; j < ptsLoL[i].Count; j++)
                        {
                            if (Math.Abs(ptsLoL[i][j].Y - minY) < 0.01 * GH_Component.DocumentTolerance())
                            {
                                lowestPointIndexs.Add(j);
                                lowestPoints.Add(ptsLoL[i][j]);
                            }
                        }

                        int shift;
                        if (lowestPointIndexs.Count == 1)
                        {
                            shift = lowestPointIndexs[0];
                        }
                        else
                        {
                            List<double> xlist = new List<double>();
                            for (int j = 0; j < lowestPoints.Count; j++)
                            {
                                xlist.Add(lowestPoints[i].X);
                            }
                            int index = xlist.IndexOf(xlist.Min());

                            shift = lowestPointIndexs[index];
                        }

                        List<Point3d> sortedPts = UtilityFunctions.Shift<Point3d>(ptsLoL[i], shift, true);
                        //sortedPts.Add(sortedPts[0]);
                        sortedPtsLoL.Add(new List<Point3d>());
                        sortedPtsLoL[i].AddRange(sortedPts);
                        //polylines.Add(new Polyline(sortedPts));
                    }

                    //return polylines;
                }
                else
                {
                    sortedPtsLoL = null;
                }
            }
            else
            {
                sortedPtsLoL = null;
            }
            #endregion

            List<Polyline> polylines = new List<Polyline>();
            if (sortedPtsLoL != null)
            {
                #region 得到能够代表每条边的直线
                List<List<Line>> unsortedLineForEachEdgeLoL = new List<List<Line>>();
                for (int i = 0; i < sortedPtsLoL.Count; i++)
                {
                    List<Line> unsortedLineForEachEdge = new List<Line>();
                    int currentIndex = 0;
                    do
                    {
                        Point3d pt0 = sortedPtsLoL[i][currentIndex];
                        Point3d pt1 = sortedPtsLoL[i][(currentIndex + 1) % sortedPtsLoL[i].Count];

                        Line line = new Line(pt0, pt1);
                        unsortedLineForEachEdge.Add(line);
                        currentIndex = sortedPtsLoL[i].IndexOf(pt1);
                    } while (currentIndex != 0);

                    unsortedLineForEachEdgeLoL.Add(new List<Line>());
                    unsortedLineForEachEdgeLoL[i].AddRange(unsortedLineForEachEdge);
                }
                #endregion

                #region 判断sortedPts[0]这个点左侧的线(lineForEdge[-1])还是右侧的线(lineForEdge[0])是southBaseLine
                for (int i = 0; i < sortedPtsLoL.Count; i++)
                {
                    Vector3d vector0 = sortedPtsLoL[i][sortedPtsLoL[i].Count - 1] - sortedPtsLoL[i][0];
                    Vector3d vector1 = sortedPtsLoL[i][1] - sortedPtsLoL[i][0];

                    double rad0 = Vector3d.VectorAngle(vector0, -Vector3d.XAxis);
                    double rad1 = Vector3d.VectorAngle(vector1, Vector3d.XAxis);

                    if (rad0 > Math.PI / 4 && rad1 > Math.PI / 4)
                    {
                        lineForEachEdgeLoL.Add(new List<Line>());
                        lineForEachEdgeLoL[i].AddRange(UtilityFunctions.Shift<Line>(unsortedLineForEachEdgeLoL[i], 1, false));
                    }
                    else
                    {
                        if (rad0 < rad1)
                        {
                            lineForEachEdgeLoL.Add(new List<Line>());
                            lineForEachEdgeLoL[i].AddRange(UtilityFunctions.Shift<Line>(unsortedLineForEachEdgeLoL[i], 1, false));
                        }
                        else if (rad0 == rad1)
                        {
                            double length0 = vector0.Length;
                            double length1 = vector1.Length;
                            if (length0 >= length1)
                            {
                                lineForEachEdgeLoL.Add(new List<Line>());
                                lineForEachEdgeLoL[i].AddRange(UtilityFunctions.Shift<Line>(unsortedLineForEachEdgeLoL[i], 1, false));
                            }
                            else
                            {
                                lineForEachEdgeLoL.Add(new List<Line>());
                                lineForEachEdgeLoL[i].AddRange(unsortedLineForEachEdgeLoL[i]);
                            }
                        }
                        else
                        {
                            lineForEachEdgeLoL.Add(new List<Line>());
                            lineForEachEdgeLoL[i].AddRange(unsortedLineForEachEdgeLoL[i]);
                        }
                    }
                }
                #endregion


                #region isGenerateable 控制
                for (int i = 0; i < lineForEachEdgeLoL.Count; i++)
                {
                    List<bool> isGenerateable = new List<bool>();
                    for (int j = 0; j < lineForEachEdgeLoL[i].Count; j++)
                    {
                        if ((lineForEachEdgeLoL[i][j].From == publicPoint && lineForEachEdgeLoL[i][j].To == turningPt)
                        || (lineForEachEdgeLoL[i][j].To == publicPoint && lineForEachEdgeLoL[i][j].From == turningPt))
                        {
                            isGenerateable.Add(false);
                        }
                        else
                        {
                            isGenerateable.Add(true);
                        }
                    }
                    isGenerateablesLoL.Add(new List<bool>());
                    isGenerateablesLoL[i].AddRange(isGenerateable);
                }
                #endregion

                #region 基于lineForEachEdge，得到sortedBoundaryPolyline
                for (int i = 0; i < sortedPtsLoL.Count; i++)
                {
                    List<Point3d> lineStartAndEnd = new List<Point3d>();
                    lineStartAndEnd.Add(lineForEachEdgeLoL[i][0].From);
                    for (int j = 0; j < lineForEachEdgeLoL[i].Count; j++)
                    {
                        lineStartAndEnd.Add(lineForEachEdgeLoL[i][j].To);
                    }
                    polylines.Add(new Polyline(lineStartAndEnd));
                }
                #endregion


                List<Brep> sortedBreps = new List<Brep>();
                for (int i = 0; i < polylines.Count; i++)
                {
                    Point3d center = new Point3d(0, 0, 0);
                    for (int j = 0; j < polylines[i].Count-1; j++)
                    {
                        center += polylines[i][j];
                    }
                    center /= polylines[i].Count - 1;

                    //List<Curve> splitedBrepEdges = new List<Curve>();
                    for (int j = 0; j < splitedBreps.Count; j++)
                    {
                        Curve crv = Curve.JoinCurves(splitedBreps[j].DuplicateNakedEdgeCurves(true, false))[0];
                        //splitedBrepEdges.Add(crv);
                        PointContainment contain = crv.Contains(center, Plane.WorldXY,GH_Component.DocumentTolerance());
                        if (contain == PointContainment.Inside)
                        {
                            sortedBreps.Add(splitedBreps[j]);
                        }
                    }

                    //for (int j = 0; j < splitedBrepEdges.Count; j++)
                    //{
                    //    splitedBrepEdges[j].co
                    //}
                }

                sortedBrepAreas = new List<double>();
                for (int i = 0; i < sortedBreps.Count; i++)
                {
                    sortedBrepAreas.Add(sortedBreps[i].GetArea());
                }
            }
            else
            {
                lineForEachEdgeLoL = null;
                polylines = null;
                sortedBrepAreas = null;
            }



            return polylines;
        }

        private List<List<Cell>> GenerateVolumeBaseLine(double w, 
                                                              double W, 
                                                              double lMin, 
                                                              double lMax, 
                                                              double d, 
                                                              bool isJitter, 

                                                              bool isSmallerBlock,
                                                              Point3d turningPoint,

                                                              out List<List<List<Curve>>> allHCurvesLoL, 
                                                              out List<List<List<Curve>>> allVCurvesLoL, 
                                                              out List<List<List<Polyline>>> contourLoL, 
                                                              out List<List<List<Polyline>>> holeLoL, 

                                                              out Curve southBaseLineForShow,
                                                              out Curve northBaseLineForShow,
                                                              out Curve eastBaseLineForShow,
                                                              out Curve westBaseLineForShow,

                                                              out Curve southBoundaryLineForShow,
                                                              out Curve northBoundaryLineForShow,
                                                              out Curve eastBoundaryLineForShow,
                                                              out Curve westBoundaryLineForShow,

                                                              out List<Curve> gap,

                                                              out List<List<Curve>> cellBoundaryLoL,

                                                              //out DataTree<Curve> horizontalCellBoundary,
                                                              //out DataTree<Curve> verticalCellBoundary,

                                                              Polyline sortedBoundaryPolyline, 
                                                              List<Line> lineForEachEdge, 
                                                              List<bool> isGenerateables)
        {
            //List<List<List<Brep>>> breps = new List<List<List<Brep>>>();

            DataTree<Curve> hCurvesDT;
            DataTree<Curve> vCurvesDT;

            DataTree<Curve> horizontalCellBoundary;
            DataTree<Curve> verticalCellBoundary;


            Curve newBoundary;
            Curve singleRegion = GenerateLayoutLinesOnQuadBlock(sortedBoundaryPolyline,
                                                                lineForEachEdge,
                                                                isGenerateables,
                                                                w,
                                                                W,
                                                                lMin,
                                                                lMax,
                                                                d,
                                                                false,

                                                                isSmallerBlock,
                                                                turningPoint,

                                                                out southBaseLineForShow,
                                                                out northBaseLineForShow,
                                                                out eastBaseLineForShow,
                                                                out westBaseLineForShow,

                                                                out southBoundaryLineForShow,
                                                                out northBoundaryLineForShow,
                                                                out eastBoundaryLineForShow,
                                                                out westBoundaryLineForShow,

                                                                out gap,

                                                                out newBoundary,

                                                                out hCurvesDT,
                                                                out vCurvesDT,
                                                                
                                                                out horizontalCellBoundary,
                                                                out verticalCellBoundary);
            int turningPointIndex = sortedBoundaryPolyline.IndexOf(turningPoint);

            // out
            allHCurvesLoL = new List<List<List<Curve>>>();
            allVCurvesLoL = new List<List<List<Curve>>>();
            contourLoL = new List<List<List<Polyline>>>();
            holeLoL = new List<List<List<Polyline>>>();

            List<List<Cell>> cellLoL;
            if (singleRegion == null)
            {
                #region 由得到的baseLine构造生成初步的Cell对象，并且对初步生成的Cell对象进行随机修改，形成工，口，C
                cellLoL = GenerateBasicCells(hCurvesDT, 
                                        vCurvesDT, 

                                        horizontalCellBoundary,
                                        verticalCellBoundary,

                                        lineForEachEdge,
                                        turningPointIndex,

                                        isJitter,
                                        lMin, lMax, d, w);
                #endregion
            }
            else
            {
                if (singleRegion.IsClosed)
                {
                    #region 构造Cell对象
                    Curve[] segment = singleRegion.DuplicateSegments();

                    int pulicBaseLineIndex = isGenerateables.IndexOf(false);
                    int anotherBaseLineIndexRelatedToCutPoint = -1;
                    if (pulicBaseLineIndex != -1)
                    {
                        Line publicLine = lineForEachEdge[pulicBaseLineIndex];
                        if (publicLine.From.EpsilonEquals(turningPoint, GH_Component.DocumentTolerance()))
                        {
                            anotherBaseLineIndexRelatedToCutPoint = (pulicBaseLineIndex + 1) % lineForEachEdge.Count;
                        }
                        else
                        {
                            anotherBaseLineIndexRelatedToCutPoint = ((pulicBaseLineIndex - 1) + lineForEachEdge.Count) % lineForEachEdge.Count;
                        }
                    }

                    BilinearCell bCell = new BilinearCell(singleRegion, segment[0], segment[2], segment[3], segment[1], pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, -1, -1, turningPoint);
                    cellLoL = new List<List<Cell>>();
                    cellLoL.Add(new List<Cell>());
                    cellLoL[0].Add(bCell);
                    #endregion
                }
                else
                {
                    Curve[] crvs = singleRegion.DuplicateSegments();
                    Curve[] result = new Curve[1];

                    BilinearCell bCell;
                    int pulicBaseLineIndex = isGenerateables.IndexOf(false);
                    Line publicLine = lineForEachEdge[pulicBaseLineIndex];
                    if (publicLine.From.EpsilonEquals(turningPoint, GH_Component.DocumentTolerance()))
                    {
                        int anotherBaseLineIndexRelatedToCutPoint = ((pulicBaseLineIndex - 1) + lineForEachEdge.Count) % lineForEachEdge.Count;
                        int prev = ((pulicBaseLineIndex - 2) + lineForEachEdge.Count) % lineForEachEdge.Count;

                        Curve crv = lineForEachEdge[((pulicBaseLineIndex - 1) + lineForEachEdge.Count) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve prevCrv = lineForEachEdge[((pulicBaseLineIndex - 2) + lineForEachEdge.Count) % lineForEachEdge.Count].ToNurbsCurve();

                        Curve offsetCrv = crv.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        Curve offsetPrevCrv = prevCrv.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                        CurveIntersections event0 = Intersection.CurveCurve(offsetCrv, offsetPrevCrv, GH_Component.DocumentTolerance(), GH_Component.DocumentTolerance());
                        Curve newCrv = new Line(event0[0].PointA, offsetCrv.PointAtEnd).ToNurbsCurve();
                        Curve newPrevCrv = new Line(offsetPrevCrv.PointAtStart, event0[0].PointA).ToNurbsCurve();
                        
                        if (anotherBaseLineIndexRelatedToCutPoint == 0)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), newCrv, null, newPrevCrv, null, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, prev, -1, turningPoint);
                        }
                        else if (anotherBaseLineIndexRelatedToCutPoint == 1)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), newPrevCrv, null, null, newCrv, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, prev, -1, turningPoint);
                        }
                        else if (anotherBaseLineIndexRelatedToCutPoint == 2)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), null, newCrv, null, newPrevCrv, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, prev, -1, turningPoint);
                        }
                        else
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), null, newPrevCrv, newCrv, null, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, prev, -1, turningPoint);
                        }
                    }
                    else
                    {
                        int anotherBaseLineIndexRelatedToCutPoint = (pulicBaseLineIndex + 1) % lineForEachEdge.Count;
                        int next = (pulicBaseLineIndex + 2) % lineForEachEdge.Count;

                        Curve crv = lineForEachEdge[(pulicBaseLineIndex + 1) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve nextCrv = lineForEachEdge[(pulicBaseLineIndex + 2) % lineForEachEdge.Count].ToNurbsCurve();

                        Curve offsetCrv = crv.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        Curve offsetNextCrv = nextCrv.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                        CurveIntersections event0 = Intersection.CurveCurve(offsetCrv, offsetNextCrv, GH_Component.DocumentTolerance(), GH_Component.DocumentTolerance());
                        Curve newCrv = new Line(offsetCrv.PointAtStart, event0[0].PointA).ToNurbsCurve();
                        Curve newNextCrv = new Line(event0[0].PointA, offsetNextCrv.PointAtEnd).ToNurbsCurve();

                        if (anotherBaseLineIndexRelatedToCutPoint == 0)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), newCrv, null, null, newNextCrv, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, -1, next, turningPoint);
                        }
                        else if (anotherBaseLineIndexRelatedToCutPoint == 1)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), null, newNextCrv, null, newCrv, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, -1, next, turningPoint);
                        }
                        else if (anotherBaseLineIndexRelatedToCutPoint == 2)
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), null, newCrv, newNextCrv, null, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, -1, next, turningPoint);
                        }
                        else
                        {
                            bCell = new BilinearCell(sortedBoundaryPolyline.ToNurbsCurve(), newNextCrv, null, newCrv, null, pulicBaseLineIndex, anotherBaseLineIndexRelatedToCutPoint, -1, next, turningPoint);
                        }

                    }

                    cellLoL = new List<List<Cell>>();
                    cellLoL.Add(new List<Cell>());
                    cellLoL[0].Add(bCell);
                }
            }

            cellBoundaryLoL = new List<List<Curve>>();
            for (int i = 0; i < cellLoL.Count; i++)
            {
                cellBoundaryLoL.Add(new List<Curve>());
                for (int j = 0; j < cellLoL[i].Count; j++)
                {
                    cellBoundaryLoL[i].Add(cellLoL[i][j].CellBoundary);
                }
            }

            return cellLoL;
        }

        private List<int> RandomGetNum(int countForThree, int countForTwo)
        {
            List<int> numList = new List<int>();
            int beenChosenCountForThree = 0;
            int beenChosenCountForTwo = 0;

            for (int i = 0; i < countForThree+countForTwo; i++)
            {
                int n = m_random.Next(0, 2);
                if (n == 0)
                {
                    beenChosenCountForThree++;
                    if (beenChosenCountForThree > countForThree)
                    {
                        numList.Add(2);
                    }
                    else
                    {
                        numList.Add(3);
                    }
                }
                else
                {
                    beenChosenCountForTwo++;
                    if (beenChosenCountForTwo > countForTwo)
                    {
                        numList.Add(3);
                    }
                    else
                    {
                        numList.Add(2);
                    }
                }
            }

            return numList;
        }

        private List<List<Cell>> GenerateBasicCells(DataTree<Curve> hCurvesDT, 
                                               DataTree<Curve> vCurvesDT, 

                                               DataTree<Curve> horizontalCellBoundary,
                                               DataTree<Curve> verticalCellBoundary,

                                               List<Line> lineForEachEdge,
                                               int turningPointIndex,

                                               bool isJitter,
                                               double lMin, 
                                               double lMax, 
                                               double d, 
                                               double w)
        {
            
            List<List<Cell>> cellLoL = new List<List<Cell>>();
            #region 由得到的baseLine构造生成初步的Cell对象
            // 如果有两排以上的HCurve
            if (hCurvesDT.BranchCount > 1)
            {
                #region 得到3或2的排列组合
                List<List<int>> randomNumLoL = new List<List<int>>();
                int countOnPrevBranch = hCurvesDT.Branch(0).Count;

                List<int[]> countPairs = UtilityFunctions.GetAllPositiveIntegerSolution(3, 2, hCurvesDT.Branch(0).Count);
                int randomIndex = m_random.Next(countPairs.Count);
                int countForThree = countPairs[randomIndex][0];
                int countForTwo = countPairs[randomIndex][1];
                List<int> randomNum = RandomGetNum(countForThree, countForTwo);

                for (int i = 0; i < hCurvesDT.BranchCount; i++)
                {
                    randomNumLoL.Add(new List<int>());

                    //bool isSame = true;

                    //List<List<int>> randomNumLoL = new List<List<int>>();
                    if (hCurvesDT.Branch(i).Count == countOnPrevBranch)
                    {
                        if (isJitter)
                        {
                            randomNumLoL[i].AddRange(RandomGetNum(countForThree, countForTwo));
                        }
                        else
                        {
                            randomNumLoL[i].AddRange(randomNum);
                        }
                    }
                    else
                    {
                        List<int[]> countPairsNew = UtilityFunctions.GetAllPositiveIntegerSolution(3, 2, hCurvesDT.Branch(i).Count);
                        int randomIndexNew = m_random.Next(countPairs.Count);
                        int countForThreeNew = countPairs[randomIndex][0];
                        int countForTwoNew = countPairs[randomIndex][1];

                        randomNumLoL[i].AddRange(RandomGetNum(countForThreeNew, countForTwoNew));
                    }
                }
                

                
                #endregion

                for (int i = 0; i < hCurvesDT.BranchCount; i++)
                {
                    cellLoL.Add(new List<Cell>());

                    int randomNumListIndex = 0;
                    int j = 0;
                    while (randomNumListIndex != randomNumLoL[i].Count)
                    {
                        int current = randomNumLoL[i][randomNumListIndex];

                        if (randomNumLoL[i][randomNumListIndex] == 2)
                        {
                            // 两条线的情况
                            if (HBranchType == HorizontalVolumeBranchType.OnlyOneCenter)
                            {
                                Curve centerH = hCurvesDT.Branch(i)[j];
                                Curve westLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                {
                                    westLine = vCurvesDT.Branch(j)[2 * i];
                                }
                                Curve eastLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                {
                                    eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                }

                                //Curve southLine = new Line(westLine.PointAtStart, eastLine.PointAtStart).ToNurbsCurve();
                                //Curve northLine = new Line(westLine.PointAtEnd, eastLine.PointAtEnd).ToNurbsCurve();
                                Curve horizontalCellBoundary_South = lineForEachEdge[0].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                Curve horizontalCellBoundary_North = lineForEachEdge[2].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                Curve verticalCellBoundary_West = verticalCellBoundary.Branch(0)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                Curve verticalCellBoundary_East = verticalCellBoundary.Branch(0)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                BilinearCell bilinearCell = new BilinearCell(centerH, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                cellLoL[i].Add(bilinearCell);
                            }
                            else
                            {
                                Curve southLine = hCurvesDT.Branch(i)[j];
                                Curve northLine = null;
                                bool flag = false;
                                if (j+1 < hCurvesDT.Branch(i).Count)
                                {
                                    if (hCurvesDT.Branch(i)[j + 1] != null)
                                    {
                                        northLine = hCurvesDT.Branch(i)[j + 1];
                                        flag = false;
                                    }
                                    else
                                    {
                                        flag = true;
                                    }
                                }
                                else
                                {
                                    flag = true;
                                }

                                Curve westLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                {
                                    westLine = vCurvesDT.Branch(j)[2 * i];
                                }
                                Curve eastLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                {
                                    eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                }


                                Curve horizontalCellBoundary_South = horizontalCellBoundary.Branch(i)[2 * j].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                Curve horizontalCellBoundary_North;
                                if (flag)
                                {
                                    //northLine = new Line(westLine.PointAtEnd, eastLine.PointAtEnd).ToNurbsCurve();
                                    horizontalCellBoundary_North = lineForEachEdge[2].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                }
                                else
                                {
                                    horizontalCellBoundary_North = horizontalCellBoundary.Branch(i)[2 * (j + 1) + 1].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                }
                                

                                Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                                BilinearCell bilinearCell = new BilinearCell(southLine, northLine, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                cellLoL[i].Add(bilinearCell);
                            }
                        }
                        else
                        {
                            // 三条线的情况
                            Curve southLine = hCurvesDT.Branch(i)[j];
                            Curve middleLine = hCurvesDT.Branch(i)[j + 1];
                            Curve northLine = hCurvesDT.Branch(i)[j + 2];

                            Curve westLine = null;
                            if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                            {
                                westLine = vCurvesDT.Branch(j)[2 * i];
                            }
                            Curve eastLine = null;
                            if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                            {
                                eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                            }
                            Curve westLine1 = null;
                            if (vCurvesDT.ItemExists(new GH_Path(j + 1), 2 * i))
                            {
                                westLine1 = vCurvesDT.Branch(j + 1)[2 * i];
                            }
                            Curve eastLine1 = null;
                            if (vCurvesDT.ItemExists(new GH_Path(j + 1), 2 * i + 1))
                            {
                                eastLine1 = vCurvesDT.Branch(j + 1)[2 * i + 1];
                            }

                            Curve horizontalCellBoundary_South = horizontalCellBoundary.Branch(i)[2 * j].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                            Curve horizontalCellBoundary_North = horizontalCellBoundary.Branch(i)[2 * (j + 2) + 1].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                            Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i];
                            Curve verticalCellBoundary_West1 = verticalCellBoundary.Branch(j + 1)[4 * i];
                            Curve joined_West = Curve.JoinCurves(new Curve[2] { verticalCellBoundary_West, verticalCellBoundary_West1 })[0];
                            joined_West = joined_West.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                            Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3];
                            Curve verticalCellBoundary_East1 = verticalCellBoundary.Branch(j + 1)[4 * i + 3];
                            Curve joined_East = Curve.JoinCurves(new Curve[2] { verticalCellBoundary_East, verticalCellBoundary_East1 })[0];
                            joined_East = joined_East.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                            TrilinearCell trilinearCell = new TrilinearCell(southLine, middleLine, northLine, westLine, westLine1, eastLine, eastLine1,
                                                                            horizontalCellBoundary_South,
                                                                            horizontalCellBoundary_North,
                                                                            joined_West,
                                                                            joined_East,
                                                                            turningPointIndex);
                            cellLoL[i].Add(trilinearCell);
                        }

                        j += current;
                        randomNumListIndex++;
                    }
                }
            }
            // HCurve为1排或以下
            else
            {
                if (HBranchType == HorizontalVolumeBranchType.OnlyOneCenter)
                {
                    cellLoL.Add(new List<Cell>());

                    Curve centerH = hCurvesDT.Branch(0)[0];
                    Curve westLine = null;
                    if (vCurvesDT.ItemExists(new GH_Path(0), 0))
                    {
                        westLine = vCurvesDT.Branch(0)[0];
                    }
                    Curve eastLine = null;
                    if (vCurvesDT.ItemExists(new GH_Path(0), 1))
                    {
                        eastLine = vCurvesDT.Branch(0)[1];
                    }

                    Curve horizontalCellBoundary_South = lineForEachEdge[0].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                    Curve horizontalCellBoundary_North = lineForEachEdge[2].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                    Curve verticalCellBoundary_West = verticalCellBoundary.Branch(0)[0].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                    Curve verticalCellBoundary_East = verticalCellBoundary.Branch(0)[3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                    BilinearCell cell = new BilinearCell(centerH, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                    cellLoL[0].Add(cell);
                }
                else
                {
                    #region 得到3或2的排列组合
                    List<int[]> countPairs = UtilityFunctions.GetAllPositiveIntegerSolution(3, 2, hCurvesDT.Branch(0).Count);

                    int randomIndex = m_random.Next(countPairs.Count);
                    int countForThree = countPairs[randomIndex][0];
                    int countForTwo = countPairs[randomIndex][1];

                    //bool isSame = true;

                    List<List<int>> randomNumLoL = new List<List<int>>();
                    if (isJitter)
                    {
                        for (int i = 0; i < hCurvesDT.BranchCount; i++)
                        {
                            randomNumLoL.Add(new List<int>());
                            randomNumLoL[i].AddRange(RandomGetNum(countForThree, countForTwo));
                        }
                    }
                    else
                    {
                        List<int> randomNum = RandomGetNum(countForThree, countForTwo);
                        for (int i = 0; i < hCurvesDT.BranchCount; i++)
                        {
                            randomNumLoL.Add(new List<int>());
                            randomNumLoL[i].AddRange(randomNum);
                        }
                    }
                    #endregion

                    for (int i = 0; i < hCurvesDT.BranchCount; i++)
                    {
                        cellLoL.Add(new List<Cell>());

                        int randomNumListIndex = 0;
                        int j = 0;
                        while (randomNumListIndex != randomNumLoL[i].Count)
                        {
                            int current = randomNumLoL[i][randomNumListIndex];

                            if (randomNumLoL[i][randomNumListIndex] == 2)
                            {
                                // 两条线的情况
                                if (HBranchType == HorizontalVolumeBranchType.OnlyOneCenter)
                                {
                                    Curve centerH = hCurvesDT.Branch(i)[j];
                                    Curve westLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                    {
                                        westLine = vCurvesDT.Branch(j)[2 * i];
                                    }
                                    Curve eastLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                    {
                                        eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                    }

                                    Curve horizontalCellBoundary_South = lineForEachEdge[0].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve horizontalCellBoundary_North = lineForEachEdge[2].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    BilinearCell bilinearCell = new BilinearCell(centerH, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                    cellLoL[i].Add(bilinearCell);
                                }
                                else if (HBranchType == HorizontalVolumeBranchType.OnlyOneSouth)
                                {
                                    Curve southLine = hCurvesDT.Branch(i)[j];
                                    Curve northLine = null;
                                    Curve westLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                    {
                                        westLine = vCurvesDT.Branch(j)[2 * i];
                                    }
                                    Curve eastLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                    {
                                        eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                    }

                                    Curve horizontalCellBoundary_North = lineForEachEdge[2].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    Curve horizontalCellBoundary_South = horizontalCellBoundary.Branch(i)[2 * j].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    BilinearCell bilinearCell = new BilinearCell(southLine, northLine, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                    cellLoL[i].Add(bilinearCell);
                                }
                                else if (HBranchType == HorizontalVolumeBranchType.OnlyOneNorth)
                                {
                                    Curve southLine = null;
                                    Curve northLine = hCurvesDT.Branch(i)[j];
                                    Curve westLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                    {
                                        westLine = vCurvesDT.Branch(j)[2 * i];
                                    }
                                    Curve eastLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                    {
                                        eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                    }

                                    Curve horizontalCellBoundary_South = lineForEachEdge[0].ToNurbsCurve().Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    Curve horizontalCellBoundary_North = horizontalCellBoundary.Branch(i)[2 * j + 1].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    BilinearCell bilinearCell = new BilinearCell(southLine, northLine, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                    cellLoL[i].Add(bilinearCell);
                                }
                                else
                                {
                                    Curve southLine = hCurvesDT.Branch(i)[j];
                                    Curve northLine = null;
                                    if (hCurvesDT.ItemExists(new GH_Path(i), j + 1))
                                    {
                                        northLine = hCurvesDT.Branch(i)[j + 1];
                                    }
                                    Curve westLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                    {
                                        westLine = vCurvesDT.Branch(j)[2 * i];
                                    }
                                    Curve eastLine = null;
                                    if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                    {
                                        eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                    }

                                    Curve horizontalCellBoundary_South = horizontalCellBoundary.Branch(i)[2 * j].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve horizontalCellBoundary_North = horizontalCellBoundary.Branch(i)[2 * (j + 1) + 1].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                    Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                    Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                                    BilinearCell bilinearCell = new BilinearCell(southLine, northLine, westLine, eastLine, horizontalCellBoundary_South, horizontalCellBoundary_North, verticalCellBoundary_West, verticalCellBoundary_East, turningPointIndex);
                                    cellLoL[i].Add(bilinearCell);
                                }
                            }
                            else
                            {
                                // 三条线的情况
                                Curve southLine = hCurvesDT.Branch(i)[j];
                                Curve middleLine = hCurvesDT.Branch(i)[j + 1];
                                Curve northLine = hCurvesDT.Branch(i)[j + 2];

                                Curve westLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i))
                                {
                                    westLine = vCurvesDT.Branch(j)[2 * i];
                                }
                                Curve eastLine = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j), 2 * i + 1))
                                {
                                    eastLine = vCurvesDT.Branch(j)[2 * i + 1];
                                }
                                Curve westLine1 = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j + 1), 2 * i))
                                {
                                    westLine1 = vCurvesDT.Branch(j + 1)[2 * i];
                                }
                                Curve eastLine1 = null;
                                if (vCurvesDT.ItemExists(new GH_Path(j + 1), 2 * i + 1))
                                {
                                    eastLine1 = vCurvesDT.Branch(j + 1)[2 * i + 1];
                                }

                                Curve horizontalCellBoundary_South = horizontalCellBoundary.Branch(i)[2 * j].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                                Curve horizontalCellBoundary_North = horizontalCellBoundary.Branch(i)[2 * (j + 2) + 1].Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                                Curve verticalCellBoundary_West = verticalCellBoundary.Branch(j)[4 * i];
                                Curve verticalCellBoundary_West1 = verticalCellBoundary.Branch(j + 1)[4 * i];
                                Curve joined_West = Curve.JoinCurves(new Curve[2] { verticalCellBoundary_West, verticalCellBoundary_West1 })[0];
                                joined_West = joined_West.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);


                                Curve verticalCellBoundary_East = verticalCellBoundary.Branch(j)[4 * i + 3];
                                Curve verticalCellBoundary_East1 = verticalCellBoundary.Branch(j + 1)[4 * i + 3];
                                Curve joined_East = Curve.JoinCurves(new Curve[2] { verticalCellBoundary_East, verticalCellBoundary_East1 })[0];
                                joined_East = joined_East.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);

                                TrilinearCell trilinearCell = new TrilinearCell(southLine, middleLine, northLine, westLine, westLine1, eastLine, eastLine1,
                                                                            horizontalCellBoundary_South,
                                                                            horizontalCellBoundary_North,
                                                                            joined_West,
                                                                            joined_East,
                                                                            turningPointIndex);
                                cellLoL[i].Add(trilinearCell);
                            }

                            j += current;
                            randomNumListIndex++;
                        }
                    }
                }
                
            }
            #endregion

            #region 把面宽方向上过长（即长度大于lMax）的体量线，分解为多个
            for (int i = 0; i < cellLoL.Count; i++)
            {
                for (int j = 0; j < cellLoL[i].Count; j++)
                {
                    if (cellLoL[i][j].IsBilinearCell())
                    {
                        BilinearCell bCell = cellLoL[i][j] as BilinearCell;
                        bCell.LengthConstraint(lMin, lMax, d, w, m_random);
                    }
                    else
                    {
                        TrilinearCell tCell = cellLoL[i][j] as TrilinearCell;
                        tCell.LengthConstraint(lMin, lMax, d, w, m_random);
                    }
                }
            }
            #endregion

            return cellLoL;
        } 


        private List<List<BoundarySegment>> OffsetBS(List<List<BoundarySegment>> allFaceBS, 
                                                     List<Polyline> facePolylines,
                                                     out List<Polyline> innerResultPolylines)
                                                     //out DataTree<Curve> cuts,
                                                     //out DataTree<Curve> publicSides,
                                                     //out DataTree<Brep> cutsBrepDT,
                                                     //out List<Brep> resultBreps)
        {
            double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;

            //cuts = new DataTree<Curve>();
            //publicSides = new DataTree<Curve>();
            //cutsBrepDT = new DataTree<Brep>();
            //resultBreps = new List<Brep>();

            innerResultPolylines = new List<Polyline>();
            List<List<BoundarySegment>> newAllFaceBS = new List<List<BoundarySegment>>();
            for (int i = 0; i < allFaceBS.Count; i++)
            {
                BoundingBox box = new BoundingBox(facePolylines[i].ToList());
                double dx = box.Max.X - box.Min.X;
                double dy = box.Max.Y - box.Min.Y;
                double extendedValue = Math.Max(dx, dy) / 2;
                
                newAllFaceBS.Add(new List<BoundarySegment>());

                // 创建每条边的切割区域
                List<Curve> cutCurves = new List<Curve>();
                for (int j = 0; j < allFaceBS[i].Count; j++)
                {
                    List<Point3d> newPoints = new List<Point3d>();

                    if (allFaceBS[i][j].setback == 0)
                    {
                        // 如果退线距离是0而导致没有cutCurves的时候，cutCurves设置为null
                        cutCurves.Add(null);
                    }
                    else
                    {
                        if (allFaceBS[i][j].Lines.Count > 1)
                        {
                            // 当前BS非直线
                            List<Point3d> points = allFaceBS[i][j].Points;

                            List<Vector3d[]> directionPairs = new List<Vector3d[]>();
                            for (int k = 0; k < points.Count; k++)
                            {
                                int prev = ((k - 1) + points.Count) % points.Count;
                                int curr = k;
                                int next = (k + 1) % points.Count;
                                Vector3d[] directionPair = null;
                                if (curr == 0)
                                {
                                    Vector3d vector = new Vector3d(points[next] - points[curr]);
                                    vector.Unitize();
                                    directionPair = new Vector3d[2] { new Vector3d(0, 0, 0), vector };
                                }
                                else if (curr == points.Count - 1)
                                {
                                    Vector3d vector = new Vector3d(points[curr] - points[prev]);
                                    vector.Unitize();
                                    directionPair = new Vector3d[2] { vector, new Vector3d(0, 0, 0) };
                                }
                                else
                                {
                                    Vector3d vector0 = new Vector3d(points[curr] - points[prev]);
                                    Vector3d vector1 = new Vector3d(points[next] - points[curr]);
                                    vector0.Unitize();
                                    vector1.Unitize();
                                    directionPair = new Vector3d[2] { vector0, vector1 };
                                }
                                directionPairs.Add(directionPair);
                            }

                            List<Vector3d> moveDirections = new List<Vector3d>();
                            for (int k = 0; k < directionPairs.Count; k++)
                            {
                                Vector3d moveDirection;
                                if (directionPairs[k][0] == new Vector3d(0, 0, 0))
                                {
                                    moveDirection = Vector3d.CrossProduct(directionPairs[k][1], new Vector3d(0, 0, -1));
                                    moveDirection.Unitize();
                                    moveDirection = moveDirection * allFaceBS[i][j].setback;
                                    // 在这里添加零曲率的判断
                                    if (!allFaceBS[i][j].IsZeroCurvatureAtStart)
                                    {
                                        moveDirection = moveDirection - directionPairs[k][1] * extendedValue;
                                    }

                                    moveDirections.Add(moveDirection);
                                }
                                else if (directionPairs[k][1] == new Vector3d(0, 0, 0))
                                {
                                    moveDirection = Vector3d.CrossProduct(directionPairs[k][0], new Vector3d(0, 0, -1));
                                    moveDirection.Unitize();
                                    moveDirection = moveDirection * allFaceBS[i][j].setback;
                                    if (!allFaceBS[i][j].IsZeroCurvatureAtEnd)
                                    {
                                        moveDirection = moveDirection + directionPairs[k][0] * extendedValue;
                                    }

                                    moveDirections.Add(moveDirection);
                                }
                                else
                                {
                                    // 求角平分线
                                    if (Vector3d.CrossProduct(directionPairs[k][0], directionPairs[k][1]).Z < 0)
                                    {
                                        directionPairs[k][0] = directionPairs[k][0] * (-1);
                                        moveDirection = -(directionPairs[k][0] + directionPairs[k][1]);
                                        moveDirection.Unitize();

                                        directionPairs[k][0].Unitize();
                                        directionPairs[k][1].Unitize();

                                        double cosTheta = directionPairs[k][0] * directionPairs[k][1];

                                        double sinAlpha = Math.Sqrt((1 - cosTheta) / 2);

                                        moveDirection = moveDirection * allFaceBS[i][j].setback / sinAlpha;
                                    }
                                    else
                                    {
                                        directionPairs[k][0] = directionPairs[k][0] * (-1);
                                        moveDirection = directionPairs[k][0] + directionPairs[k][1];
                                        moveDirection.Unitize();

                                        directionPairs[k][0].Unitize();
                                        directionPairs[k][1].Unitize();

                                        double cosTheta = directionPairs[k][1] * directionPairs[k][0];

                                        double sinAlpha = Math.Sqrt((1 - cosTheta) / 2);

                                        moveDirection = moveDirection * allFaceBS[i][j].setback / sinAlpha;
                                    }
                                    moveDirections.Add(moveDirection);
                                }
                            }

                            List<Point3d> movedPoints = new List<Point3d>();
                            for (int k = 0; k < points.Count; k++)
                            {
                                Point3d movedPoint = points[k] + moveDirections[k];
                                movedPoints.Add(movedPoint);
                            }
                            movedPoints.Reverse();

                            newPoints.AddRange(points);
                            newPoints.AddRange(movedPoints);
                        }
                        else
                        {
                            // 当前BS是直线
                            List<Point3d> points = allFaceBS[i][j].Points;

                            Vector3d direction = points.Last() - points.First();
                            direction.Unitize();

                            Vector3d moveDirection = Vector3d.CrossProduct(direction, new Vector3d(0, 0, -1));
                            moveDirection.Unitize();
                            moveDirection = moveDirection * allFaceBS[i][j].setback;

                            Point3d newPoint0;
                            if (allFaceBS[i][j].IsZeroCurvatureAtStart)
                            {
                                newPoint0 = points[0] + moveDirection;
                            }
                            else
                            {
                                newPoint0 = points[0] + moveDirection - direction * extendedValue;
                            }

                            Point3d newPoint1;
                            if (allFaceBS[i][j].IsZeroCurvatureAtEnd)
                            {
                                newPoint1 = points[1] + moveDirection;
                            }
                            else
                            {
                                newPoint1 = points[1] + moveDirection + direction * extendedValue;
                            }

                            newPoints.AddRange(points);
                            newPoints.Add(newPoint1);
                            newPoints.Add(newPoint0);
                        }

                        newPoints.Add(newPoints[0]);
                        Curve cutCurve = Curve.CreateInterpolatedCurve(newPoints, 1);
                        cutCurves.Add(cutCurve);
                    }
                }

                //cuts.EnsurePath(i);
                //cuts.Branch(i).AddRange(cutCurves);

                facePolylines[i].Add(facePolylines[i][0]);
                Curve CurveToCut = Curve.CreateInterpolatedCurve(facePolylines[i], 1);
                CurveSets curveSets = new CurveSets(CurveToCut);
                foreach (Curve curve in cutCurves)
                {
                    if (curve != null)
                    {
                        foreach (Curve curve2 in curveSets.Generation0)
                        {
                            Curve[] crvs = Curve.CreateBooleanDifference(curve2, curve, GH_Component.DocumentTolerance());
                            curveSets.AppendNewCurves(crvs);
                        }
                        curveSets.NextGeneration();
                    }
                    
                }

                List<Curve> innerResultCurve = new List<Curve>();
                innerResultCurve.AddRange(curveSets.Generation0);
                Polyline resultPolyline;
                innerResultCurve[0].TryGetPolyline(out resultPolyline);
                innerResultPolylines.Add(resultPolyline);

                // 用求surface与surface之间公共边的方法，来得到最后的BS
                Brep innerResultBrep = BoundarySurfaces(innerResultCurve[0])[0];
                //resultBreps.Add(innerResultBrep);

                List<Curve> intersections = new List<Curve>();
                for (int j = 0; j < cutCurves.Count; j++)
                {
                    BrepCurveIntersectionSolveResults solveResults = null;
                    if (cutCurves[j] == null)
                    {
                        // 如果退线距离是0而导致没有cutCurves的时候，把原来的polyline作为intersection
                        intersections.Add(allFaceBS[i][j].PolylineCurve);
                    }
                    else
                    {
                        solveResults = ComputeBrepCurveIntersection(innerResultBrep, cutCurves[j]);

                        Curve[] curves = solveResults.Curves;
                        Curve[] joinedCurves = Curve.JoinCurves(curves);

                        if (joinedCurves.Length == 1)
                        {
                            intersections.Add(joinedCurves[0]);
                        }
                        else
                        {
                            // 不用null，避免下一步报错
                            intersections.Add(allFaceBS[i][j].PolylineCurve);
                        }
                    }
                    

                    
                }
                //publicSides.EnsurePath(i);
                //publicSides.Branch(i).AddRange(intersections);

                // 得到最后的BS
                for (int j = 0; j < intersections.Count; j++)
                {
                    if (intersections[j] == null)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Intersection failed");
                        return null;
                    }
                    else
                    {
                        Polyline polyline;
                        intersections[j].TryGetPolyline(out polyline);
                        BoundarySegment newBoundarySegment = new BoundarySegment(polyline.ToList(), allFaceBS[i][j].HIndex, allFaceBS[i][j].FIndex, allFaceBS[i][j].AdjacentFIndex);
                        newBoundarySegment.LocationValue = allFaceBS[i][j].LocationValue;
                        newBoundarySegment.setback = allFaceBS[i][j].setback;
                        newBoundarySegment.IsZeroCurvatureAtStart = allFaceBS[i][j].IsZeroCurvatureAtStart;
                        newBoundarySegment.IsZeroCurvatureAtEnd = allFaceBS[i][j].IsZeroCurvatureAtEnd;

                        newAllFaceBS[i].Add(newBoundarySegment);
                    }
                }
            }
            return newAllFaceBS;
        }

        /// <summary>
        /// 三点是否共线
        /// </summary>
        /// <param name="pt0"></param>
        /// <param name="pt1"></param>
        /// <param name="pt2"></param>
        /// <returns></returns>
        private bool IsThreePointsCollinear(Point3d pt0, Point3d pt1, Point3d pt2)
        {
            Vector3d dxdy01 = new Vector3d(pt1 - pt0);
            Vector3d dxdy12 = new Vector3d(pt2 - pt1);

            double k = dxdy01.X * dxdy12.Y - dxdy01.Y * dxdy12.X;
            if (Math.Abs(k) - 0 < 0.5 * GH_Component.DocumentTolerance())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public class BrepCurveIntersectionSolveResults
        {
            public Curve[] Curves { get; set; }
            public Point3d[] Points { get; set; }
        }

        public BrepCurveIntersectionSolveResults ComputeBrepCurveIntersection(Brep brp, Curve crv)
        {
            Curve[] curves = null;
            Point3d[] points = null;
            if (!Intersection.CurveBrep(crv, brp, GH_Component.DocumentTolerance(), out curves, out points))
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Intersection failed");
            }
            return new BrepCurveIntersectionSolveResults
            {
                Curves = curves,
                Points = points
            };
        }

        public Brep[] BoundarySurfaces(Curve curInput)
        {
            Brep[] result = null;

            Polyline polyline = ToPolyline(curInput);
            if (polyline != null)
            {
                
                int segmentCount = polyline.SegmentCount;
                if (segmentCount <= 4)
                {
                    if (segmentCount != 3)
                    {
                        if (segmentCount == 4)
                        {
                            Brep brep = Brep.CreateFromCornerPoints(polyline[0], polyline[1], polyline[2], polyline[3], GH_Component.DocumentTolerance());
                            if (brep != null)
                            {
                                return new Brep[] { brep };
                            }
                        }
                    }
                    else
                    {
                        Brep brep2 = Brep.CreateFromCornerPoints(polyline[0], polyline[1], polyline[2], GH_Component.DocumentTolerance());
                        if (brep2 != null)
                        {
                            return new Brep[] { brep2 };
                        }
                    }
                }
                else
                {
                    result = Brep.CreatePlanarBreps(curInput, GH_Component.DocumentTolerance());
                }
                
            }
            return result;
        }

        public Brep[] BoundarySurfaces(List<Curve> allCurves)
        {
            CurveList edges = new CurveList();
            foreach (var crv in allCurves)
            {
                if (crv != null)
                {
                    edges.Add(crv);
                }
            }

            Brep[] result;
            if (edges.Count == 0)
            {
                result = null;
            }
            else
            {
                if (edges.Count == 1)
                {
                    Polyline polyline = ToPolyline(edges[0]);
                    if (polyline != null)
                    {
                        int segmentCount = polyline.SegmentCount;
                        if (segmentCount != 3)
                        {
                            if (segmentCount == 4)
                            {
                                Brep brep = Brep.CreateFromCornerPoints(polyline[0], polyline[1], polyline[2], polyline[3], GH_Component.DocumentTolerance());
                                if (brep != null)
                                {
                                    return new Brep[] { brep };
                                }
                            }
                        }
                        else
                        {
                            Brep brep2 = Brep.CreateFromCornerPoints(polyline[0], polyline[1], polyline[2], GH_Component.DocumentTolerance());
                            if (brep2 != null)
                            {
                                return new Brep[] { brep2 };
                            }
                        }
                    }
                }
                result = Brep.CreatePlanarBreps(edges, GH_Component.DocumentTolerance());
            }
            return result;
        }

        private Polyline ToPolyline(Curve curve)
        {
            Polyline result;
            if (curve == null)
            {
                result = null;
            }
            else
            {
                Rhino.Geometry.Plane worldXY = Rhino.Geometry.Plane.WorldXY;
                if (!curve.TryGetPlane(out worldXY, 2.0 * Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance))
                {
                    result = null;
                }
                else
                {
                    Polyline polyline = null;
                    if (!curve.TryGetPolyline(out polyline))
                    {
                        result = null;
                    }
                    else if (!polyline.IsClosed)
                    {
                        result = null;
                    }
                    else
                    {
                        polyline.Transform(Transform.PlanarProjection(worldXY));
                        polyline.ReduceSegments(0.1 * Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);

                        if (polyline.SegmentCount < 3)
                        {
                            result = null;
                        }

                        result = polyline;
                        //else if (polyline.SegmentCount > 4)
                        //{
                        //    result = null;
                        //}
                        //else
                        //{
                        //    List<int> list = new List<int>();
                        //    Node2List node2List = new Node2List();
                        //    int num = polyline.Count - 1;
                        //    for (int i = 0; i < num; i++)
                        //    {
                        //        double nx;
                        //        double ny;
                        //        worldXY.ClosestParameter(polyline[i], out nx, out ny);
                        //        node2List.Append(new Node2(nx, ny));
                        //    }
                        //    if (!Solver.Compute(node2List, list))
                        //    {
                        //        result = null;
                        //    }
                        //    else if (list.Count != polyline.SegmentCount)
                        //    {
                        //        result = null;
                        //    }
                        //    else
                        //    {
                        //        result = polyline;
                        //    }
                        //}
                    }
                }
            }

            return result;
        }

        private List<BoundarySegment> SortBoundarySegment(List<BoundarySegment> bs, 
                                                          Polyline boundaryPolyline, 
                                                          out bool isHorizontalLayout,
                                                          out List<Line> lineForEachEdge,
                                                          out Polyline sortedBoundaryPolyline)
        {
            #region 得到从最南侧的点开始排序的sortedBoundaryPolyline
            List<Point3d> pts = boundaryPolyline.ToList();
            pts.RemoveAt(pts.Count - 1);
            BoundingBox boundingBox = boundaryPolyline.BoundingBox;
            double minY = boundingBox.Min.Y;
            List<int> lowestPointIndexs = new List<int>();
            List<Point3d> lowestPoints = new List<Point3d>();
            for (int i = 0; i < pts.Count; i++)
            {
                if (Math.Abs(pts[i].Y - minY) < 0.01 * GH_Component.DocumentTolerance())
                {
                    lowestPointIndexs.Add(i);
                    lowestPoints.Add(pts[i]);
                }
            }

            int shift;
            if (lowestPointIndexs.Count == 1)
            {
                // 只有一个点是最南侧的点时
                shift = lowestPointIndexs[0];
            }
            else
            {
                // 有两个点是最南侧的点时
                List<double> xlist = new List<double>();
                for (int i = 0; i < lowestPoints.Count; i++)
                {
                    xlist.Add(lowestPoints[i].X);
                }
                int index = xlist.IndexOf(xlist.Min());

                shift = lowestPointIndexs[index];
            }

            List<Point3d> sortedPts = UtilityFunctions.Shift<Point3d>(pts, shift, true);
            
            #endregion

            #region 得到能够代表每条边的直线
            List<Line> unsortedLineForEachEdge = new List<Line>();
            int currentIndex = 0;
            do
            {
                Point3d pt0 = sortedPts[currentIndex];
                Point3d pt1 = sortedPts[(currentIndex + 1) % sortedPts.Count];
                Point3d pt2 = sortedPts[(currentIndex + 2) % sortedPts.Count];
                if (IsThreePointsCollinear(pt0, pt1, pt2))
                {
                    Line line = new Line(pt0, pt2);
                    unsortedLineForEachEdge.Add(line);

                    currentIndex = sortedPts.IndexOf(pt2);
                }
                else
                {
                    Line line = new Line(pt0, pt1);
                    unsortedLineForEachEdge.Add(line);

                    currentIndex = sortedPts.IndexOf(pt1);
                }
            } while (currentIndex != 0);
            #endregion

            #region 判断sortedPts[0]这个点左侧的线(lineForEdge[-1])还是右侧的线(lineForEdge[0])是southBaseLine
            Vector3d vector0 = sortedPts[sortedPts.Count - 1] - sortedPts[0];
            Vector3d vector1 = sortedPts[1] - sortedPts[0];

            double rad0 = Vector3d.VectorAngle(vector0, -Vector3d.XAxis);
            double rad1 = Vector3d.VectorAngle(vector1, Vector3d.XAxis);

            isHorizontalLayout = false;
            if (rad0 > Math.PI/4 && rad1 > Math.PI/4)
            {
                isHorizontalLayout = true;
                lineForEachEdge = UtilityFunctions.Shift<Line>(unsortedLineForEachEdge, 1, false);
            }
            else
            {
                if (rad0 < rad1)
                {
                    isHorizontalLayout = false;
                    lineForEachEdge = UtilityFunctions.Shift<Line>(unsortedLineForEachEdge, 1, false);
                }
                else if (rad0 == rad1)
                {
                    isHorizontalLayout = false;

                    double length0 = vector0.Length;
                    double length1 = vector1.Length;
                    if (length0 >= length1)
                    {
                        lineForEachEdge = UtilityFunctions.Shift<Line>(unsortedLineForEachEdge, 1, false);
                    }
                    else
                    {
                        lineForEachEdge = unsortedLineForEachEdge;
                    }
                }
                else
                {
                    isHorizontalLayout = false;
                    lineForEachEdge = unsortedLineForEachEdge;
                }
            }
            #endregion

            #region 基于lineForEachEdge，得到sortedBoundaryPolyline
            List<Point3d> lineStartAndEnd = new List<Point3d>();
            lineStartAndEnd.Add(lineForEachEdge[0].From);
            for (int i = 0; i < lineForEachEdge.Count; i++)
            {
                lineStartAndEnd.Add(lineForEachEdge[i].To);
            }
            sortedBoundaryPolyline = new Polyline(lineStartAndEnd);
            #endregion


            #region 对BS进行排序
            int bsShift;
            List<int> jList = new List<int>();
            for (int i = 0; i < bs.Count; i++)
            {
                for (int j = 0; j < bs[i].Points.Count; j++)
                {
                    if (bs[i].Points[j].EpsilonEquals(lineForEachEdge[0].From, GH_Component.DocumentTolerance()))
                    {
                        jList.Add(j);
                    }
                }
            }
            bsShift = jList.IndexOf(jList.Min());
            List<BoundarySegment> sortedBS = UtilityFunctions.Shift<BoundarySegment>(bs, shift, true);

            #endregion

            return sortedBS;
        }

        private Curve GenerateLayoutLinesOnQuadBlock(//List<BoundarySegment> sortedBS, 
                                                     Polyline sortedBoundaryPolyline,
                                                     List<Line> lineForEachEdge,
                                                     List<bool> isGenerateables,
                                                     //bool isHorizontalLayout,
                                                     double w,
                                                     double W,
                                                     double lMin,
                                                     double lMax,
                                                     double d,
                                                     bool isParallelGeneration,

                                                     bool isSmallerBlock,
                                                     Point3d turningPoint,

                                                     out Curve southBaseLineForShow,
                                                     out Curve northBaseLineForShow,
                                                     out Curve eastBaseLineForShow,
                                                     out Curve westBaseLineForShow,

                                                     out Curve southBoundaryLineForShow,
                                                     out Curve northBoudnaryLineForShow,
                                                     out Curve eastBoundaryLineForShow,
                                                     out Curve westBoundaryLineForShow,

                                                     out List<Curve> gap,

                                                     out Curve boundary,
                                                     //out List<Curve> gapCenterLineForShow,
                                                     //out List<Curve> hCurves,
                                                     //out List<Curve> vCurves)
                                                     out DataTree<Curve> cuttedHorizontalCenterLines,
                                                     out DataTree<Curve> cuttedVerticalCenterLines,

                                                     out DataTree<Curve> horizontalCellBoundary,
                                                     out DataTree<Curve> verticalCellBoundary)
        {
            // 从确定南侧基线开始处理
            // 判断是选择起始点左侧的边还是右侧的边作为排布的基准线
            
            boundary = Curve.CreateInterpolatedCurve(sortedBoundaryPolyline, 1);

            #region 确定东南西北侧BoundaryLine
            #region 先计算东西侧，以此结果来修改 Boundary
            // 东侧
            int eastSegmentIndex = 1;
            Curve eastBoundaryLine;
            if (!isGenerateables[eastSegmentIndex])
            {
                eastBoundaryLine = lineForEachEdge[eastSegmentIndex].ToNurbsCurve();
                eastBoundaryLine = eastBoundaryLine.Offset(Plane.WorldXY, -(w+d), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                BoundingBox box = boundary.GetBoundingBox(Plane.WorldXY);
                double boxLength = box.Max.X - box.Min.X;
                eastBoundaryLine = eastBoundaryLine.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                eastBoundaryLine = TrimByBoundaryWithCurve(eastBoundaryLine, boundary);
                eastBoundaryLine.Domain = new Interval(0, 1);
            }
            else
            {
                eastBoundaryLine = lineForEachEdge[eastSegmentIndex].ToNurbsCurve();
            }

            // 西侧
            int westSegmentIndex = 3;
            Curve westBoundaryLine;
            if (!isGenerateables[westSegmentIndex])
            {
                westBoundaryLine = lineForEachEdge[westSegmentIndex].ToNurbsCurve();
                westBoundaryLine = westBoundaryLine.Offset(Plane.WorldXY, -(w+d), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                BoundingBox box = boundary.GetBoundingBox(Plane.WorldXY);
                double boxLength = box.Max.X - box.Min.X;
                westBoundaryLine = westBoundaryLine.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                westBoundaryLine = TrimByBoundaryWithCurve(westBoundaryLine, boundary);
                westBoundaryLine.Domain = new Interval(0, 1);
            }
            else
            {
                westBoundaryLine = lineForEachEdge[westSegmentIndex].ToNurbsCurve();
            }
            #endregion

            Curve southBoundaryLine = new Line(westBoundaryLine.PointAtEnd, eastBoundaryLine.PointAtStart).ToNurbsCurve();
            Curve northBoundaryLine = new Line(eastBoundaryLine.PointAtEnd, westBoundaryLine.PointAtStart).ToNurbsCurve();

            #region 再计算南北侧，一次结果来修改Boundary
            // 南侧
            int southSegmentIndex = 0;
            if (!isGenerateables[southSegmentIndex])
            {
                southBoundaryLine = southBoundaryLine.Offset(Plane.WorldXY, -W, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                BoundingBox box = boundary.GetBoundingBox(Plane.WorldXY);
                double boxLength = box.Max.X - box.Min.X;
                southBoundaryLine = southBoundaryLine.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                southBoundaryLine = TrimByBoundaryWithCurve(southBoundaryLine, boundary);
                southBoundaryLine.Domain = new Interval(0, 1);
            }

            // 北侧
            int northSegmentIndex = 2;
            if (!isGenerateables[northSegmentIndex])
            {
                northBoundaryLine = northBoundaryLine.ToNurbsCurve().Offset(Plane.WorldXY, -W, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];

                BoundingBox box = boundary.GetBoundingBox(Plane.WorldXY);
                double boxLength = box.Max.X - box.Min.X;
                northBoundaryLine = northBoundaryLine.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                northBoundaryLine = TrimByBoundaryWithCurve(northBoundaryLine, boundary);
                northBoundaryLine.Domain = new Interval(0, 1);
            }

            #endregion
            eastBoundaryLine = new Line(southBoundaryLine.PointAtEnd, northBoundaryLine.PointAtStart).ToNurbsCurve();
            westBoundaryLine = new Line(northBoundaryLine.PointAtEnd, southBoundaryLine.PointAtStart).ToNurbsCurve();

            List<Curve> segments = new List<Curve>();
            segments.Add(southBoundaryLine);
            segments.Add(eastBoundaryLine);
            segments.Add(northBoundaryLine);
            segments.Add(westBoundaryLine);

            boundary = Curve.JoinCurves(segments)[0];
            #endregion

            Vector3d directionVector;
            Vector3d verticalVector;

            #region 确定南北侧基线
            // 确定南侧基线
            directionVector = southBoundaryLine.TangentAtStart;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve southCenterLine = TrimBothEnd(boundary, verticalVector, southBoundaryLine, 1, w, W);
            southCenterLine.Domain = new Interval(0, 1);

            Curve southBaseLine = southCenterLine.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //sCrvForCalDistance = TrimBothEnd(boundary, verticalVector, sCrvForCalDistance, true, w, W);
            southBaseLine.Domain = new Interval(0, 1);

            // 确定北侧基线
            directionVector = northBoundaryLine.TangentAtStart;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve northCenterLine = TrimBothEnd(boundary, verticalVector, northBoundaryLine, 1, w, W);
            northCenterLine.Domain = new Interval(0, 1);
            Curve northCenterLineReverse = northCenterLine.DuplicateCurve();
            northCenterLineReverse.Reverse();

            Curve northBaseLine = northCenterLine.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //nCrvForCalDistance = TrimBothEnd(boundary, verticalVector, nCrvForCalDistance, true, w, W);
            northBaseLine.Domain = new Interval(0, 1);
            // 判断南北侧基线是否小于 lMin
            double sBLLength = southBaseLine.GetLength();
            double nBLLength = northBaseLine.GetLength();

            bool needToGenerateBaseLine = true;
            if (sBLLength < lMin && nBLLength < lMin)
            {
                // 如果南北两基线，都小于 lMin 那么直接让生成的形体，占满当前整个地块
                needToGenerateBaseLine = false;
            }
            else if (sBLLength < lMin && nBLLength >= lMin)
            {
                // 如果南基线的长度小于 lMin ，那么就让南基线进行退线，直到长度等于 lMin
                needToGenerateBaseLine = true;

                // todo
            }
            else if (sBLLength >= lMin && nBLLength < lMin)
            {
                // 如果北基线的长度小于 lMin ，那么就让北基线进行退线，直到长度等于 lMin
                needToGenerateBaseLine = true;

                // todo
            }
            #endregion

            #region 生成HorizontalCenterLines
            List<Curve> horizontalCenterLines = new List<Curve>();
            if (needToGenerateBaseLine)
            {
                // 确定是否需要在南北基线之间添加新的基线
                /* 此处不应该是 southCenterLine 和 northCenterLine 而应该是 偏移了 0.5*w 后的 */
                Point3d a1 = southBaseLine.PointAtStart;
                Point3d a2 = southBaseLine.PointAtEnd;
                Point3d b1 = northBaseLine.PointAtStart;
                Point3d b2 = northBaseLine.PointAtEnd;
                int indexBetween12BaseLines = 0;
                double minDistanceBetween12BaseLines = MinDistanceBetweenTwoLineSegment(a1, a2, b1, b2, out indexBetween12BaseLines);

                if (minDistanceBetween12BaseLines <= 0)
                {
                    // 如果两线段直接相交时，或者a1a2在b1b2上面时

                    b1 = northBoundaryLine.PointAtStart;
                    b2 = northBoundaryLine.PointAtEnd;

                    int indexBetweenBaseLineAndBoundaryLine = 0;
                    double minDistanceBetweenBaseLineAndBoundaryLine = MinDistanceBetweenTwoLineSegment(a1, a2, b1, b2, out indexBetweenBaseLineAndBoundaryLine);

                    if (minDistanceBetweenBaseLineAndBoundaryLine < W * 0.5)
                    {
                        // 如果southBaseLine到northBoundaryLine的距离，小于W后，就直接让体量充满地块
                        //horizontalCenterLines.Add(null);
                        HBranchType = HorizontalVolumeBranchType.OnlyOneFull;
                    }
                    else
                    {
                        // 此时场地内只有一个体量，因此要考虑在南侧，中间，还是北侧？
                        //Curve northCenterLineDP = northCenterLine.DuplicateCurve();
                        //northCenterLineDP.Reverse();

                        Surface surface = NurbsSurface.CreateRuledSurface(southCenterLine, northCenterLineReverse);
                        //List<double> pList = new List<double>() { 0.0, 0.5, 1.0 };
                        int pIndex = m_random.Next(3);
                        if (pIndex == 1)
                        {
                            double t = surface.Domain(1).ParameterAt(0.5);
                            Curve curve = surface.IsoCurve(0, t);
                            Vector3d vec = curve.TangentAtStart;
                            verticalVector = new Vector3d(-vec.Y, vec.X, vec.Z);
                            curve = TrimBothEnd(boundary, verticalVector, curve, 2, w, W);

                            horizontalCenterLines.Add(curve);
                            HBranchType = HorizontalVolumeBranchType.OnlyOneCenter;
                        }
                        else if (pIndex == 0)
                        {
                            horizontalCenterLines.Add(southCenterLine);
                            HBranchType = HorizontalVolumeBranchType.OnlyOneSouth;
                        }
                        else
                        {
                            horizontalCenterLines.Add(northCenterLineReverse);
                            HBranchType = HorizontalVolumeBranchType.OnlyOneNorth;
                        }
                    }

                }
                else if (minDistanceBetween12BaseLines > 0 && minDistanceBetween12BaseLines < W)
                {
                    // 两线段之间的距离不足以放置新的基线（进行TweenCurve）
                    // 尝试将NorthBaseLine缩短到使得 minDistanceBetweenSouthAndNorth = W

                    #region 计算setback
                    double setback;
                    if (indexBetween12BaseLines % 2 == 0)
                    {
                        // index为偶数，右端最近
                        if (indexBetween12BaseLines == 0)
                        {
                            // 与b1b2垂直
                            setback = Math.Sqrt(Dis2(a2, b1) - Math.Pow(minDistanceBetween12BaseLines, 2)) + Math.Sqrt(Math.Pow(W, 2) - Math.Pow(minDistanceBetween12BaseLines, 2));

                        }
                        else if (indexBetween12BaseLines == 2)
                        {
                            // 与a1a2垂直
                            double alpha = Vector3d.VectorAngle(-southBaseLine.TangentAtStart, northBaseLine.TangentAtStart);
                            setback = (W - minDistanceBetween12BaseLines) / Math.Sin(alpha);
                        }
                        else
                        {
                            // 都不垂直
                            double D = Dis2PointToStraightLine(b1, b2, a2);
                            setback = Math.Sqrt(Math.Pow(W, 2) - Math.Pow(D, 2)) - Math.Sqrt(Math.Pow(minDistanceBetween12BaseLines, 2) - Math.Pow(D, 2));
                        }
                    }
                    else
                    {
                        // index为奇数，左端最近
                        if (indexBetween12BaseLines == 1)
                        {
                            // 与b1b2垂直
                            setback = Math.Sqrt(Dis2(a1, b2) - Math.Pow(minDistanceBetween12BaseLines, 2)) + Math.Sqrt(Math.Pow(W, 2) - Math.Pow(minDistanceBetween12BaseLines, 2));
                        }
                        else if (indexBetween12BaseLines == 3)
                        {
                            // 与a1a2垂直
                            double alpha = Vector3d.VectorAngle(southBaseLine.TangentAtStart, -northBaseLine.TangentAtStart);
                            setback = (W - minDistanceBetween12BaseLines) / Math.Sin(alpha);
                        }
                        else
                        {
                            // 都不垂直
                            double D = Dis2PointToStraightLine(b1, b2, a1);
                            setback = Math.Sqrt(Math.Pow(W, 2) - Math.Pow(D, 2)) - Math.Sqrt(Math.Pow(minDistanceBetween12BaseLines, 2) - Math.Pow(D, 2));
                        }
                    }
                    #endregion

                    double l = northBaseLine.GetLength() - setback;
                    if (l < lMin)
                    {
                        // 此时体量的剩余距离过短，north体量被排除
                        // 此时场地内只有一个体量，因此要考虑在南侧，中间，还是北侧？
                        Curve northCenterLineDP = northCenterLine.DuplicateCurve();
                        northCenterLineDP.Reverse();

                        Surface surface = NurbsSurface.CreateRuledSurface(southCenterLine, northCenterLineDP);
                        //List<double> pList = new List<double>() { 0.0, 0.5, 1.0 };
                        int turningPtIndex = sortedBoundaryPolyline.IndexOf(turningPoint);
                        if (turningPtIndex != -1)
                        {
                            if (turningPtIndex == 0 || turningPtIndex == 1)
                            {
                                horizontalCenterLines.Add(southCenterLine);
                                HBranchType = HorizontalVolumeBranchType.OnlyOneSouth;
                            }
                            else
                            {
                                horizontalCenterLines.Add(northCenterLineDP);
                                HBranchType = HorizontalVolumeBranchType.OnlyOneNorth;
                            }
                        }
                        else
                        {
                            int pIndex = m_random.Next(3);
                            if (pIndex == 1)
                            {
                                double t = surface.Domain(1).ParameterAt(0.5);
                                Curve curve = surface.IsoCurve(0, t);
                                Vector3d vec = curve.TangentAtStart;
                                verticalVector = new Vector3d(-vec.Y, vec.X, vec.Z);
                                curve = TrimBothEnd(boundary, verticalVector, curve, 2, w, W);

                                horizontalCenterLines.Add(curve);
                                HBranchType = HorizontalVolumeBranchType.OnlyOneCenter;
                            }
                            else if (pIndex == 0)
                            {
                                horizontalCenterLines.Add(southCenterLine);
                                HBranchType = HorizontalVolumeBranchType.OnlyOneSouth;
                            }
                            else
                            {
                                horizontalCenterLines.Add(northCenterLineDP);
                                HBranchType = HorizontalVolumeBranchType.OnlyOneNorth;
                            }
                        }
                    }
                    else
                    {
                        // 此时场地内有两个体量
                        Curve newNorthCenterLine;
                        if (indexBetween12BaseLines % 2 == 0)
                        {
                            newNorthCenterLine = northCenterLine.Trim(CurveEnd.Start, setback);
                        }
                        else
                        {
                            newNorthCenterLine = northCenterLine.Trim(CurveEnd.End, setback);
                        }
                        newNorthCenterLine.Reverse();
                        newNorthCenterLine.Domain = new Interval(0, 1);

                        horizontalCenterLines.Add(southCenterLine);
                        horizontalCenterLines.Add(newNorthCenterLine);
                        HBranchType = HorizontalVolumeBranchType.TwoNorthShorten;
                    }
                }
                else
                {
                    // 两线段之间的距离足以放置新的基线（进行TweenCurve）
                    minDistanceBetween12BaseLines -= W;

                    int insertCount = (int)(Math.Floor(minDistanceBetween12BaseLines / (W + w)));

                    if (insertCount != 0)
                    {
                        // 南北基线中间插入新的基线
                        double sum = (insertCount + 1) * W + insertCount * w;
                        List<double> eachLengths = new List<double>();
                        eachLengths.Add(W + 0.5 * w);
                        for (int i = 1; i < insertCount; i++)
                        {
                            double num = eachLengths[i - 1] + W + w;
                            eachLengths.Add(num);
                        }

                        List<double> factors = new List<double>();
                        for (int i = 0; i < insertCount; i++)
                        {
                            double factor = eachLengths[i] / sum;
                            factors.Add(factor);
                        }

                        Curve trimedSouthBaseLine;
                        Curve trimedNorthBaseLine;
                        // 按照最短距离计算的结果，对进行TweenCurve的南北基线进行修正
                        // (a1, a2, b1):右端，与a1a2垂直，即index = 0
                        // (a1, a2, b2):左端，与a1a2垂直，即index = 1
                        // (b1, b2, a2):右端，与b1b2垂直，即index = 2
                        // (b1, b2, a1):左端，与b1b2垂直，即index = 3
                        // index = 4:右端，都不垂直
                        // index = 5:左端，都不垂直
                        if (indexBetween12BaseLines % 2 == 0)
                        {
                            // index为偶数，右端最近
                            if (indexBetween12BaseLines == 0)
                            {
                                // 与a1a2垂直
                                double t;
                                southBaseLine.ClosestPoint(b1, out t);
                                trimedNorthBaseLine = northBaseLine;
                                trimedSouthBaseLine = southBaseLine.Trim(0.0, t);
                            }
                            else if (indexBetween12BaseLines == 2)
                            {
                                // 与b1b2垂直
                                double t;
                                northBaseLine.ClosestPoint(a2, out t);
                                trimedSouthBaseLine = southBaseLine;
                                trimedNorthBaseLine = northBaseLine.Trim(t, 1.0);
                            }
                            else
                            {
                                // 都不垂直
                                trimedSouthBaseLine = southBaseLine;
                                trimedNorthBaseLine = northBaseLine;
                            }
                        }
                        else
                        {
                            // index为奇数，左端最近
                            if (indexBetween12BaseLines == 1)
                            {
                                // 与a1a2垂直
                                double t;
                                southBaseLine.ClosestPoint(b2, out t);
                                trimedNorthBaseLine = northBaseLine;
                                trimedSouthBaseLine = southBaseLine.Trim(t, 1.0);
                            }
                            else if (indexBetween12BaseLines == 3)
                            {
                                // 与b1b2垂直
                                double t;
                                northBaseLine.ClosestPoint(a1, out t);
                                trimedSouthBaseLine = southBaseLine;
                                trimedNorthBaseLine = northBaseLine.Trim(0.0, t);
                            }
                            else
                            {
                                // 都不垂直
                                trimedSouthBaseLine = southBaseLine;
                                trimedNorthBaseLine = northBaseLine;
                            }
                        }

                        if (isParallelGeneration)
                        {
                            // southCrv不需要修改赋值，重新对northCrv进行赋值即可
                            Vector3d motion;
                            switch (indexBetween12BaseLines)
                            {
                                case 0:
                                    motion = new Vector3d(northBaseLine.PointAtStart - southBaseLine.PointAtEnd);
                                    break;
                                case 1:
                                    motion = new Vector3d(northBaseLine.PointAtEnd - southBaseLine.PointAtStart);
                                    break;
                                case 2:
                                    motion = new Vector3d(northBaseLine.PointAtEnd - southBaseLine.PointAtStart);
                                    break;
                                case 3:
                                    motion = new Vector3d(northBaseLine.PointAtStart - southBaseLine.PointAtEnd);
                                    break;
                                default:
                                    motion = new Vector3d(0, 0, 0);
                                    break;
                            }

                            Transform move = Transform.Translation(motion);
                            Curve newNorthBaseLine = southBaseLine.DuplicateCurve();
                            newNorthBaseLine.Transform(move);
                            directionVector = newNorthBaseLine.TangentAtStart;
                            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                            trimedNorthBaseLine = TrimBothEnd(boundary, verticalVector, newNorthBaseLine, 1, w, W);
                        }

                        trimedNorthBaseLine.Reverse();
                        Surface surface = NurbsSurface.CreateRuledSurface(trimedSouthBaseLine, trimedNorthBaseLine);
                        List<Curve> insectCurves = new List<Curve>();
                        for (int i = 0; i < factors.Count; i++)
                        {
                            double t = surface.Domain(1).ParameterAt(factors[i]);
                            Curve curve = surface.IsoCurve(0, t);
                            insectCurves.Add(curve);
                        }
                        surface.Dispose();

                        horizontalCenterLines.Add(southCenterLine);
                        for (int i = 0; i < insectCurves.Count; i++)
                        {
                            directionVector = insectCurves[i].TangentAtStart;
                            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                            Curve insectCurve = TrimBothEnd(boundary, verticalVector, insectCurves[i], 2, w, W);
                            horizontalCenterLines.Add(insectCurve);
                        }
                        horizontalCenterLines.Add(northCenterLineReverse);
                        HBranchType = HorizontalVolumeBranchType.AboveThree;
                    }
                    else
                    {
                        // 南北基线中间不插入新的基线
                        horizontalCenterLines.Add(southCenterLine);
                        horizontalCenterLines.Add(northCenterLineReverse);
                        HBranchType = HorizontalVolumeBranchType.Two;
                    }
                }
            }
            #endregion

            #region 确定东西侧基线
            // 确定东侧基线
            directionVector = eastBoundaryLine.TangentAtStart;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve eastCenterLine = TrimBothEnd(boundary, verticalVector, eastBoundaryLine, 1, w, W);
            eastCenterLine.Domain = new Interval(0, 1);

            Curve eastBaseLine = eastBoundaryLine.Offset(Plane.WorldXY, -(lMin + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //if (!isGenerateables[eastSegmentIndex])
            //{
            //    eastBaseLine = eastBoundaryLine.Offset(Plane.WorldXY, -(lMin + w + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //}
            //else
            //{
            //    eastBaseLine = eastBoundaryLine.Offset(Plane.WorldXY, -(lMin + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //}
            
            eastBaseLine = eastBaseLine.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
            eastBaseLine = TrimByBoundaryWithCurve(eastBaseLine, boundary);
            eastBaseLine.Domain = new Interval(0, 1);

            // 确定西侧基线
            directionVector = westBoundaryLine.TangentAtStart;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve westCenterLine = TrimBothEnd(boundary, verticalVector, westBoundaryLine, 1, w, W);
            westCenterLine.Domain = new Interval(0, 1);

            Curve westBaseLine = westBoundaryLine.Offset(Plane.WorldXY, -(lMin + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //if (!isGenerateables[westSegmentIndex])
            //{
            //    westBaseLine = westBoundaryLine.Offset(Plane.WorldXY, -(lMin + w + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //}
            //else
            //{
            //    westBaseLine = westBoundaryLine.Offset(Plane.WorldXY, -(lMin + w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //}
            westBaseLine = westBaseLine.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
            westBaseLine = TrimByBoundaryWithCurve(westBaseLine, boundary);
            westBaseLine.Domain = new Interval(0, 1);
            #endregion

            #region 生成VerticalGapLines
            List<Curve> verticalCenterLinesForGap = new List<Curve>();
            // 确定是否需要在南北基线之间添加新的基线
            /* 此处不应该是 eastCenterLine 和 westCenterLine 而应该是 偏移了 0.5*w 后的 */
            Point3d a3 = eastBaseLine.PointAtStart;
            Point3d a4 = eastBaseLine.PointAtEnd;
            Point3d b3 = westBaseLine.PointAtStart;
            Point3d b4 = westBaseLine.PointAtEnd;
            int indexBetween34BaseLines = 0;
            double minDistanceBetween34BaseLines = MinDistanceBetweenTwoLineSegment(a3, a4, b3, b4, out indexBetween34BaseLines);
            // dw用于求Gap偏移后形成的CenterLines
            double dw = d + w;

            // 这一部分的目标就是保证lMin
            int indexBetween34BoundaryLines = 0;
            double minDistanceBetween34BoundaryLines = MinDistanceBetweenTwoLineSegment(eastBoundaryLine.PointAtStart, eastBoundaryLine.PointAtEnd, westBoundaryLine.PointAtStart, westBoundaryLine.PointAtEnd, out indexBetween34BoundaryLines);

            if (minDistanceBetween34BaseLines < dw)
            {
                /* 已检查完成 */
                // 此时需要需要检查 minDistanceBetween34BoundaryLines的大小
                if (minDistanceBetween34BoundaryLines > lMax)
                {
                    // 即此时已经长过了lMax，需要再判断，能否至少分成 2 个 lMin 体量
                    if (minDistanceBetween34BoundaryLines < 2 * lMin + dw)
                    {
                        // verticalCenterLinesForGap 不用添加新的 GapCurves
                        // 等待下一步真正构造体量时，进行长度检查即可，长度缩减直到满足lMax的约束即可
                        VBranchType = VerticalGapBranchType.ZeroGapNeedShorten;
                    }
                    else
                    {
                        // verticalCenterLinesForGap 添加 1 条新的 GapCurves
                        double tMin = (lMin + 0.5 * dw) / minDistanceBetween34BoundaryLines;
                        double tMax = 1 - (lMin + 0.5 * dw) / minDistanceBetween34BoundaryLines;
                        double t = m_random.NextDouble() * (tMax - tMin) + tMin;
                        //double t = 0.5;

                        Curve westBoundaryLineDP = westBoundaryLine.DuplicateCurve();
                        westBoundaryLineDP.Reverse();
                        Surface surface = NurbsSurface.CreateRuledSurface(eastBoundaryLine, westBoundaryLineDP);
                        Curve gapCrv = surface.IsoCurve(0, t);

                        verticalCenterLinesForGap.Add(gapCrv);
                        surface.Dispose();
                        VBranchType = VerticalGapBranchType.OneGap;
                    }
                }
                else
                {
                    // verticalCenterLinesForGap 不用添加新的 GapCurves
                    // 等待下一步真正构造体量时，进行长度检查即可，长度不用缩减，能够直接满足约束
                    VBranchType = VerticalGapBranchType.ZeroGap;
                }
            }
            // 两线段之间的距离不足以放置新的基线（进行TweenCurve）
            // minDistanceBetween34BaseLines的长度，小于放置1个lMin体量及其中间的间距，还有整体与baseLine之间的两个 0.5 * dw 间距
            // minDistanceBetween34BaseLines的长度，大于等于 dw ，即等于 dw 时只有两个体量，每部分距离为 lMin + w，dw，lMin + w
            else if (minDistanceBetween34BaseLines >= dw && minDistanceBetween34BaseLines < (lMin + dw + dw))
            {
                // minDistanceBetween34BaseLines的距离不足以放置 1 个 lMin
                if (minDistanceBetween34BaseLines > lMax)
                {
                    // minDistanceBetween34BoundaryLines的长度，一定也大于lMax，所以要看minDistanceBetween34BoundaryLines能放几个lMax，几个lMin
                    int volumeMinCount = (int)Math.Floor((minDistanceBetween34BoundaryLines - lMax) / (lMax + 0.5 * dw)) + 1;
                    int volumeMaxCount = (int)Math.Floor((minDistanceBetween34BoundaryLines - lMin) / (lMin + 0.5 * dw)) + 1;
                    int volumeCount = 0;
                    if (volumeMinCount != volumeMaxCount)
                    {
                        volumeCount = m_random.Next(volumeMinCount, volumeMaxCount);
                    }
                    else
                    {
                        volumeCount = volumeMinCount;
                    }

                    int gapCount = volumeCount - 1;

                    List<Curve> insectCurves = ConstructTweenCurve(gapCount, lMin, lMax, dw, w, eastBoundaryLine, westBoundaryLine, Point3d.Unset, Point3d.Unset, false, false, false);
                    verticalCenterLinesForGap.AddRange(insectCurves);
                    VBranchType = VerticalGapBranchType.UncertainGap;
                }
                else if (minDistanceBetween34BoundaryLines < lMin)
                {
                    // 如果minDistanceBetween34BoundaryLines < lMin，那么A4B3要进行setback，直到其长度等于 2 * lMin + dw 或 lMax 这两者中的较小值
                    // 如果是lMax更小，那么不用添加新的 GapCurves
                    // 如果是2 * lMin + dw更小，那么添加新的GapCurves
                    // lMax 与 2 * lMin + dw 同样小时，那么选择 2 * lMin + dw
                    double targetLength = 0;
                    bool needToAddGap = false;
                    if (2 * lMin + dw <= lMax)
                    {
                        targetLength = 2 * lMin + dw;
                        needToAddGap = true;
                    }
                    else
                    {
                        targetLength = lMax;
                        needToAddGap = false;
                    }

                    Vector3d a3a4 = new Vector3d(eastBoundaryLine.PointAtEnd - eastBoundaryLine.PointAtStart);
                    Vector3d a4b3 = new Vector3d(westBoundaryLine.PointAtStart - eastBoundaryLine.PointAtEnd);
                    Vector3d b3b4 = new Vector3d(westBoundaryLine.PointAtEnd - westBoundaryLine.PointAtStart);
                    Vector3d b4a3 = new Vector3d(eastBoundaryLine.PointAtStart - westBoundaryLine.PointAtEnd);
                    if (a4b3.Length < b4a3.Length)
                    {
                        if (targetLength > a4b3.Length && targetLength < b4a3.Length)
                        {
                            if (needToAddGap)
                            {
                                // 计算两线段之间的最短距离那一侧的边，要向最长距离那一侧的边 setback 多少距离
                                Curve trimedEastBoundaryLine;
                                Curve trimedWestBoundaryLine;
                                Point3d ePoint = Point3d.Unset;
                                Point3d wPoint = Point3d.Unset;

                                double setbackE = 0;
                                double setbackW = CalSetback(a3a4, a4b3, b3b4, a4b3.Length, targetLength, out setbackE);

                                trimedEastBoundaryLine = eastBoundaryLine.Trim(CurveEnd.End, setbackE);
                                trimedWestBoundaryLine = westBoundaryLine.Trim(CurveEnd.Start, setbackW);

                                if (trimedEastBoundaryLine != null)
                                {
                                    trimedEastBoundaryLine.Domain = new Interval(0, 1);
                                }
                                else
                                {
                                    ePoint = eastBoundaryLine.PointAtStart;
                                }
                                if (trimedWestBoundaryLine != null)
                                {
                                    trimedWestBoundaryLine.Domain = new Interval(0, 1);
                                }
                                else
                                {
                                    wPoint = westBoundaryLine.PointAtEnd;
                                }

                                // 基于能够放置基线的位置，来对原有的westBaseLine进行缩短，进而形成TweenCurve式的基线
                                List<Curve> insectGapCurves = ConstructTweenCurve(1, lMin, lMax, dw, w, trimedEastBoundaryLine, trimedWestBoundaryLine, ePoint, wPoint, false, false, false);
                                verticalCenterLinesForGap.AddRange(insectGapCurves);
                                VBranchType = VerticalGapBranchType.OneGap;
                            }
                            else
                            {
                                // verticalCenterLinesForGap 不用添加新的 GapCurves
                                // 等待下一步真正构造体量时，进行长度检查即可，A4B3会被剔除掉，以满足约束条件
                                VBranchType = VerticalGapBranchType.ZeroGapNorthCull;
                            }
                        }
                        else
                        {
                            // verticalCenterLinesForGap 不用添加新的 GapCurves
                            // 等待下一步真正构造体量时，进行长度检查即可，A4B3会被剔除掉，以满足约束条件
                            VBranchType = VerticalGapBranchType.ZeroGapNorthCull;
                        }
                    }
                    else if (a4b3.Length > b4a3.Length)
                    {
                        if (targetLength < a4b3.Length && targetLength > b4a3.Length)
                        {
                            if (needToAddGap)
                            {
                                // 计算两线段之间的最短距离那一侧的边，要向最长距离那一侧的边 setback 多少距离
                                Curve trimedEastBoundaryLine;
                                Curve trimedWestBoundaryLine;
                                Point3d ePoint = Point3d.Unset;
                                Point3d wPoint = Point3d.Unset;

                                double setbackE = 0;
                                double setbackW = CalSetback(a3a4, a4b3, b3b4, a4b3.Length, targetLength, out setbackE);

                                trimedEastBoundaryLine = eastBoundaryLine.Trim(CurveEnd.End, setbackE);
                                trimedWestBoundaryLine = westBoundaryLine.Trim(CurveEnd.Start, setbackW);

                                if (trimedEastBoundaryLine != null)
                                {
                                    trimedEastBoundaryLine.Domain = new Interval(0, 1);
                                }
                                else
                                {
                                    ePoint = eastBoundaryLine.PointAtStart;
                                }
                                if (trimedWestBoundaryLine != null)
                                {
                                    trimedWestBoundaryLine.Domain = new Interval(0, 1);
                                }
                                else
                                {
                                    wPoint = westBoundaryLine.PointAtEnd;
                                }

                                // 基于能够放置基线的位置，来对原有的westBaseLine进行缩短，进而形成TweenCurve式的基线
                                List<Curve> insectGapCurves = ConstructTweenCurve(1, lMin, lMax, dw, w, trimedEastBoundaryLine, trimedWestBoundaryLine, ePoint, wPoint, false, false, false);
                                verticalCenterLinesForGap.AddRange(insectGapCurves);
                                VBranchType = VerticalGapBranchType.OneGap;
                            }
                            else
                            {
                                // verticalCenterLinesForGap 不用添加新的 GapCurves
                                // 等待下一步真正构造体量时，进行长度检查即可，A4B3会被剔除掉，以满足约束条件
                                VBranchType = VerticalGapBranchType.ZeroGapNorthCull;
                            }
                        }
                        else
                        {
                            // verticalCenterLinesForGap 不用添加新的 GapCurves
                            // 等待下一步真正构造体量时，进行长度检查即可，A4B3会被剔除掉，以满足约束条件
                            VBranchType = VerticalGapBranchType.ZeroGapNorthCull;
                        }
                    }
                    else
                    {
                        // a4b3.Length == b4a3.Length
                        // verticalCenterLinesForGap 不用添加新的 GapCurves
                        // 等待下一步真正构造体量时，进行长度检查即可，A4B3会被剔除掉，以满足约束条件
                        VBranchType = VerticalGapBranchType.ZeroGapNorthCull;
                    }
                }
                else
                {
                    // verticalCenterLinesForGap 不用添加新的 GapCurves
                    // 等待下一步真正构造体量时，进行长度检查即可，长度不用缩减，能够直接满足约束
                    VBranchType = VerticalGapBranchType.ZeroGap;
                }
            }
            // 两线段之间的距离足以放置新的基线（进行TweenCurve）
            else
            {
                /* 已检查完成 */
                minDistanceBetween34BaseLines -= dw;
                int insertLMinCount = (int)(Math.Floor(minDistanceBetween34BaseLines / (lMin + dw)));
                // 因为是求gap的位置，所以是体量的数量 - 1
                int insertGapCount = insertLMinCount - 1;

                #region 修正BoundaryLine
                // insertLMinCount 必定 > 0
                // 东西基线中间插入新的分割线
                Curve trimedEastBoundaryLine;
                Curve trimedWestBoundaryLine;
                Point3d ePoint = Point3d.Unset;
                Point3d wPoint = Point3d.Unset;
                // 按照最短距离计算的结果，对进行TweenCurve的东西基线进行修正
                // (a3,a4,b3):上端，与a3a4垂直，即index = 0
                // (a3,a4,b4):下端，与a3a4垂直，即index = 1
                // (b3,b4,a3):上端，与b3b4垂直，即index = 2
                // (b3,b4,a4):下端，与b3b4垂直，即index = 3
                // index = 4:上端，都不垂直
                // index = 5:左端，都不垂直
                if (indexBetween34BaseLines % 2 == 0)
                {
                    // index为偶数，上端最近
                    if (indexBetween34BaseLines == 0)
                    {
                        // 与a3a4垂直
                        double t;
                        eastBoundaryLine.ClosestPoint(b3, out t);
                        trimedWestBoundaryLine = westBoundaryLine;
                        trimedEastBoundaryLine = eastBoundaryLine.Trim(0.0, t);
                    }
                    else if (indexBetween34BaseLines == 2)
                    {
                        // 与b3b4垂直
                        double t;
                        westBoundaryLine.ClosestPoint(a4, out t);
                        trimedEastBoundaryLine = eastBoundaryLine;
                        trimedWestBoundaryLine = westBoundaryLine.Trim(t, 1.0);
                    }
                    else
                    {
                        // 都不垂直
                        trimedEastBoundaryLine = eastBoundaryLine;
                        trimedWestBoundaryLine = westBoundaryLine;
                    }
                }
                else
                {
                    // index为奇数，左端最近
                    if (indexBetween34BaseLines == 1)
                    {
                        // 与a3a4垂直
                        double t;
                        eastBoundaryLine.ClosestPoint(b4, out t);
                        trimedWestBoundaryLine = westBoundaryLine;
                        trimedEastBoundaryLine = eastBoundaryLine.Trim(t, 1.0);
                    }
                    else if (indexBetween34BaseLines == 3)
                    {
                        // 与b3b4垂直
                        double t;
                        westBoundaryLine.ClosestPoint(a3, out t);
                        trimedEastBoundaryLine = eastBoundaryLine;
                        trimedWestBoundaryLine = westBoundaryLine.Trim(0.0, t);
                    }
                    else
                    {
                        // 都不垂直
                        trimedEastBoundaryLine = eastBoundaryLine;
                        trimedWestBoundaryLine = westBoundaryLine;
                    }
                }
                #endregion

                List<Curve> insectCurves = ConstructTweenCurve(insertGapCount, lMin, lMax, dw, w, trimedEastBoundaryLine, trimedWestBoundaryLine, ePoint, wPoint, true, false, false);
                verticalCenterLinesForGap.AddRange(insectCurves);
                VBranchType = VerticalGapBranchType.AboveOneGap;
            }
            #endregion

            #region 基于Gap线对horizontalCenterLines做切割，生成cuttedHorizontalCenterLines，另外生成垂直方向上的建筑布局线cuttedVerticalCenterLines
            #region 调整GapLines的顺序，并且补充东西两侧边缘处的切割线
            // 将GapLine的顺序，由自东向西排列，改为自西向东排列
            verticalCenterLinesForGap.Reverse();
            Curve westBoundaryLineReverse = westBoundaryLine.DuplicateCurve();
            westBoundaryLineReverse.Reverse();
            // 向verticalCenterLinesForGap中添加东西两侧的边界
            verticalCenterLinesForGap.Insert(0, westBoundaryLineReverse);
            verticalCenterLinesForGap.Add(eastBoundaryLine);

            List<Curve[]> verticalCenterLinePairs = new List<Curve[]>();
            if (horizontalCenterLines.Count == 1)
            {
                for (int i = 0; i < verticalCenterLinesForGap.Count; i++)
                {
                    Curve[] centerLinePairs;
                    if (i == 0)
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(false, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        //directionVector = centerLinePairs[1].TangentAtStart;
                        //verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                        centerLinePairs[1] = TrimBothEnd(boundary, centerLinePairs[1], 0, w);
                    }
                    else if (i == verticalCenterLinesForGap.Count - 1)
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(false, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        //directionVector = centerLinePairs[0].TangentAtStart;
                        //verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                        centerLinePairs[0] = TrimBothEnd(boundary, centerLinePairs[0], 0, w);
                        //centerLinePairs[0] = TrimBothEnd(boundary, verticalVector, centerLinePairs[0], 2, w, W);
                    }
                    else
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(true, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        for (int j = 0; j < centerLinePairs.Length; j++)
                        {
                            //directionVector = centerLinePairs[j].TangentAtStart;
                            //verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                            centerLinePairs[j] = TrimBothEnd(boundary, centerLinePairs[j], 0, w);
                        }
                    }

                    verticalCenterLinePairs.Add(centerLinePairs);
                }
            }
            else
            {
                for (int i = 0; i < verticalCenterLinesForGap.Count; i++)
                {
                    Curve[] centerLinePairs;
                    if (i == 0)
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(false, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        //directionVector = centerLinePairs[1].TangentAtStart;
                        //verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                        centerLinePairs[1] = TrimBothEnd(boundary, centerLinePairs[1], 0, w);
                    }
                    else if (i == verticalCenterLinesForGap.Count - 1)
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(false, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        //verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                        //centerLinePairs[0] = TrimBothEnd(boundary, verticalVector, centerLinePairs[0], 2, w, W);
                        centerLinePairs[0] = TrimBothEnd(boundary, centerLinePairs[0], 0, w);
                    }
                    else
                    {
                        centerLinePairs = GenerateCenterLineNearGapCurve(true, verticalCenterLinesForGap[i], d, w, W, boundary);
                        centerLinePairs[0] = new Line(centerLinePairs[0].PointAtStart, centerLinePairs[0].PointAtEnd).ToNurbsCurve();
                        centerLinePairs[1] = new Line(centerLinePairs[1].PointAtStart, centerLinePairs[1].PointAtEnd).ToNurbsCurve();
                        for (int j = 0; j < centerLinePairs.Length; j++)
                        {
                            directionVector = centerLinePairs[j].TangentAtStart;
                            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
                            centerLinePairs[j] = TrimBothEnd(boundary, centerLinePairs[j], 0, w);
                        }
                    }

                    verticalCenterLinePairs.Add(centerLinePairs);
                }
            }
            #endregion

            // Temp
            List<Curve> gapPairs = new List<Curve>();
            for (int i = 0; i < verticalCenterLinePairs.Count; i++)
            {
                for (int j = 0; j < verticalCenterLinePairs[i].Length; j++)
                {
                    gapPairs.Add(verticalCenterLinePairs[i][j]);
                }
            }
            gap = gapPairs;

            #region 计算t值，并且得到水平方向上被切割后的centerLines
            cuttedHorizontalCenterLines = new DataTree<Curve>();
            List<List<double>> verticalTLoL = new List<List<double>>();
            for (int i = 0; i < horizontalCenterLines.Count; i++)
            {
                cuttedHorizontalCenterLines.EnsurePath(i);
                verticalTLoL.Add(new List<double>());

                //List<double[]> horizontalTPairs = new List<double[]>();
                //List<double[]> verticalTPairs = new List<double[]>();
                List<double> horizontalTPairs = new List<double>();
                List<double> verticalTPairs = new List<double>();
                for (int j = 0; j < verticalCenterLinesForGap.Count; j++)
                {
                    Curve[] centerLinePairs = verticalCenterLinePairs[j];
                    double[] horizontalTPair = new double[2] { -1, -1 };
                    double[] verticalTPair = new double[2] { -1, -1 };
                    for (int k = 0; k < centerLinePairs.Length; k++)
                    {
                        CurveIntersections curveIntersections = Intersection.CurveCurve(horizontalCenterLines[i], centerLinePairs[k], GH_Component.DocumentTolerance(), GH_Component.DocumentTolerance());

                        if (curveIntersections.Count == 1)
                        {
                            Interval overlapA = curveIntersections[0].OverlapA;
                            Interval overlapB = curveIntersections[0].OverlapB;
                            horizontalTPair[k] = overlapA.T0;
                            verticalTPair[k] = overlapB.T0;
                        }
                    }
                    //List<double> horizontalTPairList = horizontalTPair.ToList();
                    //horizontalTPairList.Sort();
                    //List<double> verticalTPairList = verticalTPair.ToList();
                    //verticalTPairList.Sort();

                    //horizontalTPairs.Add(horizontalTPairList.ToArray());
                    //verticalTPairs.Add(verticalTPairList.ToArray());
                    horizontalTPairs.AddRange(horizontalTPair);
                    verticalTPairs.AddRange(verticalTPair);
                }
                verticalTLoL[i].AddRange(verticalTPairs);

                #region 得到水平方向上的建筑布局线
                double t0 = horizontalCenterLines[i].Domain.T0;
                double t1 = horizontalCenterLines[i].Domain.T1;

                List<double> horizontalT = new List<double>(horizontalTPairs);
                horizontalT.RemoveAll(a => a == -1);
                if (horizontalT.Count == 0)
                {
                    horizontalT.Add(t0);
                    horizontalT.Add(t1);
                }

                if (horizontalT.Count % 2 != 0)
                {
                    //horizontalT.RemoveAt(horizontalT.Count - 1);
                    horizontalT.Add(t1);
                }

                int index = 0;
                while (index != horizontalT.Count)
                {
                    Curve cut = horizontalCenterLines[i].Trim(horizontalT[index], horizontalT[index + 1]);
                    if (cut == null)
                    {
                        cut = horizontalCenterLines[i].Trim(horizontalT[index + 1], horizontalT[index]);
                        if (cut != null)
                        {
                            cuttedHorizontalCenterLines.Branch(i).Add(cut);
                        }
                        else
                        {
                            // 不向cuttedHorizontalCenterLines中添加
                        }
                    }
                    else
                    {
                        cuttedHorizontalCenterLines.Branch(i).Add(cut);
                    }

                    index += 2;
                }

                #endregion
            }
            #endregion

            #region 得到垂直方向上的建筑布局线
            cuttedVerticalCenterLines = new DataTree<Curve>();
            if (verticalTLoL.Count != 0)
            {
                List<List<double>> verticalTLoLTranspose = UtilityFunctions.Transpose<double>(verticalTLoL);
                // 从西到东排列的 VerticalCenterLines
                List<Curve> verticalCenterLines = new List<Curve>();
                for (int i = 0; i < verticalCenterLinePairs.Count; i++)
                {
                    verticalCenterLines.AddRange(verticalCenterLinePairs[i]);
                }

                for (int i = 0; i < verticalCenterLines.Count; i++)
                {
                    cuttedVerticalCenterLines.EnsurePath(i);

                    double t0 = verticalCenterLines[i].Domain.T0;
                    double t1 = verticalCenterLines[i].Domain.T1;
                    List<double> verticalT = verticalTLoLTranspose[i];
                    verticalT.RemoveAll(a => a == -1);

                    if (verticalT.Count != 0)
                    {
                        if (verticalT.Count == 1)
                        {
                            verticalT.Insert(0, t0);
                            verticalT.Add(t1);
                        }

                        if (HBranchType == HorizontalVolumeBranchType.OnlyOneSouth)
                        {
                            int index = 1;
                            while (index + 1 != verticalT.Count)
                            {
                                Curve cut = verticalCenterLines[i].Trim(verticalT[index], verticalT[index + 1]);
                                cuttedVerticalCenterLines.Branch(i).Add(cut);
                                index += 1;
                            }
                        }
                        else if (HBranchType == HorizontalVolumeBranchType.OnlyOneNorth)
                        {
                            int index = 0;
                            while (index + 1 != verticalT.Count - 1)
                            {
                                Curve cut = verticalCenterLines[i].Trim(verticalT[index], verticalT[index + 1]);
                                cuttedVerticalCenterLines.Branch(i).Add(cut);
                                index += 1;
                            }
                        }
                        else if (HBranchType == HorizontalVolumeBranchType.TwoNorthShorten)
                        {
                            if (verticalT.Count == 3)
                            {
                                int index = 1;
                                while (index + 1 != verticalT.Count)
                                {
                                    Curve cut = verticalCenterLines[i].Trim(verticalT[index], verticalT[index + 1]);
                                    cuttedVerticalCenterLines.Branch(i).Add(cut);
                                    index += 1;
                                }
                            }
                            else
                            {
                                int index = 0;
                                while (index + 1 != verticalT.Count)
                                {
                                    Curve cut = verticalCenterLines[i].Trim(verticalT[index], verticalT[index + 1]);
                                    cuttedVerticalCenterLines.Branch(i).Add(cut);
                                    index += 1;
                                }
                            }
                        }
                        else
                        {
                            int index = 0;
                            while (index + 1 != verticalT.Count)
                            {
                                Curve cut = verticalCenterLines[i].Trim(verticalT[index], verticalT[index + 1]);
                                cuttedVerticalCenterLines.Branch(i).Add(cut);
                                index += 1;
                            }
                        }
                    }
                }
            }
            else
            {
                // 不向cuttedVerticalCenterLines中添加
            }


            if (HBranchType == HorizontalVolumeBranchType.OnlyOneCenter)
            {
                for (int i = 0; i < cuttedVerticalCenterLines.BranchCount; i++)
                {
                    Curve[] crvs = Curve.JoinCurves(cuttedVerticalCenterLines.Branch(i).ToArray());
                    cuttedVerticalCenterLines.Branch(i).Clear();
                    cuttedVerticalCenterLines.Branch(i).AddRange(crvs);
                }
            }

            // cuttedHorizontalCenterLines需要进行行列转换
            cuttedHorizontalCenterLines = UtilityFunctions.FlipMatrix<Curve>(cuttedHorizontalCenterLines);

            // cuttedVerticalCenterLines需要进行行列转换
            List<int> emptyBranchIndex = new List<int>();
            for (int i = 0; i < cuttedVerticalCenterLines.BranchCount; i++)
            {
                if (cuttedVerticalCenterLines.Branch(i).Count == 0)
                {
                    emptyBranchIndex.Add(i);
                }
            }
            for (int i = 0; i < emptyBranchIndex.Count; i++)
            {
                cuttedVerticalCenterLines.RemovePath(new GH_Path(emptyBranchIndex[i]));
            }
            cuttedVerticalCenterLines.RenumberPaths();
            cuttedVerticalCenterLines = UtilityFunctions.FlipMatrix<Curve>(cuttedVerticalCenterLines);

            #endregion
            #endregion

            #region 得到水平方向上centerLines偏移后的，成对的horizontalCellBoundary
            horizontalCellBoundary = new DataTree<Curve>();
            for (int i = 0; i < cuttedHorizontalCenterLines.BranchCount; i++)
            {
                horizontalCellBoundary.EnsurePath(i);
                for (int j = 0; j < cuttedHorizontalCenterLines.Branch(i).Count; j++)
                {
                    if (cuttedHorizontalCenterLines.Branch(i)[j] != null)
                    {
                        Curve crv0 = cuttedHorizontalCenterLines.Branch(i)[j].Offset(Plane.WorldXY, 0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        horizontalCellBoundary.Branch(i).Add(crv0);
                    }
                    if (cuttedHorizontalCenterLines.Branch(i)[j] != null)
                    {
                        Curve crv1 = cuttedHorizontalCenterLines.Branch(i)[j].Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        horizontalCellBoundary.Branch(i).Add(crv1);
                    }

                    //crv0 = crv0.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                    //crv1 = crv1.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                }
            }
            #endregion

            #region 得到垂直方向上centerLines偏移后的，成对的VerticalCellBoundary
            verticalCellBoundary = new DataTree<Curve>();
            for (int i = 0; i < cuttedVerticalCenterLines.BranchCount; i++)
            {
                verticalCellBoundary.EnsurePath(i);
                for (int j = 0; j < cuttedVerticalCenterLines.Branch(i).Count; j++)
                {
                    if (cuttedVerticalCenterLines.Branch(i)[j] != null)
                    {
                        Curve crv0 = cuttedVerticalCenterLines.Branch(i)[j].Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        verticalCellBoundary.Branch(i).Add(crv0);
                    }
                    if (cuttedVerticalCenterLines.Branch(i)[j] !=null)
                    {
                        Curve crv1 = cuttedVerticalCenterLines.Branch(i)[j].Offset(Plane.WorldXY, 0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        verticalCellBoundary.Branch(i).Add(crv1);
                    }
                    

                    //crv0 = crv0.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                    //crv1 = crv1.Extend(CurveEnd.Both, 5 * w, CurveExtensionStyle.Line);
                }
            }
            #endregion


            // Temp
            southBaseLineForShow = southBaseLine;
            northBaseLineForShow = northBaseLine;
            eastBaseLineForShow = eastBaseLine;
            westBaseLineForShow = westBaseLine;
            //gapCenterLineForShow = verticalCenterLinesForGap;

            southBoundaryLineForShow = southBoundaryLine;
            northBoudnaryLineForShow = northBoundaryLine;
            eastBoundaryLineForShow = eastBoundaryLine;
            westBoundaryLineForShow = westBoundaryLine;


            // 输出
            if (cuttedHorizontalCenterLines.DataCount == 0)
            {
                if (isSmallerBlock)
                {
                    int index = isGenerateables.IndexOf(false);
                    Line publicLine = lineForEachEdge[index];
                    if (publicLine.From.EpsilonEquals(turningPoint,GH_Component.DocumentTolerance()))
                    {
                        Curve crv = lineForEachEdge[(index + 1) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve nextCrv = lineForEachEdge[(index + 2) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve join = Curve.JoinCurves(new Curve[2] { crv, nextCrv })[0];
                        join = join.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        return join;
                    }
                    else
                    {
                        Curve crv = lineForEachEdge[((index - 1)+ lineForEachEdge.Count) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve prevCrv = lineForEachEdge[((index - 2) + lineForEachEdge.Count) % lineForEachEdge.Count].ToNurbsCurve();
                        Curve join = Curve.JoinCurves(new Curve[2] { crv, prevCrv })[0];
                        join = join.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                        return join;
                    }
                }
                else
                {
                    return Curve.CreateInterpolatedCurve(sortedBoundaryPolyline, 1);
                }
            }
            else
            {
                return null;
            }
        }

        private void GenerateBuilding(Curve[] curves,double w, out List<Polyline> outContour, out List<Polyline> outHoles)
        {
            if (curves.Length == 0)
            {
                outContour = new List<Polyline>();
                outHoles = new List<Polyline>();
                return;
            }
            List<Polyline> polylines = Polyline3D.ConvertCurvesToPolyline(curves).ToList();
            Plane pln = polylines.First().FitPlane();

            //List<Polyline3D.ClosedFilletType> closedType = new List<Polyline3D.ClosedFilletType>() { Polyline3D.ClosedFilletType.Miter };
            //List<Polyline3D.OpenFilletType> openType = new List<Polyline3D.OpenFilletType>() { Polyline3D.OpenFilletType.Butt };

            double tolerance = GH_Component.DocumentTolerance();
            //List<Polyline> outContourLoL;
            //List<Polyline> outHolesLoL;
            //Polyline3D.Offset(polylines, openType, closedType, pln, tolerance, new List<double> { 0.5 * w }, 2, 0.25, out outContourLoL, out outHolesLoL);
            Polyline3D.Offset(polylines, Polyline3D.OpenFilletType.Butt, Polyline3D.ClosedFilletType.Miter, 0.5 * w, pln, tolerance, out outContour, out outHoles);

            //outContour = outContourLoL.First();
            //outHoles = outHolesLoL.First();
            return;
        }

        private Curve TrimBothEnd(Curve boundary, Vector3d depthDirection, Curve needToTrim, int type, double w, double W)
        {
            double tolerance = GH_Component.DocumentTolerance();

            BoundingBox box = boundary.GetBoundingBox(Plane.WorldXY);
            double boxLength = box.Max.X - box.Min.X;

            Plane localCooridinate = new Plane(new Point3d(0, 0, 0), needToTrim.TangentAtStart, depthDirection);
            Plane transformedCoordinate = Plane.WorldXY;
            Transform transform = Transform.PlaneToPlane(transformedCoordinate, localCooridinate);

            Curve current = needToTrim;
            Curve new0;
            Curve new1;
            Curve newCurrent;
            Curve trimedCurve;
            if (type == 1)
            {
                Curve offset0 = current.Offset(Plane.WorldXY, -w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, -w, tolerance, CurveOffsetCornerStyle.None)[0];

                new0 = current.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                newCurrent = offset0.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                new1 = offset1.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                newCurrent.Domain = new Interval(0, 1);
            }
            else if (type == 0)
            {
                Curve offset0 = current.Offset(Plane.WorldXY, -W, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, -(w * 0.5 + W), tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset2 = current.Offset(Plane.WorldXY, -(w + W), tolerance, CurveOffsetCornerStyle.None)[0];

                new0 = offset0.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                newCurrent = offset1.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                new1 = offset2.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                newCurrent.Domain = new Interval(0, 1);
            }
            else
            {
                Curve offset0 = current.Offset(Plane.WorldXY, -w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];

                new0 = offset0.Extend(CurveEnd.Both, boxLength, CurveExtensionStyle.Line);
                newCurrent = current.Extend(CurveEnd.Both, boxLength, CurveExtensionStyle.Line);
                new1 = offset1.Extend(CurveEnd.Both, boxLength, CurveExtensionStyle.Line);
                newCurrent.Domain = new Interval(0, 1);
            }
            

            CurveIntersections curveIntersections0 = Intersection.CurveCurve(new0, boundary, tolerance, tolerance);
            CurveIntersections curveIntersections1 = Intersection.CurveCurve(new1, boundary, tolerance, tolerance);

            Point3d point00;
            Point3d point01;
            Point3d point10;
            Point3d point11;
            if (curveIntersections0.Count == 2 && curveIntersections1.Count == 2)
            {
                // 在对于边缘baseline进行处理时，只有 isGenerateable == false 的时候，是这种情况
                point00 = curveIntersections0[0].PointA;
                point01 = curveIntersections0[1].PointA;
                point10 = curveIntersections1[0].PointA;
                point11 = curveIntersections1[1].PointA;
            }
            else if (curveIntersections0.Count == 1 && curveIntersections0[0].IsOverlap && curveIntersections1.Count == 2)
            {
                // todo：curveIntersections0.Count == 1时，即overlap时要怎么处理
                // 在对于边缘baseline进行处理时，只有 isGenerateable == true 的时候，是这种情况
                if (curveIntersections0[0].IsOverlap)
                {
                    point00 = curveIntersections0[0].PointA;
                    point01 = curveIntersections0[0].PointA2;
                    point10 = curveIntersections1[0].PointA;
                    point11 = curveIntersections1[1].PointA;
                }
                else
                {
                    return newCurrent;
                }
                
            }
            else if (curveIntersections1.Count == 1 && curveIntersections1[0].IsOverlap && curveIntersections0.Count == 2)
            {
                if (curveIntersections1[0].IsOverlap)
                {
                    point00 = curveIntersections0[0].PointA;
                    point01 = curveIntersections0[1].PointA;
                    point10 = curveIntersections1[0].PointA;
                    point11 = curveIntersections1[0].PointA2;
                }
                else
                {
                    return newCurrent;
                }
            }
            else
            {
                return newCurrent;
            }

            Point3d point00Transformed = new Point3d(point00);
            Point3d point01Transformed = new Point3d(point01);
            Point3d point10Transformed = new Point3d(point10);
            Point3d point11Transformed = new Point3d(point11);

            point00Transformed.Transform(transform);
            point01Transformed.Transform(transform);
            point10Transformed.Transform(transform);
            point11Transformed.Transform(transform);

            Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
            lineTransformed.Transform(transform);

            // 投影后，求t值，并比较
            double t00 = lineTransformed.ClosestParameter(point00Transformed);
            double t01 = lineTransformed.ClosestParameter(point01Transformed);
            double t10 = lineTransformed.ClosestParameter(point10Transformed);
            double t11 = lineTransformed.ClosestParameter(point11Transformed);

            double tStart;
            double tEnd;
            if (t00 >= t10)
            {
                //newCurrent.ClosestPoint(point00, out tStart);
                tStart = t00;
            }
            else
            {
                //newCurrent.ClosestPoint(point10, out tStart);
                tStart = t10;
            }

            if (t01 <= t11)
            {
                //newCurrent.ClosestPoint(point01, out tEnd);
                tEnd = t01;
            }
            else
            {
                //newCurrent.ClosestPoint(point11, out tEnd);
                tEnd = t11;
            }

            Interval interval = new Interval(tStart, tEnd);
            trimedCurve = newCurrent.Trim(interval);
            if (trimedCurve == null)
            {
                return newCurrent;
            }

            return trimedCurve;
        }

        private Curve TrimBothEnd(Curve boundary, Curve needToTrim, int type,double w)
        {
            double tolerance = GH_Component.DocumentTolerance();

            bool isOnlyOneSegment = false;
            Point3d pointOnStart = Point3d.Unset;
            Point3d pointOnEnd = Point3d.Unset;
            Curve[] segments;
            if (needToTrim.SpanCount > 1)
            {
                segments = needToTrim.DuplicateSegments();
                isOnlyOneSegment = false;
            }
            else
            {
                segments = new Curve[1] { needToTrim };
                isOnlyOneSegment = true;
            }

            Curve[] resultSegments;
            if (needToTrim.SpanCount > 1)
            {
                resultSegments = needToTrim.DuplicateSegments();
            }
            else
            {
                resultSegments = new Curve[1] { needToTrim };
            }

            Curve new0;
            Curve new1;
            Curve newCurrent;
            Curve trimedSegment;

            CurveIntersections curveIntersections0;
            CurveIntersections curveIntersections1;

            if (type == 0 || type == 2)
            {
                Vector3d vec = segments.First().TangentAtStart;
                Vector3d verticalVec = new Vector3d(-vec.Y, vec.X, vec.Z);

                Plane localCooridinate = new Plane(new Point3d(0, 0, 0), vec, verticalVec);
                Plane transformedCoordinate = Plane.WorldXY;
                Transform transform = Transform.PlaneToPlane(transformedCoordinate, localCooridinate);

                Curve current = segments.First();
                //Curve trimedCurve;
                Curve offset0 = current.Offset(Plane.WorldXY, -w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                new0 = offset0.Extend(CurveEnd.Start, w, CurveExtensionStyle.Line);
                newCurrent = current.Extend(CurveEnd.Start, w, CurveExtensionStyle.Line);
                new1 = offset1.Extend(CurveEnd.Start, w, CurveExtensionStyle.Line);
                new0.Domain = new Interval(0, 1);
                newCurrent.Domain = new Interval(0, 1);
                new1.Domain = new Interval(0, 1);

                curveIntersections0 = Intersection.CurveCurve(new0, boundary, tolerance, tolerance);
                curveIntersections1 = Intersection.CurveCurve(new1, boundary, tolerance, tolerance);

                if (curveIntersections0.Count == 0 && curveIntersections1.Count == 0)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnStart = current.PointAtStart;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if (curveIntersections0.Count == 0 && curveIntersections1.Count == 1 && curveIntersections1[0].IsOverlap)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnStart = current.PointAtStart;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if (curveIntersections0.Count == 1 && curveIntersections1.Count == 0 && curveIntersections0[0].IsOverlap)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnStart = current.PointAtStart;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if(curveIntersections0.Count == 0 && curveIntersections1.Count != 0)
                {
                    IntersectionEvent event1;
                    if (curveIntersections1.Count > 1)
                    {
                        int minIndex = 0;
                        double min = curveIntersections1[0].ParameterA - 0;
                        for (int i = 0; i < curveIntersections1.Count; i++)
                        {
                            double pA = curveIntersections1[i].ParameterA;
                            double delta = pA - 0;

                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event1 = curveIntersections1[minIndex];
                    }
                    else
                    {
                        // curveIntersections1.Count == 1
                        event1 = curveIntersections1[0];
                    }

                    if (!event1.IsOverlap)
                    {
                        if (event1.ParameterA > 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnStart = current.PointAtStart;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else
                        {
                            Point3d point;
                            if (!event1.IsOverlap)
                            {
                                point = event1.PointA;
                            }
                            else
                            {
                                point = event1.PointA;
                            }
                            Point3d pointTransformed = new Point3d(point);
                            pointTransformed.Transform(transform);

                            Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                            lineTransformed.Transform(transform);

                            double t = lineTransformed.ClosestParameter(pointTransformed);

                            Interval interval = new Interval(t, 1);
                            trimedSegment = newCurrent.Trim(interval);
                            if (trimedSegment != null)
                            {
                                resultSegments[0] = trimedSegment;

                                if (isOnlyOneSegment)
                                {
                                    pointOnStart = trimedSegment.PointAtStart;
                                }
                            }
                            else
                            {
                                resultSegments[0] = newCurrent;

                                if (isOnlyOneSegment)
                                {
                                    pointOnStart = newCurrent.PointAtStart;
                                }
                            }
                        }
                    }
                    else
                    {
                        Point3d point;
                        if (!event1.IsOverlap)
                        {
                            point = event1.PointA;
                        }
                        else
                        {
                            point = event1.PointA;
                        }
                        Point3d pointTransformed = new Point3d(point);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(t, 1);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[0] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                        else
                        {
                            resultSegments[0] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = newCurrent.PointAtStart;
                            }
                        }
                    }
                    
                }
                else if (curveIntersections0.Count != 0 && curveIntersections1.Count == 0)
                {
                    IntersectionEvent event0;
                    if (curveIntersections0.Count > 1)
                    {
                        int minIndex = 0;
                        double min = curveIntersections0[0].ParameterA - 0;
                        for (int i = 0; i < curveIntersections0.Count; i++)
                        {
                            double pA = curveIntersections0[i].ParameterA;
                            double delta = pA - 0;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event0 = curveIntersections0[minIndex];
                    }
                    else
                    {
                        event0 = curveIntersections0[0];
                    }

                    if (!event0.IsOverlap)
                    {
                        if (event0.ParameterA > 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnStart = current.PointAtStart;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else
                        {
                            Point3d point;
                            if (!event0.IsOverlap)
                            {
                                point = event0.PointA;
                            }
                            else
                            {
                                point = event0.PointA;
                            }
                            Point3d pointTransformed = new Point3d(point);
                            pointTransformed.Transform(transform);

                            Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                            lineTransformed.Transform(transform);

                            double t = lineTransformed.ClosestParameter(pointTransformed);

                            Interval interval = new Interval(t, 1);
                            trimedSegment = newCurrent.Trim(interval);
                            if (trimedSegment != null)
                            {
                                resultSegments[0] = trimedSegment;

                                if (isOnlyOneSegment)
                                {
                                    pointOnStart = trimedSegment.PointAtStart;
                                }
                            }
                            else
                            {
                                resultSegments[0] = newCurrent;

                                if (isOnlyOneSegment)
                                {
                                    pointOnStart = newCurrent.PointAtStart;
                                }
                            }
                        }
                    }
                    else
                    {
                        Point3d point;
                        if (!event0.IsOverlap)
                        {
                            point = event0.PointA;
                        }
                        else
                        {
                            point = event0.PointA;
                        }
                        Point3d pointTransformed = new Point3d(point);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(t, 1);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[0] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                        else
                        {
                            resultSegments[0] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = newCurrent.PointAtStart;
                            }
                        }
                    }
                    
                }
                else
                {
                    IntersectionEvent event0;
                    IntersectionEvent event1;
                    if (curveIntersections0.Count > 1)
                    {
                        int minIndex = 0;
                        double min = curveIntersections0[0].ParameterA - 0;

                        for (int i = 0; i < curveIntersections0.Count; i++)
                        {
                            double pA = curveIntersections0[i].ParameterA;
                            double delta = pA - 0;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event0 = curveIntersections0[minIndex];
                    }
                    else
                    {
                        //curveIntersections0 == 1
                        event0 = curveIntersections0[0];
                    }
                    if (curveIntersections1.Count > 1)
                    {
                        int minIndex = 0;
                        double min = curveIntersections1[0].ParameterA - 0;

                        for (int i = 0; i < curveIntersections1.Count; i++)
                        {
                            double pA = curveIntersections1[i].ParameterA;
                            double delta = pA - 0;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event1 = curveIntersections1[minIndex];
                    }
                    else
                    {
                        event1 = curveIntersections1[0];
                    }

                    Point3d point00 = Point3d.Unset;
                    Point3d point10 = Point3d.Unset;
                    if (!event0.IsOverlap && !event1.IsOverlap)
                    {
                        // 即不贴边的状态
                        if (event0.ParameterA > 0.5 && event1.ParameterA > 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnStart = current.PointAtStart;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else if (event0.ParameterA > 0.5 && event1.ParameterA <= 0.5)
                        {
                            point10 = event1.PointA;
                        }
                        else if (event0.ParameterA <= 0.5 && event1.ParameterA > 0.5)
                        {
                            point00 = event0.PointA;
                        }
                        else
                        {
                            point00 = event0.PointA;
                            point10 = event1.PointA;
                        }
                    }
                    else if (event0.IsOverlap && !event1.IsOverlap)
                    {
                        // new0贴边
                        point00 = event0.PointA;
                        if (event1.ParameterA > 0.5)
                        {
                            point10 = Point3d.Unset;
                        }
                        else
                        {
                            point10 = event1.PointA;
                        }
                    }
                    else if (!event0.IsOverlap && event1.IsOverlap)
                    {
                        // new1贴边
                        //point00 = event0.PointA;
                        point10 = event1.PointA;
                        if (event0.ParameterA > 0.5)
                        {
                            point00 = Point3d.Unset;
                        }
                        else
                        {
                            point00 = event0.PointA;
                        }
                    }
                    else
                    {
                        point00 = event0.PointA;
                        point10 = event1.PointA;
                    }

                    if (point00 != Point3d.Unset && point10 != Point3d.Unset)
                    {
                        Point3d point00Transformed = new Point3d(point00);
                        Point3d point10Transformed = new Point3d(point10);
                        point00Transformed.Transform(transform);
                        point10Transformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        // 投影后，求t值，并比较
                        double t00 = lineTransformed.ClosestParameter(point00Transformed);
                        double t10 = lineTransformed.ClosestParameter(point10Transformed);

                        double tStart;

                        if (t00 >= t10)
                        {
                            //newCurrent.ClosestPoint(point00, out tStart);
                            tStart = t00;
                        }
                        else
                        {
                            //newCurrent.ClosestPoint(point10, out tStart);
                            tStart = t10;
                        }

                        Interval interval = new Interval(tStart, 1);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[0] = trimedSegment;
                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                        else
                        {
                            resultSegments[0] = newCurrent;
                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                    }
                    else if (point00 == Point3d.Unset && point10 != Point3d.Unset)
                    {
                        Point3d pointTransformed = new Point3d(point10);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(t, 1);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[0] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                        else
                        {
                            resultSegments[0] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = newCurrent.PointAtStart;
                            }
                        }
                    }
                    else if (point00 != Point3d.Unset && point10 == Point3d.Unset)
                    {
                        Point3d pointTransformed = new Point3d(point00);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(t, 1);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[0] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = trimedSegment.PointAtStart;
                            }
                        }
                        else
                        {
                            resultSegments[0] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnStart = newCurrent.PointAtStart;
                            }
                        }
                    }
                    else
                    {
                        if (isOnlyOneSegment)
                        {
                            pointOnStart = current.PointAtStart;
                        }
                        else
                        {
                            // 不对resultSegment做任何处理
                        }
                    }
                    
                }
            }

            if (type == 0 || type == 1)
            {
                Vector3d vec = segments.Last().TangentAtStart;
                Vector3d verticalVec = new Vector3d(-vec.Y, vec.X, vec.Z);

                Plane localCooridinate = new Plane(new Point3d(0, 0, 0), vec, verticalVec);
                Plane transformedCoordinate = Plane.WorldXY;
                Transform transform = Transform.PlaneToPlane(transformedCoordinate, localCooridinate);

                Curve current = segments.Last();
                Curve offset0 = current.Offset(Plane.WorldXY, -w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                new0 = offset0.Extend(CurveEnd.End, w, CurveExtensionStyle.Line);
                newCurrent = current.Extend(CurveEnd.End, w, CurveExtensionStyle.Line);
                new1 = offset1.Extend(CurveEnd.End, w, CurveExtensionStyle.Line);
                new0.Domain = new Interval(0, 1);
                newCurrent.Domain = new Interval(0, 1);
                new1.Domain = new Interval(0, 1);

                curveIntersections0 = Intersection.CurveCurve(new0, boundary, tolerance, tolerance);
                curveIntersections1 = Intersection.CurveCurve(new1, boundary, tolerance, tolerance);

                if (curveIntersections0.Count == 0 && curveIntersections1.Count == 0)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnEnd = current.PointAtEnd;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if (curveIntersections0.Count == 0 && curveIntersections1.Count == 1 && curveIntersections1[0].IsOverlap)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnEnd = current.PointAtEnd;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if (curveIntersections0.Count == 1 && curveIntersections1.Count == 0 && curveIntersections0[0].IsOverlap)
                {
                    if (isOnlyOneSegment)
                    {
                        pointOnEnd = current.PointAtEnd;
                    }
                    else
                    {
                        // 不对resultSegment做任何处理
                    }
                }
                else if (curveIntersections0.Count == 0 && curveIntersections1.Count != 0)
                {
                    IntersectionEvent event1;
                    if (curveIntersections1.Count > 1)
                    {
                        int minIndex = 0;
                        double min = 1 - curveIntersections1[0].ParameterA;
                        for (int i = 0; i < curveIntersections1.Count; i++)
                        {
                            double pA = curveIntersections1[i].ParameterA;
                            double delta = 1 - pA;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event1 = curveIntersections1[minIndex];
                    }
                    else
                    {
                        // curveIntersections1.Count == 1
                        event1 = curveIntersections1[0];
                    }

                    if (!event1.IsOverlap)
                    {
                        if (event1.ParameterA < 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = current.PointAtEnd;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else
                        {
                            Point3d point;
                            if (!event1.IsOverlap)
                            {
                                point = event1.PointA2;
                            }
                            else
                            {
                                point = event1.PointA;
                            }

                            Point3d pointTransformed = new Point3d(point);
                            pointTransformed.Transform(transform);

                            Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                            lineTransformed.Transform(transform);

                            double t = lineTransformed.ClosestParameter(pointTransformed);

                            Interval interval = new Interval(0, t);
                            trimedSegment = newCurrent.Trim(interval);
                            if (trimedSegment != null)
                            {
                                resultSegments[resultSegments.Length - 1] = trimedSegment;

                                if (isOnlyOneSegment)
                                {
                                    pointOnEnd = trimedSegment.PointAtEnd;
                                }
                            }
                            else
                            {
                                resultSegments[resultSegments.Length - 1] = newCurrent;

                                if (isOnlyOneSegment)
                                {
                                    pointOnEnd = newCurrent.PointAtEnd;
                                }
                            }
                        }
                    }
                    else
                    {
                        Point3d point;
                        if (!event1.IsOverlap)
                        {
                            point = event1.PointA2;
                        }
                        else
                        {
                            point = event1.PointA;
                        }

                        Point3d pointTransformed = new Point3d(point);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(0, t);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[resultSegments.Length - 1] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = trimedSegment.PointAtEnd;
                            }
                        }
                        else
                        {
                            resultSegments[resultSegments.Length - 1] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = newCurrent.PointAtEnd;
                            }
                        }
                    }
                }
                else if (curveIntersections0.Count != 0 && curveIntersections1.Count == 0)
                {
                    IntersectionEvent event0;
                    if (curveIntersections0.Count > 1)
                    {
                        int minIndex = 0;
                        double min = 1 - curveIntersections0[0].ParameterA;
                        for (int i = 0; i < curveIntersections0.Count; i++)
                        {
                            double pA = curveIntersections0[i].ParameterA;
                            double delta = 1 - pA;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }
                        event0 = curveIntersections0[minIndex];
                    }
                    else
                    {
                        // curveIntersections0.Count == 1
                        event0 = curveIntersections0[0];
                    }

                    if (!event0.IsOverlap)
                    {
                        if (event0.ParameterA < 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = current.PointAtEnd;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else
                        {
                            Point3d point;
                            if (!event0.IsOverlap)
                            {
                                point = event0.PointA2;
                            }
                            else
                            {
                                point = event0.PointA;
                            }

                            Point3d pointTransformed = new Point3d(point);
                            pointTransformed.Transform(transform);

                            Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                            lineTransformed.Transform(transform);

                            double t = lineTransformed.ClosestParameter(pointTransformed);

                            Interval interval = new Interval(0, t);
                            trimedSegment = newCurrent.Trim(interval);
                            if (trimedSegment != null)
                            {
                                resultSegments[resultSegments.Length - 1] = trimedSegment;

                                if (isOnlyOneSegment)
                                {
                                    pointOnEnd = trimedSegment.PointAtEnd;
                                }
                            }
                            else
                            {
                                resultSegments[resultSegments.Length - 1] = newCurrent;

                                if (isOnlyOneSegment)
                                {
                                    pointOnEnd = newCurrent.PointAtEnd;
                                }
                            }
                        }
                    }
                    else
                    {
                        Point3d point;
                        if (!event0.IsOverlap)
                        {
                            point = event0.PointA2;
                        }
                        else
                        {
                            point = event0.PointA;
                        }

                        Point3d pointTransformed = new Point3d(point);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(0, t);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[resultSegments.Length - 1] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = trimedSegment.PointAtEnd;
                            }
                        }
                        else
                        {
                            resultSegments[resultSegments.Length - 1] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = newCurrent.PointAtEnd;
                            }
                        }
                    }
                }
                else
                {
                    //todo
                    IntersectionEvent event0;
                    IntersectionEvent event1;
                    if (curveIntersections0.Count > 1)
                    {
                        int minIndex = 0;
                        double min = 1 - curveIntersections0[0].ParameterA;

                        for (int i = 0; i < curveIntersections0.Count; i++)
                        {
                            double pA = curveIntersections0[i].ParameterA;
                            double delta = 1 - pA;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }

                        event0 = curveIntersections0[minIndex];
                    }
                    else
                    {
                        event0 = curveIntersections0[0];
                    }
                    if (curveIntersections1.Count > 1)
                    {
                        int minIndex = 0;
                        double min = 1 - curveIntersections1[0].ParameterA;

                        for (int i = 0; i < curveIntersections1.Count; i++)
                        {
                            double pA = curveIntersections1[i].ParameterA;
                            double delta = 1 - pA;
                            if (delta < 0)
                            {
                                continue;
                            }
                            if (delta < min)
                            {
                                min = delta;
                                minIndex = i;
                            }
                        }

                        event1 = curveIntersections1[minIndex];
                    }
                    else
                    {
                        event1 = curveIntersections1[0];
                    }

                    Point3d point01 = Point3d.Unset;
                    Point3d point11 = Point3d.Unset;

                    if (!event0.IsOverlap && !event1.IsOverlap)
                    {
                        // 即不贴边的状态
                        if (event0.ParameterA < 0.5 && event1.ParameterA < 0.5)
                        {
                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = current.PointAtEnd;
                            }
                            else
                            {
                                // 不对resultSegment做任何处理
                            }
                        }
                        else if (event0.ParameterA < 0.5 && event1.ParameterA >= 0.5)
                        {
                            point11 = event1.PointA;
                        }
                        else if (event0.ParameterA >= 0.5 && event1.ParameterA < 0.5)
                        {
                            point01 = event0.PointA;
                        }
                        else
                        {
                            point01 = event0.PointA;
                            point11 = event1.PointA;
                        }
                    }
                    else if (event0.IsOverlap && !event1.IsOverlap)
                    {
                        // new0贴边
                        point01 = event0.PointA2;
                        if (event1.ParameterA < 0.5)
                        {
                            point11 = Point3d.Unset;
                        }
                        else
                        {
                            point11 = event1.PointA;
                        }
                    }
                    else if (!event0.IsOverlap && event1.IsOverlap)
                    {
                        point11 = event1.PointA2;
                        if (event0.ParameterA < 0.5)
                        {
                            point01 = Point3d.Unset;
                        }
                        else
                        {
                            point01 = event0.PointA;
                        }
                    }
                    else
                    {
                        point01 = event0.PointA;
                        point11 = event1.PointA;
                    }


                    if (point01 != Point3d.Unset && point11 != Point3d.Unset)
                    {
                        Point3d point01Transformed = new Point3d(point01);
                        Point3d point11Transformed = new Point3d(point11);

                        point01Transformed.Transform(transform);
                        point11Transformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t01 = lineTransformed.ClosestParameter(point01Transformed);
                        double t11 = lineTransformed.ClosestParameter(point11Transformed);

                        double tEnd;
                        if (t01 <= t11)
                        {
                            //newCurrent.ClosestPoint(point01, out tEnd);
                            tEnd = t01;
                        }
                        else
                        {
                            //newCurrent.ClosestPoint(point11, out tEnd);
                            tEnd = t11;
                        }

                        Interval interval = new Interval(0, tEnd);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[resultSegments.Length - 1] = trimedSegment;
                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = trimedSegment.PointAtEnd;
                            }
                        }
                        else
                        {
                            resultSegments[resultSegments.Length - 1] = newCurrent;
                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = newCurrent.PointAtEnd;
                            }
                        }
                    }
                    else if (point01 == Point3d.Unset && point11 != Point3d.Unset)
                    {
                        Point3d pointTransformed = new Point3d(point11);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(0, t);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[resultSegments.Length - 1] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = trimedSegment.PointAtEnd;
                            }
                        }
                        else
                        {
                            resultSegments[resultSegments.Length - 1] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = newCurrent.PointAtEnd;
                            }
                        }
                    }
                    else if (point01 != Point3d.Unset && point11 == Point3d.Unset)
                    {
                        Point3d pointTransformed = new Point3d(point01);
                        pointTransformed.Transform(transform);

                        Line lineTransformed = new Line(newCurrent.PointAtStart, newCurrent.PointAtEnd);
                        lineTransformed.Transform(transform);

                        double t = lineTransformed.ClosestParameter(pointTransformed);

                        Interval interval = new Interval(0, t);
                        trimedSegment = newCurrent.Trim(interval);
                        if (trimedSegment != null)
                        {
                            resultSegments[resultSegments.Length - 1] = trimedSegment;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = trimedSegment.PointAtEnd;
                            }
                        }
                        else
                        {
                            resultSegments[resultSegments.Length - 1] = newCurrent;

                            if (isOnlyOneSegment)
                            {
                                pointOnEnd = newCurrent.PointAtEnd;
                            }
                        }
                    }
                    else
                    {
                        if (isOnlyOneSegment)
                        {
                            pointOnEnd = current.PointAtEnd;
                        }
                        else
                        {
                            // 不对resultSegment做任何处理
                        }
                    }

                }
            }

            if (isOnlyOneSegment)
            {
                resultSegments[0] = new Line(pointOnStart, pointOnEnd).ToNurbsCurve();
                return resultSegments[0];
            }
            else
            {
                Curve[] crvs = Curve.JoinCurves(resultSegments);
                return crvs[0];
            }
        }


        private Curve TrimByBoundaryWithCurve(Curve crv, Curve boundary)
        {
            Brep plan = Brep.CreatePlanarBreps(boundary, GH_Component.DocumentTolerance())[0];
            BrepFace planFace = plan.Faces[0]; ;

            Curve[] overlapCurves = null;
            Point3d[] intersectionPoints = null;
            Intersection.CurveBrepFace(crv, planFace, GH_Component.DocumentTolerance(), out overlapCurves, out intersectionPoints);

            if (overlapCurves.Length != 0)
            {
                return overlapCurves[0];
            }
            else
            {
                return crv;
            }
        }

        private List<Curve> ConstructTweenCurve(int insertGapCount, 
                                                double lMin,
                                                double lMax, 
                                                double dw,
                                                double w, 
                                                Curve baseLine0, 
                                                Curve baseLine1,
                                                Point3d point0, 
                                                Point3d point1, 
                                                bool isGap,
                                                bool isLMax,
                                                bool NeedRandom)
        {
            if (baseLine0 == null && baseLine1 == null)
            {
                return new List<Curve>();
            }

            Curve baseLine1Reverse = baseLine1.DuplicateCurve();
            baseLine1Reverse.Reverse();

            List<double> factors = new List<double>();
            if (isGap)
            {
                // 如果是进行GapCurve的计算
                double sum = 0.0;
                if (isLMax)
                {
                    sum = 2 * (lMin + w) + dw + (insertGapCount + 1) * (lMax + dw);
                }
                else
                {
                    sum = 2 * (lMin + w) + dw + (insertGapCount + 1) * (lMin + dw);
                }

                //// 在计算完总长度 sum 后，要将输入的 insertGapCount + 2，因为要重新生成新的eastBaseLine和westBaseLine，也作为GapCurve
                // insertGapCount = insertGapCount + 2;

                List<double> eachLengths = new List<double>();
                eachLengths.Add(lMin + w + dw * 0.5);
                for (int i = 0; i < insertGapCount + 1; i++)
                {
                    if (isLMax)
                    {
                        eachLengths.Add(lMax + dw);
                    }
                    else
                    {
                        eachLengths.Add(lMin + dw);
                    }
                }
                eachLengths.Add(lMin + w + dw * 0.5);

                if (NeedRandom)
                {
                    eachLengths = RandomSort(eachLengths);
                }
                //List<double> randomEachLengths = RandomSort(eachLengths);
                List<double> cumulatives = CalCumulative(eachLengths);

                //List<double> factors = new List<double>();
                for (int i = 0; i < insertGapCount + 2; i++)
                {
                    double factor = cumulatives[i] / sum;
                    factors.Add(factor);
                }
            }
            else
            {
                double sum = 0.0;
                if (isLMax)
                {
                    sum = (insertGapCount + 1) * (lMax + dw) - dw;
                }
                else
                {
                    sum = (insertGapCount + 1) * (lMin + dw) - dw;
                }
                List<double> eachLengths = new List<double>();
                eachLengths.Add(lMin + 0.5 * dw);
                for (int i = 1; i < insertGapCount; i++)
                {
                    if (isLMax)
                    {
                        eachLengths.Add(lMax + dw);
                    }
                    else
                    {
                        eachLengths.Add(lMin + dw);
                    }
                }
                eachLengths.Add(lMin + 0.5 * dw);

                if (NeedRandom)
                {
                    eachLengths = RandomSort(eachLengths);
                }
                List<double> cumulatives = CalCumulative(eachLengths);
                for (int i = 0; i < insertGapCount; i++)
                {
                    double factor = cumulatives[i] / sum;
                    factors.Add(factor);
                }
            }
            
            ////double sum = (insertLMinOrLMaxCount * 1) * lMinOrlMax + insertLMinOrLMaxCount * d;
            //double sum = (insertGapCount - 1) * (lMinOrlMax + dw) + dw;
            //List<double> eachLengths = new List<double>();
            //eachLengths.Add(0.5 * dw);
            //for (int i = 1; i < insertGapCount; i++)
            //{
            //    double num = eachLengths[i - 1] + lMinOrlMax + dw;
            //    eachLengths.Add(num);
            //}

            //List<double> factors = new List<double>();
            //for (int i = 0; i < insertGapCount; i++)
            //{
            //    double factor = eachLengths[i] / sum;
            //    factors.Add(factor);
            //}

            Surface surface;
            if (baseLine0 == null && baseLine1 != null)
            {
                // baseLine1.Reverse();
                surface = Surface.CreateExtrusionToPoint(baseLine1Reverse, point0);
                surface = surface.Transpose();
                surface = surface.Reverse(0);
            }
            else if (baseLine0 != null && baseLine1 == null)
            {
                surface = Surface.CreateExtrusionToPoint(baseLine0, point1);
                surface = surface.Transpose();
            }
            else
            {
                // baseLine1.Reverse();
                surface = NurbsSurface.CreateRuledSurface(baseLine0, baseLine1Reverse);
            }

            List<Curve> insectCurves = new List<Curve>();
            for (int i = 0; i < factors.Count; i++)
            {
                double t = surface.Domain(1).ParameterAt(factors[i]);
                Curve curve = surface.IsoCurve(0, t);
                //if (isGap)
                //{
                //    curve.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
                //}
                insectCurves.Add(curve);
            }
            surface.Dispose();

            return insectCurves;
        }

        private double CalSetback(Vector3d a1a2, Vector3d a2b1, Vector3d b1b2, double a2b1Length, double targetLength, out double setbackE)
        {
            double angle0 = Vector3d.VectorAngle(a1a2, a2b1);
            double angle1 = Vector3d.VectorAngle(a2b1, b1b2);
            double a = Rhino.RhinoMath.ToDegrees(angle0);
            double b = Rhino.RhinoMath.ToDegrees(angle1);

            double x = 0;
            if (angle0 < Math.PI * 0.5 && angle1 < Math.PI * 0.5)
            {
                double tan0 = Math.Tan(angle0);
                double tan1 = Math.Tan(angle1);

                x = (targetLength - a2b1Length) * tan0 * tan1 / (tan0 + tan1);
            }
            else if (angle0 < Math.PI * 0.5 && angle1 == Math.PI * 0.5)
            {
                double tan0 = Math.Tan(angle0);

                x = (targetLength - a2b1Length) * tan0;
            }
            else if (angle0 < Math.PI * 0.5 && angle1 > Math.PI * 0.5)
            {
                double tan0 = Math.Tan(angle0);
                double tan1 = Math.Tan(Math.PI - angle1);

                if (tan0 == tan1)
                {
                    x = -1;
                }
                else
                {
                    x = (targetLength - a2b1Length) * tan0 * tan1 / (tan1 - tan0);
                }
            }
            else if (angle0 == Math.PI * 0.5 && angle1 < Math.PI * 0.5)
            {
                double tan1 = Math.Tan(angle1);

                x = (targetLength - a2b1Length) * tan1;
            }
            //else if (angle0 == Math.PI * 0.5 && angle1 == Math.PI * 0.5)
            //{
            //    // x = 任意
            //}
            else if (angle0 == Math.PI * 0.5 && angle1 > Math.PI * 0.5)
            {
                double tan1 = Math.Tan(Math.PI - angle1);
                x = (a2b1Length - targetLength) * tan1;
            }
            else if (angle0 > Math.PI * 0.5 && angle1 < Math.PI * 0.5)
            {
                double tan0 = Math.Tan(Math.PI - angle0);
                double tan1 = Math.Tan(angle1);

                if (tan0 == tan1)
                {
                    x = -1;
                }
                else
                {
                    x = (targetLength - a2b1Length) * tan0 * tan1 / (tan0 - tan1);
                }
            }
            else if (angle0 > Math.PI * 0.5 && angle1 == Math.PI * 0.5)
            {
                double tan0 = Math.Tan(Math.PI - angle0);

                x = (a2b1Length - targetLength) * tan0;
            }
            else if (angle0 > Math.PI * 0.5 && angle1 > Math.PI * 0.5)
            {
                double tan0 = Math.Tan(Math.PI - angle0);
                double tan1 = Math.Tan(Math.PI - angle1);
                x = (a2b1Length - targetLength) * tan0 * tan1 / (tan0 + tan1);
            }
            else
            {
                // x = 任意
                x = -1;
            }

            double setbackW;
            if (angle1 < Math.PI * 0.5)
            {
                setbackW = x / Math.Sin(angle1);
            }
            else if (angle1 > Math.PI * 0.5)
            {
                setbackW = x / Math.Sin(Math.PI  - angle1);
            }
            else
            {
                setbackW = x;
            }

            if (angle0 < Math.PI * 0.5)
            {
                setbackE = x / Math.Sin(angle0);
            }
            else if (angle0 > Math.PI * 0.5)
            {
                setbackE = x / Math.Sin(Math.PI - angle0);
            }
            else
            {
                setbackE = x;
            }

            return setbackW;
        }

        private Curve[] GenerateCenterLineNearGapCurve(bool isBasedGapCurve,Curve gapCurve, double d, double w, double W, Curve boundary)
        {
            Curve offset0;
            Curve offset1;
            if (isBasedGapCurve)
            {
                offset0 = gapCurve.Offset(Plane.WorldXY, -(0.5 * d + 0.5 * w), GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                offset1 = gapCurve.Offset(Plane.WorldXY, 0.5 * d + 0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            }
            else
            {
                offset0 = gapCurve.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
                offset1 = gapCurve.Offset(Plane.WorldXY, 0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            }
            
            offset0 = offset0.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
            offset1 = offset1.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
            Curve trimed0 = TrimByBoundaryWithCurve(offset0, boundary);
            Curve trimed1 = TrimByBoundaryWithCurve(offset1, boundary);
            trimed0.Domain = new Interval(0, 1);
            trimed1.Domain = new Interval(0, 1);
            return new Curve[2] { trimed0, trimed1 };
        }
        
        /// <summary>
        /// 向量AB与向量AP的外积（叉积）
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private double Cross(Point3d a, Point3d b,Point3d p)
        {
            Vector3d ab = new Vector3d(b - a);
            Vector3d ap = new Vector3d(p - a);
            return ab.X * ap.Y - ab.Y * ap.X;
        }

        /// <summary>
        /// 向量AB与向量AP的内积（点积）
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private double Dot(Point3d a, Point3d b, Point3d p)
        {
            Vector3d ab = new Vector3d(b - a);
            Vector3d ap = new Vector3d(p - a);
            return ab.X * ap.X + ab.Y * ap.Y;
        }

        /// <summary>
        /// 点P与线段AB位置关系(判断点P在线段AB的哪个方向上)
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        private int Dir(Point3d a, Point3d b, Point3d p)
        {
            double cross = Cross(a, b, p);
            double dot = Dot(a, b, p);
            if (cross < 0) return -1;       // 逆时针
            else if (cross > 0) return 1;   // 顺时针
            else if (dot < 0) return -2;    // 反延长线
            else
            {
                if (Dis2(a, b) < Dis2(a, p)) return 2;  // 延长线
                return 0;                               // p在线段AB上
            }
        }

        /// <summary>
        /// 点a、b距离的平方
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        private double Dis2(Point3d a,Point3d b)
        {
            return Math.Pow((a.X - b.X), 2) + Math.Pow((a.Y - b.Y), 2);
        }

        /// <summary>
        /// 点b到直线a1a2的距离
        /// </summary>
        /// <param name="a1"></param>
        /// <param name="a2"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        private double Dis2PointToStraightLine(Point3d a1, Point3d a2, Point3d b)
        {
            Vector3d a2a1 = a1 - a2;
            Vector3d a2b = b - a2;
            Vector3d normal = Vector3d.CrossProduct(a2a1, a2b);
            normal.Unitize();

            return Math.Abs(a2b * normal);
        }

        /// <summary>
        /// 点P到线段AB的最短距离
        /// </summary>
        /// <returns></returns>
        private double DisMin(Point3d a, Point3d b, Point3d p)
        {
            double r = ((p.X - a.X) * (b.X - a.X) + (p.Y - a.Y) * (b.Y - a.Y)) / Dis2(a, b);
            if (r <= 0)//第一种情况，P对AB的垂点在AB反向延长线上, 返回AP的长
            {
                return Math.Sqrt(Dis2(a, p));
            }
            else if (r >= 1)//第二种情况, 返回BP的长度
            {
                return Math.Sqrt(Dis2(b, p));
            }
            else//第三种情况, 返回PC的长度
            {
                double AC = r * Math.Sqrt(Dis2(a, b));//先求AC的长度,(AC=r*|AB|)
                return Math.Sqrt(Dis2(a, p) - AC * AC);//再勾股定理返回PC的长度
            }
        }

        /// <summary>
        /// 两个线段之间的最短距离
        /// </summary>
        /// <param name="a1"></param>
        /// <param name="a2"></param>
        /// <param name="b1"></param>
        /// <param name="b2"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        private double MinDistanceBetweenTwoLineSegment(Point3d a1, Point3d a2, Point3d b1, Point3d b2, out int index)
        {
            
            if (Dir(a1, a2, b1)*Dir(a1, a2, b2) <= 0 && Dir(b1, b2, a1)*Dir(b1, b2, a2) <= 0)
            {
                // 考虑相交的情况，两线段相交，距离为0
                index = -1;
                return 0;
            }
            else
            {
                // 不相交，则最短距离为每个端点到另一条线段距离的最小值

                // 为了区分最后返回距离的正负(规定a1a1在b1b2下方时，距离值为正，反之为负)
                bool isA1A2UnderB1B2 = true;
                // 按照a1a2b1b2的顺序，将这四个点连接形成四边形，如果是逆时针的方向，则a1a2一定在b1b2下方，否则就在上方
                List<Point3d> pts = new List<Point3d>() { a1, a2, b1, b2 };
                int num = GetPolygonDirection(pts);
                if (num == 0)
                {
                    isA1A2UnderB1B2 = true;
                }
                else if (num == 1)
                {
                    isA1A2UnderB1B2 = false;
                }
                else
                {
                    // 四点共线
                    index = -1;
                    return 0;
                }

                List<double> distance = new List<double>();
                distance.Add(DisMin(a1, a2, b1));
                distance.Add(DisMin(a1, a2, b2));
                distance.Add(DisMin(b1, b2, a2));
                distance.Add(DisMin(b1, b2, a1));
                
                double min = distance.Min();
                // (a1, a2, b1):右端，与a1a2垂直，即index = 0
                // (a1, a1, b2):左端，与a1a2垂直，即index = 1
                // (b1, b2, a2):右端，与b1b2垂直，即index = 2
                // (b1, b2, a1):左端，与b1b2垂直，即index = 3
                index = distance.IndexOf(min);
                List<double> distanceForSecondMin = new List<double>();
                distanceForSecondMin.AddRange(distance);
                double secondMin = distanceForSecondMin.OrderBy(z => z).Skip(1).First();

                if (min - secondMin > -0.001 && min - secondMin < 0.001)
                {
                    if (index % 2 == 0)
                    {
                        // 右端，都不垂直
                        index = 4;
                    }
                    else
                    {
                        // 左端，都不垂直
                        index = 5;
                    }
                }

                if (isA1A2UnderB1B2)
                {
                    return min;
                }
                else
                {
                    return -min;
                }
            }
        }

        /// <summary>
        /// 判断多边形是顺时针还是逆时针, -1表示可能是不可计算的图形，比如多点共线；0表示逆时针；1表示顺时针
        /// </summary>
        /// <param name="pts"></param>
        /// <returns></returns>
        private int GetPolygonDirection(List<Point3d> pts)
        {
            int j, k;
            int count = 0;
            double z;
            if (pts == null || pts.Count < 3)
            {
                // 无.可能是不可计算的图形，比如多点共线
                return -1;
            }
            int n = pts.Count;
            for (int i = 0; i < n; i++)
            {
                j = (i + 1) % n;
                k = (i + 2) % n;
                z = (pts[j].X - pts[i].X) * (pts[k].Y - pts[j].Y);
                z -= (pts[j].Y - pts[i].Y) * (pts[k].X - pts[j].X);
                if (z < 0)
                {
                    count--;
                }
                else if (z > 0)
                {
                    count++;
                }
            }

            if (count > 0)
            {
                // 逆时针
                return 0;
            }
            else if (count < 0)
            {
                // 顺时针
                return 1;
            }
            else
            {
                // 无.可能是不可计算的图形，比如多点共线
                return -1;
            }
        }

        private List<T> RandomSort<T>(List<T> list)
        {
            //var random = new Random();
            var newList = new List<T>();
            foreach (var item in list)
            {
                newList.Insert(m_random.Next(newList.Count), item);
            }
            return newList;
        }

        private List<double> CalCumulative(List<double> list)
        {
            List<double> cumulatives = new List<double>();
            cumulatives.Add(list[0]);
            for (int i = 1; i < list.Count; i++)
            {
                cumulatives.Add(cumulatives[i - 1] + list[i]);
            }

            return cumulatives;
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            //base.DrawViewportMeshes(args);
            for (int i = 0; i < InnerResultPolyline.Count; i++)
            {
                args.Display.DrawPolyline(InnerResultPolyline[i], Color.DarkRed, Thickness);
            }
            for (int i = 0; i < InnerNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            // Temp
            //for (int i = 0; i < SCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(SCrvForShowList[i], Color.Red, 2);
            //}
            //for (int i = 0; i < NCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(NCrvForShowList[i], Color.Green, 2);
            //}
            //for (int i = 0; i < ECrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(ECrvForShowList[i], Color.Blue, 2);
            //}
            //for (int i = 0; i < WCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(WCrvForShowList[i], Color.Yellow, 2);
            //}
            //for (int i = 0; i < GapCenterLineList.Count; i++)
            //{
            //    args.Display.DrawCurve(GapCenterLineList[i], Color.Black, 3);
            //}

            //for (int i = 0; i < SBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(SBCrvForShowList[i], Color.Red, 2);
            //}
            //for (int i = 0; i < NBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(NBCrvForShowList[i], Color.Green, 2);
            //}
            //for (int i = 0; i < EBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(EBCrvForShowList[i], Color.Blue, 2);
            //}
            //for (int i = 0; i < WBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(WBCrvForShowList[i], Color.Yellow, 2);
            //}

            for (int i = 0; i < AllBlockCellBoundaryDT.BranchCount; i++)
            {
                for (int j = 0; j < AllBlockCellBoundaryDT.Branch(i).Count; j++)
                {
                    args.Display.DrawCurve(AllBlockCellBoundaryDT.Branch(i)[j], Color.DeepPink, 4);
                }
            }
            //for (int i = 0; i < AllBlockCellBoundaryDT2.BranchCount; i++)
            //{
            //    for (int j = 0; j < AllBlockCellBoundaryDT2.Branch(i).Count; j++)
            //    {
            //        args.Display.DrawCurve(AllBlockCellBoundaryDT2.Branch(i)[j], Color.DeepSkyBlue, 4);
            //    }
            //}
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            //base.DrawViewportMeshes(args);
            for (int i = 0; i < InnerResultPolyline.Count; i++)
            {
                args.Display.DrawPolyline(InnerResultPolyline[i], Color.DarkRed, Thickness);
            }
            for (int i = 0; i < InnerNodeTextDots.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerNodeTextDots[i], Color.ForestGreen, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            // Temp
            //for (int i = 0; i < SCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(SCrvForShowList[i], Color.Red, 2);
            //}
            //for (int i = 0; i < NCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(NCrvForShowList[i], Color.Green, 2);
            //}
            //for (int i = 0; i < ECrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(ECrvForShowList[i], Color.Blue, 2);
            //}
            //for (int i = 0; i < WCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(WCrvForShowList[i], Color.Yellow, 2);
            //}
            //for (int i = 0; i < GapCenterLineList.Count; i++)
            //{
            //    args.Display.DrawCurve(GapCenterLineList[i], Color.Black, 3);
            //}

            //for (int i = 0; i < SBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(SBCrvForShowList[i], Color.Red, 2);
            //}
            //for (int i = 0; i < NBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(NBCrvForShowList[i], Color.Green, 2);
            //}
            //for (int i = 0; i < EBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(EBCrvForShowList[i], Color.Blue, 2);
            //}
            //for (int i = 0; i < WBCrvForShowList.Count; i++)
            //{
            //    args.Display.DrawCurve(WBCrvForShowList[i], Color.Yellow, 2);
            //}

            for (int i = 0; i < AllBlockCellBoundaryDT.BranchCount; i++)
            {
                for (int j = 0; j < AllBlockCellBoundaryDT.Branch(i).Count; j++)
                {
                    args.Display.DrawCurve(AllBlockCellBoundaryDT.Branch(i)[j], Color.DeepPink, 4);
                }
            }

            //for (int i = 0; i < AllBlockCellBoundaryDT2.BranchCount; i++)
            //{
            //    for (int j = 0; j < AllBlockCellBoundaryDT2.Branch(i).Count; j++)
            //    {
            //        args.Display.DrawCurve(AllBlockCellBoundaryDT2.Branch(i)[j], Color.DeepSkyBlue, 4);
            //    }
            //}
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
            get { return new Guid("4b850a3f-93e2-44a8-b44a-79f109cd61c3"); }
        }

        public override void CreateAttributes()/* 重写CreateAttribute方法以启用自定义电池外观 */
        {
            Attributes = new ShowAttribute(this);
        }

        public class ShowAttribute: GH_ComponentAttributes
        {
            public ShowAttribute(GhcVolumeBaseOnTopo component) : base(component) { }

            protected override void Layout()
            {
                base.Layout();
                /* 先执行base.Layout()，可以按GH电池默认方式计算电池的出/入口需要的高度，我们在下面基于这个高度进行更改 */
                Bounds = new RectangleF(Bounds.X, Bounds.Y, Bounds.Width, Bounds.Height + 42.0f);
            }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);/* 执行基本的电池渲染 */

                /* 额外的电池渲染，仅在“Objects”这个渲染轨道绘制 */
                if (channel == GH_CanvasChannel.Objects)
                {
                    RectangleF buttonRect1 = /* 按钮的位置 */ new RectangleF(Bounds.X, Bounds.Bottom - 42, Bounds.Width, 20.0f);

                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect1.Inflate(-2.0f, -2.0f);

                    using (GH_Capsule capsule1 = GH_Capsule.CreateCapsule(buttonRect1, GH_Palette.Normal))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule1.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("Horizontal:{0}", ((GhcVolumeBaseOnTopo)Owner).HBranchType.ToString()),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.Black,
                                        buttonRect1,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });

                    RectangleF buttonRect2 = new RectangleF(Bounds.X, Bounds.Bottom - 20, Bounds.Width, 20.0f);
                    /* 在X、Y方向分别留出2px的空隙，以免button贴住电池边 */
                    buttonRect2.Inflate(-2.0f, -2.0f);

                    using (GH_Capsule capsule2 = GH_Capsule.CreateCapsule(buttonRect2, GH_Palette.Black))
                    {
                        /* 按照该电池的“是否被选中”、“是否被锁定”、“是否隐藏”三个属性来决定渲染的按钮样式 */
                        /* 这样可以使得我们的按钮更加贴合GH原生的样式 */
                        /* 也可以自己换用其他的capsule.Render()重载，渲染不同样式电池 */
                        capsule2.Render(graphics, Selected, Owner.Locked, Owner.Hidden);
                    }

                    graphics.DrawString(string.Format("Vertical:{0}", ((GhcVolumeBaseOnTopo)Owner).VBranchType.ToString()),
                                        new Font(GH_FontServer.ConsoleSmall, FontStyle.Bold),
                                        Brushes.White,
                                        buttonRect2,
                                        new StringFormat()
                                        {
                                            Alignment = StringAlignment.Center,
                                            LineAlignment = StringAlignment.Center
                                        });
                }
            }
        }
    }
}