using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
//using Grasshopper.Kernel.Geometry;
//using Grasshopper.Kernel.Geometry.ConvexHull;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Linq;
using System.Drawing;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

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

            pManager.AddIntegerParameter("Seed", "", "", GH_ParamAccess.item);
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

            pManager.AddCurveParameter("H", "", "", GH_ParamAccess.tree);
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

            int seed = 0;

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetData<PlanktonMesh>("newDualPlanktonMesh", ref dual)
                && DA.GetDataList("SetbackList", setbackList)
                && DA.GetData("w", ref w)
                && DA.GetData("W", ref W)
                && DA.GetData("lMin", ref lMin)
                && DA.GetData("lMax", ref lMax)
                && DA.GetData("d", ref d)
                && DA.GetData<int>("Seed", ref seed))
            {
                DA.GetDataTree<GH_Integer>("WhichPairSetbackNeedToChange", out pairNeedToChange);
                DA.GetDataTree<GH_Number>("SetbackValueNeedToChange", out setbackValueNeedToChange);

                // 固定种子随机数生成器
                m_random = new Random(seed);

                Thickness = 2;

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
                //List<List<double>> allFaceBSSetBack = new List<List<double>>();
                //for (int i = 0; i < allFaceBS.Count; i++)
                //{
                //    allFaceBSSetBack.Add(new List<double>());
                //}

                for (int i = 0; i < innerNodeCount; i++)
                {
                    int dFIndex = dFIndex_PVIndex.IndexOf(i);

                    //allFaceBSSetBack.Add(new List<double>());

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
                //DataTree<Curve> cuts;
                //DataTree<Curve> publicSide;
                //DataTree<Brep> cutsBrepDT;
                //List<Brep> resultBreps;
                List<List<BoundarySegment>> newAllFaceBS = OffsetBS(allFaceBS,
                                                                    facePolylines,
                                                                    out innerResultPolylines); 
                                                                    //out cuts, 
                                                                    //out publicSide,
                                                                    //out cutsBrepDT,
                                                                    //out resultBreps);
                #endregion
                DA.SetDataList("SetBack Polyline", innerResultPolylines);
                //DA.SetDataTree(1, cuts);
                //DA.SetDataTree(2, publicSide);
                //DA.SetDataTree(3, cutsBrepDT);
                //DA.SetDataList("ResultBreps", resultBreps);

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



                //for (int i = 0; i < allFaceBS.Count; i++)
                //{
                //    Curve joinedCurve = GetJoinedCurveFromBS(allFaceBS[i]);

                //    CalEachBSDomain(allFaceBS[i], joinedCurve);
                //}

                List<List<Curve>> curveUnionLoL = new List<List<Curve>>();

                SCrvForShowList.Clear();
                NCrvForShowList.Clear();
                ECrvForShowList.Clear();
                WCrvForShowList.Clear();
                GapCenterLineList.Clear();

                for (int i = 0; i < newAllFaceBS.Count; i++)
                {
                    curveUnionLoL.Add(new List<Curve>());

                    bool isHorizontalLayout;
                    Polyline sortedBoundaryPolyline;
                    List<Line> lineForEachEdge;
                    List<BoundarySegment> sortedBS = SortBoundarySegment(newAllFaceBS[i], innerResultPolylines[i], out isHorizontalLayout, out lineForEachEdge, out sortedBoundaryPolyline);

                    Curve[] curveUnion = new Curve[] { };
                    List<bool> isGenerateables = new List<bool>() { true, true, true, true };
                    // todo:根据BS的属性（两个Face相连的BS的setback值），来生成 List<bool> isGenerateables

                    Curve sCrvForShow = null;
                    Curve nCrvForShow = null;
                    Curve eCrvForShow = null;
                    Curve wCrvForShow = null;
                    List<Curve> gap = new List<Curve>();
                    if (sortedBoundaryPolyline.Count == 5)
                    {
                        curveUnion = GenerateLayoutLinesOnQuadBlock(sortedBoundaryPolyline,
                                                                    lineForEachEdge,
                                                                    isGenerateables,
                                                                    w,
                                                                    W,
                                                                    lMin,
                                                                    lMax,
                                                                    d,
                                                                    false,
                                                                    out sCrvForShow,
                                                                    out nCrvForShow,
                                                                    out eCrvForShow,
                                                                    out wCrvForShow,
                                                                    out gap);
                    }

                    if (sCrvForShow != null)
                    {
                        SCrvForShowList.Add(sCrvForShow);
                    }
                    if (nCrvForShow != null)
                    {
                        NCrvForShowList.Add(nCrvForShow);
                    }
                    if (eCrvForShow != null)
                    {
                        ECrvForShowList.Add(eCrvForShow);
                    }
                    if (wCrvForShow != null)
                    {
                        WCrvForShowList.Add(wCrvForShow);
                    }
                    GapCenterLineList.AddRange(gap);

                    //List<Curve> regions = new List<Curve>();
                    //for (int j = 0; j < horizontalLayoutLine.Length; j++)
                    //{
                    //    if (horizontalLayoutLine[j] != null)
                    //    {
                    //        regions.Add(GenerateSingleLayoutRegion(horizontalLayoutLine[j], w));
                    //    }
                    //    else
                    //    {
                    //        // 如果southBaseLine到northBoundaryLine的距离，小于W后，就直接让体量充满地块
                    //        regions.Add(innerResultPolylines[i].ToNurbsCurve());
                    //    }

                    //}

                    curveUnionLoL[i].AddRange(curveUnion);
                }

                DataTree<GH_Curve> ghcurves = new DataTree<GH_Curve>();
                for (int i = 0; i < curveUnionLoL.Count; i++)
                {
                    ghcurves.EnsurePath(i);
                    for (int j = 0; j < curveUnionLoL[i].Count; j++)
                    {
                        GH_Curve ghcurve = new GH_Curve(curveUnionLoL[i][j]);
                        ghcurves.Branch(i).Add(ghcurve);
                    }
                }

                //DataTree<Curve> tree = UtilityFunctions.LoLToDataTree<Curve>(horizontalLayoutLineLoL);
                //GH_Curve
                DA.SetDataTree(2, ghcurves);







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
                    SolveResults solveResults = null;
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

        public class SolveResults
        {
            public Curve[] Curves { get; set; }
            public Point3d[] Points { get; set; }
        }

        public SolveResults ComputeBrepCurveIntersection(Brep brp, Curve crv)
        {
            Curve[] curves = null;
            Point3d[] points = null;
            if (!Intersection.CurveBrep(crv, brp, GH_Component.DocumentTolerance(), out curves, out points))
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Intersection failed");
            }
            return new SolveResults
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


        //public Curve GetJoinedCurveFromBS(List<BoundarySegment> bs)
        //{
        //    List<Curve> boundaryCurves = new List<Curve>();
        //    for (int i = 0; i < bs.Count; i++)
        //    {
        //        boundaryCurves.Add(bs[i].PolylineCurve);
        //    }
        //    Curve joinedBoundaryCurve = Curve.JoinCurves(boundaryCurves)[0];

        //    return joinedBoundaryCurve;
        //}


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
                if (pts[i].Y == minY)
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
            Vector3d vector1 = sortedPts[0] - sortedPts[1];

            double rad0 = Vector3d.VectorAngle(vector0, Vector3d.XAxis);
            double rad1 = Vector3d.VectorAngle(vector1, Vector3d.YAxis);

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

        //private List<BoundarySegment> decomposeQuadBlock(List<BoundarySegment> sortedBS, 
        //                                                 Polyline sortedBoundaryPolyline, 
        //                                                 out List<Polyline> decomposedBoundaryPolyline,
        //                                                 double lMin,
        //                                                 double lMax)
        //{
        //    int sIndex = 0;
        //    int eIndex = 1;
        //    int nIndex = 2;
        //    int wIndex = 3;

            
        //    //BoundarySegment sBS = sortedBS[sIndex];

        //    // 由于输入的bs中可能含有折线，所以要重新生成不包含折线的多边形各个边
        //    List<Line> sortedLines = new List<Line>();
        //    int currentIndex = 0;
        //    while (currentIndex != sortedBoundaryPolyline.Count - 1)
        //    {
        //        Point3d pt0 = sortedBoundaryPolyline[currentIndex];
        //        Point3d pt1 = sortedBoundaryPolyline[(currentIndex + 1) % sortedBoundaryPolyline.Count];
        //        Point3d pt2 = sortedBoundaryPolyline[(currentIndex + 2) % sortedBoundaryPolyline.Count];
        //        if (IsThreePointsCollinear(pt0, pt1, pt2))
        //        {
        //            Line line = new Line(pt0, pt2);
        //            sortedLines.Add(line);

        //            currentIndex = sortedBoundaryPolyline.IndexOf(pt2);
        //        }
        //        else
        //        {
        //            Line line = new Line(pt0, pt1);
        //            sortedLines.Add(line);

        //            currentIndex = sortedBoundaryPolyline.IndexOf(pt1);
        //        }
        //    }

        //    //// 水平方向上的分解
        //    //double sL = sortedLines[sIndex].Length;
        //    //double nL = sortedLines[nIndex].Length;

        //    // 判断是选择起始点左侧的边还是右侧的边作为排布的基准线
        //    Vector3d vector0 = sortedLines.First().Direction;
        //    Vector3d vector1 = sortedLines.Last().Direction * (-1);
        //    bool isFirst;
        //    Vector3d verticalVector;
        //    Vector3d layoutDirection = GetLayoutDirection(vector0, vector1,out isFirst, out verticalVector);

        //    Curve boundary = Curve.CreateInterpolatedCurve(sortedBoundaryPolyline, 1);
        //    Curve needToTrim;
        //    if (isFirst)
        //    {
        //        needToTrim = sortedLines.First().ToNurbsCurve();
        //    }
        //    else
        //    {
        //        needToTrim = sortedLines.Last().ToNurbsCurve();
        //    }
        //    Curve newCurve = TrimBothEnd(boundary, verticalVector, needToTrim, isGenerateable, double w);



        //}

        private Curve[] GenerateLayoutLinesOnQuadBlock(//List<BoundarySegment> sortedBS, 
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
                                                       out Curve southBaseLineForShow,
                                                       out Curve northBaseLineForShow,
                                                       out Curve eastBaseLineForShow,
                                                       out Curve westBaseLineForShow,
                                                       out List<Curve> gapCenterLineForShow)
        {
            // 从确定南侧基线开始处理
            // 判断是选择起始点左侧的边还是右侧的边作为排布的基准线
            //Vector3d vector0 = lineForEachEdge.First().Direction;
            //Vector3d vector1 = lineForEachEdge.Last().Direction * (-1);
            //bool isFirst;
            Vector3d directionVector;
            Vector3d verticalVector;
            
            //Vector3d layoutDirection = GetLayoutDirection(vector0, vector1, out isFirst, out verticalVector);

            Curve boundary = Curve.CreateInterpolatedCurve(sortedBoundaryPolyline, 1);
            // 确定南侧基线
            int southSegmentIndex = 0;

            Curve southBoundaryLine = lineForEachEdge.First().ToNurbsCurve();
            directionVector = lineForEachEdge[southSegmentIndex].Direction;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve southCenterLine = TrimBothEnd(boundary, verticalVector, southBoundaryLine, Convert.ToInt32(isGenerateables[southSegmentIndex]), w, W);
            southCenterLine.Domain = new Interval(0, 1);

            Curve southBaseLine = southCenterLine.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            //sCrvForCalDistance = TrimBothEnd(boundary, verticalVector, sCrvForCalDistance, true, w, W);
            southBaseLine.Domain = new Interval(0, 1);

            // 确定北侧基线
            int northSegmentIndex = 2;

            Curve northBoundaryLine = lineForEachEdge[northSegmentIndex].ToNurbsCurve();
            directionVector = lineForEachEdge[northSegmentIndex].Direction;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve northCenterLine = TrimBothEnd(boundary, verticalVector, northBoundaryLine, Convert.ToInt32(isGenerateables[northSegmentIndex]), w, W);
            northCenterLine.Domain = new Interval(0, 1);

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
                    }
                    else
                    {
                        // 此时场地内只有一个体量，因此要考虑在南侧，中间，还是北侧？
                        Curve northCenterLineDP = northCenterLine.DuplicateCurve();
                        northCenterLineDP.Reverse();

                        Surface surface = NurbsSurface.CreateRuledSurface(southCenterLine, northCenterLineDP);
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
                        }
                        else if (pIndex == 0)
                        {
                            horizontalCenterLines.Add(southCenterLine);
                        }
                        else
                        {
                            horizontalCenterLines.Add(northCenterLine);
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
                        int pIndex = m_random.Next(3);
                        if (pIndex == 1)
                        {
                            double t = surface.Domain(1).ParameterAt(0.5);
                            Curve curve = surface.IsoCurve(0, t);
                            Vector3d vec = curve.TangentAtStart;
                            verticalVector = new Vector3d(-vec.Y, vec.X, vec.Z);
                            curve = TrimBothEnd(boundary, verticalVector, curve, 2, w, W);

                            horizontalCenterLines.Add(curve);
                        }
                        else if (pIndex == 0)
                        {
                            horizontalCenterLines.Add(southCenterLine);
                        }
                        else
                        {
                            horizontalCenterLines.Add(northCenterLine);
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

                        horizontalCenterLines.Add(southCenterLine);
                        horizontalCenterLines.Add(newNorthCenterLine);
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
                        horizontalCenterLines.Add(northCenterLine);
                    }
                    else
                    {
                        // 南北基线中间不插入新的基线
                        horizontalCenterLines.Add(southCenterLine);
                        horizontalCenterLines.Add(northCenterLine);
                    }
                }
            }



            // todo：处理东西侧方向
            // 确定东侧基线
            int eastSegmentIndex = 1;
            Curve eastBoundaryLine = lineForEachEdge[eastSegmentIndex].ToNurbsCurve();
            directionVector = lineForEachEdge[eastSegmentIndex].Direction;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve eastCenterLine = TrimBothEnd(boundary, verticalVector, eastBoundaryLine, Convert.ToInt32(isGenerateables[eastSegmentIndex]), w, W);
            eastCenterLine.Domain = new Interval(0, 1);

            Curve eastBaseLine = eastCenterLine.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            eastBaseLine.Domain = new Interval(0, 1);

            // 确定西侧基线
            int westSegmentIndex = 3;
            Curve westBoundaryLine = lineForEachEdge[westSegmentIndex].ToNurbsCurve();
            directionVector = lineForEachEdge[westSegmentIndex].Direction;
            verticalVector = new Vector3d(-directionVector.Y, directionVector.X, directionVector.Z);
            Curve westCenterLine = TrimBothEnd(boundary, verticalVector, westBoundaryLine, Convert.ToInt32(isGenerateables[westSegmentIndex]), w, W);
            westCenterLine.Domain = new Interval(0, 1);

            Curve westBaseLine = westCenterLine.Offset(Plane.WorldXY, -0.5 * w, GH_Component.DocumentTolerance(), CurveOffsetCornerStyle.None)[0];
            westBaseLine.Domain = new Interval(0, 1);

            List<Curve> verticalCenterLines = new List<Curve>();
            List<Curve> verticalCenterLinesForGap = new List<Curve>();
            // 确定是否需要在南北基线之间添加新的基线
            /* 此处不应该是 eastCenterLine 和 westCenterLine 而应该是 偏移了 0.5*w 后的 */
            Point3d a3 = eastBaseLine.PointAtStart;
            Point3d a4 = eastBaseLine.PointAtEnd;
            Point3d b3 = westBaseLine.PointAtStart;
            Point3d b4 = westBaseLine.PointAtEnd;
            int indexBetween34BaseLines = 0;
            double minDistanceBetween34BaseLines = MinDistanceBetweenTwoLineSegment(a3, a4, b3, b4, out indexBetween34BaseLines);

            if (minDistanceBetween34BaseLines <= 0)
            {
                // 无法插入Gap，verticalCenterLinesForGap 不添加新的 insectCurves

            }
            // 两线段之间的距离不足以放置新的基线（进行TweenCurve）
            else if (minDistanceBetween34BaseLines > 0 && minDistanceBetween34BaseLines < (lMin + 2 * d)) 
            {
                // 两线段之间的距离不足以放置新的基线（进行TweenCurve）
                // 就要考虑两线段之间的最长距离，能够放置几根基线
                // 基于能够放置基线的位置，来对原有的westBaseLine进行缩短，进而形成TweenCurve式的基线

                // 计算两线段之间的最长距离
                int indexBetween34BaseLinesMax = 0;
                double maxDistanceBetween34BaseLine = MaxDistanceBetweenTwoLineSegment(a3, a4, b3, b4, out indexBetween34BaseLinesMax);

                if (minDistanceBetween34BaseLines - maxDistanceBetween34BaseLine < 0.001
                    && minDistanceBetween34BaseLines - maxDistanceBetween34BaseLine > -0.001)
                {
                    // EW两个边平行
                    // 并且都小于 lMin * 2 + d ，无法插入Gap
                    // verticalCenterLinesForGap 不添加新的 insectCurves
                }
                else
                {
                    // EW两个边不平行
                    // 计算最长距离可以塞下几个 lMax ,以能塞下几个 lMax 为基准
                    maxDistanceBetween34BaseLine -= d;
                    int insertCountForLMax = (int)(Math.Floor(maxDistanceBetween34BaseLine / (lMax + d)));
                    // 因为是求gap的位置，所以+1
                    int insertGapCount = insertCountForLMax + 1;

                    if (insertCountForLMax > 0)
                    {
                        // 如果最长距离能够塞下 lMax
                        // 计算两线段之间的最短距离那一侧的边，要向最长距离那一侧的边 setback 多少距离
                        Vector3d a3a4 = new Vector3d(a4 - a3);
                        Vector3d a4b3 = new Vector3d(b3 - a4);
                        Vector3d b3b4 = new Vector3d(b4 - b3);

                        Curve trimedEastBaseLine;
                        Curve trimedWestBaseLine;
                        if (a4b3.Length >= insertCountForLMax * (lMin + d) + d)
                        {
                            trimedEastBaseLine = eastBaseLine;
                            trimedWestBaseLine = westBaseLine;
                        }
                        else
                        {
                            double setbackE = 0;
                            double setbackW = CalSetback(a3a4, a4b3, b3b4, a4b3.Length, insertCountForLMax * (lMin + d) + d, out setbackE);

                            trimedEastBaseLine = eastBaseLine.Trim(CurveEnd.End, setbackE);
                            trimedWestBaseLine = westBaseLine.Trim(CurveEnd.Start, setbackW);
                            trimedEastBaseLine.Domain = new Interval(0, 1);
                            trimedWestBaseLine.Domain = new Interval(0, 1);
                        }
                        

                        // 基于能够放置基线的位置，来对原有的westBaseLine进行缩短，进而形成TweenCurve式的基线
                        List<Curve> insectGapCurves = ConstructTweenCurve(insertGapCount, lMax, d, trimedEastBaseLine, trimedWestBaseLine, true, W);
                        verticalCenterLinesForGap.AddRange(insectGapCurves);
                    }
                    else
                    {
                        // 如果最长距离不能塞下 lMax，就看能塞下几个 lMin
                        //int insertGapCount;
                        int insertCountForLMin = (int)(Math.Floor(maxDistanceBetween34BaseLine / (lMin + d)));
                        maxDistanceBetween34BaseLine += d;

                        // 通过在最大插入数量和最小插入数量之间取随机值，得到最终插入Gap的数量 insertCount
                        //int insertCount;
                        if (insertCountForLMin != 0)
                        {
                            // 此时最大可以插入 insertCountForLMin + 1 个Gap ，即分成 insertCountForLMin 个体量
                            int maxInsertGapCount = insertCountForLMin + 1;
                            // 最小可以插入 minInsertCount 个Gap
                            int minInsertGapCount;
                            if (maxDistanceBetween34BaseLine > lMax)
                            {
                                minInsertGapCount = 1;
                            }
                            else
                            {
                                minInsertGapCount = 0;
                            }

                            if (minInsertGapCount == maxInsertGapCount)
                            {
                                insertGapCount = minInsertGapCount;
                            }
                            else
                            {
                                insertGapCount = m_random.Next(minInsertGapCount, maxInsertGapCount + 1);
                            }
                        }
                        else
                        {
                            if (maxDistanceBetween34BaseLine > lMax)
                            {
                                insertGapCount = 1;
                            }
                            else
                            {
                                insertGapCount = 0;
                            }
                        }

                        if (insertGapCount == 0)
                        {
                            // 最终判断为在最长距离一侧的边，无法插入Gap

                            // verticalCenterLinesForGap 不添加新的 insectCurves
                        }
                        else
                        {
                            // 最终判断为在最长距离一侧的边，可以插入 insertCount 个 Gap
                            // 计算两线段之间的最短距离那一侧的边，要向最长距离那一侧的边 setback 多少距离
                            Vector3d a3a4 = new Vector3d(a4 - a3);
                            Vector3d a4b3 = new Vector3d(b3 - a4);
                            Vector3d b3b4 = new Vector3d(b4 - b3);
                            Curve trimedEastBaseLine;
                            Curve trimedWestBaseLine;
                            if (a4b3.Length >= (insertGapCount - 1) * (lMin + d) + d)
                            {
                                trimedEastBaseLine = eastBaseLine;
                                trimedWestBaseLine = westBaseLine;
                            }
                            else
                            {
                                double setbackE = 0;
                                double setbackW = CalSetback(a3a4, a4b3, b3b4, a4b3.Length, (insertGapCount - 1) * (lMin + d) + d, out setbackE);

                                trimedEastBaseLine = eastBaseLine.Trim(CurveEnd.End, setbackE);
                                trimedWestBaseLine = westBaseLine.Trim(CurveEnd.Start, setbackW);
                                trimedEastBaseLine.Domain = new Interval(0, 1);
                                trimedWestBaseLine.Domain = new Interval(0, 1);
                            }

                            // 基于能够放置基线的位置，来对原有的westBaseLine进行缩短，进而形成TweenCurve式的基线
                            List<Curve> insectCurves = ConstructTweenCurve(insertGapCount, lMin, d, trimedEastBaseLine, trimedWestBaseLine, true, W);
                            verticalCenterLinesForGap.AddRange(insectCurves);
                        }
                    }
                }
            }
            // 两线段之间的距离足以放置新的基线（进行TweenCurve）
            else
            {
                minDistanceBetween34BaseLines -= d;
                int insertLMinCount = (int)(Math.Floor(minDistanceBetween34BaseLines / (lMin + d)));
                // 因为是求gap的位置，所以+1
                int insertGapCount = insertLMinCount + 1;

                if (insertLMinCount > 0)
                {
                    // 东西基线中间插入新的分割线
                    double sum = insertLMinCount * (lMin + d) + d;
                    List<double> eachLengths = new List<double>();
                    eachLengths.Add(0.5 * d);
                    for (int i = 1; i < insertGapCount; i++)
                    {
                        double num = eachLengths[i - 1] + lMin + d;
                        eachLengths.Add(num);
                    }

                    List<double> factors = new List<double>();
                    for (int i = 0; i < insertGapCount; i++)
                    {
                        double factor = eachLengths[i] / sum;
                        factors.Add(factor);
                    }

                    Curve trimedEastBaseLine;
                    Curve trimedWestBaseLine;
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
                            eastBaseLine.ClosestPoint(b3, out t);
                            trimedWestBaseLine = westBaseLine;
                            trimedEastBaseLine = eastBaseLine.Trim(0.0, t);
                        }
                        else if (indexBetween34BaseLines == 2)
                        {
                            // 与b3b4垂直
                            double t;
                            westBaseLine.ClosestPoint(a4, out t);
                            trimedEastBaseLine = eastBaseLine;
                            trimedWestBaseLine = westBaseLine.Trim(t, 1.0);
                        }
                        else
                        {
                            // 都不垂直
                            trimedEastBaseLine = eastBaseLine;
                            trimedWestBaseLine = westBaseLine;
                        }
                    }
                    else
                    {
                        // index为奇数，左端最近
                        if (indexBetween34BaseLines == 1)
                        {
                            // 与a3a4垂直
                            double t;
                            eastBaseLine.ClosestPoint(b4, out t);
                            trimedWestBaseLine = westBaseLine;
                            trimedEastBaseLine = eastBaseLine.Trim(t, 1.0);
                        }
                        else if (indexBetween34BaseLines == 3)
                        {
                            // 与b3b4垂直
                            double t;
                            westBaseLine.ClosestPoint(a3, out t);
                            trimedEastBaseLine = eastBaseLine;
                            trimedWestBaseLine = westBaseLine.Trim(0.0, t);
                        }
                        else
                        {
                            // 都不垂直
                            trimedEastBaseLine = eastBaseLine;
                            trimedWestBaseLine = westBaseLine;
                        }
                    }

                    trimedWestBaseLine.Reverse();
                    Surface surface = NurbsSurface.CreateRuledSurface(trimedEastBaseLine, trimedWestBaseLine);
                    List<Curve> insectCurves = new List<Curve>();
                    for (int i = 0; i < factors.Count; i++)
                    {
                        double t = surface.Domain(1).ParameterAt(factors[i]);
                        Curve curve = surface.IsoCurve(0, t);
                        curve.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
                        insectCurves.Add(curve);
                    }
                    surface.Dispose();

                    verticalCenterLinesForGap.AddRange(insectCurves);
                }
                else
                {
                    // verticalCenterLinesForGap 不添加新的 insectCurves
                }
            }







            // Temp
            southBaseLineForShow = southBaseLine;
            northBaseLineForShow = northBaseLine;
            eastBaseLineForShow = eastBaseLine;
            westBaseLineForShow = westBaseLine;
            gapCenterLineForShow = verticalCenterLinesForGap;



            List<Curve> singleRegions = new List<Curve>();
            if (horizontalCenterLines.Count != 0)
            {
                for (int i = 0; i < horizontalCenterLines.Count; i++)
                {
                    singleRegions.Add(GenerateSingleLayoutRegion(horizontalCenterLines[i], w));
                }
            }
            else
            {
                // 如果southBaseLine到northBoundaryLine的距离，小于W后，就直接让体量充满地块
                // 此时horizontalCenterLines这个列表中没有元素，singleRegions列表只有一个元素
                singleRegions.Add(sortedBoundaryPolyline.ToNurbsCurve());
            }

            // Temp添加GapCurve
            if (verticalCenterLinesForGap.Count != 0)
            {
                for (int i = 0; i < verticalCenterLinesForGap.Count; i++)
                {
                    singleRegions.Add(GenerateSingleLayoutRegion(verticalCenterLinesForGap[i], w));
                }
            }


            Curve[] curveUnion;
            if (singleRegions.Count > 1)
            {
                curveUnion = Curve.CreateBooleanUnion(singleRegions, GH_Component.DocumentTolerance());
                return curveUnion;
            }
            else
            {
                return singleRegions.ToArray();
            }
        }

        private Curve GenerateSingleLayoutRegion(Curve layoutLine, double w)
        {
            Vector3d direction = layoutLine.TangentAtStart;
            direction.Unitize();

            Vector3d moveDirection0 = Vector3d.CrossProduct(direction, new Vector3d(0, 0, -1));
            moveDirection0.Unitize();
            moveDirection0 = moveDirection0 * (w * 0.5);

            Vector3d moveDirection1 = new Vector3d(moveDirection0);
            moveDirection1 *= -1;

            Point3d start = layoutLine.PointAtStart;
            Point3d end = layoutLine.PointAtEnd;

            List<Point3d> regionPts = new List<Point3d>();
            regionPts.Add(start + moveDirection0);
            regionPts.Add(start + moveDirection1);
            regionPts.Add(end + moveDirection1);
            regionPts.Add(end + moveDirection0);

            regionPts.Add(start + moveDirection0);

            Curve crv = Curve.CreateInterpolatedCurve(regionPts, 1);

            return crv;
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
                Curve offset0 = current.Offset(Plane.WorldXY, w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];
                Curve offset1 = current.Offset(Plane.WorldXY, -w * 0.5, tolerance, CurveOffsetCornerStyle.None)[0];

                new0 = offset0.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                newCurrent = current.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
                new1 = offset1.Extend(CurveEnd.Both, boxLength * 0.5, CurveExtensionStyle.Line);
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
                    return null;
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
                    return null;
                }
            }
            else
            {
                return null;
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

            return trimedCurve;
        }

        private List<Curve> ConstructTweenCurve(int insertGapCount, double lMinOrlMax, double d, Curve baseLine0, Curve baseLine1, bool isGap, double W)
        {
            //double sum = (insertLMinOrLMaxCount * 1) * lMinOrlMax + insertLMinOrLMaxCount * d;
            double sum = (insertGapCount - 1) * (lMinOrlMax + d) + d;
            List<double> eachLengths = new List<double>();
            eachLengths.Add(0.5 * d);
            for (int i = 1; i < insertGapCount; i++)
            {
                double num = eachLengths[i - 1] + lMinOrlMax + d;
                eachLengths.Add(num);
            }

            List<double> factors = new List<double>();
            for (int i = 0; i < insertGapCount; i++)
            {
                double factor = eachLengths[i] / sum;
                factors.Add(factor);
            }

            baseLine1.Reverse();
            Surface surface = NurbsSurface.CreateRuledSurface(baseLine0, baseLine1);
            List<Curve> insectCurves = new List<Curve>();
            for (int i = 0; i < factors.Count; i++)
            {
                double t = surface.Domain(1).ParameterAt(factors[i]);
                Curve curve = surface.IsoCurve(0, t);
                if (isGap)
                {
                    curve.Extend(CurveEnd.Both, W, CurveExtensionStyle.Line);
                }
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
        /// 两个线段之间的最大距离
        /// </summary>
        /// <param name="a1"></param>
        /// <param name="a2"></param>
        /// <param name="b1"></param>
        /// <param name="b2"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        private double MaxDistanceBetweenTwoLineSegment(Point3d a1,Point3d a2,Point3d b1,Point3d b2,out int index)
        {
            // 因为这里用到的两条线段是反向的，所以是a1b2,a2b1
            double a1b2Dis2 = Dis2(a1, b2);
            double a2b1Dis2 = Dis2(a2, b1);

            if (a1b2Dis2 - a2b1Dis2 > 0.5 * GH_Component.DocumentTolerance())
            {
                // 最大值的位置在左端
                index = 0;
                return Math.Sqrt(a1b2Dis2);
            }
            else if (a1b2Dis2 - a2b1Dis2 <= 0.0001 && a1b2Dis2 - a2b1Dis2 >= -0.0001)
            {
                // 两线段平行
                index = 2;
                return Math.Sqrt(a1b2Dis2);
            }
            else
            {
                // 最大值的位置在右端
                index = 1;
                return Math.Sqrt(a2b1Dis2);
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

        //private List<Point3d> GenerateSplitPointOnHorizontalDirection(double lMin, double lMax, double w, Line line)
        //{
        //    List<Point3d> splitPoints = new List<Point3d>();
        //    double L = line.Length;
        //    if (L < lMin)
        //    {
        //        // L∈(0,lMin)，需要退线后，再次执行分解操作

        //        return null;
        //    }
        //    else if (L >= lMin && L <= lMax)
        //    {
        //        // L∈[lMin,lMax]，不需要分解
        //        return new List<Point3d>();
        //    }
        //    else if (L >= 2 * lMin + w && L <= 2 * lMax + w)
        //    {
        //        // L∈[2 * lMin + w,2 * lMax + w]，分解为两个
        //        double[] lengths = new double[2] { lMin + w / 2, lMin + w / 2 };

        //        double remaining = L - lengths[0] - lengths[1];

        //        int modulus = (int)Math.Truncate(remaining / 3);
        //        double remainder = remaining % 3;

        //        int[] moduluses = new int[2] { 0, 0 };
        //        double[] additions = new double[2] { 0, 0 };
        //        if (modulus == 0 || modulus == 1)
        //        {

        //            int indexForModulus = m_random.Next(0, 2);
        //            int indexForRemainder = m_random.Next(0, 2);
        //            additions[indexForModulus] += modulus * 3;
        //            additions[indexForRemainder] += remainder;
        //        }
        //        else
        //        {

        //            moduluses[0] = m_random.Next(0, modulus);
        //            moduluses[1] = modulus - moduluses[0];

        //            int indexForModulus0 = m_random.Next(0, 2);
        //            int indexForModulus1 = 1 - indexForModulus0;
        //            int indexForRemainder = m_random.Next(0, 2);

        //            additions[indexForModulus0] += moduluses[0] * 3;
        //            additions[indexForModulus1] += moduluses[1] * 3;
        //            additions[indexForRemainder] += remainder;
        //        }
        //        for (int i = 0; i < lengths.Length; i++)
        //        {
        //            lengths[i] += additions[i];
        //        }

        //        double t = lengths[0] / line.Length;
        //        Point3d point = line.PointAt(t);

        //        splitPoints.Add(point);
        //    }
        //    else if (L >= 3 * lMin + 2 * w && L <= 3 * lMax + 2 * w)
        //    {
        //        // L∈[3 * lMin + 2 * w,3 * lMax + 2 * w]，分解为三个

        //        double[] lengths = new double[3] { lMin + w / 2, lMin + w, lMin + w / 2 };

        //        double remaining = L - lengths[0] - lengths[1] - lengths[2];
        //        int modulus = (int)Math.Truncate(remaining / 3);
        //        double remainder = remaining % 3;

        //        int[] moduluses = new int[3] { 0, 0, 0 };
        //        double[] additions = new double[3] { 0, 0, 0 };
        //        if (modulus == 0 || modulus == 1)
        //        {
        //            int indexForModulus = m_random.Next(0, 3);
        //            int indexForRemainder = m_random.Next(0, 3);
        //            additions[indexForModulus] += modulus * 3;
        //            additions[indexForRemainder] += remainder;
        //        }
        //        else
        //        {
        //            moduluses[0] = m_random.Next(0, modulus);
        //            moduluses[1] = m_random.Next(moduluses[0], modulus);
        //            moduluses[2] = modulus - moduluses[0] - moduluses[1];

        //            List<int> indexForModuluses = new List<int>();
        //            List<int> list = new List<int>() { 0, 1, 2 };
        //            int n0 = m_random.Next(list.Count - 1);
        //            indexForModuluses.Add(list[n0]);
        //            list.RemoveAt(n0);
        //            int n1 = m_random.Next(list.Count - 1);
        //            indexForModuluses.Add(list[n1]);
        //            list.RemoveAt(n1);
        //            indexForModuluses.Add(list[0]);

        //            int indexForRemainder = m_random.Next(0, 3);

        //            for (int i = 0; i < indexForModuluses.Count; i++)
        //            {
        //                additions[indexForModuluses[i]] += moduluses[i] * 3;
        //            }
        //            additions[indexForRemainder] += remainder;
        //        }
        //        for (int i = 0; i < lengths.Length; i++)
        //        {
        //            lengths[i] += additions[i];
        //        }

        //        double t0 = lengths[0] / line.Length;
        //        double t1 = lengths[1] / line.Length;
        //        Point3d pt0 = line.PointAt(t0);
        //        Point3d pt1 = line.PointAt(t1);

        //        splitPoints.Add(pt0);
        //        splitPoints.Add(pt1);
        //    }
        //    else
        //    {
        //        // 先进行二分，然后迭代直到符合上面几个条件中的一个
        //        double sample = m_random.NextDouble();
        //        // 0.4~0.6
        //        double scaled = sample * 0.2 + 0.4;
        //        Point3d pt = line.PointAt(scaled);
        //        List<Line> newLines = new List<Line>();
        //        newLines.Add(new Line(line.From, pt));
        //        newLines.Add(new Line(pt, line.To));

        //        List<List<Point3d>> splitPts = new List<List<Point3d>>();
        //        for (int i = 0; i < newLines.Count; i++)
        //        {
        //            splitPts.Add(new List<Point3d>());
        //            splitPts[i].AddRange(GenerateSplitPointOnHorizontalDirection(lMin, lMax, w, newLines[i]));
        //        }

        //        splitPoints.AddRange(splitPts[0]);
        //        splitPoints.Add(pt);
        //        splitPoints.AddRange(splitPts[1]);
        //    }

        //    return splitPoints;
        //}

        //public void CalEachBSDomain(List<BoundarySegment> bs, Curve joinedCurve)
        //{
        //    Interval allDomain = joinedCurve.Domain;

        //    for (int i = 0; i < bs.Count; i++)
        //    {
        //        Interval spanDomain = joinedCurve.SpanDomain(i);

        //        bs[i].SubDomain = new Interval(spanDomain);
        //    }
        //}

        //public Polyline GenerateInnerZone(List<BoundarySegment> bs, int bsIndex)
        //{

        //}

        //public Polyline GenerateInnerOffset(List<BoundarySegment> bs)
        //{
        //    List<Curve> boundaryCurves = new List<Curve>();
        //    for (int i = 0; i < bs.Count; i++)
        //    {
        //        boundaryCurves.Add(bs[i].PolylineCurve);
        //    }
        //    Curve joinedBoundaryCurve = Curve.JoinCurves(boundaryCurves)[0];
        //    BoundingBox box = joinedBoundaryCurve.GetBoundingBox(Rhino.Geometry.Plane.WorldXY);
        //    double dx = box.Max.X - box.Min.X;
        //    double dy = box.Max.Y - box.Min.Y;
        //    double extendedValue = Math.Max(dx, dy) / 2;

        //    // 创建offset的退线区域
        //    List<Curve> cutCurves = new List<Curve>();
        //    for (int i = 0; i < length; i++)
        //    {

        //    }

        //}

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
            for (int i = 0; i < SCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(SCrvForShowList[i], Color.Red, 2);
            }
            for (int i = 0; i < NCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(NCrvForShowList[i], Color.Green, 2);
            }
            for (int i = 0; i < ECrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(ECrvForShowList[i], Color.Blue, 2);
            }
            for (int i = 0; i < WCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(WCrvForShowList[i], Color.Yellow, 2);
            }
            for (int i = 0; i < GapCenterLineList.Count; i++)
            {
                args.Display.DrawCurve(GapCenterLineList[i], Color.Black, 3);
            }
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
            for (int i = 0; i < SCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(SCrvForShowList[i], Color.Red, 2);
            }
            for (int i = 0; i < NCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(NCrvForShowList[i], Color.Green, 2);
            }
            for (int i = 0; i < ECrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(ECrvForShowList[i], Color.Blue, 2);
            }
            for (int i = 0; i < WCrvForShowList.Count; i++)
            {
                args.Display.DrawCurve(WCrvForShowList[i], Color.Yellow, 2);
            }
            for (int i = 0; i < GapCenterLineList.Count; i++)
            {
                args.Display.DrawCurve(GapCenterLineList[i], Color.Black, 3);
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
            get { return new Guid("4b850a3f-93e2-44a8-b44a-79f109cd61c3"); }
        }
    }
}