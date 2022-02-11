using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel.Geometry;
using Grasshopper.Kernel.Geometry.ConvexHull;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Linq;
using System.Drawing;
using System.Collections.Generic;

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
        }

        private int Thickness;
        private List<Polyline> InnerResultPolyline;

        private List<TextDot> InnerNodeTextDots;

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

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetData<PlanktonMesh>("newDualPlanktonMesh", ref dual)
                && DA.GetDataList("SetbackList", setbackList))
            {
                DA.GetDataTree<GH_Integer>("WhichPairSetbackNeedToChange", out pairNeedToChange);
                DA.GetDataTree<GH_Number>("SetbackValueNeedToChange", out setbackValueNeedToChange);

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

                        // 判断首尾两个点位置处的曲率，即是否在首尾处转折
                        double first, last;
                        GetSlopeAtStartAndEnd(dual, allFaceHIndex, i, j, currPoints, out first, out last);

                        BoundarySegment bs = new BoundarySegment(currPoints, allFaceHIndex[i][j].ToArray(), allFaceFIndex[i][j], allFaceAdjacentFIndex[i][j]);
                        if (first - 0 < 0.5 * Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
                        {
                            bs.IsZeroCurvatureAtStart = true;
                        }
                        else
                        {
                            bs.IsZeroCurvatureAtStart = false;
                        }
                        if (last - 0 < 0.5 * Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
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

        private void GetSlopeAtStartAndEnd(PlanktonMesh dual, List<List<List<int>>> allFaceHIndex, int i, int j, List<Point3d> currPoints, out double first, out double last)
        {
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

            Vector3d dxdy01 = new Vector3d(pt1 - pt0);
            Vector3d dxdy12 = new Vector3d(pt2 - pt1);
            Vector3d dxdy34 = new Vector3d(pt4 - pt3);
            Vector3d dxdy45 = new Vector3d(pt5 - pt4);

            first = dxdy01.X * dxdy12.Y - dxdy01.Y * dxdy12.X;
            last = dxdy34.X * dxdy45.Y - dxdy34.Y * dxdy45.X;
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