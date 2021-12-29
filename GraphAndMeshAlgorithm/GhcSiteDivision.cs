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
using System.Text;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcSiteDivision : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcSiteDivision class.
        /// </summary>
        public GhcSiteDivision()
          : base("GhcSiteDivision", "SiteDivision",
              "进行场地划分",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            RealisticSiteFaceCenter = new List<TextDot>();
            DebugPolylinesForGreenColor = new List<Polyline>();
            DebugPolylinesForYellowColor = new List<Polyline>();

            DebugPolylinesForGreenColorHatch = new List<Hatch>();
            DebugPolylinesForYellowColorHatch = new List<Hatch>();
        }

        private List<TextDot> RealisticSiteFaceCenter;

        private List<Polyline> DebugPolylinesForGreenColor;
        private List<Polyline> DebugPolylinesForYellowColor;

        private List<Hatch> DebugPolylinesForGreenColorHatch;
        private List<Hatch> DebugPolylinesForYellowColorHatch;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DGHM", "生成的对偶图", GH_ParamAccess.item);
            pManager.AddCurveParameter("ReorganizedBoundary", "RB", "经过重新组织的边界polyline", GH_ParamAccess.item);
            pManager.AddGenericParameter("SubBoundarySegments", "SS", "分割过后形成的子BoundarySegment对象", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("DeBug Polyline", "DP", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            //Curve reorganizedBoundaryCurve = null;

            //// 因为BoundarySegment没有实现IGH_Goo接口，所以作为GH_Structure无法传递
            //// 直接传递List<List<BoundarySegment>>也是不行的，所以用以下方式进行传递
            //List<GH_ObjectWrapper> objList = new List<GH_ObjectWrapper>();
            //DA.GetDataList(2, objList);

            //List<List<FaceEdgeSegment>> subBoundarySegments = new List<List<FaceEdgeSegment>>();
            //for (int i = 0; i < objList.Count; i++)
            //{
            //    subBoundarySegments.Add(new List<FaceEdgeSegment>());
            //    if (objList[i] != null)
            //    {
            //        List<FaceEdgeSegment> currentSubBoundarySegment = objList[i].Value as List<FaceEdgeSegment>;
            //        subBoundarySegments[i].AddRange(currentSubBoundarySegment);
            //    }
            //}

            ////GH_Structure<GH_Integer> gh_structure = new GH_Structure<GH_Integer>();
            ////DataTree<BoundarySegment> sortedSubBoundarySegmentDT = new DataTree<BoundarySegment>();
            ////List<List<BoundarySegment>> subBoundarySegments = new List<List<BoundarySegment>>();

            //if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
            //    && DA.GetData<Curve>("ReorganizedBoundary", ref reorganizedBoundaryCurve)
            //    && subBoundarySegments.Count > 0)
            //{
            //    DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);
            //    PlanktonMesh D = dualGraphWithHMDP.DualPlanktonMesh;
            //    // 深拷贝对象，以便对对象进行修改
            //    PlanktonMesh DDeepCopy = new PlanktonMesh(D);

            //    // List<string> debug = UtilityFunctions.PrintFacesHalfedges(D);
            //    // List<string> debug2 = UtilityFunctions.PrintHalfedgeStartAndEnd(D);

            //    // 对输入的Curve类型的reorganizedBoundary进行类型转换，转换成Curve类的子类Polyline
            //    Polyline reorganizedBoundary = null;
            //    reorganizedBoundaryCurve.TryGetPolyline(out reorganizedBoundary);

            //    // 记录所有边缘第一层级的Face
            //    // 记录新生成的面，以及新生成的面和之前的dFace的关系
            //    List<List<int>> realisticSiteFaceVerticesForAllFaceForGreenColor = new List<List<int>>();
            //    List<List<Point3d>> realisticSiteFacePointsForAllFaceForGreenColor = new List<List<Point3d>>();

            //    List<int> realisticSiteFaceIndexForGreenColor = new List<int>();


            //    // 记录所有第二层级的Face
            //    List<List<int>> realisticSiteFaceVerticesForAllFaceForYellowColor = new List<List<int>>();
            //    List<List<Point3d>> realisticSiteFacePointsForAllFaceForYellowColor = new List<List<Point3d>>();

            //    List<int> realisticSiteFaceIndexForYellowColor = new List<int>();



            //    // 保存已经被绘制过的halfedge所对应的FE
            //    List<FaceEdgeSegment> drawnFE = new List<FaceEdgeSegment>();

            //    // key：最开始的D中halfedge的首尾点index；value：对应的已经绘制过segment的halfedge的首尾点index
            //    List<int[]> drawnDHV = new List<int[]>();

            //    // 废弃掉的半边
            //    List<int[]> abandonedDHV = new List<int[]>();

            //    Dictionary<int, Point3d> vi_Point3d = new Dictionary<int, Point3d>();

            //    #region 找到边缘的面（有连续的三边可以确定）
            //    // 找到边缘的面（有连续的三边可以确定），因为都是四边面，所以但凡是有三边的，都是连续的三边
            //    List<int> bFaceIndexThreeSideDetermined = new List<int>();
            //    for (int i = 0; i < DDeepCopy.Faces.Count; i++)
            //    {
            //        if (DDeepCopy.Faces.NakedEdgeCount(i) >= 3)
            //        {
            //            bFaceIndexThreeSideDetermined.Add(i);
            //        }
            //    }
            //    #endregion
            //    #region 找到边缘的面（有连续的两边可以确定）
            //    // 找到边缘的面（有连续的两边可以确定）
            //    List<int> bFaceIndexTwoSidedDetermined = new List<int>();
            //    for (int i = 0; i < DDeepCopy.Faces.Count; i++)
            //    {
            //        if (DDeepCopy.Faces.NakedEdgeCount(i) == 2)
            //        {
            //            List<int> nakedEdges = GetNakedEdges(DDeepCopy, i);
            //            if (IsNakedEdgesContinuous(DDeepCopy, nakedEdges))
            //            {
            //                bFaceIndexTwoSidedDetermined.Add(i);
            //            }
            //        }
            //    }
            //    #endregion

            //    #region 找到与边缘的面（有连续的三边可以确定）邻接的面，及公共边
            //    List<int> aFaceIndexThreeSideDetermined = new List<int>();
            //    // List<int> pHalfedgeThreeSideDetermined = new List<int>();
            //    List<int[]> publicHVThreeSideDetermined = new List<int[]>();
            //    for (int i = 0; i < bFaceIndexThreeSideDetermined.Count; i++)
            //    {
            //        // int publicHalfedge = -1;
            //        int adjacentFace = -1;
            //        int[] hs = DDeepCopy.Faces.GetHalfedges(bFaceIndexThreeSideDetermined[i]);
            //        for (int j = 0; j < hs.Length; j++)
            //        {
            //            if (DDeepCopy.Halfedges[DDeepCopy.Halfedges.GetPairHalfedge(hs[j])].AdjacentFace != -1)
            //            {
            //                // publicHalfedge = hs[j];
            //                int[] publicHV = new int[2] { DDeepCopy.Halfedges[hs[j]].StartVertex, DDeepCopy.Halfedges.EndVertex(hs[j]) };
            //                adjacentFace = DDeepCopy.Halfedges[DDeepCopy.Halfedges.GetPairHalfedge(hs[j])].AdjacentFace;

            //                aFaceIndexThreeSideDetermined.Add(adjacentFace);
            //                // pHalfedgeThreeSideDetermined.Add(publicHalfedge);
            //                publicHVThreeSideDetermined.Add(publicHV);
            //            }
            //        }
            //    }
            //    #endregion

            //    #region 找到与边缘的面（有连续的两边可以确定）邻接的面，及公共边
            //    List<List<int>> aFaceIndexTwoSideDetermined = new List<List<int>>();
            //    // List<List<int>> pHalfedgeTwoSideDetermined = new List<List<int>>();
            //    List<List<int[]>> publicHVTwoSideDetermined = new List<List<int[]>>();
            //    for (int i = 0; i < bFaceIndexTwoSidedDetermined.Count; i++)
            //    {
            //        aFaceIndexTwoSideDetermined.Add(new List<int>());
            //        // pHalfedgeTwoSideDetermined.Add(new List<int>());
            //        publicHVTwoSideDetermined.Add(new List<int[]>());

            //        // List<int> publicHalfedges = new List<int>();
            //        List<int> adjacentFaces = new List<int>();
            //        int[] hs = DDeepCopy.Faces.GetHalfedges(bFaceIndexTwoSidedDetermined[i]);
            //        List<int[]> publicHVs = new List<int[]>();
            //        for (int j = 0; j < hs.Length; j++)
            //        {
            //            if (DDeepCopy.Halfedges[DDeepCopy.Halfedges.GetPairHalfedge(hs[j])].AdjacentFace != -1)
            //            {
            //                // publicHalfedges.Add(hs[j]);
            //                int[] publicHV = new int[2] { DDeepCopy.Halfedges[hs[j]].StartVertex, DDeepCopy.Halfedges.EndVertex(hs[j]) };
            //                adjacentFaces.Add(DDeepCopy.Halfedges[DDeepCopy.Halfedges.GetPairHalfedge(hs[j])].AdjacentFace);

            //                publicHVs.Add(publicHV);
            //            }
            //        }

            //        aFaceIndexTwoSideDetermined[i].AddRange(adjacentFaces);
            //        // pHalfedgeTwoSideDetermined[i].AddRange(publicHalfedges);
            //        publicHVTwoSideDetermined[i].AddRange(publicHVs);
            //    }
            //    #endregion

            //    #region 清除邻接面中的重复值，清除邻接面中与边缘面的重复
            //    // 找到与边缘面邻接的面
            //    List<int> aFaceIndex = new List<int>();
            //    // List<int> aFacePublicHalfedges = new List<int>();
            //    List<List<int[]>> aFacePublicHVs = new List<List<int[]>>();
            //    // 边缘的面（有连续的三边可以确定）邻接的面一定不会重复，故直接添加
            //    aFaceIndex.AddRange(aFaceIndexThreeSideDetermined);
            //    // aFacePublicHalfedges.AddRange(pHalfedgeThreeSideDetermined);
            //    //aFacePublicHVs.AddRange(publicHVThreeSideDetermined);
            //    for (int i = 0; i < aFaceIndex.Count; i++)
            //    {
            //        aFacePublicHVs.Add(new List<int[]>());
            //        aFacePublicHVs[i].Add(publicHVThreeSideDetermined[i]);
            //    }

            //    List<int> bFace = bFaceIndexThreeSideDetermined.Union(bFaceIndexTwoSidedDetermined).ToList<int>();
            //    List<int> cleanedAFaceIndexTwoSideDetermined = new List<int>();
            //    // List<int> cleanedPHalfedgeTwoSideDetermined = new List<int>();
            //    List<int[]> cleanedPublicHVTwoSideDetermined = new List<int[]>();
            //    for (int i = 0; i < aFaceIndexTwoSideDetermined.Count; i++)
            //    {
            //        for (int j = 0; j < aFaceIndexTwoSideDetermined[i].Count; j++)
            //        {
            //            // 如果找到与边缘的面（有连续的两边可以确定）邻接的面，没有在所有的边缘面中出现，就添加
            //            if (bFace.Contains(aFaceIndexTwoSideDetermined[i][j]))
            //            {
            //                continue;
            //            }
            //            else
            //            {
            //                cleanedAFaceIndexTwoSideDetermined.Add(aFaceIndexTwoSideDetermined[i][j]);
            //                // cleanedPHalfedgeTwoSideDetermined.Add(pHalfedgeTwoSideDetermined[i][j]);
            //                cleanedPublicHVTwoSideDetermined.Add(publicHVTwoSideDetermined[i][j]);
            //            }
            //        }
            //    }

            //    // 向aFace列表中添加不重复的元素
            //    for (int i = 0; i < cleanedAFaceIndexTwoSideDetermined.Count; i++)
            //    {
            //        if (!aFaceIndexThreeSideDetermined.Contains(cleanedAFaceIndexTwoSideDetermined[i]))
            //        {
            //            aFaceIndex.Add(cleanedAFaceIndexTwoSideDetermined[i]);
            //            aFacePublicHVs.Add(new List<int[]>());
            //            aFacePublicHVs.Last().Add(cleanedPublicHVTwoSideDetermined[i]);
            //        }
            //        else
            //        {
            //            int index = aFaceIndexThreeSideDetermined.IndexOf(cleanedAFaceIndexTwoSideDetermined[i]);
            //            aFacePublicHVs[index].Add(cleanedPublicHVTwoSideDetermined[i]);
            //        }
            //    }

            //    //aFaceIndex.AddRange(cleanedAFaceIndexTwoSideDetermined);
            //    //// aFacePublicHalfedges.AddRange(cleanedPHalfedgeTwoSideDetermined);
            //    //aFacePublicHVs.AddRange(cleanedPublicHVTwoSideDetermined);
            //    #endregion

            //    #region 处理边缘的面（有连续的三边可以确定），以及边缘的面（有连续的两边可以确定）

            //    // 对于每个边缘的有三边可以确定的面，找到公共边的对边
            //    for (int i = 0; i < bFaceIndexThreeSideDetermined.Count; i++)
            //    {
            //        List<Point3d> realisticSiteFacePoints;
            //        List<int> realisticSiteFaceVertice;
            //        DDeepCopy = NakedEdge_3_Operation(DDeepCopy,
            //                                          subBoundarySegments,
            //                                          bFaceIndexThreeSideDetermined[i],
            //                                          publicHVThreeSideDetermined[i],
            //                                          out realisticSiteFaceVertice,
            //                                          out realisticSiteFacePoints,
            //                                          drawnFE,
            //                                          drawnDHV,
            //                                          abandonedDHV,
            //                                          vi_Point3d);

            //        realisticSiteFaceVerticesForAllFaceForGreenColor.Add(new List<int>());
            //        realisticSiteFacePointsForAllFaceForGreenColor.Add(new List<Point3d>());
            //        realisticSiteFaceVerticesForAllFaceForGreenColor[i].AddRange(realisticSiteFaceVertice);
            //        realisticSiteFacePointsForAllFaceForGreenColor[i].AddRange(realisticSiteFacePoints);

            //        realisticSiteFaceIndexForGreenColor.Add(bFaceIndexThreeSideDetermined[i]);
            //    }

            //    for (int i = 0; i < bFaceIndexTwoSidedDetermined.Count; i++)
            //    {
            //        List<Point3d> realisticSiteFacePoints;
            //        List<int> realisticSiteFaceVertice;
            //        DDeepCopy = NakedEdge_2_Operation(DDeepCopy,
            //                                         subBoundarySegments,
            //                                         bFaceIndexTwoSidedDetermined[i],
            //                                         publicHVTwoSideDetermined[i],
            //                                         out realisticSiteFaceVertice,
            //                                         out realisticSiteFacePoints,
            //                                         drawnFE,
            //                                         drawnDHV,
            //                                         vi_Point3d);

            //        realisticSiteFaceVerticesForAllFaceForGreenColor.Add(new List<int>());
            //        realisticSiteFacePointsForAllFaceForGreenColor.Add(new List<Point3d>());
            //        realisticSiteFaceVerticesForAllFaceForGreenColor[i + bFaceIndexThreeSideDetermined.Count].AddRange(realisticSiteFaceVertice);
            //        realisticSiteFacePointsForAllFaceForGreenColor[i + bFaceIndexThreeSideDetermined.Count].AddRange(realisticSiteFacePoints);

            //        realisticSiteFaceIndexForGreenColor.Add(bFaceIndexTwoSidedDetermined[i]);
            //    }
            //    #endregion

            //    #region 处理边缘面的邻接面
            //    Dictionary<int, List<int>> fIndex_SplitedIndex = new Dictionary<int, List<int>>();
            //    List<int[]> secondIndexPairs = new List<int[]>();
            //    for (int i = 0; i < aFaceIndex.Count; i++)
            //    {
            //        List<Point3d> realisticSiteFacePoints;
            //        List<int> realisticSiteFaceVertice;
            //        int[] secondIndexPair;
            //        // aFacePublicHVs[i][0]一定是跟每个边缘的有三边可以确定的面所对应的公共边
            //        DDeepCopy = AFace_Operation1(DDeepCopy,
            //                                     subBoundarySegments,
            //                                     aFaceIndex[i],
            //                                     aFacePublicHVs[i][0],
            //                                     abandonedDHV,
            //                                     out realisticSiteFaceVertice,
            //                                     out realisticSiteFacePoints,
            //                                     drawnFE,
            //                                     drawnDHV,
            //                                     vi_Point3d,
            //                                     out fIndex_SplitedIndex,
            //                                     out secondIndexPair);

            //        secondIndexPairs.Add(secondIndexPair);

            //        realisticSiteFaceVerticesForAllFaceForYellowColor.Add(new List<int>());
            //        realisticSiteFacePointsForAllFaceForYellowColor.Add(new List<Point3d>());
            //        realisticSiteFaceVerticesForAllFaceForYellowColor[i].AddRange(realisticSiteFaceVertice);
            //        realisticSiteFacePointsForAllFaceForYellowColor[i].AddRange(realisticSiteFacePoints);

            //        realisticSiteFaceIndexForYellowColor.Add(aFaceIndex[i]);
            //    }

            //    List<string> debug3 = UtilityFunctions.PrintFacesHalfedges(DDeepCopy);
            //    List<string> debug4 = UtilityFunctions.PrintHalfedgeStartAndEnd(DDeepCopy);
            //    List<string> debug5 = UtilityFunctions.PrintFacesVertices(DDeepCopy);

            //    for (int i = 0; i < aFaceIndex.Count; i++)
            //    {
            //        //if (i == 1)
            //        //{
            //        //    realisticSiteFaceVerticesForAllFaceForYellowColor.Add(new List<int>());
            //        //    realisticSiteFacePointsForAllFaceForYellowColor.Add(new List<Point3d>());
            //        //    continue;
            //        //}

            //        List<Point3d> realisticSiteFacePoints;
            //        List<int> realisticSiteFaceVertice;

            //        DDeepCopy = AFace_Operation2(DDeepCopy,
            //                                     subBoundarySegments,
            //                                     aFaceIndex[i],
            //                                     secondIndexPairs[i],
            //                                     out realisticSiteFaceVertice,
            //                                     out realisticSiteFacePoints,
            //                                     drawnFE,
            //                                     drawnDHV,
            //                                     abandonedDHV,
            //                                     vi_Point3d,
            //                                     fIndex_SplitedIndex);

            //        realisticSiteFaceVerticesForAllFaceForYellowColor.Add(new List<int>());
            //        realisticSiteFacePointsForAllFaceForYellowColor.Add(new List<Point3d>());
            //        realisticSiteFaceVerticesForAllFaceForYellowColor[aFaceIndex.Count + i].AddRange(realisticSiteFaceVertice);
            //        realisticSiteFacePointsForAllFaceForYellowColor[aFaceIndex.Count + i].AddRange(realisticSiteFacePoints);

            //        realisticSiteFaceIndexForYellowColor.Add(aFaceIndex[i]);


            //    }

            //    List<string> debug7 = UtilityFunctions.PrintFacesHalfedges(DDeepCopy);
            //    List<string> debug8 = UtilityFunctions.PrintHalfedgeStartAndEnd(DDeepCopy);
            //    List<string> debug9 = UtilityFunctions.PrintFacesVertices(DDeepCopy);



            //    #endregion


            //    #region 结果输出部分
            //    DebugPolylinesForGreenColor.Clear();
            //    DebugPolylinesForYellowColor.Clear();
            //    DebugPolylinesForGreenColorHatch.Clear();
            //    DebugPolylinesForYellowColorHatch.Clear();
            //    RealisticSiteFaceCenter.Clear();

            //    // 结果输出部分
            //    List<Polyline> debugPolylinesForGreenColor = new List<Polyline>();
            //    for (int i = 0; i < realisticSiteFacePointsForAllFaceForGreenColor.Count; i++)
            //    {
            //        List<Point3d> list = new List<Point3d>();
            //        list.AddRange(realisticSiteFacePointsForAllFaceForGreenColor[i]);
            //        list.Add(realisticSiteFacePointsForAllFaceForGreenColor[i][0]);
            //        Polyline polyline = new Polyline(list);
            //        debugPolylinesForGreenColor.Add(polyline);
            //    }
            //    DebugPolylinesForGreenColor.AddRange(debugPolylinesForGreenColor);

            //    List<Polyline> debugPolylinesForYellowColor = new List<Polyline>();
            //    for (int i = 0; i < realisticSiteFacePointsForAllFaceForYellowColor.Count; i++)
            //    {
            //        //if (realisticSiteFacePointsForAllFaceForYellowColor[i].Count == 0)
            //        //{
            //        //    continue;
            //        //}

            //        List<Point3d> list = new List<Point3d>();
            //        list.AddRange(realisticSiteFacePointsForAllFaceForYellowColor[i]);
            //        list.Add(realisticSiteFacePointsForAllFaceForYellowColor[i][0]);
            //        Polyline polyline = new Polyline(list);
            //        debugPolylinesForYellowColor.Add(polyline);
            //    }
            //    DebugPolylinesForYellowColor.AddRange(debugPolylinesForYellowColor);


            //    List<Polyline> allPolyline = new List<Polyline>();
            //    allPolyline.AddRange(debugPolylinesForGreenColor);
            //    allPolyline.AddRange(debugPolylinesForYellowColor);
            //    DA.SetDataList("DeBug Polyline", allPolyline);


            //    // Hatch
            //    List<Hatch> greenHatches = new List<Hatch>();
            //    for (int i = 0; i < debugPolylinesForGreenColor.Count; i++)
            //    {
            //        Hatch[] hatch = Hatch.Create(debugPolylinesForGreenColor[i].ToNurbsCurve(), 1, 1, 5, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            //        greenHatches.Add(hatch[0]);
            //    }
            //    DebugPolylinesForGreenColorHatch = greenHatches;
            //    List<Hatch> yellowHatches = new List<Hatch>();
            //    for (int i = 0; i < debugPolylinesForYellowColor.Count; i++)
            //    {
            //        Hatch[] hatch = Hatch.Create(debugPolylinesForYellowColor[i].ToNurbsCurve(), 1, 1, 5, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            //        yellowHatches.Add(hatch[0]);
            //    }
            //    DebugPolylinesForYellowColorHatch = yellowHatches;





            //    List<int> allFaceIndex = new List<int>();
            //    allFaceIndex.AddRange(realisticSiteFaceIndexForGreenColor);
            //    allFaceIndex.AddRange(realisticSiteFaceIndexForYellowColor);

            //    for (int i = 0; i < allPolyline.Count; i++)
            //    {
            //        Point3d center = allPolyline[i].CenterPoint();
            //        TextDot fCenterTextDot = new TextDot(allFaceIndex[i].ToString(), center);
            //        RealisticSiteFaceCenter.Add(fCenterTextDot);
            //    }
            //    #endregion
            //}
        }

        //private List<int> GetNakedEdges(PlanktonMesh D, int currentFaceIndex)
        //{
        //    List<int> nakedEdges = new List<int>();
        //    for (int i = 0; i < D.Faces.GetHalfedges(currentFaceIndex).Length; i++)
        //    {
        //        int pairH = D.Halfedges.GetPairHalfedge(D.Faces.GetHalfedges(currentFaceIndex)[i]);
        //        if (D.Halfedges[pairH].AdjacentFace == -1)
        //        {
        //            nakedEdges.Add(D.Faces.GetHalfedges(currentFaceIndex)[i]);
        //        }
        //    }
        //    return nakedEdges;
        //}

        //private bool IsNakedEdgesContinuous(PlanktonMesh D, List<int> nakedEdges)
        //{
        //    List<int> startVertice = new List<int>();
        //    List<int> endVertice = new List<int>();
        //    int nakedEdgesCount = nakedEdges.Count;

        //    for (int i = 0; i < nakedEdges.Count; i++)
        //    {
        //        int start = D.Halfedges[nakedEdges[i]].StartVertex;
        //        int end = D.Halfedges.EndVertex(nakedEdges[i]);

        //        startVertice.Add(start);
        //        endVertice.Add(end);
        //    }

        //    List<int> intersect = startVertice.Intersect(endVertice).ToList();
        //    if (intersect.Count == nakedEdgesCount - 1)
        //    {
        //        return true;
        //    }
        //    else
        //    {
        //        return false;
        //    }
        //}

        //private PlanktonMesh NakedEdge_3_Operation(PlanktonMesh D,
        //                                           List<List<FaceEdgeSegment>> subBoundarySegments,
        //                                           int currentBFaceIndexThreeSideDetermined,
        //                                           int[] publicHV,
        //                                           out List<int> realisticSiteFaceVertice,
        //                                           out List<Point3d> realisticSiteFacePoints,
        //                                           List<FaceEdgeSegment> drawnFE,
        //                                           List<int[]> drawnDHV,
        //                                           List<int[]> abandonedDHV,
        //                                           Dictionary<int, Point3d> vi_Point3d)
        //{
        //    #region 根据publicHV找到publicHalfedge
        //    int publicHalfedge = -1;
        //    int[] hs = D.Faces.GetHalfedges(currentBFaceIndexThreeSideDetermined);
        //    for (int i = 0; i < D.Halfedges.Count; i++)
        //    {
        //        int[] hV = new int[2] { D.Halfedges[i].StartVertex, D.Halfedges.EndVertex(i) };
        //        if (publicHV.Except(hV).ToList().Count == 0
        //            && hs.Contains(i))
        //        {
        //            publicHalfedge = i;
        //            break;
        //        }
        //    }
        //    #endregion

        //    // adjacentFace的序号
        //    int adjacentFaceIndex = D.Halfedges[D.Halfedges.GetPairHalfedge(publicHalfedge)].AdjacentFace;

        //    #region 获取目前对应的三个relatedSubBoundarySegments，进而得到目前的四个边的向量vectorsForCurrentFace
        //    #region 从公共边的nextHalfedge开始，依次添加这个面的所有Halfedge
        //    // 按顺序存储这个面的currentFaceHIndexList
        //    List<int> currentFaceHIndexList = new List<int>();
        //    List<int[]> currentFaceHVList = new List<int[]>();
        //    int iter = 0;
        //    int currentEdge = publicHalfedge;
        //    do
        //    {
        //        currentFaceHIndexList.Add(D.Halfedges[currentEdge].NextHalfedge);
        //        int[] faceHV = new int[2] { D.Halfedges[D.Halfedges[currentEdge].NextHalfedge].StartVertex, D.Halfedges.EndVertex(D.Halfedges[currentEdge].NextHalfedge) };
        //        currentFaceHVList.Add(faceHV);
        //        currentEdge = currentFaceHIndexList.Last();
        //        iter++;
        //    } while (iter < D.Faces.GetHalfedges(currentBFaceIndexThreeSideDetermined).Length);
        //    #endregion

        //    #region 按照currentFaceHIndexList中的顺序，添加现有的relatedSubBoundarySegments
        //    // 按顺序存储这个面目前有的BoundarySegment
        //    // 注意，此时公共边没有对应的FaceEdgeSegment
        //    List<FaceEdgeSegment> relatedSubBoundarySegments = new List<FaceEdgeSegment>();
        //    for (int i = 0; i < currentFaceHVList.Count - 1; i++)
        //    {
        //        for (int j = 0; j < subBoundarySegments.Count; j++)
        //        {
        //            for (int k = 0; k < subBoundarySegments[j].Count; k++)
        //            {
        //                if (subBoundarySegments[j][k].IncludedDVertice.Except(currentFaceHVList[i]).ToList().Count == 0)
        //                {
        //                    relatedSubBoundarySegments.Add(subBoundarySegments[j][k]);
        //                }
        //            }
        //        }
        //    }
        //    #endregion

        //    #region 按照currentFaceHIndexList中的顺序（即现有的relatedSubBoundarySegments的顺序），添加Vertex对应的Point位置
        //    // 注意，此时的relatedSubBoundarySegmentPoints并不是最终的结果
        //    // 注意，segment方向与半边方向是反的，所以先添加To，然后全部添加From
        //    List<Point3d> relatedSubBoundarySegmentPoints = new List<Point3d>();
        //    relatedSubBoundarySegmentPoints.Add(relatedSubBoundarySegments[0].To);
        //    for (int i = 0; i < relatedSubBoundarySegments.Count; i++)
        //    {
        //        relatedSubBoundarySegmentPoints.Add(relatedSubBoundarySegments[i].From);
        //    }
        //    #endregion

        //    #region 根据relatedSubBoundarySegmentPoints，来构造代表每条边的vector
        //    List<Vector3d> vectorsForCurrentFace = new List<Vector3d>();
        //    for (int i = 0; i < relatedSubBoundarySegmentPoints.Count - 1; i++)
        //    {
        //        Vector3d vector = new Vector3d(relatedSubBoundarySegmentPoints[i + 1] - relatedSubBoundarySegmentPoints[i]);
        //        vectorsForCurrentFace.Add(vector);
        //    }
        //    vectorsForCurrentFace.Add(new Vector3d(relatedSubBoundarySegmentPoints[0] - relatedSubBoundarySegmentPoints.Last()));
        //    #endregion
        //    #endregion

        //    #region 判断是否需要生成新的顶点
        //    List<bool> flagList = new List<bool>();
        //    for (int i = 0; i < relatedSubBoundarySegmentPoints.Count; i++)
        //    {
        //        Vector3d prevV = vectorsForCurrentFace[((i - 1) + vectorsForCurrentFace.Count) % vectorsForCurrentFace.Count];
        //        Vector3d nextV = vectorsForCurrentFace[i];
        //        Vector3d crossProduct = Vector3d.CrossProduct(prevV, nextV);

        //        if (crossProduct.Z < 0)
        //        {
        //            flagList.Add(false);
        //        }
        //        else
        //        {
        //            flagList.Add(true);
        //        }
        //    }

        //    int index = -1;
        //    bool needToModifyD = false;
        //    if (flagList.Contains(true))
        //    {
        //        needToModifyD = true;
        //        index = flagList.IndexOf(true);

        //        // 添加废弃半边
        //        int[] abandoned = new int[2] { D.Halfedges[currentFaceHIndexList[index]].StartVertex, D.Halfedges.EndVertex(currentFaceHIndexList[index]) };
        //        abandonedDHV.Add(abandoned);
        //    }
        //    #endregion

        //    #region 把目前正确的realisticSiteFaceVertice和realisticSiteFacePoints的值都存储好
        //    // 目前realisticSiteFaceVertice有四个，realisticSiteFacePoints有四个
        //    realisticSiteFaceVertice = new List<int>();
        //    for (int i = 0; i < currentFaceHIndexList.Count; i++)
        //    {
        //        realisticSiteFaceVertice.Add(D.Halfedges[currentFaceHIndexList[i]].StartVertex);
        //    }

        //    realisticSiteFacePoints = new List<Point3d>();
        //    realisticSiteFacePoints.AddRange(relatedSubBoundarySegmentPoints);
        //    // realisticSiteFaceVertice和realisticSiteFacePoints，在需要修改D时，需要对其中某个值做替换

        //    for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //    {
        //        vi_Point3d.Add(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //    }
        //    #endregion

        //    int count = realisticSiteFaceVertice.Count;
        //    PlanktonMesh newD = new PlanktonMesh();
        //    if (needToModifyD)
        //    {
        //        /* 需要修改D时 */
        //        // 修改D
        //        newD = RegeneratePMesh_3_Operation(D, currentBFaceIndexThreeSideDetermined, adjacentFaceIndex, index, currentFaceHIndexList);

        //        #region 修改替换正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //        // 要修改index点的nextPoint
        //        Point3d prevPrevPoint = relatedSubBoundarySegmentPoints[((index - 2) + count) % count];
        //        Point3d newNextPoint = prevPrevPoint + vectorsForCurrentFace[((index - 1) + count) % count];

        //        realisticSiteFaceVertice[(index + 1) % count] = newD.Vertices.Count - 1;
        //        realisticSiteFacePoints[(index + 1) % count] = newNextPoint;

        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //        }
        //        #endregion

        //        #region 生成正确的nextFE和nextNextFE
        //        List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //            string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            List<int> includeDV = new List<int>();
        //            includeDV.Add(realisticSiteFaceVertice[i]);
        //            includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //            realisticSiteFaceEdgeSegments.Add(FE);
        //        }

        //        //// 注意，segment方向与半边方向是反的
        //        //Line newNextLine = new Line(newNextPoint, realisticSiteFacePoints[index]);
        //        //string label1 = string.Format("v:{0},{1}", realisticSiteFaceVertice[index], realisticSiteFaceVertice[(index + 1) % count]);
        //        //List<int> includeDV1 = new List<int>();
        //        //includeDV1.Add(realisticSiteFaceVertice[index]);
        //        //includeDV1.Add(newD.Vertices.Count - 1);
        //        //FaceEdgeSegment newNextFE = new FaceEdgeSegment(newNextLine, label1, includeDV1, -1);

        //        //// 实际上此时realisticSiteFaceVertice.Last()与realisticSiteFacePoints[index + 1]等效
        //        //// 0与index+2等效
        //        //// 注意，segment方向与半边方向是反的
        //        //Line newNextNextLine = new Line(realisticSiteFacePoints[(index + 2) % count], newNextPoint);
        //        //string label2 = string.Format("v:{0},{1}", newD.Vertices.Count - 1, realisticSiteFaceVertice[(index + 2) % count]);
        //        //List<int> includeDV2 = new List<int>();
        //        //includeDV2.Add(newD.Vertices.Count - 1);
        //        //includeDV2.Add(realisticSiteFaceVertice[(index + 2) % count]);
        //        //FaceEdgeSegment newNextNextFE = new FaceEdgeSegment(newNextNextLine, label2, includeDV2, -1);
        //        #endregion
        //        //// hIndex已经没用了，暂时先设为-1了

        //        #region 其他需要更新的参数
        //        ///* 更新relatedSubBoundarySegments，输出最终正确的realisticSiteFaceEdgeSegments */
        //        //relatedSubBoundarySegments[index] = newNextFE;
        //        //List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        //realisticSiteFaceEdgeSegments.AddRange(relatedSubBoundarySegments);
        //        //realisticSiteFaceEdgeSegments.Add(newNextNextFE);

        //        /* 输出最终正确的drawnDHV */
        //        for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //        {
        //            int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //            drawnDHV.Add(dHV);
        //        }
        //        /* 输出最终正确的feHasH */
        //        drawnFE.AddRange(realisticSiteFaceEdgeSegments);
        //        #endregion
        //    }
        //    else
        //    {
        //        /* 不需要修改D时 */
        //        // 不需要修改D
        //        newD = new PlanktonMesh(D);
        //        // 不需要修改替换正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //        }

        //        #region 生成正确的newFE

        //        List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //            string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            List<int> includeDV = new List<int>();
        //            includeDV.Add(realisticSiteFaceVertice[i]);
        //            includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //            realisticSiteFaceEdgeSegments.Add(FE);
        //        }


        //        //// 注意，segment方向与半边方向是反的
        //        //Line newLine = new Line(realisticSiteFacePoints[0], realisticSiteFacePoints.Last());
        //        //string label = string.Format("h:{0},v:{1},{2}", currentFaceHIndexList.Last(), realisticSiteFaceVertice.Last(), realisticSiteFaceVertice[0]);
        //        //List<int> includeDV = new List<int>();
        //        //includeDV.Add(realisticSiteFaceVertice.Last());
        //        //includeDV.Add(realisticSiteFaceVertice[0]);
        //        //FaceEdgeSegment newFE = new FaceEdgeSegment(newLine, label, includeDV, -1);
        //        #endregion
        //        //// hIndex已经没用了，暂时先设为-1了

        //        #region 其他需要更新的参数
        //        ///* 更新relatedSubBoundarySegments，输出最终正确的realisticSiteFaceEdgeSegments */
        //        //List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        //realisticSiteFaceEdgeSegments.AddRange(relatedSubBoundarySegments);
        //        //realisticSiteFaceEdgeSegments.Add(newFE);

        //        /* 输出最终正确的drawnDHV */
        //        for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //        {
        //            int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //            drawnDHV.Add(dHV);
        //        }
        //        /* 输出最终正确的feHasH */
        //        drawnFE.AddRange(realisticSiteFaceEdgeSegments);
        //        #endregion
        //    }

        //    return newD;
        //}

        //private PlanktonMesh NakedEdge_2_Operation(PlanktonMesh D,
        //                                          List<List<FaceEdgeSegment>> subBoundarySegments,
        //                                          int currentBFaceIndexTwoSideDetermined,
        //                                          List<int[]> publicHVs,
        //                                          out List<int> realisticSiteFaceVertice,
        //                                          out List<Point3d> realisticSiteFacePoints,
        //                                          List<FaceEdgeSegment> drawnFE,
        //                                          List<int[]> drawnDHV,
        //                                          Dictionary<int, Point3d> vi_Point3d)
        //{
        //    #region 根据publicHVs找到publicHalfedges
        //    List<int> publicHalfedges = new List<int>();
        //    int[] hs = D.Faces.GetHalfedges(currentBFaceIndexTwoSideDetermined);
        //    for (int i = 0; i < publicHVs.Count; i++)
        //    {
        //        for (int j = 0; j < D.Halfedges.Count; j++)
        //        {
        //            int[] hV = new int[2] { D.Halfedges[j].StartVertex, D.Halfedges.EndVertex(j) };
        //            if (publicHVs[i].Except(hV).ToList().Count == 0
        //                && hs.Contains(j))
        //            {
        //                publicHalfedges.Add(j);
        //                break;
        //            }
        //        }
        //    }
        //    #endregion

        //    // adjacentFace的序号
        //    List<int> adjacentFaceIndexs = new List<int>();
        //    for (int i = 0; i < publicHalfedges.Count; i++)
        //    {
        //        adjacentFaceIndexs.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(publicHalfedges[i])].AdjacentFace);
        //    }

        //    #region 从最后一个公共边的nextHalfedge开始，依次添加这个面的所有Halfedge
        //    // 按顺序存储这个面的halfedge
        //    List<int> currentFaceHIndexList = new List<int>();
        //    List<int[]> currentFaceHVList = new List<int[]>();
        //    int iter = 0;
        //    int currentEdge;
        //    // 找到最后一个公共边的index，作为currentEdge
        //    if (D.Halfedges[publicHalfedges[0]].StartVertex == D.Halfedges.EndVertex(publicHalfedges[1]))
        //    {
        //        currentEdge = publicHalfedges[0];
        //    }
        //    else
        //    {
        //        currentEdge = publicHalfedges[1];
        //    }

        //    do
        //    {
        //        currentFaceHIndexList.Add(D.Halfedges[currentEdge].NextHalfedge);
        //        int[] faceHV = new int[2] { D.Halfedges[D.Halfedges[currentEdge].NextHalfedge].StartVertex, D.Halfedges.EndVertex(D.Halfedges[currentEdge].NextHalfedge) };
        //        currentFaceHVList.Add(faceHV);
        //        currentEdge = currentFaceHIndexList.Last();
        //        iter++;
        //    } while (iter < D.Faces.GetHalfedges(currentBFaceIndexTwoSideDetermined).Length);
        //    #endregion

        //    #region 按照currentFaceHIndexList中的顺序，添加现有的relatedSubBoundarySegments
        //    // 按顺序存储这个面目前有的BoundarySegment
        //    // 注意，此时公共边没有对应的FaceEdgeSegment
        //    List<FaceEdgeSegment> relatedSubBoundarySegments = new List<FaceEdgeSegment>();
        //    for (int i = 0; i < currentFaceHIndexList.Count - 1; i++)
        //    {
        //        for (int j = 0; j < subBoundarySegments.Count; j++)
        //        {
        //            for (int k = 0; k < subBoundarySegments[j].Count; k++)
        //            {
        //                if (subBoundarySegments[j][k].IncludedDVertice.Except(currentFaceHVList[i]).ToList().Count == 0)
        //                {
        //                    relatedSubBoundarySegments.Add(subBoundarySegments[j][k]);
        //                }
        //            }
        //        }
        //    }
        //    #endregion

        //    #region 按照currentHIndexList中的顺序（即现有的relatedSubBoundarySegments的顺序），添加Vertex对应的Point位置
        //    // 注意，此时的relatedSubBoundarySegmentPoints并不是最终的结果
        //    // 注意，segment方向与半边方向是反的，所以先添加了To，然后全部添加From
        //    List<Point3d> relatedSubBoundarySegmentPoints = new List<Point3d>();
        //    relatedSubBoundarySegmentPoints.Add(relatedSubBoundarySegments[0].To);
        //    for (int i = 0; i < relatedSubBoundarySegments.Count; i++)
        //    {
        //        relatedSubBoundarySegmentPoints.Add(relatedSubBoundarySegments[i].From);
        //    }
        //    #endregion

        //    #region 判断哪条公共边已经绘制过了
        //    // index实际上是2
        //    int index = relatedSubBoundarySegments.Count;

        //    // 目前relatedSubBoundarySegments只有两个元素。currentFaceHIndexList[relatedSubBoundarySegments.Count]是它的下一条边
        //    int[] indexPair0 = new int[2] { D.Halfedges[currentFaceHIndexList[index]].StartVertex, D.Halfedges.EndVertex(currentFaceHIndexList[index]) };
        //    int[] indexPair1 = new int[2] { D.Halfedges[currentFaceHIndexList[index + 1]].StartVertex, D.Halfedges.EndVertex(currentFaceHIndexList[index + 1]) };
        //    bool flag0 = false;
        //    bool flag1 = false;
        //    FaceEdgeSegment alreadyExistFE0 = null;
        //    FaceEdgeSegment alreadyExistFE1 = null;
        //    for (int i = 0; i < drawnDHV.Count; i++)
        //    {
        //        if (drawnDHV[i].Except(indexPair0).ToList().Count == 0)
        //        {
        //            flag0 = true;
        //            alreadyExistFE0 = drawnFE[i];
        //            continue;
        //        }
        //        if (drawnDHV[i].Except(indexPair1).ToList().Count == 0)
        //        {
        //            flag1 = true;
        //            alreadyExistFE1 = drawnFE[i];
        //            continue;
        //        }
        //    }
        //    #endregion

        //    #region 把目前正确的realisticSiteFaceVertice和realisticSiteFacePoints的值都存储好
        //    // 目前realisticSiteFaceVertice有四个，realisticSiteFacePoints有三个
        //    realisticSiteFaceVertice = new List<int>();
        //    for (int i = 0; i < currentFaceHIndexList.Count; i++)
        //    {
        //        realisticSiteFaceVertice.Add(D.Halfedges[currentFaceHIndexList[i]].StartVertex);
        //    }

        //    realisticSiteFacePoints = new List<Point3d>();
        //    realisticSiteFacePoints.AddRange(relatedSubBoundarySegmentPoints);
        //    // realisticSiteFaceVertice和realisticSiteFacePoints，在需要修改D时，需要对其中某个值做替换

        //    for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //    {
        //        vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //    }
        //    #endregion

        //    PlanktonMesh newD = new PlanktonMesh();
        //    int count = realisticSiteFaceVertice.Count;
        //    if (!flag0 && !flag1)
        //    {
        //        /* 如果两条公共边的对边都没有在hHasCorrespondFE中出现过的话（表示这两条公共边都没有绘制过对应的FaceEdgeSegment，即这两条半边的对边都没有生成过对应的FE） */
        //        // 不用修改D
        //        newD = new PlanktonMesh(D);

        //        #region realisticSiteFaceVertice不用修改，添加正确的realisticSiteFacePoints
        //        Point3d currPoint = realisticSiteFacePoints[index];
        //        Vector3d move = new Vector3d(realisticSiteFacePoints[((index - 2) + count) % count] - realisticSiteFacePoints[((index - 1) + count) % count]);
        //        Point3d nextPoint = currPoint + move;

        //        // realisticSiteFaceVertice不用修改
        //        realisticSiteFacePoints.Add(nextPoint);

        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //        }
        //        #endregion

        //        #region 生成正确的nextFE和nextNextFE

        //        List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //            string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            List<int> includeDV = new List<int>();
        //            includeDV.Add(realisticSiteFaceVertice[i]);
        //            includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //            realisticSiteFaceEdgeSegments.Add(FE);
        //        }

        //        //Line nextLine = new Line(realisticSiteFacePoints[(index + 1) % count], realisticSiteFacePoints[index]);
        //        //string label1 = string.Format("v:{0},{1}", realisticSiteFaceVertice[index], realisticSiteFaceVertice[(index + 1) % count]);
        //        //List<int> includedDV1 = new List<int>();
        //        //includedDV1.Add(realisticSiteFaceVertice[index]);
        //        //includedDV1.Add(realisticSiteFaceVertice[(index + 1) % count]);
        //        //FaceEdgeSegment nextFE = new FaceEdgeSegment(nextLine, label1, includedDV1, -1);

        //        //Line nextNextLine = new Line(realisticSiteFacePoints[(index + 2) % count], realisticSiteFacePoints[(index + 1) % count]);
        //        //string label2 = string.Format("v:{0},{1}", realisticSiteFaceVertice[(index + 1) % count], realisticSiteFaceVertice[(index + 2) % count]);
        //        //List<int> includedDV2 = new List<int>();
        //        //includedDV2.Add(realisticSiteFaceVertice[realisticSiteFaceVertice[(index + 1) % count]]);
        //        //includedDV2.Add(realisticSiteFaceVertice[(index + 2) % count]);
        //        //FaceEdgeSegment nextNextFE = new FaceEdgeSegment(nextNextLine, label2, includedDV2, -1);
        //        #endregion
        //        // hIndex已经没用了，暂时先设为-1了

        //        #region 其他需要更新的参数
        //        ///* 更新relatedSubBoundarySegments，输出最终正确的realisticSiteFaceEdgeSegments */
        //        //relatedSubBoundarySegments.Add(nextFE);
        //        //relatedSubBoundarySegments.Add(nextNextFE);
        //        //List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        //realisticSiteFaceEdgeSegments.AddRange(relatedSubBoundarySegments);

        //        /* 输出最终正确的drawnDHV */
        //        for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //        {
        //            int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //            drawnDHV.Add(dHV);
        //        }
        //        /* 输出最终正确的feHasH */
        //        drawnFE.AddRange(realisticSiteFaceEdgeSegments);
        //        #endregion
        //    }
        //    else if ((flag0 && !flag1) || (!flag0 && flag1))
        //    {
        //        /* 如果两条中有一条出现，另一条没出现，那么求交点 */

        //        // 需要修改D
        //        List<int> sortedAdjacentFace = new List<int>();
        //        List<int> sortedPublicHalfedge = new List<int>();

        //        #region 确定正确的PublicHalfedge顺序和AdjacentFace顺序

        //        FaceEdgeSegment alreadyExistFE;
        //        FaceEdgeSegment feToMove;
        //        Line lineToMove;
        //        int vIndex1;
        //        Point3d point1;
        //        int vIndex2;
        //        Point3d point2;
        //        Transform move;
        //        if (flag0 && !flag1)
        //        {
        //            int hIndex = D.Halfedges.GetPairHalfedge(currentFaceHIndexList[index]);
        //            //alreadyExistFE = feHasH[hHasCorrespondFE.IndexOf(hIndex)];
        //            alreadyExistFE = alreadyExistFE0;

        //            feToMove = relatedSubBoundarySegments[1];
        //            lineToMove = feToMove.Line;
        //            vIndex1 = realisticSiteFaceVertice[2];
        //            point1 = realisticSiteFacePoints[2];
        //            vIndex2 = realisticSiteFaceVertice[0];
        //            point2 = realisticSiteFacePoints[0];
        //            move = Transform.Translation(new Vector3d(realisticSiteFacePoints[0] - realisticSiteFacePoints[1]));

        //            // 保证sortedAdjacentFace中第一个adjacentFace是那个已经绘制过segment的adjacentFace
        //            sortedAdjacentFace.Add(D.Halfedges[hIndex].AdjacentFace);
        //            sortedAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(currentFaceHIndexList[index + 1])].AdjacentFace);

        //            sortedPublicHalfedge.Add(currentFaceHIndexList[index]);
        //            sortedPublicHalfedge.Add(currentFaceHIndexList[index + 1]);
        //        }
        //        else
        //        {
        //            int hIndex = D.Halfedges.GetPairHalfedge(currentFaceHIndexList[index + 1]);
        //            //alreadyExistFE = feHasH[hHasCorrespondFE.IndexOf(hIndex)];
        //            alreadyExistFE = alreadyExistFE1;

        //            feToMove = relatedSubBoundarySegments[0];
        //            lineToMove = feToMove.Line;
        //            vIndex1 = realisticSiteFaceVertice[0];
        //            point1 = realisticSiteFacePoints[0];
        //            vIndex2 = realisticSiteFaceVertice[2];
        //            point2 = realisticSiteFacePoints[2];
        //            move = Transform.Translation(new Vector3d(realisticSiteFacePoints[2] - realisticSiteFacePoints[1]));

        //            // 保证sortedAdjacentFace中第一个adjacentFace是那个已经绘制过segment的adjacentFace
        //            sortedAdjacentFace.Add(D.Halfedges[hIndex].AdjacentFace);
        //            sortedAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(currentFaceHIndexList[index])].AdjacentFace);

        //            sortedPublicHalfedge.Add(currentFaceHIndexList[index + 1]);
        //            sortedPublicHalfedge.Add(currentFaceHIndexList[index]);
        //        }
        //        #endregion


        //        #region 求交点
        //        Line newLine = new Line(lineToMove.From, lineToMove.To);
        //        newLine.Transform(move);

        //        bool isOutOfRange = false;
        //        double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
        //        double t1;
        //        double t2;
        //        if (Intersection.LineLine(alreadyExistFE.Line, newLine, out t1, out t2, tolerance, false))
        //        {
        //            if (t1 > 1 || t1 < 0)
        //            {
        //                isOutOfRange = true;
        //            }

        //            // 修改D
        //            newD = RegeneratePMesh_2_Operation(D, currentBFaceIndexTwoSideDetermined, sortedAdjacentFace, sortedPublicHalfedge, isOutOfRange, flag0, flag1);
        //            //List<string> debug1 = UtilityFunctions.PrintFacesHalfedges(newD);
        //            //List<string> debug2 = UtilityFunctions.PrintFacesVertices(newD);
        //            //List<string> debug3 = UtilityFunctions.PrintHalfedgeStartAndEnd(newD);

        //            #region 修改替换正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //            Point3d intersectPoint1 = alreadyExistFE.Line.PointAt(t1);
        //            Point3d intersectPoint2 = newLine.PointAt(t2);

        //            realisticSiteFaceVertice[(index + 1) % count] = newD.Vertices.Count - 1;
        //            realisticSiteFacePoints.Add(intersectPoint1);

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //            }
        //            #endregion

        //            #region 生成正确的newFE1和newFE2

        //            List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                List<int> includeDV = new List<int>();
        //                includeDV.Add(realisticSiteFaceVertice[i]);
        //                includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                realisticSiteFaceEdgeSegments.Add(FE);
        //            }

        //            //Line lineCoincideWithAlreadyExistFE = new Line(point1, intersectPoint1);
        //            //string label1 = string.Format("v:{0},{1} OriginLabel:{2}", vIndex1, realisticSiteFaceVertice[(index + 1) % count], alreadyExistFE.Label);
        //            //List<int> includedDV1 = new List<int>();
        //            //includedDV1.Add(vIndex1);
        //            //includedDV1.Add(realisticSiteFaceVertice[(index + 1) % count]);

        //            //FaceEdgeSegment newFE1 = new FaceEdgeSegment(lineCoincideWithAlreadyExistFE, label1, includedDV1, -1);

        //            //Line lineForNewFE = new Line(intersectPoint2, point2);
        //            //string label2 = string.Format("v:{0},{1} OriginLabel:{2}", realisticSiteFaceVertice[(index + 1) % count], vIndex2, alreadyExistFE.Label);
        //            //List<int> includedDV2 = new List<int>();
        //            //includedDV2.Add(realisticSiteFaceVertice[(index + 1) % count]);
        //            //includedDV2.Add(vIndex2);

        //            //FaceEdgeSegment newFE2 = new FaceEdgeSegment(lineForNewFE, label2, includedDV2, -1);
        //            #endregion
        //            // hIndex已经没用了，暂时先设为-1了

        //            #region 其他需要更新的参数
        //            ///* 更新relatedSubBoundarySegments，输出最终正确的realisticSiteFaceEdgeSegments */
        //            //List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //            //realisticSiteFaceEdgeSegments.AddRange(relatedSubBoundarySegments);
        //            //realisticSiteFaceEdgeSegments.Add(newFE1);
        //            //realisticSiteFaceEdgeSegments.Add(newFE2);

        //            /* 输出最终正确的drawnDHV */
        //            for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //            {
        //                int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //                drawnDHV.Add(dHV);
        //            }
        //            /* 输出最终正确的feHasH */
        //            drawnFE.AddRange(realisticSiteFaceEdgeSegments);
        //            #endregion

        //        }
        //        #endregion
        //    }
        //    else
        //    {
        //        /* 如果都出现了，那么直接找两个FE的公共点即可 */

        //        // 不需要修改D
        //        newD = new PlanktonMesh(D);

        //        #region 修改替换正确的realisticSiteFaceVertice和realistcSiteFacePoints
        //        Point3d commonPoint;
        //        if (alreadyExistFE0.From == alreadyExistFE1.To)
        //        {
        //            commonPoint = alreadyExistFE0.From;
        //            realisticSiteFaceVertice[(index + 1) % count] = alreadyExistFE0.IncludedDVertice[1];
        //            realisticSiteFacePoints.Add(commonPoint);
        //        }
        //        if (alreadyExistFE0.To == alreadyExistFE1.From)
        //        {
        //            commonPoint = alreadyExistFE0.To;
        //            realisticSiteFaceVertice[(index + 1) % count] = alreadyExistFE0.IncludedDVertice[0];
        //            realisticSiteFacePoints.Add(commonPoint);
        //        }

        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //        }
        //        #endregion

        //        #region 生成正确的newFE1和newFE2

        //        List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //        {
        //            Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //            string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            List<int> includeDV = new List<int>();
        //            includeDV.Add(realisticSiteFaceVertice[i]);
        //            includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //            FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //            realisticSiteFaceEdgeSegments.Add(FE);
        //        }

        //        //FaceEdgeSegment newFE1 = new FaceEdgeSegment(alreadyExistFE0);
        //        //FaceEdgeSegment newFE2 = new FaceEdgeSegment(alreadyExistFE1);
        //        #endregion

        //        #region 其他需要更新的参数
        //        ///* 更新relatedSubBoundarySegments，输出最终正确的realisticSiteFaceEdgeSegments */
        //        //List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //        //realisticSiteFaceEdgeSegments.AddRange(relatedSubBoundarySegments);
        //        //realisticSiteFaceEdgeSegments.Add(newFE1);
        //        //realisticSiteFaceEdgeSegments.Add(newFE2);

        //        /* 输出最终正确的drawnDHV */
        //        for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //        {
        //            int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //            drawnDHV.Add(dHV);
        //        }
        //        /* 输出最终正确的feHasH */
        //        drawnFE.AddRange(realisticSiteFaceEdgeSegments);
        //        #endregion
        //    }

        //    return newD;
        //}



        //private PlanktonMesh AFace_Operation1(PlanktonMesh D,
        //                                      List<List<FaceEdgeSegment>> subBoundarySegments,
        //                                      int currentFaceIndex,
        //                                      int[] publicWithThreeNakedEdgesHV,
        //                                      List<int[]> abandonedDHV,
        //                                      out List<int> realisticSiteFaceVertice,
        //                                      out List<Point3d> realisticSiteFacePoints,
        //                                      List<FaceEdgeSegment> drawnFE,
        //                                      List<int[]> drawnDHV,
        //                                      Dictionary<int, Point3d> vi_Point3d,
        //                                      out Dictionary<int, List<int>> fIndex_SplitedIndex,
        //                                      out int[] secondIndexPair)
        //{
        //    #region 得到排序好的firstHalfedge
        //    List<int> arrangedFaceHalfedges = ArrangeFaceBoundaryHalfedges(D, currentFaceIndex, publicWithThreeNakedEdgesHV, abandonedDHV);
        //    int firstHalfedge = arrangedFaceHalfedges[0];
        //    int secondHalfedge = arrangedFaceHalfedges[1];
        //    #endregion

        //    #region 找到currH，prevH，nextH
        //    int currH = firstHalfedge;
        //    int prevH = D.Halfedges[firstHalfedge].PrevHalfedge;
        //    int nextH = D.Halfedges[firstHalfedge].NextHalfedge;
        //    int preStartVertex = D.Halfedges[prevH].StartVertex;
        //    int nextEndVertex = D.Halfedges.EndVertex(nextH);

        //    int[] prevIndexPair = new int[2] { D.Halfedges[prevH].StartVertex, D.Halfedges.EndVertex(prevH) };
        //    int[] nextIndexPair = new int[2] { D.Halfedges[nextH].StartVertex, D.Halfedges.EndVertex(nextH) };
        //    int[] currIndexPair = new int[2] { D.Halfedges[currH].StartVertex, D.Halfedges.EndVertex(currH) };
        //    #endregion

        //    #region 找到firstHalfedge和secondHalfedge各自对应的FaceEdgeSegment
        //    FaceEdgeSegment firstH_Segment = null;
        //    FaceEdgeSegment secondH_Segment = null;
        //    // int[] firstIndexPair = new int[2] { D.Halfedges[firstHalfedge].StartVertex, D.Halfedges.EndVertex(firstHalfedge) };
        //    secondIndexPair = new int[2] { D.Halfedges[secondHalfedge].StartVertex, D.Halfedges.EndVertex(secondHalfedge) };
        //    for (int i = 0; i < subBoundarySegments.Count; i++)
        //    {
        //        for (int j = 0; j < subBoundarySegments[i].Count; j++)
        //        {
        //            if (subBoundarySegments[i][j].IncludedDVertice.Except(currIndexPair).ToList().Count == 0)
        //            {
        //                firstH_Segment = subBoundarySegments[i][j];
        //                break;
        //            }
        //            if (subBoundarySegments[i][j].IncludedDVertice.Except(secondIndexPair).ToList().Count == 0)
        //            {
        //                secondH_Segment = subBoundarySegments[i][j];
        //                break;
        //            }
        //        }
        //    }
        //    #region 查找firstNextHalfedge对应的Segment
        //    // int firstNextHalfedge = D.Halfedges[firstHalfedge].NextHalfedge;
        //    // int[] firstNextHV = new int[2] { firstIndexPair[1], nextEndVertex };

        //    FaceEdgeSegment nextFE = null;
        //    for (int i = 0; i < drawnFE.Count; i++)
        //    {
        //        if (drawnFE[i].IncludedDVertice.Except(nextIndexPair).ToList().Count == 0)
        //        {
        //            nextFE = drawnFE[i];
        //        }
        //    }
        //    #endregion
        //    #endregion

        //    /* 这里需要进行判断，prevFE有没有绘制过 */
        //    // int[] firstPrevHV = new int[2] { D.Halfedges[D.Halfedges[firstHalfedge].PrevHalfedge].StartVertex, D.Halfedges[firstHalfedge].StartVertex };
        //    bool prevHasDrawnFlag = false;
        //    int drawnPrevH = -1;
        //    for (int i = 0; i < drawnDHV.Count; i++)
        //    {
        //        if (drawnDHV[i].Except(prevIndexPair).ToList().Count == 0)
        //        {
        //            prevHasDrawnFlag = true;
        //            drawnPrevH = i;
        //        }
        //    }

        //    fIndex_SplitedIndex = new Dictionary<int, List<int>>();
        //    PlanktonMesh newD = new PlanktonMesh();

        //    // 修改D
        //    List<int> sortedAdjacentFace = new List<int>();
        //    // List<int> sortedPublicHalfedge = new List<int>();
        //    // 确定正确的PublicHalfedge顺序和AdjacentFace顺序
        //    sortedAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(prevH)].AdjacentFace);
        //    sortedAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(nextH)].AdjacentFace);

        //    int prevMidVertex = -1;
        //    int nextMidVertex = -1;
        //    newD = RegeneratePMesh_AFace_Operation1(D,
        //                                            currentFaceIndex,
        //                                            sortedAdjacentFace,
        //                                            preStartVertex,
        //                                            currIndexPair[0],
        //                                            currIndexPair[1],
        //                                            nextEndVertex,
        //                                            fIndex_SplitedIndex,
        //                                            out prevMidVertex,
        //                                            out nextMidVertex);

        //    double t = firstH_Segment.Line.Length / (firstH_Segment.Line.Length + secondH_Segment.Line.Length);
        //    // 注意，segment方向与半边方向是反的，所以t的分子是firstH_Segment.Line.Length，而不是secondH_Segment.Line.Length
        //    Point3d pointOnNextHalfedge = nextFE.Line.PointAt(t);

        //    Vector3d moveVector = new Vector3d(pointOnNextHalfedge - firstH_Segment.From);

        //    Point3d pointOnPrevHalfedge;

        //    bool tBelow0 = false;
        //    bool tAbove1 = false;
        //    if (prevHasDrawnFlag)
        //    {
        //        FaceEdgeSegment alreadyExistFE = drawnFE[drawnPrevH];

        //        Line lineToMove = firstH_Segment.Line;
        //        Transform move = Transform.Translation(moveVector);
        //        Line newLine = new Line(lineToMove.From, lineToMove.To);
        //        newLine.Transform(move);

        //        double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
        //        double t1;
        //        double t2;
                
        //        // 求交点
        //        if (Intersection.LineLine(alreadyExistFE.Line, newLine, out t1, out t2, tolerance, false))
        //        {
        //            double t1Round = Math.Round(t1, 2);
        //            if (Math.Round(t1, 2) < 0)
        //            {
        //                tBelow0 = true;
        //            }
        //            else if (Math.Round(t1, 2) > 1)
        //            {
        //                tAbove1 = true;
        //            }

        //            pointOnPrevHalfedge = alreadyExistFE.Line.PointAt(t1);
        //        }
        //        else
        //        {
        //            pointOnPrevHalfedge = firstH_Segment.To + moveVector;
        //        }
        //    }
        //    else
        //    {
        //        pointOnPrevHalfedge = firstH_Segment.To + moveVector;
        //    }


        //    #region 构造和更新relatedVertice和relatedPoints列表
        //    // 注意，此时的relatedVertice和relatedPoints并不是最终的结果
        //    // 注意，segment方向与半边方向是反的，所以先添加To，然后添加From
        //    List<int> relatedVertice = new List<int>();
        //    relatedVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));

        //    List<Point3d> relatedPoints = new List<Point3d>();
        //    relatedPoints.Add(firstH_Segment.To);
        //    relatedPoints.Add(firstH_Segment.From);

        //    relatedPoints.Insert(0, pointOnPrevHalfedge);
        //    relatedPoints.Add(pointOnNextHalfedge);
        //    #endregion

        //    #region 生成正确的realistcSiteFaceVertice和realisticSiteFacePoints
        //    realisticSiteFaceVertice = new List<int>();
        //    realisticSiteFaceVertice.AddRange(relatedVertice);
        //    realisticSiteFacePoints = new List<Point3d>();
        //    realisticSiteFacePoints.AddRange(relatedPoints);

        //    for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //    {
        //        vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //    }
        //    #endregion

        //    #region 生成正确的FE
        //    List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //    for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //    {
        //        Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //        string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //        List<int> includeDV = new List<int>();
        //        includeDV.Add(realisticSiteFaceVertice[i]);
        //        includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //        FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //        realisticSiteFaceEdgeSegments.Add(FE);
        //    }
        //    #endregion

        //    // 分割产生的子段
        //    if (prevHasDrawnFlag)
        //    {
        //        if (tBelow0)
        //        {
        //            // 不用构造没有包含在realisticSiteFaceEdgeSegments的子段
        //        }
        //        else if (tAbove1)
        //        {
        //            // 不用构造没有包含在realisticSiteFaceEdgeSegments的子段
        //        }
        //        else
        //        {
        //            #region 构造没有包含在realisticSiteFaceEdgeSegments的子段
        //            // 如果prevH已经绘制过，就既要新增prevH的子段，也要新增nextH的子段
        //            // v1,v2的构造不同于上一个分支
        //            int v1 = prevIndexPair[0];
        //            // int v1 = realisticSiteFacePoints.IndexOf(firstH_Segment.To);
        //            int v2 = realisticSiteFaceVertice[realisticSiteFacePoints.IndexOf(pointOnPrevHalfedge)];

        //            Line newSubPrevLine = new Line(pointOnPrevHalfedge, drawnFE[drawnPrevH].To);
        //            string newSubPrevLabel = string.Format("v:{0},{1}", v1, v2);
        //            List<int> newSubPrevIncludeDV = new List<int>();
        //            newSubPrevIncludeDV.Add(v1);
        //            newSubPrevIncludeDV.Add(v2);
        //            FaceEdgeSegment newSubPrevFE = new FaceEdgeSegment(newSubPrevLine, newSubPrevLabel, newSubPrevIncludeDV, -1);

        //            int v3 = nextIndexPair[0];
        //            //int v3 = realisticSiteFacePoints.IndexOf(nextFE.To);
        //            int v4 = realisticSiteFaceVertice[realisticSiteFacePoints.IndexOf(pointOnNextHalfedge)];

        //            Line newSubNextLine = new Line(pointOnNextHalfedge, nextFE.To);
        //            string newSubNextLabel = string.Format("v:{0},{1}", v3, v4);
        //            List<int> newSubNextIncludeDV = new List<int>();
        //            newSubNextIncludeDV.Add(v3);
        //            newSubNextIncludeDV.Add(v4);
        //            FaceEdgeSegment newSubNextFE = new FaceEdgeSegment(newSubNextLine, newSubNextLabel, newSubNextIncludeDV, -1);
        //            #endregion

        //            #region 其他需要更新的参数
        //            /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //             * 目的是为了让drawn列表中没有重复的对象
        //             */
        //            AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);

        //            /* 添加新的子段，移除旧的父段 */
        //            drawnDHV.Add(newSubPrevIncludeDV.ToArray());
        //            drawnDHV.Add(newSubNextIncludeDV.ToArray());
        //            drawnFE.Add(newSubPrevFE);
        //            drawnFE.Add(newSubNextFE);

        //            RemoveOldFE(prevIndexPair, drawnDHV, drawnFE);
        //            RemoveOldFE(nextIndexPair, drawnDHV, drawnFE);
        //            #endregion
        //        }
        //    }
        //    else
        //    {
        //        #region 构造没有包含在realisticSiteFaceEdgeSegments的子段
        //        // 如果prevH没有绘制过，就只要新增nextH的子段即可
        //        int v1 = nextFE.IncludedDVertice[0];
        //        int v2 = realisticSiteFaceVertice[realisticSiteFacePoints.IndexOf(pointOnNextHalfedge)];

        //        Line newLine = new Line(pointOnNextHalfedge, nextFE.To);
        //        string newLabel = string.Format("v:{0},{1}", v1, v2);
        //        List<int> newIncludeDV = new List<int>();
        //        newIncludeDV.Add(v1);
        //        newIncludeDV.Add(v2);
        //        FaceEdgeSegment newSubNextFE = new FaceEdgeSegment(newLine, newLabel, newIncludeDV, -1);
        //        #endregion

        //        #region 其他需要更新的参数
        //        /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //         * 目的是为了让drawn列表中没有重复的对象
        //         */
        //        AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);

        //        /* 添加新的子段，移除旧的父段 */
        //        drawnDHV.Add(newIncludeDV.ToArray());
        //        drawnFE.Add(newSubNextFE);

        //        RemoveOldFE(nextIndexPair, drawnDHV, drawnFE);
        //        #endregion
        //    }
        //    return newD;
        //}

        //private PlanktonMesh AFace_Operation2(PlanktonMesh D,
        //                                      List<List<FaceEdgeSegment>> subBoundarySegments,
        //                                      int currentFaceIndex,
        //                                      int[] secondIndexPair,
        //                                      out List<int> realisticSiteFaceVertice,
        //                                      out List<Point3d> realisticSiteFacePoints,
        //                                      List<FaceEdgeSegment> drawnFE,
        //                                      List<int[]> drawnDHV,
        //                                      List<int[]> abandonedDHV,
        //                                      Dictionary<int, Point3d> vi_Point3d,
        //                                      Dictionary<int, List<int>> fIndex_SplitedIndex)
        //{
        //    #region 找到secondIndexPair对应的currH
        //    int currH = -1;
        //    int[] hs = D.Faces.GetHalfedges(currentFaceIndex);
        //    for (int i = 0; i < hs.Length; i++)
        //    {
        //        int[] hV = new int[2] { D.Halfedges[hs[i]].StartVertex, D.Halfedges.EndVertex(hs[i]) };
        //        if (hV.Except(secondIndexPair).ToList().Count == 0)
        //        {
        //            currH = hs[i];
        //            break;
        //        }
        //    }
        //    #endregion

        //    #region 找到currH的prevprevH和nextH
        //    //int prev = D.Halfedges[currH].PrevHalfedge;
        //    //int prevprevH = D.Halfedges[D.Halfedges[currH].PrevHalfedge].PrevHalfedge;
        //    int nextH = D.Halfedges[currH].NextHalfedge;

        //    // 跳过abandonedDHV，找到真正的prevprevH
        //    bool hasAbandonedFlag = false;
        //    int h = currH;
        //    int iter = 0;
        //    do
        //    {
        //        int[] indexPair = new int[2] { D.Halfedges[h].StartVertex, D.Halfedges.EndVertex(h) };
        //        bool flag = false;
        //        for (int i = 0; i < abandonedDHV.Count; i++)
        //        {
        //            if (abandonedDHV[i].Except(indexPair).ToList().Count == 0)
        //            {
        //                flag = true;
        //                hasAbandonedFlag = true;
        //                break;
        //            }
        //        }

        //        if (!flag)
        //        {
        //            h = D.Halfedges[h].PrevHalfedge;
        //            iter++;
        //        }
        //        else
        //        {
        //            // abandonedDHV[i]和它的前一条边，本质上应该是都不存在的
        //            // 故这里选择了前一条边的前一条边
        //            h = D.Halfedges[D.Halfedges[h].PrevHalfedge].PrevHalfedge;
        //        }
        //    } while (iter < 2);

        //    int prevprevH = h;
        //    //int[] prevIndexPair = new int[2] { D.Halfedges[prev].StartVertex, D.Halfedges.EndVertex(prev) };
        //    int[] prevPrevIndexPair = new int[2] { D.Halfedges[prevprevH].StartVertex, D.Halfedges.EndVertex(prevprevH) };
        //    int[] nextIndexPair = new int[2] { D.Halfedges[nextH].StartVertex, D.Halfedges.EndVertex(nextH) };
        //    int[] currIndexPair = new int[2] { D.Halfedges[currH].StartVertex, D.Halfedges.EndVertex(currH) };
        //    #endregion
        //    #region 找到currH,prevprevH和nextH所对应的FE
        //    FaceEdgeSegment prevPrevFE = null;
        //    FaceEdgeSegment nextFE = null;

        //    //FaceEdgeSegment prevFE = null;
        //    for (int i = 0; i < drawnDHV.Count; i++)
        //    {
        //        if (drawnDHV[i].Except(prevPrevIndexPair).ToList().Count == 0)
        //        {
        //            prevPrevFE = drawnFE[i];
        //            continue;
        //        }
        //        if (drawnDHV[i].Except(nextIndexPair).ToList().Count == 0)
        //        {
        //            nextFE = drawnFE[i];
        //            continue;
        //        }
        //        //if (drawnDHV[i].Except(prevIndexPair).ToList().Count == 0)
        //        //{
        //        //    prevFE = drawnFE[i];
        //        //    continue;
        //        //}
        //    }
        //    FaceEdgeSegment currFE = null;
        //    for (int i = 0; i < subBoundarySegments.Count; i++)
        //    {
        //        for (int j = 0; j < subBoundarySegments[i].Count; j++)
        //        {
        //            if (subBoundarySegments[i][j].IncludedDVertice.Except(currIndexPair).ToList().Count == 0)
        //            {
        //                currFE = subBoundarySegments[i][j];
        //                break;
        //            }
        //        }
        //    }
        //    #endregion

        //    PlanktonMesh newD = new PlanktonMesh();
        //    realisticSiteFaceVertice = new List<int>();
        //    realisticSiteFacePoints = new List<Point3d>();

        //    // 对于包含abandonedDHV的特殊情况
        //    if (hasAbandonedFlag)
        //    {
        //        int prevH = D.Halfedges[prevprevH].NextHalfedge;
        //        int[] prevIndexPair = new int[2] { D.Halfedges[prevH].StartVertex, D.Halfedges.EndVertex(prevH) };
        //        int nextAnotherH = D.Halfedges[D.Halfedges.GetPairHalfedge(nextH)].PrevHalfedge;
        //        int[] nextAnotherIndexPair = new int[2] { D.Halfedges[nextAnotherH].StartVertex, D.Halfedges.EndVertex(nextAnotherH) };

        //        #region 找到prevH和nextAnotherH所对应的FE
        //        FaceEdgeSegment prevFE = null;
        //        FaceEdgeSegment nextAnotherFE = null;
        //        for (int i = 0; i < drawnDHV.Count; i++)
        //        {
        //            if (drawnDHV[i].Except(prevIndexPair).ToList().Count == 0)
        //            {
        //                prevFE = drawnFE[i];
        //                continue;
        //            }
        //            if (drawnDHV[i].Except(nextAnotherIndexPair).ToList().Count == 0)
        //            {
        //                nextAnotherFE = drawnFE[i];
        //                continue;
        //            }
        //        }
        //        #endregion

        //        bool isOutOfPrevFERange = false;
        //        bool isOutOfNextAnotherFERange = false;
        //        double tolerance1 = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
        //        double t3;
        //        double t4;
        //        if (Intersection.LineLine(prevFE.Line, nextAnotherFE.Line, out t3, out t4, tolerance1, false))
        //        {
        //            double t3Round = Math.Round(t3, 2);
        //            double t4Round = Math.Round(t4, 2);

        //            if (t3Round > 1 || t3Round < 0)
        //            {
        //                isOutOfPrevFERange = true;
        //            }
        //            if (t4Round > 1 || t4Round < 0)
        //            {
        //                isOutOfNextAnotherFERange = true;
        //            }

        //            #region 修改D
        //            int midVertex;

        //            #region 根据sortedAdjacentFace的顺序，得到对应的viAroundAdjacentFace列表
        //            List<int> viAroundFace = new List<int>();
        //            // 包围currFace的顶点
        //            viAroundFace = D.Faces.GetFaceVertices(currentFaceIndex).ToList();
        //            #endregion

        //            List<PlanktonXYZ> newPlanktonXYZ = new List<PlanktonXYZ>();
        //            List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //            newPlanktonXYZ = ModifyVertice_AFace_Operation2_1(D, viAroundFace, currentFaceIndex, D.Halfedges.EndVertex(prevH), D.Halfedges.EndVertex(nextAnotherH), out newFaceVertexOrder, out midVertex);

        //            newD = new PlanktonMesh();
        //            newD.Vertices.AddVertices(newPlanktonXYZ);
        //            newD.Faces.AddFaces(newFaceVertexOrder);
        //            #endregion

        //            Point3d intersectPoint = prevFE.Line.PointAt(t3);

        //            if (!isOutOfPrevFERange && !isOutOfNextAnotherFERange)
        //            {
        //                // 理论上不存在这种情况，因为在进行AFace_Operation1的操作时，就已经避免了这种情况的出现

        //                realisticSiteFaceVertice = null;
        //                realisticSiteFacePoints = null;
        //            }
        //            else if (!isOutOfPrevFERange && isOutOfNextAnotherFERange)
        //            {
        //                #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //                realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                    {
        //                        realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                    }
        //                    else
        //                    {
        //                        realisticSiteFacePoints.Add(intersectPoint);
        //                    }
        //                }

        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //                }
        //                #endregion

        //                #region 生成正确的FE
        //                List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                    string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    List<int> includeDV = new List<int>();
        //                    includeDV.Add(realisticSiteFaceVertice[i]);
        //                    includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                    realisticSiteFaceEdgeSegments.Add(FE);
        //                }
        //                #endregion

        //                #region 构造没有包含在realisticSiteFaceEdgeSegments的子段

        //                #endregion
        //            }
        //            else if (isOutOfPrevFERange && !isOutOfNextAnotherFERange)
        //            {
        //                #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //                realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                    {
        //                        realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                    }
        //                    else
        //                    {
        //                        realisticSiteFacePoints.Add(intersectPoint);
        //                    }
        //                }

        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //                }
        //                #endregion

        //                #region 生成正确的FE
        //                List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                    string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    List<int> includeDV = new List<int>();
        //                    includeDV.Add(realisticSiteFaceVertice[i]);
        //                    includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                    realisticSiteFaceEdgeSegments.Add(FE);
        //                }
        //                #endregion

        //                #region 构造没有包含在realisticSiteFaceEdgeSegments的子段

        //                #endregion
        //            }
        //            else
        //            {
        //                #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //                realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                    {
        //                        realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                    }
        //                    else
        //                    {
        //                        realisticSiteFacePoints.Add(intersectPoint);
        //                    }
        //                }

        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //                }
        //                #endregion

        //                #region 生成正确的FE
        //                List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //                for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //                {
        //                    Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                    string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    List<int> includeDV = new List<int>();
        //                    includeDV.Add(realisticSiteFaceVertice[i]);
        //                    includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                    FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                    realisticSiteFaceEdgeSegments.Add(FE);
        //                }
        //                #endregion

        //                // 因为都超出了，所以不用构造没有包含在realisticSiteFaceEdgeSegments的子段

        //                #region 其他要更新的参数
        //                /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //                 * 目的是为了让drawn列表中没有重复的对象
        //                 */
        //                AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);

        //                // 也不需要添加新的子段，移除旧的父段
        //                #endregion
        //            }

        //            return newD;
        //        }
        //    }


        //    // 对于不包含abandonedDHV的正常情况
        //    // 求交点
        //    bool isOutOfPrevPrevFERange = false;
        //    bool isOutOfNextFERange = false;
        //    double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
        //    double t1;
        //    double t2;
        //    if (Intersection.LineLine(prevPrevFE.Line, nextFE.Line, out t1, out t2, tolerance, false))
        //    {
        //        double t1Round = Math.Round(t1, 2);
        //        double t2Round = Math.Round(t2, 2);

        //        if (t1Round > 1 || t1Round < 0)
        //        {
        //            isOutOfPrevPrevFERange = true;
        //        }
        //        if (t2Round > 1 || t2Round < 0)
        //        {
        //            isOutOfNextFERange = true;
        //        }

        //        #region 修改D
        //        List<int> sortedPairAdjacentFace = new List<int>();
        //        List<int> sortedPublicHalfedge = new List<int>();
        //        // 确定正确的PublicHalfedge顺序和AdjacentFace顺序
        //        sortedPairAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(prevprevH)].AdjacentFace);
        //        sortedPairAdjacentFace.Add(D.Halfedges[D.Halfedges.GetPairHalfedge(nextH)].AdjacentFace);
        //        sortedPublicHalfedge.Add(prevprevH);
        //        sortedPublicHalfedge.Add(nextH);

        //        int midVertex;
        //        newD = RegeneratePMesh_AFace_Operation2(D, currentFaceIndex, sortedPairAdjacentFace, sortedPublicHalfedge, isOutOfPrevPrevFERange, isOutOfNextFERange, fIndex_SplitedIndex, out midVertex);
        //        #endregion

        //        List<string> debug3 = UtilityFunctions.PrintFacesHalfedges(newD);
        //        List<string> debug4 = UtilityFunctions.PrintHalfedgeStartAndEnd(newD);
        //        List<string> debug5 = UtilityFunctions.PrintFacesVertices(newD);

        //        Point3d intersectPoint = prevPrevFE.Line.PointAt(t1);

        //        if (!isOutOfPrevPrevFERange && !isOutOfNextFERange)
        //        {
        //            // 理论上不存在这种情况，因为在进行AFace_Operation1的操作时，就已经避免了这种情况的出现

        //            realisticSiteFaceVertice = null;
        //            realisticSiteFacePoints = null;

        //        }
        //        else if (!isOutOfPrevPrevFERange && isOutOfNextFERange)
        //        {
        //            #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //            realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                {
        //                    realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                }
        //                else
        //                {
        //                    realisticSiteFacePoints.Add(intersectPoint);
        //                }
        //            }

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //            }
        //            #endregion

        //            #region 生成正确的FE
        //            List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                List<int> includeDV = new List<int>();
        //                includeDV.Add(realisticSiteFaceVertice[i]);
        //                includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                realisticSiteFaceEdgeSegments.Add(FE);
        //            }
        //            #endregion

        //            #region 构造没有包含在realisticSiteFaceEdgeSegments的子段
        //            int v1 = prevPrevIndexPair[0];
        //            int v2 = realisticSiteFaceVertice[realisticSiteFacePoints.IndexOf(intersectPoint)];
        //            Line newSubPrevPrevLine = new Line(intersectPoint, vi_Point3d[prevPrevIndexPair[0]]);
        //            string newSubPrevPrevLabel = string.Format("v:{0},{1}", v1, v2);
        //            List<int> newSubPrevPrevIncludeDV = new List<int>();
        //            newSubPrevPrevIncludeDV.Add(v1);
        //            newSubPrevPrevIncludeDV.Add(v2);
        //            FaceEdgeSegment newSubPrevPrevFE = new FaceEdgeSegment(newSubPrevPrevLine, newSubPrevPrevLabel, newSubPrevPrevIncludeDV, -1);
        //            #endregion

        //            #region 其他要更新的参数
        //            /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //             * 目的是为了让drawn列表中没有重复的对象
        //             */
        //            AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);

        //            /* 添加新的子段，移除旧的父段 */
        //            drawnDHV.Add(newSubPrevPrevFE.IncludedDVertice.ToArray());
        //            drawnFE.Add(newSubPrevPrevFE);

        //            RemoveOldFE(prevPrevIndexPair, drawnDHV, drawnFE);
        //            #endregion

        //        }
        //        else if (isOutOfPrevPrevFERange && !isOutOfNextFERange)
        //        {
        //            #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //            realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                {
        //                    realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                }
        //                else
        //                {
        //                    realisticSiteFacePoints.Add(intersectPoint);
        //                }
        //            }

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //            }
        //            #endregion

        //            #region 生成正确的FE
        //            List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                List<int> includeDV = new List<int>();
        //                includeDV.Add(realisticSiteFaceVertice[i]);
        //                includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                realisticSiteFaceEdgeSegments.Add(FE);
        //            }
        //            #endregion

        //            #region 构造没有包含在realisticSiteFaceEdgeSegments的子段
        //            int v1 = realisticSiteFaceVertice[realisticSiteFacePoints.IndexOf(intersectPoint)];
        //            int v2 = nextIndexPair[1];
        //            Line newSubNextLine = new Line(vi_Point3d[nextIndexPair[1]], intersectPoint);
        //            string newSubNextLabel = string.Format("v:{0},{1}", v1, v2);
        //            List<int> newSubNextIncludeDV = new List<int>();
        //            newSubNextIncludeDV.Add(v1);
        //            newSubNextIncludeDV.Add(v2);
        //            FaceEdgeSegment newSubNextFE = new FaceEdgeSegment(newSubNextLine, newSubNextLabel, newSubNextIncludeDV, -1);
        //            #endregion

        //            #region 其他要更新的参数
        //            /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //             * 目的是为了让drawn列表中没有重复的对象
        //             */
        //            AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);

        //            /* 添加新的子段，移除旧的父段 */
        //            drawnDHV.Add(newSubNextFE.IncludedDVertice.ToArray());
        //            drawnFE.Add(newSubNextFE);

        //            RemoveOldFE(nextIndexPair, drawnDHV, drawnFE);
        //            #endregion
        //        }
        //        else
        //        {
        //            #region 生成正确的realisticSiteFaceVertice和realisticSiteFacePoints
        //            realisticSiteFaceVertice.AddRange(newD.Faces.GetFaceVertices(newD.Faces.Count - 1));

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                if (vi_Point3d.ContainsKey(realisticSiteFaceVertice[i]))
        //                {
        //                    realisticSiteFacePoints.Add(vi_Point3d[realisticSiteFaceVertice[i]]);
        //                }
        //                else
        //                {
        //                    realisticSiteFacePoints.Add(intersectPoint);
        //                }
        //            }

        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                vi_Point3d.TryAdd(realisticSiteFaceVertice[i], realisticSiteFacePoints[i]);
        //            }
        //            #endregion

        //            // 因为都超出了，所以不用构造没有包含在realisticSiteFaceEdgeSegments的子段

        //            #region 生成正确的FE
        //            List<FaceEdgeSegment> realisticSiteFaceEdgeSegments = new List<FaceEdgeSegment>();
        //            for (int i = 0; i < realisticSiteFaceVertice.Count; i++)
        //            {
        //                Line line = new Line(realisticSiteFacePoints[(i + 1) % realisticSiteFaceVertice.Count], realisticSiteFacePoints[i]);
        //                string label = string.Format("v:{0},{1}", realisticSiteFaceVertice[i], realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                List<int> includeDV = new List<int>();
        //                includeDV.Add(realisticSiteFaceVertice[i]);
        //                includeDV.Add(realisticSiteFaceVertice[(i + 1) % realisticSiteFaceVertice.Count]);
        //                FaceEdgeSegment FE = new FaceEdgeSegment(line, label, includeDV, -1);

        //                realisticSiteFaceEdgeSegments.Add(FE);
        //            }
        //            #endregion

        //            #region 其他要更新的参数
        //            /* 输出最终正确的drawnDHV，以及输出最终正确的feHasH
        //             * 目的是为了让drawn列表中没有重复的对象
        //             */
        //            AddrealisticSiteFaceEdgeSegments(realisticSiteFaceEdgeSegments, drawnDHV, drawnFE);
                    
        //            // 也不需要添加新的子段，移除旧的父段
        //            #endregion
        //        }

        //    }
        //    else
        //    {
        //        newD = new PlanktonMesh(D);
        //    }

        //    return newD;
        //}

       

        //private PlanktonMesh RegeneratePMesh_3_Operation(PlanktonMesh D,
        //                                                 int currentBFaceIndexThreeSideDetermined,
        //                                                 int adjacentFaceIndex,
        //                                                 int index,
        //                                                 List<int> currentFaceHIndexList)
        //{
        //    #region 生成新的点
        //    List<int> viAroundBFace = new List<int>();
        //    // 包围AdjacentFace的顶点
        //    viAroundBFace = D.Faces.GetFaceVertices(currentBFaceIndexThreeSideDetermined).ToList();
        //    List<int> viAroundPairBFace = new List<int>();
        //    viAroundPairBFace = D.Faces.GetFaceVertices(adjacentFaceIndex).ToList();

        //    // 即公共边的起点与终点
        //    int startVertexIndex = D.Halfedges[currentFaceHIndexList[index + 1]].StartVertex;
        //    int endVertexIndex = D.Halfedges.EndVertex(currentFaceHIndexList[index + 1]);

        //    // 构造图中新的顶点
        //    Point3d startVertex = D.Vertices[startVertexIndex].ToPoint3d();
        //    Point3d endVertex = D.Vertices[endVertexIndex].ToPoint3d();
        //    Point3d midVertex = (startVertex + endVertex) / 2;
        //    D.Vertices.Add(midVertex);
        //    int midVertexIndex = D.Vertices.Count - 1;
        //    int newestVertexIndex = midVertexIndex;

        //    #region 构造添加顶点后的新Face
        //    for (int i = 0; i < viAroundBFace.Count; i++)
        //    {
        //        if (viAroundBFace[i] == startVertexIndex)
        //        {
        //            viAroundBFace.Insert(i + 1, midVertexIndex);
        //            break;
        //        }
        //    }
        //    for (int i = 0; i < viAroundPairBFace.Count; i++)
        //    {
        //        if (viAroundPairBFace[i] == endVertexIndex)
        //        {
        //            viAroundPairBFace.Insert(i + 1, midVertexIndex);
        //            break;
        //        }
        //    }
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    #region 生成新的Face
        //    // 新增新增顶点与endVertex前一个顶点的联系
        //    for (int i = 0; i < viAroundPairBFace.Count; i++)
        //    {
        //        if (viAroundPairBFace[i] == startVertexIndex)
        //        {
        //            // 向viAroundPairAdjacentFace中添加startVertex的后一个点
        //            viAroundPairBFace.Insert(i, viAroundBFace[((viAroundBFace.IndexOf(startVertexIndex) - 1) + viAroundBFace.Count) % viAroundBFace.Count]);
        //            break;
        //        }
        //    }
        //    // 移除新增顶点与endVertex之间的联系
        //    viAroundBFace.Remove(startVertexIndex);

        //    #endregion
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonVertex = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonVertex.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来D的Face属性，同时替换掉两个应该删除的面
        //    List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentBFaceIndexThreeSideDetermined)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundBFace);
        //        }
        //        else if (i == adjacentFaceIndex)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundPairBFace);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }

        //    D = new PlanktonMesh();
        //    D.Vertices.AddVertices(newPlanktonVertex);
        //    D.Faces.AddFaces(newFaceVertexOrder);

        //    // debug
        //    List<string> debugPrint1 = UtilityFunctions.PrintFacesHalfedges(D);
        //    List<string> debugPrint2 = UtilityFunctions.PrintHalfedgeStartAndEnd(D);
        //    #endregion

        //    return D;
        //}

        //private PlanktonMesh RegeneratePMesh_2_Operation(PlanktonMesh D,
        //                                                 int currentBFaceIndexTwoSideDetermined,
        //                                                 List<int> sortedAdjacentFace,
        //                                                 List<int> sortedPublicHalfedge,
        //                                                 bool isOutOfRange,
        //                                                 bool flag0,
        //                                                 bool flag1)
        //{
        //    List<int> viAroundBFace = new List<int>();
        //    // 包围AdjacentFace的顶点
        //    viAroundBFace = D.Faces.GetFaceVertices(currentBFaceIndexTwoSideDetermined).ToList();
        //    List<List<int>> viAroundPairBFace = new List<List<int>>();
        //    for (int i = 0; i < sortedAdjacentFace.Count; i++)
        //    {
        //        viAroundPairBFace.Add(new List<int>());
        //        viAroundPairBFace[i].AddRange(D.Faces.GetFaceVertices(sortedAdjacentFace[i]).ToList());
        //    }

        //    if (isOutOfRange)
        //    {
        //        // 即公共边的起点与终点
        //        int startVertexIndex0 = D.Halfedges.EndVertex(sortedPublicHalfedge[0]);
        //        // int endVertexIndex0 = D.Halfedges.EndVertex(sortedPublicHalfedge[0]);
        //        int anotherEndVertexIndex0 = -1;
        //        for (int i = 0; i < D.Faces.GetHalfedges(sortedAdjacentFace[1]).Length; i++)
        //        {
        //            if (D.Halfedges[D.Faces.GetHalfedges(sortedAdjacentFace[1])[i]].StartVertex == startVertexIndex0)
        //            {
        //                anotherEndVertexIndex0 = D.Halfedges.EndVertex(D.Faces.GetHalfedges(sortedAdjacentFace[1])[i]);
        //            }
        //        }

        //        int startVertexIndex1 = D.Halfedges[sortedPublicHalfedge[1]].StartVertex;
        //        // int endVertexIndex1 = D.Halfedges.EndVertex(sortedPublicHalfedge[1]);
        //        int anotherEndVertexIndex1 = -1;
        //        for (int i = 0; i < D.Faces.GetHalfedges(sortedAdjacentFace[0]).Length; i++)
        //        {
        //            if (D.Halfedges.EndVertex(D.Faces.GetHalfedges(sortedAdjacentFace[0])[i]) == startVertexIndex1)
        //            {
        //                anotherEndVertexIndex1 = D.Halfedges[D.Faces.GetHalfedges(sortedAdjacentFace[0])[i]].StartVertex;
        //            }
        //        }

        //        if (flag0 && !flag1)
        //        {
        //            // 全部选startVertexIndex0，endVertexIndex0
        //            List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //            List<PlanktonXYZ> newPlanktonVertex = ModifyFaceVertice_2_Operation_OutOfRange(D,
        //                                                                                           viAroundBFace,
        //                                                                                           viAroundPairBFace,
        //                                                                                           currentBFaceIndexTwoSideDetermined,
        //                                                                                           sortedAdjacentFace,
        //                                                                                           startVertexIndex0,
        //                                                                                           anotherEndVertexIndex0,
        //                                                                                           out newFaceVertexOrder);
        //            D = new PlanktonMesh();
        //            D.Vertices.AddVertices(newPlanktonVertex);
        //            D.Faces.AddFaces(newFaceVertexOrder);
        //        }
        //        else
        //        {
        //            // 全部选startVertexIndex1，endVertexIndex1
        //            List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //            List<PlanktonXYZ> newPlanktonVertex = ModifyFaceVertice_2_Operation_OutOfRange(D,
        //                                                                                           viAroundBFace,
        //                                                                                           viAroundPairBFace,
        //                                                                                           currentBFaceIndexTwoSideDetermined,
        //                                                                                           sortedAdjacentFace,
        //                                                                                           startVertexIndex1,
        //                                                                                           anotherEndVertexIndex1,
        //                                                                                           out newFaceVertexOrder);
        //            D = new PlanktonMesh();
        //            D.Vertices.AddVertices(newPlanktonVertex);
        //            D.Faces.AddFaces(newFaceVertexOrder);
        //        }

        //    }
        //    else
        //    {
        //        // 即公共边的起点与终点
        //        int startVertexIndex0 = D.Halfedges[sortedPublicHalfedge[0]].StartVertex;
        //        int endVertexIndex0 = D.Halfedges.EndVertex(sortedPublicHalfedge[0]);

        //        int startVertexIndex1 = D.Halfedges[sortedPublicHalfedge[1]].StartVertex;
        //        int endVertexIndex1 = D.Halfedges.EndVertex(sortedPublicHalfedge[1]);

        //        if (flag0 && !flag1)
        //        {
        //            // 全部选startVertexIndex0，endVertexIndex0
        //            List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //            List<PlanktonXYZ> newPlanktonVertex = ModifyFaceVertice_2_Operation_NotOutOfRange(D,
        //                                                                                              viAroundBFace,
        //                                                                                              viAroundPairBFace,
        //                                                                                              currentBFaceIndexTwoSideDetermined,
        //                                                                                              sortedAdjacentFace,
        //                                                                                              startVertexIndex0,
        //                                                                                              endVertexIndex0,
        //                                                                                              out newFaceVertexOrder);

        //            D = new PlanktonMesh();
        //            D.Vertices.AddVertices(newPlanktonVertex);
        //            D.Faces.AddFaces(newFaceVertexOrder);
        //        }
        //        else
        //        {
        //            // 全部选startVertexIndex1，endVertexIndex1
        //            List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //            List<PlanktonXYZ> newPlanktonVertex = ModifyFaceVertice_2_Operation_NotOutOfRange(D,
        //                                                                                              viAroundBFace,
        //                                                                                              viAroundPairBFace,
        //                                                                                              currentBFaceIndexTwoSideDetermined,
        //                                                                                              sortedAdjacentFace,
        //                                                                                              startVertexIndex1,
        //                                                                                              endVertexIndex1,
        //                                                                                              out newFaceVertexOrder);

        //            D = new PlanktonMesh();
        //            D.Vertices.AddVertices(newPlanktonVertex);
        //            D.Faces.AddFaces(newFaceVertexOrder);
        //        }
        //    }

        //    return D;
        //}

        //private PlanktonMesh RegeneratePMesh_AFace_Operation1(PlanktonMesh D,
        //                                                      int currentFaceIndex,
        //                                                      List<int> sortedAdjacentFace,
        //                                                      //List<int> sortedAdjacentHalfedge,
        //                                                      int preStartVertex,
        //                                                      int startVertex,
        //                                                      int endVertex,
        //                                                      int nextEndVertex,
        //                                                      Dictionary<int, List<int>> fIndex_SplitedIndex,
        //                                                      out int prevMidVertex,
        //                                                      out int nextMidVertex)
        //{
        //    List<int> viAroundFace = new List<int>();
        //    // 包围currFace的顶点
        //    viAroundFace = D.Faces.GetFaceVertices(currentFaceIndex).ToList();
        //    List<List<int>> viAroundAdjacentFace = new List<List<int>>();
        //    for (int i = 0; i < sortedAdjacentFace.Count; i++)
        //    {
        //        viAroundAdjacentFace.Add(new List<int>());
        //        viAroundAdjacentFace[i].AddRange(D.Faces.GetFaceVertices(sortedAdjacentFace[i]).ToList());
        //    }

        //    List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //    List<PlanktonXYZ> newPlanktonXYZ = ModifyFaceVertice_AFace_Operation1(D,viAroundFace,
        //                                                                         viAroundAdjacentFace[0],
        //                                                                         viAroundAdjacentFace[1],
        //                                                                         currentFaceIndex,
        //                                                                         sortedAdjacentFace[0],
        //                                                                         sortedAdjacentFace[1],preStartVertex,
        //                                                                         startVertex,
        //                                                                         endVertex,
        //                                                                         nextEndVertex,
        //                                                                         out newFaceVertexOrder,
        //                                                                         out prevMidVertex,
        //                                                                         out nextMidVertex);
        //    D = new PlanktonMesh();
        //    D.Vertices.AddVertices(newPlanktonXYZ);
        //    D.Faces.AddFaces(newFaceVertexOrder);

        //    List<int> value = new List<int>();
        //    value.Add(currentFaceIndex);
        //    value.Add(D.Faces.Count - 1);

        //    fIndex_SplitedIndex.Add(currentFaceIndex, value);

        //    return D;
        //}

        //private PlanktonMesh RegeneratePMesh_AFace_Operation2(PlanktonMesh D,
        //                                                      int currentFaceIndex,
        //                                                      List<int> sortedAdjacentFace,
        //                                                      List<int> sortedPublicHalfedge,
        //                                                      bool isOutOfPrevPrevFERange,
        //                                                      bool isOutOfNextFERange,
        //                                                      Dictionary<int, List<int>> fIndex_SplitedIndex,
        //                                                      out int midVertex)
        //{
        //    #region 根据sortedAdjacentFace的顺序，得到对应的viAroundAdjacentFace列表
        //    List<int> viAroundFace = new List<int>();
        //    // 包围currFace的顶点
        //    viAroundFace = D.Faces.GetFaceVertices(currentFaceIndex).ToList();
        //    List<List<int>> viAroundAdjacentFace = new List<List<int>>();
        //    for (int i = 0; i < sortedAdjacentFace.Count; i++)
        //    {
        //        viAroundAdjacentFace.Add(new List<int>());
        //        viAroundAdjacentFace[i].AddRange(D.Faces.GetFaceVertices(sortedAdjacentFace[i]).ToList());
        //    }
        //    #endregion

        //    // startVertexOnPrevPrevH
        //    int startVertex0 = D.Halfedges[sortedPublicHalfedge[0]].StartVertex;
        //    // endVertexOnPrevPrevH
        //    int endVertex0 = D.Halfedges.EndVertex(sortedPublicHalfedge[0]);
        //    // startVertexOnNextH
        //    int startVertex1 = D.Halfedges[sortedPublicHalfedge[1]].StartVertex;
        //    // endVertexOnNextH
        //    int endVertex1 = D.Halfedges.EndVertex(sortedPublicHalfedge[1]);

        //    List<List<int>> newFaceVertexOrder = new List<List<int>>();
        //    List<PlanktonXYZ> newPlanktonXYZ = ModifyFaceVertice_AFace_Operation2(D,
        //                                                                          viAroundFace,
        //                                                                          viAroundAdjacentFace[0],
        //                                                                          viAroundAdjacentFace[1],
        //                                                                          currentFaceIndex,
        //                                                                          sortedAdjacentFace[0],
        //                                                                          sortedAdjacentFace[1],
        //                                                                          startVertex0,
        //                                                                          endVertex0,
        //                                                                          startVertex1,
        //                                                                          endVertex1,
        //                                                                          isOutOfPrevPrevFERange,
        //                                                                          isOutOfNextFERange,
        //                                                                          out newFaceVertexOrder,
        //                                                                          out midVertex);

        //    D = new PlanktonMesh();
        //    D.Vertices.AddVertices(newPlanktonXYZ);
        //    D.Faces.AddFaces(newFaceVertexOrder);

        //    List<int> value = new List<int>();
        //    value.Add(currentFaceIndex);
        //    value.Add(D.Faces.Count - 1);

        //    /* 需要进行添加的修正 */
        //    fIndex_SplitedIndex.TryAdd(currentFaceIndex, value);

        //    return D;
        //}

        //private List<PlanktonXYZ> ModifyFaceVertice_2_Operation_OutOfRange(PlanktonMesh D,
        //                                                                   List<int> viAroundBFace,
        //                                                                   List<List<int>> viAroundPairBFace,
        //                                                                   int currentBFaceIndexTwoSideDetermined,
        //                                                                   List<int> sortedAdjacentFace,
        //                                                                   int startVertexIndex,
        //                                                                   int anotherEndVertexIndex,
        //                                                                   out List<List<int>> newFaceVertexOrder)
        //{
        //    #region 生成新的点
        //    Point3d startVertex = D.Vertices[startVertexIndex].ToPoint3d();
        //    Point3d anotherEndVertex = D.Vertices[anotherEndVertexIndex].ToPoint3d();
        //    Point3d midVertex = (startVertex + anotherEndVertex) / 2;
        //    D.Vertices.Add(midVertex);
        //    int midVertexIndex = D.Vertices.Count - 1;
        //    int newestVertexIndex = midVertexIndex;

        //    #region 构造添加顶点后的新的viAroundBFace和viAroundPairBFace[0]
        //    for (int i = 0; i < viAroundBFace.Count; i++)
        //    {
        //        if (viAroundBFace[i] == startVertexIndex)
        //        {
        //            viAroundBFace.Insert(i + 1, midVertexIndex);
        //            break;
        //        }
        //    }
        //    // viAroundPairBFace[0]不用修改
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    #region 生成新的viAroundPairBFace[1]
        //    for (int i = 0; i < viAroundPairBFace[1].Count; i++)
        //    {
        //        if (viAroundPairBFace[1][i] == startVertexIndex)
        //        {
        //            // 向viAroundPairBFace[1]中把startVertexIndex替换成midVertexIndex点
        //            viAroundPairBFace[1][i] = midVertexIndex;
        //            break;
        //        }
        //    }
        //    // 不需要在viAroundBFace中移除startVertex
        //    #endregion
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonVertex = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonVertex.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来D的Face属性，同时替换掉应该修改的面
        //    newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentBFaceIndexTwoSideDetermined)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundBFace);
        //        }
        //        else if (i == sortedAdjacentFace[1])
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundPairBFace[1]);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }
        //    #endregion
        //    return newPlanktonVertex;

            
        //}

        //private List<PlanktonXYZ> ModifyFaceVertice_2_Operation_NotOutOfRange(PlanktonMesh D, 
        //                                                                      List<int> viAroundBFace,
        //                                                                      List<List<int>> viAroundPairBFace,
        //                                                                      int currentBFaceIndexTwoSideDetermined,
        //                                                                      List<int> sortedAdjacentFace,
        //                                                                      int startVertexIndex,
        //                                                                      int endVertexIndex,
        //                                                                      out List<List<int>> newFaceVertexOrder)
        //{
        //    #region 生成新的点
        //    Point3d startVertex = D.Vertices[startVertexIndex].ToPoint3d();
        //    Point3d endVertex = D.Vertices[endVertexIndex].ToPoint3d();
        //    Point3d midVertex = (startVertex + endVertex) / 2;
        //    D.Vertices.Add(midVertex);
        //    int midVertexIndex = D.Vertices.Count - 1;
        //    int newestVertexIndex = midVertexIndex;

        //    #region 构造添加顶点后的新的viAroundBFace和viAroundPairBFace[0]
        //    for (int i = 0; i < viAroundBFace.Count; i++)
        //    {
        //        if (viAroundBFace[i] == startVertexIndex)
        //        {
        //            viAroundBFace.Insert(i + 1, midVertexIndex);
        //            break;
        //        }
        //    }
        //    for (int i = 0; i < viAroundPairBFace[0].Count; i++)
        //    {
        //        if (viAroundPairBFace[0][i] == startVertexIndex)
        //        {
        //            viAroundPairBFace[0].Insert(i, midVertexIndex);
        //            break;
        //        }
        //    }
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    #region 生成新的viAroundPairBFace[1]
        //    for (int i = 0; i < viAroundPairBFace[1].Count; i++)
        //    {
        //        if (viAroundPairBFace[1][i] == endVertexIndex)
        //        {
        //            // 向viAroundPairBFace[1]中添加midVertex点
        //            viAroundPairBFace[1].Insert(i, midVertexIndex);
        //            break;
        //        }
        //    }
        //    // 在viAroundBFace中移除startVertexIndex点
        //    viAroundBFace.Remove(endVertexIndex);
        //    #endregion
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonVertex = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonVertex.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来D的Face属性，同时替换掉三个应该修改的面
        //    newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentBFaceIndexTwoSideDetermined)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundBFace);
        //        }
        //        else if (i == sortedAdjacentFace[0])
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundPairBFace[0]);
        //        }
        //        else if (i == sortedAdjacentFace[1])
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundPairBFace[1]);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }
        //    #endregion
        //    return newPlanktonVertex;
        //}

        //private List<PlanktonXYZ> ModifyFaceVertice_AFace_Operation1(PlanktonMesh D,
        //                                                             List<int> viAroundFace,
        //                                                             List<int> viAroundAdjacentFace0,
        //                                                             List<int> viAroundAdjacentFace1,
        //                                                             int currentFaceIndex,
        //                                                             int adjacentFace0Index,
        //                                                             int adjacentFace1Index,
        //                                                             int startVertex0,
        //                                                             int endVertex0,
        //                                                             int startVertex1,
        //                                                             int endVertex1,
        //                                                             out List<List<int>> newFaceVertexOrder,
        //                                                             out int midVertex0,
        //                                                             out int midVertex1)
        //{
        //    #region 生成新的点
        //    Point3d startPoint0 = D.Vertices[startVertex0].ToPoint3d();
        //    Point3d endPoint0 = D.Vertices[endVertex0].ToPoint3d();
        //    Point3d startPoint1 = D.Vertices[startVertex1].ToPoint3d();
        //    Point3d endPoint1 = D.Vertices[endVertex1].ToPoint3d();

        //    Point3d midPoint0 = (startPoint0 + endPoint0) / 2;
        //    Point3d midPoint1 = (startPoint1 + endPoint1) / 2;

        //    D.Vertices.Add(midPoint0);
        //    midVertex0 = D.Vertices.Count - 1;
        //    D.Vertices.Add(midPoint1);
        //    midVertex1 = D.Vertices.Count - 1;
        //    #endregion

        //    #region 生成三个Face的Vertex列表
        //    List<int> allVerticeForCurrentFace = new List<int>();
        //    allVerticeForCurrentFace.AddRange(viAroundFace);

        //    // 这里必须要取一次index，插一次值
        //    int endVertex0InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(endVertex0);
        //    allVerticeForCurrentFace.Insert(endVertex0InAllVerticeIndex, midVertex0);
        //    int endVertex1InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(endVertex1);
        //    allVerticeForCurrentFace.Insert(endVertex1InAllVerticeIndex, midVertex1);

        //    int midVertex0InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(midVertex0);
        //    int midVertex1InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(midVertex1);

        //    #region 找到allVerticeForCurrentFace的两个切片
        //    List<int> viFace = new List<int>();
        //    int iter0 = midVertex0InAllVerticeIndex;
        //    while (allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != midVertex1)
        //    {
        //        iter0--;
        //        viFace.Add(allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //    }

        //    List<int> viNewFace = new List<int>();
        //    int iter1 = midVertex0InAllVerticeIndex;
        //    while (allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != midVertex1)
        //    {
        //        iter1++;
        //        viNewFace.Add(allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //    }

        //    viFace.Insert(0, midVertex0);
        //    viNewFace.Insert(0, midVertex0);

        //    viFace.Reverse();
        //    #endregion

        //    #region 处理viAroundPairBFace0和viAroundPairBFace1
        //    //int endVertex0OnPairBFace0 = viAroundAdjacentFace0.IndexOf(endVertex0);
        //    //viAroundAdjacentFace0.Insert(endVertex0OnPairBFace0, midVertex0);

        //    //int endVertex1OnPairBFace1 = viAroundAdjacentFace1.IndexOf(endVertex1);
        //    //viAroundAdjacentFace1.Insert(endVertex1OnPairBFace1, midVertex1);

        //    int startVertex0OnPairBFace0 = viAroundAdjacentFace0.IndexOf(startVertex0);
        //    viAroundAdjacentFace0.Insert(startVertex0OnPairBFace0, midVertex0);

        //    int startVertex1OnPairBFace1 = viAroundAdjacentFace1.IndexOf(startVertex1);
        //    viAroundAdjacentFace1.Insert(startVertex1OnPairBFace1, midVertex1);
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonXYZ = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonXYZ.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来D的Face属性，同时替换viAroundFace，并且添加viAroundPairFace
        //    // viAroundPairFace作为原来的面的替代，viAroundFace作为新生成的面
        //    newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentFaceIndex)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viFace);
        //        }
        //        else if (i == adjacentFace0Index)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundAdjacentFace0);
        //        }
        //        else if (i == adjacentFace1Index)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundAdjacentFace1);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }
        //    newFaceVertexOrder.Add(viNewFace);

        //    return newPlanktonXYZ;
        //    #endregion
        //}

        //private List<PlanktonXYZ> ModifyFaceVertice_AFace_Operation2(PlanktonMesh D, 
        //                                                             List<int> viAroundFace,
        //                                                             List<int> viAroundAdjacentFace0,
        //                                                             List<int> viAroundAdjacentFace1, 
        //                                                             int currentFaceIndex,
        //                                                             int adjacentFace0Index,
        //                                                             int adjacentFace1Index,
        //                                                             int startVertex0,
        //                                                             int endVertex0,
        //                                                             int startVertex1,
        //                                                             int endVertex1,
        //                                                             bool isOutOfPrevPrevFERange,
        //                                                             bool isOutOfNextFERange,
        //                                                             out List<List<int>> newFaceVertexOrder,
        //                                                             out int midVertex)
        //{
        //    List<PlanktonXYZ> newPlanktonXYZ = new List<PlanktonXYZ>();
        //    if (!isOutOfPrevPrevFERange && !isOutOfNextFERange)
        //    {
        //        // 理论上不存在这种情况，因为在进行AFace_Operation1的操作时，就已经避免了这种情况的出现
        //        newPlanktonXYZ = null;
        //        newFaceVertexOrder = null;
        //        midVertex = -1;
        //    }
        //    else if(!isOutOfPrevPrevFERange && isOutOfNextFERange)
        //    {
        //        newPlanktonXYZ = ModifyVertice_AFace_Operation2_2(D, viAroundFace, viAroundAdjacentFace0, currentFaceIndex, adjacentFace0Index, startVertex0, endVertex0, endVertex1, out newFaceVertexOrder, out midVertex);
        //    }
        //    else if (isOutOfPrevPrevFERange && !isOutOfNextFERange)
        //    {
        //        newPlanktonXYZ = ModifyVertice_AFace_Operation2_2(D, viAroundFace, viAroundAdjacentFace1, currentFaceIndex, adjacentFace1Index, startVertex1, endVertex1, startVertex0, out newFaceVertexOrder, out midVertex);
        //    }
        //    else
        //    {
        //        newPlanktonXYZ = ModifyVertice_AFace_Operation2_1(D, viAroundFace, currentFaceIndex, startVertex0, endVertex1, out newFaceVertexOrder, out midVertex);
        //    }

        //    return newPlanktonXYZ;
        //}

        //private List<PlanktonXYZ> ModifyVertice_AFace_Operation2_1(PlanktonMesh D,
        //                                                           List<int> viAroundFace,
        //                                                           int currentFaceIndex,
        //                                                           int startVertex0,
        //                                                           int endVertex1,
        //                                                           out List<List<int>> newFaceVertexOrder,
        //                                                           out int midVertex)
        //{
        //    #region 生成新的点
        //    Point3d startPoint = D.Vertices[startVertex0].ToPoint3d();
        //    Point3d endPoint = D.Vertices[endVertex1].ToPoint3d();

        //    Point3d midPoint = (startPoint + endPoint) / 2;

        //    D.Vertices.Add(midPoint);
        //    midVertex = D.Vertices.Count - 1;
        //    #endregion

        //    #region 生成需要修改的Face的Vertex列表
        //    List<int> allVerticeForCurrentFace = new List<int>();
        //    allVerticeForCurrentFace.AddRange(viAroundFace);

        //    int startVertex0InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(startVertex0);
        //    int endVertex1InAllVerticeIndex = allVerticeForCurrentFace.IndexOf(endVertex1);

        //    #region 找到allVerticeForCurrentFace的两个切片
        //    List<int> viFace = new List<int>();
        //    int iter0 = startVertex0InAllVerticeIndex;
        //    while (allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != endVertex1)
        //    {
        //        iter0 = (iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count - 1;
        //        viFace.Add(allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //    }

        //    List<int> viNewFace = new List<int>();
        //    int iter1 = startVertex0InAllVerticeIndex;
        //    while (allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != endVertex1)
        //    {
        //        iter1 = (iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count + 1;
        //        viNewFace.Add(allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //    }

            
        //    viFace.Add(midVertex);
        //    viNewFace.Add(midVertex);
        //    viFace.Insert(0, startVertex0);
        //    viNewFace.Insert(0, startVertex0);
        //    viFace.Reverse();
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonXYZ = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonXYZ.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来的D的Face属性
        //    newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentFaceIndex)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viFace);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }
        //    newFaceVertexOrder.Add(viNewFace);

        //    return newPlanktonXYZ;
        //    #endregion
        //}

        //private List<PlanktonXYZ> ModifyVertice_AFace_Operation2_2(PlanktonMesh D, 
        //                                                           List<int> viAroundFace, 
        //                                                           List<int> viAroundFaceNeedToModify, 
        //                                                           int currentFaceIndex, 
        //                                                           int faceNeedToModify, 
        //                                                           int startVertexOnPublicHForFaceNeedToModify, 
        //                                                           int endVertexOnPublicHForFaceNeedToModify, 
        //                                                           int anotherEndVertex, 
        //                                                           out List<List<int>> newFaceVertexOrder,
        //                                                           out int midVertex)
        //{
        //    #region 生成新的点
        //    Point3d startPoint = D.Vertices[startVertexOnPublicHForFaceNeedToModify].ToPoint3d();
        //    Point3d endPoint = D.Vertices[endVertexOnPublicHForFaceNeedToModify].ToPoint3d();

        //    Point3d midPoint = (startPoint + endPoint) / 2;

        //    D.Vertices.Add(midPoint);
        //    midVertex = D.Vertices.Count - 1;
        //    #endregion

        //    #region 生成需要修改的Face的Vertex列表
        //    List<int> allVerticeForCurrentFace = new List<int>();
        //    allVerticeForCurrentFace.AddRange(viAroundFace);

        //    int endVerticeForCurrentFace = allVerticeForCurrentFace.IndexOf(endVertexOnPublicHForFaceNeedToModify);
        //    allVerticeForCurrentFace.Insert(endVerticeForCurrentFace, midVertex);

        //    int midVertexInAllVerticeIndex = allVerticeForCurrentFace.IndexOf(midVertex);

        //    #region 判断midVertex与anotherEndVertex的关系，由此决定viFace与viNewFace的内容
        //    // true表示anotherEndVertex在midVertex右边（即由midVertex到anotherEndVertex逆时针的间隔数小于顺时针的间隔数）
        //    // false表示anotherEndVertex在midVertex左边（即由midVertex到anotherEndVertex逆时针的间隔数大于顺时针的间隔数）
        //    bool flag = RelationOfTwoVerticesInFaceVertexList(midVertex, anotherEndVertex, allVerticeForCurrentFace);
        //    List<int> viFace = new List<int>();
        //    List<int> viNewFace = new List<int>();
        //    if (flag)
        //    {
        //        #region 找到allVerticeForCurrentFace的两个切片
        //        int iter0 = midVertexInAllVerticeIndex;
        //        while (allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != anotherEndVertex)
        //        {
        //            iter0 = (iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count - 1;
        //            viNewFace.Add(allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //        }

        //        int iter1 = midVertexInAllVerticeIndex;
        //        while (allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != anotherEndVertex)
        //        {
        //            iter1 = (iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count + 1;
        //            viFace.Add(allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //        }

        //        viFace.Insert(0, midVertex);
        //        viNewFace.Insert(0, midVertex);

        //        viNewFace.Reverse();
        //        #endregion
        //    }
        //    else
        //    {
        //        #region 找到allVerticeForCurrentFace的两个切片
        //        int iter0 = midVertexInAllVerticeIndex;
        //        while (allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != anotherEndVertex)
        //        {
        //            iter0 = (iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count - 1;
        //            viFace.Add(allVerticeForCurrentFace[(iter0 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //        }

        //        int iter1 = midVertexInAllVerticeIndex;
        //        while (allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count] != anotherEndVertex)
        //        {
        //            iter1 = (iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count + 1;
        //            viNewFace.Add(allVerticeForCurrentFace[(iter1 + allVerticeForCurrentFace.Count) % allVerticeForCurrentFace.Count]);
        //        }

        //        viFace.Insert(0, midVertex);
        //        viNewFace.Insert(0, midVertex);

        //        viFace.Reverse();
        //        #endregion
        //    }
        //    #endregion

        //    #region 处理viAroundPairBFace0
        //    int startVertex0OnPairBFace0 = viAroundFaceNeedToModify.IndexOf(startVertexOnPublicHForFaceNeedToModify);
        //    viAroundFaceNeedToModify.Insert(startVertex0OnPairBFace0, midVertex);
        //    #endregion
        //    #endregion

        //    #region 生成新的PlanktonMesh
        //    // 转移原来D的Vertex属性
        //    List<PlanktonXYZ> newPlanktonXYZ = new List<PlanktonXYZ>();
        //    for (int i = 0; i < D.Vertices.Count; i++)
        //    {
        //        newPlanktonXYZ.Add(D.Vertices[i].ToXYZ());
        //    }
        //    // 转移原来D的Face属性，同时替换viAroundFace
        //    // viAroundPairFace作为原来的面的替代，viAroundFace作为新生成的面
        //    newFaceVertexOrder = new List<List<int>>();
        //    for (int i = 0; i < D.Faces.Count; i++)
        //    {
        //        if (i == currentFaceIndex)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viFace);
        //        }
        //        else if (i == faceNeedToModify)
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            newFaceVertexOrder[i].AddRange(viAroundFaceNeedToModify);
        //        }
        //        else
        //        {
        //            newFaceVertexOrder.Add(new List<int>());
        //            int[] faceVertexOrder = D.Faces.GetFaceVertices(i);
        //            newFaceVertexOrder[i].AddRange(faceVertexOrder);
        //        }
        //    }
        //    newFaceVertexOrder.Add(viNewFace);
        //    #endregion

        //    return newPlanktonXYZ;
        //}

        //private List<int> ArrangeFaceBoundaryHalfedges(PlanktonMesh D,
        //                                               int currentAFaceIndex,
        //                                               int[] publicHV,
        //                                               List<int[]> abandonedDHV)
        //{
        //    #region 根据publicHV找到publicHalfedge
        //    int publicHalfedge = -1;
        //    int[] hs = D.Faces.GetHalfedges(currentAFaceIndex);
        //    for (int i = 0; i < D.Halfedges.Count; i++)
        //    {
        //        int[] hV = new int[2] { D.Halfedges[i].StartVertex, D.Halfedges.EndVertex(i) };
        //        if (publicHV.Except(hV).ToList().Count == 0
        //            && hs.Contains(i))
        //        {
        //            publicHalfedge = i;
        //            break;
        //        }
        //    }
        //    #endregion

        //    #region 找到adjacentFace的boundaryHalfedges
        //    List<int> boundaryHalfedges = new List<int>();
        //    for (int i = 0; i < hs.Length; i++)
        //    {
        //        if (D.Halfedges[D.Halfedges.GetPairHalfedge(hs[i])].AdjacentFace == -1)
        //        {
        //            boundaryHalfedges.Add(hs[i]);
        //        }
        //    }
        //    #endregion

        //    #region 移除boundaryHalfedges中标记为废弃的halfedge
        //    List<int[]> boundaryHalfedgesV = new List<int[]>();
        //    for (int i = 0; i < boundaryHalfedges.Count; i++)
        //    {
        //        int[] indexPair = new int[2] { D.Halfedges[boundaryHalfedges[i]].StartVertex, D.Halfedges.EndVertex(boundaryHalfedges[i]) };
        //        boundaryHalfedgesV.Add(indexPair);
        //    }

        //    // for倒序遍历删除
        //    for (int i = boundaryHalfedgesV.Count - 1; i >= 0; i--)
        //    {
        //        for (int j = 0; j < abandonedDHV.Count; j++)
        //        {
        //            if (boundaryHalfedgesV[i].Except(abandonedDHV[j]).ToList().Count == 0)
        //            {
        //                boundaryHalfedgesV.RemoveAt(i);
        //                boundaryHalfedges.RemoveAt(i);
        //            }
        //        }
        //    }
        //    #endregion

        //    // 找到nextHalfedge是publicHalfedge的那条半边
        //    int index = -1;
        //    for (int i = 0; i < boundaryHalfedges.Count; i++)
        //    {
        //        if (D.Halfedges[boundaryHalfedges[i]].NextHalfedge == publicHalfedge)
        //        {
        //            index = i;
        //            break;
        //        }
        //    }
        //    // 将第index个元素放在列表首位，其他元素依次排在它后面
        //    int iter = 0;
        //    while (iter < index)
        //    {
        //        int item = boundaryHalfedges[0];
        //        boundaryHalfedges.Add(item);
        //        boundaryHalfedges.RemoveAt(0);
        //        iter++;
        //    }

        //    return boundaryHalfedges;
        //}

        //public bool RelationOfTwoVerticesInFaceVertexList(int vertexToSplitIndex, int targetVertexIndex, List<int> viAroundFaceList)
        //{
        //    int counterClockwiseCount;
        //    int clockwiseCount;
        //    int index0 = 0;
        //    int index1 = 0;
        //    bool flag = false;

        //    for (int i = 0; i < viAroundFaceList.Count; i++)
        //    {

        //        if (viAroundFaceList[i] == vertexToSplitIndex)
        //        {
        //            index0 = i;
        //        }
        //        if (viAroundFaceList[i] == targetVertexIndex)
        //        {
        //            index1 = i;
        //        }
        //    }

        //    if (index0 < index1)
        //    {
        //        counterClockwiseCount = index1 - index0;
        //        clockwiseCount = viAroundFaceList.Count - index1 + index0;
        //    }
        //    else
        //    {
        //        counterClockwiseCount = viAroundFaceList.Count - index0 + index1;
        //        clockwiseCount = index0 - index1;
        //    }

        //    // 如果VertexToSplit和TargetVertex逆时针数的间隔大于顺时针数的间隔，那么VertexToSplit在TargetVertex的左边
        //    // 判断VertexToSplit和TargetVertex谁在谁左边，true表示t在v右边；false表示t在v左边
        //    if (counterClockwiseCount < clockwiseCount)
        //    {
        //        flag = true;
        //    }
        //    else
        //    {
        //        flag = false;
        //    }

        //    return flag;
        //}

        ///// <summary>
        ///// 向drawnFE添加不重复的realisticSiteFaceEdgeSegments
        ///// </summary>
        ///// <param name="realisticSiteFaceEdgeSegments"></param>
        ///// <param name="drawnDHV"></param>
        ///// <param name="drawnFE"></param>
        //private void AddrealisticSiteFaceEdgeSegments(List<FaceEdgeSegment> realisticSiteFaceEdgeSegments, List<int[]> drawnDHV, List<FaceEdgeSegment> drawnFE)
        //{
        //    for (int i = 0; i < realisticSiteFaceEdgeSegments.Count; i++)
        //    {
        //        int[] dHV = realisticSiteFaceEdgeSegments[i].IncludedDVertice.ToArray();
        //        bool hasInDrawn = false;
        //        for (int j = 0; j < drawnDHV.Count; j++)
        //        {
        //            if (drawnDHV[j].Except(dHV).ToList().Count == 0)
        //            {
        //                hasInDrawn = true;
        //                break;
        //            }
        //        }

        //        if (!hasInDrawn)
        //        {
        //            drawnDHV.Add(dHV);
        //            drawnFE.Add(realisticSiteFaceEdgeSegments[i]);
        //        }
        //        else
        //        {
        //            continue;
        //        }
        //    }
        //}

        //private void RemoveOldFE(int[] indexPairToRemove, List<int[]> drawnDHV, List<FaceEdgeSegment> drawnFE)
        //{
        //    int removeIndex = -1;
        //    for (int i = 0; i < drawnDHV.Count; i++)
        //    {
        //        if (drawnDHV[i].Except(indexPairToRemove).ToList().Count == 0)
        //        {
        //            removeIndex = i;
        //            break;
        //        }
        //    }
        //    drawnDHV.RemoveAt(removeIndex);
        //    drawnFE.RemoveAt(removeIndex);
        //}

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            for (int i = 0; i < DebugPolylinesForGreenColorHatch.Count; i++)
            {
                //args.Display.EnableDepthTesting(false);
                args.Display.DrawHatch(DebugPolylinesForGreenColorHatch[i], Color.DarkGreen, Color.Black);
                //args.Display.EnableDepthTesting(true);
            }
            for (int i = 0; i < DebugPolylinesForYellowColorHatch.Count; i++)
            {
                //args.Display.EnableDepthTesting(false);
                args.Display.DrawHatch(DebugPolylinesForYellowColorHatch[i], Color.DarkOrange, Color.Black);
                //args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < RealisticSiteFaceCenter.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(RealisticSiteFaceCenter[i], Color.Black, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            base.DrawViewportMeshes(args);

            for (int i = 0; i < RealisticSiteFaceCenter.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(RealisticSiteFaceCenter[i], Color.Black, Color.White, Color.White);
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
            get { return new Guid("3b884dcf-666e-4908-9a87-5e719ab579d8"); }
        }
    }
}