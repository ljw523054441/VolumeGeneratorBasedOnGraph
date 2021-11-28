using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Plankton;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using VolumeGeneratorBasedOnGraph.Class;

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
        }

        private List<TextDot> RealisticSiteFaceCenter;

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
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            Curve reorganizedBoundaryCurve = null;

            // 因为BoundarySegment没有实现IGH_Goo接口，所以作为GH_Structure无法传递
            // 直接传递List<List<BoundarySegment>>也是不行的，所以用以下方式进行传递
            List<GH_ObjectWrapper> objList = new List<GH_ObjectWrapper>();
            DA.GetDataList(2, objList);

            List<List<FaceEdgeSegment>> subBoundarySegments = new List<List<FaceEdgeSegment>>();
            for (int i = 0; i < objList.Count; i++)
            {
                subBoundarySegments.Add(new List<FaceEdgeSegment>());
                if (objList[i] != null)
                {
                    List<FaceEdgeSegment> currentSubBoundarySegment = objList[i].Value as List<FaceEdgeSegment>;
                    subBoundarySegments[i].AddRange(currentSubBoundarySegment);
                }
            }

            //GH_Structure<GH_Integer> gh_structure = new GH_Structure<GH_Integer>();
            //DataTree<BoundarySegment> sortedSubBoundarySegmentDT = new DataTree<BoundarySegment>();
            //List<List<BoundarySegment>> subBoundarySegments = new List<List<BoundarySegment>>();

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetData<Curve>("ReorganizedBoundary",ref reorganizedBoundaryCurve)
                && subBoundarySegments.Count > 0)
            {
                DualGraphWithHM dualGraphWithHMDP = new DualGraphWithHM(dualGraphWithHM);
                PlanktonMesh D = dualGraphWithHMDP.DualPlanktonMesh;

                List<string> debug = UtilityFunctions.PrintFacesHalfedges(D);
                List<string> debug2 = UtilityFunctions.PrintHalfedgeStartAndEnd(D);

                // 对输入的Curve类型的reorganizedBoundary进行类型转换，转换成Curve类的子类Polyline
                Polyline reorganizedBoundary = null;
                reorganizedBoundaryCurve.TryGetPolyline(out reorganizedBoundary);

                // 记录新生成的面，以及新生成的面和之前的dFace的关系
                List<List<int>> realisticSiteFaceVerticesForAllFace = new List<List<int>>();
                List<List<Point3d>> realisticSiteFacePointsForAllFace = new List<List<Point3d>>();

                List<int> realisticSiteFaceIndex = new List<int>();

                SortedDictionary<int, int> dFace_realisticFace = new SortedDictionary<int, int>();

                // 记录当前最新的顶点的序号
                int newestVertexIndex = D.Vertices.Count - 1;

                int newestHIndex = D.Halfedges.Count - 1;

                List<int> hHasCorrespondFE = new List<int>();
                List<FaceEdgeSegment> feHasH = new List<FaceEdgeSegment>();

                // 找到边缘的面（有连续的三边可以确定），因为都是四边面，所以但凡是有三边的，都是连续的三边
                List<int> bFaceIndexThreeSideDetermined = new List<int>();
                for (int i = 0; i < D.Faces.Count; i++)
                {
                    if (D.Faces.NakedEdgeCount(i) >= 3)
                    {
                        bFaceIndexThreeSideDetermined.Add(i);
                    }
                }
                // 找到边缘的面（有连续的两边可以确定）
                List<int> bFaceIndexTwoSidedDetermined = new List<int>();
                for (int i = 0; i < D.Faces.Count; i++)
                {
                    if (D.Faces.NakedEdgeCount(i) == 2)
                    {
                        List<int> nakedEdges = GetNakedEdges(D,i);
                        if (IsNakedEdgesContinuous(D, nakedEdges))
                        {
                            bFaceIndexTwoSidedDetermined.Add(i);
                        }
                    }
                }

                // 对于每个边缘的有三边可以确定的面，找到公共边的对边
                for (int i = 0; i < bFaceIndexThreeSideDetermined.Count; i++)
                {
                    List<Point3d> realisticSiteFacePoints;
                    List<int> realisticSiteFaceVertice = NakedEdge_3_Operation(D, 
                                                                               subBoundarySegments, 
                                                                               bFaceIndexThreeSideDetermined[i], 
                                                                               newestVertexIndex,
                                                                               newestHIndex,
                                                                               out realisticSiteFacePoints,
                                                                               hHasCorrespondFE,
                                                                               feHasH);

                    realisticSiteFaceVerticesForAllFace.Add(new List<int>());
                    realisticSiteFacePointsForAllFace.Add(new List<Point3d>());
                    realisticSiteFaceVerticesForAllFace[i].AddRange(realisticSiteFaceVertice);
                    realisticSiteFacePointsForAllFace[i].AddRange(realisticSiteFacePoints);

                    realisticSiteFaceIndex.Add(bFaceIndexThreeSideDetermined[i]);
                }

                for (int i = 0; i < bFaceIndexTwoSidedDetermined.Count; i++)
                {
                    List<Point3d> realisticSiteFacePoints;
                    List<int> realisticSiteFaceVertice = NakeEdge_2_Operation(D, 
                                                                              subBoundarySegments, 
                                                                              bFaceIndexTwoSidedDetermined[i],
                                                                              newestVertexIndex,
                                                                              newestHIndex,
                                                                              out realisticSiteFacePoints,
                                                                              hHasCorrespondFE,
                                                                              feHasH);

                    realisticSiteFaceVerticesForAllFace.Add(new List<int>());
                    realisticSiteFacePointsForAllFace.Add(new List<Point3d>());
                    realisticSiteFaceVerticesForAllFace[i + bFaceIndexThreeSideDetermined.Count].AddRange(realisticSiteFaceVertice);
                    realisticSiteFacePointsForAllFace[i + bFaceIndexThreeSideDetermined.Count].AddRange(realisticSiteFacePoints);

                    realisticSiteFaceIndex.Add(bFaceIndexTwoSidedDetermined[i]);
                }






                // 结果输出部分
                List<Polyline> debugPolylines = new List<Polyline>();
                for (int i = 0; i < realisticSiteFacePointsForAllFace.Count; i++)
                {
                    List<Point3d> list = new List<Point3d>();
                    list.AddRange(realisticSiteFacePointsForAllFace[i]);
                    list.Add(realisticSiteFacePointsForAllFace[i][0]);
                    Polyline polyline = new Polyline(list);
                    debugPolylines.Add(polyline);
                }
                DA.SetDataList("DeBug Polyline", debugPolylines);

                RealisticSiteFaceCenter.Clear();
                for (int i = 0; i < debugPolylines.Count; i++)
                {
                    Point3d center = debugPolylines[i].CenterPoint();
                    TextDot fCenterTextDot = new TextDot(realisticSiteFaceIndex[i].ToString(), center);
                    RealisticSiteFaceCenter.Add(fCenterTextDot);
                }


            }
        }

        private List<int> GetNakedEdges(PlanktonMesh D, int currentFaceIndex)
        {
            List<int> nakedEdges = new List<int>();
            for (int i = 0; i < D.Faces.GetHalfedges(currentFaceIndex).Length; i++)
            {
                int pairH = D.Halfedges.GetPairHalfedge(D.Faces.GetHalfedges(currentFaceIndex)[i]);
                if (D.Halfedges[pairH].AdjacentFace == -1)
                {
                    nakedEdges.Add(D.Faces.GetHalfedges(currentFaceIndex)[i]);
                }
            }
            return nakedEdges;
        }

        private bool IsNakedEdgesContinuous(PlanktonMesh D, List<int> nakedEdges)
        {
            List<int> startVertice = new List<int>();
            List<int> endVertice = new List<int>();
            int nakedEdgesCount = nakedEdges.Count;

            for (int i = 0; i < nakedEdges.Count; i++)
            {
                int start = D.Halfedges[nakedEdges[i]].StartVertex;
                int end = D.Halfedges.EndVertex(nakedEdges[i]);

                startVertice.Add(start);
                endVertice.Add(end);
            }

            List<int> intersect = startVertice.Intersect(endVertice).ToList();
            if (intersect.Count == nakedEdgesCount - 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        private List<int> NakedEdge_3_Operation(PlanktonMesh D,
                                                List<List<FaceEdgeSegment>> subBoundarySegments,
                                                int currentBFaceIndexThreeSideDetermined,
                                                int newestVertexIndex,
                                                int newestHIndex,
                                                // out List<int> realisticSiteFaceVertices,
                                                // out int realisticSiteFaceIndex,
                                                out List<Point3d> realisticSiteFacePoints,
                                                List<int> hHasCorrespondFE,
                                                List<FaceEdgeSegment> feHasH)
        {
            List<int> realisticSiteFaceVertices = new List<int>();
            realisticSiteFacePoints = new List<Point3d>();
            List<FaceEdgeSegment> realisticSiteFE = new List<FaceEdgeSegment>();

            #region 找到公共边的index
            int publicHalfedge = -1;
            int[] hs = D.Faces.GetHalfedges(currentBFaceIndexThreeSideDetermined);
            for (int i = 0; i < hs.Length; i++)
            {
                if (D.Halfedges[D.Halfedges.GetPairHalfedge(hs[i])].AdjacentFace != -1)
                {
                    publicHalfedge = hs[i];
                }
            }
            #endregion

            #region 从公共边的nextHalfedge开始，依次添加这个面的所有Halfedge
            // 按顺序存储这个面的halfedge
            List<int> hIndexList = new List<int>();
            int iter = 0;
            int currentEdge = publicHalfedge;
            do
            {
                hIndexList.Add(D.Halfedges[currentEdge].NextHalfedge);
                currentEdge = hIndexList.Last();
                iter++;
            } while (iter < D.Faces.GetHalfedges(currentBFaceIndexThreeSideDetermined).Length);
            #endregion


            #region 按照hIndexList中的顺序，添加现有的relatedSubBoundarySegments
            // 按顺序存储这个面目前有的BoundarySegment
            // 注意，此时公共边没有对应的FaceEdgeSegment
            List<FaceEdgeSegment> relatedSubBoundarySegments = new List<FaceEdgeSegment>();
            for (int i = 0; i < hIndexList.Count - 1; i++)
            {
                for (int j = 0; j < subBoundarySegments.Count; j++)
                {
                    for (int k = 0; k < subBoundarySegments[j].Count; k++)
                    {
                        if (subBoundarySegments[j][k].HIndex == hIndexList[i])
                        {
                            relatedSubBoundarySegments.Add(subBoundarySegments[j][k]);
                        }
                    }
                }
            }
            #endregion

            #region 按照hIndexList中的顺序（即现有的relatedSubBoundarySegments的顺序），添加Vertex对应的Point位置
            // 注意，此时的realisticSiteFacePoints并不是最终的结果
            // 注意，segment方向与半边方向是反的
            realisticSiteFacePoints.Add(relatedSubBoundarySegments[0].To);
            for (int i = 0; i < relatedSubBoundarySegments.Count; i++)
            {
                realisticSiteFacePoints.Add(relatedSubBoundarySegments[i].From);
            }
            #endregion

            #region 判断每次faceEdge转折的方向是否是我们想要的（右手定则）

            // 计算时[0]都被直接添加，每次循环都是添加next
            realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[0]].StartVertex);
            realisticSiteFE.Add(relatedSubBoundarySegments[0]);
            hHasCorrespondFE.Add(hIndexList[0]);
            // 所以i从1开始
            for (int i = 1; i < hIndexList.Count; i++)
            {
                // currPoint是realisticSiteFacePoints[i]
                Vector3d prevV = new Vector3d(realisticSiteFacePoints[i] - realisticSiteFacePoints[i - 1]);
                Vector3d nextV = new Vector3d(realisticSiteFacePoints[(i + 1) % realisticSiteFacePoints.Count] - realisticSiteFacePoints[i]);
                Vector3d crossProduct = Vector3d.CrossProduct(prevV, nextV);

                // 如果是，那么
                if (crossProduct.Z < 0)
                {
                    // 如果当前转折方向是我们想要的（右手定则），并且当前这个半边的序号，在hHasCorrespondFE中出现过的话（即该半边生成过对应的FE）
                    // 判断已经有对应FaceEdgeSegment的半边序号列表中，是否存在该半边的序号，即，判断该半边是否已经绘制过对应的FaceEdgeSegment
                    if (hHasCorrespondFE.Contains(hIndexList[i]))
                    {
                        FaceEdgeSegment alreadyExistFE = feHasH[hHasCorrespondFE.IndexOf(hIndexList[i])];

                        FaceEdgeSegment feToMove = relatedSubBoundarySegments[i - 1];
                        Line lineToMove = feToMove.Line;

                        Transform move = Transform.Translation(new Vector3d(realisticSiteFacePoints[((i - 2) + realisticSiteFacePoints.Count) % realisticSiteFacePoints.Count]
                                                                            - realisticSiteFacePoints[((i - 3) + realisticSiteFacePoints.Count) % realisticSiteFacePoints.Count]));

                        // 求交点
                        Line newLine = new Line(lineToMove.From, lineToMove.To);
                        newLine.Transform(move);

                        double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                        double t1;
                        double t2;
                        if (Intersection.LineLine(alreadyExistFE.Line, newLine, out t1, out t2, tolerance, false))
                        {
                            Point3d intersectPoint1 = alreadyExistFE.Line.PointAt(t1);
                            Point3d intersectPoint2 = newLine.PointAt(t2);

                            Line lineCoincideWithAlreadyExistFE = new Line(realisticSiteFacePoints[i - 1], intersectPoint1);
                            // string label1 = alreadyExistFE.Label + ":" + realisticSiteFaceVertices[i - 1].ToString() + "-" + (newestVertexIndex + 1).ToString();
                            string label1 = string.Format("h:{0},v:{1},{2}", hIndexList[i], realisticSiteFaceVertices[i - 1], newestVertexIndex + 1);
                            List<int> includedDVertice1 = new List<int>();
                            includedDVertice1.Add(realisticSiteFaceVertices[i - 1]);
                            includedDVertice1.Add(newestVertexIndex + 1);

                            FaceEdgeSegment newFE1 = new FaceEdgeSegment(lineCoincideWithAlreadyExistFE, label1, includedDVertice1, hIndexList[i]);

                            Line lineForNewFE = new Line(intersectPoint2, realisticSiteFacePoints[i - 3]);
                            // string label2 = (newestVertexIndex + 1).ToString() + "-" + realisticSiteFaceVertices[i - 3].ToString();
                            string label2 = string.Format("h:{0},v:{1},{2}", newestHIndex + 1, newestVertexIndex + 1, realisticSiteFaceVertices[i - 3]);
                            List<int> includedDVertice2 = new List<int>();
                            includedDVertice2.Add(newestVertexIndex + 1);
                            includedDVertice2.Add(realisticSiteFaceVertices[i - 3]);

                            FaceEdgeSegment newFE2 = new FaceEdgeSegment(lineForNewFE, label2, includedDVertice2, newestHIndex + 1);

                            realisticSiteFE.Add(newFE1);
                            // realisticSiteFaceBS.Add(newFE2);
                            // 不用再向hHasCorrespondFE中添加当前的hIndex了
                            //hHasCorrespondFE.Add()

                            newestVertexIndex = newestVertexIndex + 1;
                        }
                        else
                        {
                            throw new Exception(string.Format("无法找到交点，因为{0}与另一条边两条线是平行的", alreadyExistFE));
                        }
                    }
                    // 如果当前转折方向是我们想要的（右手定则），并且当前这个半边的序号，没有在hHasCorrespondFE中出现过的话（即该半边没有生成过对应的FE）
                    else
                    {
                        // 向FaceVertices中添加当前halfedge的起点
                        realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[i]].StartVertex);

                        // 向FaceBS中添加当前的faceEdgeSegment
                        if (i < relatedSubBoundarySegments.Count)
                        {
                            realisticSiteFE.Add(relatedSubBoundarySegments[i]);
                        }
                        else
                        {
                            Line newLine = new Line(realisticSiteFacePoints[i], realisticSiteFacePoints[(i + 1) % realisticSiteFacePoints.Count]);
                            string label = string.Format("h:{0},v:{1},{2}", hIndexList[i], realisticSiteFaceVertices[i], realisticSiteFaceVertices[(i + 1) % realisticSiteFaceVertices.Count]);
                            List<int> includeDV = new List<int>();
                            includeDV.Add(realisticSiteFaceVertices[i]);
                            includeDV.Add(realisticSiteFaceVertices[(i + 1) % realisticSiteFaceVertices.Count]);

                            FaceEdgeSegment newFE = new FaceEdgeSegment(newLine, label, includeDV, hIndexList[i]);
                            realisticSiteFE.Add(newFE);
                        }
                        
                        // 向已经有对应FaceEdgeSegment的半边序号列表中，添加该半边的序号
                        hHasCorrespondFE.Add(hIndexList[i]);

                        // 没有新的Vertex生成，没有新的HIndex生成
                    }
                
                }
                // 如果转折不是我们想要的（右手定则），那么
                else
                {
                    // 如果当前转折方向不是我们想要的（右手定则），并且当前这个半边的序号，在hHasCorrespondFE中出现过的话（即该半边生成过对应的FE）
                    // 判断已经有对应FaceEdgeSegment的半边序号列表中，是否存在该半边的序号，即，判断该半边是否已经绘制过对应的FaceEdgeSegment
                    if (hHasCorrespondFE.Contains(hIndexList[i]))
                    {

                    }
                    // 如果当前转折方向不是我们想要的（右手定则），并且当前这个半边的序号，没有在hHasCorrespondFE中出现过的话（即该半边没有生成过对应的FE）
                    else
                    {
                        realisticSiteFaceVertices.Add(newestVertexIndex + 1);

                        Point3d prevPrevPoint = realisticSiteFacePoints[i - 2];
                        Point3d newNextPoint = prevPrevPoint + prevV;

                        realisticSiteFacePoints[(i + 1) % realisticSiteFacePoints.Count] = newNextPoint;

                        Line nextLine = new Line(newNextPoint, realisticSiteFacePoints[i]);
                        // string label = currentBFaceIndexThreeSideDetermined.ToString() + "-" + hIndexList[i];
                        string label = string.Format("h:{0},v:{1},{2}", hIndexList[i], realisticSiteFaceVertices[i], newestVertexIndex + 1);
                        List<int> includedDVertice = new List<int>();
                        includedDVertice.Add(realisticSiteFaceVertices[i]);
                        includedDVertice.Add(newestVertexIndex + 1);

                        FaceEdgeSegment nextFE = new FaceEdgeSegment(nextLine, label, includedDVertice, hIndexList[i]);

                        realisticSiteFE.Add(nextFE);
                        hHasCorrespondFE.Add(hIndexList[i]);

                        // 有新的Vertex生成，故newestVertexIndex++
                        newestVertexIndex++;
                        //没有新的HIndex生成
                    }
                }
            }

            #endregion



            feHasH.AddRange(realisticSiteFE);

            return realisticSiteFaceVertices;
            // 先暂时用i作为realisticFace的序号
            // dFace_realisticFace.Add(bFaceIndexThreeSideDetermined[i], i);
            // realisticSiteFaceIndex = i
        }

        private List<int> NakeEdge_2_Operation(PlanktonMesh D,
                                               List<List<FaceEdgeSegment>> subBoundarySegments,
                                               int currentBFaceIndexTwoSideDetermined,
                                               int newestVertexIndex,
                                               int newestHIndex,
                                               out List<Point3d> realisticSiteFacePoints,
                                               List<int> hHasCorrespondFE,
                                               List<FaceEdgeSegment> feHasH)
        {
            List<int> realisticSiteFaceVertices = new List<int>();
            realisticSiteFacePoints = new List<Point3d>();
            List<FaceEdgeSegment> realisticSiteFE = new List<FaceEdgeSegment>();

            #region 找到公共边的index
            List<int> publicHalfedges = new List<int>();
            int[] hs = D.Faces.GetHalfedges(currentBFaceIndexTwoSideDetermined);
            for (int i = 0; i < hs.Length; i++)
            {
                if (D.Halfedges[D.Halfedges.GetPairHalfedge(hs[i])].AdjacentFace != -1)
                {
                    publicHalfedges.Add(hs[i]);
                }
            }
            #endregion

            #region 从最后一个公共边的nextHalfedge开始，依次添加这个面的所有Halfedge
            // 按顺序存储这个面的halfedge
            List<int> hIndexList = new List<int>();
            int iter = 0;
            int currentEdge;
            // 找到最后一个公共边的index，作为currentEdge
            if (D.Halfedges[publicHalfedges[0]].StartVertex == D.Halfedges.EndVertex(publicHalfedges[1]))
            {
                currentEdge = publicHalfedges[0];
            }
            else
            {
                currentEdge = publicHalfedges[1];
            }

            do
            {
                hIndexList.Add(D.Halfedges[currentEdge].NextHalfedge);
                currentEdge = hIndexList.Last();
                iter++;
            } while (iter < D.Faces.GetHalfedges(currentBFaceIndexTwoSideDetermined).Length);
            #endregion

            #region 按照hIndexList中的顺序，添加现有的relatedSubBoundarySegments
            // 按顺序存储这个面目前有的BoundarySegment
            // 注意，此时公共边没有对应的FaceEdgeSegment
            List<FaceEdgeSegment> relatedSubBoundarySegments = new List<FaceEdgeSegment>();
            for (int i = 0; i < hIndexList.Count - 1; i++)
            {
                for (int j = 0; j < subBoundarySegments.Count; j++)
                {
                    for (int k = 0; k < subBoundarySegments[j].Count; k++)
                    {
                        if (subBoundarySegments[j][k].HIndex == hIndexList[i])
                        {
                            relatedSubBoundarySegments.Add(subBoundarySegments[j][k]);
                        }
                    }
                }
            }
            #endregion

            #region 按照hIndexList中的顺序（即现有的relatedSubBoundarySegments的顺序），添加Vertex对应的Point位置
            realisticSiteFacePoints.Add(relatedSubBoundarySegments[0].To);
            for (int i = 0; i < relatedSubBoundarySegments.Count; i++)
            {
                realisticSiteFacePoints.Add(relatedSubBoundarySegments[i].From);
            }
            #endregion

            #region 将能够确定的FaceEdge的相关信息，储存好
            realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[0]].StartVertex);
            for (int i = 0; i < relatedSubBoundarySegments.Count; i++)
            {
                realisticSiteFaceVertices.Add(D.Halfedges.EndVertex(hIndexList[i]));
                realisticSiteFE.Add(relatedSubBoundarySegments[i]);
                hHasCorrespondFE.Add(hIndexList[i]);
            }
            #endregion

            // 目前relatedSubBoundarySegments只有两个元素。hIndexList[relatedSubBoundarySegments.Count]是它的下一条边
            bool flag1 = hHasCorrespondFE.Contains(hIndexList[relatedSubBoundarySegments.Count]);
            bool flag3 = hHasCorrespondFE.Contains(D.Halfedges.GetPairHalfedge(hIndexList[relatedSubBoundarySegments.Count]));
            bool flag2 = hHasCorrespondFE.Contains(hIndexList[relatedSubBoundarySegments.Count + 1]);
            bool flag4 = hHasCorrespondFE.Contains(D.Halfedges.GetPairHalfedge(hIndexList[relatedSubBoundarySegments.Count + 1]));

            //如果两条公共边都没有在hHasCorrespondFE中出现过的话（表示这两条公共边都没有绘制过对应的FaceEdgeSegment，即这两条半边都没有生成过对应的FE）
            if (!(flag1||flag3) && !(flag2||flag4))
            {
                Point3d point = realisticSiteFacePoints[2];
                Vector3d move = new Vector3d(realisticSiteFacePoints[0] - realisticSiteFacePoints[1]);
                Point3d newPoint = point + move;

                realisticSiteFacePoints.Add(newPoint);
                realisticSiteFaceVertices.Add(D.Halfedges.EndVertex(hIndexList[2]));

                // 注意构造FE的Line，首尾点顺序跟halfedge顺序应该是反的
                Line prev = new Line(newPoint, realisticSiteFacePoints[2]);
                string prevLabel = string.Format("h:{0},v:{1},{2}", hIndexList[2], realisticSiteFaceVertices[2], realisticSiteFaceVertices[3]);
                List<int> includeDV1 = new List<int>();
                includeDV1.Add(realisticSiteFaceVertices[2]);
                includeDV1.Add(realisticSiteFaceVertices[3]);
                FaceEdgeSegment prevFE = new FaceEdgeSegment(prev, prevLabel, includeDV1, hIndexList[2]);
                // 注意构造FE的Line，首尾点顺序跟halfedge顺序应该是反的
                Line next = new Line(realisticSiteFacePoints[0], newPoint);
                string nextLabel = string.Format("h:{0},v:{1},{2}", hIndexList[3], realisticSiteFaceVertices[3], realisticSiteFaceVertices[0]);
                List<int> includeDV2 = new List<int>();
                includeDV2.Add(realisticSiteFaceVertices[3]);
                includeDV2.Add(realisticSiteFaceVertices[0]);
                FaceEdgeSegment nextFE = new FaceEdgeSegment(next, nextLabel, includeDV2, hIndexList[3]);

                realisticSiteFE.Add(prevFE);
                realisticSiteFE.Add(nextFE);

                hHasCorrespondFE.Add(hIndexList[2]);
                hHasCorrespondFE.Add(hIndexList[3]);
            }
            // 如果两条中有一条出现，另一条没出现，那么求交点
            else if (((flag1||flag3) && !(flag2||flag4)) || (!(flag1||flag3) && (flag2||flag4)))
            {
                FaceEdgeSegment alreadyExistFE;
                FaceEdgeSegment feToMove;
                Line lineToMove;
                int vIndex1;
                Point3d point1;
                int vIndex2;
                Point3d point2;
                Transform move;
                if (((flag1 || flag3) && !(flag2 || flag4)))
                {
                    int hIndex = -1;
                    if (flag1)
                    {
                        hIndex = hIndexList[relatedSubBoundarySegments.Count];
                    }
                    if (flag3)
                    {
                        hIndex = D.Halfedges.GetPairHalfedge(hIndexList[relatedSubBoundarySegments.Count]);
                    }

                    alreadyExistFE = feHasH[hHasCorrespondFE.IndexOf(hIndex)];
                    feToMove = relatedSubBoundarySegments[1];
                    lineToMove = feToMove.Line;
                    vIndex1 = realisticSiteFaceVertices[2];
                    point1 = realisticSiteFacePoints[2];
                    vIndex2 = realisticSiteFaceVertices[0];
                    point2 = realisticSiteFacePoints[0];
                    move = Transform.Translation(new Vector3d(realisticSiteFacePoints[0] - realisticSiteFacePoints[1]));
                }
                else
                {
                    int hIndex = -1;
                    if (flag2)
                    {
                        hIndex = hIndexList[relatedSubBoundarySegments.Count + 1];
                    }
                    if (flag4)
                    {
                        hIndex = D.Halfedges.GetPairHalfedge(hIndexList[relatedSubBoundarySegments.Count + 1]);
                    }

                    alreadyExistFE = feHasH[hHasCorrespondFE.IndexOf(hIndex)];
                    feToMove = relatedSubBoundarySegments[0];
                    lineToMove = feToMove.Line;
                    vIndex1 = realisticSiteFaceVertices[0];
                    point1 = realisticSiteFacePoints[0];
                    vIndex2 = realisticSiteFaceVertices[2];
                    point2 = realisticSiteFacePoints[2];
                    move = Transform.Translation(new Vector3d(realisticSiteFacePoints[2] - realisticSiteFacePoints[1]));
                }

                // 求交点
                Line newLine = new Line(lineToMove.From, lineToMove.To);
                newLine.Transform(move);

                double tolerance = Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                double t1;
                double t2;
                if (Intersection.LineLine(alreadyExistFE.Line, newLine ,out t1 ,out t2 , tolerance, false))
                {
                    Point3d intersectPoint1 = alreadyExistFE.Line.PointAt(t1);
                    Point3d intersectPoint2 = alreadyExistFE.Line.PointAt(t2);

                    Line lineCoincideWithAlreadyExistFE = new Line(point1, intersectPoint1);
                    string label1 = string.Format("h:{0},v:{1},{2} OriginLabel:{3}", hIndexList[relatedSubBoundarySegments.Count + 1], vIndex1, newestVertexIndex + 1, alreadyExistFE.Label);
                    List<int> includedDV1 = new List<int>();
                    includedDV1.Add(vIndex1);
                    includedDV1.Add(newestVertexIndex + 1);

                    FaceEdgeSegment newFE1 = new FaceEdgeSegment(lineCoincideWithAlreadyExistFE, label1, includedDV1, alreadyExistFE.HIndex);

                    Line lineForNewFE = new Line(intersectPoint2, point2);
                    string label2 = string.Format("h:{0},v:{1},{2} OriginLabel:{3}", newestHIndex + 1, newestVertexIndex + 1, vIndex2, alreadyExistFE.Label);
                    List<int> includedDV2 = new List<int>();
                    includedDV2.Add(newestVertexIndex + 1);
                    includedDV2.Add(vIndex2);

                    FaceEdgeSegment newFE2 = new FaceEdgeSegment(lineForNewFE, label2, includedDV2, newestHIndex + 1);

                    realisticSiteFacePoints.Add(intersectPoint1);
                    realisticSiteFaceVertices.Add(newestVertexIndex + 1);
                    realisticSiteFE.Add(newFE1);
                    realisticSiteFE.Add(newFE2);

                    hHasCorrespondFE.Add(newestHIndex + 1);


                    newestVertexIndex++;
                    newestHIndex++;
                }

                
            }
            // 如果都出现了，那么直接找两个FE的公共点即可
            else
            {

            }

            feHasH.AddRange(realisticSiteFE);

            return realisticSiteFaceVertices;
            //Point3d pointToMove = realisticSiteFacePoints.Last();
            //Vector3d move = new Vector3d(realisticSiteFacePoints[realisticSiteFacePoints.Count - 3] - realisticSiteFacePoints[realisticSiteFacePoints.Count - 2]);

            //Point3d newNextPoint = pointToMove + move;
            //realisticSiteFacePoints.Add(newNextPoint);


            // 计算时[0]都被直接添加，每次循环都是添加next
            //realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[0]].StartVertex);
            //realisticSiteFaceFE.Add(relatedSubBoundarySegments[0]);
            //hHasCorrespondFE.Add(hIndexList[0]);
            //// 所以i从1开始
            //for (int i = 1; i < hIndexList.Count; i++)





            //realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[0]].StartVertex);
            //for (int i = 1; i < realisticSiteFacePoints.Count; i++)
            //{
            //    Vector3d curr = new Vector3d(realisticSiteFacePoints[i] - realisticSiteFacePoints[i - 1]);
            //    Vector3d next = new Vector3d(realisticSiteFacePoints[(i + 1) % realisticSiteFacePoints.Count] - realisticSiteFacePoints[i]);
            //    Vector3d crossProduct = Vector3d.CrossProduct(curr, next);

            //    if (crossProduct.Z < 0)
            //    {
            //        realisticSiteFaceVertices.Add(D.Halfedges[hIndexList[i]].StartVertex);
            //    }
            //    else
            //    {
            //        realisticSiteFaceVertices.Add(newestVertexIndex + 1);

            //        //Point3d prevPrevPoint = realisticSiteFacePoints[i - 2];
            //        //Point3d newNextPoint = prevPrevPoint + curr;

            //        //realisticSiteFacePoints[()]
            //    }
            //}
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            for (int i = 0; i < RealisticSiteFaceCenter.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(RealisticSiteFaceCenter[i], Color.ForestGreen, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            base.DrawViewportMeshes(args);

            for (int i = 0; i < RealisticSiteFaceCenter.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(RealisticSiteFaceCenter[i], Color.ForestGreen, Color.White, Color.White);
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