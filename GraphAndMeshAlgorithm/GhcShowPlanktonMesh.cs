using Grasshopper.Kernel;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcShowPlanktonMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ShowPlanktonMesh class.
        /// </summary>
        public GhcShowPlanktonMesh()
          : base("ShowPlanktonMesh", "ShowPlanktonMesh",
              "显示PlanktonMesh",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
            InnerNodeTextDot = new List<TextDot>();
            OuterNodeTextDot = new List<TextDot>();

            FaceTextDot = new List<TextDot>();

            LineForHalfEdges = new List<List<Line>>();
            TextDotForHalfEdges = new List<List<TextDot>>();
        }

        private List<TextDot> InnerNodeTextDot;
        private List<TextDot> OuterNodeTextDot;

        private List<TextDot> FaceTextDot;

        private List<List<Line>> LineForHalfEdges;
        private List<List<TextDot>> TextDotForHalfEdges;

        private double ScreenSize;
        private double RelativeSize;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            // pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);
            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

            pManager.AddNumberParameter("OffsetDistance", "OD", "Offset的距离", GH_ParamAccess.item);
            pManager.AddNumberParameter("ScreenSize", "SS", "", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("RelativeSize", "RS", "", GH_ParamAccess.item, 1.0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //// 全局参数传递
            //GlobalParameter globalParameter = new GlobalParameter();
            //DA.GetData("GlobalParameter", ref globalParameter);
            //int innerNodeCount = globalParameter.VolumeNodeCount;
            //int outerNodeCount = globalParameter.BoundaryNodeCount;

            List<Point3d> innerPoints = new List<Point3d>();
            List<Point3d> outerPoints = new List<Point3d>();
            List<int> innerPointsIndex = new List<int>();
            List<int> outerPointsIndex = new List<int>();

            

            PlanktonMesh P = new PlanktonMesh();

            double distance = 0.0;

            double screenSize = 0.0;
            double relativeSize = 0.0;

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh", ref P)
                && DA.GetData<double>("OffsetDistance", ref distance))
            {
                PlanktonMesh PDeepCopy = new PlanktonMesh(P);
                
                InnerNodeTextDot.Clear();
                OuterNodeTextDot.Clear();
                FaceTextDot.Clear();
                LineForHalfEdges.Clear();
                TextDotForHalfEdges.Clear();

                DA.GetData("ScreenSize", ref screenSize);
                DA.GetData("RelativeSize", ref relativeSize);

                ScreenSize = screenSize;
                RelativeSize = relativeSize;


                // 区分和记录innerPoints和outerPoints的位置和序号
                for (int i = 0; i < PDeepCopy.Vertices.Count; i++)
                {
                    if (PDeepCopy.Vertices.IsBoundary(i))
                    {
                        outerPoints.Add(RhinoSupport.ToPoint3d(PDeepCopy.Vertices[i]));
                        outerPointsIndex.Add(i);
                    }
                    else
                    {
                        innerPoints.Add(RhinoSupport.ToPoint3d(PDeepCopy.Vertices[i]));
                        innerPointsIndex.Add(i);
                    }
                }

                // 对innerNode设置用于可视化的TextDot
                for (int i = 0; i < innerPoints.Count; i++)
                {
                    TextDot textDot = new TextDot(string.Format("{0}", innerPointsIndex[i]), innerPoints[i]);
                    InnerNodeTextDot.Add(textDot);
                }
                // 对outerNode设置用于可视化的TextDot
                for (int i = 0; i < outerPoints.Count; i++)
                {
                    TextDot textDot = new TextDot(string.Format("{0}", outerPointsIndex[i]), outerPoints[i]);
                    OuterNodeTextDot.Add(textDot);
                }

                // 记录和存储face的中心，并设置用于可视化的TextDot
                for (int i = 0; i < PDeepCopy.Faces.Count; i++)
                {
                    Point3d center = RhinoSupport.ToPoint3d(PDeepCopy.Faces.GetFaceCenter(i));
                    TextDot textDot = new TextDot(string.Format("{0}", i), center);
                    FaceTextDot.Add(textDot);
                }

                // 将半边转化为箭头
                //for (int i = 0; i < P.Halfedges.Count; i++)
                //{
                //    Point3d startPoint = RhinoSupport.ToPoint3d(P.Vertices[P.Halfedges[i].StartVertex]);
                //    Point3d endPoint = RhinoSupport.ToPoint3d(P.Vertices[P.Halfedges[P.Halfedges[i].NextHalfedge].StartVertex]);

                //    Point3d midPoint = (startPoint + endPoint) / 2;

                //    Vector3d v = new Vector3d(endPoint - startPoint);
                //    Vector3d moveV = Vector3d.CrossProduct(v, new Vector3d(0, 0, -1));
                //    moveV.Unitize();
                //    moveV = moveV * 0.1;

                //    v = v * 0.3;
                //}

                // 将原来的PlanktonMesh变为Polyline，并记录每个面上的半边关系
                // List<List<int>> allHalfEdgesIndex;
                // List<Polyline> polylines = ToPolylines(P, out allHalfEdgesIndex);
                List<List<int>> allHalfEdgesIndex = GetHalfedgeIndexs(PDeepCopy);
                List<Polyline> polylines = RhinoSupport.ToPolylines(PDeepCopy).ToList();

                // offset
                List<Polyline> offsetPolylines = new List<Polyline>();
                for (int i = 0; i < polylines.Count; i++)
                {
                    PolylineCurve polylineCurve = new PolylineCurve(polylines[i]);
                    Curve[] offsets = polylineCurve.Offset(Plane.WorldXY, - distance, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance, CurveOffsetCornerStyle.None);
                    Curve offset = offsets[0];
                    Polyline offsetPolyline;
                    offset.TryGetPolyline(out offsetPolyline);
                    offsetPolylines.Add(offsetPolyline);
                }

                // 生成每个halfEdge所代表的Line
                List<List<Line>> lineForHalfEdges = new List<List<Line>>();
                for (int i = 0; i < offsetPolylines.Count; i++)
                {
                    lineForHalfEdges.Add(new List<Line>());
                    lineForHalfEdges[i].AddRange(offsetPolylines[i].GetSegments());
                }
                LineForHalfEdges = lineForHalfEdges;

                // 生成能够代表每个halfEdge的TextDot
                List<List<TextDot>> textDotForHalfEdges = new List<List<TextDot>>();
                for (int i = 0; i < lineForHalfEdges.Count; i++)
                {
                    textDotForHalfEdges.Add(new List<TextDot>());
                    for (int j = 0; j < lineForHalfEdges[i].Count; j++)
                    {
                        TextDot textDot = new TextDot(string.Format("{0}", allHalfEdgesIndex[i][j]), lineForHalfEdges[i][j].PointAt(0.5));
                        textDotForHalfEdges[i].Add(textDot);
                    }
                }
                TextDotForHalfEdges = textDotForHalfEdges;
            }
        }

        public List<List<int>> GetHalfedgeIndexs(PlanktonMesh P)
        {
            List<List<int>> allHalfEdgesIndex = new List<List<int>>();

            for (int i = 0; i < P.Faces.Count; i++)
            {
                allHalfEdgesIndex.Add(new List<int>());
                int[] halfEdgesIndex = P.Faces.GetHalfedges(i);
                allHalfEdgesIndex[i].AddRange(halfEdgesIndex);
            }

            return allHalfEdgesIndex;
        }

        public List<Polyline> ToPolylines(PlanktonMesh P, out List<List<int>> allHalfEdgesIndex)
        {
            List<Polyline> polylines = new List<Polyline>();
            allHalfEdgesIndex = new List<List<int>>();

            for (int i = 0; i < P.Faces.Count; i++)
            {
                allHalfEdgesIndex.Add(new List<int>());
                int[] halfEdgesIndex = P.Faces.GetHalfedges(i);
                allHalfEdgesIndex[i].AddRange(halfEdgesIndex);

                Polyline polyline = new Polyline();
                for (int j = 0; j < halfEdgesIndex.Length; j++)
                {
                    Point3d point = RhinoSupport.ToPoint3d(P.Vertices[P.Halfedges[halfEdgesIndex[j]].StartVertex]);
                    polyline.Add(point);
                }
                polyline.Add(polyline[0]);

                polylines.Add(polyline);
            }

            return polylines;
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // 屏蔽掉电池原本的预览
            // base.DrawViewportWires(args);

            for (int i = 0; i < LineForHalfEdges.Count; i++)
            {
                for (int j = 0; j < LineForHalfEdges[i].Count; j++)
                {
                    args.Display.DrawArrow(LineForHalfEdges[i][j], Color.DodgerBlue, ScreenSize, RelativeSize);
                }
            }

            for (int i = 0; i < InnerNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(InnerNodeTextDot[i], Color.Black, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < OuterNodeTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(OuterNodeTextDot[i], Color.Gray, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < FaceTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(FaceTextDot[i], Color.Red, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }

            for (int i = 0; i < TextDotForHalfEdges.Count; i++)
            {
                for (int j = 0; j < TextDotForHalfEdges[i].Count; j++)
                {
                    args.Display.EnableDepthTesting(false);
                    args.Display.DrawDot(TextDotForHalfEdges[i][j], Color.DodgerBlue, Color.White, Color.White);
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
            get { return new Guid("f562e1a1-8bf4-44d6-ac8b-e2a1b79cfdca"); }
        }
    }
}