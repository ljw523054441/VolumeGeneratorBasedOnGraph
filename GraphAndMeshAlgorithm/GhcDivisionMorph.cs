using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Linq;
using System.Drawing;
using System.Collections.Generic;
using VolumeGeneratorBasedOnGraph.Class;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDivisionMorph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDivisionMorph class.
        /// </summary>
        public GhcDivisionMorph()
          : base("DivisionMorph", "DivisionMorph",
              "Description",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            PolylinePoints = new List<List<Point3d>>();
            //InnerNodeTextDots = new List<TextDot>();
        }

        private int Thickness;
        private List<List<Point3d>> PolylinePoints;
        //private List<TextDot> InnerNodeTextDots;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualPlanktonMesh", "DualPlanktonMesh", "", GH_ParamAccess.item);
            pManager.AddBoxParameter("Reference", "R", "Reference box to map from", GH_ParamAccess.item);
            pManager.AddSurfaceParameter("Surface", "S", "Surface to map onto", GH_ParamAccess.item);
            pManager.AddIntervalParameter("U Domain", "U", "Surface space U extents", GH_ParamAccess.item);
            pManager.AddIntervalParameter("V Domain", "V", "Surface space V extents", GH_ParamAccess.item);
            pManager.AddIntervalParameter("W Domain", "W", "Surface space W extents", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("DualPlanktonMesh", "DualPlanktonMesh", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Thickness = 2;
            // DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();

            PlanktonMesh dualPlanktonMesh = new PlanktonMesh();
            
            Box unset = Box.Unset;
            Brep brep = null;
            Interval unset2 = Interval.Unset;
            Interval unset3 = Interval.Unset;
            Interval unset4 = Interval.Unset;

            if (!DA.GetData("DualPlanktonMesh",ref dualPlanktonMesh))
            {
                return;
            }
            if (!DA.GetData<Box>("Reference", ref unset))
            {
                return;
            }
            if (!DA.GetData<Brep>("Surface", ref brep))
            {
                return;
            }
            if (!DA.GetData<Interval>("U Domain", ref unset2))
            {
                return;
            }
            if (!DA.GetData<Interval>("V Domain", ref unset3))
            {
                return;
            }
            if (!DA.GetData<Interval>("W Domain", ref unset4))
            {
                return;
            }

            if (!unset.IsValid)
            {
                return;
            }
            if (!unset2.IsValid)
            {
                return;
            }
            if (!unset3.IsValid)
            {
                return;
            }
            if (!unset4.IsValid)
            {
                return;
            }
            if (unset2.IsSingleton)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Surface space U extents cannot be zero");
                return;
            }
            if (unset3.IsSingleton)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Surface space U extents cannot be zero");
                return;
            }
            if (unset4.IsSingleton)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Surface space U extents cannot be zero");
                return;
            }
            if (unset.Volume == 0.0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Reference space volume cannot be zero");
                return;
            }

            SurfaceUVWSpaceMorph surfaceUVWSpaceMorph = new SurfaceUVWSpaceMorph(brep.Faces[0], unset, unset2, unset3, unset4);
            surfaceUVWSpaceMorph.PreserveStructure = false;
            surfaceUVWSpaceMorph.QuickPreview = false;

            List<GH_Point> vertices = new List<GH_Point>();
            for (int i = 0; i < dualPlanktonMesh.Vertices.Count; i++)
            {
                Point3d point = dualPlanktonMesh.Vertices[i].ToPoint3d();
                vertices.Add(new GH_Point(point));
            }

            List<Point3d> newVertices = new List<Point3d>();
            for (int i = 0; i < vertices.Count; i++)
            {
                GH_Point ghPoint = vertices[i];
                GH_Point ghPointDP = ghPoint.DuplicatePoint();
                ghPointDP = (GH_Point)ghPointDP.Morph(surfaceUVWSpaceMorph);

                newVertices.Add(ghPointDP.Value);
            }

            List<PlanktonXYZ> newXYZ = new List<PlanktonXYZ>();
            for (int i = 0; i < newVertices.Count; i++)
            {
                newXYZ.Add(new PlanktonXYZ((float)newVertices[i].X, (float)newVertices[i].Y, (float)newVertices[i].Z));
            }

            List<List<int>> faces = new List<List<int>>();
            for (int i = 0; i < dualPlanktonMesh.Faces.Count; i++)
            {
                faces.Add(new List<int>());
                faces[i].AddRange(dualPlanktonMesh.Faces.GetFaceVertices(i));
            }

            PlanktonMesh newDualPlanktonMesh = new PlanktonMesh();
            newDualPlanktonMesh.Vertices.AddVertices(newXYZ);
            newDualPlanktonMesh.Faces.AddFaces(faces);

            DA.SetData("DualPlanktonMesh", newDualPlanktonMesh);

            #region 可视化
            PolylinePoints.Clear();
            PolylinePoints = new List<List<Point3d>>();
            for (int i = 0; i < faces.Count; i++)
            {
                PolylinePoints.Add(new List<Point3d>());
                int[] vi = newDualPlanktonMesh.Faces.GetFaceVertices(i);
                for (int j = 0; j < vi.Length; j++)
                {
                    Point3d point = newDualPlanktonMesh.Vertices[vi[j]].ToPoint3d();
                    PolylinePoints[i].Add(point);
                }
                PolylinePoints[i].Add(PolylinePoints[i].First());
            }
            #endregion
        }


        public class SurfaceUVWSpaceMorph : SpaceMorph
        {
            protected Line m_x;
            protected Line m_y;
            protected Line m_z;
            protected Surface m_srf;
            protected Interval m_u;
            protected Interval m_v;
            protected Interval m_w;

            public SurfaceUVWSpaceMorph(Surface surface, Box reference, Interval u, Interval v, Interval w)
            {
                this.m_x = new Line(reference.PointAt(0.0, 0.0, 0.0), reference.PointAt(1.0, 0.0, 0.0));
                this.m_y = new Line(reference.PointAt(0.0, 0.0, 0.0), reference.PointAt(0.0, 1.0, 0.0));
                this.m_z = new Line(reference.PointAt(0.0, 0.0, 0.0), reference.PointAt(0.0, 0.0, 1.0));
                this.m_srf = surface;
                this.m_u = u;
                this.m_v = v;
                this.m_w = w;
            }

            public override Point3d MorphPoint(Point3d point)
            {
                double num = this.m_x.ClosestParameter(point);
                double num2 = this.m_y.ClosestParameter(point);
                double num3 = this.m_z.ClosestParameter(point);
                num = this.m_u.ParameterAt(num);
                num2 = this.m_v.ParameterAt(num2);
                num3 = this.m_w.ParameterAt(num3);
                Point3d point3d = this.m_srf.PointAt(num, num2);
                Vector3d vector3d = this.m_srf.NormalAt(num, num2);
                return point3d + vector3d * num3;
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            // base.DrawViewportWires(args);
            for (int i = 0; i < PolylinePoints.Count; i++)
            {
                args.Display.DrawPolyline(PolylinePoints[i], Color.DarkOrange, Thickness);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            // base.DrawViewportMeshes(args);
            for (int i = 0; i < PolylinePoints.Count; i++)
            {
                args.Display.DrawPolyline(PolylinePoints[i], Color.DarkOrange, Thickness);
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
            get { return new Guid("6ff4c6c6-ca59-4a27-a3a2-ddbd0e00452d"); }
        }
    }
}