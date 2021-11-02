using Grasshopper.Kernel;
using Rhino.Geometry;
using Plankton;
using PlanktonGh;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcDualPolygonSplit : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcDualPolygonSplit class.
        /// </summary>
        public GhcDualPolygonSplit()
          : base("DualPolygonSplit", "DualPolygonSplit",
              "对对偶图中的奇异多边形进行拆分",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {

        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("DualHalfedgeMesh", "DHM", "生成的对偶图（半边数据结构）", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("RegeneratedDualHalfedgeMesh", "RDHM", "重新生成的对偶图（半边数据结构）", GH_ParamAccess.item);

            // 看是否需要调整原来的图结构（半边数据结构）


            pManager.AddGenericParameter("DebugVerticesOutput", "DebugV", "Debug结果顶点", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugHalfedgesOutput", "DebugH", "Debug结果半边", GH_ParamAccess.list);
            pManager.AddGenericParameter("DebugFacesOutput", "DebugF", "Debug结果面", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            DA.GetData("GlobalParameter", ref globalParameter);
            int innerNodeCount = globalParameter.VolumeNodeCount;
            int outerNodeCount = globalParameter.BoundaryNodeCount;

            PlanktonMesh D = new PlanktonMesh();

            if (DA.GetData<PlanktonMesh>("DualHalfedgeMesh", ref D))
            {
                /* 将对偶图中的奇异多边形进行拆点
                 * 1.首先找到奇异多边形
                 * 2.对关键边进行拆分
                 * 3.构造新的四边形
                 */

                List<int> faceEdgeCounts = new List<int>();
                int[] halfedgeindex;
                for (int i = 0; i < D.Faces.Count; i++)
                {
                    halfedgeindex = D.Faces.GetHalfedges(i);
                    faceEdgeCounts.Add(halfedgeindex.Length);
                }

                for (int i = 0; i < D.Faces.Count; i++)
                {
                    if (faceEdgeCounts[i] == 5)
                    {
                        Split5gon(ref D, i);

                        break;
                    }
                    // else
                    // {
                    //     AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "对偶图中存在边数大于9的多边形");
                    // }
                }

                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(D);
                DA.SetDataList("DebugVerticesOutput", printVertices);

                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(D);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);

                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFaces(D);
                DA.SetDataList("DebugFacesOutput", printFaces);
            }
        }

        public void Split5gon(ref PlanktonMesh planktonMesh, int splitEdgeIndex)
        {
            planktonMesh.Halfedges.SplitEdge(splitEdgeIndex);

            // return planktonMesh;
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
            get { return new Guid("8e5f5b1c-1e0a-4980-9b8e-d3628ff6ae32"); }
        }
    }
}