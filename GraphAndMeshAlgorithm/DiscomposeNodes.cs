using Grasshopper.Kernel;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Plankton;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class DiscomposeNodes : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DiscomposeNodes class.
        /// </summary>
        public DiscomposeNodes()
          : base("DiscomposeNodes", "DiscomposeNodes",
              "将度大于4的顶点按照规则进行拆解",
              "VolumeGeneratorBasedOnGraph", "GraphEmbeding")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GlobalParameter", "全局参数传递", GH_ParamAccess.item);
            pManager.AddGenericParameter("TheChosenTriangleHalfedgeMesh", "THMesh", "所选择的那个三角形剖分结果(半边数据结构)", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("NewTriangleHalfedgeMesh", "NTHMesh", "经过顶点拆解后的三角形剖分结果(半边数据结构)", GH_ParamAccess.item);

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

            PlanktonMesh P = new PlanktonMesh();

            if (DA.GetData<PlanktonMesh>("TheChosenTriangleHalfedgeMesh", ref P))
            {
                /* 获取度大于4的顶点的序号 */
                List<int> degreeOver5Indexs = new List<int>();
                for (int i = 0; i < P.Vertices.Count; i++)
                {
                    int degree = P.Vertices.GetValence(i);
                    if (degree > 4)
                    {
                        degreeOver5Indexs.Add(i);
                    }
                }

                //for (int i = 0; i < degreeOver5Indexs.Count; i++)
                //{
                //    switch (degreeOver5Indexs[i])
                //    {
                //        case 5:

                //            break;

                //        case 6:
                //            break;

                //        case 7:
                //            break;

                //        case 8:
                //            break;

                //        default:
                //            break;
                //    }
                //}

                // 要把下面的switch结构套进这个for循环中
                //for (int i = 0; i < degreeOver5Indexs.Count; i++)
                //{

                //}

                switch (P.Vertices.GetValence(degreeOver5Indexs[1]))
                    {
                        case 5:
                            
                            P.Vertices.SplitVertex(33, 15);
                            break;

                        case 6:
                            break;

                        case 7:
                            break;

                        case 8:
                            break;

                        default:
                            break;
                    }

                DA.SetData("NewTriangleHalfedgeMesh", P);

                List<string> printVertices = new List<string>();
                printVertices = UtilityFunctions.PrintVertices(P);
                DA.SetDataList("DebugVerticesOutput", printVertices);

                List<string> printHalfedges = new List<string>();
                printHalfedges = UtilityFunctions.PrintHalfedges(P);
                DA.SetDataList("DebugHalfedgesOutput", printHalfedges);


                List<string> printFaces = new List<string>();
                printFaces = UtilityFunctions.PrintFaces(P);
                DA.SetDataList("DebugFacesOutput", printFaces);

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
            get { return new Guid("167d61e3-6ad9-4512-bda1-f3518b7d9f41"); }
        }
    }
}