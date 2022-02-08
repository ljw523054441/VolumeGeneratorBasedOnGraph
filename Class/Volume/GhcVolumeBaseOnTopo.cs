using Grasshopper.Kernel;
using Plankton;
using PlanktonGh;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class.Volume
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
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("DualGraphWithHM", "DualGraphWithHM", "", GH_ParamAccess.item);
            pManager.AddGenericParameter("newDualPlanktonMesh", "newDualPlanktonMesh", "", GH_ParamAccess.item);

            pManager.AddNumberParameter("SetbackList", "SetbackList", "", GH_ParamAccess.list);


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
            DualGraphWithHM dualGraphWithHM = new DualGraphWithHM();
            PlanktonMesh dual = new PlanktonMesh();

            PlanktonMesh graph = new PlanktonMesh();
            List<GraphNode> graphNodes = new List<GraphNode>();
            List<List<int>> graphTable = new List<List<int>>();

            List<double> setbackList = new List<double>();

            if (DA.GetData<DualGraphWithHM>("DualGraphWithHM", ref dualGraphWithHM)
                && DA.GetData<PlanktonMesh>("newDualPlanktonMesh", ref dual)
                && DA.GetDataList("SetbackList", setbackList))
            {
                dualGraphWithHM.DualPlanktonMesh = new PlanktonMesh(dual);
                
                graph = dualGraphWithHM.PlanktonMesh;
                graphNodes = dualGraphWithHM.GraphNodes;
                graphTable = dualGraphWithHM.GraphTables;

                // 各边退界
                double setbackN = setbackList[0];
                double setbackW = setbackList[1];
                double setbackS = setbackList[2];
                double setbackE = setbackList[3];
                double setback = setbackList[4];

                #region 得到所有的innerNode的序号和outerNode的序号
                List<int> innerNodeIndexs = dualGraphWithHM.InnerNodeIndexList;
                List<int> outerNodeIndexs = dualGraphWithHM.OuterNodeIndexList;
                int innerNodeCount = dualGraphWithHM.InnerNodeCount;
                int outerNodeCount = dualGraphWithHM.OuterNodeCount;
                #endregion

                //#region 计算faceIndexsAroundOuterNodes，找到每个outerNode相关的面和发出的半边
                //// 对于每个outerNode找到由它发出的半边
                //List<List<int>> halfedgeIndexsStartFromOuterNodes = new List<List<int>>();
                //// 对于每个outerNode找到与它相关的面
                //List<List<int>> faceIndexsAroundOuterNodes = new List<List<int>>();

                //for (int i = 0; i < outerNodeIndexs.Count; i++)
                //{
                //    // 找到每个outerNode所发出的halfedge的index
                //    int[] halfedgeIndexsStartFromOuterNode = graph.Vertices.GetHalfedges(outerNodeIndexs[i]);
                //    halfedgeIndexsStartFromOuterNodes.Add(new List<int>());
                //    halfedgeIndexsStartFromOuterNodes[i].AddRange(halfedgeIndexsStartFromOuterNode);
                //    // 找到每个outerNode所邻接的Face的index，-1表示邻接外界
                //    int[] faceIndexsAroundOuterNode = graph.Vertices.GetVertexFaces(outerNodeIndexs[i]);
                //    faceIndexsAroundOuterNodes.Add(new List<int>());
                //    faceIndexsAroundOuterNodes[i].AddRange(faceIndexsAroundOuterNode);
                //}
                //#endregion


                for (int i = 0; i < dual.Faces.Count; i++)
                {

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
            get { return new Guid("4b850a3f-93e2-44a8-b44a-79f109cd61c3"); }
        }
    }
}