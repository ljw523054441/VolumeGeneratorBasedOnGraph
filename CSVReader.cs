using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph
{
    public class CSVReader : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CSVReader class.
        /// </summary>
        public CSVReader()
          : base("CSVReader", "CSVReader",
              "读取CSV文件",
              "VolumeGeneratorBasedOnGraph", "CSV Import")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("CSVFilePath", "Path", "CSV文件的路径", GH_ParamAccess.item);
            pManager.AddIntegerParameter("NumberOfCustomBoundaryNode", "Number", "自定义的边界邻接点的数量", GH_ParamAccess.item, 0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("ConnectivityArributeTree", "ConnectivityTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "BoundaryAdjacencyTree", "体量与边界的邻接关系树", GH_ParamAccess.tree);
            pManager.AddGenericParameter("LabelList", "LabelList", "体量标签列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("AttributeList", "AttributeList", "体量属性列表", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Receive input 接收输入
            string csvPath = null;
            DA.GetData("CSVFilePath", ref csvPath);

            int number = 0;
            DA.GetData("NumberOfCustomBoundaryNode", ref number);

            // Initialize variables 初始化变量
            List<string> nodeLabelList = new List<string>();
            //List<double> volumeAreaAttributeList = new List<double>();
            List<NodeAttribute> nodeAttributesList = new List<NodeAttribute>();
            DataTree<int> volumeConnectivityArributeDataTree = new DataTree<int>();
            DataTree<int> volumeBoundaryAdjacencyDataTree = new DataTree<int>();


            // Read the data of the CSV files as individual lines 按行读取csv文件中的数据
            string[] csvLines = System.IO.File.ReadAllLines(csvPath);
            // Parse all lines 解析所有行的数据
            for (int i = 1; i < csvLines.Length; i++)
            {
                // Split rowData with "," 按逗号分隔每一行的字符串
                string[] rowData = csvLines[i].Split(',');

                if (rowData[0] == "")
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Node must have a Label, please check csv file.");
                    return;
                }
                if (rowData[1] == "")
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Node must have attribute, please check csv file.");
                    return;
                }
                // Add first and second colume data into corresponding list 将得到的第一列和第二列数据分别添加到对应的列表中
                nodeLabelList.Add(rowData[0]);
                nodeAttributesList.Add(new NodeAttribute(rowData[0], Convert.ToDouble(rowData[1])));
                //volumeAreaAttributeList.Add(Convert.ToDouble(rowData[1]));

                // Split ConnectivityArribute with "-" and add it to corresponding datatree with correct path 用-号将第三列数据分隔，并按照对应的路径添加到树形结构中
                string[] connectivity = rowData[2].Split('-');
                GH_Path path = new GH_Path(0, i - 1);
                for (int j = 0; j < connectivity.Length; j++)
                {
                    if (connectivity[j] == "")
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Node must have connectivity, please check csv file.");
                        return;
                    }
                    if (Convert.ToInt32(connectivity[j]) >= csvLines.Length - 1)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, connectivity[j].ToString()+"节点的连接对象的序号不在节点列表中");
                    }
                    else
                    {
                        volumeConnectivityArributeDataTree.Add(Convert.ToInt32(connectivity[j]), path);
                    }
                }

                // 边界邻接点的数量，先设置成4，可以继续增加，比如8
                // NodeAttribute.NEWSCount
                NodeAttribute.BoundaryNodeCount = 4 + number;
                NodeAttribute.VolumeNodeCount = csvLines.Length - 1;

                string[] adjacency = rowData[3].Split('-');
                for (int j = 0; j < adjacency.Length; j++)
                {
                    if (adjacency[j] == "")
                    {
                        volumeBoundaryAdjacencyDataTree.Add(-1 - NodeAttribute.BoundaryNodeCount, path);
                        continue;
                    }
                    else if (Convert.ToInt32(adjacency[j]) >= NodeAttribute.BoundaryNodeCount)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查adjacency在csv文件中的输入。");
                    }
                    //else if (adjacency[j] == "N")
                    //{
                    //    volumeBoundaryAdjacencyDataTree.Add(0 - NodeAttribute.NEWSCount, path);
                    //    continue;
                    //}
                    //else if (adjacency[j] == "E")
                    //{
                    //    volumeBoundaryAdjacencyDataTree.Add(1 - NodeAttribute.NEWSCount, path);
                    //    continue;
                    //}
                    //else if (adjacency[j] == "W")
                    //{
                    //    volumeBoundaryAdjacencyDataTree.Add(2 - NodeAttribute.NEWSCount, path);
                    //    continue;
                    //}
                    //else if (adjacency[j] == "S")
                    //{
                    //    volumeBoundaryAdjacencyDataTree.Add(3 - NodeAttribute.NEWSCount, path);
                    //    continue;
                    //}
                    else
                    {
                        volumeBoundaryAdjacencyDataTree.Add(Convert.ToInt32(adjacency[j]) - NodeAttribute.BoundaryNodeCount, path);
                    }
                    // volumeBoundaryAdjacencyDataTree.Add(Convert.ToInt32(adjacency[j]), path);
                }
            }



            // Output 输出
            DA.SetDataTree(0, volumeConnectivityArributeDataTree);
            DA.SetDataTree(1, volumeBoundaryAdjacencyDataTree);
            DA.SetDataList("LabelList", nodeLabelList);
            DA.SetDataList("AttributeList", nodeAttributesList);
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
            get { return new Guid("a08c725c-f6d4-4fda-a553-a11c8ec3d3d8"); }
        }
    }
}