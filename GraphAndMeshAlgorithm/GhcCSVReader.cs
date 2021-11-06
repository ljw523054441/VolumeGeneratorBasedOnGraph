﻿using Grasshopper;
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
              "VolumeGeneratorBasedOnGraph", "Construct Graph")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("GlobalParameter", "GP", "全局参数传递", GH_ParamAccess.item);

            pManager.AddGenericParameter("CSVFilePath", "Path", "CSV文件的路径", GH_ParamAccess.item);
            // pManager.AddIntegerParameter("NumberOfCustomBoundaryNode", "Number", "自定义的边界邻接点的数量", GH_ParamAccess.item, 0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("VolumeConnectivityTree", "CTree", "体量连接关系树", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("BoundaryAdjacencyTree", "ATree", "体量与边界的邻接关系树", GH_ParamAccess.tree);
            // pManager.AddGenericParameter("LabelList", "LabelList", "体量标签列表", GH_ParamAccess.list);
            pManager.AddGenericParameter("VolumeNodeAttributes", "VNAttributes", "体量属性列表", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // 全局参数传递
            GlobalParameter globalParameter = new GlobalParameter();
            

            if (DA.GetData("GlobalParameter", ref globalParameter))
            {
                int volumeNodeCount = globalParameter.VolumeNodeCount;
                int boundaryNodeCount = globalParameter.BoundaryNodeCount;

                // Receive input 接收输入
                string csvPath = null;
                DA.GetData("CSVFilePath", ref csvPath);

                #region Initialize variables 初始化变量
                List<string> nodeLabelList = new List<string>();
                List<NodeAttribute> nodeAttributesList = new List<NodeAttribute>();
                DataTree<int> volumeConnectivityArributeDataTree = new DataTree<int>();
                DataTree<int> volumeBoundaryAdjacencyDataTree = new DataTree<int>();
                #endregion

                #region Read the data of the CSV files as individual lines 按行读取csv文件中的数据
                string[] csvLines = System.IO.File.ReadAllLines(csvPath);
                #endregion

                #region Parse all lines 解析所有行的数据
                for (int i = 1; i < csvLines.Length; i++)
                {
                    #region Split rowData with "," 按逗号分隔每一行的字符串
                    string[] rowData = csvLines[i].Split(',');
                    #endregion

                    #region 读取时可能出现的报错
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
                    #endregion

                    #region Add first and second colume data into corresponding list 将得到的第一列和第二列数据分别添加到对应的列表中
                    nodeLabelList.Add(rowData[0]);
                    nodeAttributesList.Add(new NodeAttribute(rowData[0],
                                                             Convert.ToDouble(rowData[1]),
                                                             Convert.ToInt32(rowData[4]),
                                                             Convert.ToInt32(rowData[5]),
                                                             Convert.ToInt32(rowData[6]),
                                                             Convert.ToDouble(rowData[7]),
                                                             Convert.ToDouble(rowData[8]),
                                                             Convert.ToDouble(rowData[9])
                                                             ));
                    //volumeAreaAttributeList.Add(Convert.ToDouble(rowData[1]));
                    #endregion

                    #region Split ConnectivityArribute with "-" and add it to corresponding datatree with correct path 用-号将第三列数据分隔，并按照对应的路径添加到树形结构中
                    string[] connectivity = rowData[2].Split('-');
                    GH_Path path = new GH_Path(0, i - 1);
                    for (int j = 0; j < connectivity.Length; j++)
                    {
                        volumeConnectivityArributeDataTree.EnsurePath(path);

                        if (connectivity[j] == "")
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "节点必须有连接关系, 请检查CSV文件。");
                            return;
                        }
                        if (Convert.ToInt32(connectivity[j]) >= csvLines.Length - 1)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, connectivity[j].ToString() + "节点的连接对象的序号不在节点列表中。");
                        }
                        else
                        {
                            volumeConnectivityArributeDataTree.Add(Convert.ToInt32(connectivity[j]), path);
                        }
                    }
                    nodeAttributesList[i - 1].ConnectivityTable = volumeConnectivityArributeDataTree.Branch(path).ToArray();
                    #endregion

                    #region Split AjacencyArribute with "-" and add it to corresponding datatree with correct path 用-号将第四列数据分隔，并按照对应的路径添加到树形结构中
                    string[] adjacency = rowData[3].Split('-');
                    for (int j = 0; j < adjacency.Length; j++)
                    {
                        volumeBoundaryAdjacencyDataTree.EnsurePath(path);

                        // 如果是空，就跳过它继续后面的循环
                        if (adjacency[j] == "")
                        {
                            volumeBoundaryAdjacencyDataTree.AddRange(new List<int>(), path);
                            continue;
                        }
                        else if (Convert.ToInt32(adjacency[j]) >= boundaryNodeCount)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查adjacency中序号的输入，是否超过了BoundaryNode数量的上限。");
                        }
                        else if (Convert.ToInt32(adjacency[j]) < 0)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "请检查adjacency中序号的输入，是否超过了BoundaryNode数量的下限。");
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
                            volumeBoundaryAdjacencyDataTree.Add(Convert.ToInt32(adjacency[j]), path);
                        }
                        // volumeBoundaryAdjacencyDataTree.Add(Convert.ToInt32(adjacency[j]), path);
                    }
                    nodeAttributesList[i - 1].AdjacencyTable = volumeBoundaryAdjacencyDataTree.Branch(path).ToArray();
                    #endregion
                }
                #endregion



                #region Output 输出
                DA.SetDataTree(0, volumeConnectivityArributeDataTree);
                DA.SetDataTree(1, volumeBoundaryAdjacencyDataTree);
                DA.SetDataList("VolumeNodeAttributes", nodeAttributesList);
                #endregion
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
            get { return new Guid("a08c725c-f6d4-4fda-a553-a11c8ec3d3d8"); }
        }
    }
}