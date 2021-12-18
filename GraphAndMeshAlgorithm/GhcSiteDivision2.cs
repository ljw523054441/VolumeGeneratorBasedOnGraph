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

namespace VolumeGeneratorBasedOnGraph.GraphAndMeshAlgorithm
{
    public class GhcSiteDivision2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcSiteDivision2 class.
        /// </summary>
        public GhcSiteDivision2()
          : base("GhcSiteDivision2", "SiteDivision2",
              "进行场地划分",
              "VolumeGeneratorBasedOnGraph", "CreateVolume")
        {
            VolumeJunctionTextDot = new List<TextDot>();
        }

        private List<TextDot> VolumeJunctionTextDot;

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SortedBoundarySegments", "SBS", "经过排序后的BoundarySegment", GH_ParamAccess.list);
            pManager.AddIntegerParameter("BSIndexContainVolumeJunctions", "BSI", "包含边界分裂点的BoundarySegment序号", GH_ParamAccess.list);

            pManager.AddTextParameter("VolumeJunctionsTexts", "VTD", "表示边界分裂点的Text", GH_ParamAccess.list);
            // pManager.AddGenericParameter("PairVolumeJunctionsIndex", "PJI", "作为分界线的一对分裂点的序号", GH_ParamAccess.list);
            pManager.AddBooleanParameter("NeedToConnect", "NTC", "的这一对分裂点是否作为分界线", GH_ParamAccess.list);

            pManager.AddNumberParameter("tList", "ts", "每个边界分裂点对应的t值", GH_ParamAccess.list);

            pManager.AddNumberParameter("AdditionalTXList", "atXs", "绘制内部分界线需要的tX值", GH_ParamAccess.list);
            pManager.AddNumberParameter("AdditionalTYList", "atYs", "绘制内部分界线需要的tY值", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("DeBug Polyline", "DP", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<BoundarySegment> sortedBoundarySegments = new List<BoundarySegment>(); 
            List<int> bSIndexContainVolumeJunctions = new List<int>();

            List<string> volumeJunctionsTexts = new List<string>();
            //List<int[]> pairVolumeJunctionsIndex = new List<int[]>();
            List<bool> needToConnect = new List<bool>();

            List<double> tList = new List<double>();
            List<double> additionalTXList = new List<double>();
            List<double> additionalTYList = new List<double>();

            if (DA.GetDataList("SortedBoundarySegments", sortedBoundarySegments) 
                && DA.GetDataList("BSIndexContainVolumeJunctions", bSIndexContainVolumeJunctions)
                && DA.GetDataList("VolumeJunctionsTexts", volumeJunctionsTexts)
                // && DA.GetDataList("PairVolumeJunctionsIndex", pairVolumeJunctionsIndex)
                && DA.GetDataList("NeedToConnect", needToConnect)
                && DA.GetDataList("tList", tList)
                && DA.GetDataList("AdditionalTXList",additionalTXList)
                && DA.GetDataList("AdditionalTYList",additionalTYList))
            {
                //// 对输入的Curve类型的reorganizedBoundary进行类型转换，转换成Curve类的子类Polyline
                //Polyline reorganizedBoundary = null;
                //reorganizedBoundaryCurve.TryGetPolyline(out reorganizedBoundary);

                //// gh_Structure转化为dataTree
                //DataTree<Point3d> volumeJunctionPointsDT = new DataTree<Point3d>();
                //foreach (GH_Path gh_Path in gh_Structure.Paths)
                //{
                //    List<Point3d> branchItems = new List<Point3d>();
                //    foreach (GH_Point element in gh_Structure.get_Branch(gh_Path))
                //    {
                //        branchItems.Add(element.Value);
                //    }
                //    volumeJunctionPointsDT.AddRange(branchItems, gh_Path);
                //}

                List<Point3d> cornerPoints = new List<Point3d>();
                for (int i = 0; i < sortedBoundarySegments.Count; i++)
                {
                    cornerPoints.Add(sortedBoundarySegments[i].From);
                }

                BoundingBox boundingBox = new BoundingBox(cornerPoints);
                

                List<BoundarySegment> bsContainVolumeJunctions = new List<BoundarySegment>();
                for (int i = 0; i < bSIndexContainVolumeJunctions.Count; i++)
                {
                    bsContainVolumeJunctions.Add(sortedBoundarySegments[bSIndexContainVolumeJunctions[i]]);
                }

                List<double> tForBSLine = new List<double>();
                int iter = 0;
                while (tForBSLine.Count < bsContainVolumeJunctions.Count)
                {
                    if (iter > tList.Count - 1)
                    {
                        throw new Exception("tList的数量不足");
                    }

                    tForBSLine.Add(tList[iter]);
                    iter++;
                }

                // 全局的addtionalTX的使用量
                int additionalTXCount = 0;
                // 全局的addtionalTY的使用量
                int additionalTYCount = 0;

                List<Point3d> volumeJunctions = new List<Point3d>();
                for (int i = 0; i < bsContainVolumeJunctions.Count; i++)
                {
                    volumeJunctions.Add(bsContainVolumeJunctions[i].Line.PointAt(tForBSLine[i]));
                }


                #region 构造折线，先做正X正Y
                // 相邻的两个volumeJunction所在的BS的序号
                List<int[]> bSIndexPair = new List<int[]>();
                for (int i = 0; i < bSIndexContainVolumeJunctions.Count; i++)
                {
                    int[] pair = new int[2] { bSIndexContainVolumeJunctions[i], bSIndexContainVolumeJunctions[(i + 1) % bSIndexContainVolumeJunctions.Count] };
                    bSIndexPair.Add(pair);
                }

                // 对于每个bsIndexPair，进行连接线的绘制，后面的会在绘制时，参照前面已经绘制过的轨迹
                bool isFirstDrawn = true;
                for (int i = 0; i < bSIndexPair.Count; i++)
                {
                    // 如果不需要作为分界线，直接当前的计算
                    if (!needToConnect[i])
                    {
                        continue;
                    }
                    
                    int interval = CalInterval(bSIndexPair[i], sortedBoundarySegments.Count);
                    // 如果相邻的两个volumeJunction属于同一个volume，那么不用做折线
                    if (interval == 0)
                    {
                        continue;
                    }

                    // 如果间隔是1，那么直接求交点，构成折线
                    if (interval == 1)
                    {
                        
                        Point3d point1 = new Point3d(volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][0])].X, volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][1])].Y, 0);
                        Point3d point2 = new Point3d(volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][1])].X, volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][0])].Y, 0);

                        Point3d corner = bsContainVolumeJunctions[bSIndexPair[i][0]].To;

                        if (isFirstDrawn)
                        {
                            // 如果之前没有其他线被绘制过
                            List<Point3d> pointList = new List<Point3d>();
                            pointList.Add(volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][0])]);
                            pointList.Add(corner);
                            pointList.Add(volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][1])]);
                            pointList.Add(point1);

                            if (IsConvex(pointList))
                            {
                                // 选point1
                                Point3d[] line1 = new Point3d[2] { volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][0])], point1 };
                                Point3d[] line2 = new Point3d[2] { point1, volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][1])] };
                            }
                            else
                            {
                                // 选point2
                                Point3d[] line1 = new Point3d[2] { volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][0])], point2 };
                                Point3d[] line2 = new Point3d[2] { point2, volumeJunctions[bSIndexContainVolumeJunctions.IndexOf(bSIndexPair[i][1])] };
                            }
                        }
                        else
                        {
                            // 如果之前有其他线被绘制过
                        }
                        
                    }

                }

                #endregion

                #region 可视化
                VolumeJunctionTextDot.Clear();

                for (int i = 0; i < volumeJunctions.Count; i++)
                {
                    TextDot textDot = new TextDot(volumeJunctionsTexts[i], volumeJunctions[i]);
                    VolumeJunctionTextDot.Add(textDot);
                }
                #endregion


                volumeJunctions.Add(volumeJunctions[0]);
                Polyline debugPolyline = new Polyline(volumeJunctions);

                DA.SetData("DeBug Polyline", debugPolyline);
            }
        }

        /// <summary>
        /// 计算两个相邻的volumeJunction点之间间隔了几个转角
        /// </summary>
        /// <param name="bsIndexPair"></param>
        /// <param name="bsCount"></param>
        /// <returns></returns>
        private int CalInterval(int[] bsIndexPair,int bsCount)
        {
            int counterClockwiseCount;
            int clockwiseCount;
            int index0 = bsIndexPair[0];
            int index1 = bsIndexPair[1];

            if (index0 < index1)
            {
                counterClockwiseCount = index1 - index0;
                clockwiseCount = bsCount - index1 + index0;
            }
            else
            {
                counterClockwiseCount = bsCount - index0 + index1;
                clockwiseCount = index0 - index1;
            }

            if (counterClockwiseCount < clockwiseCount)
            {
                return counterClockwiseCount;
            }
            else
            {
                return clockwiseCount;
            }


        }

        private List<Point3d> MakeStepLine(Point3d point1,
                                           Point3d point2,
                                           int interval,
                                           BoundarySegment bsForPoint1,
                                           BoundarySegment bsForPoint2,
                                           List<Point3d[]> linesPassingThroughPoint1,
                                           BoundingBox boundingBox,
                                           List<double> additionalTXList, 
                                           List<double> additionalTYList, 
                                           ref int additionalTXCount,
                                           ref int additionalTYCount)
        {
            double totalDX = point2.X - point1.X;
            double totalDY = point2.Y - point1.Y;

            if (interval % 2 == 1)
            {
                // 如果转折数是奇数
                #region 基础计算
                List<Point3d> turningPoints = new List<Point3d>();

                int quadrant = WhichQuadrant(point1, point2);

                #region 求point1和point2所对应的BS的法线（逆时针90度）
                // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
                Vector3d vectorForBS1 = new Vector3d(bsForPoint1.To - bsForPoint1.From);
                Vector3d vectorForBS2 = new Vector3d(bsForPoint2.To - bsForPoint2.From);
                Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.From, bsForPoint1.To);
                Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.From, bsForPoint2.To);

                // 求在x轴正方向或y轴正方向的投影
                Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
                Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
                // 投影后向量的normal，为正X正Y方向
                Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
                Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
                #endregion

                #region 判断两个法向量的交点，是否都在两个向量的正方向上
                bool isIntersectOnPositiveDirection = true;
                // 即由交点到向量原点所构成的向量是否与原来的法向量方向相同
                // 1.求交点
                Point3d intersectPoint = GetLineIntersection(point1, nVectorForBS1, point2, nVectorForBS2);
                if (intersectPoint == Point3d.Unset)
                {
                    throw new Exception("interval为奇数时划线，point1和point2的法向量平行");
                }
                // 2.构造新向量
                Vector3d newVector1 = new Vector3d(point1 - intersectPoint);
                Vector3d newVector2 = new Vector3d(point2 - intersectPoint);
                // 3.判断新向量与原向量是否同向
                if (newVector1 * nVectorForBS1 > 0 && newVector2 * nVectorForBS2 > 0)
                {
                    isIntersectOnPositiveDirection = true;
                }
                else
                {
                    isIntersectOnPositiveDirection = false;
                }
                #endregion

                #region 判断由nVectorForBS1到nVectorForBS2是顺时针还是逆时针
                bool isCounterClockwise = true;
                Vector3d crossProduct = Vector3d.CrossProduct(nVectorForBS1, nVectorForBS2);
                if (crossProduct.Z < 0)
                {
                    // 逆时针
                    isCounterClockwise = true;
                }
                else
                {
                    // 顺时针
                    isCounterClockwise = false;
                }
                #endregion
                #endregion

                if (linesPassingThroughPoint1 != null)
                {
                    // 如果point1已经被绘制过
                    #region 找到与point1相关的endPoints列表
                    List<Point3d> endPoints = new List<Point3d>();
                    for (int i = 0; i < linesPassingThroughPoint1.Count; i++)
                    {
                        endPoints.Add(linesPassingThroughPoint1[i][1]);
                    }
                    #endregion

                    if (isIntersectOnPositiveDirection)
                    {
                        // 如果两个法向量的交点，在两个向量的正方向上

                        #region 对所有的endPoint进行排序
                        List<Point3d> sortedEndPoints = new List<Point3d>();
                        List<Point3d> pointsAbove = new List<Point3d>();
                        List<Point3d> pointsBelow = new List<Point3d>();
                        switch (quadrant)
                        {
                            case 1:
                                // 第一象限，比较Y值，同时优先找below Point2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].Y <= point2.Y)
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 3:
                                // 第三象限，比较Y值，同时优先找above Point2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].Y >= point2.Y)
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 2:
                                // 第二象限，比较X值，同时优先找above Point2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].X >= point2.X)
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 4:
                                // 第四象限，比较X值，同时优先找below Point2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].X <= point2.X)
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                }
                                break;
                            default:
                                break;
                        }

                        List<Point3d> sortedPointsAbove = new List<Point3d>();
                        List<Point3d> sortedPointsBelow = new List<Point3d>();
                        #endregion

                        switch (quadrant)
                        {
                            case 1:
                                int iter = 0;
                                #region 确定第一个turningPoint的位置
                                // 1象限时，根据Y值进行排序
                                sortedPointsBelow = pointsBelow;
                                BubbleSortPoint3d(sortedPointsBelow, true);

                                if (sortedPointsBelow.Count != 0)
                                {
                                    // 从below Point2.Y的点中，找到有最大Y值的那个点，添加到turningPoint中
                                    turningPoints.Add(sortedPointsBelow.Last());
                                }
                                else
                                {
                                    // 如果没有below Point2.Y的点时，直接构造新的点（即Y值等于point2.Y的点），添加到turningPoint中
                                    turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                }
                                #endregion
                                iter++;

                                bool flag;
                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 因为已经添加了第一个turningPoint
                                    // 所以要紧接着先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 因为已经添加了第一个turningPoint
                                    // 所以要紧接着先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion
                                break;
                            case 2:
                                iter = 0;
                                #region 确定第一个turningPoint的位置
                                // 2象限时，根据X值进行排序
                                sortedPointsAbove = pointsAbove;
                                BubbleSortPoint3d(sortedPointsAbove, false);

                                if (sortedPointsAbove.Count != 0)
                                {
                                    // 从above Point2.X的点中，找到有最小X值的那个点，添加到turningPoint中
                                    turningPoints.Add(sortedPointsAbove[0]);
                                }
                                else
                                {
                                    // 如果没有above Point2.X的点时，直接构造新的点（即X值等于point2.X的点），添加到turningPoint中
                                    turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                }
                                #endregion
                                iter++;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion
                                break;
                            case 3:
                                iter = 0;
                                #region 确定第一个turningPoint的位置
                                // 3象限时，根据Y值进行排序
                                sortedPointsAbove = pointsAbove;
                                BubbleSortPoint3d(sortedPointsAbove, true);

                                if (sortedPointsAbove.Count != 0)
                                {
                                    // 从above Point2.Y的点中，找到有最小Y值的那个点，添加到turningPoint中
                                    turningPoints.Add(sortedPointsAbove[0]);
                                }
                                else
                                {
                                    // 如果没有above Point2.Y的点时，直接构造新的点（即Y值等于point2.Y的点），添加到turningPoint中
                                    turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                }
                                #endregion
                                iter++;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion
                                break;
                            case 4:
                                iter = 0;
                                #region 确定第一个turningPoint的位置
                                // 4象限时，根据X值进行排序
                                sortedPointsBelow = pointsBelow;
                                BubbleSortPoint3d(sortedPointsBelow, false);
                                //sortedPointsBelow.Reverse();

                                if (sortedPointsBelow.Count != 0)
                                {
                                    // 从below Point2.X的点中，找到有最小X值的那个点，添加到turningPoint中
                                    turningPoints.Add(sortedPointsBelow.Last());
                                }
                                else
                                {
                                    // 如果没有below Point2.X的点时，直接构造新的点（即X值等于point2.X的点），添加到turningPoint中
                                    turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                }
                                #endregion
                                iter++;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion
                                break;
                            default:
                                break;
                        }
                    }
                    else
                    {
                        // 如果两个法向量的交点，没有在两个向量的正方向是
                        if (isCounterClockwise)
                        {
                            // 由nVectorForBS1到nVectorForBS2是逆时针
                            interval = interval + 2;

                            // 添加point1
                            turningPoints.Add(point1);

                            #region 对所有的endPoint进行排序
                            List<Point3d> newEndPoints = new List<Point3d>();
                            switch (quadrant)
                            {
                                case 1:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(point1.X + (point1.X - endPoints[i].X), endPoints[i].Y, endPoints[i].Z));
                                    }
                                    break;
                                case 2:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(endPoints[i].X , point1.Y + (point1.Y - endPoints[i].Y), endPoints[i].Z));
                                    }
                                    break;
                                case 3:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(point1.X - (point1.X - endPoints[i].X), endPoints[i].Y, endPoints[i].Z));
                                    }
                                    break;
                                case 4:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(endPoints[i].X, point1.Y - (point1.Y - endPoints[i].Y), endPoints[i].Z));
                                    }
                                    break;
                                default:
                                    break;
                            }
                            
                            List<Point3d> sortedEndPoints = new List<Point3d>();
                            List<Point3d> pointsAbove = new List<Point3d>();
                            List<Point3d> pointsBelow = new List<Point3d>();
                            switch (quadrant)
                            {
                                case 1:
                                    // 第一象限，比较X值，同时优先找below Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].X <= point2.X)
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 3:
                                    // 第三象限，比较X值，同时优先找above Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].X >= point2.X)
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 2:
                                    // 第二象限，比较Y值，同时优先找below Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].Y <= point2.Y)
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 4:
                                    // 第四象限，比较Y值，同时优先找above Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].Y >= point2.Y)
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                default:
                                    break;
                            }

                            List<Point3d> sortedPointsAbove = new List<Point3d>();
                            List<Point3d> sortedPointsBelow = new List<Point3d>();
                            #endregion

                            switch (quadrant)
                            {
                                case 1:
                                    int iter = 0;
                                    #region 添加首个turningPoint
                                    // 1象限时，根据X值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, false);

                                    if (sortedPointsBelow.Count != 0)
                                    {
                                        // 从below Point2.X的点中，找到有最大X值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsBelow.Last());
                                    }
                                    else
                                    {
                                        // 如果没有below Point2.X的点时，直接构造新的点（即X值等于Point2.X的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    double yMax = boundingBox.Max.Y;
                                    double dySecond = yMax - point2.Y;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYSecond = additionalTYList[additionalTYCount - 1];
                                    double addYSecond = tYSecond * dySecond;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y + addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    bool flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 2:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 2象限时，根据Y值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, true);
                                    if (sortedPointsBelow.Count != 0)
                                    {
                                        // 从below Point2.Y的点中，找到有最大Y值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsBelow.Last());
                                    }
                                    else
                                    {
                                        // 如果没有below Point2.Y的点时，直接构造新的点（即Y值等于Point2.Y的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    double xMin = boundingBox.Min.X;
                                    double dxSecond = point2.X - xMin;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXSecond = additionalTXList[additionalTXCount - 1];
                                    double addXSecond = tXSecond * dxSecond;

                                    turningPoints.Add(new Point3d(point2.X - addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 3:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 3象限时，根据X值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, false);
                                    if (sortedPointsAbove.Count != 0)
                                    {
                                        // 从above Point2.X的点中，找到有最小X值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsAbove[0]);
                                    }
                                    else
                                    {
                                        // 如果没有above Point2.X的点时，直接构造新的点（即X值等于Point2.X的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    double yMin = boundingBox.Min.Y;
                                    dySecond = point2.Y - yMin;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYSecond = additionalTYList[additionalTYCount - 1];
                                    addYSecond = tYSecond * dySecond;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y - addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 4:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 4象限时，根据Y值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, true);
                                    if (sortedPointsAbove.Count != 0)
                                    {
                                        // 从above Point2.Y的点中，找到有最小Y值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsAbove[0]);
                                    }
                                    else
                                    {
                                        // 如果没有above Point2.Y的点时，直接构造新的点（即Y值等于Point2.Y的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    double xMax = boundingBox.Max.X;
                                    dxSecond = xMax - point2.X;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXSecond = additionalTXList[additionalTXCount - 1];
                                    addXSecond = tXSecond * dxSecond;

                                    turningPoints.Add(new Point3d(point2.X + addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                default:
                                    break;
                            }
                        }
                        else
                        {
                            // 由nVectorForBS1到nVectorForBS2是顺时针

                            // 添加point1
                            turningPoints.Add(point1);

                            #region 对所有的endPoint进行排序
                            List<Point3d> newEndPoints = new List<Point3d>();
                            switch (quadrant)
                            {
                                case 1:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(point1.X + (point1.X - endPoints[i].X), endPoints[i].Y, endPoints[i].Z));
                                    }
                                    break;
                                case 2:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(endPoints[i].X, point1.Y + (point1.Y - endPoints[i].Y), endPoints[i].Z));
                                    }
                                    break;
                                case 3:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(point1.X - (point1.X - endPoints[i].X), endPoints[i].Y, endPoints[i].Z));
                                    }
                                    break;
                                case 4:
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        newEndPoints.Add(new Point3d(endPoints[i].X, point1.Y - (point1.Y - endPoints[i].Y), endPoints[i].Z));
                                    }
                                    break;
                                default:
                                    break;
                            }

                            List<Point3d> sortedEndPoints = new List<Point3d>();
                            List<Point3d> pointsAbove = new List<Point3d>();
                            List<Point3d> pointsBelow = new List<Point3d>();
                            switch (quadrant)
                            {
                                case 1:
                                    // 第一象限，比较X值，同时优先找below Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].X <= point2.X)
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 3:
                                    // 第三象限，比较X值，同时优先找above Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].X >= point2.X)
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 2:
                                    // 第二象限，比较Y值，同时优先找below Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].Y <= point2.Y)
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                case 4:
                                    // 第四象限，比较Y值，同时优先找above Point2的
                                    for (int i = 0; i < newEndPoints.Count; i++)
                                    {
                                        if (newEndPoints[i].Y >= point2.Y)
                                        {
                                            pointsAbove.Add(newEndPoints[i]);
                                        }
                                        else
                                        {
                                            pointsBelow.Add(newEndPoints[i]);
                                        }
                                    }
                                    break;
                                default:
                                    break;
                            }

                            List<Point3d> sortedPointsAbove = new List<Point3d>();
                            List<Point3d> sortedPointsBelow = new List<Point3d>();
                            #endregion

                            switch (quadrant)
                            {
                                case 1:
                                    int iter = 0;
                                    #region 添加首个turningPoint
                                    // 1象限时，根据X值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, false);

                                    if (sortedPointsBelow.Count != 0)
                                    {
                                        // 从below Point2.X的点中，找到有最大X值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsBelow.Last());
                                    }
                                    else
                                    {
                                        // 如果没有below Point2.X的点时，直接构造新的点（即X值等于Point2.X的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    double dyMax = point2.Y - point1.Y;
                                    //double dySecond = point2.Y
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYSecond = additionalTYList[additionalTYCount - 1];
                                    double addYSecond = tYSecond * dyMax;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y - addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    bool flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 2:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 2象限时，根据Y值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, true);
                                    if (sortedPointsBelow.Count != 0)
                                    {
                                        // 从below Point2.Y的点中，找到有最大Y值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsBelow.Last());
                                    }
                                    else
                                    {
                                        // 如果没有below Point2.Y的点时，直接构造新的点（即Y值等于Point2.Y的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    double dxMax = point1.X - point2.X;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXSecond = additionalTXList[additionalTXCount - 1];
                                    double addXSecond = tXSecond * dxMax;

                                    turningPoints.Add(new Point3d(point2.X + addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 3:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 3象限时，根据X值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, false);
                                    if (sortedPointsAbove.Count != 0)
                                    {
                                        // 从above Point2.X的点中，找到有最小X值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsAbove[0]);
                                    }
                                    else
                                    {
                                        // 如果没有above Point2.X的点时，直接构造新的点（即X值等于Point2.X的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point2.X, point1.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    dyMax = point1.Y - point2.Y;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYSecond = additionalTYList[additionalTYCount - 1];
                                    addYSecond = tYSecond * dyMax;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y + addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 4:
                                    iter = 0;
                                    #region 确定第一个turningPoint的计算方式
                                    // 4象限时，根据Y值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, true);
                                    if (sortedPointsAbove.Count != 0)
                                    {
                                        // 从above Point2.Y的点中，找到有最小Y值的那个点，添加到turningPoint中
                                        turningPoints.Add(sortedPointsAbove[0]);
                                    }
                                    else
                                    {
                                        // 如果没有above Point2.Y的点时，直接构造新的点（即Y值等于Point2.Y的点），添加到turningPoint中
                                        turningPoints.Add(new Point3d(point1.X, point2.Y, point1.Z));
                                    }
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    dxMax = point2.X - point1.X;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXSecond = additionalTXList[additionalTXCount - 1];
                                    addXSecond = tXSecond * dxMax;

                                    turningPoints.Add(new Point3d(point2.X - addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                    return turningPoints;
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制

                    if (isIntersectOnPositiveDirection)
                    {
                        // 如果两个法向量的交点，在两个向量的正方向上
                        switch (quadrant)
                        {
                            case 1:
                                int iter = 0;
                                bool flag;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 2:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 3:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 4:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            default:
                                break;
                        }
                    }
                    else
                    {
                        // 如果两个法向量的交点，没有在两个向量的正方向上
                        if (isCounterClockwise)
                        {
                            // 由nVectorForBS1到nVectorForBS2是逆时针
                            interval = interval + 2;

                            // 添加point1
                            turningPoints.Add(point1);

                            switch (quadrant)
                            {
                                case 1:
                                    int iter = 0;
                                    #region 添加首个turningPoint
                                    // 找boundingBox的左界
                                    double xMin = boundingBox.Min.X;
                                    double dxFirst;
                                    double addXFirst;

                                    // point1一定更靠近左边
                                    dxFirst = point1.X - xMin;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXFirst = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tXFirst * dxFirst;

                                    additionalTXCount++;

                                    //if (point1.X <= point2.X)
                                    //{
                                    //    // point1更靠近左边
                                    //    dxFirst = point1.X - xMin;

                                    //    // 动用一个tX
                                    //    if (additionalTXCount > additionalTXList.Count)
                                    //    {
                                    //        throw new Exception("additionalTXList的数量不足");
                                    //    }
                                    //    double tX = additionalTXList[additionalTXCount - 1];
                                    //    addXFirst = tX * dxFirst;

                                    //    additionalTXCount++;
                                    //    //point1IsLefter = true;
                                    //}
                                    //else
                                    //{
                                    //    // point2更靠近左边
                                    //    dxFirst = point2.X - xMin;

                                    //    // 动用一个tX
                                    //    if (additionalTXCount > additionalTXList.Count)
                                    //    {
                                    //        throw new Exception("additionalTXList的数量不足");
                                    //    }
                                    //    double tX = additionalTXList[additionalTXCount - 1];
                                    //    addXFirst = tX * dxFirst + (point1.X - point2.X);

                                    //    additionalTXCount++;
                                    //    //point1IsLefter = false;
                                    //}

                                    turningPoints.Add(new Point3d(turningPoints.Last().X - addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    double yMax = boundingBox.Max.Y;
                                    double dySecond = yMax - point2.Y;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYSecond = additionalTYList[additionalTYCount - 1];
                                    double addYSecond = tYSecond * dySecond;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y + addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    bool flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 2:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    // 找boundingBox的下界
                                    double yMin = boundingBox.Min.Y;
                                    double dyFirst;
                                    double addYFirst;

                                    // point1一定更靠近下边
                                    dyFirst = point1.Y - yMin;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYFirst = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tYFirst * dyFirst;

                                    additionalTYCount++;

                                    //if (point1.Y <= point2.Y)
                                    //{
                                    //    // point1更靠近下边
                                    //    dyFirst = point1.Y - yMin;

                                    //    // 动用一个tY
                                    //    if (additionalTYCount > additionalTYList.Count)
                                    //    {
                                    //        throw new Exception("additionalTYList的数量不足");
                                    //    }
                                    //    double tY = additionalTYList[additionalTYCount - 1];
                                    //    addYFirst = tY * dyFirst;

                                    //    additionalTYCount++;
                                    //    //point1IsLower = true;
                                    //}
                                    //else
                                    //{
                                    //    // point2更靠近下边
                                    //    dyFirst = point2.Y - yMin;

                                    //    // 动用一个tY
                                    //    if (additionalTYCount > additionalTYList.Count)
                                    //    {
                                    //        throw new Exception("additionalTYList的数量不足");
                                    //    }
                                    //    double tY = additionalTYList[additionalTYCount - 1];
                                    //    addYFirst = tY * dyFirst + (point1.Y - point2.Y);

                                    //    additionalTYCount++;
                                    //    //point1IsLower = false;
                                    //}

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addYFirst, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    xMin = boundingBox.Min.X;
                                    double dxSecond = point2.X - xMin;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXSecond = additionalTXList[additionalTXCount - 1];
                                    double addXSecond = tXSecond * dxSecond;

                                    turningPoints.Add(new Point3d(point2.X - addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 3:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    // 找boundingBox的右界
                                    double xMax = boundingBox.Max.X;
                                    // point1更靠近右边
                                    dxFirst = xMax - point1.X;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXFirst = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tXFirst * dxFirst;

                                    additionalTXCount++;

                                    //if (point1.Y >= point2.Y)
                                    //{
                                    //    // point1更靠近上边
                                    //    dyFirst = yMax - point1.Y;

                                    //    // 动用一个tY
                                    //    if (additionalTYCount > additionalTYList.Count)
                                    //    {
                                    //        throw new Exception("additionalTYList的数量不足");
                                    //    }
                                    //    double tY = additionalTYList[additionalTYCount - 1];
                                    //    addYFirst = tY * dyFirst;

                                    //    additionalTYCount++;
                                    //    //point1IsUpper = true;
                                    //}
                                    //else
                                    //{
                                    //    // point2更靠近上边
                                    //    dyFirst = yMax - point2.Y;

                                    //    // 动用一个tY
                                    //    if (additionalTYCount > additionalTYList.Count)
                                    //    {
                                    //        throw new Exception("additionalTYList的数量不足");
                                    //    }
                                    //    double tY = additionalTYList[additionalTYCount - 1];
                                    //    addYFirst = tY * dyFirst + (point2.Y - point1.Y);

                                    //    additionalTYCount++;
                                    //    //point1IsUpper = false;
                                    //}

                                    turningPoints.Add(new Point3d(turningPoints.Last().X + addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    yMin = boundingBox.Min.Y;
                                    dySecond = point2.Y - yMin;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYSecond = additionalTYList[additionalTYCount - 1];
                                    addYSecond = tYSecond * dySecond;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y - addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 4:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    // 找boundingBox的上界
                                    yMax = boundingBox.Max.Y;

                                    // point1更靠近上边
                                    dyFirst = yMax - point1.Y;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYFirst = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tYFirst * dyFirst;

                                    additionalTYCount++;

                                    //if (point1.X >= point2.X)
                                    //{
                                    //    // point1更靠近右边
                                    //    dxFirst = xMax - point1.X;

                                    //    // 动用一个tX
                                    //    if (additionalTXCount > additionalTXList.Count)
                                    //    {
                                    //        throw new Exception("additionalTXList的数量不足");
                                    //    }
                                    //    double tX = additionalTXList[additionalTXCount - 1];
                                    //    addXFirst = tX * dxFirst;

                                    //    additionalTXCount++;
                                    //    point1IsRighter = true;
                                    //}
                                    //else
                                    //{
                                    //    // point2更靠近右边
                                    //    dxFirst = xMax - point2.X;
                                    //    // 动用一个tX
                                    //    if (additionalTXCount > additionalTXList.Count)
                                    //    {
                                    //        throw new Exception("additionalTXList的数量不足");
                                    //    }
                                    //    double tX = additionalTXList[additionalTXCount - 1];
                                    //    addXFirst = tX * dxFirst + (point2.X - point1.X);

                                    //    additionalTXCount++;
                                    //    point1IsRighter = false;
                                    //}

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addYFirst, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint点
                                    xMax = boundingBox.Max.X;
                                    dxSecond = xMax - point2.X;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXSecond = additionalTXList[additionalTXCount - 1];
                                    addXSecond = tXSecond * dxSecond;

                                    turningPoints.Add(new Point3d(point2.X + addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余的turningPoint点
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                default:
                                    break;
                            }

                            // 移除point1
                            turningPoints.RemoveAt(0);
                        }
                        else
                        {
                            // 由nVectorForBS1到nVectorForBS2是顺时针

                            // 添加point1
                            turningPoints.Add(point1);

                            switch (quadrant)
                            {
                                case 1:
                                    int iter = 0;
                                    #region 添加首个turningPoint
                                    double xMin = boundingBox.Min.X;
                                    double dxFirst = point1.X - xMin;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXFirst = additionalTXList[additionalTXCount - 1];
                                    double addXFirst = tXFirst * dxFirst;

                                    additionalTXCount++;
                                    turningPoints.Add(new Point3d(turningPoints.Last().X - addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    double dyMax = point2.Y - point1.Y;
                                    //double dySecond = point2.Y
                                    // 动用一个tY
                                    if (additionalTYCount>additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYSecond = additionalTYList[additionalTYCount - 1];
                                    double addYSecond = tYSecond * dyMax;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y - addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    bool flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 2:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    double yMin = boundingBox.Min.Y;
                                    double dyFirst = point1.X - yMin;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tYFirst = additionalTYList[additionalTYCount - 1];
                                    double addYFirst = tYFirst * dyFirst;

                                    additionalTYCount++;
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addYFirst, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    double dxMax = point1.X - point2.X;
                                    // 动用一个tX
                                    if (additionalTXCount>additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tXSecond = additionalTXList[additionalTXCount - 1];
                                    double addXSecond = tXSecond * dxMax;

                                    turningPoints.Add(new Point3d(point2.X + addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = point2.Y - turningPoints.Last().Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 3:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    double xMax = boundingBox.Max.X;
                                    dxFirst = xMax - point1.X;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXFirst = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tXFirst * dxFirst;

                                    additionalTXCount++;
                                    turningPoints.Add(new Point3d(turningPoints.Last().X + addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    dyMax = point1.Y - point2.Y;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYSecond = additionalTYList[additionalTYCount - 1];
                                    addYSecond = tYSecond * dyMax;

                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y + addYSecond, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = true;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = turningPoints.Last().X - point2.X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                case 4:
                                    iter = 0;
                                    #region 添加首个turningPoint
                                    double yMax = boundingBox.Max.Y;
                                    dyFirst = yMax - point1.Y;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    tYFirst = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tYFirst * dyFirst;

                                    additionalTYCount++;
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addYFirst, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    #region 构造第二个turningPoint
                                    dxMax = point2.X - point1.X;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    tXSecond = additionalTXList[additionalTXCount - 1];
                                    addXSecond = tXSecond * dxMax;

                                    turningPoints.Add(new Point3d(point2.X - addXSecond, turningPoints.Last().Y, turningPoints.Last().Z));
                                    #endregion
                                    iter++;

                                    #region 添加剩余turningPoint
                                    flag = false;
                                    while (iter < interval - 1)
                                    {
                                        if (flag)
                                        {
                                            double dx = point2.X - turningPoints.Last().X;

                                            // 动用一个tX
                                            if (additionalTXCount > additionalTXList.Count)
                                            {
                                                throw new Exception("additionalTXList的数量不足");
                                            }
                                            double tX = additionalTXList[additionalTXCount - 1];
                                            double addX = tX * dx;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                            additionalTXCount++;
                                            iter++;
                                            flag = false;
                                        }
                                        else
                                        {
                                            double dy = turningPoints.Last().Y - point2.Y;

                                            // 动用一个tY
                                            if (additionalTYCount > additionalTYList.Count)
                                            {
                                                throw new Exception("additionalTYList的数量不足");
                                            }
                                            double tY = additionalTYList[additionalTYCount - 1];
                                            double addY = tY * dy;

                                            turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                            additionalTYCount++;
                                            iter++;
                                            flag = true;
                                        }
                                    }
                                    #endregion

                                    #region 添加最后一个turningPoint点，并移除point1
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, point2.Z));
                                    turningPoints.RemoveAt(0);
                                    #endregion
                                    #endregion
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                    return turningPoints;
                }
            }
            else
            {
                // 如果转折数是偶数
                #region 基础计算
                List<Point3d> turningPoints = new List<Point3d>();

                int quadrant = WhichQuadrant(point1, point2);

                #region 求point1和point2所对应的BS的法线（逆时针90度）
                // 这里不能用bsForPoint，必须用bsForPoint所对应的x轴方向或者y轴方向
                Vector3d vectorForBS1 = new Vector3d(bsForPoint1.To - bsForPoint1.From);
                Vector3d vectorForBS2 = new Vector3d(bsForPoint2.To - bsForPoint2.From);
                Vector3d projectVectorBS1 = CalProjectVector(bsForPoint1.From, bsForPoint1.To);
                Vector3d projectVectorBS2 = CalProjectVector(bsForPoint2.From, bsForPoint2.To);

                // 求在x轴正方向或y轴正方向的投影
                Vector3d vectorForBS1OnProjectVectorBS1 = vectorForBS1 * projectVectorBS1 * projectVectorBS1 / Math.Sqrt(projectVectorBS1.Length);
                Vector3d vectorForBS2OnProjectVectorBS2 = vectorForBS2 * projectVectorBS2 * projectVectorBS2 / Math.Sqrt(projectVectorBS2.Length);
                // 投影后向量的normal，为正X正Y方向
                Vector3d nVectorForBS1 = Normal(vectorForBS1OnProjectVectorBS1);
                Vector3d nVectorForBS2 = Normal(vectorForBS2OnProjectVectorBS2);
                #endregion
                #endregion

                if (linesPassingThroughPoint1 != null)
                {
                    // 如果point1已经被绘制过

                    #region 找到与point1相关的endPoints列表
                    List<Point3d> endPoints = new List<Point3d>();
                    for (int i = 0; i < linesPassingThroughPoint1.Count; i++)
                    {
                        endPoints.Add(linesPassingThroughPoint1[i][1]);
                    }
                    #endregion

                    // 判断两向量是否同向
                    double dot = nVectorForBS1 * nVectorForBS2;
                    if (dot > 0)
                    {
                        // 同向
                        if (nVectorForBS1.Y == 0)
                        {
                            // 如果point1的法向量是x轴方向
                            if (nVectorForBS1.X > 0)
                            {
                                // 法向量为x轴正方向
                                // 添加point
                                turningPoints.Add(point1);

                                #region 对所有的endPoint进行排序(找到resultEndPoint)
                                List<Point3d> sortedEndPoints = endPoints;
                                Point3d resultEndPoint = Point3d.Unset;
                                #endregion

                                bool point1IsLefter;
                                #region 添加首个turningPoint
                                if (point1.X <= point2.X)
                                {
                                    // point1更靠近左边(第四象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 4);

                                    point1IsLefter = true;
                                }
                                else
                                {
                                    // point2更靠近左边(第三象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 3);

                                    point1IsLefter = false;
                                }

                                turningPoints.Add(resultEndPoint);
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = false;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx;
                                        if (point1IsLefter)
                                        {
                                            // point1更靠近左边(第四象限)
                                            dx = turningPoints.Last().X - point2.X;
                                        }
                                        else
                                        {
                                            // point2更靠近左边(第三象限)
                                            dx = turningPoints.Last().X - point1.X;
                                        }
                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                            else
                            {
                                // 法向量为x轴负方向
                                // 添加point
                                turningPoints.Add(point1);

                                #region 对所有的endPoint进行排序(找到resultEndPoint)
                                List<Point3d> sortedEndPoints = endPoints;
                                Point3d resultEndPoint = Point3d.Unset;
                                #endregion

                                bool point1IsLefter;
                                #region 添加首个turningPoint
                                if (point1.X <= point2.X)
                                {
                                    // point1更靠近左边(第一象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 1);

                                    point1IsLefter = true;
                                }
                                else
                                {
                                    // point2更靠近左边(第二象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 2);

                                    point1IsLefter = false;
                                }

                                turningPoints.Add(resultEndPoint);
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx;
                                        if (point1IsLefter)
                                        {
                                            // point1更靠近左边(第一象限)
                                            dx = point1.X - turningPoints.Last().X;
                                        }
                                        else
                                        {
                                            // point2更靠近左边(第二象限)
                                            dx = point2.X - turningPoints.Last().X;
                                        }
                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                        }
                        else
                        {
                            // 如果point1的法向量是y轴方向
                            if (nVectorForBS1.Y > 0)
                            {
                                // 法向量为y轴正方向
                                // 添加point1
                                turningPoints.Add(point1);

                                #region 对所有的endPoint进行排序(找到resultEndPoint)
                                List<Point3d> sortedEndPoints = endPoints;
                                Point3d resultEndPoint = Point3d.Unset;
                                #endregion

                                bool point1IsLower;
                                #region 添加首个turningPoint
                                if (point1.Y <= point2.Y)
                                {
                                    // point1更靠近顶边(第四象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 4);

                                    point1IsLower = true;
                                }
                                else
                                {
                                    // point2更靠近顶边(第一象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 1);

                                    point1IsLower = false;
                                }
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy;
                                        if (point1IsLower)
                                        {
                                            // point1更靠近顶边(第四象限)
                                            dy = turningPoints.Last().Y - point2.Y;
                                        }
                                        else
                                        {
                                            // point2更靠近顶边(第一象限)
                                            dy = turningPoints.Last().Y - point1.Y;
                                        }
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                            else
                            {
                                // 法向量为y轴负方向
                                // 添加point1
                                turningPoints.Add(point1);

                                #region 对所有的endPoint进行排序(找到resultEndPoint)
                                List<Point3d> sortedEndPoints = endPoints;
                                Point3d resultEndPoint = Point3d.Unset;
                                #endregion

                                bool point1IsLower;
                                #region 添加首个turningPoint
                                if (point1.Y <= point2.Y)
                                {
                                    // point1更靠近顶边(第二象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 2);

                                    point1IsLower = true;
                                }
                                else
                                {
                                    // point2更靠近顶边(第三象限)
                                    resultEndPoint = FindExtremeValue(sortedEndPoints, point1, 3);

                                    point1IsLower = false;
                                }
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy;
                                        if (point1IsLower)
                                        {
                                            // point1更靠近顶边(第二象限)
                                            dy = point1.Y - turningPoints.Last().Y;
                                        }
                                        else
                                        {
                                            // point2更靠近顶边(第三象限)
                                            dy = point2.Y - turningPoints.Last().Y;
                                        }
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                        }
                    }
                    else if (dot < 0)
                    {
                        #region 对所有的endPoint进行排序
                        List<Point3d> sortedEndPoints = new List<Point3d>();
                        List<Point3d> pointsAbove = new List<Point3d>();
                        List<Point3d> pointsBelow = new List<Point3d>();
                        switch (quadrant)
                        {
                            case 1:
                                // 第一象限，比较Y值，同时优先找below (Point1+Point2)/2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].Y <= ((point1+point2) / 2).Y)
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 3:
                                // 第三象限，比较Y值，同时优先找above (Point1+Point2)/2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].Y >= ((point1 + point2) / 2).Y)
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 2:
                                // 第二象限，比较X值，同时优先找above (Point1+Point2)/2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].X >= ((point1 + point2) / 2).X)
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                }
                                break;
                            case 4:
                                // 第四象限，比较X值，同时优先找below (Point1+Point2)/2的
                                for (int i = 0; i < endPoints.Count; i++)
                                {
                                    if (endPoints[i].X <= ((point1 + point2) / 2).X)
                                    {
                                        pointsBelow.Add(endPoints[i]);
                                    }
                                    else
                                    {
                                        pointsAbove.Add(endPoints[i]);
                                    }
                                }
                                break;
                            default:
                                break;
                        }
                        #endregion

                        List<Point3d> sortedPointsAbove = new List<Point3d>();
                        List<Point3d> sortedPointsBelow = new List<Point3d>();

                        switch (quadrant)
                        {
                            case 1:
                                int iter = 0;
                                bool flag;
                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsBelow = new List<Point3d>();
                                    // 第一象限，比较Y值，同时优先找below (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].Y <= ((point1 + point2) / 2).Y)
                                        {
                                            pointsBelow.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;

                                    // 1象限时，根据Y值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, true);

                                    if (nVectorForBS1.Y > 0)
                                    {
                                        // nVectorForBS1.Y > 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsBelow.Last());
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = point2.X - turningPoints.Last().X;

                                                // 动用一个tX
                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = point2.Y - turningPoints.Last().Y;

                                                // 动用一个tY
                                                if (additionalTYCount > additionalTYList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在x方向上有投影
                                        // 就要point2的x值，和目前最后一个turningPoint的y值
                                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // nVectorForBS1.Y < 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsBelow.Last().X, 2 * point1.Y - sortedPointsBelow.Last().Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即X值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        // 去参考同向的计算方式
                                        throw new Exception("偶数间隔，反向，case1，vectorForBS1OnProjectVectorBS1.Y == 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        #endregion
                                    }
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsBelow = new List<Point3d>();
                                    // 第一象限，比较X值，同时优先找below (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].X <= ((point1 + point2) / 2).X)
                                        {
                                            pointsBelow.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;

                                    // 1象限时，根据X值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, false);

                                    if (nVectorForBS1.X > 0)
                                    {
                                        // 法向量为x轴正方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.X的点中，找到最大X值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsBelow.Last());
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.X的点，直接构造新的点（即X值等于(point1+point2)/2.X的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(((point1 + point2) / 2).X, point1.Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = point2.X - turningPoints.Last().X;

                                                // 动用一个tX
                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = point2.Y - turningPoints.Last().Y;

                                                // 动用一个tY
                                                if (additionalTYCount > additionalTYList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在y方向上有投影
                                        // 就要目前最后一个turningPoint的x值，和point2的y值
                                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                        #endregion

                                    }
                                    else
                                    {
                                        // 法向量为x轴负方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.X的点中，找到最大X值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsBelow.Last().X, 2 * point1.Y - sortedPointsBelow.Last().Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.X的点，直接构造新的点（即X值等于(point1+point2)/2.X的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        // 去参考同向的计算方式
                                        throw new Exception("偶数间隔，反向，case1，vectorForBS1OnProjectVectorBS1.Y != 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }
                                }
                                #endregion
                                break;
                            case 2:
                                iter = 0;
                                //flag;
                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsBelow = new List<Point3d>();
                                    // 第一象限，比较Y值，同时优先找below (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].Y <= ((point1 + point2) / 2).Y)
                                        {
                                            pointsBelow.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;

                                    // 2象限时，根据Y值进行排序
                                    sortedPointsBelow = pointsBelow;
                                    BubbleSortPoint3d(sortedPointsBelow, true);

                                    if (nVectorForBS1.Y > 0)
                                    {
                                        // nVectorForBS1.Y > 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsBelow.Last());
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = turningPoints.Last().X - point2.X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = point2.Y - turningPoints.Last().Y;

                                                if (additionalTYCount > additionalTYList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在x方向上有投影
                                        // 就要point2的x值，和目前最后一个turningPoint的y值
                                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // nVectorForBS1.Y < 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsBelow.Last().X, 2 * point1.Y - sortedPointsBelow.Last().Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即X值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        throw new Exception("偶数间隔，反向，case2，vectorForBS1OnProjectVectorBS1.Y == 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }

                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsAbove = new List<Point3d>();
                                    // 第一象限，比较X值，同时优先找above (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].X >= ((point1 + point2) / 2).X)
                                        {
                                            pointsAbove.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;

                                    // 2象限时，根据X值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, true);

                                    if (nVectorForBS1.X > 0)
                                    {
                                        // 法向量为x轴正方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsAbove[0]);
                                        }
                                        else
                                        {
                                            // 如果没有above (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = turningPoints.Last().X - point2.X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = point2.Y - turningPoints.Last().Y;

                                                if (additionalTYCount > additionalTYList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在y方向上有投影
                                        // 就要目前最后一个turningPoint的x值，和point2的y值
                                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // 法向量为x轴负方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsAbove[0].X, 2 * point1.Y - sortedPointsAbove[0].Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有above (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        // 去参考同向的计算方式
                                        throw new Exception("偶数间隔，反向，case2，vectorForBS1OnProjectVectorBS1.Y != 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }
                                }
                                #endregion

                                break;
                            case 3:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsAbove = new List<Point3d>();
                                    // 第一象限，比较Y值，同时优先找above (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].Y >= ((point1 + point2) / 2).Y)
                                        {
                                            pointsAbove.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;

                                    // 3象限时，根据Y值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, true);

                                    if (nVectorForBS1.Y > 0)
                                    {
                                        // nVectorForBS1.Y > 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsAbove[0]);
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = turningPoints.Last().X - point2.X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = turningPoints.Last().Y - point2.Y;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在x方向上有投影
                                        // 就要point2的x值，和目前最后一个turningPoint的y值
                                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // nVectorForBS1.Y < 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsAbove[0].X, 2 * point1.Y - sortedPointsAbove[0].Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即X值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        throw new Exception("偶数间隔，反向，case3，vectorForBS1OnProjectVectorBS1.Y == 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }

                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsAbove = new List<Point3d>();
                                    // 第一象限，比较X值，同时优先找above (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].X >= ((point1 + point2) / 2).X)
                                        {
                                            pointsAbove.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                    if (nVectorForBS1.X > 0)
                                    {
                                        // 法向量为x轴正方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsAbove[0]);
                                        }
                                        else
                                        {
                                            // 如果没有above (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = turningPoints.Last().X - point2.X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = turningPoints.Last().Y - point2.Y;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在y方向上有投影
                                        // 就要目前最后一个turningPoint的x值，和point2的y值
                                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // 法向量为x轴负方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsAbove[0].X, 2 * point1.Y - sortedPointsAbove[0].Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有above (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        // 去参考同向的计算方式
                                        throw new Exception("偶数间隔，反向，case3，vectorForBS1OnProjectVectorBS1.Y != 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }
                                }
                                #endregion
                                break;
                            case 4:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsAbove = new List<Point3d>();
                                    // 第一象限，比较Y值，同时优先找above (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].Y >= ((point1 + point2) / 2).Y)
                                        {
                                            pointsAbove.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;

                                    // 4象限时，根据Y值进行排序
                                    sortedPointsAbove = pointsAbove;
                                    BubbleSortPoint3d(sortedPointsAbove, true);

                                    if (nVectorForBS1.Y > 0)
                                    {
                                        // nVectorForBS1.Y > 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从above (point1+point2)/2.Y的点中，找到最小Y值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsAbove[0]);
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即Y值等于(point1+point2)/2.Y的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(point1.X, ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = point2.X - turningPoints.Last().X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = turningPoints.Last().Y - point2.Y;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在x方向上有投影
                                        // 就要point2的x值，和目前最后一个turningPoint的y值
                                        turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // nVectorForBS1.Y < 0的情况
                                        #region 添加第一个turningPoint
                                        if (sortedPointsAbove.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.Y的点中，找到最大Y值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsAbove[0].X, 2 * point1.Y - sortedPointsAbove[0].Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.Y的点，直接构造新的点（即X值等于(point1+point2)/2.Y的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        throw new Exception("偶数间隔，反向，case4，vectorForBS1OnProjectVectorBS1.Y == 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }

                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    #region 对所有的endPoint进行排序
                                    pointsBelow = new List<Point3d>();
                                    // 第一象限，比较X值，同时优先找below (Point1+Point2)/2的
                                    for (int i = 0; i < endPoints.Count; i++)
                                    {
                                        if (endPoints[i].X <= ((point1 + point2) / 2).X)
                                        {
                                            pointsBelow.Add(endPoints[i]);
                                        }
                                    }
                                    #endregion

                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                    if (nVectorForBS1.X > 0)
                                    {
                                        // 法向量为x轴正方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.X的点中，找到最大X值的那个点，添加到turningPoint中
                                            turningPoints.Add(sortedPointsBelow.Last());
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.X的点，直接构造新的点（即X值等于(point1+point2)/2.X的点），添加到turningPoint中
                                            turningPoints.Add(new Point3d(((point1 + point2) / 2).X, point1.Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        do
                                        {
                                            if (flag)
                                            {
                                                double dx = point2.X - turningPoints.Last().X;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTXList的数量不足");
                                                }
                                                double tX = additionalTXList[additionalTXCount - 1];
                                                double addX = tX * dx;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                                additionalTXCount++;
                                                flag = false;
                                            }
                                            else
                                            {
                                                double dy = turningPoints.Last().Y - point2.Y;

                                                if (additionalTXCount > additionalTXList.Count)
                                                {
                                                    throw new Exception("additionalTYList的数量不足");
                                                }
                                                double tY = additionalTYList[additionalTYCount - 1];
                                                double addY = tY * dy;

                                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                                additionalTYCount++;
                                                flag = true;
                                            }

                                            iter++;
                                        } while (iter < interval - 1);
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式
                                        // 如果原来的bs是在y方向上有投影
                                        // 就要目前最后一个turningPoint的x值，和point2的y值
                                        turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                        #endregion
                                    }
                                    else
                                    {
                                        // 法向量为x轴负方向
                                        #region 添加第一个turningPoint
                                        if (sortedPointsBelow.Count != 0)
                                        {
                                            // 从below (point1+point2)/2.X的点中，找到最大X值的那个点，做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - sortedPointsBelow.Last().X, 2 * point1.Y - sortedPointsBelow.Last().Y, point1.Z));
                                        }
                                        else
                                        {
                                            // 如果没有below (point1+point2)/2.X的点，直接构造新的点（即X值等于(point1+point2)/2.X的点），做镜像点，添加到turningPoint中
                                            turningPoints.Add(new Point3d(2 * point1.X - ((point1 + point2) / 2).X, 2 * point1.Y - ((point1 + point2) / 2).Y, point1.Z));
                                        }
                                        #endregion
                                        iter++;

                                        #region 循环计算turningPoint
                                        // 去参考同向的计算方式
                                        throw new Exception("偶数间隔，反向，case2，vectorForBS1OnProjectVectorBS1.Y != 0，nVectorForBS1.Y < 0，代码未完成");
                                        #endregion

                                        #region 选择最后一个turningPoint的计算方式

                                        #endregion
                                    }
                                }
                                #endregion
                                break;
                            default:
                                break;
                        }
                    }
                    else
                    {
                        return null;
                        throw new Exception("转折数为偶数时，两向量不平行(即垂直)");
                    }

                    return turningPoints;
                }
                else
                {
                    // 如果point1没有被绘制过，即这是全局的第一次绘制
                    
                    // 判断两向量是否同向
                    double dot = nVectorForBS1 * nVectorForBS2;
                    if (dot > 0)
                    {
                        // 同向

                        if (nVectorForBS1.Y == 0)
                        {
                            // 如果point1的法向量是x轴方向
                            if (nVectorForBS1.X > 0)
                            {
                                // 法向量为x轴正方向
                                // 添加point
                                turningPoints.Add(point1);

                                #region 添加首个turningPoint
                                // 找boundingbox的右界
                                double xMax = boundingBox.Max.X;
                                double dxFirst;
                                double addXFirst;
                                bool point1IsLefter;
                                if (point1.X <= point2.X)
                                {
                                    // point1更靠近左边
                                    dxFirst = xMax - point2.X;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tX = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tX * dxFirst + (point2.X - point1.X);

                                    point1IsLefter = true;
                                }
                                else
                                {
                                    // point2更靠近左边
                                    dxFirst = xMax - point1.X;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tX = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tX * dxFirst;

                                    point1IsLefter = false;
                                }

                                turningPoints.Add(new Point3d(turningPoints.Last().X + addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx;
                                        if (point1IsLefter)
                                        {
                                            dx = turningPoints.Last().X - point2.X;
                                        }
                                        else
                                        {
                                            dx = turningPoints.Last().X - point1.X;
                                        }
                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                            else
                            {
                                // 法向量为x轴负方向
                                // 添加point
                                turningPoints.Add(point1);

                                #region 添加首个turningPoint
                                // 找boundingbox的左界
                                double xMin = boundingBox.Min.X;
                                double dxFirst;
                                double addXFirst;
                                bool point1IsLefter;
                                if (point1.X <= point2.X)
                                {
                                    // point1更靠近左边
                                    dxFirst = point1.X - xMin;

                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tX = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tX * dxFirst;

                                    point1IsLefter = true;
                                }
                                else
                                {
                                    // point2更靠近左边
                                    dxFirst = point2.X - xMin;
                                    // 动用一个tX
                                    if (additionalTXCount > additionalTXList.Count)
                                    {
                                        throw new Exception("additionalTXList的数量不足");
                                    }
                                    double tX = additionalTXList[additionalTXCount - 1];
                                    addXFirst = tX * dxFirst + (point1.X - point2.X);

                                    point1IsLefter = false;
                                }

                                turningPoints.Add(new Point3d(turningPoints.Last().X - addXFirst, turningPoints.Last().Y, turningPoints.Last().Z));
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx;
                                        if (point1IsLefter)
                                        {
                                            dx = point1.X - turningPoints.Last().X;
                                        }
                                        else
                                        {
                                            dx = point2.X - turningPoints.Last().X;
                                        }
                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                        }
                        else
                        {
                            // 如果point1的法向量是y轴方向
                            if (nVectorForBS1.Y > 0)
                            {
                                // 法向量为y轴正方向
                                // 添加point1
                                turningPoints.Add(point1);

                                #region 添加首个turningPoint
                                // 找boundingbox的上界
                                double yMax = boundingBox.Max.Y;
                                double dyFirst;
                                double addYFirst;
                                bool point1IsLower;
                                if (point1.Y <= point2.Y)
                                {
                                    // point2更靠近顶边
                                    dyFirst = yMax - point2.Y;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tY = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tY * dyFirst + (point2.Y - point1.Y);

                                    point1IsLower = true;
                                }
                                else
                                {
                                    // point1更靠近顶边
                                    dyFirst = yMax - point1.Y;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tY = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tY * dyFirst;

                                    point1IsLower = false;
                                }

                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addYFirst, turningPoints.Last().Z));
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy;
                                        if (point1IsLower)
                                        {
                                            dy = turningPoints.Last().Y - point2.Y;
                                        }
                                        else
                                        {
                                            dy = turningPoints.Last().Y - point1.Y;
                                        }
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                            else
                            {
                                // 法向量为y轴负方向
                                // 添加point1
                                turningPoints.Add(point1);

                                #region 添加首个turningPoint
                                // 找boundingbox的下界
                                double yMin = boundingBox.Min.Y;
                                double dyFirst;
                                double addYFirst;
                                bool point1IsLower;
                                if (point1.Y <= point2.Y)
                                {
                                    // point1更靠近底边 
                                    dyFirst = point1.Y - yMin;

                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tY = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tY * dyFirst;

                                    point1IsLower = true;
                                }
                                else
                                {
                                    // point2更靠近底边
                                    dyFirst = point2.Y - yMin;
                                    // 动用一个tY
                                    if (additionalTYCount > additionalTYList.Count)
                                    {
                                        throw new Exception("additionalTYList的数量不足");
                                    }
                                    double tY = additionalTYList[additionalTYCount - 1];
                                    addYFirst = tY * dyFirst + (point1.Y - point2.Y);

                                    point1IsLower = false;
                                }

                                turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addYFirst, turningPoints.Last().Z));
                                #endregion

                                #region 添加剩余turningPoint
                                int iter = 0;
                                bool flag = true;
                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy;
                                        if (point1IsLower)
                                        {
                                            dy = point1.Y - turningPoints.Last().Y;
                                        }
                                        else
                                        {
                                            dy = point2.Y - turningPoints.Last().Y;
                                        }
                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);
                            }
                        }
                    }
                    else if (dot < 0)
                    {
                        // 反向
                        switch (quadrant)
                        {
                            case 1:
                                int iter = 0;
                                bool flag;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        // 动用一个tX
                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        // 动用一个tY
                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 2:
                                iter = 0;
                                //flag;
                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = point2.Y - turningPoints.Last().Y;

                                        if (additionalTYCount > additionalTYList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y + addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 3:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = turningPoints.Last().X - point2.X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X - addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            case 4:
                                iter = 0;

                                #region 选择bool flag的初始值
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了bool flag的初始值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要先添加y值，来得到turningPoint
                                    flag = false;
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要先添加x值，来得到turningPoint
                                    flag = true;
                                }
                                #endregion

                                #region 循环计算turningPoint
                                // 将point1添加到turningPoints列表中，方便循环计算的开始
                                turningPoints.Add(point1);

                                do
                                {
                                    if (flag)
                                    {
                                        double dx = point2.X - turningPoints.Last().X;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTXList的数量不足");
                                        }
                                        double tX = additionalTXList[additionalTXCount - 1];
                                        double addX = tX * dx;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X + addX, turningPoints.Last().Y, turningPoints.Last().Z));

                                        additionalTXCount++;
                                        flag = false;
                                    }
                                    else
                                    {
                                        double dy = turningPoints.Last().Y - point2.Y;

                                        if (additionalTXCount > additionalTXList.Count)
                                        {
                                            throw new Exception("additionalTYList的数量不足");
                                        }
                                        double tY = additionalTYList[additionalTYCount - 1];
                                        double addY = tY * dy;

                                        turningPoints.Add(new Point3d(turningPoints.Last().X, turningPoints.Last().Y - addY, turningPoints.Last().Z));

                                        additionalTYCount++;
                                        flag = true;
                                    }

                                    iter++;
                                } while (iter < interval - 1);
                                #endregion

                                #region 选择最后一个turningPoint的计算方式
                                // 要考虑是在x轴方向的投影，还是在y轴方向上的投影，这决定了最后一个turningPoint的值
                                if (vectorForBS1OnProjectVectorBS1.Y == 0)
                                {
                                    // 如果原来的bs是在x方向上有投影
                                    // 就要point2的x值，和目前最后一个turningPoint的y值
                                    turningPoints.Add(new Point3d(point2.X, turningPoints.Last().Y, turningPoints.Last().Z));
                                }
                                else
                                {
                                    // 如果原来的bs是在y方向上有投影
                                    // 就要目前最后一个turningPoint的x值，和point2的y值
                                    turningPoints.Add(new Point3d(turningPoints.Last().X, point2.Y, turningPoints.Last().Z));
                                }
                                #endregion

                                // 移除point1
                                turningPoints.RemoveAt(0);

                                break;
                            default:
                                break;
                        }
                    }
                    else
                    {
                        return null;
                        throw new Exception("转折数为偶数时，两向量不平行(即垂直)");
                    }

                    return turningPoints;
                }
            }
        }

        private Vector3d CalProjectVector(Point3d point1, Point3d point2)
        {
            Vector3d projectVector;

            // 要考虑原来的向量是属于哪个象限的问题
            int quadrant = WhichQuadrant(point1, point2);

            if (point2.X - point1.X != 0)
            {
                double k = (point2.Y - point1.Y) / (point2.X - point1.X);
                if (k > 0 && k <= 1)
                {
                    // 斜率小于45度时
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k > 1)
                {
                    if (quadrant == 1)
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                    else//quadrant == 3
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                }
                else if (k == 0)
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else if (k >= -1 && k < 0)
                {
                    // 斜率小于45度时
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(1, 0, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(-1, 0, 0);
                    }
                }
                else
                {
                    if (quadrant == 4)
                    {
                        projectVector = new Vector3d(0, -1, 0);
                    }
                    else//quadrant == 2
                    {
                        projectVector = new Vector3d(0, 1, 0);
                    }
                }
            }
            else
            {
                // 即y轴正方向
                if (quadrant == 1)
                {
                    projectVector = new Vector3d(0, 1, 0);
                }
                else//quadrant == 3
                {
                    projectVector = new Vector3d(0, -1, 0);
                }
            }

            return projectVector;
        }

        private int WhichQuadrant(Point3d point1,Point3d point2)
        {
            double dx = point2.X - point1.X;
            double dy = point2.Y - point1.Y;

            if (dx == 0)
            {
                // y轴
                if (dy > 0)
                {
                    // y轴正方向，并入第一象限的情况计算
                    return 1;
                }
                else
                {
                    // y轴反方向，并入第三象限的情况计算
                    return 3;
                }
            }
            if (dy == 0)
            {
                // x轴
                if (dx > 0)
                {
                    // x轴正方向，并入第四象限计算
                    return 4;
                }
                else
                {
                    // x轴反方向，并入第二象限计算
                    return 2;
                }
            }

            if (dy / dx > 0 && dx > 0)
            {
                return 1;
            }
            else if (dy / dx > 0 && dx < 0)
            {
                return 3;
            }
            else if (dy / dx < 0 && dx < 0)
            {
                return 2;
            }
            else
            {
                return 4;
            }
        }

        /// <summary>
        /// 从小到大进行冒泡
        /// </summary>
        /// <param name="needToSort"></param>
        /// <param name="flag">true是按Y排序，false是按X排序</param>
        /// <returns></returns>
        private bool BubbleSortPoint3d(List<Point3d> needToSort, bool flag)
        {
            if (needToSort.Count < 2)
            {
                return false;
            }

            if (flag)
            {
                for (int i = 0; i < needToSort.Count - 1; i++)
                {
                    for (int j = 0; j < needToSort.Count - 1; j++)
                    {
                        if (needToSort[j].Y > needToSort[j + 1].Y)
                        {
                            Point3d temp = needToSort[j + 1];
                            needToSort[j + 1] = needToSort[j];
                            needToSort[j] = temp;
                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < needToSort.Count - 1; i++)
                {
                    for (int j = 0; j < needToSort.Count - 1; j++)
                    {
                        if (needToSort[j].X > needToSort[j + 1].X)
                        {
                            Point3d temp = needToSort[j + 1];
                            needToSort[j + 1] = needToSort[j];
                            needToSort[j] = temp;
                        }
                    }
                }
            }
            return true;
        }

        private Point3d FindExtremeValue(List<Point3d> points, Point3d point1, int options)
        {
            Point3d resultPoint = Point3d.Unset;
            switch (options)
            {
                case 1:
                    // 第一象限，找最靠上的或者最靠左的
                    resultPoint = points[0];
                    Vector3d resultVector = new Vector3d(points[0] - point1);
                    double resultDistance = resultVector.Length;
                    List<Vector3d> vectors = new List<Vector3d>();
                    for (int i = 1; i < points.Count; i++)
                    {
                        Vector3d vector = new Vector3d(points[i] - point1);
                        double distance = vector.Length;
                        if (distance > resultDistance)
                        {
                            resultDistance = distance;
                            resultPoint = points[i];
                        }
                    }
                    break;
                case 2:
                    // 第二象限，找最靠下的或者最靠左的
                    resultPoint = points[0];
                    resultVector = new Vector3d(points[0] - point1);
                    resultDistance = resultVector.Length;
                    vectors = new List<Vector3d>();
                    for (int i = 1; i < points.Count; i++)
                    {
                        Vector3d vector = new Vector3d(points[i] - point1);
                        double distance = vector.Length;
                        if (distance > resultDistance)
                        {
                            resultDistance = distance;
                            resultPoint = points[i];
                        }
                    }
                    break;
                case 3:
                    // 第三象限，找最靠下的或者最靠右的
                    resultPoint = points[0];
                    resultVector = new Vector3d(points[0] - point1);
                    resultDistance = resultVector.Length;
                    vectors = new List<Vector3d>();
                    for (int i = 1; i < points.Count; i++)
                    {
                        Vector3d vector = new Vector3d(points[i] - point1);
                        double distance = vector.Length;
                        if (distance > resultDistance)
                        {
                            resultDistance = distance;
                            resultPoint = points[i];
                        }
                    }
                    break;
                case 4:
                    // 第四象限，找最靠上的或者最靠右的
                    resultPoint = points[0];
                    resultVector = new Vector3d(points[0] - point1);
                    resultDistance = resultVector.Length;
                    vectors = new List<Vector3d>();
                    for (int i = 1; i < points.Count; i++)
                    {
                        Vector3d vector = new Vector3d(points[i] - point1);
                        double distance = vector.Length;
                        if (distance > resultDistance)
                        {
                            resultDistance = distance;
                            resultPoint = points[i];
                        }
                    }
                    break;
                default:
                    break;
            }

            return resultPoint;
        }

        private bool IsConvex(List<Point3d> pointList)
        {
            if (pointList.Count<3)
            {
                return false;
            }

            for (int i = 0; i < pointList.Count; i++)
            {
                double dx0 = pointList[(i + 1) % pointList.Count].X - pointList[i].X;
                double dy0 = pointList[(i + 1) % pointList.Count].Y - pointList[i].Y;
                double dx1 = pointList[(i + 2) % pointList.Count].X - pointList[i].X;
                double dy1 = pointList[(i + 2) % pointList.Count].Y - pointList[i].Y;
                // crossProduct > 0 表示points是按逆时针输出的；< 0，顺时针
                double crossProduct = dx0 * dy1 - dx1 * dy0;

                double flag = 0;

                if (crossProduct != 0)
                {
                    // 说明异号，有个角大于180度
                    if (crossProduct * flag < 0)
                    {
                        return false;
                    }
                    else
                    {
                        flag = crossProduct;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// 计算逆时针旋转后的向量
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="rad"></param>
        /// <returns></returns>
        private Vector3d Rotate(Vector3d vector, double rad)
        {
            double x = vector.X * Math.Cos(rad) - vector.Y * Math.Sin(rad);
            double y = vector.X * Math.Sin(rad) + vector.Y * Math.Cos(rad);
            return new Vector3d(x, y, vector.Z);
        }

        /// <summary>
        /// 计算向量逆时针旋转九十度的单位法向量
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        private Vector3d Normal(Vector3d vector)
        {
            double L = Math.Sqrt(vector * vector);
            return new Vector3d(-vector.Y / L, vector.X / L, vector.Z);
        }

        /// <summary>
        /// 找到两个射线的交点
        /// </summary>
        /// <param name="p"></param>
        /// <param name="v"></param>
        /// <param name="q"></param>
        /// <param name="w"></param>
        /// <returns></returns>
        private Point3d GetLineIntersection(Point3d p,Vector3d v, Point3d q, Vector3d w)
        {
            Vector3d u = p - q;
            if (CrossProduct(v, w) == 0)
            {
                return Point3d.Unset;
            }
            else
            {
                double t = CrossProduct(w, u) / CrossProduct(v, w);
                return p + v * t;
            }
        }

        /// <summary>
        /// 不考虑Z值的叉积
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        private double CrossProduct(Vector3d a, Vector3d b)
        {
            return a.X * b.Y - a.Y * b.X;
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            //base.DrawViewportWires(args);
            for (int i = 0; i < VolumeJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.Red, Color.White, Color.White);
                args.Display.EnableDepthTesting(true);
            }
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            //base.DrawViewportWires(args);
            for (int i = 0; i < VolumeJunctionTextDot.Count; i++)
            {
                args.Display.EnableDepthTesting(false);
                args.Display.DrawDot(VolumeJunctionTextDot[i], Color.Red, Color.White, Color.White);
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
            get { return new Guid("bb14004b-a69c-402c-ab3c-b36d211bdce5"); }
        }
    }
}