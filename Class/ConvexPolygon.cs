using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace VolumeGeneratorBasedOnGraph.Class
{
    /// <summary>
    /// 凸包类
    /// </summary>
    public class ConvexPolygon : List<int>
    {
        /// <summary>
        /// 用一个index的列表来存储和表示凸包ConvexPolygon，它不像Polyline那样需要在末尾再存上[0]
        /// </summary>
        private List<int> indices;

        /// <summary>
        /// Debug用
        /// </summary>
        public string Description
        {
            get
            {
                return string.Join<int>(",", indices);
            }
        }

        /// <summary>
        /// 判断凸包ConvexPolygon是否已经是三角形了（三角形也是一种凸包）
        /// </summary>
        public bool IsTriangle
        {
            get
            {
                return indices.Count == 3;
            }
        }

        /// <summary>
        /// 判断凸包ConvexPolygon是否已经是四边形了（四边形也是一种凸包）
        /// </summary>
        public bool IsQuadrangle
        {
            get
            {
                return indices.Count == 4;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="vertex_indices"></param>
        public ConvexPolygon(List<int> vertex_indices)
        {
            indices = new List<int>();
            indices = vertex_indices;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="PolygonalPolyline"></param>
        public ConvexPolygon(Polyline PolygonalPolyline)
        {
            if (PolygonalPolyline.IsClosed)
            {
                indices = new List<int>();
                // 因为闭合的polyline会有首尾重复的一个点
                for (int i = 0; i < PolygonalPolyline.Count - 1; i++)
                {
                    indices.Add(i);
                }
            }
        }

        /// <summary>
        /// 将凸包ConvexPolygon转化为Polyline类，凸包ConvexPolygon只有点序号，通过输入实际的点，来形成闭合的Polyline
        /// </summary>
        /// <param name="Vertices"></param>
        /// <returns></returns>
        public Polyline ToPolyline(IEnumerable<Point3d> Vertices)
        {
            List<Point3d> list = new List<Point3d>();

            foreach (int num in indices)
            {
                list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, num));
            }

            // 保证生成的Polyline的首尾闭合，将起始点又放在队尾
            list.Add(Enumerable.ElementAtOrDefault<Point3d>(Vertices, indices[0]));

            return new Polyline(list);
        }

        /// <summary>
        /// 将凸包ConvexPolygon三角形化
        /// </summary>
        /// <returns></returns>
        public List<List<MeshFace>> TriangulateConvexPolygon()
        {
            if (this == null)
            {
                return null;
            }


            List<List<MeshFace>> LoLConvexPolygonMeshFace = new List<List<MeshFace>>();

            // 如果对象是三角化的
            if (IsTriangle)
            {
                LoLConvexPolygonMeshFace.Add(new List<MeshFace>());

                MeshFace currentConvexPolygon = new MeshFace(indices[0], indices[1], indices[2]);
                LoLConvexPolygonMeshFace[0].Add(currentConvexPolygon);

                return LoLConvexPolygonMeshFace;
            }
            // 如果对象不是三角化的
            else
            {
                // 就从第二个顶点开始，对对象进行分割
                for (int i = 2; i < indices.Count; i++)
                {
                    List<ConvexPolygon> subConvexPolygonList = DivideConvexPolygon(i);

                    // 对于分割后生成的子对象列表来说，如果只有一个元素
                    if (subConvexPolygonList.Count == 1)
                    {
                        ConvexPolygon subConvexPolygon = subConvexPolygonList[0];
                        List<List<MeshFace>> LoLSubConvexPolygonMeshFace = subConvexPolygon.TriangulateConvexPolygon();

                        foreach (List<MeshFace> meshFace in LoLSubConvexPolygonMeshFace)
                        {
                            List<MeshFace> list = new List<MeshFace>();
                            MeshFace currentConvexPolygonMeshFace = new MeshFace(indices[0], indices[1], indices[i]);

                            list.Add(currentConvexPolygonMeshFace);
                            list.AddRange(meshFace);

                            LoLConvexPolygonMeshFace.Add(list);
                        }
                    }
                    else
                    {
                        ConvexPolygon behindDivideConvexPolygon = subConvexPolygonList[0];
                        ConvexPolygon frontDivideConvexPolygon = subConvexPolygonList[1];

                        List<List<MeshFace>> LoLBehindDivideConvexPolygonMeshFace = behindDivideConvexPolygon.TriangulateConvexPolygon();
                        List<List<MeshFace>> LoLFrontDivideConvexPolygonMeshFace = frontDivideConvexPolygon.TriangulateConvexPolygon();

                        foreach (List<MeshFace> behindMeshFace in LoLBehindDivideConvexPolygonMeshFace)
                        {
                            foreach (List<MeshFace> frontMeshFace in LoLFrontDivideConvexPolygonMeshFace)
                            {
                                List<MeshFace> list = new List<MeshFace>();
                                MeshFace currentConvexPolygonMeshFace = new MeshFace(indices[0], indices[1], indices[i]);
                                list.Add(currentConvexPolygonMeshFace);
                                list.AddRange(behindMeshFace);
                                list.AddRange(frontMeshFace);
                                LoLConvexPolygonMeshFace.Add(list);
                            }
                        }
                    }
                }
            }
            return LoLConvexPolygonMeshFace;
        }

        /// <summary>
        /// 在要分割位置的点处，分割整个大的Polygon，形成一个或两个子Polygon(behindConvexPolygon和frontConvexPolygon)
        /// </summary>
        /// <param name="divideIndex">要分割位置的点的序号</param>
        /// <returns></returns>
        public List<ConvexPolygon> DivideConvexPolygon(int divideIndex)
        {
            // 分割后得到的两个ConvexPolygon的列表
            List<ConvexPolygon> dividedConvexPolygonList;

            // 当分割位置的点的序号是0（表示polyline的第一个点），1（polyline的第二个点）或者indices.Count - 1（表示polyline的最后一个点）时，不能分割出两个新的ConvexPolygon
            if (divideIndex <= 1 | divideIndex > indices.Count - 1)
            {
                dividedConvexPolygonList = null;
            }
            else
            {
                // 从分割点往后，包括分割点和整个大的Polygon起点的这些点，所形成的子Polygon
                List<int> behindDividedIndexList = new List<int>();
                // 从分割点往前，包括分割点
                List<int> frontDividedIndexList = new List<int>();

                // 包括分割点并且往后
                for (int i = divideIndex; i < indices.Count; i++)
                {
                    behindDividedIndexList.Add(indices[i]);
                }
                // 包括整个大的Polygon起点
                behindDividedIndexList.Add(indices[0]);

                // 从第二个点开始，到分割点（包括分割点），所形成的的子Polygon
                for (int i = 1; i <= divideIndex; i++)
                {
                    frontDividedIndexList.Add(indices[i]);
                }

                ConvexPolygon behindDivideConvexPolygon = null;
                ConvexPolygon frontDivideConvexPolygon = null;

                // 大于三个点，才形成ConvexPolygon
                if (behindDividedIndexList.Count >= 3)
                {
                    // polygonFace -> behindDivideConvexPolygon
                    behindDivideConvexPolygon = new ConvexPolygon(behindDividedIndexList);
                }
                if (frontDividedIndexList.Count >= 3)
                {
                    // polygonFace2 -> frontDivideConvexPolygon
                    frontDivideConvexPolygon = new ConvexPolygon(frontDividedIndexList);
                }

                // 输出最后分割出来符合顶点数>=3的子ConvexPolygon
                if (behindDivideConvexPolygon != null & frontDivideConvexPolygon != null)
                {
                    dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                    {
                            behindDivideConvexPolygon,
                            frontDivideConvexPolygon
                    });
                }
                else
                {
                    if (behindDivideConvexPolygon == null)
                    {
                        dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                        {
                                frontDivideConvexPolygon
                        });
                    }
                    else if (frontDivideConvexPolygon == null)
                    {
                        dividedConvexPolygonList = Enumerable.ToList<ConvexPolygon>(new ConvexPolygon[]
                        {
                                behindDivideConvexPolygon
                        });
                    }
                    else
                    {
                        dividedConvexPolygonList = null;
                    }
                }
            }
            return dividedConvexPolygonList;
        }
    }


}