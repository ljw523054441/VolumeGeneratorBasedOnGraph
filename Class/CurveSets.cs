using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Linq;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    public class CurveSets
    {
		private List<Curve> m_set_a;

		private List<Curve> m_set_b;

		public CurveSets(Curve crv)
		{
			this.m_set_a = new List<Curve>();
			this.m_set_b = new List<Curve>();
			this.m_set_a.Add(crv);
		}

		public CurveSets(IEnumerable<Curve> crvs)
		{
			this.m_set_a = new List<Curve>();
			this.m_set_b = new List<Curve>();
			this.m_set_a.AddRange(crvs);
		}

		public void AppendNewCurve(Curve crv)
		{
			if (crv != null)
			{
				this.m_set_b.Add(crv);
			}
		}

		public void AppendNewCurves(IEnumerable<Curve> crvs)
		{
			if (crvs != null)
			{
				foreach (Curve crv in crvs)
				{
					this.AppendNewCurve(crv);
				}
			}
		}

		public IList<Curve> Generation0
		{
			get
			{
				return this.m_set_a;
			}
		}

		public IList<Curve> Generation1
		{
			get
			{
				return this.m_set_b;
			}
		}

		public void NextGeneration()
		{
			this.m_set_a.Clear();
			this.m_set_a.AddRange(this.m_set_b);
			this.m_set_b.Clear();
		}
	}
}
