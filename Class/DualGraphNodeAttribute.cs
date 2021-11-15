using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace VolumeGeneratorBasedOnGraph.Class
{
    class DualGraphNodeAttribute : ICloneable
    {
        public DualGraphNodeAttribute()
        {

        }

        public DualGraphNodeAttribute(DualGraphNodeAttribute nodeAttribute)
        {

        }

        public object Clone()
        {
            DualGraphNodeAttribute nodeAttribute = (DualGraphNodeAttribute)MemberwiseClone();
            return nodeAttribute;
        }
    }
}
