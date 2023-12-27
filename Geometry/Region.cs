using TransferMatrix;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BasicLib;


namespace Geometry
{
    /// <summary>
    /// Class, which represents a certain area with special preferences
    /// </summary>
    public class Region
    {
        /// <summary>
        /// index of the region (telling which region override the other) (lowest index means highest hierarchy order)
        /// </summary>
        public int index { get; set; }

        /// <summary>
        /// list of orderd corner points
        /// </summary>
        public List<ContourJunction> orderedPoints { get; set; }
        /// <summary>
        /// list of linesegments, which connet the corner points (same order as points)
        /// </summary>
        public List<ContourSegment> orderedSegments { get; set; }

        /// <summary>
        /// type of the region
        /// </summary>
        public pointType type { get; set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of a region (needs Initialize afterwards)
        /// </summary>
        public Region()
        {
        }

        // Initialize Region ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor, which sets the corner points and the type. Moreover it is check, if the region has intersections of borders.
        /// </summary>
        /// <param name="index">index of region</param>
        /// <param name="orderedSegments">contour segments ordered clock or counter clock wise</param>
        /// <param name="regionType">type of region (semiconductor, cell)</param>
        public void Initialize(int index, List<ContourSegment> orderedSegments, pointType regionType)
        {
            this.index = index;
            this.orderedSegments = orderedSegments;
            this.type = regionType;

            orderedPoints = new List<ContourJunction>();
            orderedPoints.Add(orderedSegments[0].firstAdjacentContourJunction);
            orderedPoints.Add(orderedSegments[0].secondAdjacentContourJunction);

            if (orderedSegments.Count > 1)
                if (orderedSegments[1].firstAdjacentContourJunction.position.SameAs(orderedPoints[0].position, 0) || orderedSegments[1].secondAdjacentContourJunction.position.SameAs(orderedPoints[0].position, 0))
                    orderedPoints.Reverse();

            for (int i = 1; i < orderedSegments.Count; i++)
            {
                if (!orderedPoints.Select(p => p.position).Any(c => c.SameAs(orderedSegments[i].firstAdjacentContourJunction.position, 0)))
                    orderedPoints.Add(orderedSegments[i].firstAdjacentContourJunction);
                if (!orderedPoints.Select(p => p.position).Any(c => c.SameAs(orderedSegments[i].secondAdjacentContourJunction.position, 0)))
                    orderedPoints.Add(orderedSegments[i].secondAdjacentContourJunction);
            }
        }

        // set all specific properties of a region ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// set all specific properties of a region
        /// </summary>
        /// <returns></returns>
        public virtual void SetProperties(double[] preferencesArray, pointType regionType)
        {
            Console.WriteLine("Override the function \"" + System.Reflection.MethodBase.GetCurrentMethod().Name
                + "()\" in the class \"" + this.GetType().Name + "\"!");
        }
        public virtual void SetOpticalModel((double roughnessFrontGrid, double roughnessFrontContact, double roughnessAbsorber, double roughnessBackContact, double roughnessBackGrid, int IDbefore, (int ID, double roughness) behind, List<(int ID, double thickness)> incoherent,
            List<(int ID, double thickness, double roughness)> aboveFrontGrid, List<(int ID, double thickness, double roughness)> aboveAbsorber,
            List<(int ID, double thickness, double roughness)> belowAbsorber, List<(int ID, double thickness, double roughness)> belowBackGrid) opticalModel)
        {
            Console.WriteLine("Override the function \"" + System.Reflection.MethodBase.GetCurrentMethod().Name
                + "()\" in the class \"" + this.GetType().Name + "\"!");
        }

        // Calculate area of the region ████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the 2D-area of the region
        /// </summary>
        public double GetArea()
        {
            // Shoelace formula (Gaussian trapezoidal formula)
            // https://en.wikipedia.org/wiki/Shoelace_formula
            List<Position> cornerPoints = orderedPoints.Select(p => p.position).ToList();
            double sum = 0;
            for (int k = 0; k < cornerPoints.Count; k++)
                sum += (cornerPoints[k].x + cornerPoints[(k + 1) % cornerPoints.Count].x) * (cornerPoints[(k + 1)
                    % cornerPoints.Count].y - cornerPoints[k].y);
            return sum / 2;
        }

        // print console ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// prints region to console
        /// </summary>
        public void Print()
        {
            Console.WriteLine("Region " + type.ToString());
            Console.WriteLine("corner points:");
            foreach (Position pos in orderedPoints.Select(p => p.position).ToList())
            {
                Console.Write("\t");
                pos.Print();
            }
        }
    }
}