using System;
using System.Collections.Generic;
using System.Linq;

namespace Geometry
{
    /// <summary>
    /// Class for a branching-junction of the contour
    /// </summary>
    public class ContourJunction
    {
        /// <summary>
        /// index of this contour junction
        /// </summary>
        public int index { get; set; }

        /// <summary>
        /// position of this contour junction
        /// </summary>
        public Position position { get; set; }

        /// <summary>
        /// list of neighbor contour junctions (have a contour segment in between)
        /// </summary>
        public List<ContourJunction> neighbors { get; set; }
        /// <summary>
        /// list of adjacent contour segments (same order as neighbors contour junctions)
        /// </summary>
        public List<ContourSegment> adjacentContourSegments { get; set; }
        /// <summary>
        /// list of adjacent regions
        /// </summary>
        public List<Region> adjacentRegions { get; set; }

        /// <summary>
        /// radius, within no other meshpoint is allowed to be placed
        /// </summary>
        public double radiusForbiddenArea { get; set; } = -1;

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor, which sets all values
        /// </summary>
        public ContourJunction(int index, Position position)
        {
            this.index = index;
            this.position = position;
            neighbors = new List<ContourJunction>();
            adjacentContourSegments = new List<ContourSegment>();
            adjacentRegions = new List<Region>();
        }

        // Set the list of neighboring contour junctions ████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Sets the list of neighboring contour junctions
        /// </summary>
        public void SetNeighborJunctionsAndSegments(List<ContourJunction> neighborContourJunctions, List<ContourSegment> adjacentContourSegments)
        {
            this.neighbors = neighborContourJunctions;
            this.adjacentContourSegments = adjacentContourSegments;
        }

        // Calculates the plumb point positions █████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculates the plumb point positions and puts them on the adjacent segments
        /// </summary>
        public void SetDefaultPositionTriplets(double dPointGlobalJunctions, double dSegmentGlobal)
        {
            for (int i = 0; i < neighbors.Count; i++)
                adjacentContourSegments[i].AddPlumbPosition(position.GetPositionMovedDistanceToOtherPosition(
                    neighbors[i].position, dPointGlobalJunctions), dSegmentGlobal);
        }

        // Sets the forbidden area ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Sets the forbidden area, where no other meshpoint can be placed to
        /// </summary>
        public void SetForbiddenArea<T, R>(Meshing2D_DelaunayVoronoi<T, R> delaunayVoronoi)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            radiusForbiddenArea = Math.Sqrt(delaunayVoronoi.dPointToJunctionGlobal * delaunayVoronoi.dPointToJunctionGlobal
                + delaunayVoronoi.dPointToSegmentGlobal * delaunayVoronoi.dPointToSegmentGlobal);
        }

        // Checks whether a position is within the forbidden area ███████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Checks whether a position is within the forbidden area
        /// </summary>
        public bool IsInForbiddenArea(Position positionToCheck)
        {
            return position.DistanceTo(positionToCheck) < radiusForbiddenArea;
        }

        // Write to console █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Write contour junction to Console
        /// </summary>
        public void Print()
        {
            Console.Write("Contour junction " + index + " at ");
            position.Print();
        }
    }
}