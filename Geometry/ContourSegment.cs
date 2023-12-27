using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;

namespace Geometry
{
    /// <summary>
    /// Class for a contour segment between two contour junctions
    /// </summary>
    public class ContourSegment
    {
        /// <summary>
        /// index of this contour segment
        /// </summary>
        public int index { get; set; }

        /// <summary>
        /// size of this segment, i.e. its length
        /// </summary>
        public double Size { get; private set; }

        /// <summary>
        /// frist adjacent contour junction of this contour segment
        /// </summary>
        public ContourJunction firstAdjacentContourJunction { get; private set; }
        /// <summary>
        /// second adjacent contour junction of this contour segment
        /// </summary>
        public ContourJunction secondAdjacentContourJunction { get; private set; }
        /// <summary>
        /// list of adjacent regions
        /// </summary>
        public List<Region> adjacentRegions { get; set; }

        /// <summary>
        /// line segment of the contour segment from the first junction to the second junction
        /// </summary>
        public LineSegment lineSegment { get; private set; }

        /// <summary>
        /// default distance of the plumb position to each other (along the segment)
        /// </summary>
        public double dPoint { get; private set; }

        /// <summary>
        /// list of position triplets including plumbpoints (auxiliary positions which are directly on the contour) and mesh points (will be placed orthogonally to the plumb position with distnace dSegment) (ordered from frist contour junction to second one)
        /// </summary>
        public List<(Position plumbPos, Position firstMeshPos, Position secondMeshPos)> pointTriplets { get; set; }
            = new List<(Position plumbPos, Position firstMeshPos, Position secondMeshPos)>();

        /// <summary>
        /// list with plumb points and construction rule (position + plumbPointShift). polygonIndex: polygon where the shifted plumb point will be in, plumbPointShift: direction of plumb point shift, positions: list with plumb points on segment
        /// </summary>
        public List<(int polygonIndex, double[] plumbPointShift, List<Position> positions)> PlumbPoints = new List<(int, double[], List<Position>)>();

        /// <summary>
        /// minimal distance allowed for a meshpoint to this segment
        /// </summary>
        public double MinimalMeshPointDistanceToSegmentAllowed = 0;

        /// <summary>
        /// corner positions of the area, where no other mesh point is allowed to put to (otherwise, the contour of the segment would be destroyed)
        /// </summary>
        public Position[] forbiddenArea { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor, which sets all values
        /// </summary>
        public ContourSegment(int index, ContourJunction firstAdjacentContourJunction, ContourJunction secondAdjacentContourJunction)
        {
            this.index = index;
            this.firstAdjacentContourJunction = firstAdjacentContourJunction;
            this.secondAdjacentContourJunction = secondAdjacentContourJunction;
            lineSegment = new LineSegment(firstAdjacentContourJunction.position, secondAdjacentContourJunction.position);
            adjacentRegions = new List<Region>();

            // length of this segment
            Size = firstAdjacentContourJunction.position.DistanceTo(secondAdjacentContourJunction.position);
        }

        // Adds a single Plumbposition and directly calculates its corresponding mesh positions █████████████████████████████████████████████████████
        /// <summary>
        /// Adds a single Plumbposition and directly calculates its corresponding mesh positions
        /// </summary>
        public void AddPlumbPosition(Position plumbPosition, double dSegmentGlobal)
        {
            Position meshpoint1 = plumbPosition.GetPositionMovedDistancePerpendicularToOtherPosition(firstAdjacentContourJunction.position,
                    dSegmentGlobal, false);

            Position meshpoint2 = plumbPosition.GetPositionMovedDistancePerpendicularToOtherPosition(firstAdjacentContourJunction.position,
                    dSegmentGlobal, true);

            pointTriplets.Add((plumbPosition, meshpoint1, meshpoint2));
        }

        // Calculates the plumb point positions █████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculates the plumb point positions
        /// </summary>
        public void SetDefaultPositionTriplets(double dPointGlobal, double dPointGlobalJunctions, double dSegmentGlobal)
        {
            Position plumbPositionJunction1 = firstAdjacentContourJunction.position.GetPositionMovedDistanceToOtherPosition(
                secondAdjacentContourJunction.position, dPointGlobalJunctions);
            Position plumbPositionJunction2 = secondAdjacentContourJunction.position.GetPositionMovedDistanceToOtherPosition(
                firstAdjacentContourJunction.position, dPointGlobalJunctions);

            int amountOfDivisions = Convert.ToInt32(Math.Ceiling(plumbPositionJunction1.DistanceTo(plumbPositionJunction2)
                / dPointGlobal));
            int amountOfPoints = amountOfDivisions - 1;
            Position vectorFrom1to2 = plumbPositionJunction2 - plumbPositionJunction1;
            for (int i = 0; i < amountOfPoints; i++)
                AddPlumbPosition(plumbPositionJunction1 + Convert.ToDouble(i + 1) / Convert.ToDouble(amountOfDivisions) * vectorFrom1to2,
                    dSegmentGlobal);

            this.dPoint = plumbPositionJunction1.DistanceTo(plumbPositionJunction2) / Convert.ToDouble(amountOfDivisions);
        }

        // Orders all point triplets ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Orders all point triplets from the first adjacent contour junction to the second one
        /// </summary>
        public void OrderPositionTriplets()
        {
            pointTriplets = pointTriplets.OrderBy(t => t.plumbPos.DistanceTo(firstAdjacentContourJunction.position)).ToList();
        }

        // Sets the forbidden area ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Sets the forbidden area, where no other meshpoint can be placed to
        /// </summary>
        public void SetForbiddenArea<T, R>(Meshing2D_DelaunayVoronoi<T, R> delaunayVoronoi)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            Position plumbPositionJunction1 = firstAdjacentContourJunction.position.GetPositionMovedDistanceToOtherPosition(
                secondAdjacentContourJunction.position, delaunayVoronoi.dPointToJunctionGlobal);
            Position plumbPositionJunction2 = secondAdjacentContourJunction.position.GetPositionMovedDistanceToOtherPosition(
                firstAdjacentContourJunction.position, delaunayVoronoi.dPointToJunctionGlobal);

            double minDistance = 0.5 * Math.Sqrt(dPoint * dPoint + 4 * delaunayVoronoi.dPointToSegmentGlobal * delaunayVoronoi.dPointToSegmentGlobal);
            forbiddenArea = new Position[4];
            forbiddenArea[0] = plumbPositionJunction1.GetPositionMovedDistancePerpendicularToOtherPosition(plumbPositionJunction2,
                minDistance, true);
            forbiddenArea[1] = plumbPositionJunction1.GetPositionMovedDistancePerpendicularToOtherPosition(plumbPositionJunction2,
                minDistance, false);
            forbiddenArea[2] = plumbPositionJunction2.GetPositionMovedDistancePerpendicularToOtherPosition(plumbPositionJunction1,
                minDistance, true);
            forbiddenArea[3] = plumbPositionJunction2.GetPositionMovedDistancePerpendicularToOtherPosition(plumbPositionJunction1,
                minDistance, false);
        }

        // Checks whether a position is within the forbidden area ███████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Checks whether a position is within the forbidden area
        /// </summary>
        public bool IsInForbiddenArea(Position positionToCheck)
        {
            return positionToCheck.InPolygon(forbiddenArea.ToList(), true);
        }

        // Write to console █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Write contour segment to Console
        /// </summary>
        public void Print()
        {
            Console.Write("Contour segment " + index + " from ");
            firstAdjacentContourJunction.position.Print(false);
            Console.Write(" to ");
            secondAdjacentContourJunction.position.Print();
        }
    }
}