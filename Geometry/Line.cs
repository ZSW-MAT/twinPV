using BasicLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Geometry
{
    /// <summary>
    /// Class, which represents a line
    /// </summary>
    public class Line
    {
        /// <summary>
        /// position vector
        /// </summary>
        public double[] positionVector { get; private set; }
        /// <summary>
        /// direction vector
        /// </summary>
        public double[] directionVector { get; private set; }

        // Constructor from 2 vectors ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor from 2 vectors
        /// </summary>
        /// <param name="positionVector">position Vector of the new line</param>
        /// <param name="directionVector">direction Vector of the new line</param>
        public Line(double[] positionVector, double[] directionVector)
        {
            this.positionVector = positionVector;
            this.directionVector = directionVector;
        }

        // Constructor from 1 Position and 2 saclars ████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor from 1 Position and 2 saclars
        /// </summary>
        /// <param name="positionVector">position Vector of the new line</param>
        /// <param name="directionVectorX">x-component of the direction vector</param>
        /// <param name="directionVectorY">y-component of the direction vector</param>
        public Line(Position positionVector, double directionVectorX, double directionVectorY)
            : this(new double[] { positionVector.x, positionVector.y }, new double[] { directionVectorX, directionVectorY })
        { }

        // Constructor from 2 positions █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor from 2 positions
        /// </summary>
        /// <param name="startingPosition">starting positon</param>
        /// <param name="endPosition">endposition</param>
        public Line(Position startingPosition, Position endPosition)
            : this(new double[] { startingPosition.x, startingPosition.y, startingPosition.z },
                  new double[] { endPosition.x - startingPosition.x, endPosition.y - startingPosition.y, endPosition.z - startingPosition.z })
        { }

        // mirrow point on line █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Return the mirrored position on this line
        /// </summary>
        /// <param name="pos">position which will be mirrored</param>
        /// <returns></returns>
        public Position Mirrorpoint(Position pos)
        {
            // get perpendicular position
            Position perpendicularPosition = pos.PerpendicularPoint(this);

            return new Position(2 * perpendicularPosition.x - pos.x, 2 * perpendicularPosition.y - pos.y);
        }

        // intersection with other line █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the position which is the intersection of both lines
        /// </summary>
        /// <param name="secondLine">second line, where the point should be on</param>
        /// <returns></returns>
        public Position Intersection(Line secondLine)
        {
            // set both lines in one equation (vectorial), calculate line parameter of 2nd line, put into first line
            double numerator = directionVector[0] * secondLine.positionVector[1] - directionVector[0] * positionVector[1]
                + directionVector[1] * positionVector[0] - directionVector[1] * secondLine.positionVector[0];
            double denominator = directionVector[1] * secondLine.directionVector[0] - directionVector[0] * secondLine.directionVector[1];
            double lineParameter = numerator / denominator;

            double x = secondLine.positionVector[0] + lineParameter * secondLine.directionVector[0];
            double y = secondLine.positionVector[1] + lineParameter * secondLine.directionVector[1];

            return new Position(x, y);
        }

        // Checks, wheather a second line is parallel to this one ███████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Checks, wheather a second line is parallel to this one
        /// </summary>
        /// <param name="secondLine">second line, which will be checked</param>
        /// <returns></returns>
        public bool IsParallelTo(Line secondLine)
        {
            double scalarProduct = directionVector[0] * secondLine.directionVector[0] + directionVector[1] * secondLine.directionVector[1];
            double squareOfAbsoluteValues = Misc.GetPNormOfVector(directionVector, 2) * Misc.GetPNormOfVector(secondLine.directionVector, 2);

            if (squareOfAbsoluteValues == 0)
                return false;

            double division = Math.Abs(scalarProduct / squareOfAbsoluteValues) - 1;

            return Math.Abs(division) < 1e-12;
        }

        // Checks, wheather a second line is perpendicular to this one ██████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Checks, wheather a second line is perpendicular to this one
        /// </summary>
        /// <param name="secondLine">second line, which will be checked</param>
        /// <returns></returns>
        public bool IsPerpendicularTo(Line secondLine)
        {
            double scalarProduct = directionVector[0] * secondLine.directionVector[0] + directionVector[1] * secondLine.directionVector[1];
            double squareOfAbsoluteValues = Misc.GetPNormOfVector(directionVector, 2) * Misc.GetPNormOfVector(secondLine.directionVector, 2);

            if (squareOfAbsoluteValues == 0)
                return false;

            double division = Math.Abs(scalarProduct / squareOfAbsoluteValues);

            return Math.Abs(division) < 1e-12;
        }

        // Returns distance to parallel line ████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns distance to parallel line
        /// </summary>
        /// <param name="parallelLine">parallel line, where the distance will be calculate to (NEEDS TO BE PARALLEL TO THE FRIST LINE! Otherwise returns NaN)</param>
        /// <returns></returns>
        public double DistanceToParallelLine(Line parallelLine)
        {
            if (!IsParallelTo(parallelLine))
                return double.NaN;

            return new Position(positionVector).DistanceTo(parallelLine);
        }

        // Returns the bisectric line ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns the bisectric line (the one, which has the smaller angle to both input lines)
        /// </summary>
        /// <param name="secondLine">line, to which the bisectrix is constructed</param>
        /// <returns></returns>
        public Line GetBisectrixWithSmallerAngle(Line secondLine)
        {
            Position intersection = Intersection(secondLine);

            // Get normalized direction vectors
            Position direction1 = new Position(directionVector.Select(c => c / Misc.GetPNormOfVector(directionVector, 2)).ToArray());
            Position direction2 = new Position(secondLine.directionVector.Select(c => c / Misc.GetPNormOfVector(secondLine.directionVector, 2)).ToArray());

            // Get two possible second positions (one is for the bisectrix with the larger, one for the smaller angle)
            Position secondpoint_a = intersection + direction1 + direction2;
            Position secondpoint_b = intersection + direction1 - direction2;

            // pic the point for the line with the smaller angle => point is farer away from intersection
            Position secondpoint = new List<Position> { secondpoint_a, secondpoint_b }.
                OrderBy(p => p.DistanceTo(intersection)).Last();

            return new Line(intersection, secondpoint);
        }

        // Print point in console ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// write line to console
        /// </summary>
        public void Print()
        {
            Console.WriteLine("x = (" + positionVector[0] + "|" + positionVector[1] + ") + p * ("
                + directionVector[0] + "|" + directionVector[1] + ")");
        }
    }

    /// <summary>
    /// Class, which represents a linesegment (subcalss of linesegment)
    /// </summary>
    public class LineSegment : Line
    {
        /// <summary>
        /// starting Position of linesegment
        /// </summary>
        public Position startingPoint { get; private set; }
        /// <summary>
        /// end position of linesegment
        /// </summary>
        public Position endPoint { get; private set; }

        // Constuctor from 2 positions ██████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor from 2 positions
        /// </summary>
        /// <param name="startingPoint">starting point of the line segment</param>
        /// <param name="endPoint">end point of the line segment</param>
        public LineSegment(Position startingPoint, Position endPoint) : base(startingPoint, endPoint)
        {
            this.startingPoint = startingPoint;
            this.endPoint = endPoint;
        }
        /// <summary>
        /// Constructor from geometrySegment
        /// </summary>
        public LineSegment(ContourSegment segment) : base(segment.firstAdjacentContourJunction.position, segment.secondAdjacentContourJunction.position)
        {
            startingPoint = segment.firstAdjacentContourJunction.position;
            endPoint = segment.secondAdjacentContourJunction.position;
        }

        // Return the start and end position of the common parallel section █████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns the start and end position of the common parallel section on the first line segment (returns ((NaN, NaN), (NaN, NaN)) if there is no common section)
        /// </summary>
        /// <param name="parallelLinesegment">parallel Linesegment, which will be checked for a parallel section with the first one</param>
        /// <returns></returns>
        public (Position positionFrom, Position positionTo) GetSectionOfParallelLinesegments(LineSegment parallelLinesegment)
        {
            // all involved positions
            Position secondStartingPoint = parallelLinesegment.startingPoint.PerpendicularPoint(this);
            Position secondEndPoint = parallelLinesegment.endPoint.PerpendicularPoint(this);
            (Position position, int origin)[] allPositions =
                new (Position position, int origin)[] { (startingPoint, 1), (endPoint, 1), (secondStartingPoint, -1), (secondEndPoint, -1) };

            // reference, where all lengths will be meassured to -> one of the outer points -> point with the largest distance to any other point
            (Position position, int index) refPos = (new Position(0, 0), 0);
            double currentMaxDistance = 0;

            for (int i = 0; i < allPositions.Length; i++)
                for (int j = i + 1; j < allPositions.Length; j++)
                {
                    double distance = allPositions[i].position.DistanceTo(allPositions[j].position);
                    if (distance > currentMaxDistance)
                    {
                        currentMaxDistance = distance;
                        refPos = allPositions[i];
                    }
                }
            (Position position, int origin, double maximumDistance) referencePosition = (refPos.position, refPos.index, currentMaxDistance);

            // order all positions according to distance to referecn point (might be that same distance appreas twice)
            List<(Position position, int origin)> orderedAllPositions =
                allPositions.OrderBy(p => p.position.DistanceTo(referencePosition.position)).ToList();

            // same list of same lenght as top one, but if one position is in twice, then both their origin indexes will be 0 (instead of -1 or 1)
            List<(Position position, int origin)> orderedPositions = new List<(Position position, int origin)>();
            orderedPositions.Add(orderedAllPositions[0]);
            for (int i = 1; i < 4; i++)
            {
                if (orderedPositions.Last().position.SameAs(orderedAllPositions[i].position, 1e-14 * referencePosition.maximumDistance))
                {
                    orderedPositions.RemoveAt(orderedPositions.Count - 1);
                    orderedPositions.Add((orderedAllPositions[i].position, 0));
                    orderedPositions.Add((orderedAllPositions[i].position, 0));
                }
                else
                    orderedPositions.Add(orderedAllPositions[i]);
            }

            // determine in which order the positions appear and return borders, if they have a common parallel section
            switch (orderedPositions[0].origin * orderedPositions[1].origin)
            {
                case 1: // first 2 points are from same segment (++--) or (--++) => segments do not have common parallel section
                    break;

                case -1: // first 2 points are from different segments (+-  ) or (-+  ) => common parallel section: second point to third point (even if last 2 origins are 0)
                    return (orderedPositions[1].position, orderedPositions[2].position);

                case 0: // 3 possibilities: 1. (0+  )/(0-  ) cannot happen! 2. (+00-)/(-00+) no common section, only endpoints overlay 3. (00+-)/(00-+) or (0000) common section goes from second to third
                    if (orderedPositions[0].origin != 0)
                        break;
                    else
                        return (orderedPositions[1].position, orderedPositions[2].position);
            }

            // segments do not have common parallel section
            return (new Position(double.NaN, double.NaN), new Position(double.NaN, double.NaN));
        }

        // Check, if 2 linesegments do intersect ████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Check, if 2 linesegments do intersect
        /// </summary>
        /// <param name="secondLinesegment">second Linesegment, which will be checked for an intersection</param>
        /// <param name="includeEdgesOfFirstSegment">do endpoints of first segment also count as intersection?</param>
        /// <param name="includeEdgesOfFirstSegment">do endpoints of second segment also count as intersection?</param>
        /// <param name="firstSegementIntersects">output parameter telling that the endpoints of the first segment intersects</param>
        /// <param name="secondSegementIntersects">output parameter telling that the endpoints of the second segment intersects</param>
        /// <returns></returns>
        public bool DoesIntersect(LineSegment secondLinesegment, bool includeEdgesOfFirstSegment, bool includeEdgesOfSecondSegment,
            out bool firstSegementIntersects, out bool secondSegementIntersects)
        {
            // https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

            // Find the four orientations needed for general and special cases 
            int o1 = Orientation(startingPoint, endPoint, secondLinesegment.startingPoint);
            int o2 = Orientation(startingPoint, endPoint, secondLinesegment.endPoint);
            int o3 = Orientation(secondLinesegment.startingPoint, secondLinesegment.endPoint, startingPoint);
            int o4 = Orientation(secondLinesegment.startingPoint, secondLinesegment.endPoint, endPoint);

            firstSegementIntersects = false;
            secondSegementIntersects = false;

            // if only edges are on segment
            if (!includeEdgesOfFirstSegment)
                if (o3 == 0 || o4 == 0)
                    return false;
            if (!includeEdgesOfSecondSegment)
                if (o1 == 0 || o2 == 0)
                    return false;

            // Special Cases
            // if the 3 points a colinear and the center point lies in between the outside points
            if (o1 == 0 && OnSegment(startingPoint, secondLinesegment.startingPoint, endPoint))
            {
                secondSegementIntersects = true;
                return true;
            }
            if (o2 == 0 && OnSegment(startingPoint, secondLinesegment.endPoint, endPoint))
            {
                secondSegementIntersects = true;
                return true;
            }
            if (o3 == 0 && OnSegment(secondLinesegment.startingPoint, startingPoint, secondLinesegment.endPoint))
            {
                firstSegementIntersects = true;
                return true;
            }
            if (o4 == 0 && OnSegment(secondLinesegment.startingPoint, endPoint, secondLinesegment.endPoint))
            {
                firstSegementIntersects = true;
                return true;
            }

            // General case 
            if (o1 != o2 && o3 != o4)
                return true;

            // Doesn't fall in any of the above cases
            return false;
        }
        /// <summary>
        /// Check, if 2 linesegments do intersect
        /// </summary>
        /// <param name="secondLinesegment">second Linesegment, which will be checked for an intersection</param>
        /// <param name="includeEdges">do endpoints also count as intersection?</param>
        public bool DoesIntersect(LineSegment secondLinesegment, bool includeEdges = true)
        {
            return DoesIntersect(secondLinesegment, includeEdges, includeEdges, out bool firstSegementIntersects, out bool secondSegementIntersects);
        }
        /// <summary>
        /// checks, if Position q lies between Positions p and r, if all 3 are colinear to each other
        /// </summary>
        /// <param name="p">left point</param>
        /// <param name="q">middle point</param>
        /// <param name="r">right point</param>
        /// <returns></returns>
        bool OnSegment(Position p, Position q, Position r)
        {
            // Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr' 
            if (q.x <= Math.Max(p.x, r.x) && q.x >= Math.Min(p.x, r.x) &&
                q.y <= Math.Max(p.y, r.y) && q.y >= Math.Min(p.y, r.y))
                return true;

            return false;
        }
        /// <summary>
        /// returns the local orientation of the 3 points ("0": colinear, "1": clockwise, "2": counterclockwise)
        /// </summary>
        /// <param name="p">frist point</param>
        /// <param name="q">second point</param>
        /// <param name="r">third point</param>
        /// <returns></returns>
        int Orientation(Position p, Position q, Position r)
        {
            // To find orientation of ordered triplet (p, q, r). The function returns following values:
            // 0 --> p, q and r are colinear 
            // 1 --> Clockwise 
            // 2 --> Counterclockwise 
            // https://www.geeksforgeeks.org/orientation-3-ordered-points/

            double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

            if (val == 0) return 0; // colinear

            return (val > 0) ? 1 : 2; // clock or counterclock wise 
        }

        // distance to line segment █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return sthe euclidian distance to another finite long line segment (either distance to infinite long line or to ond endpoint)
        /// </summary>
        /// <param name="linesegment">other linesegment, to which the distance will be measured</param>
        /// <returns></returns>
        public double DistanceToLinesegment(LineSegment linesegment)
        {
            // calculate all distances of the four end-points to the other linesegment and return minimum
            double distance1 = startingPoint.DistanceTo(linesegment);
            double distance2 = endPoint.DistanceTo(linesegment);
            double distance3 = linesegment.startingPoint.DistanceTo(this);
            double distance4 = linesegment.endPoint.DistanceTo(this);

            return Math.Min(Math.Min(distance1, distance2), Math.Min(distance3, distance4));
        }
    }
}