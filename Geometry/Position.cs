using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using BasicLib;

namespace Geometry
{
    /// <summary>
    /// Class, that contains a tuple of x and y values and is able to perform several geometrical operations
    /// </summary>
    [DataContract]
    public class Position
    {
        /// <summary>
        /// x-coordinate
        /// </summary>
        [DataMember]
        public double x { get; set; }
        /// <summary>
        /// y-coordinate
        /// </summary>
        [DataMember]
        public double y { get; set; }
        /// <summary>
        /// z-coordinate
        /// </summary>
        [DataMember]
        public double z { get; set; }

        // Constructor from to double numbers ███████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// parameterless constructor
        /// </summary>
        public Position()
        {
        }
        /// <summary>
        /// Constructor (1d) from to double numbers
        /// </summary>
        /// <param name="x">x-coordinate of the position</param>
        public Position(double x)
        {
            this.x = x;
            this.y = 0;
            this.z = 0;
        }
        /// <summary>
        /// Constructor (2d) from to double numbers
        /// </summary>
        /// <param name="x">x-coordinate of the position</param>
        /// <param name="y">y-coordinate of the position</param>
        public Position(double x, double y)
        {
            this.x = x;
            this.y = y;
            this.z = 0;
        }
        /// <summary>
        /// Constructor (3d) from double numbers
        /// </summary>
        /// <param name="x">x-coordinate of the position</param>
        /// <param name="y">y-coordinate of the position</param>
        /// <param name="z">z-coordinate of the position</param>
        public Position(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        // Constructor from vector ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor from vector (2d)
        /// </summary>
        /// <param name="positionVector">position vector of the position</param>
        public Position(double[] positionVector)
        {
            x = positionVector[0];
            y = positionVector[1];

            if (positionVector.Length > 2) z = positionVector[2];
        }

        // aux/help/support ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns x, y, z as array
        /// </summary>
        /// <returns></returns>
        public double[] PositionToArray()
        {
            return new double[] { x, y, z };
        }

        // add two positions ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// adds two positions component-wise: x_new = x1 + x2, y_new = y1 + y2, and z_new = z1 + z2
        /// </summary>
        /// <param name="firstPosition">frist position</param>
        /// <param name="secondPosition">second position</param>
        /// <returns></returns>
        public static Position operator +(Position firstPosition, Position secondPosition)
        {
            return new Position(firstPosition.x + secondPosition.x, firstPosition.y + secondPosition.y, firstPosition.z + secondPosition.z);
        }

        // subtract two positions ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// subtracts two positions component-wise: x_new = x1 + x2, y_new = y1 + y2, and z_new = z1 + z2
        /// </summary>
        /// <param name="firstPosition">positive weighted position</param>
        /// <param name="secondPosition">negative weighted position</param>
        /// <returns></returns>
        public static Position operator -(Position firstPosition, Position secondPosition)
        {
            return new Position(firstPosition.x - secondPosition.x, firstPosition.y - secondPosition.y, firstPosition.z - secondPosition.z);
        }

        // multiply position component-wise with a scalar factor ████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// multiplies the position component-wise with a scalar factor 
        /// </summary>
        /// <param name="position">Position</param>
        /// <param name="multiplicator">scalar, whith which the position will be multiplied component-wise</param>
        /// <returns></returns>
        public static Position operator *(Position position, double multiplicator)
        {
            return new Position(position.x * multiplicator, position.y * multiplicator, position.z * multiplicator);
        }
        /// <summary>
        /// multiplies the position component-wise with a scalar factor 
        /// </summary>
        /// <param name="multiplicator">scalar, whith which the position will be multiplied component-wise</param>
        /// <param name="position">Position</param>
        /// <returns></returns>
        public static Position operator *(double multiplicator, Position position)
        {
            return position * multiplicator;
        }

        // divide position component-wise by a scalar ███████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// divides position component-wise by a scalar
        /// </summary>
        /// <param name="pos">Position</param>
        /// <param name="divisor">scalar, which the position componentes will be divided by</param>
        /// <returns></returns>
        public static Position operator /(Position position, double divisor)
        {
            return new Position(position.x / divisor, position.y / divisor, position.z / divisor);
        }

        // checks, if two positions are the same ████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if two positions are the same (with tolerance)
        /// </summary>
        /// <param name="secondPosition">positon, which will be checked to be the same</param>
        /// <param name="tolerance">tolerance, which will be accepted to be the same positon</param>
        /// <returns></returns>
        public bool SameAs(Position secondPosition, double tolerance)
        {
            return DistanceTo(secondPosition) <= tolerance;
        }

        // checks, if two positions are the same ████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if two vectors (positions) show in the same direction (1) or opposite to each other (-1) or anywhere else (according to scalar product)
        /// </summary>
        /// <param name="secondPosition">positon, which will be checked for its direction</param>
        /// <returns></returns>
        public double CheckVectorDirection(Position secondPosition)
        {
            double scalar = x * secondPosition.x + y * secondPosition.y + z * secondPosition.z;
            double norm = Math.Sqrt((x * x + y * y + z * z) * (secondPosition.x * secondPosition.x + secondPosition.y * secondPosition.y + secondPosition.z * secondPosition.z));
            double divided = scalar / norm;

            // round for small deviations
            if (Math.Abs(divided - 1) < 1e-12)
                divided = 1;
            if (Math.Abs(divided + 1) < 1e-12)
                divided = -1;

            return divided;
        }

        // center position with other position ██████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return the position in the middle of both positions
        /// </summary>
        /// <param name="secondPosition">second position</param>
        /// <returns></returns>
        public Position CenterWith(Position secondPosition)
        {
            return new Position((x + secondPosition.x) / 2, (y + secondPosition.y) / 2, (z + secondPosition.z) / 2);
        }

        // point on ratio to the second point ███████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Return the point on the line beween the two points (also outside the two positions possible)
        /// </summary>
        /// <param name="secondPosition">second position</param>
        /// <param name="ratioTofirstPoint">ratio, where the new point will be
        /// /// (ratio = 0 -> point1, ratio = 0.5 -> center, ratio = 1 -> point2)</param>
        /// <returns></returns>
        public Position PercentalCenterWith(Position secondPosition, double ratioTofirstPoint)
        {
            double x_new = x + (secondPosition.x - x) * ratioTofirstPoint;
            double y_new = y + (secondPosition.y - y) * ratioTofirstPoint;
            double z_new = z + (secondPosition.z - z) * ratioTofirstPoint;
            return new Position(x_new, y_new, z_new);
        }

        // distance to other position ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Return the squared of the euclidian distance between two positions
        /// </summary>
        /// <param name="secondPosition">second position</param>
        /// <returns></returns>
        public double DistanceSquaredTo(Position secondPosition)
        {
            double deltaX = x - secondPosition.x;
            double deltaY = y - secondPosition.y;
            double deltaZ = z - secondPosition.z;
            return deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
        }
        /// <summary>
        /// Return the euclidian distance between two positions
        /// </summary>
        /// <param name="secondPosition">second position</param>
        /// <returns></returns>
        public double DistanceTo(Position secondPosition)
        {
            double deltaX = x - secondPosition.x;
            double deltaY = y - secondPosition.y;
            double deltaZ = z - secondPosition.z;
            return Math.Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        }

        /// <summary>
        /// Return the euclidian distance between two positions
        /// </summary>
        /// <param name="secondPosition">second position (must be in 3d!)</param>
        /// <returns>distance</returns>
        public double DistanceTo(double[] secondPosition)
        {
            double deltaX = x - secondPosition[0];
            double deltaY = y - secondPosition[1];
            double deltaZ = z - secondPosition[2];
            return Math.Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        }

        /// <summary>
        /// Returns the distance to the next region
        /// </summary>
        /// <param name="regions">list of regions, whcere the nearest distance to, is returned</param>
        /// <returns></returns>
        public double DistanceToNextRegion(List<Region> regions)
        {
            double minimumDistance = 1e50;
            foreach (Region region in regions)
                if (DistanceToRegion(region) < minimumDistance)
                    minimumDistance = DistanceToRegion(region);
            return minimumDistance;
        }
        /// <summary>
        /// Returns the distance to a certain region
        /// </summary>
        /// <param name="region">region, whereto the distance is calculated</param>
        /// <returns></returns>
        public double DistanceToRegion(Region region)
        {
            // return 0 if point in in region
            if (InPolygon(region, false))
                return 0;

            // otherwise: calculate distance to all edges and choose minimum
            List<Position> cornerPoints = region.orderedPoints.Select(p => p.position).ToList();
            double minimumDistance = double.PositiveInfinity;
            for (int i = 0; i < cornerPoints.Count; i++)
            {
                LineSegment linesegment = new LineSegment(cornerPoints[i], cornerPoints[Misc.mod(i + 1, cornerPoints.Count)]);
                if (DistanceTo(linesegment) < minimumDistance)
                    minimumDistance = DistanceTo(linesegment);
            }
            return minimumDistance;
        }

        /// <summary>
        /// returns the euclidian distance to an (infinite long) line
        /// </summary>
        /// <param name="line">line, to which the distance is measured</param>
        /// <returns></returns>
        public double DistanceTo(Line line)
        {
            // calculate perpendicular point
            Position perpendicularPoint = PerpendicularPoint(line);

            // return distance to perpendicular point
            return DistanceTo(perpendicularPoint);
        }

        // distance to line segment █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return sthe euclidian distance to a finite long line segment (either distance to infinite long line or to ond endpoint)
        /// </summary>
        /// <param name="linesegment">linesegment, to which the distance will be measured</param>
        /// <returns></returns>
        public double DistanceTo(LineSegment linesegment)
        {
            // perpendicular point on infinite long line
            Position perpendicularPosition = PerpendicularPoint(linesegment);

            // create line from perpendicular point to start positoin of linesegment
            Line g1 = new Line(perpendicularPosition, linesegment.startingPoint.x - perpendicularPosition.x,
                linesegment.startingPoint.y - perpendicularPosition.y);

            // calculate line parameter for starting and endpoint (for starting it is always 1)
            double para2;
            // calculation only with one coordinate (larger one for numerical stability), cuz perpendicular point is for sure on line
            if (Math.Abs(g1.directionVector[0]) > Math.Abs(g1.directionVector[1]))
                para2 = (linesegment.endPoint.x - g1.positionVector[0]) / g1.directionVector[0];
            else
                para2 = (linesegment.endPoint.y - g1.positionVector[1]) / g1.directionVector[1];

            // if para1 * para2 = para2 < 0, then perpendicular point is between starting and end point => within line segment
            if (para2 < 0)
                return DistanceTo(perpendicularPosition);
            // if not, then perpendicular point it outside the segment and the distance will be the distance to the nearer endpoint
            else
                return Math.Min(DistanceTo(linesegment.startingPoint), DistanceTo(linesegment.endPoint));
        }

        // get perpendicular point on line ██████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// get perpendicular point on line (intersection of perpendicular with line)
        /// </summary>
        /// <param name="line">line, on which the perpendicular will be</param>
        /// <returns></returns>
        public Position PerpendicularPoint(Line line)
        {
            // create perpendicular line, (position vector of point and direction vector 90° rotated)
            Line perpendicularLine;
            if (line.directionVector[0] < line.directionVector[1])
                perpendicularLine = new Line(this, -line.directionVector[1], line.directionVector[0]);
            else
                perpendicularLine = new Line(this, line.directionVector[1], -line.directionVector[0]);

            // return intersetcion of normal line with perpendicular line
            return line.Intersection(perpendicularLine);
        }

        // Returns a position, which was moved towards another position █████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// takes the position and moves a certain distance towards another point (distance can also be larger than the distance of the two points)
        /// </summary>
        /// <param name="positionToMoveTo">position, in which direction the returned point will be</param>
        /// <param name="distanceToMove">distance how far the returned point will be set to</param>
        /// <returns></returns>
        public Position GetPositionMovedDistanceToOtherPosition(Position positionToMoveTo, double distanceToMove)
        {
            Position vector = positionToMoveTo - this;
            double normalizeddistance = distanceToMove / Math.Sqrt(vector.x * vector.x + vector.y * vector.y);
            return this + normalizeddistance * vector;
        }

        // Returns a position, which was moved perpendicular to another position ████████████████████████████████████████████████████████████████████
        /// <summary>
        /// takes the position and moves a certain distance perpendicular to another point
        /// </summary>
        /// <param name="positionToMovePerpendicularTo">position, in which direction the returned point will be</param>
        /// <param name="distanceToMove">distance how far the returned point will be set to</param>
        /// <returns></returns>
        public Position GetPositionMovedDistancePerpendicularToOtherPosition(Position positionToMovePerpendicularTo, double distanceToMove,
            bool rotateClockwise)
        {
            Position vectorParallel = positionToMovePerpendicularTo - this;

            double multiplicator = 1;
            if (rotateClockwise)
                multiplicator = -1;

            Position vector = new Position(vectorParallel.y * multiplicator, -vectorParallel.x * multiplicator);
            double normalizeddistance = distanceToMove / Math.Sqrt(vector.x * vector.x + vector.y * vector.y);
            return this + normalizeddistance * vector;
        }

        // Vector to another position ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns an MathNET vector to another position
        /// </summary>
        /// <param name="secondPosition">second Position</param>
        /// <returns></returns>
        public double[] VectorTo(Position secondPosition)
        {
            double[] vec = new double[2];
            vec[0] = secondPosition.x - x;
            vec[1] = secondPosition.y - y;
            return vec;
        }

        // checks, if a position is within a polygon ████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if a position is within a polygon (given by corner positions) (Ray cast to the right)
        /// </summary>
        /// <param name="positions">list of corners of the polygon</param>
        /// <param name="edgeCountsAsInside">if point is on edge/corner, does it count as inside?</param>
        /// <returns></returns>
        public bool InPolygon(List<Position> positions, bool edgeCountsAsInside)
        {
            // Ray-Cast to the Right
            // "<=" und "=>" means return true, if point is on edge or corner
            // "<" und ">" means return false, if point is on edge or corner

            int i, j;
            bool c = false;

            if (edgeCountsAsInside)
            {
                for (i = 0, j = positions.Count() - 1; i < positions.Count(); j = i++) // i start from 0, j starts from top, j gets old i, i is increased by 1 => i and j are always two consecutive numbers, which both are smaller than the total number (incl first with last (will be executed as first iteration))
                {
                    if (((positions[i].y >= y) != (positions[j].y >= y)) &&
                        (x <= (positions[j].x - positions[i].x) * (y - positions[i].y)
                        / (positions[j].y - positions[i].y) + positions[i].x))
                        c = !c;
                }
            }
            else
            {
                for (i = 0, j = positions.Count() - 1; i < positions.Count(); j = i++)
                {
                    if (((positions[i].y > y) != (positions[j].y > y)) &&
                        (x < (positions[j].x - positions[i].x) * (y - positions[i].y)
                        / (positions[j].y - positions[i].y) + positions[i].x))
                        c = !c;
                }
            }
            return c;

            // several methods here:
            // http://geomalgorithms.com/a03-_inclusion.html
        }
        /// <summary>
        /// checks, if a position is within a polygon (given by corner positions) (Ray cast to the right)
        /// </summary>
        /// <param name="positions">list of corners of the polygon</param>
        /// <param name="edgeCountsAsInside">if point is on edge/corner, does it count as inside? (takes a lot more time if true or false!!!)</param>
        /// <param name="absoluteToleranceForEdge">how far away (absolute value) can the point be away to count as "on edge"</param>
        /// <returns></returns>
        public bool InPolygonWinding(List<Position> positions, bool? edgeCountsAsInside = null, double absoluteToleranceForEdge = 1e-12)
        {
            /*
            
            Position2D poly1 = new Position2D(0, 0);
            Position2D poly2 = new Position2D(1, 1);
            Position2D poly3 = new Position2D(3, 2);
            Position2D poly4 = new Position2D(2, 4);
            Position2D poly5 = new Position2D(0, 4);
            List<Position2D> polygon = new List<Position2D>() { poly1, poly2, poly3, poly4, poly5 };
            Point point1 = new PointCell(); point1.position = poly1;
            Point point3 = new PointCell(); point3.position = poly3;
            Point point5 = new PointCell(); point5.position = poly5;
            Triangle triangle = new Triangle(0, point1, point3, point5, new List<int>());

            Position2D p1 = new Position2D(-1, 5);
            Position2D p2 = new Position2D(-1, 4);
            Position2D p3 = new Position2D(1, 3);
            Position2D p4 = new Position2D(0, 3);
            Position2D p5 = new Position2D(2, 2);
            Position2D p6 = new Position2D(3, 1);
            Position2D p7 = new Position2D(1, 1);
            Position2D p8 = new Position2D(-1, 1);
            Position2D p9 = new Position2D(-100, 0);
            Position2D p10 = new Position2D(2, 4);
            Position2D p11 = new Position2D(0, 4);

            Console.WriteLine(p1.InPolygon(triangle));

            Console.WriteLine("1\t" + p1.InPolygonWinding(polygon, false) + "\t" + p1.InPolygonWinding(polygon, true));
            Console.WriteLine("2\t" + p2.InPolygonWinding(polygon, false) + "\t" + p2.InPolygonWinding(polygon, true));
            Console.WriteLine("3\t" + p3.InPolygonWinding(polygon, false) + "\t" + p3.InPolygonWinding(polygon, true));
            Console.WriteLine("4\t" + p4.InPolygonWinding(polygon, false) + "\t" + p4.InPolygonWinding(polygon, true));
            Console.WriteLine("5\t" + p5.InPolygonWinding(polygon, false) + "\t" + p5.InPolygonWinding(polygon, true));
            Console.WriteLine("6\t" + p6.InPolygonWinding(polygon, false) + "\t" + p6.InPolygonWinding(polygon, true));
            Console.WriteLine("7\t" + p7.InPolygonWinding(polygon, false) + "\t" + p7.InPolygonWinding(polygon, true));
            Console.WriteLine("8\t" + p8.InPolygonWinding(polygon, false) + "\t" + p8.InPolygonWinding(polygon, true));
            Console.WriteLine("9\t" + p9.InPolygonWinding(polygon, false) + "\t" + p9.InPolygonWinding(polygon, true));
            Console.WriteLine("10\t" + p10.InPolygonWinding(polygon, false) + "\t" + p10.InPolygonWinding(polygon, true));
            Console.WriteLine("11\t" + p11.InPolygonWinding(polygon, false) + "\t" + p11.InPolygonWinding(polygon, true));

            */
            int amountOfCorners = positions.Count;

            if (edgeCountsAsInside != null)
            {
                // Check if position is equal to any corner
                for (int i = 0; i < amountOfCorners; i++)
                    if (SameAs(positions[i], absoluteToleranceForEdge))
                        if (edgeCountsAsInside ?? false)
                            return true;
                        else
                            return false;

                // Check if position is on any edge
                for (int i = 0; i < amountOfCorners; i++)
                    if (OnLineSegmentWithTolerance(new LineSegment(positions[i], positions[(i + 1) % amountOfCorners])))
                        if (edgeCountsAsInside ?? false)
                            return true;
                        else
                            return false;
            }

            // If position is not on edge (or it doesn't matter -> edgeCountsAsInside is null) then do the winding number algorithm

            // from http://geomalgorithms.com/a03-_inclusion.html
            // on edge or not: https://stackoverflow.com/questions/37703202/winding-number-algorithm-and-point-on-boundary-edge-of-convex
            int amountofWindings = 0;

            // loop through all edges of the polygon
            for (int i = 0; i < amountOfCorners; i++)
            {   // edge from V[i] to  V[i+1]
                if (positions[i].y <= y) // start y <= P->y
                {
                    if (positions[(i + 1) % amountOfCorners].y > y)      // an upward crossing
                        if (isLeft(positions[i], positions[(i + 1) % amountOfCorners], this) > 0)  // P left of  edge
                            amountofWindings++;            // have  a valid up intersect
                }
                else// start y > P->y (no test needed)
                {
                    if (positions[(i + 1) % amountOfCorners].y <= y)     // a downward crossing
                        if (isLeft(positions[i], positions[(i + 1) % amountOfCorners], this) < 0)  // P right of  edge
                            amountofWindings--;            // have  a valid down intersect
                }
            }

            return amountofWindings != 0;

            double isLeft(Position P0, Position P1, Position P2)
            {
                return ((P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y));
            }
        }
        /// <summary>
        /// checks, if position is inside the voronoi shell of a meshpoint
        /// </summary>
        /// <param name="point">Meshpoint</param>
        /// <param name="edgeCountsAsInside">if point is on edge/corner, does it count as inside?</param>
        /// <returns></returns>
        public bool InPolygon(FiniteElement point, bool edgeCountsAsInside)
        {
            if (point.corners.Count == 2) // 1D Simulation
            {
                if (edgeCountsAsInside)
                {
                    if (x >= point.corners.Min(c => c.position.x) && x <= point.corners.Max(c => c.position.x))
                        return true;
                    else
                        return false;
                }
                else
                {
                    if (x > point.corners.Min(c => c.position.x) && x < point.corners.Max(c => c.position.x))
                        return true;
                    else
                        return false;
                }
            }

            return InPolygon(point.corners.Select(c => c.position).ToList(), edgeCountsAsInside);
        }
        /// <summary>
        /// checks, if position is inside a region
        /// </summary>
        /// <param name="region">region, which will be checked</param>
        /// <param name="edgeCountsAsInside">if point is on edge/corner, does it count as inside?</param>
        /// <returns></returns>
        public bool InPolygon(Region region, bool edgeCountsAsInside)
        {
            return InPolygon(region.orderedPoints.Select(p => p.position).ToList(), edgeCountsAsInside);
        }
        /// <summary>
        /// checks, if position is inside a triangle (edge counts always as inside)
        /// </summary>
        /// <param name="triangle">triangle</param>
        /// <returns></returns>
        public bool InPolygon(PlatformTriangle triangle)
        {
            double p0x = triangle.corners[0].position.x;
            double p0y = triangle.corners[0].position.y;
            double p1x = triangle.corners[1].position.x;
            double p1y = triangle.corners[1].position.y;
            double p2x = triangle.corners[2].position.x;
            double p2y = triangle.corners[2].position.y;

            double dX = x - p2x;
            double dY = y - p2y;
            double dX21 = p2x - p1x;
            double dY12 = p1y - p2y;
            double D = dY12 * (p0x - p2x) + dX21 * (p0y - p2y);
            double s = dY12 * dX + dX21 * dY;
            double t = (p2y - p0y) * dX + (p0x - p2x) * dY;
            if (D < 0) return s <= 0 && t <= 0 && s + t >= D;
            return s >= 0 && t >= 0 && s + t <= D;
        }

        // check, if point lies on linesegment ██████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if point lies on linesegment (with tolerance)
        /// </summary>
        /// <param name="linesegment">line segment, which will be checked</param>
        /// <returns></returns>
        public bool OnLineSegmentWithTolerance(LineSegment linesegment)
        {
            /*
            
            Check, if point is within small rectangle around the line segment --> gets tolerance

            ┌─────┐
            │  ●  |
            │  |  |
            │  |  |
            │  |  |
            │  |  |
            │  |  |
            │  |○ |
            │  |  |
            │  |  |
            │  |  |
            │  ●  |
            └─────┘

            */

            // vector of linesegment from starting point to endpoint
            double[] v = new double[] { linesegment.endPoint.x - linesegment.startingPoint.x, linesegment.endPoint.y - linesegment.startingPoint.y };

            // craete tiny vectors along and perpendiculatr to linesegment
            double[] v_p = new double[] { 1e-12 * v[0], 1e-12 * v[1] };
            double[] v_s = new double[] { v_p[1], -v_p[0] };

            // create corner of rectangle
            Position p1 = new Position(linesegment.startingPoint.x - v_p[0] - v_s[0], linesegment.startingPoint.y - v_p[1] - v_s[1]);
            Position p2 = new Position(linesegment.startingPoint.x - v_p[0] + v_s[0], linesegment.startingPoint.y - v_p[1] + v_s[1]);
            Position p3 = new Position(linesegment.endPoint.x + v_p[0] + v_s[0], linesegment.endPoint.y + v_p[1] + v_s[1]);
            Position p4 = new Position(linesegment.endPoint.x + v_p[0] - v_s[0], linesegment.endPoint.y + v_p[1] - v_s[1]);

            // now check if point in rectangle
            return InPolygon(new List<Position>() { p1, p2, p3, p4 }, true);
        }

        // get vector from position █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return the position vector of the position
        /// </summary>
        /// <returns></returns>
        public double[] ToVector()
        {
            double[] positionVector = new double[2];
            positionVector[0] = x;
            positionVector[1] = y;
            return positionVector;
        }

        // Determine Line segment bisector ██████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Return the Line, which is the Line segment bisector of the two positions
        /// </summary>
        public Line LineSegmentBisector(Position secondPosition)
        {
            // connection vector
            double[] connectionVector = VectorTo(secondPosition);

            // rotate connection vector 90°
            double[] directionVector = new double[2];
            directionVector[0] = -connectionVector[1];
            directionVector[1] = connectionVector[0];

            // center as position vector
            double[] positionVector = CenterWith(secondPosition).ToVector();

            // return Line Segment Bisector
            return new Line(positionVector, directionVector);
        }

        // Returns the quadrant index of a position █████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns the quadrant index of a position
        /// </summary>
        /// <returns></returns>
        public int GetQuadrantIndex(((double min, double max) x, (double min, double max) y) minMaxValues, int amountQuadrantsX, int amountQuadrantsY)
        {
            // structure
            //  8 /  9 / 10 / 11  ...
            //  4 /  5 /  6 /  7
            //  0 /  1 /  2 /  3
            double quadrantDistanceX = (minMaxValues.x.max - minMaxValues.x.min) / amountQuadrantsX;
            double quadrantDistanceY = (minMaxValues.y.max - minMaxValues.y.min) / amountQuadrantsY;

            int indexX = (int)((x - minMaxValues.x.min) / quadrantDistanceX);
            int indexY = (int)((y - minMaxValues.y.min) / quadrantDistanceY);

            int index = indexY * amountQuadrantsX + indexX;

            return (index >= 0 && index < amountQuadrantsX * amountQuadrantsY) ? index : amountQuadrantsX * amountQuadrantsY;
        }
        /// <summary>
        /// Returns the quadrant index of a position
        /// </summary>
        /// <returns>quadrant id</returns>
        public int GetQuadrantIndex(((double min, double max) x, (double min, double max) y, (double min, double max) z) minMaxValues, int amountQuadrantsX, int amountQuadrantsY, int amountQuadrantsZ)
        {
            // scheme
            //
            ///  24  25
            //  20  21 .......
            // 16  17  18  19
            // |
            // |  12  13  14  15
            // | 8   9   10  11
            // |4   5   6   7
            // 0   1   2   3

            double quadrantDistanceX = (minMaxValues.x.max - minMaxValues.x.min) / amountQuadrantsX;
            double quadrantDistanceY = (minMaxValues.y.max - minMaxValues.y.min) / amountQuadrantsY;
            double quadrantDistanceZ = (minMaxValues.z.max - minMaxValues.z.min) / amountQuadrantsZ;

            int indexX = (int)((x - minMaxValues.x.min) / quadrantDistanceX);
            int indexY = (int)((y - minMaxValues.y.min) / quadrantDistanceY);
            int indexZ = (int)((z - minMaxValues.z.min) / quadrantDistanceZ);

            int index = amountQuadrantsX * amountQuadrantsY * indexZ + amountQuadrantsX * indexY + indexX;

            return (index >= 0 && index < amountQuadrantsX * amountQuadrantsY * amountQuadrantsZ) ? index : amountQuadrantsX * amountQuadrantsY * amountQuadrantsZ;
        }

        // Print Position to console ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Write Position to console
        /// </summary>
        public void Print(bool newLine = true)
        {
            if (newLine)
                Console.WriteLine("(" + x + "|" + y + ")");
            else
                Console.Write("(" + x + "|" + y + ")");
        }
    }
}