using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using BasicLib;

namespace Geometry
{
    /// <summary>
    /// Class, which represents a triangle made out of 3 points
    /// </summary>
    public class PlatformTriangle
    {
        /// <summary>
        /// index of the triangle
        /// </summary>
        public int index { get; private set; }

        /// <summary>
        /// number telling in which quadrant the center of the circum circle of this triangle is in (for finding the triangle faster)
        /// </summary>
        public int quadrantIndex { get; private set; }

        /// <summary>
        /// array with all three corner points
        /// </summary>
        public FiniteElement[] corners { get; private set; } = new FiniteElement[3];

        /// <summary>
        /// adjacent triangles to this triangle
        /// </summary>
        public List<int> adjacentTriangles { get; private set; } = new List<int>();

        /// <summary>
        /// Center of the circum circle
        /// </summary>
        public Position circumCircleCenter;

        /// <summary>
        /// squared radius of the circum circle of this triangle
        /// </summary>
        public double squaredCircumcircleRadius { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor of a triangle
        /// </summary>
        /// <param name="index">index of the new triangle</param>
        /// <param name="point1">first corner point</param>
        /// <param name="point2">second corner point</param>
        /// <param name="point3">third corner point</param>
        /// <param name="adjacentTriangles">list of adjacent triangles</param>
        public PlatformTriangle(int index, FiniteElement point1, FiniteElement point2, FiniteElement point3, List<int> adjacentTriangles,
            ((double min, double max) x, (double min, double max) y) minMaxValues, int amountQuadrantsX, int amountQuadrantsY)
        {
            // set index
            this.index = index;

            // set neighbors
            this.adjacentTriangles = adjacentTriangles;

            // set corner points
            corners[0] = point1;
            corners[1] = point2;
            corners[2] = point3;

            // add corner points as neighbors of each other
            point1.neighbors.Add((point2.index, 0, 0));
            point1.neighbors.Add((point3.index, 0, 0));
            point1.neighbors = point1.neighbors.Distinct().ToList();
            point2.neighbors.Add((point1.index, 0, 0));
            point2.neighbors.Add((point3.index, 0, 0));
            point2.neighbors = point2.neighbors.Distinct().ToList();
            point3.neighbors.Add((point1.index, 0, 0));
            point3.neighbors.Add((point2.index, 0, 0));
            point3.neighbors = point3.neighbors.Distinct().ToList();

            // update circumcircle
            UpdateCircumCircle();

            quadrantIndex = circumCircleCenter.GetQuadrantIndex(minMaxValues, amountQuadrantsX, amountQuadrantsY);
        }

        // update the circumcircle ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// updates the radius and the position of the circum circle
        /// </summary>
        private void UpdateCircumCircle()
        {
            // https://codefound.wordpress.com/2013/02/21/how-to-compute-a-circumcircle/#more-58
            // https://en.wikipedia.org/wiki/Circumscribed_circle

            // auxiliary variables
            double A0 = corners[0].position.x * corners[0].position.x + corners[0].position.y * corners[0].position.y;
            double A1 = corners[1].position.x * corners[1].position.x + corners[1].position.y * corners[1].position.y;
            double A2 = corners[2].position.x * corners[2].position.x + corners[2].position.y * corners[2].position.y;

            double aux1 = A0 * (corners[2].position.y - corners[1].position.y) + A1 * (corners[0].position.y - corners[2].position.y)
                + A2 * (corners[1].position.y - corners[0].position.y);
            double aux2 = -(A0 * (corners[2].position.x - corners[1].position.x) + A1 * (corners[0].position.x - corners[2].position.x)
                + A2 * (corners[1].position.x - corners[0].position.x));
            double div = 2 * (corners[0].position.x * (corners[2].position.y - corners[1].position.y)
                            + corners[1].position.x * (corners[0].position.y - corners[2].position.y)
                            + corners[2].position.x * (corners[1].position.y - corners[0].position.y));

            // if denominator == 0, triangle is colinear
            if (div == 0)
            {
                throw new Exception("Corners of the triangle:\n" + corners[0].position.x + " | " + corners[0].position.y + "\n"
                     + corners[1].position.x + " | " + corners[1].position.y + "\n" + corners[2].position.x + " | " + corners[2].position.y);
            }

            // calculate center of circum circle (link from Wikipedia)
            circumCircleCenter = new Position(aux1 / div, aux2 / div);

            // calculate squared radius, (distance from center to any corner)
            squaredCircumcircleRadius = (circumCircleCenter.x - corners[0].position.x) * (circumCircleCenter.x - corners[0].position.x)
                + (circumCircleCenter.y - corners[0].position.y) * (circumCircleCenter.y - corners[0].position.y);
        }

        // check if position is inside circumcircle █████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// check if a position is inside the circumcircle
        /// </summary>
        /// <param name="position">position which will be checked</param>
        /// <returns></returns>
        public bool PositionWithinCircumcircle(Position position)
        {
            // calculate squared radius
            double deltaX = position.x - circumCircleCenter.x;
            double deltaY = position.y - circumCircleCenter.y;

            var squaredDistanceFromPositionToCircumCircleCenter = deltaX * deltaX + deltaY * deltaY;

            // return, if distance is smaller than radius
            return squaredDistanceFromPositionToCircumCircleCenter < squaredCircumcircleRadius;
        }

        // check, if triangle is completely in a polygon ████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if the triangle is completely in a polygon
        /// </summary>
        /// <param name="positions">list of positions of the polygon</param>
        /// <returns></returns>
        public bool CompletelyInPolygon(List<Position> positions)
        {
            // create linesegments of consecute positions of the triangle
            LineSegment edge1 = new LineSegment(corners[0].position, corners[1].position);
            LineSegment edge2 = new LineSegment(corners[1].position, corners[2].position);
            LineSegment edge3 = new LineSegment(corners[2].position, corners[0].position);

            for (int p = 0; p < positions.Count; p++)
            {
                // create linesegment of consecutive positions of the polygon
                LineSegment contour = new LineSegment(positions[p], positions[Misc.mod(p + 1, positions.Count)]);

                // if there are intersections, return false
                if (contour.DoesIntersect(edge1, false) || contour.DoesIntersect(edge2, false) || contour.DoesIntersect(edge3, false))
                    return false;

                // now, triangle is completely out or completely in the polygon
                // -> check if ANY point of each linesegment is inside the polygon
                if (!edge1.startingPoint.CenterWith(edge1.endPoint).InPolygon(positions, false)
                    || !edge2.startingPoint.CenterWith(edge2.endPoint).InPolygon(positions, false)
                    || !edge3.startingPoint.CenterWith(edge3.endPoint).InPolygon(positions, false))
                    return false;
            }
            return true;
        }

        // check if, triangle is partly in polygon ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// checks, if the triangle is partly in a polygon
        /// </summary>
        /// <param name="positionen">list of positions of the polygon</param>
        /// <returns></returns>
        public bool PartlyInPolygon(List<Position> positionen)
        {
            // create linesegments of consecute positions of the triangle
            LineSegment edge1 = new LineSegment(corners[0].position, corners[1].position);
            LineSegment edge2 = new LineSegment(corners[1].position, corners[2].position);
            LineSegment edge3 = new LineSegment(corners[2].position, corners[0].position);

            for (int p = 0; p < positionen.Count; p++)
            {
                // create linesegment of consecutive positions of the polygon
                LineSegment contour = new LineSegment(positionen[p], positionen[Misc.mod(p + 1, positionen.Count)]);

                // if there are intersections, return false
                if (contour.DoesIntersect(edge1, false) || contour.DoesIntersect(edge2, false) || contour.DoesIntersect(edge3, false))
                    return true;

                // now, triangle is completely out or completely in the polygon
                // -> check if ANY point of each linesegment is inside the polygon
                if (edge1.startingPoint.CenterWith(edge1.endPoint).InPolygon(positionen, false)
                    || edge2.startingPoint.CenterWith(edge2.endPoint).InPolygon(positionen, false)
                    || edge3.startingPoint.CenterWith(edge3.endPoint).InPolygon(positionen, false))
                    return true;
            }
            return false;
        }

        // return triangle index, which borders to this triangle and two given points ███████████████████████████████████████████████████████████████
        /// <summary>
        /// return triangle index, which borders to this triangle and two given points
        /// returns "-1", if corner points do no border to this triangle or if there is no opposite triangle (e.g. triangle at the rim)
        /// </summary>
        /// <param name="cornerPoint1">first point, which needs to be common of both triangles</param>
        /// <param name="cornerPoint2">second point, which needs to be common of both triangles</param>
        /// <returns></returns>
        public int TriangleWithCommonBorders(FiniteElement cornerPoint1, FiniteElement cornerPoint2, Dictionary<int, List<PlatformTriangle>> adjacentTriangles)
        {
            // return opposite triangle
            foreach (PlatformTriangle triangleOfPoint1 in adjacentTriangles[cornerPoint1.index])
                if (triangleOfPoint1.index != index)
                    if (adjacentTriangles[cornerPoint2.index].Any(t => t.index == triangleOfPoint1.index))
                        return triangleOfPoint1.index;

            // return "-1" of opposite triangle does not exist
            return -1;
        }

        // print to console █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// print to console
        /// </summary>
        public void Print()
        {
            Console.WriteLine("triangle " + index);

            Console.WriteLine("corners:");
            foreach (FiniteElement eck in corners)
            {
                Console.Write("\t");
                eck.position.Print();
            }

            Console.Write("adjavent triangles:");
            foreach (int indexTriangles in adjacentTriangles)
                Console.Write(" " + indexTriangles);
            Console.WriteLine();

            Console.Write("Center of circumcircle: ");
            circumCircleCenter.Print();

            Console.WriteLine("Squared radius = " + squaredCircumcircleRadius + "\n");
        }
    }
}