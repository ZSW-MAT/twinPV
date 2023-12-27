using MoreLinq;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using BasicLib;

namespace Geometry
{
    /// <summary>
    /// Class for creating a Delaunay- and Voronoi-Mesh in 2D
    /// </summary>
    public class Meshing2D_DelaunayVoronoi<T, R> : MeshingAlgorithm<T, R>
        where T : FiniteElement, new()
        where R : Region, new()
    {
        /// <summary>
        /// adjacent triangles, of each finite element
        /// </summary>
        public Dictionary<int, List<PlatformTriangle>> adjacentTriangles { get; private set; } = new Dictionary<int, List<PlatformTriangle>>();

        /// <summary>
        /// bool, which checks if the super triangles are already removed
        /// </summary>
        bool superTrianglesRemoved { get; set; } = false;

        /// <summary>
        /// dictionary of all meshed triangles
        /// </summary>
        public Dictionary<int, PlatformTriangle> triangles { get; private set; }
        /// <summary>
        /// array of list of quadrant of all meshed triangles (for performance reasons)
        /// </summary>
        public List<int>[] triangleQuadrants { get; private set; }
        /// <summary>
        /// highest index of all simulation triangles
        /// </summary>
        public int nextAvailableTriangleIndex { get; private set; }

        /// <summary>
        /// amount of quadrant for triangles in x direction
        /// </summary>
        int amountTriangleQuadrantsX { get; set; } = 15;
        /// <summary>
        /// amount of quadrant for triangles in y direction
        /// </summary>
        int amountTriangleQuadrantsY { get; set; } = 15;

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Creates a Delaunay triangulated mesh with 4 outer points around the simulation area and two super triangles.
        /// </summary>
        /// <param name="regions">list of regions with a certain property</param>
        /// <param name="additionalContourJunctions">list of contourjunctions, which are not part of a region</param>
        public Meshing2D_DelaunayVoronoi(out Mesh<T> mesh, GeometryFileData2D geometryFileData2D) : base(geometryFileData2D)
        {
            mesh = Initialize();
        }
        /// <summary>
        /// Creates a Delaunay triangulated mesh for a single area with the given cornerpoints (should be easy to use -> meshed points accessible via delaunayVoronoi.mesh.points)
        /// </summary>
        /// <param name="cornerpoints">list of regions with a certain property (must be ordered clockwise or counter clockwise) no connection crossing allowed!</param>
        /// <param name="dPoint">mean distance of two meshpoints</param>
        public Meshing2D_DelaunayVoronoi(out Mesh<T> mesh, List<(double x, double y)> cornerpoints, int desiredAmountOfPoints)
        {
            dimension = 2;

            // define points
            Dictionary<int, ContourJunction> points = new Dictionary<int, ContourJunction>();
            for (int i = 0; i < cornerpoints.Count; i++)
                points.Add(i, new ContourJunction(i, new Position(cornerpoints[i].x, cornerpoints[i].y)));

            // define segments
            Dictionary<int, ContourSegment> segments = new Dictionary<int, ContourSegment>();
            for (int i = 0; i < cornerpoints.Count; i++)
            {
                segments.Add(i, new ContourSegment(i, points[i], points[(i + 1) % cornerpoints.Count]));
                points[i].neighbors.Add(points[(i + 1) % cornerpoints.Count]);
                points[i].adjacentContourSegments.Add(segments.Last().Value);
                points[(i + 1) % cornerpoints.Count].neighbors.Add(points[i]);
                points[(i + 1) % cornerpoints.Count].adjacentContourSegments.Add(segments.Last().Value);
            }

            // define areas
            List<R> areas = new List<R>();
            R region = new R();
            region.Initialize(0, segments.Values.OrderBy(s => s.index).ToList(), pointType.none);
            areas.Add(region);

            // write region to adjacentRegions of points and segments
            foreach (var segment in segments.Values)
            {
                if (!segment.adjacentRegions.Any(a => a.index == areas.Last().index))
                    segment.adjacentRegions.Add(areas.Last());
                if (!segment.firstAdjacentContourJunction.adjacentRegions.Any(a => a.index == areas.Last().index))
                    segment.firstAdjacentContourJunction.adjacentRegions.Add(areas.Last());
                if (!segment.secondAdjacentContourJunction.adjacentRegions.Any(a => a.index == areas.Last().index))
                    segment.secondAdjacentContourJunction.adjacentRegions.Add(areas.Last());
            }

            // initialize mesh
            mesh = Initialize();
            // Construct contour junctions and segments
            ConstructContours(desiredAmountOfPoints);
            // Add additional regular points
            mesh = AddRegularAndRandomlyShiftedPoints(mesh);
            // Add points for outer contours, simulation holes and regions
            mesh = AddContours(mesh);
            // Create Voronoi mesh from Delaunay triangulation
            mesh = CreateVoronoiFromDelaunay(mesh);
        }
        /// <summary>
        /// Creates a Delaunay triangulated mesh with 4 outer points around the simulation area and two super triangles.
        /// </summary>
        /// <param name="regions">list of regions with a certain property</param>
        /// <param name="additionalContourJunctions">list of contourjunctions, which are not part of a region</param>
        Mesh<T> Initialize()
        {
            /*
            
            ○ = outer contour points
            ● = points of the super triangles (moved 10% to the outside)

            (-4)       (-3)
             ● ───────── ●
             │\  ○       │
             │ \      ○  │
             │  \        │
             │ ○ \       │
             │    \      │
             │     \   ○ │
             │      \    │
             │ ○     \   │
             │        \  │
             │       ○ \ │
             │  ○       \│    
             ● ───────── ●
            (-1)       (-2)

            */

            // construct 4 outer points
            T p0 = new T();
            p0.index = -1;
            p0.position = new Position(minMaxValues.x.min - (minMaxValues.x.max - minMaxValues.x.min),
                minMaxValues.y.min - 0.1 * (minMaxValues.y.max - minMaxValues.y.min));

            T p1 = new T();
            p1.index = -2;
            p1.position = new Position(minMaxValues.x.max + (minMaxValues.x.max - minMaxValues.x.min),
                minMaxValues.y.min - 0.1 * (minMaxValues.y.max - minMaxValues.y.min));

            T p2 = new T();
            p2.index = -3;
            p2.position = new Position(minMaxValues.x.max + (minMaxValues.x.max - minMaxValues.x.min),
                minMaxValues.y.max + 0.1 * (minMaxValues.y.max - minMaxValues.y.min));

            T p3 = new T();
            p3.index = -4;
            p3.position = new Position(minMaxValues.x.min - (minMaxValues.x.max - minMaxValues.x.min),
                minMaxValues.y.max + 0.1 * (minMaxValues.y.max - minMaxValues.y.min));

            // construct super triangles
            PlatformTriangle d0 = new PlatformTriangle(0, p0, p1, p3, new List<int> { 1 }, minMaxValues, amountTriangleQuadrantsX, amountTriangleQuadrantsY);
            PlatformTriangle d1 = new PlatformTriangle(1, p1, p2, p3, new List<int> { 0 }, minMaxValues, amountTriangleQuadrantsX, amountTriangleQuadrantsY);

            // create entry in adjacent triangles
            adjacentTriangles.Add(p0.index, new List<PlatformTriangle>());
            adjacentTriangles.Add(p1.index, new List<PlatformTriangle>());
            adjacentTriangles.Add(p2.index, new List<PlatformTriangle>());
            adjacentTriangles.Add(p3.index, new List<PlatformTriangle>());

            // triangle to all corner points
            adjacentTriangles[p0.index].Add(d0);
            adjacentTriangles[p1.index].Add(d0);
            adjacentTriangles[p3.index].Add(d0);

            adjacentTriangles[p1.index].Add(d1);
            adjacentTriangles[p2.index].Add(d1);
            adjacentTriangles[p3.index].Add(d1);

            // create dictionaries of points and triangles
            Dictionary<int, T> points = new Dictionary<int, T>() { { -1, p0 }, { -2, p1 }, { -3, p2 }, { -4, p3 } };
            Dictionary<int, PlatformTriangle> triangles = new Dictionary<int, PlatformTriangle>() { { 0, d0 }, { 1, d1 } };
            List<int>[] triangleQuadrants = new List<int>[amountTriangleQuadrantsX * amountTriangleQuadrantsY + 1];
            for (int i = 0; i < triangleQuadrants.Length; i++)
                triangleQuadrants[i] = new List<int>();
            triangleQuadrants[d0.quadrantIndex].Add(0); // add triangle indices to their quadrant
            triangleQuadrants[d1.quadrantIndex].Add(1); // add triangle indices to their quadrant

            // set mesh
            this.triangles = triangles;
            this.triangleQuadrants = triangleQuadrants;
            this.nextAvailableTriangleIndex = 2;
            return new Mesh<T>(points, 0);
        }

        // Contours junctions and contour segments ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Determines contours junctions and contour segments
        /// </summary>  
        /// <param name="dPoint">mean distance of two meshpoints</param>
        public void ConstructContours(int desiredAmountOfPoints)
        {
            if (desiredAmountOfPoints < 1)
                desiredAmountOfPoints = 1;
            double dPoint = Math.Sqrt(outerContour.GetArea() / (double)desiredAmountOfPoints);

            // determine distances to mesh
            SetDistances(dPoint);

            // set positions of junctions
            foreach (ContourJunction contourJunction in contourJunctions)
            {
                contourJunction.SetDefaultPositionTriplets(dPointToJunctionGlobal, dPointToSegmentGlobal);
                contourJunction.SetForbiddenArea(this);
            }

            // set positions of segments
            foreach (ContourSegment contourSegment in contourSegments)
            {
                contourSegment.SetDefaultPositionTriplets(dPointToPointGlobal, dPointToJunctionGlobal, dPointToSegmentGlobal);
                contourSegment.SetForbiddenArea(this);
            }

            // reflect interference meshpoints on other segments
            bool pointAdded = true;
            while (pointAdded)
            {
                pointAdded = false;
                for (int cs = 0; cs < contourSegments.Count; cs++) // Check all segments
                    for (int ns = 0; ns < contourSegments.Count; ns++) // with all other segments
                        if (cs != ns) // except itself
                        {
                            // check all pairs in mesh points
                            foreach (var triplet in contourSegments[ns].pointTriplets)
                                if (triplet.firstMeshPos.InPolygon(contourSegments[cs].forbiddenArea.ToList(), true) // check if one or both points of the pair lay in the forbidden area
                                    || triplet.secondMeshPos.InPolygon(contourSegments[cs].forbiddenArea.ToList(), true))
                                {
                                    // Construct plumb point (parallel -> perpendicular intersection point, not parallel -> mirror point at bisectrix line)
                                    Position plumbPosition;
                                    Position interferingPosition = triplet.plumbPos;
                                    if (contourSegments[cs].lineSegment.IsParallelTo(contourSegments[ns].lineSegment)
                                        || contourSegments[cs].lineSegment.IsPerpendicularTo(contourSegments[ns].lineSegment)) // segments are parallel
                                    {
                                        plumbPosition = interferingPosition.PerpendicularPoint(contourSegments[cs].lineSegment);
                                    }
                                    else // segments are not parallel
                                    {
                                        Line mirrorLine = contourSegments[cs].lineSegment.GetBisectrixWithSmallerAngle(contourSegments[ns].lineSegment);
                                        plumbPosition = mirrorLine.Mirrorpoint(interferingPosition);
                                    }

                                    if (plumbPosition.InPolygon(contourSegments[cs].forbiddenArea.ToList(), true)) // if plumb position is in its own forbidden area
                                        if (!contourSegments[cs].pointTriplets.Any(t => t.plumbPos.SameAs(plumbPosition, defaultTolerance))) // if plumb point is not yet present in normal or additional plumb points
                                        {
                                            // do not add, if "interfering point" does not destroy the region contour
                                            Position plumbpointInPosDirection = contourSegments[cs].secondAdjacentContourJunction.position;
                                            Position plumbpointInNegDirection = contourSegments[cs].firstAdjacentContourJunction.position;
                                            foreach (var existingTriplet in contourSegments[cs].pointTriplets)
                                            {
                                                switch ((existingTriplet.plumbPos - plumbPosition).CheckVectorDirection(new Position(contourSegments[cs].lineSegment.directionVector)))
                                                {
                                                    case 1: // in positive direction
                                                        if (existingTriplet.plumbPos.DistanceTo(plumbPosition) < plumbpointInPosDirection.DistanceTo(plumbPosition))
                                                            plumbpointInPosDirection = existingTriplet.plumbPos;
                                                        break;
                                                    case -1: // in negative direction
                                                        if (existingTriplet.plumbPos.DistanceTo(plumbPosition) < plumbpointInNegDirection.DistanceTo(plumbPosition))
                                                            plumbpointInNegDirection = existingTriplet.plumbPos;
                                                        break;
                                                }
                                            }
                                            if (interferingPosition.DistanceTo(plumbpointInNegDirection) < plumbpointInPosDirection.DistanceTo(plumbpointInNegDirection)
                                                && interferingPosition.DistanceTo(plumbpointInPosDirection) < plumbpointInPosDirection.DistanceTo(plumbpointInNegDirection)) // only add if point interferes
                                            {
                                                // Add point triplet
                                                contourSegments[cs].AddPlumbPosition(plumbPosition, dPointToSegmentGlobal);

                                                // Go in while loop again
                                                pointAdded = true;
                                            }
                                        }
                                }
                        }
            }

            // order position triplets (so its easier to get voronoi edges)
            foreach (ContourSegment contourSegment in contourSegments)
                contourSegment.OrderPositionTriplets();
        }
        /// <summary>
        /// Adds meshpoints of contours junctions and contour segments
        /// </summary>  
        public Mesh<T> AddContours(Mesh<T> mesh)
        {
            foreach (var cs in contourSegments)
            {
                // corner 1 points
                mesh = AddPoint(mesh, cs.pointTriplets[0].firstMeshPos, true,
                    new List<Position> { cs.firstAdjacentContourJunction.position,
                        cs.pointTriplets[0].plumbPos.CenterWith(cs.pointTriplets[1].plumbPos) }, cs.index);
                mesh = AddPoint(mesh, cs.pointTriplets[0].secondMeshPos, true,
                    new List<Position> { cs.firstAdjacentContourJunction.position,
                        cs.pointTriplets[0].plumbPos.CenterWith(cs.pointTriplets[1].plumbPos) }, cs.index);

                // middle points
                for (int i = 1; i < cs.pointTriplets.Count - 1; i++)
                {
                    mesh = AddPoint(mesh, cs.pointTriplets[i].firstMeshPos, true,
                        new List<Position> { cs.pointTriplets[i - 1].plumbPos.CenterWith(cs.pointTriplets[i].plumbPos),
                            cs.pointTriplets[i + 1].plumbPos.CenterWith(cs.pointTriplets[i].plumbPos) }, cs.index);
                    mesh = AddPoint(mesh, cs.pointTriplets[i].secondMeshPos, true,
                        new List<Position> { cs.pointTriplets[i - 1].plumbPos.CenterWith(cs.pointTriplets[i].plumbPos),
                            cs.pointTriplets[i + 1].plumbPos.CenterWith(cs.pointTriplets[i].plumbPos) }, cs.index);
                }

                // corner 2 points
                mesh = AddPoint(mesh, cs.pointTriplets[cs.pointTriplets.Count - 1].firstMeshPos, true,
                    new List<Position> { cs.secondAdjacentContourJunction.position,
                        cs.pointTriplets[cs.pointTriplets.Count - 1].plumbPos.CenterWith(cs.pointTriplets[cs.pointTriplets.Count - 2].plumbPos) }, cs.index);
                mesh = AddPoint(mesh, cs.pointTriplets[cs.pointTriplets.Count - 1].secondMeshPos, true,
                    new List<Position> { cs.secondAdjacentContourJunction.position,
                        cs.pointTriplets[cs.pointTriplets.Count - 1].plumbPos.CenterWith(cs.pointTriplets[cs.pointTriplets.Count - 2].plumbPos) }, cs.index);
            }

            return mesh;
        }

        // Add points to mesh ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Adds regular arranged and randomly shifted points to the mesh
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <param name="randomMoveFactor">determines, how much the points are randomly moved (0: point is not moved, 1: point can be moved directly over the neighbor point)</param>
        /// <returns></returns>
        public Mesh<T> AddRegularAndRandomlyShiftedPoints(Mesh<T> mesh, double randomMoveFactor = 0.5)
        {
            // create random class
            Random random = new Random();

            // Maximum and mimimum boundaries of contour
            double MinX = outerContour.orderedPoints.Select(p => p.position).Select(e => e.x).Min();
            double MaxX = outerContour.orderedPoints.Select(p => p.position).Select(e => e.x).Max();
            double MinY = outerContour.orderedPoints.Select(p => p.position).Select(e => e.y).Min();
            double MaxY = outerContour.orderedPoints.Select(p => p.position).Select(e => e.y).Max();

            // amount of points in both directions
            int amountX = Convert.ToInt32((MaxX - MinX) / dPointToPointGlobal);
            int amountY = Convert.ToInt32((MaxY - MinY) / dPointToPointGlobal);

            // distance of points in both directions
            double deltaX = (MaxX - MinX) / amountX;
            double deltaY = (MaxY - MinY) / amountY;

            // array with all indexes from 0 to amount (random order, to get better stability, cuz no long thin triangle)
            (int x, int y)[] combinations = new (int x, int y)[amountX * amountY];
            int k = 0;
            for (int x = 0; x < amountX; x++)
                for (int y = 0; y < amountY; y++)
                    combinations[k++] = (x, y);
            combinations = combinations.OrderBy(e => Misc.random.Next()).ToArray();

            // interate through mesh in x and y direction and create regularly arranged and randomly shifted points
            foreach (var combination in combinations)
            {
                Position pos = new Position(MinX + (Convert.ToDouble(combination.x) + 0.5 + randomMoveFactor * (random.NextDouble() - 0.5))
                    * deltaX, MinY + (Convert.ToDouble(combination.y) + 0.5 + randomMoveFactor * (random.NextDouble() - 0.5)) * deltaY);
                mesh = AddPoint(mesh, pos);
            }

            return mesh;
        }

        // create Voronoi-corners from Delaunay-mesh ████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Generates a convex Voronoi-shell for each meshpoint
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <returns></returns>
        public Mesh<T> CreateVoronoiFromDelaunay(Mesh<T> mesh)
        {
            // order neighbors ——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            foreach (var point in mesh.finiteElements.Values)
                point.neighbors = point.neighbors.OrderBy(n => mesh.finiteElements[n.index].position.y).ToList();

            // remove triangles, which are more or less one line (large circumcircle) ———————————————————————————————————————————————————————————————
            for (int d = 0; d < nextAvailableTriangleIndex; d++)
                if (triangles.ContainsKey(d))
                    if (triangles[d].squaredCircumcircleRadius > 1e11)
                        RemoveTriangle(d);

            // Determine Voronoi corners ————————————————————————————————————————————————————————————————————————————————————————————————————————————
            // assign all center of the circumcirles to each points as corner
            foreach (PlatformTriangle triangle in triangles.Values)
                foreach (T point in triangle.corners)
                    point.corners.Add((triangle.circumCircleCenter, false));

            // remove triangles —————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            // remove super triangles
            if (!superTrianglesRemoved)
                mesh = RemoveSuperTriangles(mesh);

            // remove triangles outside the outer contour and outside the simulation holes
            RemoveTrianglesOutsideContour();

            foreach (T point in mesh.finiteElements.Values)
            {
                // Voronoi corners ——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
                // remove all corners outside the outer contour (including diretly ON the contour)
                point.corners.RemoveAll(e => !e.isOnContour && !e.position.InPolygon(outerContour.orderedPoints.Select(p => p.position).ToList(), false));

                // remove all corners inside simulation holes (including diretly ON the contour)
                foreach (Region hole in regions.Where(r => r.type == pointType.hole))
                    point.corners.RemoveAll(e => e.position.InPolygon(hole.orderedPoints.Select(p => p.position).ToList(), true));

                // order corners counterclockwise
                point.corners = point.corners.OrderBy(corner => Misc.Atan3(corner.position.x - point.position.x, corner.position.y - point.position.y)).ToList();
            }

            return CalculateEdgeLengthAndArea(mesh);
        }

        // refine mesh adaptive █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// refine mesh adaptive
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <param name="maximumDifference">maximum difference which is allowed between points (new point if larger)</param>
        /// <param name="between3Points">"true" = point is inserted at center of gravity of 3 points,
        /// "false" = point is inserted at center of gravity of 2 points</param>
        /// <returns></returns>
        public Mesh<T> RefineMesh(Mesh<T> mesh, Func<T, double>[] selectors, double[] maximumDifferences, bool between3Points,
            out int amountBefore, out int amountAfter)
        {
            // Selectors and maximum differences must be of same length
            if (selectors.Count() != maximumDifferences.Count())
                throw new Exception("amount of checking quantities for adaptive mesh must be same as the array of their maximum difference.");

            // amount of points before remeshing
            amountBefore = mesh.nextAvailableFiniteElementIndex;

            // list of new points, which will be added to the mesh
            List<T> newPoints = new List<T>();
            List<int[]> pointAddedBetween = new List<int[]>();

            // find all new points
            if (between3Points) // create point at center of gravity of 3 points
            {
                foreach (T point in mesh.finiteElements.Values.Where(p => p.neighbors.Count > 1))
                    for (int n = 0; n < point.neighbors.Count; n++)
                        if (point.neighbors[n].index > point.index)
                            for (int i = 0; i < selectors.Length; i++)
                            {
                                // calculate differences
                                double delta12 = Math.Abs(selectors[i](point) - selectors[i](mesh.finiteElements[point.neighbors[n].index]));
                                double delta13 = Math.Abs(selectors[i](point) - selectors[i](mesh.finiteElements[point.neighbors[(n + 1) % point.neighbors.Count].index]));
                                double delta23 = Math.Abs(selectors[i](mesh.finiteElements[point.neighbors[n].index]) - selectors[i](mesh.finiteElements[point.neighbors[(n + 1) % point.neighbors.Count].index]));

                                // get maximum difference
                                double deltaMax = Math.Max(delta12, Math.Max(delta13, delta23));

                                // add point, if maximum is too large
                                if (deltaMax > maximumDifferences[i])
                                {
                                    // postition at center of gravity of all 3
                                    double x = (point.position.x + mesh.finiteElements[point.neighbors[n].index].position.x + mesh.finiteElements[point.neighbors[(n + 1) % point.neighbors.Count].index].position.x) / 3;
                                    double y = (point.position.y + mesh.finiteElements[point.neighbors[n].index].position.y + mesh.finiteElements[point.neighbors[(n + 1) % point.neighbors.Count].index].position.y) / 3;

                                    // create new point and set mean value
                                    newPoints.Add(new T());
                                    newPoints[newPoints.Count - 1].position = new Position(x, y);
                                    for (int k = 0; k < point.GetDifferentialEquationVariables().Length; k++)
                                        newPoints[newPoints.Count - 1].SetDifferentialEquationVariable(
                                            (point.GetDifferentialEquationVariables()[k]
                                            + mesh.finiteElements[point.neighbors[n].index].GetDifferentialEquationVariables()[k]
                                            + mesh.finiteElements[point.neighbors[(n + 1) % point.neighbors.Count].index].GetDifferentialEquationVariables()[k]) / 3, k);

                                    // add point indixes of adjacent points, where the new point is placed inbetween
                                    pointAddedBetween.Add(new int[] { point.index, point.neighbors[n].index, point.neighbors[(n + 1) % point.neighbors.Count].index });
                                }
                            }
            }

            else // create point at center of gravity of 2 points
            {
                foreach (T point in mesh.finiteElements.Values)
                    foreach (var n in point.neighbors.FindAll(neighbor => neighbor.index > point.index))
                        for (int i = 0; i < selectors.Length; i++)
                            if (Math.Abs(selectors[i](point) - selectors[i](mesh.finiteElements[n.index]))
                                > maximumDifferences[i])
                            {
                                double x = (point.position.x + mesh.finiteElements[n.index].position.x) / 2;
                                double y = (point.position.y + mesh.finiteElements[n.index].position.y) / 2;

                                // create new point and set mean value
                                newPoints.Add(new T());
                                newPoints[newPoints.Count - 1].position = new Position(x, y);
                                for (int k = 0; k < point.GetDifferentialEquationVariables().Length; k++)
                                    newPoints[newPoints.Count - 1].SetDifferentialEquationVariable(
                                        (point.GetDifferentialEquationVariables()[k] + mesh.finiteElements[n.index].GetDifferentialEquationVariables()[k]) / 2, k);

                                // add point indixes of adjacent points, where the new point is placed inbetween
                                pointAddedBetween.Add(new int[] { point.index, n.index });
                            }
            }

            // add all new points
            for (int p = 0; p < newPoints.Count; p++)
            {
                int oldAmount = mesh.nextAvailableFiniteElementIndex;

                // add point
                mesh = AddPoint(mesh, newPoints[p].position);

                // if it was added (and not rejected), set parameters
                if (oldAmount + 1 == mesh.nextAvailableFiniteElementIndex)
                {
                    if (pointAddedBetween[p].Length == 2)
                        Console.WriteLine("> Point " + oldAmount + " between " + pointAddedBetween[p][0]
                            + " and " + pointAddedBetween[p][1] + " added. ("
                            + mesh.finiteElements[oldAmount].position.x + "|" + mesh.finiteElements[oldAmount].position.y + ")" + ".");
                    else
                        Console.WriteLine("> Point " + oldAmount + " between " + pointAddedBetween[p][0]
                            + ", " + pointAddedBetween[p][1] + " and " + pointAddedBetween[p][2] + " added. ("
                            + mesh.finiteElements[oldAmount].position.x + "|" + mesh.finiteElements[oldAmount].position.y + ")" + ".");

                    // Set data of new point
                    mesh.finiteElements[mesh.nextAvailableFiniteElementIndex - 1].index = mesh.nextAvailableFiniteElementIndex - 1;
                    for (int i = 0; i < mesh.finiteElements[mesh.nextAvailableFiniteElementIndex - 1].GetDifferentialEquationVariables().Length; i++)
                        mesh.finiteElements[mesh.nextAvailableFiniteElementIndex - 1].SetDifferentialEquationVariable(newPoints[p].GetDifferentialEquationVariables()[i], i);
                }
                else
                {
                    if (pointAddedBetween[p].Length == 2)
                        Console.WriteLine("> Point " + oldAmount + " between " + pointAddedBetween[p][0]
                            + " and " + pointAddedBetween[p][1] + " could not be added. ("
                            + newPoints[p].position.x + "|" + newPoints[p].position.y + ")" + ".");
                    else
                        Console.WriteLine("> Point " + oldAmount + " between " + pointAddedBetween[p][0]
                            + ", " + pointAddedBetween[p][1] + " and " + pointAddedBetween[p][2] + " could not be added. ("
                            + newPoints[p].position.x + "|" + newPoints[p].position.y + ")" + ".");
                }
            }

            amountAfter = mesh.nextAvailableFiniteElementIndex;
            return CreateVoronoiFromDelaunay(mesh);
        }

        // Calculate the edge lengths and the area from the Voronoi-corners █████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculates the edge lengths and the area from the Voronoi-corners
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <returns></returns>
        Mesh<T> CalculateEdgeLengthAndArea(Mesh<T> mesh)
        {
            foreach (T point in mesh.finiteElements.Values)
            {
                // add edge length for each neighbor
                for (int n_index = 0; n_index < point.neighbors.Count; n_index++)
                {
                    // list of triangles, which border to both points
                    List<PlatformTriangle> commonTriangles = adjacentTriangles[point.index].Intersect(adjacentTriangles[point.neighbors[n_index].index]).ToList();

                    // variable for edge length between to points
                    double edgeLength = -1;

                    // if both cells are in middle (not at the rim of the simulation area), the have 2 common triangles
                    // then their edge length is the distance between the 2 centers of the circumcircles of the triangles
                    if (commonTriangles.Count == 2)
                        edgeLength = commonTriangles.First().circumCircleCenter.DistanceTo(commonTriangles.Last().circumCircleCenter);
                    else
                    {
                        // list of common voronoi corners
                        List<Position> commonCorners = new List<Position>();
                        foreach (var corner1 in point.corners)
                            foreach (var corner2 in mesh.finiteElements[point.neighbors[n_index].index].corners)
                                if (corner1.position.SameAs(corner2.position, defaultTolerance))
                                    commonCorners.Add(corner1.position);

                        // if exactly 2 corners are the same, the edge length is the distance between those 2
                        if (commonCorners.Count == 2)
                            edgeLength = commonCorners[0].DistanceTo(commonCorners[1]);
                        else
                        {
                            // list of common voronoi corners, which are on the contour
                            List<Position> commonContourCorners = new List<Position>();
                            foreach (var corner1 in point.corners.Where(c => c.isOnContour))
                                foreach (var corner2 in mesh.finiteElements[point.neighbors[n_index].index].corners.Where(c => c.isOnContour))
                                    if (corner1.position.SameAs(corner2.position, defaultTolerance))
                                        commonContourCorners.Add(corner1.position);

                            // if both cells are on the rim of the contour, they have only one common triangle
                            // edge length is the distance from common corner on contour to center of circumcircle of common triangle
                            if (commonTriangles.Count == 1)
                            {
                                if (commonContourCorners.Count == 0)
                                {
                                    edgeLength = 1e-50;
                                    Console.WriteLine("coord. point: " + point.position.x + " " + point.position.y + " / " + mesh.finiteElements[point.neighbors[n_index].index].position.x + " " + mesh.finiteElements[point.neighbors[n_index].index].position.y + "..");
                                    Console.WriteLine("Point " + point.index + " and point " + point.neighbors[n_index].index + " do not have common contour corners! => Edge length set to 1e-50.");
                                    //throw new Exception("Point " + point.index + " and point " + n + " do not have common contour corners!");
                                }
                                else
                                    edgeLength = commonContourCorners.First().DistanceTo(commonTriangles.First().circumCircleCenter);
                            }
                        }
                    }

                    if (edgeLength == -1)
                        throw new Exception("Edge length from point " + point.index + " to point " + point.neighbors[n_index].index + " could not be set!");

                    // add edge length on same position as neighbor index in point.N[]
                    point.neighbors[n_index] = (point.neighbors[n_index].index, point.position.DistanceTo(mesh.finiteElements[point.neighbors[n_index].index].position), edgeLength);
                }

                // calculate volume with the Shoelace formula (Gaussian trapezoidal formula)
                // https://en.wikipedia.org/wiki/Shoelace_formula
                double sum = 0;
                for (int k = 0; k < point.corners.Count; k++)
                    sum += (point.corners[k].position.x + point.corners[(k + 1) % point.corners.Count].position.x) * (point.corners[(k + 1) % point.corners.Count].position.y - point.corners[k].position.y);
                point.size = sum / 2;
            }

            return mesh;
        }

        // remove all triangles █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Remove super triangles
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <returns></returns>
        Mesh<T> RemoveSuperTriangles(Mesh<T> mesh)
        {
            // set bool to true -> triangles are removed only once
            superTrianglesRemoved = true;

            // remove outer points
            for (int p = -4; p < 0; p++)
                mesh = RemovePoint(mesh, p);

            return mesh;
        }
        /// <summary>
        /// removes all triangles, which are outside the contour
        /// </summary>
        /// <param name="mesh">mesh-container with all points and triangles</param>
        /// <returns></returns>
        void RemoveTrianglesOutsideContour()
        {
            // triangles outside mesharea
            for (int d = 0; d < nextAvailableTriangleIndex; d++)
                if (triangles.ContainsKey(d))
                    if (!triangles[d].CompletelyInPolygon(outerContour.orderedPoints.Select(p => p.position).ToList()))
                        RemoveTriangle(d);

            // triangles in simulation holes
            for (int d = 0; d < nextAvailableTriangleIndex; d++)
                if (triangles.ContainsKey(d))
                    foreach (Region hole in regions.Where(r => r.type == pointType.hole))
                        if (triangles[d].PartlyInPolygon(hole.orderedPoints.Select(p => p.position).ToList())
                            || triangles[d].circumCircleCenter.InPolygon(hole.orderedPoints.Select(p => p.position).ToList(), true))
                        {
                            RemoveTriangle(d);
                            break;
                        }
        }

        List<Position> island1 = new List<Position>();
        List<Position> island2 = new List<Position>();

        // Add point to the mesh via Bowyer-Watson-Algorithm ████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Adds a new point to the mesh via Bowyer-Watson-Algorithm
        /// </summary>
        /// <param name="newPosition">position of the new points</param>
        /// <param name="allowDestroyKontur">true, if point is border-point and hence allowed to "destroy" contours</param>
        /// <param name="cornersOnContour">if non-null: list of voronoi corners, that will be added in advance</param>
        /// <returns></returns>
        public Mesh<T> AddPoint(Mesh<T> mesh, Position newPosition, bool allowDestroyKontur = false,
            List<Position> cornersOnContour = null, int indexOfSegmentCreatedFrom = -1)
        {
            // Only if position is within the outer contour and there is no other point at this position (with tolerance)
            if (newPosition.InPolygon(outerContour, false)
                && !mesh.finiteElements.Values.Any(p => p.position.SameAs(newPosition, defaultTolerance)))
            {
                if (!newPosition.InPolygon(island1, true) && !newPosition.InPolygon(island2, true))
                    // Only if new position is not in a simulation hole
                    if (!regions.Where(r => r.type == pointType.hole).Any(h => newPosition.InPolygon(h.orderedPoints.Select(p => p.position).ToList(), true)))
                    {
                        // Only if new point is allowed to "destroy" contour" (points, which are included to model countours itself)
                        // or if no contour segment and no contour junction is destroyed
                        if (allowDestroyKontur ||
                            (!contourSegments.Any(cs => cs.IsInForbiddenArea(newPosition))
                            && !contourJunctions.Where(cj => cj.adjacentContourSegments.Count > 0).Any(cj => cj.IsInForbiddenArea(newPosition))))
                        {
                            // List of triangles, in whose circumcircle the new point is (at least the triangle, in which the poit actually is)
                            List<PlatformTriangle> N0 = new List<PlatformTriangle> { triangles[GetAnyTriangleWithinCircumcircle(newPosition)] };

                            // Liste of triangles, which need to be checked (only neighbor triangles of N0 triangles) (see JF 21.10.2019)
                            List<PlatformTriangle> N1 = N0[0].adjacentTriangles.Select(n => triangles[n]).ToList();

                            // check all (planar adjacent) neighbor triangles
                            while (N1.Count > 0)
                            {
                                // check, if position is within their circumcircle of the triangle
                                if (N1[0].PositionWithinCircumcircle(newPosition))
                                {
                                    // add positiv check triangles to N0
                                    N0.Add(N1[0]);

                                    // check adjacent triangles of the added triangle, as well (if they are not already in N0 or N1)
                                    N1.AddRange(N1[0].adjacentTriangles.FindAll(n => !N0.Any(d => d.index == n) && !N1.Any(d => d.index == n))
                                        .Select(n => triangles[n]).ToList());
                                }

                                // removed checked triangle from checklist
                                N1.RemoveAt(0);
                            }

                            // create new point from position and add to mesh
                            T newPoint = new T();
                            newPoint.index = mesh.nextAvailableFiniteElementIndex;
                            newPoint.position = newPosition;
                            newPoint.corners = cornersOnContour == null ? new List<(Position position, bool inOnContour)>() : cornersOnContour.Select(c => (c, true)).ToList();
                            if (cornersOnContour != null)
                                if (cornersOnContour.Count == 2)
                                    newPoint.borderEdgeSize = cornersOnContour[0].DistanceTo(cornersOnContour[1]);
                            newPoint.indexOfBorderElementCreatedFrom = indexOfSegmentCreatedFrom;
                            mesh.AddPoint(newPoint);

                            // create entry in adjacent triangles
                            adjacentTriangles.Add(newPoint.index, new List<PlatformTriangle>());

                            // Output amount of points
                            if (mesh.nextAvailableFiniteElementIndex % 10000 == 0)
                            {
                                Console.ForegroundColor = ConsoleColor.White;
                                Console.Write("██  ");
                                Console.ForegroundColor = ConsoleColor.Gray;
                                Console.WriteLine("Current amount of points in Delaunay algorithm: " + mesh.nextAvailableFiniteElementIndex);
                            }

                            // list of edges, wich encompasses all removed N0 triangles
                            List<Edge> polygon = CutoutPolygon(N0);

                            // remove all N0 triangles from mesh container
                            foreach (PlatformTriangle n0 in N0)
                                RemoveTriangle(n0.index);

                            // add a new triangle for each edge of the removed polygon
                            int amountTrianglesBefore = nextAvailableTriangleIndex;
                            List<int> toBeRemovedTriangles = new List<int>();
                            for (int k = 0; k < polygon.Count; k++)
                            {
                                // neighbor indexes of the new triangle
                                List<int> neighboringTriangles = new List<int> { amountTrianglesBefore + Misc.mod(k - 1, polygon.Count),
                                amountTrianglesBefore + Misc.mod(k + 1, polygon.Count)};

                                // If edge of an triangle has an outer adjacent triangle, this is also new neighbor of the new triangle
                                if (polygon[k].indexTriangleOutside != -1)
                                {
                                    neighboringTriangles.Add(polygon[k].indexTriangleOutside);
                                    triangles[polygon[k].indexTriangleOutside].adjacentTriangles.Add(amountTrianglesBefore + k);
                                }

                                // create and add triangle
                                PlatformTriangle triangle = new PlatformTriangle(amountTrianglesBefore + k, polygon[k].startpoint,
                                    polygon[k].endpoint, newPoint, neighboringTriangles, minMaxValues, amountTriangleQuadrantsX, amountTriangleQuadrantsY);
                                AddTriangle(triangle);

                                // triangle to all corner points
                                adjacentTriangles[polygon[k].startpoint.index].Add(triangle);
                                adjacentTriangles[polygon[k].endpoint.index].Add(triangle);
                                adjacentTriangles[newPoint.index].Add(triangle);
                            }
                        }
                    }
            }

            return mesh;
        }

        // Removes point from mesh ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Adds a new point to the mesh via Bowyer-Watson-Algorithm
        /// </summary>
        public Mesh<T> RemovePoint(Mesh<T> mesh, int index)
        {
            // remove all adjacent triangles of point
            for (int t = 0; t < adjacentTriangles[index].Count; t++)
            {
                RemoveTriangle(adjacentTriangles[index][t].index);
                t--;
            }

            // remove point
            mesh.RemovePoint(index);

            return mesh;
        }

        // Returns index if triangle, wich includes a position ██████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Returns index if triangle, wich includes a position
        /// </summary>
        /// <param name="triangles">list of triangles, which will be scanned</param>
        /// <param name="position">Position, whose location will be checked</param>
        /// <returns></returns>
        int GetAnyTriangleWithinCircumcircle(Position position)
        {
            int quadrantIndexPosition = position.GetQuadrantIndex(minMaxValues, amountTriangleQuadrantsX, amountTriangleQuadrantsY);
            int quadrantPosX = quadrantIndexPosition % amountTriangleQuadrantsX;
            int quadrantPosY = quadrantIndexPosition / amountTriangleQuadrantsX;

            // iterate over own quadrant
            foreach (int t in triangleQuadrants[quadrantIndexPosition])
                if (triangles[t].PositionWithinCircumcircle(position))
                    return t;

            // iterate over quadrant with chebyshev distance == 2 and == 3
            for (int distance = 1; distance <= 3; distance++)
            {
                // get list with all quadrant indexes
                List<int> quadrantIndexesToCheck = new List<int>();
                for (int i = 0; i < amountTriangleQuadrantsX * amountTriangleQuadrantsY; i++)
                {
                    int quadrantX = i % amountTriangleQuadrantsX;
                    int quadrantY = i / amountTriangleQuadrantsX;
                    if (Math.Max(Math.Abs(quadrantPosX - quadrantX), Math.Abs(quadrantPosY - quadrantY)) == distance)
                        quadrantIndexesToCheck.Add(i);
                }

                // check triangles in this quadrants
                foreach (int quadrant in quadrantIndexesToCheck)
                    foreach (int t in triangleQuadrants[quadrant])
                        if (triangles[t].PositionWithinCircumcircle(position))
                            return t;
            }

            // iterate over outside quadrant
            foreach (int t in triangleQuadrants[amountTriangleQuadrantsX * amountTriangleQuadrantsY])
                if (triangles[t].PositionWithinCircumcircle(position))
                    return t;

            // iterate over all
            foreach (PlatformTriangle triangle in triangles.Values)
                if (triangle.PositionWithinCircumcircle(position))
                    return triangle.index;

            // if no triangle does exactly contain the position, search the triangle, whose edge is next to the point
            double minimumDistanceToEdge = double.PositiveInfinity;
            int triangleIndex = -5;
            foreach (PlatformTriangle triangle in triangles.Values)
            {
                LineSegment edge1 = new LineSegment(triangle.corners[0].position, triangle.corners[1].position);
                LineSegment edge2 = new LineSegment(triangle.corners[1].position, triangle.corners[2].position);
                LineSegment edge3 = new LineSegment(triangle.corners[2].position, triangle.corners[0].position);

                double distance = position.DistanceTo(edge1);
                if (distance < minimumDistanceToEdge)
                {
                    triangleIndex = triangle.index;
                    minimumDistanceToEdge = distance;
                }
            }
            return triangleIndex;
        }

        // Get border for hole, which is created by removed triangles ███████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// return list of edges, which border the hole-polygon, which was created by removing triangles
        /// </summary>
        /// <param name="removedTriangles">list of removed triangles</param>
        /// <returns></returns>
        List<Edge> CutoutPolygon(List<PlatformTriangle> removedTriangles)
        {
            // list of edges, which border the polygon
            List<Edge> edges = new List<Edge>();

            // add all edges of triangles
            // if edge is called twice, the edge is not at the outside but within the inside of the polygon
            // => if edge is already there, remove it cuz its not a outer border!
            foreach (PlatformTriangle n0 in removedTriangles.Distinct())
            {
                Edge edge0 = new Edge(n0.corners[0], n0.corners[1], n0.TriangleWithCommonBorders(n0.corners[0], n0.corners[1], adjacentTriangles));
                if (edges.FindAll(k => k.Equals(edge0)).Count() == 1)
                    edges.RemoveAt(edges.FindIndex(k => k.Equals(edge0)));
                else
                    edges.Add(edge0);

                Edge edge1 = new Edge(n0.corners[1], n0.corners[2], n0.TriangleWithCommonBorders(n0.corners[1], n0.corners[2], adjacentTriangles));
                if (edges.FindAll(k => k.Equals(edge1)).Count() == 1)
                    edges.RemoveAt(edges.FindIndex(k => k.Equals(edge1)));
                else
                    edges.Add(edge1);

                Edge edge2 = new Edge(n0.corners[2], n0.corners[0], n0.TriangleWithCommonBorders(n0.corners[2], n0.corners[0], adjacentTriangles));
                if (edges.FindAll(k => k.Equals(edge2)).Count() == 1)
                    edges.RemoveAt(edges.FindIndex(k => k.Equals(edge2)));
                else
                    edges.Add(edge2);
            }

            // order edges in a circle (so later neighbor indexes can be determined easily)
            List<Edge> orderedEdges = new List<Edge>();

            // pick a random edge and move it from the old list to the ordered list
            orderedEdges.Add(edges[0]);
            edges.RemoveAt(0);

            // move edges to ordered list, until old list ist empty
            while (edges.Count > 0)
            {
                // search for neighbor from list, which hasn't been used yet
                // from the second edge: only one neighbor is left, cuz all used elements have been removed => edges ordered in a circle
                Edge nextEdge = orderedEdges.Last().AdjacentEdges(edges).First();

                // move neighbor edge from unordered list to oredered list
                orderedEdges.Add(nextEdge);
                edges.Remove(nextEdge);
            }

            return orderedEdges;
        }

        // Sets the 3 different distances ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Sets the 3 different distances dPointsGlobal (user defined), dPointsGlobalJunctions (determined via minimal distance between junctions) and dSegmentGlobal (dPointsGlobalJunctions/60)
        /// </summary>
        /// <param name="dPoint">user defined minimum distance between two meshpoints</param>
        /// <returns></returns>
        void SetDistances(double dPoint)
        {
            dPointToPointGlobal = dPoint;

            if (contourJunctions == null || contourSegments == null)
            {
                dPointToJunctionGlobal = dPointToPointGlobal;
                dPointToSegmentGlobal = dPointToJunctionGlobal / 60;
            }
            else
            {
                double minimumDistanceBetweenAnyContourJunctions = contourJunctions.Min(cj => cj.neighbors.Where(j => j.index > cj.index)
                    .Select(n => n.position.DistanceTo(cj.position)).DefaultIfEmpty(double.PositiveInfinity).Min());
                double minimumDistanceBetweenAnyJunctionAndSegment = contourJunctions.Min(cj => contourSegments
                    .Where(cs => cs.firstAdjacentContourJunction.index != cj.index && cs.secondAdjacentContourJunction.index != cj.index)
                    .Min(cs => cj.position.DistanceTo(cs.lineSegment)));
                double minimumContours = 0.4 * Math.Min(minimumDistanceBetweenAnyContourJunctions, minimumDistanceBetweenAnyJunctionAndSegment);
                dPointToJunctionGlobal = Math.Min(dPointToPointGlobal, minimumContours);
                dPointToSegmentGlobal = dPointToJunctionGlobal / 60;
            }
        }

        // Add new triangle to mesh █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Add new triangle to mesh
        /// </summary>
        public void AddTriangle(PlatformTriangle triangle)
        {
            triangles.Add(nextAvailableTriangleIndex, triangle);
            triangleQuadrants[triangle.quadrantIndex].Add(triangle.index);
            nextAvailableTriangleIndex++;
        }

        // Remove triangle from mesh ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Removes a single triangle from the mesh-container
        /// </summary>
        /// <param name="index">index of the triangle, which shell be removed</param>
        /// <returns></returns>
        public void RemoveTriangle(int index)
        {
            // Remove triangle from the list of neighboring of all adjecent triangles
            foreach (int n in triangles[index].adjacentTriangles)
                triangles[n].adjacentTriangles.Remove(index);

            // remove points from each others list of neighbors, if their connection was only via the removed triangle
            for (int i = 0; i < 3; i++)
                if (adjacentTriangles[triangles[index].corners[i].index].Select(a => a.index)
                    .Intersect(adjacentTriangles[triangles[index].corners[Misc.mod(i + 1, 3)].index].Select(a => a.index))
                    .Where(a => a != index).Count() == 0)
                {
                    triangles[index].corners[i].neighbors.RemoveAll(e => e.index == triangles[index].corners[Misc.mod(i + 1, 3)].index);
                    triangles[index].corners[Misc.mod(i + 1, 3)].neighbors.RemoveAll(e => e.index == triangles[index].corners[i].index);
                }

            // remove triangle from list of adjacent triangles of all points
            foreach (T point in triangles[index].corners)
                adjacentTriangles[point.index].Remove(triangles[index]);

            // remove triangle from dictionary
            triangleQuadrants[triangles[index].quadrantIndex].Remove(index);
            triangles.Remove(index);
        }
    }
}