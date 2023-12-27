using BasicLib;
using MoreLinq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry
{
    public class MeshingAlgorithm<T, R>
        where T : FiniteElement, new()
        where R : Region, new()
    {
        /// <summary>
        /// outer contour of all points
        /// </summary>
        public Region outerContour { get; set; }
        /// <summary>
        /// meshed area (total area minus simulation holes) (do NOT use this area for efficiency calculations!!! there might also areas wich do NOT cound as active area)
        /// </summary>
        public double totalSize { get; set; }

        /// <summary>
        /// arrays telling the borders of the single quadrants
        /// </summary>
        public ((double min, double max) x, (double min, double max) y) minMaxValues;

        /// <summary>
        /// default tolerance, which is allowed between two positions to be considered as the same position
        /// </summary>
        public double defaultTolerance { get; set; }

        /// <summary>
        /// determines the dimension of this meshing algorithm
        /// </summary>
        public int dimension { get; set; }

        /// <summary>
        /// list of regions
        /// </summary>
        public List<R> regions { get; set; } = new List<R>();

        /// <summary>
        /// list of contour junctions, which represent the outer contours, the holes and the regions
        /// </summary>
        public List<ContourJunction> contourJunctions { get; set; }
        /// <summary>
        /// list of countour segments, which represent the outer contours, the holes and the regions
        /// </summary>
        public List<ContourSegment> contourSegments { get; set; }

        /// <summary>
        /// global distance between two finite elements
        /// </summary>
        public double dPointToPointGlobal { get; set; }
        /// <summary>
        /// global distance from finite element to contour junction
        /// </summary>
        public double dPointToJunctionGlobal { get; set; }
        /// <summary>
        /// global distance from meshpoints to contour segments
        /// </summary>
        public double dPointToSegmentGlobal { get; set; }

        /// <summary>
        /// parameterless constructor for subclasses
        /// </summary>
        public MeshingAlgorithm()
        {
        }

        /// <summary>
        /// constructor for loaded mesh
        /// </summary>
        /// <param name="geometryFileData"></param>
        public MeshingAlgorithm(GeometryFileData geometryFileData, double? dPointToPointGlobal = null)
        {
            if (dPointToPointGlobal != null)
                this.dPointToPointGlobal = dPointToPointGlobal ?? 0;

            // create contour junctions and contour segments
            (var regions, var additionalContourJunctions, var unit) = geometryFileData.GetRegionsAndAdditionalPoints<R>();
            contourJunctions = new List<ContourJunction>();
            contourSegments = new List<ContourSegment>();
            this.regions = regions;
            foreach (var r in regions)
            {
                foreach (var point in r.orderedPoints)
                    if (!contourJunctions.Any(p => p.index == point.index))
                        contourJunctions.Add(point);
                foreach (var segment in r.orderedSegments)
                    if (!contourSegments.Any(s => s.index == segment.index))
                        contourSegments.Add(segment);
            }
            foreach (var additionalCJ in additionalContourJunctions)
                contourJunctions.Add(additionalCJ);
            contourJunctions = contourJunctions.OrderBy(cj => cj.index).ToList();
            contourSegments = contourSegments.OrderBy(cs => cs.index).ToList();

            // boundaries
            switch (geometryFileData.dimension)
            {
                case 1:
                    SetOuterContour1D();
                    break;

                case 2:
                    SetOuterContour2D();
                    break;

                case 3:
                    break;
            }

            dimension = geometryFileData.dimension;
        }

        /// <summary>
        /// returns a list with all region indexes where a certain postion is in
        /// </summary>
        /// <param name="position">tested position</param>
        /// <returns></returns>
        public List<R> GetEnclosingRegions(Position position)
        {
            List<R> enclosingRegions = new List<R>();

            switch (dimension)
            {
                case 1:
                    for (int i = 0; i < regions.Count; i++)
                        if ((position.x >= regions[i].orderedPoints.Min(p => p.position.x)) && (position.x <= regions[i].orderedPoints.Max(p => p.position.x)))
                            enclosingRegions.Add(regions[i]);
                    break;
                case 2:
                    for (int i = 0; i < regions.Count; i++)
                        if (position.InPolygon(regions[i], true))
                            enclosingRegions.Add(regions[i]);
                    break;
            }
            return enclosingRegions;
        }

        /// <summary>
        /// calculates the outer region from all regions
        /// </summary>
        /// <param name="regions">all regions, where the outer region is created from</param>
        /// <returns></returns>
        public void SetOuterContour1D()
        {
            outerContour = new Region();
            outerContour.Initialize(-1, new List<ContourSegment> { new ContourSegment(-1, contourJunctions.MinBy(cj => cj.position.x).First(), contourJunctions.MaxBy(cj => cj.position.x).First()) }, pointType.none);

            minMaxValues = ((contourJunctions.Min(cj => cj.position.x), contourJunctions.Max(cj => cj.position.x)), (double.NaN, double.NaN));

            totalSize = minMaxValues.x.max - minMaxValues.x.min;

            defaultTolerance = Math.Max(minMaxValues.x.max - minMaxValues.x.min, minMaxValues.y.max - minMaxValues.y.min) / 1e9;
        }

        /// <summary>
        /// calculates the outer region from all regions
        /// </summary>
        /// <param name="regions">all regions, where the outer region is created from</param>
        /// <returns></returns>
        public void SetOuterContour2D()
        {
            List<ContourJunction> outerGeometryPoints = new List<ContourJunction>();
            List<ContourSegment> outerGeometrySegments = new List<ContourSegment>();

            ContourJunction startPoint = contourJunctions.OrderBy(p => p.position.x).ThenBy(p => p.position.y).First();

            ContourJunction beforePivotPoint = startPoint;
            outerGeometryPoints.Add(beforePivotPoint);

            ContourJunction pivotPoint = beforePivotPoint.neighbors.MinBy(p => Misc.Atan3(startPoint.position.x, startPoint.position.y, -1, 0, p.position.x, p.position.y)).First();

            while (startPoint.index != pivotPoint.index)
            {
                outerGeometryPoints.Add(pivotPoint);
                outerGeometrySegments.Add(contourSegments.Where(s => (s.firstAdjacentContourJunction.index == pivotPoint.index && s.secondAdjacentContourJunction.index == beforePivotPoint.index)
                    || (s.firstAdjacentContourJunction.index == beforePivotPoint.index && s.secondAdjacentContourJunction.index == pivotPoint.index)).First());

                // move pivot point to next point
                pivotPoint = pivotPoint.neighbors.Where(n => n.index != beforePivotPoint.index).MinBy(p =>
                    Misc.Atan3(pivotPoint.position.x, pivotPoint.position.y, beforePivotPoint.position.x, beforePivotPoint.position.y, p.position.x, p.position.y)).First();

                // move before point to previous pivot point
                beforePivotPoint = outerGeometryPoints.Last();
            }
            outerGeometrySegments.Add(contourSegments.Where(s => (s.firstAdjacentContourJunction.index == pivotPoint.index && s.secondAdjacentContourJunction.index == beforePivotPoint.index)
                        || (s.firstAdjacentContourJunction.index == beforePivotPoint.index && s.secondAdjacentContourJunction.index == pivotPoint.index)).First());

            outerContour = new Region();
            outerContour.Initialize(-1, outerGeometrySegments, pointType.none);

            minMaxValues = ((contourJunctions.Min(cj => cj.position.x), contourJunctions.Max(cj => cj.position.x)),
                (contourJunctions.Min(cj => cj.position.y), contourJunctions.Max(cj => cj.position.y)));

            totalSize = outerContour.GetArea() - regions.Where(r => r.type == pointType.hole).Sum(r => r.GetArea());

            defaultTolerance = Math.Max(minMaxValues.x.max - minMaxValues.x.min, minMaxValues.y.max - minMaxValues.y.min) / 1e9;
        }
    }
}