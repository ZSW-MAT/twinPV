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
    public class Meshing2D_Quadtree<T, R> : MeshingAlgorithm<T, R>
        where T : FiniteElement, new()
        where R : Region, new()
    {
        /// <summary>
        /// global x distance between two finite elements
        /// </summary>
        public double dPointToPointGlobalX { get; set; }
        /// <summary>
        /// global y distance between two finite elements
        /// </summary>
        public double dPointToPointGlobalY { get; set; }
        /// <summary>
        /// amount of finite elements in x direction
        /// </summary>
        public int amountElementsX { get; set; }
        /// <summary>
        /// amount of finite elements in y direction
        /// </summary>
        public int amountElementsY { get; set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Creates a regular Quadtree mesh for rectangular geometry. can only have rectangular outer region!!!
        /// </summary>
        public Meshing2D_Quadtree(out Mesh<T> mesh, GeometryFileData2D geometryFileData2D, int desiredAmountOfPoints) : base(geometryFileData2D)
        {
            // calculate distances
            dPointToPointGlobal = Math.Sqrt(totalSize / desiredAmountOfPoints);

            // amount of points in both directions
            amountElementsX = Convert.ToInt32(Math.Round((minMaxValues.x.max - minMaxValues.x.min) / dPointToPointGlobal));
            amountElementsY = Convert.ToInt32(Math.Round((minMaxValues.y.max - minMaxValues.y.min) / dPointToPointGlobal));

            // distance of points in both directions
            dPointToPointGlobalX = (minMaxValues.x.max - minMaxValues.x.min) / amountElementsX;
            dPointToPointGlobalY = (minMaxValues.y.max - minMaxValues.y.min) / amountElementsY;

            dPointToSegmentGlobal = dPointToPointGlobal / 2;
            dPointToJunctionGlobal = Math.Sqrt(2) * dPointToSegmentGlobal;

            // initialize mesh
            mesh = new Mesh<T>();

            // create finite elements and add to mesh
            double FEsize = dPointToPointGlobalX * dPointToPointGlobalY;
            for (int y = 0; y < amountElementsY; y++)
                for (int x = 0; x < amountElementsX; x++)
                {
                    T newElement = new T();
                    newElement.index = mesh.nextAvailableFiniteElementIndex;
                    newElement.position = new Position((x + 0.5) * dPointToPointGlobalX, (y + 0.5) * dPointToPointGlobalY);
                    newElement.size = FEsize;

                    // neighbors
                    bool isBottomBoundaryElement = y == 0;
                    bool isRightBoundaryElement = x == amountElementsX - 1;
                    bool isTopBoundaryElement = y == amountElementsY - 1;
                    bool isLeftBoundaryElement = x == 0;
                    if (!isBottomBoundaryElement)
                        newElement.neighbors.Add((newElement.index - amountElementsX, dPointToPointGlobalY, dPointToPointGlobalX));
                    if (!isLeftBoundaryElement)
                        newElement.neighbors.Add((newElement.index - 1, dPointToPointGlobalX, dPointToPointGlobalY));
                    if (!isRightBoundaryElement)
                        newElement.neighbors.Add((newElement.index + 1, dPointToPointGlobalX, dPointToPointGlobalY));
                    if (!isTopBoundaryElement)
                        newElement.neighbors.Add((newElement.index + amountElementsX, dPointToPointGlobalY, dPointToPointGlobalX));

                    // segment created from
                    if (isBottomBoundaryElement)
                    {
                        newElement.indexOfBorderElementCreatedFrom = 0;
                        newElement.borderEdgeSize = dPointToPointGlobalX;
                    }
                    if (isRightBoundaryElement)
                    {
                        newElement.indexOfBorderElementCreatedFrom = 1;
                        newElement.borderEdgeSize = dPointToPointGlobalY;
                    }
                    if (isTopBoundaryElement)
                    {
                        newElement.indexOfBorderElementCreatedFrom = 2;
                        newElement.borderEdgeSize = dPointToPointGlobalX;
                    }
                    if (isLeftBoundaryElement)
                    {
                        newElement.indexOfBorderElementCreatedFrom = 3;
                        newElement.borderEdgeSize = dPointToPointGlobalY;
                    }

                    // corners
                    Position toBottomLeft = new Position(-dPointToPointGlobalX / 2, -dPointToPointGlobalY / 2);
                    Position toBottomRight = new Position(dPointToPointGlobalX / 2, -dPointToPointGlobalY / 2);
                    Position toTopRight = new Position(dPointToPointGlobalX / 2, dPointToPointGlobalY / 2);
                    Position toTopLeft = new Position(-dPointToPointGlobalX / 2, dPointToPointGlobalY / 2);
                    newElement.corners.Add((newElement.position + toBottomLeft, isBottomBoundaryElement || isLeftBoundaryElement));
                    newElement.corners.Add((newElement.position + toBottomRight, isBottomBoundaryElement || isRightBoundaryElement));
                    newElement.corners.Add((newElement.position + toTopRight, isTopBoundaryElement || isLeftBoundaryElement));
                    newElement.corners.Add((newElement.position + toTopLeft, isTopBoundaryElement || isRightBoundaryElement));

                    // add element to mesh
                    mesh.AddPoint(newElement);
                }
        }
    }
}