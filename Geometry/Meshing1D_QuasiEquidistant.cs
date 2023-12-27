using BasicLib;
using MoreLinq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Geometry
{
    public class Meshing1D_QuasiEquidistant<T, R> : MeshingAlgorithm<T, R>
        where T : FiniteElement, new()
        where R : Region, new()
    {
        /// <summary>
        /// returns a 1D mesh with equidistant finite elements
        /// </summary>
        /// <param name="mesh">returned mesh with all finite elements</param>
        /// <param name="thicknessOfLayers">thickness of sinlge layers</param>
        /// <param name="desiredAmountOfPoints">approximate amount of elements for the mesh</param>
        public Meshing1D_QuasiEquidistant(out Mesh<T> mesh, GeometryFileData1D geometryFileData1D, int desiredAmountOfPoints) : base(geometryFileData1D)
        {
            // calculate distances
            dPointToPointGlobal = totalSize / desiredAmountOfPoints;
            dPointToJunctionGlobal = dPointToPointGlobal / 50;

            // initialize mesh
            mesh = new Mesh<T>();

            // calculate position of finite elements positions as plumb points in contour segments
            foreach (var cj in contourJunctions)
                cj.SetDefaultPositionTriplets(dPointToJunctionGlobal, 1);
            foreach (var cs in contourSegments)
            {
                cs.SetDefaultPositionTriplets(dPointToPointGlobal, dPointToJunctionGlobal, 1);
                cs.pointTriplets = cs.pointTriplets.OrderBy(pt => pt.plumbPos.x).ToList();
            }

            // create finite elements and add to mesh
            int totalAmountOfPoints = contourSegments.Sum(cs => cs.pointTriplets.Count);
            for (int CSindex = 0; CSindex < contourSegments.Count; CSindex++)
                for (int FEindex = 0; FEindex < contourSegments[CSindex].pointTriplets.Count; FEindex++)
                {
                    T newElement = new T();
                    newElement.index = mesh.nextAvailableFiniteElementIndex;
                    newElement.position = contourSegments[CSindex].pointTriplets[FEindex].plumbPos;

                    if (newElement.index != 0)
                        newElement.neighbors.Add((newElement.index - 1, 1, 1)); // neighbor distances in next for loop
                    if (newElement.index != totalAmountOfPoints - 1)
                        newElement.neighbors.Add((newElement.index + 1, 1, 1)); // neighbor distances in next for loop

                    if (FEindex == 0)
                    {
                        newElement.indexOfBorderElementCreatedFrom = contourSegments[CSindex].firstAdjacentContourJunction.index;
                        newElement.borderEdgeSize = 1;
                    }
                    if (FEindex == contourSegments[CSindex].pointTriplets.Count - 1)
                    {
                        newElement.indexOfBorderElementCreatedFrom = contourSegments[CSindex].secondAdjacentContourJunction.index;
                        newElement.borderEdgeSize = 1;
                    }

                    mesh.AddPoint(newElement);
                }

            // add neighbor distances
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                for (int n = 0; n < mesh.finiteElements[i].neighbors.Count; n++)
                    mesh.finiteElements[i].neighbors[n] = (mesh.finiteElements[i].neighbors[n].index,
                        mesh.finiteElements[i].position.DistanceTo(mesh.finiteElements[mesh.finiteElements[i].neighbors[n].index].position), mesh.finiteElements[i].neighbors[n].edgeSize);

            // add corners for each finite element
            mesh.finiteElements[0].corners.Add((new Position(minMaxValues.x.min), true));
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex - 1; i++)
            {
                Position corner = mesh.finiteElements[i].position.CenterWith(mesh.finiteElements[i + 1].position);

                bool isCornerOnContour = (mesh.finiteElements[i].indexOfBorderElementCreatedFrom != -1) && (mesh.finiteElements[i + 1].indexOfBorderElementCreatedFrom != -1)
                    && (mesh.finiteElements[i].indexOfBorderElementCreatedFrom == mesh.finiteElements[i + 1].indexOfBorderElementCreatedFrom);

                mesh.finiteElements[i].corners.Add((corner, isCornerOnContour));
                mesh.finiteElements[i + 1].corners.Add((corner, isCornerOnContour));
            }
            mesh.finiteElements.Last().Value.corners.Add((new Position(minMaxValues.x.max), false));

            // calculate size for each finite element
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                mesh.finiteElements[i].size = Math.Abs(mesh.finiteElements[i].corners[1].position.x - mesh.finiteElements[i].corners[0].position.x);
        }
    }
}