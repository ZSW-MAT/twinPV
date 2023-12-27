using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Geometry
{
    /// <summary>
    /// Class for an edge of a triangle incling startpoint and endpoints
    /// </summary>
    public class Edge
    {
        /// <summary>
        /// startpoint of this edge
        /// </summary>
        public FiniteElement startpoint { get; private set; }

        /// <summary>
        /// endpoint of this edge
        /// </summary>
        public FiniteElement endpoint { get; private set; }

        /// <summary>
        /// give the index of the triangle, which borders this edge and is outside the cut-out-polygon
        /// </summary>
        public int indexTriangleOutside { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="startpoint">starting point of the edge</param>
        /// <param name="endpoint">endpoint of the edge</param>
        /// <param name="indexTriangleOutside">index of the triangle, which borders this edge and is outside the cut-out-polygon</param>
        public Edge(FiniteElement startpoint, FiniteElement endpoint, int indexTriangleOutside)
        {
            this.startpoint = startpoint;
            this.endpoint = endpoint;
            this.indexTriangleOutside = indexTriangleOutside;
        }

        // check, if edge is the same as other edge █████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// check, if edge is the same as other edge (points are allowed to be switched)
        /// </summary>
        /// <param name="secondEdge">compared edge</param>
        /// <returns></returns>
        public bool Equals(Edge secondEdge)
        {
            var samePoints = startpoint == secondEdge.startpoint && endpoint == secondEdge.endpoint;
            var samePointsReversed = startpoint == secondEdge.endpoint && endpoint == secondEdge.startpoint;
            return samePoints || samePointsReversed;
        }

        // return list of adjacent egdes ████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the list of adjacent egdes
        /// </summary>
        /// <param name="edges">list of potential neighbors</param>
        /// <returns></returns>
        public List<Edge> AdjacentEdges(List<Edge> edges)
        {
            List<Edge> AdjacentEdges = new List<Edge>();
            foreach (Edge edge in edges)
            {
                if (startpoint.index == edge.startpoint.index
                    || endpoint.index == edge.startpoint.index
                    || startpoint.index == edge.endpoint.index
                    || endpoint.index == edge.endpoint.index)
                    AdjacentEdges.Add(edge);
            }
            return AdjacentEdges;
        }

        // write to console █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// write to console
        /// </summary>
        public void Print()
        {
            Console.WriteLine("Edge");
            Console.Write("Starting point: ");
            startpoint.position.Print();
            Console.Write("Endpoint: ");
            endpoint.position.Print();
            Console.Write("index of outside triangle: ");
            Console.WriteLine(indexTriangleOutside);
            Console.WriteLine();
        }
    }
}