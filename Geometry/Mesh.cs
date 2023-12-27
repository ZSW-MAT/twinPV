using BasicLib;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace Geometry
{
    /// <summary>
    /// Struct which is a container for the mesh including all finite elements and all triangles and both their amount
    /// </summary>
    [DataContract]
    public class Mesh<T>
        where T : FiniteElement
    {
        /// <summary>
        /// dictionary of all meshpoints
        /// </summary>
        [DataMember]
        public Dictionary<int, T> finiteElements { get; private set; } = new Dictionary<int, T>();
        /// <summary>
        /// next available index of for simulation points
        /// </summary>
        [DataMember]
        public int nextAvailableFiniteElementIndex { get; private set; }

        // constructors █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// parameterless constructor
        /// </summary>
        public Mesh()
        {
        }
        /// <summary>
        /// Set the processed data which contains all points and traingles
        /// </summary>
        /// <param name="points">dictionary of current meshpoints</param>
        /// <param name="triangles">dictionary of current triangles</param>
        /// <param name="triangleQuadrants">quadrant array with elements lists</param>
        /// <param name="highestPointIndex">highest index of all simulation points</param>
        /// <param name="highestTriangleIndex">highest index of all simulation triangles</param>
        public Mesh(Dictionary<int, T> points, int highestPointIndex)
        {
            this.finiteElements = points;
            this.nextAvailableFiniteElementIndex = highestPointIndex;
        }

        // Add new point to mesh ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Add new point to mesh
        /// </summary>
        /// <param name="point">point to be added</param>
        public void AddPoint(T point)
        {
            finiteElements.Add(nextAvailableFiniteElementIndex, point);
            nextAvailableFiniteElementIndex++;
        }

        // Remove point from mesh ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Remove point from mesh
        /// </summary>
        public void RemovePoint(int index)
        {
            finiteElements.Remove(index);
        }

        // Calculates the mean distance between neighbor elements ███████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculate the mean distance between neighbor elements
        /// </summary>
        public double CalculatMeanElementDistance()
        {
            double sumDistances = 0;
            double amount = 0;

            foreach (var element in finiteElements.Values)
                foreach (var n in element.neighbors)
                {
                    sumDistances += n.distance;
                    amount++;
                }

            return sumDistances / amount;
        }
    }
}