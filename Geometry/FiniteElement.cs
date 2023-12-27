using BasicLib;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace Geometry
{
    /// <summary>
    /// Class for meshpoints
    /// </summary>
    [DataContract]
    public class FiniteElement
    {
        /// <summary>
        /// index of the meshpoint
        /// </summary>
        [DataMember]
        public int index { get; set; }
        /// <summary>
        /// position of the meshpoint
        /// </summary>
        [DataMember]
        public Position position { get; set; }
        /// <summary>
        /// length / area / volume of this finite element in m / m^2 / m^3
        /// </summary>
        [DataMember]
        public double size { get; set; }
        /// <summary>
        /// list of neighbors (index, distance, edgesize) there is also a distance method in this.position.DistanceTo(AnotherPosition)
        /// </summary>
        [DataMember]
        public List<(int index, double distance, double edgeSize)> neighbors { get; set; } = new List<(int index, double distance, double edgeSize)>();
        /// <summary>
        /// list of corner positions (order counterclockwise)
        /// </summary>
        [DataMember]
        public List<(Position position, bool isOnContour)> corners { get; set; } = new List<(Position position, bool isOnContour)>();
        /// <summary>
        /// index contour segment, where this point was created from. if its not a segment point, this variable is -1
        /// </summary>
        [DataMember]
        public int indexOfBorderElementCreatedFrom { get; set; } = -1;
        /// <summary>
        /// edge length to to segment, where this point was created from (if it was not created from a segment, this variable is NaN)
        /// </summary>
        [DataMember]
        public double borderEdgeSize { get; set; } = 0;

        // constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// parameterless constructor
        /// </summary>
        public FiniteElement()
        {
        }

        // set i-th differential equation variable ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// sets the i-th differential equation variable (dependent on type of point)
        /// </summary>
        /// <param name="differentialVariable">value which will be set to the variable of the point</param>
        /// <param name="index">index of the differential variable, which will be set</param>
        public virtual void SetDifferentialEquationVariable(double differentialVariable, int index)
        {
            Console.WriteLine("Override the function \"" + System.Reflection.MethodBase.GetCurrentMethod().Name
                + "()\" in the class \"" + this.GetType().Name + "\"!");
        }

        // get all differential equation variables ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// gets all differential equation variables listed in a double array (length dependent on type of point)
        /// </summary>
        /// <returns></returns>
        public virtual double[] GetDifferentialEquationVariables()
        {
            Console.WriteLine("Override the function \"" + System.Reflection.MethodBase.GetCurrentMethod().Name
                + "()\" in the class \"" + this.GetType().Name + "\"!");
            return new double[] { 1 };
        }

        // get index in neighbor list ███████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// get index in neighbor list
        /// </summary>
        /// <param name="globalElementIndex">index of the element searching for</param>
        public int GetIndexInNeighborList(int globalElementIndex)
        {
            for (int i = 0; i < neighbors.Count; i++)
                if (globalElementIndex == neighbors[i].index)
                    return i;

            return -1;
        }
    }
}