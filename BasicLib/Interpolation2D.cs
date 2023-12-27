using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    public class Interpolation2D
    {
        public (double x, double y, double z)[] basePointsOfArea { get; private set; }

        public Interpolation2D((double x, double y, double z)[] basePointsOfArea)
        {
            this.basePointsOfArea = basePointsOfArea;
        }

        public double GetValueAt(double x, double y, int amountOfClosestPoints = 4)
        {
            (int position, double distance)[] distances = new (int position, double distance)[basePointsOfArea.Length];
            for (int i = 0; i < distances.Length; i++)
                distances[i] = (i, Math.Pow(x - basePointsOfArea[i].x, 2) + Math.Pow(y - basePointsOfArea[i].y, 2));

            distances = distances.OrderBy(d => d.distance).ToArray();

            double totalDistance = 0;
            double totalValue = 0;
            for (int i = 0; i < amountOfClosestPoints; i++)
            {
                totalDistance += distances[i].distance;
                totalValue += basePointsOfArea[distances[i].position].z * distances[i].distance;
            }
            totalValue /= totalDistance;

            return totalValue;
        }
    }
}
