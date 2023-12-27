using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Geometry;
using Database;
using BasicLib;

namespace Semiconductor
{
    public class RegionSemiconductor : Region
    {
        /// <summary>
        /// material of this semiconductor region
        /// </summary>
        public Material material { get; private set; }

        public bool isGradedEgap { get; set; } = false;

        public (double x, double y)[] gradingCoordinates;

        public bool isAbsorber {get; set;}

        public double roughnessOnTop { get; set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        public RegionSemiconductor()
        {
        }

        // Set properties for semiconductor region ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// sets all properties, a semiconductor region needs to have
        /// </summary>
        public override void SetProperties(double[] preferencesArray, pointType regionType)
        {
            material = Data.GetMaterialFromID((int)Math.Round(preferencesArray[0]));

            gradingCoordinates = new (double x, double y)[5];
            gradingCoordinates[0] = ((preferencesArray[3], preferencesArray[4]));
            gradingCoordinates[1] = ((preferencesArray[5], preferencesArray[6]));
            gradingCoordinates[2] = ((preferencesArray[7], preferencesArray[8]));
            gradingCoordinates[3] = ((preferencesArray[9], preferencesArray[10]));
            gradingCoordinates[4] = ((preferencesArray[11], preferencesArray[12]));

            if (double.IsNaN( gradingCoordinates[0].x))
                isGradedEgap = true;
            else
                isGradedEgap = false;

            isAbsorber = (int)Math.Round(preferencesArray[1]) == 1;

            roughnessOnTop = preferencesArray[2];

            type = regionType;
        }
    }
}