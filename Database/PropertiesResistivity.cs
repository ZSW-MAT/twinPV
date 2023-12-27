using BasicLib;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace Database
{
    public class PropertiesResistivity
    {
        /// <summary>
        /// electrics: determines the specific electrical resistivity of this material -> rho(z) = rho_Bulk + (rho_Thin - rho_Bulk) * exp(-decay * z),
        /// resistivity in the bulk material [Ohm*m]
        /// </summary>
        public double rhoInBulk { get; private set; }
        /// <summary>
        /// electrics: determines the specific electrical resistivity of this material -> rho(z) = rho_Bulk + (rho_Thin - rho_Bulk) * exp(-decay * z),
        /// resistivity in an imaginary infenitesimal thin layer (higher due to roughness and unever layer growth) [Ohm*m]
        /// </summary>
        public double rhoInInfinitesimalThinLayer { get; private set; }
        /// <summary>
        /// electrics: determines the specific electrical resistivity of this material -> rho(z) = rho_Bulk + (rho_Thin - rho_Bulk) * exp(-decay * z),
        /// determines how fast the resistivity goes from thin to bulk material [1/m]
        /// </summary>
        public double decayCoefficient { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// constructs a new material and sets all properties
        /// </summary>
        public PropertiesResistivity(string filepath)
        {
            if (filepath == null)
            {
                MessageBox.Show("Could not find the file " + Path.GetFullPath(filepath) + ".", "File not found exeption", MessageBoxButton.OK, MessageBoxImage.Error);
            }
            else
            {
                string[] lines = InputOutput.ReadInputFile(filepath);

                try
                {
                    rhoInBulk = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "resistivity bulk:") + 1].Trim());
                    rhoInInfinitesimalThinLayer = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "resistivity infinitesimal thin layer:") + 1].Trim());
                    decayCoefficient = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "decay coefficient:") + 1].Trim());
                }
                catch (Exception e)
                {
                    MessageBox.Show("File: " + Path.GetFullPath(filepath) + "\n" + e.Message, "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                }
            }
        }

        // Get Resistivity for given layer thickness ████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the specific resistivity for a given layer thickness z of this material via rho(z) = rho_Bulk + (rho_Thin - rho_Bulk) * exp(-decay * z)
        /// </summary>
        /// <param name="layerThickness">layer thickness z, for which the specific resistivity is calculated</param>
        public double GetSpecificResistivity(double layerThickness)
        {
            return rhoInBulk + (rhoInInfinitesimalThinLayer - rhoInBulk)
                * Math.Exp(-decayCoefficient * layerThickness);
        }
    }
}