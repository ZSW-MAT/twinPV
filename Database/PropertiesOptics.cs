using BasicLib;
using Extreme.Mathematics;
using MoreLinq;
using System;
using System.Windows;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Database
{
    public class PropertiesOptics
    {
        /// <summary>
        /// coefficienct: describes the amount of ligth, which goes through this material (0 = no light passing through, 1 = all light passing through)
        /// </summary>
        public double simpleLightTransmissionCoefficient { get; private set; }

        /// <summary>
        /// lambert beer: prefactor (nonReflectedPower -> I0) (from 0 to 1) and exponent (extinctionCoefficient -> a) in 1/meter in I(z) = I0 * exp(-a * z) which determines how much light goes into the device (0 = all light goes through)
        /// </summary>
        public double lambertBeerAbsorptionCoefficient { get; private set; }

        /// <summary>
        /// lambda-dependent complex refractive index
        /// </summary>
        public (double lambda, Complex<double> n)[] n_rawData { get; private set; }

        /// <summary>
        /// lambda-dependent complex refractive index (index in this array is index in spectrum)
        /// </summary>
        public Complex<double>[] n_toSpectrum { get; private set; }

        public (double weight1, int ID1, double weight2, int ID2) originalMaterials { get; private set; } = (-1, -1, -1, -1);

        /// <summary>
        /// Constructor from filepath
        /// </summary>
        public PropertiesOptics(string filepath)
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
                    simpleLightTransmissionCoefficient = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "simple transmission coefficient:") + 1].Trim());
                    lambertBeerAbsorptionCoefficient = InputOutput.ToDoubleWithArbitrarySeparator(lines[InputOutput.GetLineOfStringInArray(lines, "lambert beer absorption coefficient:") + 1].Trim());
                }
                catch (Exception e)
                {
                    MessageBox.Show("File: " + Path.GetFullPath(filepath) + "\n" + e.Message, "Numeric input exeption", MessageBoxButton.OK, MessageBoxImage.Error);
                }

                //  ██╗ refractive data
                //  ╚═╝
                int startRefractiveData = InputOutput.GetLineOfStringInArray(lines, "refractive data:") + 1;
                if (lines.Skip(startRefractiveData).Where(line => !string.IsNullOrWhiteSpace(line)).Count() < 1)
                {
                    n_rawData = new (double lambda, Complex<double> n)[] { (500e-9, new Complex<double>(1, 0)) };
                }
                else
                {
                    var array = InputOutput.ReadLinesTo2DArray(lines.Skip(startRefractiveData).ToArray());
                    n_rawData = new (double lambda, Complex<double> n)[array.GetLength(0)];
                    for (int i = 0; i < array.GetLength(0); i++)
                        n_rawData[i] = (array[i, 0] * 1e-9, new Complex<double>(array[i, 1], array[i, 2]));
                }
            }
        }

        /// <summary>
        /// set array for each wavelength in spectrum to array (performance reasons)
        /// </summary>
        /// <param name="spectrum">spectrum (lambda in meter, deltaLambda in meter, spectral intensity density in W/m^3)</param>
        public void InitializeNKarray((double lambda, double deltaLambda, double spectralIntensityDensity)[] spectrum)
        {
            if (originalMaterials.ID1 == -1) // normal material
            {
                n_toSpectrum = new Complex<double>[spectrum.Length];
                for (int specIndex = 0; specIndex < spectrum.Length; specIndex++)
                    n_toSpectrum[specIndex] = n_rawData.MinBy(p => Math.Abs(p.lambda - spectrum[specIndex].lambda)).First().n;
            }
            else // effective material
            {
                n_toSpectrum = new Complex<double>[Data.GetMaterialFromID(originalMaterials.ID1).propertiesOptics.n_toSpectrum.Length];
                for (int specIndex = 0; specIndex < Data.GetMaterialFromID(originalMaterials.ID1).propertiesOptics.n_toSpectrum.Length; specIndex++)
                {
                    // Bruggeman approximation
                    /*
                    double w = originalMaterials.weight1;
                    var n1 = Data.opticMaterials[originalMaterials.ID1].n_toSpectrum[specIndex];
                    var n2 = Data.opticMaterials[originalMaterials.ID2].n_toSpectrum[specIndex];

                    n_toSpectrum[specIndex] = -0.5 * Complex<double>.Sqrt(Complex<double>.Sqrt(
                        Complex<double>.Pow(-3 * n1 * n1 * w + n1 * n1 + 3 * n2 * n2 * w - 2 * n2 * n2, 2) + 8 * n1 * n1 * n2 * n2)
                        + 3 * n1 * n1 * w - n1 * n1 - 3 * n2 * n2 * w + 2 * n2 * n2);
                    */

                    // arithmetic mean
                    n_toSpectrum[specIndex] = originalMaterials.weight1 * Data.GetMaterialFromID(originalMaterials.ID1).propertiesOptics.n_toSpectrum[specIndex]
                        + originalMaterials.weight2 * Data.GetMaterialFromID(originalMaterials.ID2).propertiesOptics.n_toSpectrum[specIndex];
                }
            }
        }
    }
}