using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BasicLib
{
    /// <summary>
    /// an optical spectrum 
    /// </summary>
    public struct Spectrum
    {
        /// <summary>
        /// spectrum (wavelength in [m], wavelength-width of this intervall in [m], spectral intensity density in [ W / (m^2 * m) ])
        /// <br>lambda = center of discrete wavelength bar [m]</br>
        /// <br>deltaLambda = width of discrete wavelength bar [m]</br>
        /// <br>spectralIntensityDensity = average powerdensity in this wavelength bar [W/(s*m^3)]</br>
        /// <br>-> intensity in single wavelength bar = spectralIntensityDensity * deltaLambda [W/(s*m^2)]</br>
        /// </summary>
        public (double lambda, double deltaLambda, double spectralIntensityDensity)[] data { get; private set; }

        /// <summary>
        /// total summed power over entire spectrum in [W / m^2]
        /// </summary>
        public double totalPower { get; private set; }

        public Spectrum((double lambda, double deltaLambda, double spectralIntensityDensity)[] data)
        {
            this.data = data;
            totalPower = data.Sum(s => s.deltaLambda * s.spectralIntensityDensity);
        }

        /// <summary>
        /// returns the interpolated spectral intensity at a given wavelength of this spectrum
        /// </summary>
        /// <param name="wavelength">requested wavelength in m</param>
        /// <returns></returns>
        public (int index, double spectralIntensityDensity) SpectralIntensityDensityAtWavelength(double wavelength)
        {
            // find wavelength via bisection
            int minIndex = 0;
            int maxIndex = data.Length - 1;
            int meanIndex = maxIndex / 2;

            while (minIndex + 1 != maxIndex)
            {
                meanIndex = (maxIndex + minIndex) / 2;

                if (wavelength < data[meanIndex].lambda)
                    maxIndex = meanIndex;
                else if (wavelength > data[meanIndex].lambda)
                    minIndex = meanIndex;
                else
                    return (meanIndex, data[meanIndex].spectralIntensityDensity);
            }

            // interpolate between values
            double distanceToMin = Math.Abs(data[minIndex].lambda - wavelength);
            double distanceToMax = Math.Abs(data[maxIndex].lambda - wavelength);
            double distanceTotal = distanceToMin + distanceToMax;
            return (minIndex, (distanceToMin * data[maxIndex].spectralIntensityDensity + distanceToMax * data[minIndex].spectralIntensityDensity) / distanceTotal);
        }
    }
}