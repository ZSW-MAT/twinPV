using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Extreme.Mathematics;
using Database;
using BasicLib;

namespace TransferMatrix
{
    public class LayerTMM
    {
        /// <summary>
        /// optical material of this layer
        /// </summary>
        public Material material { get; set; }
        /// <summary>
        /// position in the whole stack, at which this layer starts
        /// </summary>
        public double layerPositionStart { get; private set; } = 0;
        /// <summary>
        /// position in the whole stack, at which this layer ends
        /// </summary>
        public double layerPositionEnd { get; private set; } = 0;
        /// <summary>
        /// lambda-dependent vector of the complex perpendicular (senkrecht) electric field at the beginning of this layer
        /// </summary>
        public Vector<Complex<double>>[] E_start_s { get; set; }
        /// <summary>
        /// lambda-dependent vector of the complex perpendicular (senkrecht) electric field at the end of this layer
        /// </summary>
        public Vector<Complex<double>>[] E_end_s { get; set; }
        /// <summary>
        /// lambda-dependent vector of the complex parallel electric field at the beginning of this layer
        /// </summary>
        public Vector<Complex<double>>[] E_start_p { get; set; }
        /// <summary>
        /// lambda-dependent vector of the complex parallel electric field at the end of this layer
        /// </summary>
        public Vector<Complex<double>>[] E_end_p { get; set; }

        /// <summary>
        /// lambda-dependent incident angle in rad in this layer within a specific layerstack
        /// </summary>
        public Complex<double>[] angleOfIncidence { get; set; }

        /// <summary>
        /// Constuctor
        /// </summary>
        /// <param name="material">material, which this layer is made of</param>
        /// <param name="layerStart">position in the whole stack, at which this layer starts</param>
        /// <param name="layerEnd">position in the whole stack, at which this layer ends</param>
        public LayerTMM(Material material, double layerStart, double layerEnd, int lenghtOfSpectrum)
        {
            this.material = material;
            layerPositionStart = layerStart;
            layerPositionEnd = layerEnd;
            E_start_s = new Vector<Complex<double>>[lenghtOfSpectrum];
            E_end_s = new Vector<Complex<double>>[lenghtOfSpectrum];
            E_start_p = new Vector<Complex<double>>[lenghtOfSpectrum];
            E_end_p = new Vector<Complex<double>>[lenghtOfSpectrum];
            angleOfIncidence = new Complex<double>[lenghtOfSpectrum];
        }

        /// <summary>
        /// Get the half-matrix D for this layer for perpendicular polarization
        /// </summary>
        /// <returns></returns>
        public Matrix<Complex<double>> SplitDiffractionMatrix_s(int specIndex)
        {
            Complex<double> n = material.propertiesOptics.n_toSpectrum[specIndex];
            return Matrix.Create(new Complex<double>[,] { { new Complex<double>(1, 0), new Complex<double>(1, 0) },
                { Complex<double>.Cos(angleOfIncidence[specIndex]) * n, -Complex<double>.Cos(angleOfIncidence[specIndex]) * n } });
        }
        /// <summary>
        /// Get the modified matrix for D for roughness effects (changes the original values for r and t in diffraction matrix to r' and t')
        /// </summary>
        /// <returns></returns>
        public Matrix<Complex<double>> ModifyingRoughnessMatrix_s(int specIndex, LayerTMM layerAfterThisLayer, double roughness, Spectrum spectrum)
        {
            // get constant values
            double n_own = material.propertiesOptics.n_toSpectrum[specIndex].Re;
            Complex<double> cosTheta_own = Complex<double>.Cos(angleOfIncidence[specIndex]);

            double n_after = layerAfterThisLayer.material.propertiesOptics.n_toSpectrum[specIndex].Re;
            Complex<double> cosTheta_after = Complex<double>.Cos(layerAfterThisLayer.angleOfIncidence[specIndex]);

            // get original r
            Complex<double> r = (n_own * cosTheta_own - n_after * cosTheta_after) / (n_own * cosTheta_own + n_after * cosTheta_after);

            // modify r (r can only be modified with a matrix) 
            // https://www.wolframalpha.com/input/?i=%7B%7B1%2Cf*r%7D%2C%7Bf*r%2C1%7D%7D+*+inverse%28%7B%7B1%2Cr%7D%2C%7Br%2C1%7D%7D%29
            Complex<double> modifyingR = Complex<double>.Exp(-2 * Math.Pow(2 * Math.PI * roughness / spectrum.data[specIndex].lambda, 2) * Math.Pow(n_own, 2));
            Complex<double> a = (modifyingR * r * r - 1) / (r * r - 1);
            Complex<double> c = (r - modifyingR * r) / (r * r - 1);

            // modify t (can be modified by factor)
            Complex<double> modifyingT = Complex<double>.Exp(-0.5 * Math.Pow(2 * Math.PI * roughness / spectrum.data[specIndex].lambda, 2) * Math.Pow(n_own - n_after, 2));

            return 1 / modifyingT * Matrix.Create(new Complex<double>[,] { { a, c }, { c, a } });
        }

        /// <summary>
        /// Get the half-matrix D for this layer for parallel polarization
        /// </summary>
        /// <returns></returns>
        public Matrix<Complex<double>> SplitDiffractionMatrix_p(int specIndex)
        {
            Complex<double> n = material.propertiesOptics.n_toSpectrum[specIndex];
            return Matrix.Create(new Complex<double>[,] { { Complex<double>.Cos(angleOfIncidence[specIndex]),
                    Complex<double>.Cos(angleOfIncidence[specIndex]) }, { n, -n } });
        }
        /// <summary>
        /// Get the modified matrix for D for roughness effects (changes the original values for r and t in diffraction matrix to r' and t')
        /// </summary>
        /// <returns></returns>
        public Matrix<Complex<double>> ModifyingRoughnessMatrix_p(int specIndex, LayerTMM layerAfterThisLayer, double roughness, Spectrum spectrum)
        {
            // get constant values
            double n_own = material.propertiesOptics.n_toSpectrum[specIndex].Re;
            Complex<double> cosTheta_own = Complex<double>.Cos(angleOfIncidence[specIndex]);

            double n_after = layerAfterThisLayer.material.propertiesOptics.n_toSpectrum[specIndex].Re;
            Complex<double> cosTheta_after = Complex<double>.Cos(layerAfterThisLayer.angleOfIncidence[specIndex]);

            // get original r
            Complex<double> r = (n_after * cosTheta_own - n_own * cosTheta_after) / (n_after * cosTheta_own + n_own * cosTheta_after);

            // modify r (r can only be modified with a matrix)
            // https://www.wolframalpha.com/input/?i=%7B%7B1%2Cf*r%7D%2C%7Bf*r%2C1%7D%7D+*+inverse%28%7B%7B1%2Cr%7D%2C%7Br%2C1%7D%7D%29
            Complex<double> modifyingR = Complex<double>.Exp(-2 * Math.Pow(2 * Math.PI * roughness / spectrum.data[specIndex].lambda, 2) * Math.Pow(n_own, 2));
            Complex<double> a = (modifyingR * r * r - 1) / (r * r - 1);
            Complex<double> c = (r - modifyingR * r) / (r * r - 1);

            // modify t (can be modified by factor)
            Complex<double> modifyingT = Complex<double>.Exp(-0.5 * Math.Pow(2 * Math.PI * roughness / spectrum.data[specIndex].lambda, 2) * Math.Pow(n_own - n_after, 2));

            return 1 / modifyingT * Matrix.Create(new Complex<double>[,] { { a, c }, { c, a } });
        }

        /// <summary>
        /// Get the Matrix P for a wave propagating in this layer: (E_toRight(x) , E_toLeft(x)) = P(x) * electricFieldEnd
        /// if position is not given, then P reaches from end to start, otherwise P reaches from end to the given (aboslute) position within the whole layerstack
        /// </summary>
        /// <param name="lambda">wavelength in m</param>
        /// <param name="position">position in the whole stack, to which the propagation matrix should reach</param>
        /// <returns></returns>
        public Matrix<Complex<double>> PropagationMatrix(double lambda, int specIndex, double position = double.NaN)
        {
            // get position within THIS SINGLE layer (input position is whithin the WHOLE stack)
            double lengthInLayer;
            if (double.IsNaN(position))
                lengthInLayer = layerPositionEnd - layerPositionStart;
            else
                lengthInLayer = layerPositionEnd - position;

            Complex<double> exponent = 2 * Math.PI * lengthInLayer / lambda * Complex<double>.Cos(angleOfIncidence[specIndex]) * material.propertiesOptics.n_toSpectrum[specIndex];

            return Matrix.Create(new Complex<double>[,] { { Complex<double>.Exp(-Complex<double>.I * exponent), 0 },
                { 0, Complex<double>.Exp(Complex<double>.I * exponent) } });
        }
    }
}