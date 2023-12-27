using BasicLib;
using Extreme.Mathematics;
using Database;
using MoreLinq;
using System;
using System.Linq;

namespace TransferMatrix
{
    public class ModelTMM
    {
        // electric fields
        /// <summary>
        /// lambda-dependent incoming perpendicular (senkrecht) amplitude of electric field before the whole layerstack
        /// </summary>
        public Complex<double>[] E_in_s { get; private set; }
        /// <summary>
        /// lambda-dependent reflected perpendicular (senkrecht) amplitude of electric field of the whole layerstack
        /// </summary>
        public Complex<double>[] E_refl_s { get; private set; }
        /// <summary>
        /// lambda-dependent transmitted perpendicular (senkrecht) amplitude of electric field through the whole layerstack
        /// </summary>
        public Complex<double>[] E_trans_s { get; private set; }
        /// <summary>
        /// lambda-dependent incoming parallel amplitude of electric field before the whole layerstack
        /// </summary>
        public Complex<double>[] E_in_p { get; private set; }
        /// <summary>
        /// lambda-dependent reflected parallel amplitude of electric field of the whole layerstack
        /// </summary>
        public Complex<double>[] E_refl_p { get; private set; }
        /// <summary>
        /// lambda-dependent transmitted parallel amplitude of electric field through the whole layerstack
        /// </summary>
        public Complex<double>[] E_trans_p { get; private set; }

        // R and T
        /// <summary>
        /// lambda-dependent reflected fraction of perpendicular (senkrecht) intensity of the whole layerstack
        /// </summary>
        public double[] R_s { get; private set; }
        /// <summary>
        /// lambda-dependent transmitted fraction of perpendicular (senkrecht) intensity through the whole layerstack
        /// </summary>
        public double[] T_s { get; private set; }
        /// <summary>
        /// lambda-dependent reflected fraction of parallel intensity of the whole layerstack
        /// </summary>
        public double[] R_p { get; private set; }
        /// <summary>
        /// lambda-dependent transmitted fraction of parallel intensity through the whole layerstack
        /// </summary>
        public double[] T_p { get; private set; }
        /// <summary>
        /// lambda-dependent reflected fraction of total intensity of the whole layerstack
        /// </summary>
        public double[] R { get; private set; }
        /// <summary>
        /// lambda-dependent transmitted fraction of total intensity through the whole layerstack
        /// </summary>
        public double[] T { get; private set; }

        /// <summary>
        /// /// fraction of parallel polarized light with respect to the total (s and p) polarized intensity
        /// </summary>
        double fractionInParallelPolarized { get; set; }

        // layer stack
        /// <summary>
        /// stack of all layers
        /// </summary>
        public LayerTMM[] layerStack { get; private set; }
        /// <summary>
        /// layer without a length before the layerstack (usually air) -> needs to be k=0 (see also https://arxiv.org/pdf/1603.02720.pdf)
        /// </summary>
        public LayerTMM layerBeforeStack { get; private set; }
        /// <summary>
        /// layer without a length behind the layerstack
        /// </summary>
        public LayerTMM layerBehindStack { get; private set; }
        /// <summary>
        /// length of the whole layerstack
        /// </summary>
        public double lengthOfStack { get; private set; }

        // spectrum
        /// <summary>
        /// <br>spectrum, which shines on this optic model</br>
        /// <br>lambda = center of discrete wavelength bar [m]</br>
        /// <br>deltaLambda = width of discrete wavelength bar [m]</br>
        /// <br>spectralIntensityDensity = average powerdensity in this wavelength bar [W/(s*m^3)]</br>
        /// <br>-> intensity in single wavelength bar = spectralIntensityDensity * deltaLambda [W/(s*m^2)]</br>
        /// </summary>
        public Spectrum spectrum { get; private set; }

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="materialBeforeStack">material, which is before the material stack (need to be k=0 in every case!!! see also https://arxiv.org/pdf/1603.02720.pdf)</param>
        /// <param name="materialBehindStack">material, which is behind the material stack</param>
        /// <param name="materialStack">array of materials (from top to bottom) with its thicknesses</param>
        /// <param name="spectrum">spectrum, which shines on this optic model</param>
        /// <param name="angleOfIncidenceInDegree">incident angle (between ray and perpendicular line) in degree</param>
        /// <param name="fractionInParallelPolarized">fraction of parallel polarized light with respect to the total (s and p) polarized intensity</param>
        public ModelTMM(Material materialBeforeStack, (Material material, double roughnessOnTop) materialBehindStack,
            (Material material, double thickness, double roughnessOnTop)[] materialStack, Spectrum spectrum,
            double multiplySpectrumByFactor = 1, double angleOfIncidenceInDegree = 0, double fractionInParallelPolarized = 0.5)
        {
            //  ██╗ 
            //  ╚██╗ initialize arrays
            //  ██╔╝
            //  ╚═╝
            if (multiplySpectrumByFactor == 1)
                this.spectrum = spectrum;
            else
                this.spectrum = new Spectrum(spectrum.data.Select(s => (s.lambda, s.deltaLambda, s.spectralIntensityDensity * multiplySpectrumByFactor)).ToArray());
            this.fractionInParallelPolarized = fractionInParallelPolarized;
            E_in_s = new Complex<double>[spectrum.data.Length];
            E_refl_s = new Complex<double>[spectrum.data.Length];
            E_trans_s = new Complex<double>[spectrum.data.Length];
            E_in_p = new Complex<double>[spectrum.data.Length];
            E_refl_p = new Complex<double>[spectrum.data.Length];
            E_trans_p = new Complex<double>[spectrum.data.Length];
            R_s = new double[spectrum.data.Length];
            T_s = new double[spectrum.data.Length];
            R_p = new double[spectrum.data.Length];
            T_p = new double[spectrum.data.Length];
            R = new double[spectrum.data.Length];
            T = new double[spectrum.data.Length];

            //  ██╗ 
            //  ╚██╗ extract n and k from data (sort by index in spectrum)
            //  ██╔╝
            //  ╚═╝
            // normal materials
            if (materialBeforeStack.propertiesOptics.originalMaterials.ID1 == -1)
                materialBeforeStack.propertiesOptics.InitializeNKarray(spectrum.data);
            if (materialBehindStack.material.propertiesOptics.originalMaterials.ID1 == -1)
                materialBehindStack.material.propertiesOptics.InitializeNKarray(spectrum.data);
            foreach (var material in materialStack.Select(m => m.material).Where(m => m.propertiesOptics.originalMaterials.ID1 == -1))
                material.propertiesOptics.InitializeNKarray(spectrum.data);

            // effective materials
            if (materialBeforeStack.propertiesOptics.originalMaterials.ID1 != -1)
                materialBeforeStack.propertiesOptics.InitializeNKarray(spectrum.data);
            if (materialBehindStack.material.propertiesOptics.originalMaterials.ID1 != -1)
                materialBehindStack.material.propertiesOptics.InitializeNKarray(spectrum.data);
            foreach (var material in materialStack.Select(m => m.material).Where(m => m.propertiesOptics.originalMaterials.ID1 != -1))
                material.propertiesOptics.InitializeNKarray(spectrum.data);

            //  ██╗ 
            //  ╚██╗ create layerstack
            //  ██╔╝
            //  ╚═╝
            layerBeforeStack = new LayerTMM(materialBeforeStack, 0, 0, spectrum.data.Length);
            layerStack = new LayerTMM[materialStack.Length];
            for (int i = 0; i < layerStack.Length; i++)
                layerStack[i] = new LayerTMM(materialStack[i].material, Enumerable.Range(0, i).Sum(k => materialStack[k].thickness), Enumerable.Range(0, i).Sum(k => materialStack[k].thickness) + materialStack[i].thickness, spectrum.data.Length);
            lengthOfStack = layerStack.Length == 0 ? 0 : layerStack.Last().layerPositionEnd;
            layerBehindStack = new LayerTMM(materialBehindStack.material, lengthOfStack, lengthOfStack, spectrum.data.Length);

            //  ██╗ 
            //  ╚██╗ calculate electric fields (TMM)
            //  ██╔╝
            //  ╚═╝
            for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
            {
                //  ██╗ set angles via Snell's law
                //  ╚═╝
                layerBeforeStack.angleOfIncidence[specIndex] = angleOfIncidenceInDegree * Math.PI / 180;
                for (int i = 0; i < layerStack.Length; i++)
                    layerStack[i].angleOfIncidence[specIndex] = Complex<double>.Asin(layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex]
                        / layerStack[i].material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Sin(layerBeforeStack.angleOfIncidence[specIndex]));
                layerBehindStack.angleOfIncidence[specIndex] = Complex<double>.Asin(layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex]
                    / layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Sin(layerBeforeStack.angleOfIncidence[specIndex]));

                //  ██╗ calculate transfer matrix
                //  ╚═╝
                var totalEfieldAmplitude = MiscTMM.GetEfieldAmplitudeFromIntensity(spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda, layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex]);
                E_in_s[specIndex] = totalEfieldAmplitude * Math.Sqrt(1 - fractionInParallelPolarized);
                E_in_p[specIndex] = totalEfieldAmplitude * Math.Sqrt(fractionInParallelPolarized);

                var M_s = MiscTMM.GetInverse2x2(layerBeforeStack.SplitDiffractionMatrix_s(specIndex));
                var M_p = MiscTMM.GetInverse2x2(layerBeforeStack.SplitDiffractionMatrix_p(specIndex));

                for (int i = 0; i < layerStack.Length; i++)
                {
                    M_s *= layerStack[i].SplitDiffractionMatrix_s(specIndex);
                    if (materialStack[i].roughnessOnTop > 0)
                    {
                        if (i == 0)
                            M_s *= layerBeforeStack.ModifyingRoughnessMatrix_s(specIndex, layerStack[i], materialStack[i].roughnessOnTop, spectrum);
                        else
                            M_s *= layerStack[i - 1].ModifyingRoughnessMatrix_s(specIndex, layerStack[i], materialStack[i].roughnessOnTop, spectrum);
                    }
                    M_s *= layerStack[i].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex);
                    M_s *= MiscTMM.GetInverse2x2(layerStack[i].SplitDiffractionMatrix_s(specIndex));

                    M_p *= layerStack[i].SplitDiffractionMatrix_p(specIndex);
                    if (materialStack[i].roughnessOnTop > 0)
                    {
                        if (i == 0)
                            M_p *= layerBeforeStack.ModifyingRoughnessMatrix_p(specIndex, layerStack[i], materialStack[i].roughnessOnTop, spectrum);
                        else
                            M_p *= layerStack[i - 1].ModifyingRoughnessMatrix_p(specIndex, layerStack[i], materialStack[i].roughnessOnTop, spectrum);
                    }
                    M_p *= layerStack[i].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex);
                    M_p *= MiscTMM.GetInverse2x2(layerStack[i].SplitDiffractionMatrix_p(specIndex));
                }

                M_s *= layerBehindStack.SplitDiffractionMatrix_s(specIndex);
                if (materialBehindStack.roughnessOnTop > 0)
                    M_s *= layerStack.Last().ModifyingRoughnessMatrix_s(specIndex, layerBehindStack, materialBehindStack.roughnessOnTop, spectrum);
                M_p *= layerBehindStack.SplitDiffractionMatrix_p(specIndex);
                if (materialBehindStack.roughnessOnTop > 0)
                    M_p *= layerStack.Last().ModifyingRoughnessMatrix_p(specIndex, layerBehindStack, materialBehindStack.roughnessOnTop, spectrum);

                //  ██╗ total reflectance
                //  ╚═╝
                // perpendicular
                Complex<double> r_s = M_s[1, 0] / M_s[0, 0];
                E_refl_s[specIndex] = r_s * E_in_s[specIndex];
                double R_s = r_s.MagnitudeSquared;
                this.R_s[specIndex] = R_s;

                // parallel
                Complex<double> r_p = M_p[1, 0] / M_p[0, 0];
                E_refl_p[specIndex] = r_p * E_in_p[specIndex];
                double R_p = r_p.MagnitudeSquared;
                this.R_p[specIndex] = R_p;

                // total
                R[specIndex] = R_s * (1 - fractionInParallelPolarized) + R_p * fractionInParallelPolarized;

                //  ██╗ total transmittance
                //  ╚═╝
                // perpendicular
                Complex<double> t_s = 1 / M_s[0, 0];
                E_trans_s[specIndex] = t_s * E_in_s[specIndex];
                double T_s = (layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Cos(layerBehindStack.angleOfIncidence[specIndex])).Re
                    / (layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex])).Re
                    * t_s.MagnitudeSquared;
                this.T_s[specIndex] = T_s;

                // parallel
                Complex<double> t_p = 1 / M_p[0, 0];
                E_trans_p[specIndex] = t_p * E_in_p[specIndex];
                double T_p = (layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Cos(layerBehindStack.angleOfIncidence[specIndex]).Conjugate()).Re
                    / (layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex] * Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]).Conjugate()).Re
                    * t_p.MagnitudeSquared;
                this.T_p[specIndex] = T_p;

                // total
                T[specIndex] = T_s * (1 - fractionInParallelPolarized) + T_p * fractionInParallelPolarized;

                //  ██╗ electric fields at start and end of each layer
                //  ╚═╝
                Vector<Complex<double>> E_behind_s = Vector.Create(new Complex<double>[] { E_trans_s[specIndex], 0 });
                Vector<Complex<double>> E_behind_p = Vector.Create(new Complex<double>[] { E_trans_p[specIndex], 0 });

                if (layerStack.Length > 0)
                {
                    // last layer
                    // perpendicular
                    layerStack.Last().E_end_s[specIndex] = MiscTMM.GetInverse2x2(layerStack.Last().SplitDiffractionMatrix_s(specIndex))
                        * layerBehindStack.SplitDiffractionMatrix_s(specIndex) /*here modified D matrix*/ * E_behind_s;
                    layerStack.Last().E_start_s[specIndex] = layerStack.Last().PropagationMatrix(spectrum.data[specIndex].lambda, specIndex) * layerStack.Last().E_end_s[specIndex];
                    // parallel
                    layerStack.Last().E_end_p[specIndex] = MiscTMM.GetInverse2x2(layerStack.Last().SplitDiffractionMatrix_p(specIndex))
                        * layerBehindStack.SplitDiffractionMatrix_p(specIndex) * E_behind_p;
                    layerStack.Last().E_start_p[specIndex] = layerStack.Last().PropagationMatrix(spectrum.data[specIndex].lambda, specIndex) * layerStack.Last().E_end_p[specIndex];

                    // all other layers
                    for (int i = layerStack.Length - 2; i >= 0; i--)
                    {
                        // perpendicular
                        layerStack[i].E_end_s[specIndex] = MiscTMM.GetInverse2x2(layerStack[i].SplitDiffractionMatrix_s(specIndex))
                            * layerStack[i + 1].SplitDiffractionMatrix_s(specIndex) * layerStack[i + 1].E_start_s[specIndex];
                        layerStack[i].E_start_s[specIndex] = layerStack[i].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex) * layerStack[i].E_end_s[specIndex];
                        // parallel
                        layerStack[i].E_end_p[specIndex] = MiscTMM.GetInverse2x2(layerStack[i].SplitDiffractionMatrix_p(specIndex))
                            * layerStack[i + 1].SplitDiffractionMatrix_p(specIndex) * layerStack[i + 1].E_start_p[specIndex];
                        layerStack[i].E_start_p[specIndex] = layerStack[i].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex) * layerStack[i].E_end_p[specIndex];
                    }
                }
            }
        }

        // public accessable methods ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// poynting vector in W/m² at a given depth within the optical stack
        /// </summary>
        /// <param name="position">depth in meter, where the field is calculated</param>
        /// <param name="fromLambda">wavelength in meter, from which the field is added up (photons with shorter wavelength / higher energy are NOT counted) (NaN -> count all photons down to lambda = 0)</param>
        /// <param name="toLambda">wavelength in meter, up to which the photons are counted (photons with longer wavelength / lower energy are NOT counted) (NaN -> count all photons up to lambda = inf)</param>
        public (double s, double p) GetPoyntingVectorAtPosition(double position, double fromLambda = double.NaN, double toLambda = double.NaN)
        {
            int lambdaIndexStart = GetIndexInSpectrumFromLambda(fromLambda, 0);
            int lambdaIndexEnd = GetIndexInSpectrumFromLambda(toLambda, spectrum.data.Length - 1);

            var lambdaDependentPoyntingVector = GetLambdaDependentPoyntingVectorAtPosition(position);

            double Ss = 0;
            double Sp = 0;

            for (int specIndex = lambdaIndexStart; specIndex <= lambdaIndexEnd; specIndex++)
            {
                Ss += lambdaDependentPoyntingVector.s[specIndex];
                Sp += lambdaDependentPoyntingVector.p[specIndex];
            }

            return (Ss, Sp);
        }

        /// <summary>
        /// returns an array of poynting vectors in W/m² with the same wavelengths as the spectrum
        /// </summary>
        /// <param name="position">position in the stack, where the intensity is calculated</param>
        public (double[] s, double[] p) GetLambdaDependentPoyntingVectorAtPosition(double position)
        {
            double[] Ss = new double[spectrum.data.Length];
            double[] Sp = new double[spectrum.data.Length];

            var Efield = GetEfieldAtPosition(position);

            // before layerstack
            if (position < 0)
            {
                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerBeforeStack.material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = cosTheta0;

                    Ss[specIndex] = (n * cosTheta * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).Conjugate()
                        * (Efield.s[specIndex][0] - Efield.s[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);

                    Sp[specIndex] = (n * cosTheta.Conjugate() * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).Conjugate()
                        * (Efield.p[specIndex][0] - Efield.p[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);
                }
            }

            // after layerstack
            else if (position > lengthOfStack)
            {
                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = Complex<double>.Cos(layerBehindStack.angleOfIncidence[specIndex]);

                    Ss[specIndex] = (n * cosTheta * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).Conjugate()
                        * (Efield.s[specIndex][0] - Efield.s[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);

                    Sp[specIndex] = (n * cosTheta.Conjugate() * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).Conjugate()
                        * (Efield.p[specIndex][0] - Efield.p[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);
                }
            }

            // within the layerstack
            else
            {
                int layerIndex = GetLayerIndexAtPosition(position);

                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerStack[layerIndex].material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = Complex<double>.Cos(layerStack[layerIndex].angleOfIncidence[specIndex]);

                    Ss[specIndex] = (n * cosTheta * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).Conjugate()
                        * (Efield.s[specIndex][0] - Efield.s[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);

                    Sp[specIndex] = (n * cosTheta.Conjugate() * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).Conjugate()
                        * (Efield.p[specIndex][0] - Efield.p[specIndex][1])).Re / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c);
                }
            }

            return (Ss, Sp);
        }

        /// <summary>
        /// returns an array of electric fields in V/m with the same wavelengths as the spectrum
        /// </summary>
        /// <param name="position">position in the stack, where the intensity is calculated</param>
        public (Vector<Complex<double>>[] s, Vector<Complex<double>>[] p) GetEfieldAtPosition(double position)
        {
            Vector<Complex<double>>[] Efield_s = new Vector<Complex<double>>[spectrum.data.Length];
            Vector<Complex<double>>[] Efield_p = new Vector<Complex<double>>[spectrum.data.Length];

            // before layerstack
            if (position < 0)
            {
                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Efield_s[specIndex] = layerBeforeStack.PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * Vector.Create(new Complex<double>[] { E_in_s[specIndex], E_refl_s[specIndex] });
                    Efield_p[specIndex] = layerBeforeStack.PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * Vector.Create(new Complex<double>[] { E_in_p[specIndex], E_refl_p[specIndex] });
                }
            }

            // after layerstack
            else if (position > lengthOfStack)
            {
                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Efield_s[specIndex] = layerBehindStack.PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * Vector.Create(new Complex<double>[] { E_trans_s[specIndex], 0 });
                    Efield_p[specIndex] = layerBehindStack.PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * Vector.Create(new Complex<double>[] { E_trans_p[specIndex], 0 });
                }
            }

            // within the layerstack
            else
            {
                int layerIndex = GetLayerIndexAtPosition(position);

                for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                {
                    Efield_s[specIndex] = layerStack[layerIndex].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * layerStack[layerIndex].E_end_s[specIndex];
                    Efield_p[specIndex] = layerStack[layerIndex].PropagationMatrix(spectrum.data[specIndex].lambda, specIndex, position)
                        * layerStack[layerIndex].E_end_p[specIndex];
                }
            }

            return (Efield_s, Efield_p);
        }

        /// <summary>
        /// returns the absolute amount of incoming photons in 1/(s*m^2) and the reflected, absorbed (in each layer) and transmitted ratio of photons (with respect to the incoming amount) within a certain wavelength range
        /// </summary>
        /// <param name="fromLambda">wavelength in meter, from which the photons are counted (photons with shorter wavelength / higher energy are NOT counted) (NaN -> count all photons down to lambda = 0)</param>
        /// <param name="toLambda">wavelength in meter, up to which the photons are counted (photons with longer wavelength / lower energy are NOT counted) (NaN -> count all photons up to lambda = inf)</param>
        public (double incomingAbsolute, double reflectedFactor, double[] absorbedFactor, double transmittedFactor) GetPhotonsInReflectionAbsorptionTransmission(
            double fromLambda = double.NaN, double toLambda = double.NaN)
        {
            int lambdaStartIndex = GetIndexInSpectrumFromLambda(fromLambda, 0);
            int lambdaStopIndex = GetIndexInSpectrumFromLambda(toLambda, spectrum.data.Length - 1);

            double incoming = 0;
            double reflected = 0;
            double[] absorbed = new double[layerStack.Length];
            double transmitted = 0;

            // #photons [1 / (s * m²)] = intensity [W/m² = J / (s * m²)] / photon energy [J]
            // #photons = I / E = I / (h * f) = I * λ / (h * c)
            // division is done at the end

            // incoming, reflected transmitted
            var poynting_start = GetLambdaDependentPoyntingVectorAtPosition(0);
            var poynting_end = GetLambdaDependentPoyntingVectorAtPosition(lengthOfStack);
            for (int specIndex = lambdaStartIndex; specIndex <= lambdaStopIndex; specIndex++)
            {
                incoming += spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda * spectrum.data[specIndex].lambda;
                reflected += (spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda - poynting_start.s[specIndex] - poynting_start.p[specIndex])
                    * spectrum.data[specIndex].lambda;
                transmitted += (poynting_end.s[specIndex] + poynting_end.p[specIndex]) * spectrum.data[specIndex].lambda;
            }

            // absorbed
            for (int layerIndex = 0; layerIndex < layerStack.Length; layerIndex++)
            {
                var poynting_layerStart = GetLambdaDependentPoyntingVectorAtPosition(layerStack[layerIndex].layerPositionStart);
                var poynting_layerEnd = GetLambdaDependentPoyntingVectorAtPosition(layerStack[layerIndex].layerPositionEnd);

                for (int specIndex = lambdaStartIndex; specIndex <= lambdaStopIndex; specIndex++)
                {
                    absorbed[layerIndex] += (poynting_layerStart.s[specIndex] - poynting_layerEnd.s[specIndex]) * spectrum.data[specIndex].lambda;
                    absorbed[layerIndex] += (poynting_layerStart.p[specIndex] - poynting_layerEnd.p[specIndex]) * spectrum.data[specIndex].lambda;
                }
            }

            // if there are no photons, return zeros
            if (incoming == 0)
                return (0, 0, new double[layerStack.Length], 0);

            reflected /= incoming;
            for (int i = 0; i < absorbed.Length; i++)
                absorbed[i] /= incoming;
            transmitted /= incoming;
            incoming /= physConstants.h * physConstants.c;

            return (incoming, reflected, absorbed, transmitted);
        }

        /// <summary>
        /// returns the local absorption of photons in 1/(s*m³) at a given position within a certain lambda range
        /// </summary>
        /// <param name="position">depth in meter, where the absoption is calculated</param>
        /// <param name="fromLambda">wavelength in meter, from which the photones are added up (photons with shorter wavelength / higher energy are NOT counted) (NaN -> count all photons down to lambda = 0)</param>
        /// <param name="toLambda">wavelength in meter, up to which the photons are counted (photons with longer wavelength / lower energy are NOT counted) (NaN -> count all photons up to lambda = inf)</param>
        public double GetLocalAbsorption(double position, double fromLambda = double.NaN, double toLambda = double.NaN)
        {
            // before layerstack
            if (position < 0)
                return 0;

            int lambdaIndexStart = GetIndexInSpectrumFromLambda(fromLambda, 0);
            int lambdaIndexEnd = GetIndexInSpectrumFromLambda(toLambda, spectrum.data.Length - 1);

            double absorption = 0;

            var Efield = GetEfieldAtPosition(position);

            // after layerstack
            if (position > lengthOfStack)
            {
                for (int specIndex = lambdaIndexStart; specIndex <= lambdaIndexEnd; specIndex++)
                {
                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = Complex<double>.Cos(layerBehindStack.angleOfIncidence[specIndex]);
                    Complex<double> k_z = 2 * Math.PI * n * cosTheta / spectrum.data[specIndex].lambda;

                    // s polarized
                    absorption += (n * cosTheta * (k_z * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.s[specIndex][0] - Efield.s[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                    // version of Steven J Byrnes (https://arxiv.org/abs/1603.02720)
                    //absorption += (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared * (n * cosTheta * k_z).Im / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum[specIndex].lambda;

                    // p polarized
                    absorption += (n * cosTheta.Conjugate() * (k_z * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.p[specIndex][0] - Efield.p[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                }

                absorption /= physConstants.h * physConstants.c;
                return absorption;
            }

            // within the layerstack
            int layerIndex = GetLayerIndexAtPosition(position);

            for (int specIndex = lambdaIndexStart; specIndex <= lambdaIndexEnd; specIndex++)
            {
                Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                Complex<double> n = layerStack[layerIndex].material.propertiesOptics.n_toSpectrum[specIndex];
                Complex<double> cosTheta = Complex<double>.Cos(layerStack[layerIndex].angleOfIncidence[specIndex]);
                Complex<double> k_z = 2 * Math.PI * n * cosTheta / spectrum.data[specIndex].lambda;

                // s polarized
                absorption += (n * cosTheta * (k_z * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared
                        - k_z.Conjugate() * (Efield.s[specIndex][0] - Efield.s[specIndex][1]).MagnitudeSquared)).Im
                        / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                // version of Steven J Byrnes (https://arxiv.org/abs/1603.02720)
                //absorption += (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared * (n * cosTheta * k_z).Im / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum[specIndex].lambda;

                // p polarized
                absorption += (n * cosTheta.Conjugate() * (k_z * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).MagnitudeSquared
                        - k_z.Conjugate() * (Efield.p[specIndex][0] - Efield.p[specIndex][1]).MagnitudeSquared)).Im
                        / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
            }

            absorption /= physConstants.h * physConstants.c;
            return absorption;
        }

        /// <summary>
        /// returns the local generation (absorption minus FCA e.g. in TCOs) of photons in 1/(s*m³) at a given position within a certain lambda range 
        /// </summary>
        /// <param name="position">depth in meter, where the absoption is calculated</param>
        /// <param name="fromLambda">wavelength in meter, from which the photones are added up (photons with shorter wavelength / higher energy are NOT counted) (NaN -> count all photons down to lambda = 0)</param>
        /// <param name="toLambda">wavelength in meter, up to which the photons are counted (photons with longer wavelength / lower energy are NOT counted) (NaN -> count all photons up to lambda = inf)</param>
        public double GetLocalGeneration(double position, double fromLambda = double.NaN, double toLambda = double.NaN)
        {
            // before layerstack
            if (position < 0)
                return 0;

            int lambdaIndexStart = GetIndexInSpectrumFromLambda(fromLambda, 0);
            int lambdaIndexEnd = GetIndexInSpectrumFromLambda(toLambda, spectrum.data.Length - 1);

            double absorption = 0;

            var Efield = GetEfieldAtPosition(position);

            // after layerstack
            if (position > lengthOfStack)
            {
                for (int specIndex = lambdaIndexStart; specIndex <= lambdaIndexEnd; specIndex++)
                {
                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerBehindStack.material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = Complex<double>.Cos(layerBehindStack.angleOfIncidence[specIndex]);
                    Complex<double> k_z = 2 * Math.PI * n * cosTheta / spectrum.data[specIndex].lambda;

                    // s polarized
                    absorption += (n * cosTheta * (k_z * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.s[specIndex][0] - Efield.s[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                    // version of Steven J Byrnes (https://arxiv.org/abs/1603.02720)
                    //absorption += (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared * (n * cosTheta * k_z).Im / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum[specIndex].lambda;

                    // p polarized
                    absorption += (n * cosTheta.Conjugate() * (k_z * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.p[specIndex][0] - Efield.p[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                }

                absorption /= physConstants.h * physConstants.c;
                return absorption;
            }

            // within the layerstack
            int layerIndex = GetLayerIndexAtPosition(position);
            double bandgapAtPosition;
            if (layerStack[layerIndex].material.propertiesSemiconductor != null)
                bandgapAtPosition = layerStack[layerIndex].material.propertiesSemiconductor.Egap;
            else
                bandgapAtPosition = 0.5;

            var bandgapInLambda = Misc.ConverteVEnergyInWavelength(bandgapAtPosition) * 1.1;

            for (int specIndex = lambdaIndexStart; specIndex <= lambdaIndexEnd; specIndex++)
            {
                if (spectrum.data[specIndex].lambda < bandgapInLambda)
                {

                    Complex<double> cosTheta0 = Complex<double>.Cos(layerBeforeStack.angleOfIncidence[specIndex]);
                    Complex<double> n = layerStack[layerIndex].material.propertiesOptics.n_toSpectrum[specIndex];
                    Complex<double> cosTheta = Complex<double>.Cos(layerStack[layerIndex].angleOfIncidence[specIndex]);
                    Complex<double> k_z = 2 * Math.PI * n * cosTheta / spectrum.data[specIndex].lambda;

                    // s polarized
                    absorption += (n * cosTheta * (k_z * (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.s[specIndex][0] - Efield.s[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                    // version of Steven J Byrnes (https://arxiv.org/abs/1603.02720)
                    //absorption += (Efield.s[specIndex][0] + Efield.s[specIndex][1]).MagnitudeSquared * (n * cosTheta * k_z).Im / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum[specIndex].lambda;

                    // p polarized
                    absorption += (n * cosTheta.Conjugate() * (k_z * (Efield.p[specIndex][0] + Efield.p[specIndex][1]).MagnitudeSquared
                            - k_z.Conjugate() * (Efield.p[specIndex][0] - Efield.p[specIndex][1]).MagnitudeSquared)).Im
                            / cosTheta0.Re / (2 * physConstants.mu0 * physConstants.c) * spectrum.data[specIndex].lambda;
                }
            }

            absorption /= physConstants.h * physConstants.c;
            return absorption;
        }


        /// <summary>
        /// returns the wavelength-dependent optical distribution into all absorbed, reflected and transmitted photons
        /// </summary>
        public (double wavelength, double reflected, double[] absorbed, double transmitted)[] GetOpticalEQE()
        {
            (double wavelength, double reflected, double[] absorbed, double transmitted)[] opticalEQE = new (double wavelength, double reflected, double[] absorbed, double transmitted)[spectrum.data.Length];

            var poynting_start = GetLambdaDependentPoyntingVectorAtPosition(0);
            var poynting_end = GetLambdaDependentPoyntingVectorAtPosition(lengthOfStack);

            ((double[] s, double[] p) start, (double[] s, double[] p) stop)[] poynting_layer = new ((double[] s, double[] p) start, (double[] s, double[] p) stop)[layerStack.Length];

            for (int layerIndex = 0; layerIndex < layerStack.Length; layerIndex++)
                poynting_layer[layerIndex] =
                    (GetLambdaDependentPoyntingVectorAtPosition(layerStack[layerIndex].layerPositionStart), GetLambdaDependentPoyntingVectorAtPosition(layerStack[layerIndex].layerPositionEnd));

            for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
            {
                double incomingAmount = spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda * spectrum.data[specIndex].lambda;


                double reflectedAmount = (spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda - poynting_start.s[specIndex] - poynting_start.p[specIndex]) * spectrum.data[specIndex].lambda;
                double transmittedAmount = (poynting_end.s[specIndex] + poynting_end.p[specIndex]) * spectrum.data[specIndex].lambda;

                double[] absorbed = new double[layerStack.Length];
                for (int layerIndex = 0; layerIndex < layerStack.Length; layerIndex++)
                {


                    absorbed[layerIndex] += (poynting_layer[layerIndex].start.s[specIndex] - poynting_layer[layerIndex].stop.s[specIndex]) * spectrum.data[specIndex].lambda;
                    absorbed[layerIndex] += (poynting_layer[layerIndex].start.p[specIndex] - poynting_layer[layerIndex].stop.p[specIndex]) * spectrum.data[specIndex].lambda;
                    absorbed[layerIndex] /= incomingAmount;
                }


                opticalEQE[specIndex] = (spectrum.data[specIndex].lambda, reflectedAmount / incomingAmount, absorbed, transmittedAmount / incomingAmount);
            }

            return opticalEQE;
        }

        // internal methods █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the index of the nearest lambda in the spectrum
        /// </summary>
        /// <param name="lambda">wavelength, which is searched for</param>
        /// <param name="defaultIndex">index, which is returned, if lambda is NaN</param>
        int GetIndexInSpectrumFromLambda(double lambda, int defaultIndex)
        {
            int index = defaultIndex;
            if (!double.IsNaN(lambda))
            {
                double deltaLambda = double.PositiveInfinity;
                for (int i = 0; i < spectrum.data.Length; i++)
                {
                    if (Math.Abs(spectrum.data[i].lambda - lambda) < deltaLambda)
                    {
                        deltaLambda = Math.Abs(spectrum.data[i].lambda - lambda);
                        index = i;
                    }
                }
            }
            return index;
        }
        /// <summary>
        /// get layerindex, where the given position is in
        /// </summary>
        /// <param name="position"></param>
        int GetLayerIndexAtPosition(double position)
        {
            for (int i = 0; i < layerStack.Length; i++)
            {
                if (position >= layerStack[i].layerPositionStart && position <= layerStack[i].layerPositionEnd)
                    return i;
            }
            return -1;
        }
    }
}