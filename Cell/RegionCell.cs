using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BasicLib;
using Geometry;
using Database;
using TransferMatrix;
using Extreme.Mathematics.Curves;

namespace Cell
{
    public class RegionCell : Region
    {
        /// <summary>
        /// material which is used as a front contact (TCO)
        /// </summary>
        public Material frontContact { get; private set; }
        /// <summary>
        /// thickness of the front contact (TCO)
        /// </summary>
        public double thicknessFrontContact { get; set; }
        /// <summary>
        /// material which is used as a front grid (null means no grid)
        /// </summary>
        public Material frontGrid { get; private set; }
        /// <summary>
        /// thickness of the front grid
        /// </summary>
        public double thicknessFrontGrid { get; set; }

        /// <summary>
        /// material which is used as a back contact
        /// </summary>
        public Material backContact { get; private set; }
        /// <summary>
        /// thickness of the back contact
        /// </summary>
        public double thicknessBackContact { get; set; }
        /// <summary>
        /// material which is used as a back grid (null means no grid)
        /// </summary>
        public Material backGrid { get; private set; }
        /// <summary>
        /// thickness of the back grid
        /// </summary>
        public double thicknessBackGrid { get; set; }

        /// <summary>
        /// pn junction which is present within this region
        /// </summary>
        public pnJunction pnJunction { get; set; }
        /// <summary>
        /// determines, whether this region counts to the active area for determining the efficiency
        /// </summary>
        public bool countsAsActiveArea { get; private set; }

        #region optics
        /// <summary>
        /// optical model in this region
        /// </summary>
        public (double roughnessFrontGrid, double roughnessFrontContact, double roughnessAbsorber, double roughnessBackContact, double roughnessBackGrid, Material materialBefore, (Material material, double roughness) behind, List<(Material material, double thickness)> incoherent,
            List<(Material material, double thickness, double roughness)> aboveFrontGrid, List<(Material material, double thickness, double roughness)> aboveAbsorber,
            List<(Material material, double thickness, double roughness)> belowAbsorber, List<(Material material, double thickness, double roughness)> belowBackGrid) opticalModel { get; private set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses in atmospheric scattering additionally to the used spectrum
        /// </summary>
        public double opticFactor_additionalAtmosphereScattering { get; set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses in shading due to clouds or laboratory arrangements
        /// </summary>
        public double opticFactor_shading { get; private set; }
        /// <summary>
        /// optical factor x*Iph: amount of effective illuminated area (cos of angle between incident sun angle and perpendicular ray on solar cell)
        /// </summary>
        public double opticFactor_effectiveArea { get; set; }
        /// <summary>
        /// optical factor x*Iph: transparency of the grid
        /// </summary>
        public double opticFactor_transmissionGrid { get; set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses due to relfection
        /// </summary>
        public double opticFactor_reflection { get; set; }
        /// <summary>
        /// optical factors (1-Σx)*Iph: losses due to parasitic absorption
        /// </summary>
        public (double factor, bool isAbsorber, int materialID, string materialName)[] opticFactor_absorption { get; set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses due to transmission
        /// </summary>
        public double opticFactor_transmission { get; set; }
        /// <summary>
        /// model for transfer matrix method
        /// </summary>
        public ModelTMM modelOpticsTMM { get; private set; }
        #endregion

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        public RegionCell()
        {
        }

        // Set properties for cell region ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// sets all properties, a cell region needs to have
        /// </summary>
        public override void SetProperties(double[] preferencesArray, pointType regionType)
        {
            frontContact = Data.GetMaterialFromID((int)Math.Round(preferencesArray[0]));
            thicknessFrontContact = preferencesArray[1];

            frontGrid = Data.GetMaterialFromID((int)Math.Round(preferencesArray[2]));
            thicknessFrontGrid = preferencesArray[3];

            backContact = Data.GetMaterialFromID((int)Math.Round(preferencesArray[4]));
            thicknessBackContact = preferencesArray[5];

            backGrid = Data.GetMaterialFromID((int)Math.Round(preferencesArray[6]));
            thicknessBackGrid = preferencesArray[7];

            pnJunction = Data.GetPNjunctionFromID((int)Math.Round(preferencesArray[8]));

            opticFactor_shading = preferencesArray[9];

            countsAsActiveArea = (int)Math.Round(preferencesArray[10]) == 1 ? true : false;

            type = regionType;

            // module materials for P1 (pn material on back) and P3 (air on front)
            if (regionType == pointType.P3)
                frontContact = Data.GetMaterialFromID(990000000);
            if (regionType == pointType.P1)
                backContact = pnJunction.absorberMaterial;
        }

        public override void SetOpticalModel((double roughnessFrontGrid, double roughnessFrontContact, double roughnessAbsorber, double roughnessBackContact, double roughnessBackGrid,
            int IDbefore, (int ID, double roughness) behind, List<(int ID, double thickness)> incoherent,
            List<(int ID, double thickness, double roughness)> aboveFrontGrid, List<(int ID, double thickness, double roughness)> aboveAbsorber,
            List<(int ID, double thickness, double roughness)> belowAbsorber, List<(int ID, double thickness, double roughness)> belowBackGrid) opticalModelInString)
        {
            Material materialBefore = Data.GetMaterialFromID(opticalModelInString.IDbefore);

            (Material material, double roughness) behind = (Data.GetMaterialFromID(opticalModelInString.behind.ID), opticalModelInString.behind.roughness);

            List<(Material material, double thickness)> incoherent = new List<(Material material, double thickness)>();
            foreach (var incoh in opticalModelInString.incoherent)
                incoherent.Add((Data.GetMaterialFromID(incoh.ID), incoh.thickness));

            List<(Material material, double thickness, double roughness)> aboveFrontGrid = new List<(Material material, double thickness, double roughness)>();
            foreach (var mat in opticalModelInString.aboveFrontGrid)
                aboveFrontGrid.Add((Data.GetMaterialFromID(mat.ID), mat.thickness, mat.roughness));

            List<(Material material, double thickness, double roughness)> aboveAbsorber = new List<(Material material, double thickness, double roughness)>();
            foreach (var mat in opticalModelInString.aboveAbsorber)
                aboveAbsorber.Add((Data.GetMaterialFromID(mat.ID), mat.thickness, mat.roughness));

            List<(Material material, double thickness, double roughness)> belowAbsorber = new List<(Material material, double thickness, double roughness)>();
            foreach (var mat in opticalModelInString.belowAbsorber)
                belowAbsorber.Add((Data.GetMaterialFromID(mat.ID), mat.thickness, mat.roughness));

            List<(Material material, double thickness, double roughness)> belowBackGrid = new List<(Material material, double thickness, double roughness)>();
            foreach (var mat in opticalModelInString.belowBackGrid)
                belowBackGrid.Add((Data.GetMaterialFromID(mat.ID), mat.thickness, mat.roughness));

            opticalModel = (opticalModelInString.roughnessFrontGrid, opticalModelInString.roughnessFrontContact, opticalModel.roughnessAbsorber, opticalModelInString.roughnessBackContact, opticalModelInString.roughnessBackGrid,
                materialBefore, behind, incoherent, aboveFrontGrid, aboveAbsorber, belowAbsorber, belowBackGrid);
        }

        public void CalculateOpticalCoefficients(OpticMode opticMode, Spectrum spectrum, double angleOfIncidence = 0, double sunElevation = 90)
        {
            switch (opticMode)
            {
                case OpticMode.coefficient:
                    opticFactor_additionalAtmosphereScattering = Misc.FactorAtmosphericScattering(sunElevation);
                    opticFactor_effectiveArea = Misc.FactorEffectiveArea(angleOfIncidence);
                    opticFactor_transmissionGrid = frontGrid.propertiesOptics.simpleLightTransmissionCoefficient;
                    opticFactor_reflection = 1 - frontContact.propertiesOptics.simpleLightTransmissionCoefficient;
                    opticFactor_absorption = new (double factor, bool isAbsorber, int materialID, string materialName)[0];
                    opticFactor_transmission = 0;
                    break;

                case OpticMode.lambertBeer:
                    //  ██╗ lambert beer
                    //  ╚═╝
                    var absorptionList = new List<(double factor, bool isAbsorber, int materialID, string materialName)>();
                    double currentLight = 1;
                    
                    // incoherent layers
                    for (int materialIndex = 0; materialIndex < opticalModel.incoherent.Count; materialIndex++)
                    {
                        var abs = currentLight - currentLight * Math.Exp(-opticalModel.incoherent[materialIndex].material.propertiesOptics.lambertBeerAbsorptionCoefficient * opticalModel.incoherent[materialIndex].thickness);
                        absorptionList.Add((abs, false, opticalModel.incoherent[materialIndex].material.ID, opticalModel.incoherent[materialIndex].material.name));
                        currentLight -= abs;
                    }
                    
                    // layers above front grid
                    foreach (var layer in opticalModel.aboveFrontGrid)
                    {
                        var abs = currentLight - currentLight * Math.Exp(-layer.material.propertiesOptics.lambertBeerAbsorptionCoefficient * layer.thickness);
                        absorptionList.Add((abs, false, layer.material.ID, layer.material.name));
                        currentLight -= abs;
                    }

                    // front grid
                    if (frontGrid.ID != 990000000)
                    {
                        var abs = currentLight - currentLight * Math.Exp(-frontGrid.propertiesOptics.lambertBeerAbsorptionCoefficient * thicknessFrontGrid);
                        absorptionList.Add((abs, false, frontGrid.ID, frontGrid.name));
                        currentLight -= abs;
                    }

                    // front contact
                    var absop = currentLight - currentLight * Math.Exp(-frontContact.propertiesOptics.lambertBeerAbsorptionCoefficient * thicknessFrontContact);
                    absorptionList.Add((absop, false, frontContact.ID, frontContact.name));
                    currentLight -= absop;

                    // layers above absorber
                    foreach (var layer in opticalModel.aboveAbsorber)
                    {
                        var abs = currentLight - currentLight * Math.Exp(-layer.material.propertiesOptics.lambertBeerAbsorptionCoefficient * layer.thickness);
                        absorptionList.Add((abs, false, layer.material.ID, layer.material.name));
                        currentLight -= abs;
                    }

                    //  ██╗ set to region
                    //  ╚═╝
                    opticFactor_additionalAtmosphereScattering = Misc.FactorAtmosphericScattering(sunElevation);
                    opticFactor_effectiveArea = Misc.FactorEffectiveArea(angleOfIncidence);
                    opticFactor_transmissionGrid = Math.Exp(-frontGrid.propertiesOptics.lambertBeerAbsorptionCoefficient * thicknessFrontGrid);
                    opticFactor_reflection = 0;
                    opticFactor_absorption = absorptionList.ToArray();
                    opticFactor_transmission = currentLight * Math.Exp(-pnJunction.absorberMaterial.propertiesOptics.lambertBeerAbsorptionCoefficient * pnJunction.thicknessAbsorberLayer);
                    break;

                case OpticMode.transferMatrixInCell:
                    //  ██╗ transfer-matrix
                    //  ╚═╝
                    List<(Material material, double thickness, double roughnessOnTop, bool isAbsorber)> materialStack = new List<(Material material, double thickness, double roughnessOnTop, bool isAbsorber)>();

                    // layers above front grid
                    foreach (var layer in opticalModel.aboveFrontGrid)
                        materialStack.Add((layer.material, layer.thickness, layer.roughness, false));

                    // front grid
                    if (frontGrid.ID != 990000000)
                        materialStack.Add((frontGrid, thicknessFrontGrid, opticalModel.roughnessFrontGrid, false));

                    // front contact
                    materialStack.Add((frontContact, thicknessFrontContact, opticalModel.roughnessFrontContact, false));

                    // layers above absorber
                    foreach (var layer in opticalModel.aboveAbsorber)
                        materialStack.Add((layer.material, layer.thickness, layer.roughness, false));

                    // absorber
                    materialStack.Add((pnJunction.absorberMaterial, pnJunction.thicknessAbsorberLayer, opticalModel.roughnessAbsorber, true));

                    // layers below absorber
                    foreach (var layer in opticalModel.belowAbsorber)
                        materialStack.Add((layer.material, layer.thickness, layer.roughness, false));

                    // back contact
                    materialStack.Add((backContact, thicknessBackContact, opticalModel.roughnessBackContact, false));

                    // back grid
                    if (backGrid.ID != 990000000)
                        materialStack.Add((backGrid, thicknessBackGrid, opticalModel.roughnessBackGrid, false));

                    // layers below back grid
                    foreach (var layer in opticalModel.belowBackGrid)
                        materialStack.Add((layer.material, layer.thickness, layer.roughness, false));

                    modelOpticsTMM = new ModelTMM(opticalModel.materialBefore, opticalModel.behind, materialStack.Select(m => (m.material, m.thickness, m.roughnessOnTop)).ToArray(), spectrum, 1, angleOfIncidence);
                    var photonsRAT = modelOpticsTMM.GetPhotonsInReflectionAbsorptionTransmission(0e-9, physConstants.h * physConstants.c / (physConstants.e * pnJunction.bandgapINeV));

                    //  ██╗ incoherent absorption
                    //  ╚═╝
                    // get absorption of all absorbers to weight absorption function of incoherent material with absorption profile of all absorbers
                    List<(double wavelength, double deltaWavelength, List<double> incoherentAbsorption, double absorperAbsorption)> incoherent = new List<(double wavelength, double deltaWavelength, List<double> incoherentAbsorption, double absorperAbsorption)>();
                    var opticalEQE = modelOpticsTMM.GetOpticalEQE();
                    for (int wavelengthIndex = 0; wavelengthIndex < opticalEQE.Length; wavelengthIndex++)
                    {
                        double abs = 0;
                        for (int layerIndex = 0; layerIndex < materialStack.Count; layerIndex++)
                            if (materialStack[layerIndex].isAbsorber)
                                abs += opticalEQE[wavelengthIndex].absorbed[layerIndex];

                        double deltaWavelength = 0;
                        if (wavelengthIndex == 0)
                            deltaWavelength = opticalEQE[1].wavelength - opticalEQE[0].wavelength;
                        else if (wavelengthIndex == opticalEQE.Length - 1)
                            deltaWavelength = opticalEQE[opticalEQE.Length - 1].wavelength - opticalEQE[opticalEQE.Length - 2].wavelength;
                        else
                            deltaWavelength = (opticalEQE[wavelengthIndex + 1].wavelength - opticalEQE[wavelengthIndex - 1].wavelength) / 2;

                        incoherent.Add((opticalEQE[wavelengthIndex].wavelength, deltaWavelength, new List<double>(), abs));
                    }

                    // get absorption of incoherent materials
                    for (int materialIndex = 0; materialIndex < opticalModel.incoherent.Count; materialIndex++)
                    {
                        List<double> wavelengths = new List<double>();
                        List<double> absorption = new List<double>();
                        if (opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData.Length < 2)
                        {
                            // duplicate point if only single wavelength data is given (otherwise no spline exists)
                            double alpha = 4.0 * Math.PI * opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData.First().n.Im
                                / opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData.First().lambda;
                            double abs = 1.0 - Math.Exp(-alpha * opticalModel.incoherent[materialIndex].thickness);
                            absorption = new List<double>() { abs, abs };
                            wavelengths = new List<double>() { 400, 1000 };
                        }
                        else
                            for (int i = 0; i < opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData.Length; i++)
                            {
                                double alpha = 4.0 * Math.PI * opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData[i].n.Im
                                    / opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData[i].lambda; // alpha = 4 * pi * k / lambda
                                absorption.Add(1.0 - Math.Exp(-alpha * opticalModel.incoherent[materialIndex].thickness));
                                wavelengths.Add(opticalModel.incoherent[materialIndex].material.propertiesOptics.n_rawData[i].lambda);
                            }
                        CubicSpline spline = new CubicSpline(wavelengths.ToArray(), absorption.ToArray());
                        for (int i = 0; i < incoherent.Count; i++)
                            incoherent[i] = (incoherent[i].wavelength, incoherent[i].deltaWavelength,
                                incoherent[i].incoherentAbsorption.Concat(new List<double>() { spline.ValueAt(incoherent[i].wavelength) }).ToList(), incoherent[i].absorperAbsorption);
                    }

                    List<(double absorption, int ID, string materialName)> incoherentAbsorption = new List<(double absorption, int ID, string materialName)>();
                    for (int materialIndex = 0; materialIndex < opticalModel.incoherent.Count; materialIndex++)
                        incoherentAbsorption.Add((incoherent.Sum(d => d.absorperAbsorption * d.incoherentAbsorption[materialIndex] * d.deltaWavelength) / incoherent.Sum(d => d.absorperAbsorption * d.deltaWavelength)
                            * (1 - photonsRAT.reflectedFactor), opticalModel.incoherent[materialIndex].material.ID, opticalModel.incoherent[materialIndex].material.name));

                    // add to absorption list
                    List<(double factor, bool isAbsorber, int materialID, string materialName)> absorbedList = new List<(double factor, bool isAbsorber, int materialID, string materialName)>();
                    double totalIncoherentAbsorption = 0;
                    for (int materialIndex = 0; materialIndex < opticalModel.incoherent.Count; materialIndex++)
                    {
                        totalIncoherentAbsorption += (1 - totalIncoherentAbsorption) * incoherentAbsorption[materialIndex].absorption;
                        absorbedList.Add(((1 - totalIncoherentAbsorption) * incoherentAbsorption[materialIndex].absorption, false, incoherentAbsorption[materialIndex].ID, incoherentAbsorption[materialIndex].materialName));
                    }
                    for (int i = 0; i < materialStack.Count; i++)
                        absorbedList.Add((photonsRAT.absorbedFactor[i] * (1 - totalIncoherentAbsorption), materialStack[i].isAbsorber, materialStack[i].material.ID, materialStack[i].material.name));

                    //  ██╗ set to region
                    //  ╚═╝
                    opticFactor_additionalAtmosphereScattering = Misc.FactorAtmosphericScattering(sunElevation);
                    opticFactor_effectiveArea = Misc.FactorEffectiveArea(angleOfIncidence);
                    opticFactor_transmissionGrid = Math.Exp(-frontGrid.propertiesOptics.lambertBeerAbsorptionCoefficient * thicknessFrontGrid);
                    opticFactor_reflection = photonsRAT.reflectedFactor;
                    opticFactor_absorption = absorbedList.ToArray();
                    opticFactor_transmission = photonsRAT.transmittedFactor * (1 - totalIncoherentAbsorption);

                    break;
            }
        }
    }
}