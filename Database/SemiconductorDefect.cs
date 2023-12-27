using BasicLib;
using Database;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Database
{
    public class SemiconductorDefect
    {
        /// <summary>
        /// gives the Charge type of the defect
        /// </summary>
        TypeOfDefect typeOfDefect { get; set; }


        public string typeOfDefectString
        {
            get
            {
                switch (typeOfDefect)
                {
                    case TypeOfDefect.Donor:
                        return "Donor";
                    case TypeOfDefect.Neutral:
                        return "Neutral";
                    case TypeOfDefect.Acceptor:
                        return "Acceptor";
                    default:
                        return "";
                }
            }
            set
            {

            }
        }

        /// <summary>
        /// specifies the energy band from which the energetic position  is measured
        /// </summary>
        DefectReferenceBand defectReferenceBand { get; set; }

        public string referenceBandString
        {
            get
            {
                switch (defectReferenceBand)
                {
                    case DefectReferenceBand.Ei:
                        return "Intrinsic Level";
                    case DefectReferenceBand.Ev:
                        return "Valence Band Maximum";
                    case DefectReferenceBand.Ec:
                        return "Conduction Band Minimum";
                    default:
                        return "";
                }
            }
            set
            {

            }
        }

        /// <summary>
        /// energetic position of the defct in SRH model, in J
        /// </summary>
        public double energeticPosition { get; private set; }
        public double materialEnergeticPosition { get; set; }

        /// <summary>
        /// Density of defects in SRH model, in m^-3
        /// </summary>
        public double defectDensity { get; set; }

        /// <summary>
        /// electron capture cross section in SRH model, in m^2
        /// </summary>
        public double electronCaptureCrosssection { get; set; }

        /// <summary>
        ///  hole capture cross section in SRH model, in m^2
        /// </summary>
        public double holeCaptureCrosssection { get; set; }

        public string energeticPositionGUI
        {
            get
            {
                return getEnergeticPosition(PropertiesSemiconductor.Egap, PropertiesSemiconductor.getMaterialIntrinsicLevel()).ToString();
            }
            set
            {

            }
        }
        public string defectDensityGUI
        {
            get
            {
                return (defectDensity * 1e-6).ToString("G4");
            }
            set
            {

            }
        }

        public string electronCaptureCrosssectionGUI
        {
            get
            {
                return (electronCaptureCrosssection * 1e4).ToString("G4");
            }
            set
            {

            }
        }

        public string holeCaptureCrosssectionGUI
        {
            get
            {
                return (holeCaptureCrosssection * 1e4).ToString("G4");
            }
            set
            {

            }
        }

        PropertiesSemiconductor PropertiesSemiconductor { get; set; }

        /// <summary>
        /// constructs a new defect (for a semiconductor material) from basic defect properties (trap density, cross sections, energetic position)
        /// Attention with the units! Construction with energy in eV, calculations with energies in J (see comment to energeticPosition above)
        /// </summary>
        /// <param name="name">name of defect</param>
        /// <param name="ID">ID of defect for database</param>
        /// <param name="typeOfDefect">charge type of defect</param>
        /// <param name="energeticPosition">in eV (current description: distance from intrinsic level)</param>
        /// <param name="defectDensity">in m^-3</param>
        /// <param name="electronCaptureCrosssection">in m^2</param>
        /// <param name="holeCaptureCrosssection">in m^2</param>
        public SemiconductorDefect(TypeOfDefect typeOfDefect, DefectReferenceBand defectReferenceBand, double energeticPosition, double defectDensity, double electronCaptureCrosssection, double holeCaptureCrosssection, PropertiesSemiconductor PropertiesSemiconductor)
        {
            this.typeOfDefect = typeOfDefect;
            this.defectReferenceBand = defectReferenceBand;
            this.energeticPosition = energeticPosition;//  * physConstants.e;
            this.defectDensity = defectDensity;
            this.electronCaptureCrosssection = electronCaptureCrosssection;
            this.holeCaptureCrosssection = holeCaptureCrosssection;
            this.PropertiesSemiconductor = PropertiesSemiconductor;


            switch (defectReferenceBand)
            {
                case DefectReferenceBand.Ei:
                    materialEnergeticPosition = energeticPosition + PropertiesSemiconductor.getMaterialIntrinsicLevel();
                    break;
                case DefectReferenceBand.Ev:
                    materialEnergeticPosition = energeticPosition + (PropertiesSemiconductor.chemicalPotential - PropertiesSemiconductor.Egap);
                    break;
                case DefectReferenceBand.Ec:
                    materialEnergeticPosition = energeticPosition + PropertiesSemiconductor.chemicalPotential ;
                    break;
                default:
                    materialEnergeticPosition = energeticPosition + PropertiesSemiconductor.getMaterialIntrinsicLevel();
                    break;
            }

        }

        /// <summary>
        /// calculates the electron lifetime for a specific defect in a semiconductor material
        /// </summary>
        /// <param name="electronThermalVelocity"> thermal velocity of electrons of the material</param>
        /// <returns></returns>
        public double tau_n(double electronThermalVelocity)
        {
            return 1 / (defectDensity * electronCaptureCrosssection * electronThermalVelocity);
        }

        /// <summary>
        /// calculates the hole lifetime for a specific defect in a semiconductor
        /// </summary>
        /// <param name="holeThermalVelocity">thermal velocity of holes of the material</param>
        /// <returns></returns>
        public double tau_p(double holeThermalVelocity)
        {
            return 1 / (defectDensity * holeCaptureCrosssection * holeThermalVelocity);
        }

        /// <summary>
        /// calculates the reference densities for electrons
        /// </summary>
        /// <param name="Nintr">intrinsic carrier density</param>
        /// <param name="temperature"> in K</param>
        /// <returns></returns>
        public double nTrapRef(double Nintr, double temperature, double localBandGap, double intrinsicLevel)
        {
            return Nintr * Math.Exp(getEnergeticPosition(localBandGap, intrinsicLevel) / (physConstants.kB * temperature));
        }

        /// <summary>
        /// calculates the reference densities for holes
        /// </summary>
        /// <param name="Nintr">intrinsic carrier density</param>
        /// <param name="temperature"> in K</param>
        /// <returns></returns>
        public double pTrapRef(double Nintr, double temperature, double localBandGap, double intrinsicLevel)
        {
            return Nintr * Math.Exp(-getEnergeticPosition(localBandGap, intrinsicLevel) / (physConstants.kB * temperature));
        }

        public double getEnergeticPosition(double localBandGap, double intrinsicLevel)
        {

            switch (defectReferenceBand)
            {
                case DefectReferenceBand.Ei:
                    return energeticPosition ;
                    break;
                case DefectReferenceBand.Ev:
                    return energeticPosition -(intrinsicLevel -( PropertiesSemiconductor.chemicalPotential - localBandGap));
                    break;
                case DefectReferenceBand.Ec:
                    return energeticPosition + (PropertiesSemiconductor.chemicalPotential - intrinsicLevel);
                    break;
                default:
                    return energeticPosition ;
                    break;
            }
        }


    }
}