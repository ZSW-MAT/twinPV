using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using BasicLib;
using Geometry;
using Database;
using Extreme.Mathematics;

namespace Cell
{
    /// <summary>
    /// Class, with respresents a subclass of Point for the Cell simulation
    /// </summary>
    public class FiniteElementCell : FiniteElement
    {
        #region potentials
        /// <summary>
        /// electrical potential on the front side of the microcell
        /// </summary>
        public double phiFront { get; set; }
        /// <summary>
        /// electrical potential on the back side of the microcell
        /// </summary>
        public double phiBack { get; set; }
        /// <summary>
        /// initial guess of the electrical potiential on the front side
        /// </summary>
        public double phiFrontInit { get; set; }
        /// <summary>
        /// initial guess of the electrical potiential on the back side
        /// </summary>
        public double phiBackInit { get; set; }
        #endregion

        #region basic properties
        /// <summary>
        /// determines, whether this meshpoint counts to the active area for determining the efficiency
        /// </summary>
        public bool countsAsActiveArea { get; private set; }
        /// <summary>
        /// pn junction in this point
        /// </summary>
        public pnJunction pnJunction { get; private set; }
        /// <summary>
        /// variable for optimization (fraction of front grid thickness: 1 -> grid is as thick as in frontThicknessGrid, 0 -> grid is not present, values between 0 and 1 are allowed)
        /// </summary>
        public double frontGridDensity { get; set; } = double.NaN;
        /// <summary>
        /// power in watt which is converted to heat via Frontgrid, FrontTCO, Backgrid, BackTCO, Diode and Shuntresistance
        /// </summary>
        public double createdHeatPower { get; private set; }
        #endregion

        #region external contacts
        /// <summary>
        /// determines, if this point is an external cell front contact
        /// </summary>
        public bool isExternalCellFrontContact { get; set; } = false;
        /// <summary>
        /// contact resistance from external cell front contact wire to top contact of this cell [Ohm]
        /// </summary>
        public double contactResistanceExternalCellFrontContact { get; set; }
        /// <summary>
        /// determines, if this point is an external cell back contact
        /// </summary>
        public bool isExternalCellBackContact { get; set; } = false;
        /// <summary>
        /// contact resistance from external cell back contact wire to bottom contact of this cell [Ohm]
        /// </summary>
        public double contactResistanceExternalCellBackContact { get; set; }
        #endregion

        #region contact layers
        /// <summary>
        /// material which is used as a front contact (TCO)
        /// </summary>
        public Material frontContact { get; private set; }
        /// <summary>
        /// thickness of the front contact layer in meter
        /// </summary>
        public double thicknessFrontContact { get; private set; }

        /// <summary>
        /// material which is used as a front grid (null means no grid)
        /// </summary>
        public Material frontGrid { get; private set; }
        /// <summary>
        /// thickness of the front grid layer in meter
        /// </summary>
        public double thicknessFrontGrid { get; private set; }

        /// <summary>
        /// material which is used as a back contact
        /// </summary>
        public Material backContact { get; private set; }
        /// <summary>
        /// thickness of the back contact layer in meter
        /// </summary>
        public double thicknessBackContact { get; private set; }

        /// <summary>
        /// material which is used as a back grid (null means no grid)
        /// </summary>
        public Material backGrid { get; private set; }
        /// <summary>
        /// thickness of the back grid layer in meter
        /// </summary>
        public double thicknessBackGrid { get; private set; }
        #endregion

        #region resistors to neighbors
        /// <summary>
        /// list of resistors to neighbors (ordered in the same order as neighbor array) on the front side in Ohm
        /// </summary>
        public List<double> resistorsFront { get; private set; } = new List<double>();
        /// <summary>
        /// list of resistors to neighbors (ordered in the same order as neighbor array) on the front side in Ohm (only up to the voronoi edge to next microcell)
        /// </summary>
        public List<double> resistorsFrontToEdge { get; private set; } = new List<double>();
        /// <summary>
        /// list of contact layer resistors to neighbors (ordered in the same order as neighbor array) on the front side in Ohm (only up to the voronoi edge to next microcell)
        /// </summary>
        public List<double> resistorsFrontToEdgeContact { get; private set; } = new List<double>();
        /// <summary>
        /// list of grid resistors to neighbors (ordered in the same order as neighbor array) on the front side in Ohm (only up to the voronoi edge to next microcell) (NaN if there is no grid in this direction)
        /// </summary>
        public List<double> resistorsFrontToEdgeGrid { get; private set; } = new List<double>();
        /// <summary>
        /// list of resistors to neighbors (ordered in the same order as neighbor array) on the back side in Ohm
        /// </summary>
        public List<double> resistorsBack { get; private set; } = new List<double>();
        /// <summary>
        /// list of resistors to neighbors (ordered in the same order as neighbor array) on the back side in Ohm (only up to the voronoi edge to next microcell)
        /// </summary>
        public List<double> resistorsBackToEdge { get; private set; } = new List<double>();
        /// <summary>
        /// list of contact layer resistors to neighbors (ordered in the same order as neighbor array) on the back side in Ohm (only up to the voronoi edge to next microcell)
        /// </summary>
        public List<double> resistorsBackToEdgeContact { get; private set; } = new List<double>();
        /// <summary>
        /// list of grid resistors to neighbors (ordered in the same order as neighbor array) on the back side in Ohm (only up to the voronoi edge to next microcell) (NaN if there is no grid in this direction)
        /// </summary>
        public List<double> resistorsBackToEdgeGrid { get; private set; } = new List<double>();
        #endregion

        #region currents to neighbors
        /// <summary>
        /// list of currents on the front side to all neighbors in Ampere (ordered in the same order as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> currentsFront { get; private set; } = new List<double>();
        /// <summary>
        /// list of current densities on the front side to all neighbors in Ampere (ordered in the same order as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> currentdensitiesFront { get; private set; } = new List<double>();
        /// <summary>
        /// list of currents on the back side to all neighbors in Ampere (ordered in the same order as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> currentsBack { get; private set; } = new List<double>();
        /// <summary>
        /// list of current densities on the back side to all neighbors in Ampere (ordered in the same order as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> currentdensitiesBack { get; private set; } = new List<double>();
        #endregion

        #region current vectors
        /// <summary>
        /// all out flowing currents on the front side added up vectorially (always positve (=outgoing) or zero)
        /// </summary>
        public double[] IvecFront { get; private set; } = new double[2];
        /// <summary>
        /// all out flowing currents on the back side added up vectorially (always positve (=outgoing) or zero)
        /// </summary>
        public double[] IvecBack { get; private set; } = new double[2];
        #endregion

        #region module preferences
        /// <summary>
        /// determines, wheter this points phiFront is free (false) or bounded to another point (true)
        /// </summary>
        public bool hasModuleConnection_isInP3 { get; set; } = false;
        /// <summary>
        /// determines, wheter this points phiBack is free (false) or bounded to another point and hence used for interpolating its phiBack (true)
        /// </summary>
        public bool hasModuleConnection_isInCell { get; set; } = false;
        /// <summary>
        /// index of the points back potential, where this points front potential is linked to
        /// </summary>
        public int indexOfModuleConnectedPoint { get; set; } = -1;
        /// <summary>
        /// determines if this point is any point of a module interconnection or not (none)
        /// </summary>
        public pointType type { get; private set; }
        /// <summary>
        /// absolute contact resistacce in Ohm going from front- to back-contact within the P2 trench
        /// </summary>
        public double contactSeriesResistanceP2 { get; set; }
        #endregion

        #region optics
        /// <summary>
        /// optical factor (1-x)*Iph: losses in atmospheric scattering additionally to the used spectrum
        /// </summary>
        public double opticFactor_additionalAtmosphereScattering { get; private set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses in shading due to clouds or laboratory arrangements
        /// </summary>
        public double opticFactor_shading { get; private set; }
        /// <summary>
        /// optical factor x*Iph: amount of effective illuminated area (cos of angle between incident sun angle and perpendicular ray on solar cell)
        /// </summary>
        public double opticFactor_effectiveArea { get; private set; }
        /// <summary>
        /// optical factor x*Iph: transparency of the grid
        /// </summary>
        public double opticFactor_transmissionGrid { get; private set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses due to relfection
        /// </summary>
        public double opticFactor_reflection { get; private set; }
        /// <summary>
        /// optical factors (1-Σx)*Iph: losses due to parasitic absorption
        /// </summary>
        public (double factor, bool isAbsorber, int materialID, string materialName)[] opticFactor_absorption { get; private set; }
        /// <summary>
        /// optical factor (1-x)*Iph: losses due to transmission
        /// </summary>
        public double opticFactor_transmission { get; private set; }
        #endregion

        // set i-th differential equation variable ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// sets the i-th differential equation variable (i=0: phiFront, i=1: phiBack)
        /// </summary>
        /// <param name="differentialVariable">value which will be set to the variable of the point</param>
        /// <param name="index">index of the differential variable, which will be set</param>
        public override void SetDifferentialEquationVariable(double differentialVariable, int index)
        {
            switch (index)
            {
                case 0:
                    phiFront = differentialVariable;
                    phiFrontInit = differentialVariable;
                    break;
                case 1:
                    phiBack = differentialVariable;
                    phiBackInit = differentialVariable;
                    break;
                default:
                    break;
            }
        }

        // get all differential equation variables ██████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// gets all differential equation variables listed in a double array (double[] { phiFront, phiBack })
        /// </summary>
        /// <returns></returns>
        public override double[] GetDifferentialEquationVariables()
        {
            return new double[] { phiFront, phiBack };
        }

        // refresh electrical parameters ████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// refreshes the electrical parameters for each single point
        /// </summary>
        public void RefreshParametersFirst(ModelCell cell, double frontGridDensity = double.NaN)
        {
            //  ██╗ 
            //  ╚██╗ Set optics and electrics
            //  ██╔╝
            //  ╚═╝
            this.frontGridDensity = frontGridDensity;

            var enclosingRegions = cell.meshingAlgorithm.GetEnclosingRegions(position);
            frontContact = enclosingRegions.Last().frontContact;
            thicknessFrontContact = enclosingRegions.Last().thicknessFrontContact;
            frontGrid = enclosingRegions.Last().frontGrid.ID == 990000000 ? null : enclosingRegions.Last().frontGrid;
            thicknessFrontGrid = enclosingRegions.Last().thicknessFrontGrid;
            backContact = enclosingRegions.Last().backContact;
            thicknessBackContact = enclosingRegions.Last().thicknessBackContact;
            backGrid = enclosingRegions.Last().backGrid.ID == 990000000 ? null : enclosingRegions.Last().backGrid;
            thicknessBackGrid = enclosingRegions.Last().thicknessBackGrid;

            pnJunction = enclosingRegions.Last().pnJunction;
            type = enclosingRegions.Last().type;
            countsAsActiveArea = enclosingRegions.Last().countsAsActiveArea;

            opticFactor_additionalAtmosphereScattering = enclosingRegions.Last().opticFactor_additionalAtmosphereScattering;
            opticFactor_shading = enclosingRegions.Last().opticFactor_shading;
            opticFactor_effectiveArea = enclosingRegions.Last().opticFactor_effectiveArea;
            opticFactor_transmissionGrid = enclosingRegions.Last().opticFactor_transmissionGrid;
            opticFactor_reflection = enclosingRegions.Last().opticFactor_reflection;
            opticFactor_absorption = enclosingRegions.Last().opticFactor_absorption;
            opticFactor_transmission = enclosingRegions.Last().opticFactor_transmission;

            //  ██╗ 
            //  ╚██╗ Set resitors to Edge
            //  ██╔╝
            //  ╚═╝
            resistorsFrontToEdgeContact.Clear();
            resistorsFrontToEdgeGrid.Clear();
            resistorsFrontToEdge.Clear();
            resistorsBackToEdgeContact.Clear();
            resistorsBackToEdgeGrid.Clear();
            resistorsBackToEdge.Clear();

            for (int n = 0; n < neighbors.Count; n++)
            {
                double distance = neighbors[n].distance;
                double edgeLength = neighbors[n].edgeSize;
                if (double.IsNaN(edgeLength) || double.IsInfinity(edgeLength) || edgeLength == 0)
                {
                    resistorsFrontToEdgeContact.Add(double.PositiveInfinity);
                    resistorsFrontToEdgeGrid.Add(double.PositiveInfinity);
                    resistorsFrontToEdge.Add(double.PositiveInfinity);
                    resistorsBackToEdgeContact.Add(double.PositiveInfinity);
                    resistorsBackToEdgeGrid.Add(double.PositiveInfinity);
                    resistorsBackToEdge.Add(double.PositiveInfinity);
                }
                else
                {
                    //  ██╗ front resistances
                    //  ╚═╝
                    // front contact resistance
                    resistorsFrontToEdgeContact.Add(frontContact.propertiesContact.GetSpecificResistivity(thicknessFrontContact) * distance / (2 * thicknessFrontContact * edgeLength));

                    // front grid resistance
                    if (frontGrid == null)
                        resistorsFrontToEdgeGrid.Add(double.NaN);
                    else
                    {
                        if (double.IsNaN(frontGridDensity)) // if in topology optimization
                            resistorsFrontToEdgeGrid.Add(frontGrid.propertiesContact.GetSpecificResistivity(thicknessFrontGrid) * distance / (2 * thicknessFrontGrid * edgeLength));
                        else
                        {
                            if (frontGridDensity == 0)
                                resistorsFrontToEdgeGrid.Add(double.NaN);
                            else
                                resistorsFrontToEdgeGrid.Add(frontGrid.propertiesContact.GetSpecificResistivity(thicknessFrontGrid) * distance / (2 * thicknessFrontGrid * edgeLength) * Misc.SIMPfunction_Conductivity(frontGridDensity));
                        }
                    }

                    // front total resistance
                    if (double.IsNaN(resistorsFrontToEdgeGrid.Last()))
                        resistorsFrontToEdge.Add(resistorsFrontToEdgeContact.Last());
                    else
                        resistorsFrontToEdge.Add(Misc.ParallelResistance(resistorsFrontToEdgeContact.Last(), resistorsFrontToEdgeGrid.Last()));

                    //  ██╗ back resistances
                    //  ╚═╝
                    // back contact resistance
                    resistorsBackToEdgeContact.Add(backContact.propertiesContact.GetSpecificResistivity(thicknessBackContact) * distance / (2 * thicknessBackContact * edgeLength));

                    // back grid and total resistance
                    if (backGrid == null)
                    {
                        resistorsBackToEdgeGrid.Add(double.NaN);
                        resistorsBackToEdge.Add(resistorsBackToEdgeContact.Last());
                    }
                    else
                    {
                        resistorsBackToEdgeGrid.Add(backGrid.propertiesContact.GetSpecificResistivity(thicknessBackGrid) * distance / (2 * thicknessBackGrid * edgeLength));
                        resistorsBackToEdge.Add(Misc.ParallelResistance(resistorsBackToEdgeContact.Last(), resistorsBackToEdgeGrid.Last()));
                    }
                }
            }
        }
        /// <summary>
        /// refreshes the electrical parameters for each single point after all points are set (e.g. total resistor list to all neighbor microcells)
        /// </summary>
        public void RefreshParametersSecond(Mesh<FiniteElementCell> mesh)
        {
            // Reset Lists
            resistorsFront.Clear();
            resistorsBack.Clear();

            // Total Resistors
            for (int n = 0; n < neighbors.Count; n++)
            {
                int indexInNeighborsArray = mesh.finiteElements[neighbors[n].index].neighbors.FindIndex(p => p.index == index);

                // Front resistors
                double Rfront = resistorsFrontToEdge[n] + mesh.finiteElements[neighbors[n].index].resistorsFrontToEdge[indexInNeighborsArray];
                resistorsFront.Add(double.IsNaN(Rfront) ? double.PositiveInfinity : Rfront);

                // Back resistors
                double Rback = resistorsBackToEdge[n] + mesh.finiteElements[neighbors[n].index].resistorsBackToEdge[indexInNeighborsArray];
                resistorsBack.Add(double.IsNaN(Rback) ? double.PositiveInfinity : Rback);
            }
        }
        /// <summary>
        /// returns the derivative of the total resistor to the edge with respect to the design variable x_e
        /// </summary>
        /// <param name="mesh">mesh of the whole cell</param>
        /// <param name="globalNeighborIndex">index of the neighbor (index in mesh, NOT index in neighbors array)</param>
        /// <returns></returns>
        public double ResistorToEdge_DerivativeXe(int globalNeighborIndex, Vector<bool> TOdensityArray_isFixed)
        {
            if (TOdensityArray_isFixed[index])
                return 0;

            int indexInNeighborList = GetIndexInNeighborList(globalNeighborIndex);

            //  ∂R      ∂R     ∂Rgrid
            // ———— = —————— * ——————
            // ∂x_e   ∂Rgrid    ∂x_e

            //   ∂R
            // ——————
            // ∂Rgrid
            if (double.IsNaN(resistorsFrontToEdgeGrid[indexInNeighborList]))
                return 0;
            double delR_delRgrid = Misc.ParallelResistanceDerivation(resistorsFrontToEdgeGrid[indexInNeighborList], resistorsFrontToEdgeContact[indexInNeighborList]);

            // ∂Rgrid
            // ——————
            //  ∂x_e
            double distance = neighbors[indexInNeighborList].distance;
            double edgeLength = neighbors[indexInNeighborList].edgeSize;
            double delRgrid_delXe = 1e10;
            if (frontGridDensity != 0)
                delRgrid_delXe = frontGrid.propertiesContact.GetSpecificResistivity(thicknessFrontGrid) * distance / (2 * thicknessFrontGrid * edgeLength) * Misc.SIMPderivative_Conductivity(frontGridDensity);

            return delR_delRgrid * delRgrid_delXe;
        }

        // Set initial guess ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Sets the initial guess for electrical potential
        /// </summary>
        public void SetInitialGuess((double lower, double upper) operatingVoltage, bool isModule, SimulationSelector simulationSelector)
        {
            if (isModule) // set initial guess for module
            {
                if (type == pointType.cell || type == pointType.P1)
                {
                    phiFrontInit = operatingVoltage.upper;
                    phiBackInit = operatingVoltage.lower;
                }
                else if (type == pointType.gap12 || type == pointType.P2 || type == pointType.gap23)
                {
                    phiFrontInit = operatingVoltage.upper;
                    phiBackInit = operatingVoltage.upper;
                }
                else // P3
                {
                    phiFrontInit = 2 * operatingVoltage.upper - operatingVoltage.lower;
                    phiBackInit = operatingVoltage.upper;
                }
            }
            else // set initial guess for cell
            {
                switch (simulationSelector)
                {
                    case SimulationSelector.bothPotentials:
                        phiFrontInit = operatingVoltage.upper;
                        phiBackInit = operatingVoltage.lower - 0.1;
                        break;

                    case SimulationSelector.frontPotential:
                        phiFrontInit = operatingVoltage.upper + 0.1;
                        phiBackInit = operatingVoltage.lower;
                        break;

                    case SimulationSelector.backPotential:
                        phiFrontInit = operatingVoltage.upper;
                        phiBackInit = operatingVoltage.lower - 0.1;
                        break;
                }
            }

            phiFront = phiFrontInit;
            phiBack = phiBackInit;
        }

        // calculate parameters after solving ███████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// calculate parameters after solving (e.g. currents and current densities to neighbors)
        /// </summary>
        /// <param name="mesh">simulation mesh container</param>
        /// <returns></returns>
        public void CalculateDerivedParameters(Mesh<FiniteElementCell> mesh)
        {
            // set currents —————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            // Reset all currents
            currentsFront.Clear();
            currentdensitiesFront.Clear();
            currentsBack.Clear();
            currentdensitiesBack.Clear();

            // set single currents and currentdensities
            for (int n = 0; n < neighbors.Count; n++)
            {
                // current in Ampere
                double Ifront = (phiFront - mesh.finiteElements[neighbors[n].index].phiFront) / ResistanceFrontTo(mesh.finiteElements[neighbors[n].index]);
                double Iback = (phiBack - mesh.finiteElements[neighbors[n].index].phiBack) / ResistanceBackTo(mesh.finiteElements[neighbors[n].index]);

                // currentdensity normalized over TCO and Grid in A/m²
                double jFront = Ifront / ((thicknessFrontContact + thicknessFrontGrid * (frontGrid == null ? 0 : 1)) * neighbors[n].edgeSize);
                double jBack = Iback / ((thicknessBackContact + thicknessBackGrid * (backGrid == null ? 0 : 1)) * neighbors[n].edgeSize);

                // fill arrays
                currentsFront.Add(Ifront);
                currentdensitiesFront.Add(jFront);
                currentsBack.Add(Iback);
                currentdensitiesBack.Add(jBack);
            }

            // vectorial sum of currents
            //         ___    
            // →       \    →
            // I_vec =  ⟩   I_out[n]
            //         /   
            //         ‾‾‾
            //        n ∈ N
            for (int n = 0; n < neighbors.Count; n++) // sum over all neighbors
            {
                if (currentsFront[n] > 0) // only outflowing currents
                {
                    // angle of the out coming current
                    double angle = Misc.Atan3(position.VectorTo(mesh.finiteElements[neighbors[n].index].position));

                    // get vector from absolut value and angle and add up all out flowing currents
                    IvecFront[0] += currentsFront[n] * Math.Cos(angle);
                    IvecFront[1] += currentsFront[n] * Math.Sin(angle);
                }

                if (currentsBack[n] > 0) // only outflowing currents
                {
                    // angle of the out coming current
                    double angle = Misc.Atan3(position.VectorTo(mesh.finiteElements[neighbors[n].index].position));

                    // get vector from absolut value and angle and add up all out flowing currents
                    IvecBack[0] += currentsBack[n] * Math.Cos(angle);
                    IvecBack[1] += currentsBack[n] * Math.Sin(angle);
                }
            }

            // created heat power ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            if (type == pointType.cell)
            {
                // front
                double powerAtFrontGrid = 0;
                if (frontGrid != null)
                    for (int n = 0; n < neighbors.Count; n++)
                        if (mesh.finiteElements[neighbors[n].index].type == pointType.cell)
                            powerAtFrontGrid += Math.Pow((mesh.finiteElements[neighbors[n].index].phiFront - phiFront) / 2, 2) / resistorsFrontToEdgeGrid[n];
                double powerAtFrontContact = 0;
                for (int n = 0; n < neighbors.Count; n++)
                    if (mesh.finiteElements[neighbors[n].index].type == pointType.cell)
                        powerAtFrontContact += Math.Pow((mesh.finiteElements[neighbors[n].index].phiFront - phiFront) / 2, 2) / resistorsFrontToEdgeContact[n];

                // back
                double powerAtBackGrid = 0;
                if (backGrid != null)
                    for (int n = 0; n < neighbors.Count; n++)
                        if (mesh.finiteElements[neighbors[n].index].type == pointType.cell)
                            powerAtBackGrid += Math.Pow((mesh.finiteElements[neighbors[n].index].phiBack - phiBack) / 2, 2) / resistorsBackToEdgeGrid[n];
                double powerAtBackContact = 0;
                for (int n = 0; n < neighbors.Count; n++)
                    if (mesh.finiteElements[neighbors[n].index].type == pointType.cell)
                        powerAtBackContact += Math.Pow((mesh.finiteElements[neighbors[n].index].phiBack - phiBack) / 2, 2) / resistorsBackToEdgeContact[n];

                // internal (diode and shunt)
                double powerInternal = Math.Abs(pnJunction.GetCurrentAtVoltage(phiFront - phiBack, 0) * (phiFront - phiBack));

                createdHeatPower = powerAtFrontGrid + powerAtFrontContact + powerAtBackGrid + powerAtBackContact + powerInternal;
            }
            else
            {
                createdHeatPower = 0;
            }
        }

        // Gets the resistances to neighbor microcell ███████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Gets the resistance on the front side to neighbor microcell in Ohm
        /// </summary>
        /// <param name="second_point">point, where the current is flowing to</param>
        /// <returns></returns>
        public double ResistanceFrontTo(FiniteElementCell second_point)
        {
            for (int i = 0; i < neighbors.Count; i++)
                if (second_point.index == neighbors[i].index)
                    return resistorsFront[i];

            return 0;
        }
        /// <summary>
        /// Gets the resistance on the back side to neighbor microcell in Ohm
        /// </summary>
        /// <param name="second_point">point, where the current is flowing to</param>
        /// <returns></returns>
        public double ResistanceBackTo(FiniteElementCell second_point)
        {
            for (int i = 0; i < neighbors.Count; i++)
                if (second_point.index == neighbors[i].index)
                    return resistorsBack[i];

            return 0;
        }

        // Calculate the current, which is produced by the microcell (minus diode and shunt) ████████████████████████████████████████████████████████
        /// <summary>
        /// Calculate the current, which is produced by the microcell (minus diode and shunt) (negative = current is generated)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentGenerated(double? voltage = null, double specialFrontGridDensity = double.NaN, double? specialRshunt = null)
        {
            if (type == pointType.P1 || type == pointType.P2 || type == pointType.P3)
                return 0;

            return pnJunction.GetCurrentAtVoltage(voltage ?? (phiFront - phiBack), OpticalPrefactor(specialFrontGridDensity) * pnJunction.characteristicCurve.currentPhoto,
                null, null, null, specialRshunt) * size;
        }
        /// <summary>
        /// Calculate the derivative according to the front potential, which is produced by the microcell (minus diode and shunt)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentGenerated_DerivativePhiFront(double? voltage = null, double specialFrontGridDensity = double.NaN)
        {
            if (type == pointType.P1 || type == pointType.P2 || type == pointType.P3)
                return 0;

            return pnJunction.GetCurrentAtVoltageDerivation(voltage ?? (phiFront - phiBack), OpticalPrefactor(specialFrontGridDensity) * pnJunction.characteristicCurve.currentPhoto) * size;
        }
        /// <summary>
        /// Calculate the derivative according to the back potential, which is produced by the microcell (minus diode and shunt)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentGenerated_DerivativePhiBack(double? voltage = null, double specialFrontGridDensity = double.NaN)
        {
            if (type == pointType.P1 || type == pointType.P2 || type == pointType.P3)
                return 0;

            return -pnJunction.GetCurrentAtVoltageDerivation(voltage ?? (phiFront - phiBack), OpticalPrefactor(specialFrontGridDensity) * pnJunction.characteristicCurve.currentPhoto) * size;
        }
        /// <summary>
        /// Calculate the derivative according to the design variable x_e, which is produced by the microcell (minus diode and shunt)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentGenerated_DerivativeXe(double? voltage = null, double specialFrontGridDensity = double.NaN)
        {
            if (type == pointType.P1 || type == pointType.P2 || type == pointType.P3)
                return 0;

            // Igen = I (Φ, Iph(x_e))
            // Iph = const * (tGrid + tTO)
            // tTO = (1 - tGrid) * (1 - x_e)^3

            // ∂Igen   ∂Igen   ∂Iph   ∂tTO
            // ————— = ————— * ———— * ————
            // ∂x_e    ∂Iph    ∂tTO   ∂x_e

            // ∂tTO/∂x_e
            double useFrontGridDensity = double.IsNaN(specialFrontGridDensity) ? frontGridDensity : specialFrontGridDensity;
            double delAdditionalTO_delXe = double.IsNaN(useFrontGridDensity) ? 0 : (1 - opticFactor_transmissionGrid) * Misc.SIMPderivative_GeneratedCurrent(frontGridDensity);

            // ∂Iph/∂tTO
            double delIph_delAdditionalTO = (1 - opticFactor_additionalAtmosphereScattering)
                * (1 - opticFactor_shading)
                * opticFactor_effectiveArea
                // * 1 ∂tTO / ∂tTO
                * (1 - opticFactor_reflection - opticFactor_absorption.Where(m => !m.isAbsorber).Sum(m => m.factor) - opticFactor_transmission)
                * pnJunction.characteristicCurve.currentPhoto;

            // ∂Igen/∂Iph
            double delIgen_delIph = pnJunction.characteristicCurve.GetCurrentAtVoltageDerivation_Iph(voltage ?? (phiFront - phiBack), OpticalPrefactor(specialFrontGridDensity) * pnJunction.characteristicCurve.currentPhoto);

            return delIgen_delIph * delIph_delAdditionalTO * delAdditionalTO_delXe * size;
        }
        /// <summary>
        /// Calculate Optical prefactor multiplied with Iph
        /// </summary>
        /// <returns></returns>
        public double OpticalPrefactor(double specialFrontGridDensity)
        {
            double useFrontGridDensity = double.IsNaN(specialFrontGridDensity) ? frontGridDensity : specialFrontGridDensity;
            double additionalTO = double.IsNaN(useFrontGridDensity) ? 0 : (1 - opticFactor_transmissionGrid) * Misc.SIMPfunction_GeneratedCurrent(frontGridDensity);

            double opticalPrefactor = 1;

            opticalPrefactor *= 1 - opticFactor_additionalAtmosphereScattering;
            opticalPrefactor *= 1 - opticFactor_shading;
            opticalPrefactor *= opticFactor_effectiveArea;
            opticalPrefactor *= opticFactor_transmissionGrid + additionalTO;
            opticalPrefactor *= 1 - opticFactor_reflection - opticFactor_absorption.Where(m => !m.isAbsorber).Sum(m => m.factor) - opticFactor_transmission;

            return opticalPrefactor;
        }

        // Calculate the current, which flows laterally through a cell at the P2-gap ████████████████████████████████████████████████████████████████
        /// <summary>
        /// Calculate the current, which flows laterally through a cell at the P2-gap (positive = current flows from front to back)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentP2(double? voltage = null)
        {
            if (type != pointType.P2)
                return 0;

            // two resitors in series: path in frontcontact material and contact resistor
            return (voltage ?? (phiFront - phiBack)) / (frontContact.propertiesContact.GetSpecificResistivity(thicknessFrontContact) * pnJunction.thicknessAbsorberLayer / size + contactSeriesResistanceP2);
        }
        /// <summary>
        /// Calculate the derivation according to the front potential of th current, which flows laterally through a cell at the P2-gap
        /// </summary>
        /// <returns></returns>
        public double GetCurrentP2_DerivativePhiFront()
        {
            return 1 / (frontContact.propertiesContact.GetSpecificResistivity(thicknessFrontContact) * pnJunction.thicknessAbsorberLayer / size + contactSeriesResistanceP2);
        }
        /// <summary>
        /// Calculate the derivation according to the back potential of th current, which flows laterally through a cell at the P2-gap
        /// </summary>
        /// <returns></returns>
        public double GetCurrentP2_DerivativePhiBack()
        {
            return -1 / (frontContact.propertiesContact.GetSpecificResistivity(thicknessFrontContact) * pnJunction.thicknessAbsorberLayer / size + contactSeriesResistanceP2);
        }

        // Calculate the sum of the currents, which flows to the next neighbors (positive = outgoing) ███████████████████████████████████████████████
        /// <summary>
        /// Calculates the sum of the currents, which flows to the next neighbors on the front side (positive = outgoing)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentToNeighborsFrontSum(Mesh<FiniteElementCell> mesh, List<(int index, double distance, double edgeSize)> neighbors = null)
        {
            // select neighbors
            List<(int index, double distance, double edgeSize)> N = neighbors ?? this.neighbors;

            // current to neighbors via TCO and grid
            //
            //                       ___
            //                       \    Φ_k - Φ_n
            //           I_k,n  =     ⟩   ——————————
            //                       /     R_{k,n}
            //                       ‾‾‾
            //                     n ∈ N(k)
            //
            double currentSum = 0;
            for (int n = 0; n < N.Count; n++)
                currentSum += (phiFront - mesh.finiteElements[N[n].index].phiFront) / ResistanceFrontTo(mesh.finiteElements[N[n].index]);
            return currentSum;
        }
        /// <summary>
        /// Calculates the sum of the currents, which flows to the next neighbors on the back side (positive = outgoing)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentToNeighborsBackSum(Mesh<FiniteElementCell> mesh, List<(int index, double distance, double edgeSize)> neighbors = null)
        {
            // select neighbors
            List<(int index, double distance, double edgeSize)> N = neighbors ?? this.neighbors;

            // current to neighbors via TCO and grid
            //
            //                       ___
            //                       \    Φ_k - Φ_n
            //           I_k,n  =     ⟩   ——————————
            //                       /     R_{k,n}
            //                       ‾‾‾
            //                     n ∈ N(k)
            //
            double currentSum = 0;
            for (int n = 0; n < N.Count; n++)
                currentSum += (phiBack - mesh.finiteElements[N[n].index].phiBack) / ResistanceBackTo(mesh.finiteElements[N[n].index]);
            return currentSum;
        }
        /// <summary>
        /// Calculates the single derivative of the sum of the currents, which flows to the next neighbors (positive = outgoing)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentToNeighborsDerivativeFrontSingle(Mesh<FiniteElementCell> mesh, int neighbor)
        {
            // current to neighbors via TCO and grid
            //
            //                          1
            //           ∂I/∂Φ_n  =  —————————
            //                        R_{k,n}
            //
            return 1 / ResistanceFrontTo(mesh.finiteElements[neighbor]);
        }
        /// <summary>
        /// Calculates the single derivative of the sum of the currents, which flows to the next neighbors (positive = outgoing)
        /// </summary>
        /// <returns></returns>
        public double GetCurrentToNeighborsDerivativeBackSingle(Mesh<FiniteElementCell> mesh, int neighbor)
        {
            // current to neighbors via TCO and grid
            //
            //                          1
            //           ∂I/∂Φ_n  =  —————————
            //                        R_{k,n}
            //
            return 1 / ResistanceBackTo(mesh.finiteElements[neighbor]);
        }

        // write point to file ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// writes point to file
        /// </summary>
        /// <param name="file">file, where the line is written to</param>
        public void OutputToFile(StreamWriter file)
        {
            StringBuilder builder = new StringBuilder();
            builder.Append(InputOutput.ToStringWithSeparator(position.x));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(position.y));
            builder.Append("\t" + index);
            builder.Append("\t" + InputOutput.ToStringWithSeparator(phiFront));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(phiBack));
            if (frontGrid == null)
                builder.Append("\t0");
            else
                builder.Append("\t1");
            if (backGrid == null)
                builder.Append("\t0");
            else
                builder.Append("\t1");
            builder.Append("\t" + InputOutput.ToStringWithSeparator(size));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(Misc.Atan3(IvecFront)));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(Misc.GetPNormOfVector(IvecFront, 2)));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(Misc.Atan3(IvecBack)));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(Misc.GetPNormOfVector(IvecBack, 2)));
            builder.Append("\t" + InputOutput.ToStringWithSeparator(createdHeatPower));

            file.WriteLine(builder.ToString());
        }

        // Print to console █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Print point to console
        /// </summary>
        /// <param name="mesh">simulation container</param>
        public void Print(Mesh<FiniteElementCell> mesh)
        {
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("point index: " + index);
            Console.ForegroundColor = ConsoleColor.Gray;

            // Geometry
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("    geometry");
            Console.ForegroundColor = ConsoleColor.Gray;

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tposition in m: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine("(" + position.x + "|" + position.y + ")");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tarea: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(size + " m²");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.WriteLine("\t" + corners.Count + " corners in m:");
            Console.ForegroundColor = ConsoleColor.Gray;
            foreach (var corner in corners)
                Console.WriteLine("\t  (" + corner.position.x + "|" + corner.position.y + ")");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.WriteLine("\t" + neighbors.Count + " neighbors, edge lengths, distances and currents on front side:");
            Console.ForegroundColor = ConsoleColor.Gray;
            for (int n = 0; n < neighbors.Count; n++)
                Console.WriteLine("\t    (" + neighbors[n] + ") -> " + neighbors[n].edgeSize + " m -> " + neighbors[n].distance + " m -> "
                    + ((phiFront - mesh.finiteElements[neighbors[n].index].phiFront) / ResistanceFrontTo(mesh.finiteElements[neighbors[n].index])) + " A");

            // Preferences
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("    preferences");
            Console.ForegroundColor = ConsoleColor.Gray;

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\ttype of point: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(type.ToString());

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tFrontGrid: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(frontGrid == null ? "no front grid" : (thicknessFrontGrid * 1e9 + "nm " + frontGrid.name));

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tFrontContact: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(thicknessFrontContact * 1e9 + "nm " + frontContact.name);

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tBackContact: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(thicknessBackContact * 1e9 + "nm " + backContact.name);

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tBackGrid: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(backGrid == null ? "no back grid" : (thicknessBackGrid * 1e9 + "nm " + backGrid.name));

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tisExternalCellFrontContact: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(isExternalCellFrontContact);

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tisExternalCellBackContact: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(isExternalCellBackContact);

            // Electrics
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("    electrics");
            Console.ForegroundColor = ConsoleColor.Gray;

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tPhiFront: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(phiFront + " V");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tPhiBack: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(phiBack + " V");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tinitial PhiFront: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(phiFrontInit + " V");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tinitial PhiBack: ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(phiBackInit + " V");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tcurrent to neighbors on front (neg = incoming): ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(currentsFront.Sum() + " A");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tcurrent to neighbors on back (neg = incoming): ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(currentsBack.Sum() + " A");

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.Write("\tgenerated current: (neg = generated): ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine(GetCurrentGenerated() + " A");

            Console.WriteLine();
        }
    }
}