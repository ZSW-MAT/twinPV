using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Geometry;
using BasicLib;
using Database;
using Extreme.Mathematics;
using MoreLinq;
using Extreme.Mathematics.Curves;

namespace Semiconductor
{
    /// <summary>
    /// Class, with respresents a subclass of Point for the Semiconductor simulation
    /// </summary>
    public class FiniteElementSemiconductor : FiniteElement
    {
        /// <summary>
        /// elmaterial.chemicalPotentialtrical potential in V
        /// </summary>
        public double phi { get; set; }
        /// <summary>
        /// initial guess of the elctrical potential in V
        /// </summary>
        public double phiInit { get; set; }
        /// <summary>
        /// elmaterial.chemicalPotentialtron quasi fermi level
        /// </summary>
        public double phi_n { get; set; }
        /// <summary>
        /// hole quasi fermi level
        /// </summary>
        public double phi_p { get; set; }

        /// <summary>
        /// initial guess of electron fermi potential (0)
        /// </summary>
        public double phi_n_init { get; set; }
        /// <summary>
        /// initial guess of hole fermi potential (0)
        /// </summary>
        public double phi_p_init { get; set; }

        public (double phi_BU, double phi_n_BU, double phi_p_BU) backUpPotentials { get; set; }


        /// <summary>
        /// electron density in thermodynamic equilibrium, is set after solving Poisson
        /// </summary>
        public double nEquilibrium { get; set; }
        /// <summary>
        /// hole density in thermodynamic equilibrium, is set after solving Poisson
        /// </summary>
        public double pEquilibrium { get; set; }

        public bool isFEinAbsorber { get; set; }

        /// <summary>
        /// material of this point
        /// </summary>
        public Material material { get; set; }

        public double localBandGap { get; set; }

        // external contacts
        /// <summary>
        /// determines, if this point is an external cell front contact
        /// </summary>
        public bool hasBoundaryCondition { get; set; } = false;
        /// <summary>
        /// determines, if this point is an external cell back contact
        /// </summary>
        public double boundaryCondition { get; set; } = 0;
        /// <summary>
        /// determines, if at this point operating voltage is applied (Vop = p contact, because positive bias voltages result in a IV curve)
        /// </summary>
        public bool hasOperatingVoltage { get; set; } = false;
        /// <summary>
        /// dcurrent boundary condition for ramp up of voltage
        /// </summary>
        public double currentBoundaryCondition { get; set; } = 0;

        public bool hasInterfaceCondition { get; set; } = false;

        public (double SRV_electrons, double SRV_holes, double contactBarrier, double interfaceTrapEnergy) contactPreferences { get; set; }
        public (double SRV_electrons, double SRV_holes, double contactBarrier, double interfaceTrapEnergy) backupContactPreferences { get; set; }

        /// <summary>
        /// gives the height of the barrier at the contacts
        /// </summary>
        public double barrierHeight { get; set; }

        /// <summary>
        /// defines the minimal difference of two neighbouring points in Ec or Ev when TE current equations should be used
        /// </summary>
        double barrierTolerancTEcurrent = 1e-2;
        /// <summary>
        /// list of type of electron current calculation (SG, TE,...)  to all neighbors (ordered in the same odrder as neighbor array)
        /// </summary>
        public List<TypeOfCurrent> typeOfElectronCurrents { get; private set; } = new List<TypeOfCurrent>();
        /// <summary>
        /// list of type of hole current calculation (SG, TE,...)  to all neighbors (ordered in the same odrder as neighbor array)
        /// </summary>
        public List<TypeOfCurrent> typeOfHoleCurrents { get; private set; } = new List<TypeOfCurrent>();

        // Currents (magnitude, no vectors) to single neighbors
        /// <summary>
        /// list of electron currents  to all neighbors in Ampere/m^2 (ordered in the same odrder as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> electronCurrent { get; private set; } = new List<double>();
        /// <summary>
        /// list of hole currents  to all neighbors in Ampere/m^2 (ordered in the same odrder as neighbor array)
        /// (positive = outgoing, negative = incoming)
        /// </summary>
        public List<double> holeCurrent { get; private set; } = new List<double>();


        // Current (absolute) vectors to single neighbors
        /// <summary>
        /// List of arrays/vectors where the entries of the list give the electron currents to the neighbors of the point.  In the same order as neighbors.
        /// </summary>
        public List<double[]> electronCurrentsToNeighbors { get; private set; } = new List<double[]>();
        /// <summary>
        /// List of arrays/vectors where the entries of the list give the hole currents to the neighbors of the point.  In the same order as neighbors.
        /// </summary>
        public List<double[]> holeCurrentsToNeighbors { get; private set; } = new List<double[]>();

        //Current Density vectors to single neighbors
        /// <summary>
        /// List of arrays/vectors where the entries of the list give the electron current densities to the neighbors of the point.  In the same order as neighbors.
        /// </summary>
        public List<double[]> electronCurrentDensitiesToNeighbors { get; private set; } = new List<double[]>();
        /// <summary>
        /// List of arrays/vectors where the entries of the list give the hole current densities to the neighbors of the point.  In the same order as neighbors.
        /// </summary>
        public List<double[]> holeCurrentDensitiesToNeighbors { get; private set; } = new List<double[]>();

        //Resulting current vectors
        /// <summary>
        /// Vector which represents the resulting electron current (after summing up all single currents to the neighbors), first component x part and second component y part
        /// </summary>
        public double[] electronCurrentVector { get; private set; } = new double[2];
        /// <summary>
        /// Vector which represents the resulting hole current (after summing up all single currents to the neighbors), first component x part and second component y part
        /// </summary>
        public double[] holeCurrentVector { get; private set; } = new double[2];

        //currents to boundaries
        /// <summary>
        /// saves the electron current to a boundary 
        /// </summary>
        public double electronCurrentToBoundary { get; set; } = 0;
        /// <summary>
        /// saves the hole current to a boundary 
        /// </summary>
        public double holeCurrentToBoundary { get; set; } = 0;

        /// <summary>
        /// Sum over all e and h currents (Sum(jn+jp)) to/from all neighbors. Should be 0 for every single point
        /// </summary>
        public double sumOfCurrents { get; set; } = 0;


        // Recombination Currents for loss analysis
        public List<(double voltage, double localGenerationCurrent)> voltageDependentLocalGenerationArray { get; set; } = new List<(double voltage, double localGenerationCurrent)>();
        public List<(double voltage, double localSRHRecombinationCurrent)> voltageDependentLocalSRHRecombinationArray { get; set; } = new List<(double voltage, double localSRHRecombinationCurrent)>();
        public List<(double voltage, double localAugerRecombinationCurrent)> voltageDependentLocalAugerRecombinationArray { get; set; } = new List<(double voltage, double localAugerRecombinationCurrent)>();
        public List<(double voltage, double localRadiativeRecombinationCurrent)> voltageDependentLocalRadiativeRecombinationArray { get; set; } = new List<(double voltage, double localRadiativeRecombinationCurrent)>();
        public List<(double voltage, double localSRVelectronVop)> voltageDependentLocalSRVelectronVop { get; set; } = new List<(double voltage, double localSRVelectronVop)>();
        public List<(double voltage, double localSRVelectronVzero)> voltageDependentLocalSRVelectronVzero { get; set; } = new List<(double voltage, double localSRVelectronVzero)>();
        public List<(double voltage, double localSRVholeVop)> voltageDependentLocalSRVholeVop { get; set; } = new List<(double voltage, double localSRVholeVop)>();
        public List<(double voltage, double localSRVholeVzero)> voltageDependentLocalSRVholeVzero { get; set; } = new List<(double voltage, double localSRVholeVzero)>();
        public List<(double voltage, double localSRVholeVzero)> voltageDependentLocalInterfaceRecombinationArray { get; set; } = new List<(double voltage, double localInterfaceRecombination)>();


        //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



        //Incomplete Ionization settings
        public double donorEnergy { get; set; } = 0.05;
        public double acceptorEnergy { get; set; } = 0.05;
        public double degeneracyFactorDonor { get; set; } = 2;
        public double degeneracyFactorAcceptor { get; set; } = 4;


        // MISC functions-----------------------------------------------------------------------------------------------
        /// <summary>
        /// writes point to file
        /// </summary>
        /// <param name="mesh">mesh container</param>
        /// <param name="file">file where the point gets written to</param>
        /// <param name="offset">Offset, if the values need to be shifted</param>
        public void OutputToFile(Mesh<FiniteElementSemiconductor> mesh, StreamWriter file, ModelSemiconductor modelSemiconductor)
        {
            // Potentials and bands
            double Ec = -phi + material.propertiesSemiconductor.chemicalPotential;

            double Ev = -phi + material.propertiesSemiconductor.chemicalPotential - localBandGap;

            double phin = phi_n;

            double phip = phi_p;


            string datastring = InputOutput.ToStringWithSeparator(position.x * 1e6) + "\t" +
                               InputOutput.ToStringWithSeparator(position.y * 1e6) + "\t" +
                               InputOutput.ToStringWithSeparator(Ec) + "\t" +
                               InputOutput.ToStringWithSeparator(Ev) + "\t" +
                               InputOutput.ToStringWithSeparator(phin) + "\t" +
                               InputOutput.ToStringWithSeparator(phip) + "\t" +
                               InputOutput.ToStringWithSeparator(-phi) + "\t" +
                               InputOutput.ToStringWithSeparator(nDensity(modelSemiconductor)) + "\t" +
                               InputOutput.ToStringWithSeparator(pDensity(modelSemiconductor)) + "\t"
                               ;

            file.WriteLine(datastring);
        }

        /// <summary>
        /// sets the i-th differential equation variable (i=0: phi)
        /// </summary>
        /// <param name="differentialVariable">value which will be set to the variable of the point</param>
        /// <param name="index">index of the differential variable, which will be set</param>
        public override void SetDifferentialEquationVariable(double differentialVariable, int index)
        {
            switch (index)
            {
                case 0:
                    phi = differentialVariable;
                    break;
                default:
                    break;
            }
        }

        /// <summary>
        /// gets all differential equation variables listed in a double array (double[] { phi })
        /// </summary>
        /// <returns></returns>
        public override double[] GetDifferentialEquationVariables()
        {
            return new double[] { phi };
        }

        /// <summary>
        /// refreshes all elmaterial.chemicalPotentialtric properties
        /// </summary>
        public void RefreshParameters(ModelSemiconductor semiconductor)
        {
            //  ██╗ 
            //  ╚██╗ Set electronic data
            //  ██╔╝
            //  ╚═╝
            var enclosingRegions = semiconductor.meshingAlgorithm.GetEnclosingRegions(position);
            for (int r = 0; r < enclosingRegions.Count; r++)
            {
                material = enclosingRegions[r].material;
                localBandGap = material.propertiesSemiconductor.Egap;

                if (!double.IsNaN( enclosingRegions[r].gradingCoordinates[0].x ))
                {
                    enclosingRegions[r].isGradedEgap = true;
                }
                else
                    enclosingRegions[r].isGradedEgap = false;

            }


            if (semiconductor.meshingAlgorithm.dimension == 2)
            { 
                if (enclosingRegions.Last().isGradedEgap)
                {
                    //TODO: Input x/ y direction

                    double minX = enclosingRegions.Last().orderedPoints.Min(p => p.position.x);
                    double maxX = enclosingRegions.Last().orderedPoints.Max(p => p.position.x);
                    double d_layer = maxX - minX;

                    // Spline Calculation:
                    var xValues = enclosingRegions.Last().gradingCoordinates.Select(d => d.x).ToArray();
                    var yValues = enclosingRegions.Last().gradingCoordinates.Select(d => d.y).ToArray();
                    var GradingSpline = new CubicSpline(xValues, yValues, CubicSplineKind.Akima);

                    localBandGap = enclosingRegions.Last().material.propertiesSemiconductor.Egap * GradingSpline.ValueAt( (position.x - minX) / d_layer ); 
                }
            }
            if (semiconductor.meshingAlgorithm.dimension == 1)
            {
                if (enclosingRegions.Last().isGradedEgap)
                {    
                    double minX = enclosingRegions.Last().orderedPoints.Min(p => p.position.x);
                    double maxX = enclosingRegions.Last().orderedPoints.Max(p => p.position.x);

                    double d_layer = maxX - minX;

                    // Spline Calculation:
                    var xValues = enclosingRegions.Last().gradingCoordinates.Select(d => d.x).ToArray();
                    var yValues = enclosingRegions.Last().gradingCoordinates.Select(d => d.y).ToArray();
                    var GradingSpline = new CubicSpline(xValues, yValues, CubicSplineKind.Akima);

                    localBandGap = enclosingRegions.Last().material.propertiesSemiconductor.Egap * GradingSpline.ValueAt((position.x - minX) / d_layer);
                    

                    //localBandGap = enclosingRegions.Last().material.propertiesSemiconductor.Egap + linearCigsGrading / d_CIGS * (position.x - minX);
                }
            }


        }

        /// <summary>
        /// Modify paramters for gradings etc.
        /// </summary>
        /// <returns></returns>
        public void ModifyParameters(ModelSemiconductor semiconductor)
        {

        }

        /// <summary>
        /// sets the type of current between a point and its neighbors ( in a list of neums in same order as the list of neighbors) 
        /// </summary>
        /// <param name="semiconductor"></param>
        public void RefreshParametersSecond(ModelSemiconductor semiconductor)
        {
            foreach (var p in neighbors)
            {
                //double diffEv = Math.Abs((semiconductor.delaunayVoronoi.mesh.points[p].material.propertiesSemiconductor.Egap - material.propertiesSemiconductor.Egap)
                //                - (semiconductor.delaunayVoronoi.mesh.points[p].material.propertiesSemiconductor.chemicalPotential - material.propertiesSemiconductor.chemicalPotential));
                //double diffEc = Math.Abs(material.propertiesSemiconductor.chemicalPotential - semiconductor.delaunayVoronoi.mesh.points[p].material.propertiesSemiconductor.chemicalPotential);

                double diffEv = Math.Abs(GetDiffEvToNeighbor(semiconductor.mesh.finiteElements[p.index]));
                double diffEc = Math.Abs(GetDiffEcToNeighbor(semiconductor.mesh.finiteElements[p.index]));
                //Console.WriteLine(" Von Punkt: " + index + " nach Punkt: " + semiconductor.delaunayVoronoi.mesh.points[p].index + " diffEv: " + diffEv);
                //Console.WriteLine(" Von Punkt: " + index + " nach Punkt: " + semiconductor.delaunayVoronoi.mesh.points[p].index + " diffEc: " + diffEc);


                if (diffEc > barrierTolerancTEcurrent)
                    typeOfElectronCurrents.Add(TypeOfCurrent.ThermionicEmission);
                else
                    typeOfElectronCurrents.Add(TypeOfCurrent.ScharfetterGummel);

                if (diffEv > barrierTolerancTEcurrent)
                    typeOfHoleCurrents.Add(TypeOfCurrent.ThermionicEmission);
                else
                    typeOfHoleCurrents.Add(TypeOfCurrent.ScharfetterGummel);


            }
        }

        /// <summary>
        /// write the point to the console
        /// </summary>
        public void Print()
        {
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("point index: " + index);
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.WriteLine("\tposition: (" + position.x + "|" + position.y + ")");
            Console.WriteLine("\tvolume: " + size);
            Console.WriteLine("\tphi: " + phi);
            Console.WriteLine("\tphiInit: " + phiInit);
            Console.WriteLine("\tphi_n: " + phi_n);
            Console.WriteLine("\tphi_p: " + phi_p);

            Console.WriteLine("\tSum of Currents: " + sumOfCurrents);
            if (hasBoundaryCondition)
                Console.WriteLine("\tboundaryCondition: " + boundaryCondition);

            Console.ForegroundColor = ConsoleColor.DarkCyan;
            Console.WriteLine("\t" + neighbors.Count + " neighbors, edge lengths, distances");
            Console.ForegroundColor = ConsoleColor.Gray;
            for (int n = 0; n < neighbors.Count; n++)
                Console.WriteLine("\t    (" + neighbors[n] + ") -> " + neighbors[n].edgeSize + " m -> " + neighbors[n].distance + " m");


            //Console.WriteLine("\t" + neighbors.Count + " neighbors and edge lengths and currents and type of currents (neighbor, edgelength, jn, jp, type jn, type jp: ");
            //for (int k = 0; k < neighbors.Count; k++)
            //Console.WriteLine("\t  (" + neighbors[k] + ") -> " + edgeLengths[k] + " -> " + electronCurrent[k] + " -> " + holeCurrent[k] + " -> " + typeOfElectronCurrents[k] + " -> " + typeOfHoleCurrents[k]);
            /*
            Console.WriteLine("\t" + corners.Count + " corners: ");
            foreach (Position2D Eck in corners)
                Console.WriteLine("\t  (" + Eck.x + "|" + Eck.y + ") ");
            Console.WriteLine();*/
        }



        // Differnet physics functions-----------------------------------------------------------------------------------------
        /// <summary>
        /// calculates the inital guess for this point (via elmaterial.chemicalPotentialtrical neutrality)
        /// </summary>
        public void SetInitialGuessFiniteElement(ModelSemiconductor modelSemiconductor)
        {
            // Set initial guess for elmaterial.chemicalPotentialtrical potential
            phiInit = material.propertiesSemiconductor.chemicalPotential - localBandGap / 2
                - (physConstants.kB * modelSemiconductor.T) / (2 * physConstants.e)
                * Math.Log(material.propertiesSemiconductor.Nc(modelSemiconductor.T) / material.propertiesSemiconductor.Nv(modelSemiconductor.T))
                    + physConstants.kB * modelSemiconductor.T / physConstants.e
                    //* Misc.Asinh(Doping(modelSemiconductor) / (2 * material.propertiesSemiconductor.Nintr(modelSemiconductor.T)));
                    * Misc.Asinh((material.propertiesSemiconductor.NDplus - material.propertiesSemiconductor.NAminus) / (2 * Nintr(modelSemiconductor.T))); // only use NaMinus and NdPlus in inital Guess
            phi = phiInit;
        }

        /// <summary>
        /// Calculates the first derivative in x-dirmaterial.chemicalPotentialtion of the potential (-> elmaterial.chemicalPotentialtric field)
        /// </summary>
        /// <param name="mesh">mesh container</param>
        /// <returns></returns>
        public double ElectricField(Mesh<FiniteElementSemiconductor> mesh)
        {
            // returns the derivation of ∂Φ/∂x in [(material.chemicalPotential - material.material.Egap) / m]

            // numerator
            double numerator = 0;
            foreach (var l in neighbors)
                numerator += (mesh.finiteElements[l.index].phi - phi) * l.edgeSize / position.DistanceTo(mesh.finiteElements[l.index].position);

            // denominator
            double denominator = 0;
            foreach (var l in neighbors)
                denominator += l.edgeSize;

            return numerator / denominator;
        }

        /// <summary>
        /// retrurns the difference of conductions band edges of two points
        /// </summary>
        /// <param name="neighbor"></param>
        /// <returns></returns>
        public double GetDiffEcToNeighbor(FiniteElementSemiconductor neighbor)
        {
            return neighbor.material.propertiesSemiconductor.chemicalPotential - material.propertiesSemiconductor.chemicalPotential;
        }

        /// <summary>
        /// retrurns the difference of valence band edges of two points
        /// </summary>
        /// <param name="neighbor"></param>
        /// <returns></returns>
        public double GetDiffEvToNeighbor(FiniteElementSemiconductor neighbor)
        {
            return neighbor.material.propertiesSemiconductor.chemicalPotential - neighbor.localBandGap - (material.propertiesSemiconductor.chemicalPotential - localBandGap);
        }




        /// <summary>
        /// Sets electron and hole currents for every point: absoulute and densities, in  total and to each neighbor
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        public void SetLocalCurrents(ModelSemiconductor modelSemiconductor)
        {
            // Reset all currents
            electronCurrent.Clear();
            holeCurrent.Clear();

            // set single currents and currentdensities —————————————————————————————————————————————————————————————————————————————————————————————
            for (int n = 0; n < neighbors.Count; n++)
            {
                electronCurrent.Add(modelSemiconductor.GetElectronCurrentToNeighbor(index, neighbors[n].index));

                holeCurrent.Add(modelSemiconductor.GetHoleCurrentToNeighbor(index, neighbors[n].index));
            }            

            for (int n = 0; n < neighbors.Count; n++)
            {
                double angle = Misc.Atan3(position.VectorTo(modelSemiconductor.mesh.finiteElements[neighbors[n].index].position));

                //absolute currents to neighbors
                electronCurrentsToNeighbors.Add(new double[] { electronCurrent[n] * Math.Cos(angle), electronCurrent[n] * Math.Sin(angle) });
                holeCurrentsToNeighbors.Add(new double[] { holeCurrent[n] * Math.Cos(angle), holeCurrent[n] * Math.Sin(angle) });

                // Current DENSITIES to neighbors (devided by edgelength)
                electronCurrentDensitiesToNeighbors.Add(new double[] { electronCurrent[n] * Math.Cos(angle) / neighbors[n].edgeSize, electronCurrent[n] * Math.Sin(angle) / neighbors[n].edgeSize });
                holeCurrentDensitiesToNeighbors.Add(new double[] { holeCurrent[n] * Math.Cos(angle) / neighbors[n].edgeSize, holeCurrent[n] * Math.Sin(angle) / neighbors[n].edgeSize });

                // Total currents as vectors
                electronCurrentVector[0] += electronCurrent[n] * Math.Cos(angle);
                electronCurrentVector[1] += electronCurrent[n] * Math.Sin(angle);
                holeCurrentVector[0] += holeCurrent[n] * Math.Cos(angle);
                holeCurrentVector[1] += holeCurrent[n] * Math.Sin(angle);

                sumOfCurrents += (electronCurrent[n] + holeCurrent[n]);

            }
            if (hasBoundaryCondition)
            {
                Position perpendicularPoint = null;
                switch (modelSemiconductor.meshingAlgorithm.dimension)
                {
                    case 1:
                        perpendicularPoint = modelSemiconductor.meshingAlgorithm.contourJunctions[indexOfBorderElementCreatedFrom].position;
                        break;
                    case 2:
                        perpendicularPoint = position.PerpendicularPoint(modelSemiconductor.meshingAlgorithm.contourSegments[indexOfBorderElementCreatedFrom].lineSegment);
                        break;
                }
                var vectorToBoundary = position.VectorTo(perpendicularPoint);
                var angle = Misc.Atan3(vectorToBoundary);

                electronCurrentToBoundary = ElectronSurfaceRecombinationCurrent(modelSemiconductor) * borderEdgeSize;
                holeCurrentToBoundary = HoleSurfaceRecombinationCurrent(modelSemiconductor) * borderEdgeSize;
                //Console.WriteLine("e " + position.x + "\t" + electronCurrentToBoundary);
                //Console.WriteLine("h " + position.x + "\t"+ holeCurrentToBoundary);
                electronCurrentVector[0] += electronCurrentToBoundary * Math.Cos(angle);
                electronCurrentVector[1] += electronCurrentToBoundary * Math.Sin(angle);
                holeCurrentVector[0] += holeCurrentToBoundary * Math.Cos(angle);
                holeCurrentVector[1] += holeCurrentToBoundary * Math.Sin(angle);

                sumOfCurrents += (electronCurrentToBoundary + holeCurrentToBoundary);
            }
        }


        // density and doping functions-----------------------------------------------------------------------------------
        /// <summary>
        /// total doping (NDplus - NAminus) with ionized impurities 
        /// </summary>
        /// <returns></returns>
        public double Doping(ModelSemiconductor modelSemiconductor)
        {

            double referenceDesnityDonor = (material.propertiesSemiconductor.Nc(modelSemiconductor.T) * Math.Exp(donorEnergy / modelSemiconductor.U_T));
            double referenceDesnityAcceptor = (material.propertiesSemiconductor.Nv(modelSemiconductor.T) * Math.Exp(-acceptorEnergy / modelSemiconductor.U_T));

            double ionizedDonorDensity = material.propertiesSemiconductor.NDplus / (1 + degeneracyFactorDonor * nDensity(modelSemiconductor) / referenceDesnityDonor);
            double ionizedAcceptorDensity = material.propertiesSemiconductor.NAminus / (1 + degeneracyFactorAcceptor * pDensity(modelSemiconductor) / referenceDesnityAcceptor);

            /*
            Console.WriteLine("\n");
            Console.WriteLine("Index: " + index);
            Console.WriteLine("NDplus alt: " + material.propertiesSemiconductor.NDplus);
            Console.WriteLine("ND neu: " + ionizedDonorDensity );
            Console.WriteLine("NAminus alt: " + material.propertiesSemiconductor.NAminus);
            Console.WriteLine("NA neu: " + ionizedAcceptorDensity);
            */

            if (modelSemiconductor.useIncompleteIonization)
                return ionizedDonorDensity - ionizedAcceptorDensity;
            else
                return material.propertiesSemiconductor.NDplus - material.propertiesSemiconductor.NAminus;
        }

        /// <summary>
        /// returns the derivation of the doping function with respect to Phi
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double DopingDerivationPhi(ModelSemiconductor modelSemiconductor)
        {
            double referenceDesnityDonor = (material.propertiesSemiconductor.Nc(modelSemiconductor.T) * Math.Exp(donorEnergy / modelSemiconductor.U_T));
            double referenceDesnityAcceptor = (material.propertiesSemiconductor.Nv(modelSemiconductor.T) * Math.Exp(-acceptorEnergy / modelSemiconductor.U_T));

            double ionizedDonorDensity = material.propertiesSemiconductor.NDplus / (1 + degeneracyFactorDonor * nDensity(modelSemiconductor) / referenceDesnityDonor);
            double ionizedAcceptorDensity = material.propertiesSemiconductor.NAminus / (1 + degeneracyFactorAcceptor * pDensity(modelSemiconductor) / referenceDesnityAcceptor);

            double derivationDonor = -ionizedDonorDensity * degeneracyFactorDonor * nDerivationPhi(modelSemiconductor) / referenceDesnityDonor;
            double derivationAcceptor = -ionizedAcceptorDensity * degeneracyFactorAcceptor * pDerivationPhi(modelSemiconductor) / referenceDesnityAcceptor;

            if (modelSemiconductor.useIncompleteIonization)
                return derivationDonor - derivationAcceptor;
            else
                return 0;

        }

        /// <summary>
        /// returns the derivation of the doping function with respect to Phi_n
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double DopingDerivationPsin(ModelSemiconductor modelSemiconductor)
        {
            double referenceDesnityDonor = (material.propertiesSemiconductor.Nc(modelSemiconductor.T) * Math.Exp(donorEnergy / modelSemiconductor.U_T));
            double referenceDesnityAcceptor = (material.propertiesSemiconductor.Nv(modelSemiconductor.T) * Math.Exp(-acceptorEnergy / modelSemiconductor.U_T));

            double ionizedDonorDensity = material.propertiesSemiconductor.NDplus / (1 + degeneracyFactorDonor * nDensity(modelSemiconductor) / referenceDesnityDonor);
            double ionizedAcceptorDensity = material.propertiesSemiconductor.NAminus / (1 + degeneracyFactorAcceptor * pDensity(modelSemiconductor) / referenceDesnityAcceptor);

            double derivationDonor = -ionizedDonorDensity * degeneracyFactorDonor * nDerivationPhi_n(modelSemiconductor) / referenceDesnityDonor;
            double derivationAcceptor = -ionizedAcceptorDensity * degeneracyFactorAcceptor * pDerivationPhi_n(modelSemiconductor) / referenceDesnityAcceptor; // =0

            if (modelSemiconductor.useIncompleteIonization)
                return derivationDonor - derivationAcceptor;
            else
                return 0;

        }

        /// <summary>
        /// returns the derivation of the doping function with respect to Phi_p
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double DopingDerivationPsip(ModelSemiconductor modelSemiconductor)
        {
            double referenceDesnityDonor = (material.propertiesSemiconductor.Nc(modelSemiconductor.T) * Math.Exp(donorEnergy / modelSemiconductor.U_T));
            double referenceDesnityAcceptor = (material.propertiesSemiconductor.Nv(modelSemiconductor.T) * Math.Exp(-acceptorEnergy / modelSemiconductor.U_T));

            double ionizedDonorDensity = material.propertiesSemiconductor.NDplus / (1 + degeneracyFactorDonor * nDensity(modelSemiconductor) / referenceDesnityDonor);
            double ionizedAcceptorDensity = material.propertiesSemiconductor.NAminus / (1 + degeneracyFactorAcceptor * pDensity(modelSemiconductor) / referenceDesnityAcceptor);

            double derivationDonor = -ionizedDonorDensity * degeneracyFactorDonor * nDerivationPhi_p(modelSemiconductor) / referenceDesnityDonor; //=0
            double derivationAcceptor = -ionizedAcceptorDensity * degeneracyFactorAcceptor * pDerivationPhi_p(modelSemiconductor) / referenceDesnityAcceptor;

            if (modelSemiconductor.useIncompleteIonization)
                return derivationDonor - derivationAcceptor;
            else
                return 0;

        }


        /// <summary>
        /// Returns the intrinsic carrier density of the material
        /// </summary>
        /// <param name="temperature"></param>
        /// <returns></returns>
        public double Nintr(double temperature)
        {
            return Math.Sqrt(material.propertiesSemiconductor.Nc(temperature) * material.propertiesSemiconductor.Nv(temperature) * Math.Exp(-localBandGap * physConstants.e / (physConstants.kB * temperature)));
        }

        /// <summary>
        /// calculates intrinsic fermi level of the Material 
        /// </summary>
        /// <param name="temperature"></param>
        /// <returns></returns>
        public double localIntrinsicLevel(double temperature = 300)
        {
            return material.propertiesSemiconductor.chemicalPotential - localBandGap / 2 - (physConstants.kB * 300) / (2 * physConstants.e) * Math.Log(material.propertiesSemiconductor.Nc(temperature) / material.propertiesSemiconductor.Nv(temperature));

        }




        public double FermiDiracDensity(ModelSemiconductor modelSemiconductor, double etaArgument)
        {
            double splineVlaue = modelSemiconductor.FermiDiracSpline.ValueAt(etaArgument);
            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * Math.Exp(splineVlaue);

        }
        public double FermiDiracDensityDerivation(ModelSemiconductor modelSemiconductor, double etaArgument)
        {
            double splineDerivationValue = modelSemiconductor.FermiDiracDerivation.ValueAt(etaArgument);
            double splineVlaue = modelSemiconductor.FermiDiracSpline.ValueAt(etaArgument);

            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * splineDerivationValue * Math.Exp(splineVlaue);
        }

        public double distributionFuction(ModelSemiconductor modelSemiconductor, double eta)
        {

            //return Math.Exp(eta);

            if (modelSemiconductor.useFermiDirac) //TODO Tim: Zusätzlich abhängig von Betrag von Eta machen
                return Math.Exp(modelSemiconductor.FermiDiracSpline.ValueAt(eta));
            else if (modelSemiconductor.useBlakemoreStatistic)
                return Math.Pow(Math.Exp(-eta) + 0.27, -1);
            else
                return Math.Exp(eta);
        }
        public double distributionFuctionDerivation(ModelSemiconductor modelSemiconductor, double eta)
        {

            //return Math.Exp(eta) ;

            if (modelSemiconductor.useFermiDirac)
                return Math.Exp(modelSemiconductor.FermiDiracDerivation.ValueAt(eta));
            else if (modelSemiconductor.useBlakemoreStatistic)
                return Math.Pow(Math.Exp(-eta) + 0.27, -2) * Math.Exp(-eta);
            else
                return Math.Exp(eta);
        }


        /// <summary>
        /// gives wether an eta of electrons or of holes is calculated
        /// </summary>
        public enum EtaChargeType { eta_Electrons, eta_Holes }

        public double etaDerivation(DerivationVariable derivationVariable, EtaChargeType chargeType, ModelSemiconductor modelSemiconductor)
        {
            if (chargeType == EtaChargeType.eta_Electrons)
            {
                switch (derivationVariable)
                {
                    case DerivationVariable.Phi:
                        return 1 / (modelSemiconductor.U_T);
                    case DerivationVariable.Phi_n:
                        return 1 / (modelSemiconductor.U_T);
                    case DerivationVariable.Phi_p:
                        return 0;
                    default:
                        break;
                }

                throw new Exception("Error in Function etaDerivation. Check Density and Current equations.");

            }
            else if (chargeType == EtaChargeType.eta_Holes)
            {
                switch (derivationVariable)
                {
                    case DerivationVariable.Phi:
                        return -1 / (modelSemiconductor.U_T);
                    case DerivationVariable.Phi_n:
                        return 0;
                    case DerivationVariable.Phi_p:
                        return -1 / (modelSemiconductor.U_T);
                    default:
                        break;
                }

                throw new Exception("Error in Function etaDerivation. Check Density and Current equations.");
            }
            else
                throw new Exception("Error in Function etaDerivation. Check Density and Current equations.");
        }

        public double eta(EtaChargeType chargeType, ModelSemiconductor modelSemiconductor)
        {
            if (chargeType == EtaChargeType.eta_Electrons)
                return (((phi + phi_n) - material.propertiesSemiconductor.chemicalPotential) / (modelSemiconductor.U_T));
            else if (chargeType == EtaChargeType.eta_Holes)
                return ((-phi_p - phi) + material.propertiesSemiconductor.chemicalPotential - localBandGap) / (modelSemiconductor.U_T);
            else
                throw new Exception("Error in eta function. Check Density and Current equations.");

        }



        /// <summary>
        /// returns the electron density at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double nDensity(ModelSemiconductor modelSemiconductor)
        {
            double etaN = eta(EtaChargeType.eta_Electrons, modelSemiconductor);
            double n = material.propertiesSemiconductor.Nc(modelSemiconductor.T) * distributionFuction(modelSemiconductor, etaN);

            //if(etaN < -80 ||etaN > 20)
            //Console.WriteLine("------------------------------------------------> EtaN out of Integration Range: " +  etaN +  "<----------------------------------------------------------------------------------");
            //Console.WriteLine("\nFermiDirac: " + n);
            //Console.WriteLine("Boltzmann: " + material.propertiesSemiconductor.Nc * Math.Exp(etaN));

            return n;
        }
        /// <summary>
        /// returns the hole density at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double pDensity(ModelSemiconductor modelSemiconductor)
        {
            double etaP = eta(EtaChargeType.eta_Holes, modelSemiconductor);
            return material.propertiesSemiconductor.Nv(modelSemiconductor.T) * distributionFuction(modelSemiconductor, etaP);

        }
        /// <summary>
        /// returns the multiplied electron and hole density at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double np_density(ModelSemiconductor modelSemiconductor)
        {
            return nDensity(modelSemiconductor) * pDensity(modelSemiconductor);
        }
        /// <summary>
        /// returns the derivation of np density with respect to Phi at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double npDerivationPhi(ModelSemiconductor modelSemiconductor)
        {
            //return 0;
            return nDensity(modelSemiconductor) * pDerivationPhi(modelSemiconductor) + nDerivationPhi(modelSemiconductor) * pDensity(modelSemiconductor);
        }
        /// <summary>
        /// returns the derivation of np density with respect to Phi_n at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double npDerivationPhi_n(ModelSemiconductor modelSemiconductor)
        {
            //return 1 / modelSemiconductor.U_T * np_density(modelSemiconductor);
            return nDensity(modelSemiconductor) * pDerivationPhi_n(modelSemiconductor) + nDerivationPhi_n(modelSemiconductor) * pDensity(modelSemiconductor);


        }
        /// <summary>
        /// returns the derivation of np density with respect to Phi_p at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double npDerivationPhi_p(ModelSemiconductor modelSemiconductor)
        {
            //return -1 / modelSemiconductor.U_T * np_density(modelSemiconductor);
            return nDensity(modelSemiconductor) * pDerivationPhi_p(modelSemiconductor) + nDerivationPhi_p(modelSemiconductor) * pDensity(modelSemiconductor);

        }
        /// <summary>
        /// returns the derivation of the electron density with respect to Phi at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double nDerivationPhi(ModelSemiconductor modelSemiconductor)
        {
            double etaN = eta(EtaChargeType.eta_Electrons, modelSemiconductor);
            double etaDerivationPhi = etaDerivation(DerivationVariable.Phi, EtaChargeType.eta_Electrons, modelSemiconductor);
            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * etaDerivationPhi * distributionFuctionDerivation(modelSemiconductor, etaN);
            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * etaDerivationPhi * Math.Exp(etaN);// distributionFuctionDerivation(modelSemiconductor, etaN) ;
        }
        /// <summary>
        /// returns the derivation of the electron density with respect to Phi_n at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double nDerivationPhi_n(ModelSemiconductor modelSemiconductor)
        {
            double etaN = eta(EtaChargeType.eta_Electrons, modelSemiconductor);
            double etaDerivationPhiN = etaDerivation(DerivationVariable.Phi_n, EtaChargeType.eta_Electrons, modelSemiconductor);
            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * etaDerivationPhiN * distributionFuctionDerivation(modelSemiconductor, etaN);
            return material.propertiesSemiconductor.Nc(modelSemiconductor.T) * etaDerivationPhiN * Math.Exp(etaN);//distributionFuctionDerivation(modelSemiconductor, etaN);
        }
        /// <summary>
        /// returns the derivation of the electron density with respect to Phi_p at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double nDerivationPhi_p(ModelSemiconductor modelSemiconductor)
        {
            return 0;
        }
        /// <summary>
        /// returns the derivation of the hole density with respect to Phi at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double pDerivationPhi(ModelSemiconductor modelSemiconductor)
        {
            double etaP = eta(EtaChargeType.eta_Holes, modelSemiconductor);
            double etaDerivationPhi = etaDerivation(DerivationVariable.Phi, EtaChargeType.eta_Holes, modelSemiconductor);

            return material.propertiesSemiconductor.Nv(modelSemiconductor.T) * etaDerivationPhi * distributionFuctionDerivation(modelSemiconductor, etaP);
            return material.propertiesSemiconductor.Nv(modelSemiconductor.T) * etaDerivationPhi * Math.Exp(etaP);//distributionFuctionDerivation(modelSemiconductor, etaP) ;

        }
        /// <summary>
        /// returns the derivation of the hole density with respect to Phi_n at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double pDerivationPhi_n(ModelSemiconductor modelSemiconductor)
        {
            return 0;
        }
        /// <summary>
        /// returns the derivation of the hole density with respect to Phi_p at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double pDerivationPhi_p(ModelSemiconductor modelSemiconductor)
        {
            double etaP = eta(EtaChargeType.eta_Holes, modelSemiconductor);
            double etaDerivationPhip = etaDerivation(DerivationVariable.Phi_p, EtaChargeType.eta_Holes, modelSemiconductor);

            return material.propertiesSemiconductor.Nv(modelSemiconductor.T) * etaDerivationPhip * distributionFuctionDerivation(modelSemiconductor, etaP);
            return material.propertiesSemiconductor.Nv(modelSemiconductor.T) * etaDerivationPhip * Math.Exp(etaP);//distributionFuctionDerivation(modelSemiconductor, etaP);
        }








        // Surface Recombination Currents---------------------------------------------------------------------------------------------
        /// <summary>
        /// returns the surface recombination velocity of electrons depending on the contact
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double SurfaceRecVelocityElectrons(ModelSemiconductor modelSemiconductor)
        {
            double electronSurfaceVelocity;
            // p contact
            if (hasOperatingVoltage)
                electronSurfaceVelocity = 1e2;
            // n contact
            else
                electronSurfaceVelocity = 1e5;

            return electronSurfaceVelocity;
        }
        /// <summary>
        /// returns the surface recombination velocity of holes depending on the contact
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double SurfaceRecVelocityHoles(ModelSemiconductor modelSemiconductor)
        {
            double holeSurfaceVelocity;
            // p contact
            if (hasOperatingVoltage)
                holeSurfaceVelocity = 1e5;
            // n contact
            else
                holeSurfaceVelocity = 1e2;

            return holeSurfaceVelocity;
        }


        /// <summary>
        /// returns the current of electrons surface recombination
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double ElectronSurfaceRecombinationCurrent(ModelSemiconductor modelSemiconductor)
        {
            //return - physConstants.e * SurfaceRecVelocityElectrons(modelSemiconductor) *  (nDensity(modelSemiconductor) - nEquilibrium);
            return -physConstants.e * contactPreferences.SRV_electrons * (nDensity(modelSemiconductor) - nEquilibrium);
        }
        /// <summary>
        /// returns the current of hole surface recombination
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double HoleSurfaceRecombinationCurrent(ModelSemiconductor modelSemiconductor)
        {

            //return physConstants.e * SurfaceRecVelocityHoles(modelSemiconductor) * (pDensity(modelSemiconductor) - pEquilibrium);
            return physConstants.e * contactPreferences.SRV_holes * (pDensity(modelSemiconductor) - pEquilibrium);
        }


        /// <summary>
        /// Returns the derivation of surface recombination currents for ELECTRONS, sign is integrated in function, use with (+)
        /// </summary>
        /// <param name="pointIndex">Index of the point</param>
        /// <param name="leftOrRigth">specification of wether left or right contact is calculated, 0 for left, 1 for right</param>
        /// <param name="derivWithRespectTo">specification of the potential wich is used for the derivation, 0 for Phi, 1 for Phi_n, 2 for Phi_p</param>
        public double surfaceRecElectronsDerivation(ModelSemiconductor modelSemiconductor, int pointIndex, int derivWithRespectTo)
        {

            //double electronSurfaceVelocity = SurfaceRecVelocityElectrons(modelSemiconductor);
            double electronSurfaceVelocity = contactPreferences.SRV_electrons;

            if (derivWithRespectTo == 0) //Phi
                return -physConstants.e * electronSurfaceVelocity * nDerivationPhi(modelSemiconductor);
            else if (derivWithRespectTo == 1) //Phi_n
                return -physConstants.e * electronSurfaceVelocity * nDerivationPhi_n(modelSemiconductor);
            else if (derivWithRespectTo == 2) //Phi_p
                return -physConstants.e * electronSurfaceVelocity * nDerivationPhi_p(modelSemiconductor);

            else
                throw new Exception("Surface Recombination Electrons:  error");
        }
        /// <summary>
        /// Returns the derivation of surface recombination currents for HOLES, sign is integrated in function, use with (+)
        /// </summary>
        /// <param name="pointIndex">Index of the point</param>
        /// <param name="leftOrRigth">specification of wether left or right contact is calculated, 0 for left, 1 for right</param>
        /// <param name="derivWithRespectTo">specification of the potential wich is used for the derivation, 0 for Phi, 1 for Phi_n, 2 for Phi_p</param>
        /// <returns></returns>
        public double surfaceRecHolesDerivation(ModelSemiconductor modelSemiconductor, int pointIndex, int derivWithRespectTo)
        {
            //double holeSurfaceVelocity = SurfaceRecVelocityHoles(modelSemiconductor);
            double holeSurfaceVelocity = contactPreferences.SRV_holes;

            if (derivWithRespectTo == 0) //Phi
                return +physConstants.e * holeSurfaceVelocity * pDerivationPhi(modelSemiconductor);
            else if (derivWithRespectTo == 1) //Phi_n
                return +physConstants.e * holeSurfaceVelocity * pDerivationPhi_n(modelSemiconductor);
            else if (derivWithRespectTo == 2) //Phi_p
                return +physConstants.e * holeSurfaceVelocity * pDerivationPhi_p(modelSemiconductor);

            else
                throw new Exception("Surface Recombination Holes:  error");
        }




        //Recombination and Generation Equations-----------------------------------------------------------------------

        /// <summary>
        /// Returns the total optical generation rate at a given mesh point
        /// </summary>
        /// <param name="modelSemiconductor"></param>
        /// <returns></returns>
        public double TotalGenerationRate(ModelSemiconductor modelSemiconductor) // R(n,p)
        {
            if (modelSemiconductor.enableGeneration)
            {
                // TMM optic calculation
                //Console.WriteLine(modelSemiconductor.generationArray[index].generationRate / modelSemiconductor.globalGenerationRampingFactor);
                return modelSemiconductor.generationArray[index].generationRate / modelSemiconductor.globalGenerationRampingFactor;

        }
            else
                return 0;
        }

        /// <summary>
        /// returns the total recombination rate at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double TotalRecombinationRate(ModelSemiconductor modelSemiconductor) // R(n,p)
        {
            double x = (TotalRecombinationCoefficient(modelSemiconductor)) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2));
            //Console.WriteLine("TotalRcoeff:........" + TotalRecombinationCoefficient(modelSemiconductor));
            //Console.WriteLine("np" + np_density(modelSemiconductor));
            //Console.WriteLine("Nintr" + Math.Pow(material.Nintr(modelSemiconductor.T), 2));
            //Console.WriteLine();
            /*
            if (hasInterfaceCondition)
            {
                Console.WriteLine("INterface Rek: " + InterfaceRecombinationRate(modelSemiconductor));
                Console.WriteLine("SRH Rek Rate: " + SRHRecombinationRate(modelSemiconductor));
                Console.WriteLine("Rad Rek Rate: " + SpontaneousRecombinationRate(modelSemiconductor));
                Console.WriteLine();
            }
            */



                return x;
        }
        /// <summary>
        /// returns the SRH recombination rate at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double SRHRecombinationRate(ModelSemiconductor modelSemiconductor) // R_SRH(n,p)
        {
            double d = (RecombinationCoefficientSRH(modelSemiconductor)) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2));
            return d;
        }
        /// <summary>
        /// returns the Auger recombination rate at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double AugerRecombinationRate(ModelSemiconductor modelSemiconductor) // R_Auger(n,p)
        {
            return (RecombinationCoefficientAuger(modelSemiconductor)) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2));
        }
        /// <summary>
        /// returns the radiative recombination rate at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double SpontaneousRecombinationRate(ModelSemiconductor modelSemiconductor) // R_Spontaneous(n,p)
        {
            if (modelSemiconductor.useRadiativeRecombination == true)
                return (material.propertiesSemiconductor.rSpontaneous) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2));
            else return 0;
        }

        // Recombination Coefficients r:
        /// <summary>
        /// returns the sum of recombination coefficients for SRH, Auger and radiative recombination
        /// </summary>
        /// <returns></returns>
        public double TotalRecombinationCoefficient(ModelSemiconductor modelSemiconductor) //r_total
        {
            if (modelSemiconductor.useRadiativeRecombination == false)
                return RecombinationCoefficientSRH(modelSemiconductor)
                    + RecombinationCoefficientAuger(modelSemiconductor)
                    + InterfaceRecombinationcoefficient(modelSemiconductor);

            return material.propertiesSemiconductor.rSpontaneous
                + RecombinationCoefficientSRH(modelSemiconductor)
                + RecombinationCoefficientAuger(modelSemiconductor)
                + InterfaceRecombinationcoefficient(modelSemiconductor);
        }
        /// <summary>
        /// returns the coefficient of SRH recombination
        /// </summary>
        /// <returns></returns>
        public double RecombinationCoefficientSRH(ModelSemiconductor modelSemiconductor) //r_SRH
        {
            if (modelSemiconductor.useSrhRecombination == true)
            {
                double recRate = 0;
                foreach (var p in material.propertiesSemiconductor.defectList)
                {
                    // R_i = 1 / (tau_p*(n + n_T) + tau_n*(p + p_T))
                    // with tau_n,p = (v_n,p * sigma_n,p * N_T)^-1
                    // and p_T or n_T = Nintr * exp(E_T / U_T)

                    double a = material.propertiesSemiconductor.holeThermalVelocity;
                    double b = p.tau_p(a);
                    double c = p.nTrapRef(Nintr(modelSemiconductor.T), modelSemiconductor.T, localBandGap, localIntrinsicLevel());
                    double e = material.propertiesSemiconductor.electronThermalVelocity;
                    double f = p.tau_n(e);
                    double g = p.pTrapRef(Nintr(modelSemiconductor.T), modelSemiconductor.T, localBandGap, localIntrinsicLevel());
                    double h = b * (nDensity(modelSemiconductor) + c) + f * (pDensity(modelSemiconductor) + g);
                    // double d = Math.Pow(h, -1);
                    double d = 1 / h;


                    //double d = Math.Pow(p.tau_p(material.propertiesSemiconductor.holeThermalVelocity) * (nDensity(modelSemiconductor) + p.nTrap(material.propertiesSemiconductor.Nintr, modelSemiconductor.T))
                    //  + p.tau_n(material.propertiesSemiconductor.electronThermalVelocity) * (pDensity(modelSemiconductor) + p.pTrap(material.propertiesSemiconductor.Nintr, modelSemiconductor.T)), -1);
                    /*
                    Console.WriteLine("Coeff 1" + p.tau_p(material.holeThermalVelocity));
     g               Console.WriteLine("Coeff 2" + p.nTrap(material.Nintr, modelSemiconductor.T));
                    Console.WriteLine("Coeff 3" + p.tau_n(material.electronThermalVelocity));
                    Console.WriteLine("Coeff 4" + p.pTrap(material.Nintr, modelSemiconductor.T));
                    */
                    recRate += d;
                }
                return recRate;
            }
            else return 0;
        }

        public double scaleCoeff()
        {
            return 1e0;
            return borderEdgeSize / size;
        }

        public double scaleCoeffSRV()
        {
            return borderEdgeSize / size;
            return 1e0;
        }

        /// <summary>
        /// returns the Interface Recombination rate at a given meshpoint
        /// </summary>
        /// <returns></returns>
        public double InterfaceRecombinationRate(ModelSemiconductor modelSemiconductor) // R_SRH(n,p)
        {
            double d = (InterfaceRecombinationcoefficient(modelSemiconductor)) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2));
            
            return d;
        }

        /// <summary>
        /// returns the SRH recombination coefficient at semiconductor interfaces
        /// </summary>
        /// <returns></returns>
        public double InterfaceRecombinationcoefficient(ModelSemiconductor modelSemiconductor)
        {
            if (hasInterfaceCondition)
            {

                double Sp = contactPreferences.SRV_holes* scaleCoeffSRV();
                double Sn = contactPreferences.SRV_electrons * scaleCoeffSRV();

                double energeticPosition = contactPreferences.interfaceTrapEnergy * physConstants.e;
                double nTrap = Nintr(modelSemiconductor.T) * Math.Exp(energeticPosition / (physConstants.kB * modelSemiconductor.T));
                double pTrap = Nintr(modelSemiconductor.T) * Math.Exp(-energeticPosition / (physConstants.kB * modelSemiconductor.T));
                double recCoeff = 1 / ((nDensity(modelSemiconductor) + nTrap) / Sp + (pDensity(modelSemiconductor) + pTrap) / Sn);
                recCoeff *= scaleCoeff();
                /*
                Console.WriteLine("Sn: " + Sn + "at Point: " + index);
                Console.WriteLine("Sp: " + Sp + "at Point: " + index);
                Console.WriteLine("IF Rec an Punkt " + index + " : " + recCoeff * borderEdgeSize / size);
                Console.WriteLine("BorderEdgeSize " + borderEdgeSize);
                Console.WriteLine("size " + size);
                Console.WriteLine(" energeticPosition " + energeticPosition);
                Console.WriteLine(" nTrap " + nTrap );
                Console.WriteLine(" pTrap " + pTrap );
                Console.WriteLine(" Nintr " + Nintr(modelSemiconductor.T));
                Console.WriteLine(" exp " + Math.Exp(energeticPosition / (physConstants.kB * modelSemiconductor.T)));
                */
                return recCoeff;
            }
            else
                return 0;
        }
        public double InterfaceRecCoeff_Diff_Phi(ModelSemiconductor modelSemiconductor)
        {
            double Sp = contactPreferences.SRV_holes * scaleCoeffSRV();
            double Sn = contactPreferences.SRV_electrons * scaleCoeffSRV();

            if (hasInterfaceCondition)
            {
                double recCoeff = -(1 / Sp * nDerivationPhi(modelSemiconductor)
                                    + 1 / Sn * pDerivationPhi(modelSemiconductor))
                                    * Math.Pow(InterfaceRecombinationcoefficient(modelSemiconductor), 2);
                return recCoeff * scaleCoeff();
            }
            else
                return 0;
        }
        public double InterfaceRecCoeff_Diff_Phin(ModelSemiconductor modelSemiconductor)
        {
            double Sp = contactPreferences.SRV_holes * scaleCoeffSRV();

            if (hasInterfaceCondition)
            {
                double recCoeff = -1 / Sp * nDerivationPhi_n(modelSemiconductor) * Math.Pow(InterfaceRecombinationcoefficient(modelSemiconductor), 2);
                return recCoeff * scaleCoeff();
            }
            else
                return 0;
        }
        public double InterfaceRecCoeff_Diff_Phip(ModelSemiconductor modelSemiconductor)
        {
            double Sn = contactPreferences.SRV_electrons * scaleCoeffSRV();

            if (hasInterfaceCondition)
            {
                double recCoeff = -1 / Sn * pDerivationPhi_p(modelSemiconductor) * Math.Pow(InterfaceRecombinationcoefficient(modelSemiconductor), 2);
                return recCoeff * scaleCoeff();
            }
            else
                return 0;
        }


        /// <summary>
        /// returns the coefficent of auger recombination
        /// </summary>
        /// <returns></returns>
        public double RecombinationCoefficientAuger(ModelSemiconductor modelSemiconductor) //r_Auger
        {
            if (modelSemiconductor.useAugerRecombination == true)
            {
                return material.propertiesSemiconductor.Cn * nDensity(modelSemiconductor) + material.propertiesSemiconductor.Cp * pDensity(modelSemiconductor);
            }
            else return 0;
        }

        //Derivations of recombination coefficients:
        /// <summary>
        /// returns the derivation of the total recombination coefficient with respect to Phi_n
        /// </summary>
        /// <returns></returns>
        public double DiffRecombCoeffPsin(ModelSemiconductor modelSemiconductor) //dr/dPsin
        {
            double SRHCoeffDiffPsin = 0;
            double AugerCoeffDiffPsin = 0;
            if (modelSemiconductor.useAugerRecombination == true)
                AugerCoeffDiffPsin = material.propertiesSemiconductor.Cn * nDerivationPhi_n(modelSemiconductor);

            if (modelSemiconductor.useSrhRecombination == true)
            {
                foreach (var p in material.propertiesSemiconductor.defectList)
                {
                    SRHCoeffDiffPsin += -p.tau_p(material.propertiesSemiconductor.holeThermalVelocity) * nDerivationPhi_n(modelSemiconductor) * Math.Pow(RecombinationCoefficientSRH(modelSemiconductor), 2);
                }
            }

            return AugerCoeffDiffPsin + SRHCoeffDiffPsin + InterfaceRecCoeff_Diff_Phin(modelSemiconductor);
            //return   (tau_p * Math.Pow(RecombinationCoefficientSRH( ),2) - Cn) * n_density( psin) / main.U_T;
        }
        /// <summary>
        /// returns the derivation of the total recombination coefficient with respect to Phi_p
        /// </summary>
        /// <returns></returns>
        public double DiffRecombCoeffPsip(ModelSemiconductor modelSemiconductor) //dr/dPsip
        {
            double SRHCoeffDiffPsip = 0;
            double AugerCoeffDiffPsip = 0;
            if (modelSemiconductor.useAugerRecombination == true)
                AugerCoeffDiffPsip = material.propertiesSemiconductor.Cp * pDerivationPhi_p(modelSemiconductor);

            if (modelSemiconductor.useSrhRecombination == true)
            {
                foreach (var p in material.propertiesSemiconductor.defectList)
                {
                    SRHCoeffDiffPsip += -p.tau_n(material.propertiesSemiconductor.electronThermalVelocity) * Math.Pow(RecombinationCoefficientSRH(modelSemiconductor), 2) * pDerivationPhi_p(modelSemiconductor);

                }
            }

            return AugerCoeffDiffPsip + SRHCoeffDiffPsip + InterfaceRecCoeff_Diff_Phip(modelSemiconductor);
            //return  ( Cp - tau_n * Math.Pow(RecombinationCoefficientSRH( ),2) ) * p_density( psip) / main.U_T;
        }
        /// <summary>
        /// returns the derivation of the total recombination coefficient with respect to Phi
        /// </summary>
        /// <returns></returns>
        public double DiffRecombCoeffPhi(ModelSemiconductor modelSemiconductor) //dr/dPhi
        {
            double SRHCoeffDiffPhi = 0;
            double AugerCoeffDiffPhi = 0;
            if (modelSemiconductor.useAugerRecombination == true)
                AugerCoeffDiffPhi = material.propertiesSemiconductor.Cn * nDerivationPhi(modelSemiconductor) + material.propertiesSemiconductor.Cp * pDerivationPhi(modelSemiconductor);

            if (modelSemiconductor.useSrhRecombination == true)
            {
                foreach (var p in material.propertiesSemiconductor.defectList)
                {
                    SRHCoeffDiffPhi += -(p.tau_p(material.propertiesSemiconductor.holeThermalVelocity) * nDerivationPhi(modelSemiconductor)
                        + p.tau_n(material.propertiesSemiconductor.electronThermalVelocity) * pDerivationPhi(modelSemiconductor)) * Math.Pow(RecombinationCoefficientSRH(modelSemiconductor), 2);// / Math.Pow(tau_p * (n_density( psin) + nTrap)

                }
            }
            return AugerCoeffDiffPhi + SRHCoeffDiffPhi + InterfaceRecCoeff_Diff_Phi(modelSemiconductor);
        }

        // Derivations of recombinations rates: 
        /// <summary>
        /// returns the derivation of the total recombination rate with respect to Phi_n
        /// </summary>
        /// <returns></returns>
        public double Diff_Rekomb_Psi_n(ModelSemiconductor modelSemiconductor) //R nach Psi_n abgeleitet
        {
            return DiffRecombCoeffPsin(modelSemiconductor) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2))
                + (TotalRecombinationCoefficient(modelSemiconductor) * npDerivationPhi_n(modelSemiconductor));
        }
        /// <summary>
        /// returns the derivation of the total recombination rate with respect to Phi_p
        /// </summary>
        /// <returns></returns>
        public double Diff_Rekomb_Psi_p(ModelSemiconductor modelSemiconductor) //R nach Psi_p abgeleitet
        {
            return DiffRecombCoeffPsip(modelSemiconductor) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2))
                + (TotalRecombinationCoefficient(modelSemiconductor) * npDerivationPhi_p(modelSemiconductor));
        }
        /// <summary>
        /// returns the derivation of the total recombination rate with respect to Phi
        /// </summary>
        /// <returns></returns>
        public double Diff_Rekomb_Phi(ModelSemiconductor modelSemiconductor) //R nach Phi abgeleitet
        {
            return DiffRecombCoeffPhi(modelSemiconductor) * (np_density(modelSemiconductor) - Math.Pow(Nintr(modelSemiconductor.T), 2))
                + (TotalRecombinationCoefficient(modelSemiconductor) * npDerivationPhi(modelSemiconductor));
        }



    }
}