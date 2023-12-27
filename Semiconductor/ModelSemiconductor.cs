using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using Extreme.Mathematics;
using Extreme.Mathematics.LinearAlgebra;
using Geometry;
using BasicLib;
using System.Threading.Tasks;
using TransferMatrix;
using Database;
using Extreme.Mathematics.Calculus;
using Extreme.Mathematics.Curves;
using MoreLinq;
using System.Windows;

namespace Semiconductor
{
    public class ModelSemiconductor
    {
        /// <summary>
        /// Name of the semiconductor model
        /// </summary>
        public string name { get; private set; }


        //todo: temperatur über GUI (dazu in Modelsemiconductor Konstruktoren ersetzen)
        public double T { get; set; }//{ get { return T; }  set{ U_T = physConstants.kB * value / physConstants.e; } }


        public double U_T { get; set; }

        /// <summary>
        /// gives the optical model of the layerstack for the TMM calculation
        /// </summary>
        public ModelTMM modelOpticsTMM { get; set; }
        /// <summary>
        /// gives the unchanged input spectrum of the model. This means th underlying spectrum with the in the GUI stated borders, 
        /// BUT WITHOUT changes due to incoherent layers or range changes in EQE or loss analysis calculations
        /// </summary>
        public Spectrum referenceSpectrum { get; set; }
        /// <summary>
        /// gives the material before the stack
        /// </summary>
        public int materialBefore { get; private set; }
        /// <summary>
        /// gives the material behind the stack
        /// </summary>
        public (int ID, double roughnessOnTop) materialBehind { get; private set; }


        /// <summary>
        /// gives the optical model /geometry of the modell: all incoherent layers (before and after) and all coherent layers (before, after and the actual semiconductor stack)
        /// </summary>
        public (List<(Material material, double thickness, double roughnessOnTop, bool isAbsorber)> SemiconductorStack, List<(Material material, double thickness)> incoherentBeforeStack,
            List<(Material material, double thickness)> coherentBeforeStack, List<(Material material, double thickness)> incoherentBehindStack,
            List<(Material material, double thickness)> coherentBehindStack) opticalGeometry { get; set; }

        /// <summary>
        /// List with absorbed parts of the spectrum per incoherent layer before the semiconductor layerstack
        /// </summary>
        public List<double> absorptanceIncohBefore { get; set; }
        /// <summary>
        /// List with absorbed power of the spectrum per incoherent layer before the semiconductor layerstack
        /// </summary>
        public List<double> absorbedPowerIncohBefore;

        /// <summary>
        /// gives the ID of the absorber (active layer) of the semiconductor modell (e.g. CIGS, Pero, ...) From this loss analysis, spectrum data or SQ calculations are dervied
        /// </summary>
        public int absorberID { get; set; }
        /// <summary>
        /// gives the thickness of the absorber (active layer) of the semiconductor modell (e.g. CIGS, Pero, ...) From this loss analysis, spectrum data or SQ calculations are dervied
        /// </summary>
        public double absorberThickness { get; set; }
        /// <summary>
        /// gives the band gap of the absorber (active layer) of the semiconductor modell (e.g. CIGS, Pero, ...) From this loss analysis, spectrum data or SQ calculations are dervied
        /// </summary>
        public double absorberBandGap { get; set; }

        /// <summary>
        /// total power of the used spectrum (before incoherent layer etc.)
        /// </summary>
        double totalInputPower;

        /// <summary>
        /// defiones the voltage steps, which are used for the ramping of the operating voltage when simulating an IV curve
        /// </summary>
        public double deltaU = 0.01;

        /// <summary>
        /// bool which defines whether the calculation of densities via the blakemore statistic is aktivated or not (for the whole device)
        /// </summary>
        public bool useBlakemoreStatistic { get; set; } = false;
        /// <summary>
        /// bool which defines whether the calculation of densities via the Fermi-Dirac statistic is aktivated or not (for the whole device)
        /// </summary>
        public bool useFermiDirac { get; set; } = false;

        /// <summary>
        /// bool which defines whether radiative recombination is aktivated or not (for the whole device)
        /// </summary>
        public bool useRadiativeRecombination { get; set; }
        /// <summary>
        /// bool which defines whether auger recombination is aktivated or not (for the whole device)
        /// </summary>
        public bool useAugerRecombination { get; set; }
        /// <summary>
        /// bool which defines whether radiative recombination is aktivated or not (for the whole device)
        /// </summary>
        public bool useSrhRecombination { get; set; }


        /// <summary>
        /// effecitve radiative coefficient
        /// </summary>
        public double radiativeCoefficientSQ { get; set; }

        public bool useFarrellBoundaryConditions { get; set; } = false;
        /// <summary>
        /// defines whether the simulation peformed under dark or light conditions
        /// </summary>
        public bool enableGeneration { get; set; }
        /// <summary>
        /// Value which is used as start for the ramping of the generation rate
        /// </summary>
        public double globalGenerationRampingFactor { get; set; }
        /// <summary>
        /// reset the generation ramping factor 
        /// </summary>
        public bool resetRampingFactor { get; set; } = false;
        /// <summary>
        /// bool whether an IV calculation is stopped after the direction of the current changed
        /// </summary>
        public bool stopAfterVoc { get; set; } = true;
        /// <summary>
        /// if stopAfterVOc is true this factor defines the stop condition (stop current = |Iph| *factor)
        /// </summary>
        public double stopAfterVocFactorIph { get; set; } = 1.2;
        /// <summary>
        /// whether temperature dependent equations for doping etc. are used
        /// </summary>
        public bool useIncompleteIonization { get; set; } = true;

        /// <summary>
        /// Barrier height at the p contact for the direct calculation (without using metal work fucntion etc.)
        /// </summary>
        public double barrierHeightpContact { get; set; } = 0;
        /// <summary>
        /// Barrier height at the n contact for the direct calculation (without using metal work fucntion etc.)
        /// </summary>
        public double barrierHeightnContact { get; set; } = 0;
        /// <summary>MO
        /// Metal work function of the p contact
        /// </summary>
        public double metalWorkfunctionPContact { get; set; } = 0;// 5.8;
        /// <summary>
        /// Metal work function of the n contact
        /// </summary>
        public double metalWorkfunctionNContact { get; set; } = 0;// 4.1;

        public (double xPosition, double generationRate)[] generationArray { get; set; }


        // Scaling quantities (not used yet)
        public bool scaleSemiconductorEquations { get; set; } = false;
        public double scaleVoltage { get; set; }
        public double scaleSpace { get; set; }
        public double scaleConcentrations { get; set; }
        public double diffusionCoefficient { get; set; }
        public double scaleMobilities { get; set; }
        public double scaleRecombination { get; set; }



        // Arrays for each recombination typ for the semiconductor Loss analysis
        /// <summary>
        /// list with SRH rec current for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double globalSRHCurrent)> voltageDependentSRHCurrentArray { get; set; } = new List<(double voltage, double globalSRHCurrent)>();
        /// <summary>
        /// list with auger rec current for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double globalAugerCurrent)> voltageDependentAugerCurrentArray { get; set; } = new List<(double voltage, double globalAugerCurrent)>();
        /// <summary>
        /// list with radiative rec current for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double globalRadiativeCurrent)> voltageDependentRadiativeCurrentArray { get; set; } = new List<(double voltage, double globalRadiativeCurrent)>();
        /// <summary>
        /// list with SR current of electreons at p contact (minorities) for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double surfaceRecElectronsVop)> voltageDependentSRVelectronVop { get; set; } = new List<(double voltage, double surfaceRecElectronsVop)>();
        /// <summary>
        /// list with SR current of electreons at n contact (majorities) for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double surfaceRecElectronsVzero)> voltageDependentSRVelectronVzero { get; set; } = new List<(double voltage, double surfaceRecElectronsVzero)>();
        /// <summary>
        /// list with SR current of holes at p contact (majorities) for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double surfaceRecHolesVop)> voltageDependentSRVholesVop { get; set; } = new List<(double voltage, double surfaceRecHolesVop)>();
        /// <summary>
        /// list with SR current of holes at n contact (minorities) for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double surfaceRecHolesVzero)> voltageDependentSRVholesVzero { get; set; } = new List<(double voltage, double surfaceRecHolesVzero)>();
        /// <summary>
        /// list with generation current for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double generation)> voltageDependentGeneration { get; set; } = new List<(double voltage, double generation)>();
        /// <summary>
        /// list with sum of all Interface recombination currents within the device for each voltage step of the IV curve
        /// </summary>
        public List<(double voltage, double InterfaceRecCurrent)> voltageDependentInterfaceRecCurrent { get; set; } = new List<(double voltage, double InterfaceRecCurrent)>();

        //public double deltaU = 0.01;


        public MeshingAlgorithm<FiniteElementSemiconductor, RegionSemiconductor> meshingAlgorithm { get; private set; }
        /// <summary>
        /// object for generatrion the Delaunay-grid and the Voronoi-mesh
        /// </summary>
        public Meshing2D_DelaunayVoronoi<FiniteElementSemiconductor, RegionSemiconductor> delaunayVoronoi { get; private set; }
        public Meshing1D_QuasiEquidistant<FiniteElementSemiconductor, RegionSemiconductor> meshing1D_QuasiEquidistant { get; private set; }

        public Mesh<FiniteElementSemiconductor> mesh;

        /// <summary>
        /// voltage, at which the pn junction is operated
        /// </summary>
        public double operatingVoltage { get; set; }

        /// <summary>
        /// characteristic curve of the pn junction
        /// </summary>
        public CharacteristicCurve semiconductorCharacteristic { get; set; }

        /// <summary>
        /// list of indexes of positions, segments and areas, which are set to a certain voltage
        /// </summary>
        public (List<(int index, int selector, List<double> conditions)> points,
            List<(int index, int selector, List<double> conditions)> segments,
            List<(int index, int selector, List<double> conditions)> regions) boundaryConditions
        { get; set; }



        //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        /// <summary>
        /// Contructor, which sets the name, control panels, unit, cornerpoints, regions points and simulation holes from the preferences file
        /// </summary>
        /// <param name="name">Name of the cell</param>
        public ModelSemiconductor(string name, double temperature, bool useSpontRecombination, bool useAugerRecombination, bool useSrhRecombination)
        {
            Console.WriteLine();
            Misc.WriteDividingLine();
            Misc.WriteFormatedLine();
            Misc.WriteFormatedLine(">>> SIMULATION OF A SEMICONDUCTOR JUNCTION <<<");
            T = temperature;
            U_T = physConstants.kB * T / physConstants.e;
            semiconductorCharacteristic = new CharacteristicCurve(T);
            this.name = name;
            this.useRadiativeRecombination = useSpontRecombination;
            this.useAugerRecombination = useAugerRecombination;
            this.useSrhRecombination = useSrhRecombination;

            FermiDiracIntegration(0.5);

        }

        //Fermi-Dirac equations-----------------------------------------------------------------------------------------------------
        public CubicSpline FermiDiracSpline { get; set; }
        public Curve FermiDiracDerivation { get; set; }
        /// <summary>
        /// Calculates the Fermi-Dirac integral of a given order (Standard is 0.5 for densities and -0.5 for derivation of densities) and at a finite number of 
        /// calculation points (eta). Afterwards a Spline with these points is constructed. For the return of values and derivation values this Spline is used. 
        /// </summary>
        /// <param name="integrationOrder"></param>
        public void FermiDiracIntegration(double integrationOrder)
        {
            double eta;
            double[] xPointsFermiDirac = new double[401];
            double[] yPointsFermiDirac = new double[401];
            double integrationResult;

            double FermiDiracIntegrand(double x)
            {
                return 2 / Math.Sqrt(Math.PI) * Math.Pow(x, integrationOrder) / (Math.Exp(x - eta) + 1);
            }
            int j = 0;
            AdaptiveIntegrator adaptiveIntegrator = new AdaptiveIntegrator();
            adaptiveIntegrator.RelativeTolerance = 1e-10;
            for (eta = -80; eta < 20.1; eta += 0.25)
            {

                integrationResult = adaptiveIntegrator.Integrate(FermiDiracIntegrand, 0, double.PositiveInfinity);

                xPointsFermiDirac[j] = eta;
                yPointsFermiDirac[j] = integrationResult;
                j++;
            }
            FermiDiracSpline = new CubicSpline(xPointsFermiDirac, yPointsFermiDirac.Select(z => Math.Log(z)).ToArray(), CubicSplineKind.Akima);
            FermiDiracDerivation = FermiDiracSpline.GetDerivative();

        }




        // Mesh and simulation preparation ------------------------------------------------------------------------------------------
        /// <summary>
        /// Create mesh (Delaunay and Voronoi)
        /// </summary>
        public void SetMesh(string[] geometryLines, int desiredAmountOfPoints, MeshingMethod meshingMethod, Mesh<FiniteElementSemiconductor> meshfromJson, string geometryPath)
        {
            if (meshfromJson == null)
            {

                switch (meshingMethod)
                {

                    case MeshingMethod.quasiEquidistant_1D:
                        GeometryFileData1D geometryFileData1D = new GeometryFileData1D(geometryLines);
                        boundaryConditions = geometryFileData1D.boundaryConditions;
                        meshing1D_QuasiEquidistant = new Meshing1D_QuasiEquidistant<FiniteElementSemiconductor, RegionSemiconductor>(out mesh, geometryFileData1D, desiredAmountOfPoints);
                        meshingAlgorithm = meshing1D_QuasiEquidistant;
                        materialBefore = geometryFileData1D.materialBefore;
                        materialBehind = geometryFileData1D.materialBehind;
                        opticalGeometry = getOpticalGeometry(geometryLines);
                       


                        break;

                    case MeshingMethod.delaunayVoronoi_2D:
                        GeometryFileData2D geometryFileData2D = new GeometryFileData2D(geometryLines);
                        boundaryConditions = geometryFileData2D.boundaryConditions;
                        delaunayVoronoi = new Meshing2D_DelaunayVoronoi<FiniteElementSemiconductor, RegionSemiconductor>(out mesh, geometryFileData2D); // Create four points, which are outside the simulation
                        delaunayVoronoi.ConstructContours(desiredAmountOfPoints); // Construct contour junctions and segments
                        mesh = delaunayVoronoi.AddRegularAndRandomlyShiftedPoints(mesh, 0.5); // Add additional regular points
                        mesh = delaunayVoronoi.AddContours(mesh); // Add points for outer contours, simulation holes and regions
                        mesh = delaunayVoronoi.CreateVoronoiFromDelaunay(mesh); // Create Voronoi mesh from Delaunay triangulation
                        meshingAlgorithm = delaunayVoronoi;
                        materialBefore = geometryFileData2D.materialBefore;
                        materialBehind = geometryFileData2D.materialBehind;
                        opticalGeometry = getOpticalGeometry(geometryLines);

                        Misc.WriteFormatedLine("Delaunay- und Voronoi-Triangulation done with " + mesh.nextAvailableFiniteElementIndex + " Points."); // Output
                        break;

                    case MeshingMethod.quadtree_2D:
                        break;

                    case MeshingMethod.delaunayVoronoi_3D:
                        break;


                }



            }
            else
            {
                mesh = meshfromJson;

                switch (Path.GetExtension(geometryPath))
                {
                    case ".1dg":
                        GeometryFileData1D geometryFileData1D = new GeometryFileData1D(geometryLines);
                        boundaryConditions = geometryFileData1D.boundaryConditions;
                        meshingAlgorithm = new MeshingAlgorithm<FiniteElementSemiconductor, RegionSemiconductor>(geometryFileData1D);
                        Misc.WriteFormatedLine("Mesh loaded with " + mesh.nextAvailableFiniteElementIndex + " Elements.");
                        break;

                    case ".2dg":
                        GeometryFileData2D geometryFileData2D = new GeometryFileData2D(geometryLines);
                        boundaryConditions = geometryFileData2D.boundaryConditions;
                        meshingAlgorithm = new MeshingAlgorithm<FiniteElementSemiconductor, RegionSemiconductor>(geometryFileData2D);
                        Misc.WriteFormatedLine("Mesh loaded with " + mesh.nextAvailableFiniteElementIndex + " Elements.");
                        break;
                }
            }

            
        }

        /// <summary>
        /// refine mesh, calculate new voronoi mesh, solve differential equation again
        /// </summary>
        /// <param name="maxAmountOfRefinements">maximum amount of refinement iterations</param>
        public void RefineMesh(bool refineMesh, int maxAmountOfRefinements, double[] maximumDifferences)
        {
            if (refineMesh)
            {
                // run, till mesh is fine enough or maximum interations are done
                for (int i = 0; i < maxAmountOfRefinements; i++)
                {
                    // Write to Console
                    Misc.WriteFormatedLine("Started refinement step " + (i + 1) + ".");

                    // break when mesh is not refined anymore
                    mesh = delaunayVoronoi.RefineMesh(mesh, new Func<FiniteElementSemiconductor, double>[] { p => p.phi }, maximumDifferences, true,
                        out int amountBefore, out int amountAfter);
                    if (amountBefore == amountAfter)
                    {
                        Misc.WriteFormatedLine("Mesh didn't change.");
                        Misc.WriteDividingLine();
                        break;
                    }

                    // create Voronoi-Mesh
                    mesh = delaunayVoronoi.CreateVoronoiFromDelaunay(mesh);

                    // Write timer to Console
                    Misc.WriteFormatedLine("Delaunay triangulation and Voronoi mesh of refinement step " + (i + 1) + " done.");

                    // set electrical parameters
                    SetBoundaryConditions();
                    foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
                    {
                        point.RefreshParameters(this);
                        point.SetInitialGuessFiniteElement(this);
                    }

                    // solve differential equation via Newtons method
                    Solve();
                }
            }

            // shift points to potentials accoring to Peter Würfel
            // Sketch
            /*
            ——————————————————————————————————————————————————— E_vak = 0
                |                                        |
                |                                        |
                | -e*Φ ≙ electrical potential            |
                |                                        |
                V                                        |
            —————————————————————————————————————————————|————— -eΦ < 0
                |                                        |
                |                                        |
                |          η ≙ electrochemical potential |
                |                                        |
                |                                        |
                | χ = μ < 0 ≙ chemical potential         |
                |                                        |
                V                                        V
            ——————————————————————————————————————————————————— E_c = χ - e*Φ < 0
                |
                |
                | - E_gap ≙ energy gap
                |
                V
            ——————————————————————————————————————————————————— E_v = χ - e*Φ - E_gap < 0
            */

            /*foreach (PointSemiconductor p in mesh.points.Values)
            {
                p.phi *= -1;
                p.phiInit *= -1;
            }
            double offset = mesh.points.Select(p => p.Value.phi).DefaultIfEmpty(1).Max();
            foreach (PointSemiconductor p in mesh.points.Values)
            {
                p.phi -= offset;
                p.phiInit -= offset;
            }*/
        }


        // Output and printing -------------------------------------------------------------------------------------------------------

        /// <summary>
        /// Write all points to file(s).
        /// </summary>
        /// <param name="filepath">filepath of the outputfile</param>
        public void OutputData(string filepath,
            out (double voltage, double current, double loadResistance, double power, double efficiency) simulationResults)
        {
            // delete old file
            if (File.Exists(filepath))
                File.Delete(filepath);

            using (StreamWriter file = new StreamWriter(filepath, true))
            {
                // Header
                file.WriteLine("x\ty\tE_c\tE_v\tphi_n\tphi_p\t-eΦ\tn\tp");
                file.WriteLine("µm\tµm\teV\teV\teV\t1/cm^3\t1/cm^3");

                // Punkte ausgeben
                foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
                    point.OutputToFile(mesh, file, this);
            }

            // print to console —————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
            Misc.WriteFormatedLine($"Simulation done.");
            Misc.WriteFormatedLine();
            Misc.WriteDividingLine();

            simulationResults = (0, 0, 0, 0, 0);
        }

        /// <summary>
        /// Prints all points with their preferences to Console
        /// </summary>
        public void PrintAllPoints()
        {
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine("\n======================== List of all points of the semiconductor junction ========================");
            Console.ForegroundColor = ConsoleColor.Gray;
            foreach (var point in mesh.finiteElements.Values)
                point.Print();
        }




        //Solving including Initial Guess ------------------------------------------------------------------------------------------
        /// <summary>
        /// set initial guess
        /// </summary>
        public void SetInitialGuess(double operatingVoltage)
        {

            this.operatingVoltage = operatingVoltage;

            // set electrical parameters
            SetBoundaryConditions();


            foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
            {
                point.RefreshParameters(this);
            }

            foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
            {
                point.RefreshParametersSecond(this);
            }
            foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
            {
                point.SetInitialGuessFiniteElement(this);
            }

            double phiInitMin = mesh.finiteElements.Values.Min(p => p.phiInit);
            foreach (var t in mesh.finiteElements.Values)
            {
                t.phiInit -= phiInitMin;
                t.phi = t.phiInit;
                t.phi_n = phiInitMin;
                t.phi_p = phiInitMin;
                t.phi_n_init = phiInitMin;
                t.phi_p_init = phiInitMin;
            }

            foreach (var t in mesh.finiteElements.Values)
            {
                //Console.WriteLine(t.index + "\t" + t.phiInit);
                //Console.WriteLine(t.index + "\t" + t.phi);

            }

            setContactBarrierHeights();

        }

        public void SetBoundaryConditions()
        {
            //  ██╗ points
            //  ╚═╝
            foreach (var pointContact in boundaryConditions.points)
            {
                // select point, where the cell contact is in
                foreach (var point in mesh.finiteElements.Values)
                {
                    bool dosomething = false;
                    int dimensions = meshingAlgorithm.dimension;

                    if (dimensions == 2)
                        if (meshingAlgorithm.contourJunctions[pointContact.index].position.InPolygon(point, true))
                            dosomething = true;

                    if (dimensions == 1)
                    {
                        if (point.indexOfBorderElementCreatedFrom == pointContact.index)
                        {
                            dosomething = true;
                        }
                    }

                    if (dosomething)    
                    {
                        point.hasBoundaryCondition = true;

                        point.boundaryCondition = pointContact.selector == 1 ? operatingVoltage : 0;
                        point.hasOperatingVoltage = pointContact.selector == 1 ? true : false;
                        point.contactPreferences = (pointContact.conditions[0] * 1e-2, pointContact.conditions[1] * 1e-2, pointContact.conditions[2], pointContact.conditions[3]);
                        if (pointContact.selector == 2)
                        {
                            point.hasBoundaryCondition = false;
                            point.hasInterfaceCondition = true;
                        }

                        if (dimensions == 2)
                            goto pointFound;
                    }
                }

                // if no point is found, take the nearest one
                int nearestPointIndex = mesh.finiteElements.Values.MinBy(p => meshingAlgorithm.contourJunctions[pointContact.index].position.DistanceTo(p.position)).First().index;
                mesh.finiteElements[nearestPointIndex].hasBoundaryCondition = true;
                mesh.finiteElements[nearestPointIndex].boundaryCondition = pointContact.selector == 1 ? operatingVoltage : 0;
                mesh.finiteElements[nearestPointIndex].hasOperatingVoltage = pointContact.selector == 1 ? true : false;
                mesh.finiteElements[nearestPointIndex].contactPreferences = (pointContact.conditions[0] * 1e-2, pointContact.conditions[1] * 1e-2, pointContact.conditions[2], pointContact.conditions[3]);
            pointFound:;
            }

            //  ██╗ segments
            //  ╚═╝
            foreach (var point in mesh.finiteElements.Values)
                foreach (var segment in this.boundaryConditions.segments)
                {
                    if (meshingAlgorithm.dimension == 2 && point.indexOfBorderElementCreatedFrom == segment.index)
                        {
                            point.hasBoundaryCondition = true;
                            point.boundaryCondition = segment.selector == 1 ? operatingVoltage : 0;
                            point.hasOperatingVoltage = segment.selector == 1 ? true : false;
                            point.contactPreferences = (segment.conditions[0] * 1e-2, segment.conditions[1] * 1e-2, segment.conditions[2], segment.conditions[3]);
                            if (segment.selector == 2)
                            {
                                Console.WriteLine(point.index);
                                point.hasBoundaryCondition = false;
                                point.hasInterfaceCondition = true;
                            }
                        }
                }

            //  ██╗ regions
            //  ╚═╝
            foreach (var point in mesh.finiteElements.Values)
                foreach (var region in this.boundaryConditions.regions)
                    if (point.position.InPolygon(meshingAlgorithm.regions[region.index].orderedPoints.Select(p => p.position).ToList(), true))
                    {
                        point.hasBoundaryCondition = true;
                        point.boundaryCondition = region.selector == 1 ? operatingVoltage : 0;
                        point.hasOperatingVoltage = region.selector == 1 ? true : false;
                        point.contactPreferences = (region.conditions[0] * 1e-2, region.conditions[1] * 1e-2, region.conditions[2], region.conditions[3]);
                        if (region.selector == 2)
                        {
                            point.hasBoundaryCondition = false;
                            point.hasInterfaceCondition = true;
                        }
                    }

        }

        public void SetInitialGuessSingleVoltage(double operatingVoltage)
        {
            this.operatingVoltage = operatingVoltage;

            // set electrical parameters
            SetBoundaryConditions();
            foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
            {
                point.RefreshParameters(this);
            }

            foreach (FiniteElementSemiconductor point in mesh.finiteElements.Values)
            {
                point.RefreshParametersSecond(this);
                point.SetInitialGuessFiniteElement(this);
            }

            double phiInitMin = mesh.finiteElements.Values.Min(p => p.phiInit);
            foreach (var t in mesh.finiteElements.Values)
            {
                t.phiInit -= phiInitMin;
                t.phi = t.phiInit;
                t.phi_n = phiInitMin;
                t.phi_p = phiInitMin;
                t.phi_n_init = phiInitMin + operatingVoltage / 2;
                t.phi_p_init = phiInitMin - operatingVoltage / 2;

                if (t.hasBoundaryCondition)
                    Console.WriteLine("SRV e: " + t.contactPreferences.SRV_electrons);
            }

            deltaU = operatingVoltage;
        }

        /// <summary>
        /// Sets the scaling factors for the sclaing of the drift diffusion equations
        /// </summary>
        public void GetScalingFactors()
        {
            if (scaleSemiconductorEquations)
            {
                scaleVoltage = 1;
                scaleSpace = 1;
                scaleConcentrations = 1;
                diffusionCoefficient = 1;
                scaleMobilities = 1;
                scaleRecombination = 1;

            }
            else
            {
                scaleVoltage = 1;
                scaleSpace = 1;
                scaleConcentrations = 1;
                diffusionCoefficient = 1;
                scaleMobilities = 1;
                scaleRecombination = 1;
            }

        }


        /// <summary>
        /// Solves the differential equation system
        /// </summary>
        public bool Solve()
        {
            // write initial guess from mesh to vector (Newtons method need an vector as input)
            Vector<double> solution = Extreme.Mathematics.Vector.Create<double>(mesh.nextAvailableFiniteElementIndex);

            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                solution[i] = mesh.finiteElements[i].phi;

            // norm of the final residual function of Newtons method
            double residualNorm;

            // Newtons method (solution is initial guess, 'getFunction' and 'getJacobi' are functions/instructions, which are given as parameter
            (solution, residualNorm) = NewtonMethod.Solve(solution, GetFunctionPoissonSync, GetJacobiPoissonSync, 1e-13, 100);

            // return false, if there is no solution found
            if (residualNorm > 1)
                return false;

            // write solution of Newtons method back to mesh
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                mesh.finiteElements[i].phi = solution[i];

            foreach (var p in mesh.finiteElements.Values)
            {
                p.nEquilibrium = p.nDensity(this);
                p.pEquilibrium = p.pDensity(this);
                /*
                if(p.hasBoundaryCondition && p.hasOperatingVoltage)
                    Console.WriteLine("p Side nEq" + p.nEquilibrium + "pEq" + p.pEquilibrium);
                if(p.hasBoundaryCondition && p.hasOperatingVoltage== false)
                    Console.WriteLine("n Side nEq" + p.nEquilibrium + "pEq" + p.pEquilibrium);
                */
            }



            return true;
        }

        public (List<(Material material, double thickness, double roughnessOnTop, bool isAbsorber)> SemiconductorStack, List<(Material material, double thickness)> incoherentBeforeStack,
            List<(Material material, double thickness)> coherentBeforeStack, List<(Material material, double thickness)> incoherentBehindStack,
            List<(Material material, double thickness)> coherentBehindStack) getOpticalGeometry(string[] geometryLines)
        {
            // transition from nm into m (*1e-9) is done here


            if (meshingAlgorithm.dimension == 2)
            {
                // 2d code is in progress
                List<(Material material, double layerthickness, double roughnessOnTop, bool isAbsorber)> layerInformationForTMM
                    = new List<(Material material, double layerthickness, double roughnessOnTop, bool isAbsorber)>();

                    double maxX = -1000;
                    double minX = 1000;
                    double maxY = -1000;
                    double minY = +1000;

                    for (int geomPointIndex = 0; geomPointIndex < delaunayVoronoi.outerContour.orderedPoints.Count; geomPointIndex++)
                    {
                        double geomPoint_X = delaunayVoronoi.outerContour.orderedPoints[geomPointIndex].position.x;
                        double geomPoint_Y = delaunayVoronoi.outerContour.orderedPoints[geomPointIndex].position.y;

                        if (geomPoint_X > maxX)
                            maxX = geomPoint_X;
                        if (geomPoint_X < minY)
                            minX = geomPoint_X;
                        if (geomPoint_Y > maxY)
                            maxY = geomPoint_Y;
                        if (geomPoint_Y < minY)
                            minY = geomPoint_Y;
                    }

                    double xDeviceExpansion = maxX - minX;
                    double yDeviceExpansion = maxY - minY;

                    //Set construction Points for new lines through the device in X and Y direction:
                    Position leftPoint = new Position(minX - 0.2 * xDeviceExpansion, minY + yDeviceExpansion / 2);
                    Position rightPoint = new Position(maxX + 0.2 * xDeviceExpansion, minY + yDeviceExpansion / 2);
                    Position upperPoint = new Position(minX + xDeviceExpansion / 2, maxY + 0.2 * yDeviceExpansion);
                    Position lowerPoint = new Position(minX + xDeviceExpansion / 2, minY - 0.2 * yDeviceExpansion);

                    //Construct Lines:
                    LineSegment LineX = new LineSegment(leftPoint, rightPoint);
                    LineSegment LineY = new LineSegment(lowerPoint, upperPoint);

                    List<Position> xIntersections = new List<Position> { };
                    List<Position> yIntersections = new List<Position> { };

                    for (int segmentIndex = 0; segmentIndex < delaunayVoronoi.contourSegments.Count; segmentIndex++)
                    {
                        if (LineX.DoesIntersect(delaunayVoronoi.contourSegments[segmentIndex].lineSegment))
                        {
                            Position xIntersection = LineX.Intersection(delaunayVoronoi.contourSegments[segmentIndex].lineSegment);
                            xIntersections.Add(xIntersection);
                        }
                        if (LineY.DoesIntersect(delaunayVoronoi.contourSegments[segmentIndex].lineSegment))
                        {
                            Position yIntersection = LineY.Intersection(delaunayVoronoi.contourSegments[segmentIndex].lineSegment);
                            yIntersections.Add(yIntersection);
                        }
                    }

                    //Set Device direction (layerstack in x or y direction)
                    bool layerstackInXdirection = true;
                    if (yIntersections.Count > xIntersections.Count)
                        layerstackInXdirection = false;

                    int amountOfLayers;


                    //Set layer thicknesses and material
                    if (layerstackInXdirection)
                    {
                        amountOfLayers = xIntersections.Count - 1;
                        xIntersections = xIntersections.OrderBy(d => d.x).ToList();

                        for (int intersectIndex = 0; intersectIndex < amountOfLayers; intersectIndex++)
                        {
                            //Set layerthicknesses
                            double layerthickness = xIntersections[intersectIndex + 1].x - xIntersections[intersectIndex].x;
                            //Set Material IDs
                            Position middleOfLayer = new Position(xIntersections[intersectIndex].x + layerthickness / 2, minY + yDeviceExpansion / 2);
                            RegionSemiconductor region = meshingAlgorithm.GetEnclosingRegions(middleOfLayer).Last();
                            layerInformationForTMM.Add((region.material, layerthickness, region.roughnessOnTop, region.isAbsorber));

                            Console.WriteLine("Material: " + layerInformationForTMM[intersectIndex].material.name + " , Tickness: " + layerInformationForTMM[intersectIndex].layerthickness);
                        }
                    }

                    if (layerstackInXdirection == false)
                    {
                        amountOfLayers = yIntersections.Count - 1;
                        yIntersections = yIntersections.OrderBy(d => d.y).ToList();

                        for (int intersectIndex = 0; intersectIndex < amountOfLayers; intersectIndex++)
                        {
                            //Set Layerthicknesses
                            double layerthickness = yIntersections[intersectIndex + 1].y - yIntersections[intersectIndex].y;
                            //Set Material IDs
                            Position middleOfLayer = new Position(yIntersections[intersectIndex].y + layerthickness / 2, minX + xDeviceExpansion / 2);
                            RegionSemiconductor region = meshingAlgorithm.GetEnclosingRegions(middleOfLayer).Last();
                            layerInformationForTMM.Add((region.material, layerthickness, region.roughnessOnTop, region.isAbsorber));

                            Console.WriteLine("Material: " + layerInformationForTMM[intersectIndex].material.name + " , Tickness: " + layerInformationForTMM[intersectIndex].layerthickness);

                        }
                    }

                //Additional layers
                List<double> thicknesses = new List<double>();
                List<(Material material, double thickness)> incoherentBeforeStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> coherentBeforeStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> coherentBehindStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> incoherentBehindStack = new List<(Material material, double thickness)>();


                GeometryFileData1D geometryFileData1D = new GeometryFileData1D(geometryLines);
                for (int i = 0; i < geometryFileData1D.points.Count - 1; i++)
                    thicknesses.Add((geometryFileData1D.points[i + 1] - geometryFileData1D.points[i]) * 1e-9);                
                for (int i = 0; i < geometryFileData1D.incoherentLayersBefore.Count; i++)
                    incoherentBeforeStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.incoherentLayersBefore[i][0]), geometryFileData1D.incoherentLayersBefore[i][1] * 1e-9));
                for (int i = 0; i < geometryFileData1D.coherentLayersBefore.Count; i++)
                {
                    coherentBeforeStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.coherentLayersBefore[i][0]), geometryFileData1D.coherentLayersBefore[i][1] * 1e-9));
                }
                for (int i = 0; i < geometryFileData1D.coherentLayersBehind.Count; i++)
                    coherentBehindStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.coherentLayersBehind[i][0]), geometryFileData1D.coherentLayersBehind[i][1] * 1e-9));
                for (int i = 0; i < geometryFileData1D.incoherentLayersBehind.Count; i++)
                    incoherentBehindStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.incoherentLayersBehind[i][0]), geometryFileData1D.incoherentLayersBehind[i][1] * 1e-9));

                return (layerInformationForTMM, incoherentBeforeStack, coherentBeforeStack, incoherentBehindStack, coherentBehindStack);
            }

            else if (meshingAlgorithm.dimension == 1)
            {
                List<double> thicknesses = new List<double>();
                List<(Material material, double layerthickness, double roughnessOnTop, bool isAbsorber)> semiconductorStack
                    = new List<(Material material, double layerthickness, double roughnessOnTop, bool isAbsorber)>();
                List<(Material material, double thickness)> incoherentBeforeStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> coherentBeforeStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> coherentBehindStack = new List<(Material material, double thickness)>();
                List<(Material material, double thickness)> incoherentBehindStack = new List<(Material material, double thickness)>();


                GeometryFileData1D geometryFileData1D = new GeometryFileData1D(geometryLines);
                for (int i = 0; i < geometryFileData1D.points.Count - 1; i++)
                    thicknesses.Add((geometryFileData1D.points[i + 1] - geometryFileData1D.points[i]) * 1e-9);
                for (int i = 0; i < geometryFileData1D.materials.Count; i++)
                    semiconductorStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.materials[i][0]), thicknesses[i], geometryFileData1D.materials[i][2],
                        (int)(geometryFileData1D.materials[i][1]) == 1 ? true : false));
                for (int i = 0; i < geometryFileData1D.incoherentLayersBefore.Count; i++)
                    incoherentBeforeStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.incoherentLayersBefore[i][0]), geometryFileData1D.incoherentLayersBefore[i][1] * 1e-9));
                for (int i = 0; i < geometryFileData1D.coherentLayersBefore.Count; i++)
                {
                    coherentBeforeStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.coherentLayersBefore[i][0]), geometryFileData1D.coherentLayersBefore[i][1] * 1e-9));
                }
                for (int i = 0; i < geometryFileData1D.coherentLayersBehind.Count; i++)
                    coherentBehindStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.coherentLayersBehind[i][0]), geometryFileData1D.coherentLayersBehind[i][1] * 1e-9));
                for (int i = 0; i < geometryFileData1D.incoherentLayersBehind.Count; i++)
                    incoherentBehindStack.Add((Data.GetMaterialFromID((int)geometryFileData1D.incoherentLayersBehind[i][0]), geometryFileData1D.incoherentLayersBehind[i][1] * 1e-9));


                return (semiconductorStack, incoherentBeforeStack, coherentBeforeStack, incoherentBehindStack, coherentBehindStack);

            }

            else
                return (null, null, null, null, null);

        }

        public void setReferenceSpectrum(double spectrumStart = 300e-9, double spectrumEnd = 1300e-9)
        {
            referenceSpectrum = new Spectrum(MiscTMM.spectrumAM15.data.Where(d => d.lambda >= spectrumStart && d.lambda < spectrumEnd).Select(d => (d.lambda, d.deltaLambda, d.spectralIntensityDensity /**10/totalPower*/)).ToArray());

        }


        /// <summary>
        /// Sets initial guess of the quasi fermi levels in case of illumination
        /// </summary>
        /// <param name="enableGeneration"></param>
        public void SetInitialGuessIllumination(bool enableGeneration, bool enableSrhRecombination, bool enableAugerRecombination,
            bool enableRadiativeRecombination, OpticModeSemiconductor opticModeSemiconductor, double constantGeneration, double spectrumStart = 300e-9, double spectrumEnd = 1300e-9)
        {
            totalInputPower = MiscTMM.spectrumAM15.data.Where(d => d.lambda >= spectrumStart && d.lambda < spectrumEnd).Sum(s => s.deltaLambda * s.spectralIntensityDensity);
            
            //Set spectrum
            Spectrum spectrum = new Spectrum(MiscTMM.spectrumAM15.data.Where(d => d.lambda >= spectrumStart && d.lambda < spectrumEnd).Select(d => (d.lambda, d.deltaLambda, d.spectralIntensityDensity /**10/totalPower*/)).ToArray());


            //Set absorber preferences
            if (meshingAlgorithm.regions.Where(r => r.isAbsorber).Count() > 0)
            {
                foreach (var layer in opticalGeometry.SemiconductorStack)
                {
                    if (meshingAlgorithm.regions.Where(r => r.isAbsorber).First().material.ID == layer.material.ID)
                    {
                        absorberID = layer.material.ID;
                        absorberThickness = layer.thickness;
                        absorberBandGap = layer.material.propertiesSemiconductor.Egap;
                    }
                }
            }
            else
            {
                MessageBox.Show("No absorber/actice layer defined. Please adapt geometry/input file.", "Geometry File", MessageBoxButton.OK, MessageBoxImage.Information);
            }



            var TMMstructure = getTMMstructure();
            var materialStack = TMMstructure.materialStack;// electircal stack together with coherent layers before and after
            double generationDisplacement = TMMstructure.coherentLayerThicknessBefore; // due to additional coherent layers before the actual semiconductor stack


            // TOdo: werden hier und in Loss Analysis immer die richtigen SPektren, Lambda grenzen usw. verwendet???
            Spectrum spectrumAfterIncoh = calculateIncoherentLayersBefore(spectrum).spectrumAfterIncohBefore; //Original SPectrum - (R +) A in incoherent layer
            absorptanceIncohBefore = calculateIncoherentLayersBefore(spectrum).absorbedPhotonsRelative;
            absorbedPowerIncohBefore = calculateIncoherentLayersBefore(spectrum).absorbedPowerPerLayer;

            calculateIncoherentLayersBehind(spectrumAfterIncoh); //in case of incoherent material behind: set new material behind + transmission from TMM stack is reduced 

            Material materialBeforeStack = Data.GetMaterialFromID(materialBefore);
            (Material material, double roughness) materialBehindStack = (Data.GetMaterialFromID(materialBehind.ID), materialBehind.roughnessOnTop * 1e-9);
            modelOpticsTMM = new ModelTMM(materialBeforeStack, materialBehindStack, materialStack, spectrumAfterIncoh);

            Console.WriteLine(" ------------------------------------------------------> material Before: " + modelOpticsTMM.layerBeforeStack.material.name);
            this.enableGeneration = enableGeneration;
            this.useSrhRecombination = enableSrhRecombination;
            this.useAugerRecombination = enableAugerRecombination;
            this.useRadiativeRecombination = enableRadiativeRecombination;

            foreach (var f in mesh.finiteElements.Values)
                if (meshingAlgorithm.GetEnclosingRegions(f.position).First().isAbsorber)
                    f.isFEinAbsorber = true;

            if (enableGeneration == true)
            {
                generationArray = new (double xPosition, double generationRate)[mesh.nextAvailableFiniteElementIndex];
                switch (opticModeSemiconductor)
                {
                    case OpticModeSemiconductor.TMM:
                        foreach (var p in mesh.finiteElements.Values)
                        {
                            generationArray[p.index] = (p.position.x, modelOpticsTMM.GetLocalGeneration(p.position.x+generationDisplacement, spectrumStart, spectrumEnd));
                        }
                        break;
                    case OpticModeSemiconductor.AbsorberConstantGeneration:
                        foreach (var p in mesh.finiteElements.Values)
                        {
                            if (p.material.ID == absorberID)
                                generationArray[p.index] = (p.position.x, constantGeneration);
                            else
                                generationArray[p.index] = (p.position.x, 1e22);

                        }
                        break;
                    case OpticModeSemiconductor.noParasiticAbsorbtion:
                        double startAbsorberLayer = meshingAlgorithm.regions.Where(d => d.isAbsorber).First().orderedPoints.MinBy(p => p.position.x).First().position.x; //meshingAlgorithm.regions.Where(d => d.material.ID == 020000000 || d.material.ID == 020005000 || d.material.ID == 020006000).First().orderedPoints.MinBy(p => p.position.x).First().position.x;
                        double localGeneration;
                        foreach (var p in mesh.finiteElements.Values)
                        {

                            if (p.material.ID == absorberID)
                                localGeneration = (1.354e28 * Math.Exp(-6167586 * (p.position.x - startAbsorberLayer)) + 2.697e26);
                            else
                                localGeneration = 1e22;
                            generationArray[p.index] = (p.position.x, localGeneration);
                        }
                        break;
                    default:
                        break;
                }



                Console.WriteLine("\tIllumination enabled.");
            }
            else
                Console.WriteLine("\tIllumination disabled.");
        }

        /// <summary>
        /// Calculate lambert beer in incoherent layers before the actual TMM layerstack and reduces the spectrum for each wavelength by this factor
        /// Set last incoherent layer (= adjacent layer to the TMM layerstack)
        /// </summary>
        public (Spectrum spectrumAfterIncohBefore, double[] transmittance, List<double> absorbedPowerPerLayer, List<double> absorbedPhotonsRelative) calculateIncoherentLayersBefore(Spectrum spectrum)
        {
            double[] lambertBeerTransmittance = new double[spectrum.data.Length];
            for (int i = 0; i < lambertBeerTransmittance.Length; i++)
                lambertBeerTransmittance[i] = 1;


            List<double> absorbedPowerIncohPerLayer = new List<double>();
            List<double> totalAbsorbedIncohRelative = new List<double>();


            // in case of additional incoherent layer behind the TMM layerstack: change material behind.
            if (opticalGeometry.incoherentBehindStack.Count > 0)
            {
                materialBehind = (opticalGeometry.incoherentBehindStack.First().material.ID, 0e-9);

            }


            if (opticalGeometry.incoherentBeforeStack.Count > 0)
            {
                //Set new material before
                materialBefore = opticalGeometry.incoherentBeforeStack.Last().material.ID;
                
                Spectrum spectrumAfterIncoh = new Spectrum();

                (double wavelgth, double deltaWaveltgh, double specIntDensity)[] changedSpectrumData = new (double wavelgth, double deltaWaveltgh, double specIntDensity)[spectrum.data.Length];
                // Calculate wavelength dependent absoprtion in each incoherent layer 
                foreach (var layer in opticalGeometry.incoherentBeforeStack)
                {

                    double incomingAbsoultePhotons = 0;
                    double absorbedIncohPhotons = 0;

                    for (int specIndex = 0; specIndex < spectrum.data.Length; specIndex++)
                    {
                        double transmittedLambertBeer = Misc.LambertBeerTransmittedRatio(spectrum.data[specIndex].lambda, layer.thickness,
                            layer.material.propertiesOptics.n_rawData.MinBy(d => Math.Abs(d.lambda - spectrum.data[specIndex].lambda)).First().n.Im);
                        lambertBeerTransmittance[specIndex] *= transmittedLambertBeer;

                        //total amount of photons i
                        incomingAbsoultePhotons += spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda * spectrum.data[specIndex].lambda;

                        //total amount of absorbed photons in incho layer 
                        absorbedIncohPhotons += spectrum.data[specIndex].spectralIntensityDensity * spectrum.data[specIndex].deltaLambda * spectrum.data[specIndex].lambda * (1-transmittedLambertBeer);
                    }

                    // ratio of absorbed and incoming photons
                    totalAbsorbedIncohRelative.Add( absorbedIncohPhotons / incomingAbsoultePhotons);

                    //power loss  = ratio * totalPower (in Spectrum)
                    double totalPower = spectrum.totalPower;

                    double absorbedIncohPower = totalPower * totalAbsorbedIncohRelative.Last();

                    absorbedPowerIncohPerLayer.Add(absorbedIncohPower);
                }

                for (int j = 0; j < changedSpectrumData.Length; j++)
                    changedSpectrumData[j] = (spectrum.data[j].lambda, spectrum.data[j].deltaLambda, spectrum.data[j].spectralIntensityDensity * lambertBeerTransmittance[j]);

                spectrumAfterIncoh = new Spectrum(changedSpectrumData);

                return (spectrumAfterIncoh, lambertBeerTransmittance, absorbedPowerIncohPerLayer, totalAbsorbedIncohRelative);
            }
            else 
                return (spectrum, lambertBeerTransmittance, absorbedPowerIncohPerLayer, totalAbsorbedIncohRelative);
        }
        /// <summary>
        /// Set the TMM structure for optical calculation of coherent layers. Consits of the actual semiconductor layerstack plus the coherent layers before and after.
        /// </summary>
        ((Material material, double thickness, double roughnessOnTop)[] materialStack, double coherentLayerThicknessBefore) getTMMstructure()
        {
                int lengthTMMstack = opticalGeometry.coherentBeforeStack.Count
                + opticalGeometry.SemiconductorStack.Count
                + opticalGeometry.coherentBehindStack.Count;
                (Material material, double thickness, double roughnessOnTop)[] materialStack =
                    new (Material material, double thickness, double roughnessOnTop)[lengthTMMstack];
                double coherentThicknessBefore = 0; // Sum of all thicknesses of coherent layers before the semiconductor stack (for the displacement of the generation rate)


                int countCoh1 = opticalGeometry.coherentBeforeStack.Count;
                int countSC = opticalGeometry.SemiconductorStack.Count;
                int countCoh2 = opticalGeometry.coherentBehindStack.Count;
                for (int i = 0; i < countCoh1; i++)
                {
                    materialStack[i] = (opticalGeometry.coherentBeforeStack[i].material, opticalGeometry.coherentBeforeStack[i].thickness, 0);
                    coherentThicknessBefore += opticalGeometry.coherentBeforeStack[i].thickness;
                }
                for (int i = 0; i < countSC; i++)
                {
                    materialStack[i + countCoh1] = (opticalGeometry.SemiconductorStack[i].material, opticalGeometry.SemiconductorStack[i].thickness, opticalGeometry.SemiconductorStack[i].roughnessOnTop * 1e-9);
                }
                for (int i = 0; i < countCoh2; i++)
                {
                    materialStack[i + countCoh1 + countSC] = (opticalGeometry.coherentBehindStack[i].material, opticalGeometry.coherentBehindStack[i].thickness, 0);
                }

                return (materialStack, coherentThicknessBefore);

            
        }
        /// <summary>
        /// Calculate the lambert beer absorption in the incoherent layers after the actual TMM layerstack
        /// Set first incoherent layer (= adjacent layer to the TMM layerstack) as new material behind
        /// </summary>
        public (Spectrum spectrumAfterIncohBefore, double[] transmittance, List<double> absorbedPowerPerLayer, List<double> absorbedPhotonsRelative) calculateIncoherentLayersBehind(Spectrum spectrumAfterIncoh)
        {
            (double lambda, double deltaLambda, double spectralPowerDensity)[] newSPectrumData = new (double lambda, double deltaLambda, double spectralPowerDensity)[spectrumAfterIncoh.data.Length];
            for (int i = 0; i > spectrumAfterIncoh.data.Length; i++)
            {
                newSPectrumData[i] = (spectrumAfterIncoh.data[i].lambda, spectrumAfterIncoh.data[i].deltaLambda, spectrumAfterIncoh.data[i].spectralIntensityDensity * modelOpticsTMM.T[i]);
            }
            Spectrum spectrumAfterTMM = new Spectrum(newSPectrumData);

            double[] lambertBeerTransmittance = new double[spectrumAfterTMM.data.Length];
            for (int i = 0; i < lambertBeerTransmittance.Length; i++)
                lambertBeerTransmittance[i] = 1;


            List<double> absorbedPowerIncohPerLayer = new List<double>();
            List<double> totalAbsorbedIncohRelative = new List<double>();

            if (opticalGeometry.incoherentBehindStack.Count > 0)
            {

                (double wavelgth, double deltaWaveltgh, double specIntDensity)[] changedSpectrumData = new (double wavelgth, double deltaWaveltgh, double specIntDensity)[spectrumAfterTMM.data.Length];
                // Calculate wavelength dependent absoprtion in each incoherent layer 
                foreach (var layer in opticalGeometry.incoherentBeforeStack)
                {

                    double incomingAbsoultePhotons = 0;
                    double absorbedIncohPhotons = 0;

                    for (int specIndex = 0; specIndex < spectrumAfterTMM.data.Length; specIndex++)
                    {
                        double transmittedLambertBeer = Misc.LambertBeerTransmittedRatio(spectrumAfterTMM.data[specIndex].lambda, layer.thickness,
                            layer.material.propertiesOptics.n_rawData.MinBy(d => Math.Abs(d.lambda - spectrumAfterTMM.data[specIndex].lambda)).First().n.Im);
                        lambertBeerTransmittance[specIndex] *= transmittedLambertBeer;

                        //total amount of photons i
                        incomingAbsoultePhotons += spectrumAfterTMM.data[specIndex].spectralIntensityDensity * spectrumAfterTMM.data[specIndex].deltaLambda * spectrumAfterTMM.data[specIndex].lambda;

                        //total amount of absorbed photons in incho layer 
                        absorbedIncohPhotons += spectrumAfterTMM.data[specIndex].spectralIntensityDensity * spectrumAfterTMM.data[specIndex].deltaLambda * spectrumAfterTMM.data[specIndex].lambda * (1 - transmittedLambertBeer);
                    }

                    // ratio of absorbed and incoming photons
                    totalAbsorbedIncohRelative.Add(absorbedIncohPhotons / incomingAbsoultePhotons);
                    Console.WriteLine("Absorbed Photons: " + absorbedIncohPhotons / physConstants.h / physConstants.c);
                    Console.WriteLine("Total Photons: " + incomingAbsoultePhotons / physConstants.h / physConstants.c);

                    //power loss  = ratio * totalPower (in Spectrum)
                    double totalPower = spectrumAfterTMM.totalPower;

                    double absorbedIncohPower = totalPower * totalAbsorbedIncohRelative.Last();

                    absorbedPowerIncohPerLayer.Add(absorbedIncohPower);
                }

                for (int j = 0; j < changedSpectrumData.Length; j++)
                    changedSpectrumData[j] = (spectrumAfterTMM.data[j].lambda, spectrumAfterTMM.data[j].deltaLambda, spectrumAfterTMM.data[j].spectralIntensityDensity * lambertBeerTransmittance[j]);

                Spectrum spectrumAfterIncohBehind = new Spectrum(changedSpectrumData);

                return (spectrumAfterIncohBehind, lambertBeerTransmittance, absorbedPowerIncohPerLayer, totalAbsorbedIncohRelative);
            }
            else
                return (spectrumAfterTMM, lambertBeerTransmittance, absorbedPowerIncohPerLayer, totalAbsorbedIncohRelative);
        }
        

        public (Position position, double phi_electron, double phi_hole, double phi_potential)[] SaveInitialGuessIllumination()
        {
            (Position position, double phi_electron, double phi_hole, double phi_potential)[] startingValuesIllumination
                = new (Position position, double phi_electron, double phi_hole, double phi_potential)[mesh.nextAvailableFiniteElementIndex];

            for (int FE = 0; FE < mesh.nextAvailableFiniteElementIndex; FE++)
                startingValuesIllumination[FE] = ((mesh.finiteElements[FE].position), mesh.finiteElements[FE].phi_n, mesh.finiteElements[FE].phi_p, mesh.finiteElements[FE].phi);

            return startingValuesIllumination;
        }

        public void SetInitialGuessIllumination((Position position, double phi_electron, double phi_hole, double phi_potential)[] startingValuesIllumination
            , ModelTMM modelOptics, bool enableGeneration, bool enableSrhRecombination, bool enableAugerRecombination, bool enableRadiativeRecombination)
        {
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("\nStarte Initial Guess Illumination");
            Console.ForegroundColor = ConsoleColor.Gray;
            this.enableGeneration = enableGeneration;
            this.useSrhRecombination = enableSrhRecombination;
            this.useAugerRecombination = enableAugerRecombination;
            this.useRadiativeRecombination = enableRadiativeRecombination;
            if (enableGeneration == true)
            {
                generationArray = new (double xPosition, double generationRate)[mesh.nextAvailableFiniteElementIndex];
                //Generate optic model, write generation values in array

                foreach (var p in mesh.finiteElements.Values)
                    generationArray[p.index] = (p.position.x, modelOptics.GetLocalAbsorption(p.position.x, 300e-9, 1080e-9));


                foreach (var b in mesh.finiteElements.Values)
                {
                    var minDistancePoint = startingValuesIllumination.MinBy(e => e.position.DistanceSquaredTo(b.position)).First();
                    b.phi_n = minDistancePoint.phi_electron;
                    b.phi_p = minDistancePoint.phi_hole;
                    //b.phi = minDistancePoint.phi_potential;
                }

                resetRampingFactor = true;
                Console.WriteLine("Illumination enabled.");

            }
            else
            {
                Console.WriteLine("Illumination disabled.");
            }
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Initial Guess Illumination beendet");
            Console.ForegroundColor = ConsoleColor.Gray;
        }

        public void setContactBarrierHeights()
        {
            foreach (var p in mesh.finiteElements.Values)
            {
                if (p.hasBoundaryCondition)
                {
                    p.barrierHeight = p.contactPreferences.contactBarrier;

                    /*
                    if (p.hasOperatingVoltage) // p-contact
                    {
                        if (Math.Abs(metalWorkfunctionPContact) > 1e-3) // Calculate barrier height   via work function of the used metal
                            p.barrierHeight = -p.phiInit - metalWorkfunctionPContact - p.phi_p_init;
                        else
                            p.barrierHeight = barrierHeightpContact; // takes given barrier height for p contact as difference to initial quasi fermi level
                        
                    }
                    else //n-contact
                    {
                        if (Math.Abs(metalWorkfunctionNContact) > 1e-3) // Calculate barrier height for n contact via work function of the used metal
                        { 
                            p.barrierHeight = -p.phiInit -  metalWorkfunctionNContact - p.phi_n_init;
                        }
                        else
                            p.barrierHeight = barrierHeightnContact;  // takes given barrier height for n contact as difference to initial quasi fermi level
                       
                    }
                */
                }
            }
        }

        /// <summary>
        /// Solves the Van Roosbroeck System
        /// </summary>
        /// <param name="operatingVoltage"></param>
        public void SolvingVrb(double operatingVoltage)
        {
            //Tolerance für das Newtonverfahren angeben:███████████████████████████████████████████████████
            double toleranceVrb = 1e-6;

            semiconductorCharacteristic = new CharacteristicCurve(T);

            int Usteps = (int)Math.Abs(operatingVoltage / deltaU); // Number of voltage steps
            if (operatingVoltage == 0)
                Usteps = 0;


            if (operatingVoltage < 0)
                deltaU *= -1;

            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("\nStarte van Roosbroeck-Newtonverfahren");
            Console.ForegroundColor = ConsoleColor.Gray;

            //SCALING
            //ScaleValues();

            double quasiJsc = 0;

            // █████████████████████████████████████████████████████████████████████████████████████████████
            //Eigentliches Newtonverfahren mit Schleife über die Spannung.
            for (int Ucounter = 0; Ucounter < Usteps + 1; Ucounter++)
            {
                //Creating vectors for solving Van Roosbrock problem ██████████████████████████████████████████
                Vector<double> phiPsinPsip = Vector.Create<double>(3 * mesh.nextAvailableFiniteElementIndex);
                Vector<double> phiPsinPsipInitial = Vector.Create<double>(3 * mesh.nextAvailableFiniteElementIndex);

                // Set boundary voltages
                foreach (var p in mesh.finiteElements.Values)
                {

                    if (operatingVoltage == 0)
                        p.currentBoundaryCondition = 0;
                    else
                        p.currentBoundaryCondition = Ucounter * deltaU * p.boundaryCondition / operatingVoltage;

                    if (p.hasBoundaryCondition)
                    {
                        //Console.WriteLine("Punkt mit Boundary Condition: " + p.index);
                        //Console.WriteLine("Boundary Condition von Punkt " + p.index + ": " + p.boundaryCondition);
                        //Console.WriteLine("CURRENT Boundary Condition von Punkt " + p.index + ": " + p.currentBoundaryCondition);
                    }
                }

                int generationRampingSteps;

                if (enableGeneration && Ucounter == 0)
                    generationRampingSteps = 60;
                else
                    generationRampingSteps = 0;

                if (resetRampingFactor && Ucounter == 0)
                    generationRampingSteps = 20;

                for (double generationCounter = 0; generationCounter < generationRampingSteps + 1; generationCounter++)
                {

                    double generationRampingFactor = Math.Pow(10, generationRampingSteps - (generationCounter));

                    globalGenerationRampingFactor = generationRampingFactor;
                    Console.WriteLine("globalGenerationRampingFactor = " + globalGenerationRampingFactor);

                    //Vektoren aus Dict befüllen
                    int k = 0;
                    for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    {
                        if (mesh.finiteElements.ContainsKey(i))
                        {
                            phiPsinPsip[k] = mesh.finiteElements[i].phi;
                            phiPsinPsip[k + mesh.nextAvailableFiniteElementIndex] = mesh.finiteElements[i].phi_n;
                            phiPsinPsip[k + 2 * mesh.nextAvailableFiniteElementIndex] = mesh.finiteElements[i].phi_p;
                            k++;
                        }
                    }


                    //Newton lösen███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
                    phiPsinPsip = NewtonMethod.Solve(phiPsinPsip, GetFunctionVanRoosbroeckSync, GetJacobiVanRoosbroeckAsync, toleranceVrb, 100).solution;

                    //Werte wieder in Dictionary schreiben█████████████████
                    for (int i = 0; i < mesh.finiteElements.Keys.Max() + 1; i++)
                    {
                        mesh.finiteElements[i].phi = phiPsinPsip[i];
                        mesh.finiteElements[i].phi_n = phiPsinPsip[i + mesh.nextAvailableFiniteElementIndex];
                        mesh.finiteElements[i].phi_p = phiPsinPsip[i + 2 * mesh.nextAvailableFiniteElementIndex];

                    }



                }
                setDeviceCurrents(Ucounter * deltaU);



                //Set Recombination Currents for every voltage step
                setLocalGenRecCurrents(this, Ucounter * deltaU);
                setGlobalLossCurrents(Ucounter * deltaU);
                SetLocalCurrentsForLastVoltage();
                writeBotPotData(Ucounter * deltaU);

                if (Ucounter == 0)
                    quasiJsc = semiconductorCharacteristic.experimentalData[0].current;
                else
                {
                    if (Math.Abs(semiconductorCharacteristic.experimentalData.Last().current) > 1000)
                        goto BreakVoltageLoop;
                    if (stopAfterVoc)
                        if (semiconductorCharacteristic.experimentalData.Last().current > stopAfterVocFactorIph * Math.Abs(quasiJsc) && enableGeneration)
                            goto BreakVoltageLoop;

                }
            }
        BreakVoltageLoop:;


            SetLocalCurrentsForLastVoltage();


            // Roosbroeck Verfahren beendet██████████████████████████████████████████████████████████████████
            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine("Roosbroeck-Newtonverfahren beendet");
            Console.ForegroundColor = ConsoleColor.Gray;

        }


        // SCaling, Simulation preparation and simulation processing--------------------------------------------------------------------------------

        public void SetLocalCurrentsForLastVoltage()
        {
            //Set currents from point to each neighbor
            foreach (var n in mesh.finiteElements.Values)
            {
                n.SetLocalCurrents(this);
            }
        }

        public void writeBotPotData(double voltage)
        {
            // Streamwriter Semiconductor Data
            string filepath = InputOutput.pathSemiconductor.output + "semiconductorBoTPotData";
            filepath += ".dat";

            if (voltage == 0)
                File.Delete(filepath);

            using (StreamWriter file = new StreamWriter(filepath, true))
            {
                string line = null;


                if (new FileInfo(filepath).Length == 0)
                {
                    file.WriteLine("voltage\tindex\tPosition\tEc\tEv\tEfn\tEfp\tjn\tjp\tjG\tjR");
                    file.WriteLine("V\tm\teV\teV\teV\teV\tA/m^2\tA/m^2\tA/m^3\tA/m^3");
                }
                foreach (var v in mesh.finiteElements.Values)
                {
                    double Ec = -v.phi + v.material.propertiesSemiconductor.chemicalPotential;
                    double Ev = -v.phi + v.material.propertiesSemiconductor.chemicalPotential - v.localBandGap;
                    double phin = v.phi_n;
                    double phip = v.phi_p;

                    string datastring = InputOutput.ToStringWithSeparator(voltage) + "\t" +
                        InputOutput.ToStringWithSeparator(v.index) + "\t"
                     + InputOutput.ToStringWithSeparator(v.position.x) + "\t"
                     + InputOutput.ToStringWithSeparator(Ec) + "\t"
                     + InputOutput.ToStringWithSeparator(Ev) + "\t"
                     + InputOutput.ToStringWithSeparator(phin) + "\t"
                     + InputOutput.ToStringWithSeparator(phip) + "\t"
                     + InputOutput.ToStringWithSeparator(v.electronCurrent.First()) + "\t"
                     + InputOutput.ToStringWithSeparator(v.holeCurrent.First()) + "\t"
                     + InputOutput.ToStringWithSeparator(v.TotalGenerationRate(this)) + "\t"
                     + InputOutput.ToStringWithSeparator(v.TotalRecombinationRate(this))
                     ;
                    file.WriteLine(datastring);
                }
            }
        }
        public void setDeviceCurrents(double voltage)
        {
            double sumCurrent = 0;
            double sumEdge = 0;
            double sumEdgeZeroVolts = 0;

            //foreach (var p in mesh.finiteElements.Values.Where(a => a.hasOperatingVoltage))
            foreach (var p in mesh.finiteElements.Values.Where(a => a.hasBoundaryCondition && a.hasOperatingVoltage == false))
            {
                foreach (var n in p.neighbors)
                {
                    sumCurrent += jElectron(p.index, n.index) * n.edgeSize + jHole(p.index, n.index) * n.edgeSize;
                }
                sumEdge += p.borderEdgeSize;

            }

            // For geometries with passivation pattern (area with contact/operting voltage is smaller than cell area)
            foreach (var p in mesh.finiteElements.Values.Where(a => a.hasBoundaryCondition))
            {

                if (p.hasOperatingVoltage == false)
                {
                    sumEdgeZeroVolts += p.borderEdgeSize;
                }
            }


            if (sumEdgeZeroVolts > sumEdge)
                sumEdge = sumEdgeZeroVolts;

            double sumCurrentDensity = -1 * sumCurrent / sumEdge;

            semiconductorCharacteristic.AddExperimentalPoint((voltage, sumCurrentDensity, voltage * sumCurrentDensity, 1, 1));
        }

        // funtions for calculation of recombination currents
        /// <summary>
        /// sets the recombination currents for each meshpoint (as a voltage dependent arrray at each point) for every recombination mechanism and generation. 
        /// </summary>
        /// <param name="voltage">voltage which is to be set</param>
        public void setLocalGenRecCurrents(ModelSemiconductor modelSemiconductor, double voltage)
        {
            if (enableGeneration)
            {
                foreach (var p in mesh.finiteElements.Values)
                {
                    p.voltageDependentLocalGenerationArray.Add((voltage, p.TotalGenerationRate(modelSemiconductor) * physConstants.e * p.size));
                    p.voltageDependentLocalSRHRecombinationArray.Add((voltage, p.SRHRecombinationRate(modelSemiconductor) * physConstants.e * p.size));
                    p.voltageDependentLocalAugerRecombinationArray.Add((voltage, p.AugerRecombinationRate(modelSemiconductor) * physConstants.e * p.size));
                    p.voltageDependentLocalRadiativeRecombinationArray.Add((voltage, p.SpontaneousRecombinationRate(modelSemiconductor) * physConstants.e * p.size));

                    if (p.hasBoundaryCondition) //Boundary Point?
                    {
                        if (p.hasOperatingVoltage)
                        {
                            p.voltageDependentLocalSRVelectronVop.Add((voltage, p.ElectronSurfaceRecombinationCurrent(modelSemiconductor) * p.borderEdgeSize));
                            p.voltageDependentLocalSRVholeVop.Add((voltage, p.HoleSurfaceRecombinationCurrent(modelSemiconductor) * p.borderEdgeSize));

                        }
                        else
                        {
                            p.voltageDependentLocalSRVelectronVzero.Add((voltage, p.ElectronSurfaceRecombinationCurrent(modelSemiconductor) * p.borderEdgeSize));
                            p.voltageDependentLocalSRVholeVzero.Add((voltage, p.HoleSurfaceRecombinationCurrent(modelSemiconductor) * p.borderEdgeSize));
                        }
                    }
                    if (p.hasInterfaceCondition)
                    {
                        p.voltageDependentLocalInterfaceRecombinationArray.Add((voltage, p.InterfaceRecombinationRate(modelSemiconductor) * physConstants.e * p.size));
                    }
                }
            }
        }

        /// <summary>
        /// sets the global loss currents for each recombination mechansim (in public voltage dependent lists)  Currents in A/m^2. Recalculation in mA/cm^2 while plotting.
        /// </summary>
        /// <param name="voltage">voltage which is to be set</param>
        public void setGlobalLossCurrents(double voltage)
        {
            if (enableGeneration)
            {
                double SRH = 0;
                double auger = 0;
                double radiative = 0;
                double SRVelectronsVop = 0;
                double SRVelectronsVzero = 0;
                double SRVholesVop = 0;
                double SRVholesVzero = 0;
                double generation = 0;
                double contactAreaVop = 0;
                double contactAreaVzero = 0;
                double interfaceRec = 0;

                foreach (var p in mesh.finiteElements.Values)
                {
                    SRH += p.voltageDependentLocalSRHRecombinationArray.Last().localSRHRecombinationCurrent;
                    auger += p.voltageDependentLocalAugerRecombinationArray.Last().localAugerRecombinationCurrent;
                    radiative += p.voltageDependentLocalRadiativeRecombinationArray.Last().localRadiativeRecombinationCurrent;
                    generation += p.voltageDependentLocalGenerationArray.Last().localGenerationCurrent;

                    if (p.hasBoundaryCondition)
                    {
                        if (p.hasOperatingVoltage)
                        {
                            SRVelectronsVop += p.voltageDependentLocalSRVelectronVop.Last().localSRVelectronVop;
                            SRVholesVop += p.voltageDependentLocalSRVholeVop.Last().localSRVholeVop;
                            contactAreaVop += p.borderEdgeSize;
                        }
                        else
                        {
                            SRVelectronsVzero += p.voltageDependentLocalSRVelectronVzero.Last().localSRVelectronVzero;
                            SRVholesVzero += p.voltageDependentLocalSRVholeVzero.Last().localSRVholeVzero;
                            contactAreaVzero += p.borderEdgeSize;

                        }
                    }
                    if (p.hasInterfaceCondition)
                    {
                        interfaceRec += p.voltageDependentLocalInterfaceRecombinationArray.Last().localSRVholeVzero;
                    }


                }

                //active area, in most cases equal to contactareas, but not for geometries with passivation features like point contacts through a passivation layer
                double activeArea = 0;
                if (contactAreaVop > contactAreaVzero)
                    activeArea = contactAreaVop;
                else
                    activeArea = contactAreaVzero;


                //Console.WriteLine("Vop area" + contactAreaVop);
                //Console.WriteLine("Vzero area" + contactAreaVzero);
                voltageDependentSRHCurrentArray.Add((voltage, SRH / activeArea));
                voltageDependentAugerCurrentArray.Add((voltage, auger / activeArea));
                voltageDependentRadiativeCurrentArray.Add((voltage, radiative / activeArea));
                voltageDependentGeneration.Add((voltage, generation / activeArea));
                voltageDependentSRVelectronVop.Add((voltage, SRVelectronsVop / activeArea));
                voltageDependentSRVelectronVzero.Add((voltage, SRVelectronsVzero / activeArea));
                voltageDependentSRVholesVop.Add((voltage, SRVholesVop / activeArea));
                voltageDependentSRVholesVzero.Add((voltage, SRVholesVzero / activeArea));
                voltageDependentInterfaceRecCurrent.Add((voltage, interfaceRec / activeArea));





            }
        }


        //F and J equations -----------------------------------------------------------------------------------------------------------------------
        // Poisson

        /// <summary>
        /// returns the residuum of the function
        /// </summary>
        /// <param name="solution">array with all electrical potentials</param>
        /// <returns></returns>
        private Vector<double> GetFunctionPoissonSync(Vector<double> solution)
        {
            // set solution vector of previous iteration to mesh
            for (int i = 0; i < solution.Length; i++)
                mesh.finiteElements[i].phi = solution[i];

            Vector<double> F = Extreme.Mathematics.Vector.Create<double>(mesh.finiteElements.Count);
            for (int k = 0; k < mesh.finiteElements.Count; k++)
            {
                //if (mesh.points[k].material == Database.Data.materials[050000000] || mesh.points[k].material == Database.Data.materials[050100000])
                //mesh.points[k].Print();
                // microcells at the edge: F = Φ_k - Φ_init_k - U_k
                if (mesh.finiteElements[k].hasBoundaryCondition)
                    F[k] = mesh.finiteElements[k].phi - mesh.finiteElements[k].barrierHeight - mesh.finiteElements[k].phiInit;
                // other microcells
                else
                {
                    // function for non-edge-microcells
                    //     ___
                    //     \             Φ_l - Φ_k     q_e · 𝜔_k
                    // F =  ⟩  L_{k,l} · ——————————  + ——————————— · (C_k + p_k - n_k)
                    //     /             |h_{k,l}|     𝜀_r · 𝜀_0
                    //     ‾‾‾
                    //   l ∈ N(k)

                    double sum = 0;
                    foreach (var l in mesh.finiteElements[k].neighbors)
                    {

                        double x = (mesh.finiteElements[l.index].phi - mesh.finiteElements[k].phi) * l.edgeSize
                            / mesh.finiteElements[k].position.DistanceTo(mesh.finiteElements[l.index].position);
                        //  Console.WriteLine("Index: " + k + ", x: " + x);
                        sum += x;
                    }

                    double y = physConstants.e * mesh.finiteElements[k].size / (physConstants.eps0 * mesh.finiteElements[k].material.propertiesSemiconductor.epsR)
                       * (mesh.finiteElements[k].Doping(this) + mesh.finiteElements[k].pDensity(this) - mesh.finiteElements[k].nDensity(this));
                    //Console.WriteLine("p an Punkt " + k + ": " + mesh.points[k].pDensity(this));
                    //Console.WriteLine("n an Punkt " + k + ": " + mesh.points[k].nDensity(this));

                    F[k] = sum + y;
                    //Console.WriteLine("Index: " + k + ", y: " + y );
                    //Console.WriteLine("Index: " + k + ", nDensity: " + mesh.points[k].nDensity(this));
                    //Console.WriteLine("Index: " + k + ", pDensity: " + mesh.points[k].pDensity(this) + "\n");

                    //Console.WriteLine("Point " + k + ": " + F[k]);
                }
            }


            return F;
        }

        /// <summary>
        /// returns the Jacobi matrix of the function
        /// </summary>
        /// <param name="solution">array with all electrical potentials</param>
        /// <returns></returns>
        private SparseMatrix<double> GetJacobiPoissonSync(Vector<double> solution)
        {
            int amountOfNonZeroElements = 0;
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (mesh.finiteElements[i].hasBoundaryCondition)
                    amountOfNonZeroElements += 1;
                else
                    amountOfNonZeroElements += 1 + mesh.finiteElements[i].neighbors.Count;
            }

            int[] rowIndex = new int[amountOfNonZeroElements];
            int[] columnIndex = new int[amountOfNonZeroElements];
            double[] values = new double[amountOfNonZeroElements];

            // set solution vector of previous iteration to mesh
            for (int i = 0; i < solution.Length; i++)
                mesh.finiteElements[i].phi = solution[i];

            // Create J matrix
            //SparseMatrix<double> J = Matrix.CreateSparse<double>(mesh.points.Count, mesh.points.Count);

            int arrayIndex = 0;
            for (int k = 0; k < mesh.finiteElements.Count; k++) // iterate over single lines of Jacobi matrix
            {
                // for edge microcells: ∂J/∂Φ_k = 0 (diagonal) and ∂J/∂Φ_l = 0 (off-diagonal)
                if (mesh.finiteElements[k].hasBoundaryCondition)
                {
                    //J[k, k] = 1;
                    rowIndex[arrayIndex] = k;
                    columnIndex[arrayIndex] = k;
                    values[arrayIndex++] = 1;
                }

                // for non-edge-microcells
                else
                {
                    // off-diagonal
                    //             L_{k,l}
                    // ∂J/∂Φ_l =  ——————————
                    //            |h_{k,l}|

                    double sumDiagonal = 0;

                    foreach (var l in mesh.finiteElements[k].neighbors)
                    {
                        //J[k, l] = mesh.points[k].EdgeLength(mesh.points[l]) / mesh.points[k].position.DistanceTo(mesh.points[l].position);
                        //sumDiagonal += J[k, l];
                        rowIndex[arrayIndex] = k;
                        columnIndex[arrayIndex] = l.index;
                        values[arrayIndex++] = l.edgeSize / mesh.finiteElements[k].position.DistanceTo(mesh.finiteElements[l.index].position);

                        sumDiagonal += values[arrayIndex - 1];
                    }

                    // diagonal
                    //             ___
                    //             \     L_{k,l}     q_e · 𝜔_k
                    // ∂J/∂Φ_k = -  ⟩   —————————— - ——————————— · (- ∂p_k/∂Φ_k + ∂n_k/∂Φ_k)
                    //             /    |h_{k,l}|    𝜀_r · 𝜀_0     
                    //             ‾‾‾                                              \_______________________/
                    //           l ∈ N(k)                                          q_e * (p_k + n_k) / (k_B·T)

                    //J[k, k] = -sumDiagonal - BasicLib.physConstants.e * mesh.points[k].area / (BasicLib.physConstants.eps0 * mesh.points[k].material.propertiesSemiconductor.epsR)
                    //    * (mesh.points[k].DopingDerivationPhi(this) +  mesh.points[k].pDensity(this) + mesh.points[k].nDensity(this)) * BasicLib.physConstants.e / (BasicLib.physConstants.kB * T);
                    rowIndex[arrayIndex] = k;
                    columnIndex[arrayIndex] = k;
                    values[arrayIndex++] = -sumDiagonal + BasicLib.physConstants.e * mesh.finiteElements[k].size / (BasicLib.physConstants.eps0 * mesh.finiteElements[k].material.propertiesSemiconductor.epsR)
                        * (mesh.finiteElements[k].DopingDerivationPhi(this) + mesh.finiteElements[k].pDerivationPhi(this) - mesh.finiteElements[k].nDerivationPhi(this));

                }
            }

            var J_poisson = Matrix.CreateSparse<double>(mesh.finiteElements.Count, mesh.finiteElements.Count, rowIndex, columnIndex, values);
            //Console.WriteLine("J: " + J_poisson );
            return J_poisson;
        }

        //Van Roosbroeck
        /// <summary>
        /// residuum F for solving the van Roosbroeck system
        /// </summary>
        /// <param name="phiPsinPsip"></param>
        /// <returns></returns>
        public Vector<double> GetFunctionVanRoosbroeckSync(Vector<double> phiPsinPsip)
        {
            double[] values = new double[3 * mesh.nextAvailableFiniteElementIndex];

            // Writing last values of the potentials in dictionary
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                mesh.finiteElements[i].phi = phiPsinPsip[i];
                mesh.finiteElements[i].phi_n = phiPsinPsip[i + mesh.nextAvailableFiniteElementIndex];
                mesh.finiteElements[i].phi_p = phiPsinPsip[i + 2 * mesh.nextAvailableFiniteElementIndex];
            }

            //Set F-vector elements in values array
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                //Boundary condition
                if (mesh.finiteElements[i].hasBoundaryCondition)
                {
                    //Dirichlet BC
                    if (useFarrellBoundaryConditions == true)
                    {
                        values[i] = (mesh.finiteElements[i].phi - mesh.finiteElements[i].phiInit - mesh.finiteElements[i].currentBoundaryCondition);
                        values[i + mesh.nextAvailableFiniteElementIndex] = (mesh.finiteElements[i].phi_n - mesh.finiteElements[i].phi_n_init + mesh.finiteElements[i].currentBoundaryCondition);
                        values[i + 2 * mesh.nextAvailableFiniteElementIndex] = (mesh.finiteElements[i].phi_p - mesh.finiteElements[i].phi_p_init + mesh.finiteElements[i].currentBoundaryCondition);
                    }
                    //Neumann BC
                    if (useFarrellBoundaryConditions == false)
                    {
                        if (mesh.finiteElements[i].hasBoundaryCondition)
                            values[i] = mesh.finiteElements[i].phi - mesh.finiteElements[i].barrierHeight - mesh.finiteElements[i].phiInit - mesh.finiteElements[i].currentBoundaryCondition;
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        double sumN = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            sumN += (GetElectronCurrentToNeighbor(i, t.index)) * t.edgeSize;
                        }
                        values[i + mesh.nextAvailableFiniteElementIndex]
                                = sumN - BasicLib.physConstants.e * (mesh.finiteElements[i].TotalRecombinationRate(this) - mesh.finiteElements[i].TotalGenerationRate(this)) * mesh.finiteElements[i].size
                                                                + mesh.finiteElements[i].ElectronSurfaceRecombinationCurrent(this) * mesh.finiteElements[i].borderEdgeSize;
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        double sumP = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            sumP += (GetHoleCurrentToNeighbor(i, t.index)) * t.edgeSize;
                        }
                        values[i + 2 * mesh.nextAvailableFiniteElementIndex]
                                = sumP + BasicLib.physConstants.e * (mesh.finiteElements[i].TotalRecombinationRate(this) - mesh.finiteElements[i].TotalGenerationRate(this)) * mesh.finiteElements[i].size
                                                                    + mesh.finiteElements[i].HoleSurfaceRecombinationCurrent(this) * mesh.finiteElements[i].borderEdgeSize;
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    }

                }

                //Points without boundary condition
                else
                {
                    // Phi, electrical potential
                    double sum = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        sum += (mesh.finiteElements[t.index].phi - mesh.finiteElements[i].phi) * t.edgeSize
                            / mesh.finiteElements[i].position.DistanceTo(mesh.finiteElements[t.index].position);
                    }
                    values[i] = sum + BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                            * (mesh.finiteElements[i].Doping(this) + mesh.finiteElements[i].pDensity(this) - mesh.finiteElements[i].nDensity(this));
                    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                    // Phi_n, e- quasi fermi niveau 
                    double sumN = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        sumN += (GetElectronCurrentToNeighbor(i, t.index) * t.edgeSize);
                    }
                    values[i + mesh.nextAvailableFiniteElementIndex]
                            = sumN - BasicLib.physConstants.e * (mesh.finiteElements[i].TotalRecombinationRate(this) - mesh.finiteElements[i].TotalGenerationRate(this)) * mesh.finiteElements[i].size;
                    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                    //Phi_p, h+ quasi fermi niveau 
                    double sumP = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        sumP += GetHoleCurrentToNeighbor(i, t.index) * t.edgeSize;
                    }
                    values[i + 2 * mesh.nextAvailableFiniteElementIndex]
                            = sumP + BasicLib.physConstants.e * (mesh.finiteElements[i].TotalRecombinationRate(this) - mesh.finiteElements[i].TotalGenerationRate(this)) * mesh.finiteElements[i].size;
                    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                }
            }

            return Vector.Create<double>(values);
        }

        /// <summary>
        /// OLD VERSION, jacobi matrix for solving the van Roosbroeck system, old version for runtime tests
        /// </summary>
        /// <param name="phiPsinPsip"></param>
        /// <returns></returns>
        public SparseMatrix<double> JvRB(Vector<double> phiPsinPsip)
        {
            //Tolerance for abrupt heterojunctions (in eV)
            double heteroTolerance = 1e-4;

            // Writing last values of the potentials in dictionary██████████████████████████████████████████
            int k = 0; // geht nur bis Number_Meshpoints
            int[] translate = new int[mesh.nextAvailableFiniteElementIndex];
            for (int i = 0; i < mesh.finiteElements.Keys.Max() + 1; i++)
            {
                if (mesh.finiteElements.ContainsKey(i))
                {
                    mesh.finiteElements[i].phi = phiPsinPsip[k];
                    mesh.finiteElements[i].phi_n = phiPsinPsip[k + mesh.nextAvailableFiniteElementIndex];
                    mesh.finiteElements[i].phi_p = phiPsinPsip[k + 2 * mesh.nextAvailableFiniteElementIndex];
                    translate[k] = i;
                    k++;
                }
            }

            bool calculationAtZeroVolts = false;
            foreach (var v in mesh.finiteElements.Values)
            {
                if (v.hasBoundaryCondition)
                    if (v.currentBoundaryCondition > 0)
                        calculationAtZeroVolts = false;
            }

            // Creating new matrix vor van Roosbroeck problem███████████████████████████████████████████████
            SparseMatrix<double> JvRB = Matrix.CreateSparse<double>(3 * mesh.nextAvailableFiniteElementIndex, 3 * mesh.nextAvailableFiniteElementIndex);

            //  ( dF1/dPhi  dF1/dPsin  dF1/dPsip )
            //  ( dF2/dPhi  dF2/dPsin  dF2/dPsip )    = J
            //  ( dF3/dPhi  dF3/dPsin  dF3/dPsip )

            k = 0;
            // Definition of single blocks of the matrix:███████████████████████████████████████████████████
            for (int i = 0; i < mesh.finiteElements.Keys.Max() + 1; i++)
            {
                if (mesh.finiteElements.ContainsKey(i))
                {
                    //Boundary Conditions
                    if (mesh.finiteElements[i].hasBoundaryCondition)
                    {
                        #region BC Farrell /  Dirichlet 
                        if (useFarrellBoundaryConditions == true)
                        {
                            //Block 1
                            JvRB[k, k] = 1;
                            //Block 2
                            JvRB[k, k + mesh.nextAvailableFiniteElementIndex] = 0;
                            //Block 3
                            JvRB[k, k + 2 * mesh.nextAvailableFiniteElementIndex] = 0;
                            //Block 4
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k] = 0;
                            //Block 5
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex] = 1;
                            //Block 6
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex] = 0;
                            //Block 7
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k] = 0;
                            //Block 8
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex] = 0;
                            //Block 9
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex] = 1;
                        }
                        #endregion

                        #region BC Surface recombination
                        if (useFarrellBoundaryConditions == false)
                        {
                            // Block 1 Block 2 and 3 equal 0
                            JvRB[k, k] = 1;

                            // Block 2 and 3 equal 0

                            //Block 4
                            #region RB Block 4

                            double sumB4 = 0;
                            foreach (var t in mesh.finiteElements[i].neighbors)
                            {
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index] = jnDerivationPhiIndex2(i, t.index) * t.edgeSize;
                                //JvRB[k + mesh.highestPointIndex, t] = -jnDerivationPhiIndex1(t, i) * mesh.points[i].EdgeLength(mesh.points[t]);
                                sumB4 += jnDerivationPhiIndex1(i, t.index) * t.edgeSize;
                            }
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k] = sumB4 - Math.Abs(mesh.finiteElements[i].size)
                                * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                                            + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                            //Block 5
                            #region RB Block 5

                            double sumB5 = 0;
                            foreach (var t in mesh.finiteElements[i].neighbors)
                            {
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index + mesh.nextAvailableFiniteElementIndex]
                                    = jnDerivationPsinIndex2(i, t.index) * t.edgeSize;
                                sumB5 += jnDerivationPsinIndex1(i, t.index) * t.edgeSize;//JvRB[k + mesh.highestPointIndex, t + mesh.highestPointIndex];
                            }
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex]
                                = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this)
                                 + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 1) * mesh.finiteElements[i].borderEdgeSize;

                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                            //Block 6
                            #region RB Block 6
                            JvRB[k + mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex]
                                = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                            //Block 7
                            #region RB Block 7

                            double sumB7 = 0;
                            foreach (var t in mesh.finiteElements[i].neighbors)
                            {
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index] = jpDerivationPhiIndex2(i, t.index) * t.edgeSize;
                                //JvRB[k + 2 * mesh.highestPointIndex, t] = -jpDerivationPhiIndex1(t, i) * mesh.points[i].EdgeLength(mesh.points[t]);
                                //sumB7 += JvRB[k + 2 * mesh.highestPointIndex, t];
                                sumB7 += jpDerivationPhiIndex1(i, t.index) * t.edgeSize;
                            }
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k] = sumB7
                                + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                                 + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                            //Block 8
                            #region RB Block 8
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex]
                                = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                            //Block 9
                            #region RB Block 9

                            double sumB9 = 0;
                            foreach (var t in mesh.finiteElements[i].neighbors)
                            {
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index + 2 * mesh.nextAvailableFiniteElementIndex]
                                    = jpDerivationPsipIndex2(i, t.index) * t.edgeSize;
                                sumB9 += jpDerivationPsipIndex1(i, t.index) * t.edgeSize;// JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                            }
                            JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex] = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this)
                                + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 2) * mesh.finiteElements[i].borderEdgeSize;

                            #endregion
                            //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        }

                        #endregion
                    }

                    // Points without boundary conditions
                    else
                    {
                        //  (x - -) Block 1█████████████████████████████████████████████████████████████████████████████
                        //  (- - -) 
                        //  (- - -) dF1/dPhi
                        #region Block 1 Mitte
                        double sumB1 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            JvRB[k, t.index] = t.edgeSize
                                / mesh.finiteElements[i].position.DistanceTo(mesh.finiteElements[t.index].position);
                            sumB1 += JvRB[k, t.index];
                            //sumB1 += JvRB[i, t]; //TODO: richtig?
                        }

                        JvRB[k, k] = -sumB1 + BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size)
                            / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                            * (mesh.finiteElements[i].DopingDerivationPhi(this) + mesh.finiteElements[i].pDerivationPhi(this) - mesh.finiteElements[i].nDerivationPhi(this));
                        /*
                        JvRB[k, k] =  (-(1 / step_width(i, mesh.points[i].neighbors[1]) + 1 / step_width(mesh.points[i].neighbors[0], i))
                            + Basic.e * Math.Abs(mesh.points[i].area) / (mesh.points[i].material.propertiesSemiconductor.epsR * Basic.eps0) 
                            * (mesh.points[i].DopingDerivationPhi() + mesh.points[i].pDerivationPhi(this) - mesh.points[i].nDerivationPhi(this)));

                        JvRB[k, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[1])] =  (1 / step_width(i, mesh.points[i].neighbors[1]));
                        JvRB[k, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[0])] =  (1 / step_width(mesh.points[i].neighbors[0], i));
                        */
                        #endregion
                        //  (- x -) Block 2█████████████████████████████████████████████████████████████████████████████
                        //  (- - -)
                        //  (- - -)  dF1/dPsin
                        #region Block 2 Mitte
                        JvRB[k, k + mesh.nextAvailableFiniteElementIndex] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                        * (mesh.finiteElements[i].DopingDerivationPsin(this) - mesh.finiteElements[i].nDerivationPhi_n(this)));
                        #endregion
                        //  (- - x) Block 3█████████████████████████████████████████████████████████████████████████████
                        //  (- - -)
                        //  (- - -)  dF1/dPsip
                        #region Block 3 Mitte
                        JvRB[k, k + 2 * mesh.nextAvailableFiniteElementIndex] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                        * (mesh.finiteElements[i].DopingDerivationPsip(this) + mesh.finiteElements[i].pDerivationPhi_p(this)));
                        #endregion
                        //  (- - -) Block 4█████████████████████████████████████████████████████████████████████████████
                        //  (x - -)
                        //  (- - -)  dF2/dPhi
                        #region Block 4 Mitte
                        double sumB4 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {

                            double diffEc = mesh.finiteElements[i].GetDiffEcToNeighbor(mesh.finiteElements[t.index]);

                            if (Math.Abs(diffEc) > heteroTolerance && calculationAtZeroVolts == false)
                            {
                                // off diagonal element
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index]
                                    = JnTE_PhiDerivationNeighbor(i, t.index, diffEc) * t.edgeSize;

                                //diagonal element
                                sumB4 += JnTE_PhiDerivationPoint(i, t.index, diffEc) * t.edgeSize;
                            }
                            else
                            {
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index] = jnDerivationPhiIndex2(i, t.index) * t.edgeSize;
                                sumB4 += jnDerivationPhiIndex1(i, t.index) * t.edgeSize;
                            }
                        }
                        JvRB[k + mesh.nextAvailableFiniteElementIndex, k] = sumB4
                            - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);


                        #endregion
                        //  (- - -) Block 5█████████████████████████████████████████████████████████████████████████████
                        //  (- x -)
                        //  (- - -)  dF2/dPsin
                        #region Block 5 Mitte
                        double sumB5 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {

                            double diffEc = mesh.finiteElements[i].GetDiffEcToNeighbor(mesh.finiteElements[t.index]);

                            if (Math.Abs(diffEc) > heteroTolerance && calculationAtZeroVolts == false)
                            {
                                //off diagonal element
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index + mesh.nextAvailableFiniteElementIndex]
                                    = JnTE_PhinDerivationNeighbor(i, t.index, diffEc) * t.edgeSize;

                                // diagonal element
                                sumB5 += JnTE_PhinDerivationPoint(i, t.index, diffEc) * t.edgeSize;
                            }
                            else
                            {
                                JvRB[k + mesh.nextAvailableFiniteElementIndex, t.index + mesh.nextAvailableFiniteElementIndex]
                                = jnDerivationPsinIndex2(i, t.index) * t.edgeSize;
                                sumB5 += jnDerivationPsinIndex1(i, t.index) * t.edgeSize;//JvRB[k + mesh.highestPointIndex, t + mesh.highestPointIndex];
                            }
                        }
                        JvRB[k + mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex]
                            = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this);



                        #endregion
                        //  (- - -) Block 6█████████████████████████████████████████████████████████████████████████████
                        //  (- - x)
                        //  (- - -)  dF2/dPsip
                        #region Block 6 Mitte
                        JvRB[k + mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex]
                            = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                        #endregion
                        //  (- - -) Block 7█████████████████████████████████████████████████████████████████████████████
                        //  (- - -)
                        //  (x - -)  dF3/dPhi
                        #region Block 7 Mitte
                        double sumB7 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            double diffEv = mesh.finiteElements[i].GetDiffEvToNeighbor(mesh.finiteElements[t.index]);

                            if (Math.Abs(diffEv) > heteroTolerance && calculationAtZeroVolts == false)
                            {
                                //off diagonal element
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index]
                                    = JpTE_PhiDerivationNeighbor(i, t.index, diffEv) * t.edgeSize;

                                //diagonal element
                                sumB7 += JpTE_PhiDerivationPoint(i, t.index, diffEv) * t.edgeSize;
                            }
                            else
                            {
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index] = jpDerivationPhiIndex2(i, t.index) * t.edgeSize;
                            }
                            //JvRB[k + 2*mesh.highestPointIndex, t] = -jpDerivationPhiIndex1(t, i) * mesh.points[i].EdgeLength(mesh.points[t]);
                            sumB7 += jpDerivationPhiIndex1(i, t.index) * t.edgeSize;
                        }
                        JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k]
                            = sumB7 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);

                        #endregion
                        //  (- - -) Block 8█████████████████████████████████████████████████████████████████████████████
                        //  (- - -)
                        //  (- x -)  dF3/dPsin
                        #region Block 8 Mitte
                        JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + mesh.nextAvailableFiniteElementIndex]
                            = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                        #endregion
                        //  (- - -) Block 9█████████████████████████████████████████████████████████████████████████████
                        //  (- - -)
                        //  (- - x)  dF3/dPsip
                        #region Block 9 Mitte
                        double sumB9 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {

                            double diffEv = mesh.finiteElements[i].GetDiffEvToNeighbor(mesh.finiteElements[t.index]);

                            if (Math.Abs(diffEv) > heteroTolerance && calculationAtZeroVolts == false)
                            {

                                //Console.WriteLine("Punkt " + i + " , Nachbarn: " + t + " ,Edgelength: " + mesh.points[i].EdgeLength(mesh.points[t]));

                                //off diagonal element
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index + 2 * mesh.nextAvailableFiniteElementIndex]
                                    = JpTE_PhipDerivationNeighbor(i, t.index, diffEv) * t.edgeSize;

                                //diagonal element
                                sumB9 += JpTE_PhipDerivationPoint(i, t.index, diffEv) * t.edgeSize;
                            }
                            else
                            {
                                JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, t.index + 2 * mesh.nextAvailableFiniteElementIndex]
                                = jpDerivationPsipIndex2(i, t.index) * t.edgeSize;
                                sumB9 += jpDerivationPsipIndex1(i, t.index) * t.edgeSize;// JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                                                                                         //sumB9 += JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                            }
                        }
                        JvRB[k + 2 * mesh.nextAvailableFiniteElementIndex, k + 2 * mesh.nextAvailableFiniteElementIndex]
                            = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this);

                        //JvRB[k + 2* mesh.highestPointIndex, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[1]) + 2 * mesh.highestPointIndex] =  (jpDerivationPsipIndex2(i, mesh.points[i].neighbors[1]));

                        #endregion
                    }

                    k++;
                }
            }
            /*
            //Block (1)
            Console.WriteLine("J_(1) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));

            // Block (2) und (3)
            Console.WriteLine("J_n_(2) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(3) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));

            // Block (4) und (7)
            Console.WriteLine("J_n_(4) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(7) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));

            // Block (5) und (9)
            Console.WriteLine("J_n_(5) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(9) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));

            // Block (6) und (8)
            Console.WriteLine("J_n_(6) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(8) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            */

            //Console.WriteLine("J Determinante: " +  JvRB.GetDeterminant());
            /*
            for (int i = 0; i < JvRB.ColumnCount; i++)
            {
                for (int j = 0; j < JvRB.ColumnCount; j++)
                    JvRB[i, j] *= 1;// 1e100;
            }
            */
            return JvRB;
        }

        /// <summary>
        /// OLD VERSION without parallelization, jacobi matrix for solving the van Roosbroeck system
        /// </summary>
        /// <param name="phiPsinPsip"></param>
        /// <returns></returns>
        public SparseMatrix<double> GetJacobiVanRoosbroeckSync(Vector<double> phiPsinPsip)
        {
            double heteroTolerance = 1e-4;
            //Create arrays for sparse marix 
            int amountOfNonZeroElements = 0;
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (mesh.finiteElements[i].hasBoundaryCondition)
                {
                    if (useFarrellBoundaryConditions == true)
                        amountOfNonZeroElements += 9;
                    else
                        amountOfNonZeroElements += 7 + 4 * mesh.finiteElements[i].neighbors.Count;
                }
                else
                    amountOfNonZeroElements += 9 + 5 * mesh.finiteElements[i].neighbors.Count;
            }
            int[] rowIndex = new int[amountOfNonZeroElements];
            int[] columnIndex = new int[amountOfNonZeroElements];
            double[] values = new double[amountOfNonZeroElements];


            // Writing last values of the potentials in dictionary
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                mesh.finiteElements[i].phi = phiPsinPsip[i];
                mesh.finiteElements[i].phi_n = phiPsinPsip[i + mesh.nextAvailableFiniteElementIndex];
                mesh.finiteElements[i].phi_p = phiPsinPsip[i + 2 * mesh.nextAvailableFiniteElementIndex];
            }

            //  ( dF1/dPhi  dF1/dPhi_n  dF1/dPhi_p )
            //  ( dF2/dPhi  dF2/dPhi_n  dF2/dPhi_p )    = J
            //  ( dF3/dPhi  dF3/dPhi_n  dF3/dPdPhi_psip )


            // Definition of single blocks of the matrix:
            int arrayIndex = 0;
            for (int i = 0; i < mesh.finiteElements.Keys.Max() + 1; i++)
            {
                //Boundary Conditions
                if (mesh.finiteElements[i].hasBoundaryCondition)
                {
                    if (useFarrellBoundaryConditions == true)
                    {
                        //Block 1
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 1;
                        //Block 2
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 3
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 4
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 0;
                        //Block 5
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 1;
                        //Block 6
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 7
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 0;
                        //Block 8
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 9
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 1;
                    }

                    if (useFarrellBoundaryConditions == false)
                    {
                        // Block 1 Block 2 and 3 equal 0
                        //JvRB[i, i] = 1;
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 1;

                        // Block 2 and 3 equal 0

                        //Block 4
                        #region RB Block 4

                        double sumB4 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            //JvRB[i + mesh.highestPointIndex, t] = jnDerivationPhiIndex2(i, t) * mesh.points[i].EdgeLength(mesh.points[t]);
                            sumB4 += jnDerivationPhiIndex1(i, t.index) * t.edgeSize;
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = jnDerivationPhiIndex2(i, t.index) * t.edgeSize;
                        }
                        //JvRB[i + mesh.highestPointIndex, i] 
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = sumB4 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                                        + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 5
                        #region RB Block 5

                        double sumB5 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            //JvRB[i + mesh.highestPointIndex, t + mesh.highestPointIndex]
                            //  = jnDerivationPsinIndex2(i, t) * mesh.points[i].EdgeLength(mesh.points[t]);
                            sumB5 += jnDerivationPsinIndex1(i, t.index) * t.edgeSize;//JvRB[k + mesh.highestPointIndex, t + mesh.highestPointIndex];
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = jnDerivationPsinIndex2(i, t.index) * t.edgeSize;
                        }
                        //JvRB[i + mesh.highestPointIndex, i + mesh.highestPointIndex]
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this)
                            + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 1) * mesh.finiteElements[i].borderEdgeSize;

                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 6
                        #region RB Block 6
                        //JvRB[i + mesh.highestPointIndex, i + 2 * mesh.highestPointIndex]
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 7
                        #region RB Block 7

                        double sumB7 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            //JvRB[i + 2 * mesh.highestPointIndex, t] 
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = jpDerivationPhiIndex2(i, t.index) * t.edgeSize;
                            sumB7 += jpDerivationPhiIndex1(i, t.index) * t.edgeSize;
                        }
                        //JvRB[i + 2 * mesh.highestPointIndex, i] 
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = sumB7 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                            + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 8
                        #region RB Block 8
                        //JvRB[i + 2 * mesh.highestPointIndex, i + mesh.highestPointIndex]
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 9
                        #region RB Block 9

                        double sumB9 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            //JvRB[i + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex]
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + 2 * mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = jpDerivationPsipIndex2(i, t.index) * t.edgeSize;
                            sumB9 += jpDerivationPsipIndex1(i, t.index) * t.edgeSize;// JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                        }
                        //JvRB[i + 2 * mesh.highestPointIndex, i + 2 * mesh.highestPointIndex] = 
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this)
                            + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 2) * mesh.finiteElements[i].borderEdgeSize;

                        #endregion
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    }

                }

                // Points without boundary conditions
                else
                {
                    //  (x - -) Block 1█████████████████████████████████████████████████████████████████████████████
                    //  (- - -) 
                    //  (- - -) dF1/dPhi
                    #region Block 1 Mitte
                    double sumB1 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        //JvRB[i, t] 
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = t.index;
                        values[arrayIndex++] = t.edgeSize
                            / mesh.finiteElements[i].position.DistanceTo(mesh.finiteElements[t.index].position);
                        sumB1 += values[arrayIndex - 1];
                    }

                    //JvRB[i, i]
                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = -sumB1 + BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size)
                        / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                        * (mesh.finiteElements[i].DopingDerivationPhi(this) + mesh.finiteElements[i].pDerivationPhi(this) - mesh.finiteElements[i].nDerivationPhi(this));
                    /*
                    JvRB[k, k] =  (-(1 / step_width(i, mesh.points[i].neighbors[1]) + 1 / step_width(mesh.points[i].neighbors[0], i))
                        + Basic.e * Math.Abs(mesh.points[i].area) / (mesh.points[i].material.propertiesSemiconductor.epsR * Basic.eps0) 
                        * (mesh.points[i].DopingDerivationPhi() + mesh.points[i].pDerivationPhi(this) - mesh.points[i].nDerivationPhi(this)));

                    JvRB[k, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[1])] =  (1 / step_width(i, mesh.points[i].neighbors[1]));
                    JvRB[k, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[0])] =  (1 / step_width(mesh.points[i].neighbors[0], i));
                    */
                    #endregion
                    //  (- x -) Block 2█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - -)  dF1/dPsin
                    #region Block 2 Mitte
                    // JvRB[i, i + mesh.highestPointIndex] 
                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                    * (mesh.finiteElements[i].DopingDerivationPsin(this) - mesh.finiteElements[i].nDerivationPhi_n(this)));
                    #endregion
                    //  (- - x) Block 3█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - -)  dF1/dPsip
                    #region Block 3 Mitte
                    //JvRB[i, i + 2 * mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                    * (mesh.finiteElements[i].DopingDerivationPsip(this) + mesh.finiteElements[i].pDerivationPhi_p(this)));
                    #endregion
                    //  (- - -) Block 4█████████████████████████████████████████████████████████████████████████████
                    //  (x - -)
                    //  (- - -)  dF2/dPhi
                    #region Block 4 Mitte
                    double sumB4 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {

                        double diffEc = mesh.finiteElements[i].GetDiffEcToNeighbor(mesh.finiteElements[t.index]);

                        if (Math.Abs(diffEc) > heteroTolerance)
                        {
                            // off diagonal element
                            //JvRB[i + mesh.highestPointIndex, t]
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = JnTE_PhiDerivationNeighbor(i, t.index, diffEc) * t.edgeSize;

                            //diagonal element
                            sumB4 += JnTE_PhiDerivationPoint(i, t.index, diffEc) * t.edgeSize;
                        }
                        else
                        {
                            //JvRB[i + mesh.highestPointIndex, t] 
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = jnDerivationPhiIndex2(i, t.index) * t.edgeSize;
                            sumB4 += jnDerivationPhiIndex1(i, t.index) * t.edgeSize;
                        }
                    }
                    //JvRB[i + mesh.highestPointIndex, i] =
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = sumB4 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);


                    #endregion
                    //  (- - -) Block 5█████████████████████████████████████████████████████████████████████████████
                    //  (- x -)
                    //  (- - -)  dF2/dPsin
                    #region Block 5 Mitte
                    double sumB5 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {

                        double diffEc = mesh.finiteElements[i].GetDiffEcToNeighbor(mesh.finiteElements[t.index]);

                        if (Math.Abs(diffEc) > heteroTolerance)
                        {
                            //off diagonal element
                            //JvRB[i + mesh.highestPointIndex, t + mesh.highestPointIndex]
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = JnTE_PhinDerivationNeighbor(i, t.index, diffEc) * t.edgeSize;

                            // diagonal element
                            sumB5 += JnTE_PhinDerivationPoint(i, t.index, diffEc) * t.edgeSize;
                        }
                        else
                        {
                            //JvRB[i + mesh.highestPointIndex, t + mesh.highestPointIndex]
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = jnDerivationPsinIndex2(i, t.index) * t.edgeSize;
                            sumB5 += jnDerivationPsinIndex1(i, t.index) * t.edgeSize;//JvRB[k + mesh.highestPointIndex, t + mesh.highestPointIndex];
                        }
                    }
                    //JvRB[i + mesh.highestPointIndex, i + mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this);



                    #endregion
                    //  (- - -) Block 6█████████████████████████████████████████████████████████████████████████████
                    //  (- - x)
                    //  (- - -)  dF2/dPsip
                    #region Block 6 Mitte
                    //JvRB[i + mesh.highestPointIndex, i + 2 * mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                    #endregion
                    //  (- - -) Block 7█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (x - -)  dF3/dPhi
                    #region Block 7 Mitte
                    double sumB7 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        double diffEv = mesh.finiteElements[i].GetDiffEvToNeighbor(mesh.finiteElements[t.index]);

                        if (Math.Abs(diffEv) > heteroTolerance)
                        {
                            //off diagonal element
                            //JvRB[i + 2 * mesh.highestPointIndex, t]
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = JpTE_PhiDerivationNeighbor(i, t.index, diffEv) * t.edgeSize;

                            //diagonal element
                            sumB7 += JpTE_PhiDerivationPoint(i, t.index, diffEv) * t.edgeSize;
                        }
                        else
                        {
                            //JvRB[i + 2 * mesh.highestPointIndex, t] 
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = jpDerivationPhiIndex2(i, t.index) * t.edgeSize;
                        }
                        //JvRB[k + 2*mesh.highestPointIndex, t] = -jpDerivationPhiIndex1(t, i) * mesh.points[i].EdgeLength(mesh.points[t]);
                        sumB7 += jpDerivationPhiIndex1(i, t.index) * t.edgeSize;
                    }
                    //JvRB[i + 2 * mesh.highestPointIndex, i]
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = sumB7 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);

                    #endregion
                    //  (- - -) Block 8█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- x -)  dF3/dPsin
                    #region Block 8 Mitte
                    //JvRB[i + 2 * mesh.highestPointIndex, i + mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                    #endregion
                    //  (- - -) Block 9█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - x)  dF3/dPsip
                    #region Block 9 Mitte
                    double sumB9 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {

                        double diffEv = mesh.finiteElements[i].GetDiffEvToNeighbor(mesh.finiteElements[t.index]);

                        if (Math.Abs(diffEv) > heteroTolerance)
                        {

                            //Console.WriteLine("Punkt " + i + " , Nachbarn: " + t + " ,Edgelength: " + mesh.points[i].EdgeLength(mesh.points[t]));

                            //off diagonal element
                            //JvRB[i + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex]
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + 2 * mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = JpTE_PhipDerivationNeighbor(i, t.index, diffEv) * t.edgeSize;

                            //diagonal element
                            sumB9 += JpTE_PhipDerivationPoint(i, t.index, diffEv) * t.edgeSize;
                        }
                        else
                        {
                            // JvRB[i + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex]
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + 2 * mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = jpDerivationPsipIndex2(i, t.index) * t.edgeSize;
                            sumB9 += jpDerivationPsipIndex1(i, t.index) * t.edgeSize;// JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                                                                                     //sumB9 += JvRB[k + 2 * mesh.highestPointIndex, t + 2 * mesh.highestPointIndex];
                        }
                    }
                    //JvRB[i + 2 * mesh.highestPointIndex, i + 2 * mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this);

                    //JvRB[k + 2* mesh.highestPointIndex, Array.FindIndex(translate, row => row == mesh.points[i].neighbors[1]) + 2 * mesh.highestPointIndex] =  (jpDerivationPsipIndex2(i, mesh.points[i].neighbors[1]));

                    #endregion
                }

            }
            /*
            //Block (1)
            Console.WriteLine("J_(1) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));

            // Block (2) und (3)
            Console.WriteLine("J_n_(2) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(3) = " + JvRB.GetSubmatrix(0, mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));

            // Block (4) und (7)
            Console.WriteLine("J_n_(4) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(7) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, 0, mesh.highestPointIndex - 1));

            // Block (5) und (9)
            Console.WriteLine("J_n_(5) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(9) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));

            // Block (6) und (8)
            Console.WriteLine("J_n_(6) = " + JvRB.GetSubmatrix(mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1, 2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1));
            Console.WriteLine("J_p_(8) = " + JvRB.GetSubmatrix(2 * mesh.highestPointIndex, 3 * mesh.highestPointIndex - 1, mesh.highestPointIndex, 2 * mesh.highestPointIndex - 1));
            */

            //Console.WriteLine("J Determinante: " +  JvRB.GetDeterminant());
            int matrixSize = 3 * mesh.nextAvailableFiniteElementIndex;
            SparseMatrix<double> J = Matrix.CreateSparse(matrixSize, matrixSize, rowIndex, columnIndex, values);
            return J;
        }

        /// <summary>
        /// returns the Jacobi matrix of the function
        /// </summary>
        /// <param name="solution">array with all electrical potentials</param>
        /// <returns></returns>
        private SparseMatrix<double> GetJacobiVanRoosbroeckAsync(Vector<double> phiPsinPsip)
        {
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                mesh.finiteElements[i].phi = phiPsinPsip[i];
                mesh.finiteElements[i].phi_n = phiPsinPsip[i + mesh.nextAvailableFiniteElementIndex];
                mesh.finiteElements[i].phi_p = phiPsinPsip[i + 2 * mesh.nextAvailableFiniteElementIndex];

            }

            // Get the amount of logical processor units and create a list of as many tasks as processors
            int amountProcessors = Environment.ProcessorCount;
            Task<(int[] rowIndexes, int[] columnIndexes, double[] values)>[] tasks
                = new Task<(int[] rowIndexes, int[] columnIndexes, double[] values)>[amountProcessors];

            // Get the amount of rows, one single processor will calculate
            // (always rounded down and the last processor will have to do the modulo-rest on top)
            int lengthOfBlock = mesh.nextAvailableFiniteElementIndex / amountProcessors;

            // Create and start all tasks
            for (int p = 0; p < amountProcessors - 1; p++)
            {
                int processorNumber = p;
                tasks[p] = Task.Run(() => GetJacobiVanRoosbroeckBlock(processorNumber * lengthOfBlock, (processorNumber + 1) * lengthOfBlock));
            }
            tasks[amountProcessors - 1] = Task.Run(() => GetJacobiVanRoosbroeckBlock((amountProcessors - 1) * lengthOfBlock, mesh.nextAvailableFiniteElementIndex));

            // Wait for all tasks to be finished => all blocks are set
            var resultListArray = Task.WhenAll(tasks).Result;

            // Create Arrays of all single arrays
            int[][] rowIndexesSingle = new int[amountProcessors][];
            int[][] columnIndexesSingle = new int[amountProcessors][];
            double[][] valuesSingle = new double[amountProcessors][];
            for (int i = 0; i < amountProcessors; i++)
            {
                rowIndexesSingle[i] = resultListArray[i].rowIndexes;
                columnIndexesSingle[i] = resultListArray[i].columnIndexes;
                valuesSingle[i] = resultListArray[i].values;
            }

            // Create and return Jacobi sparse matrix from arrays
            int matrixSize = 3 * mesh.nextAvailableFiniteElementIndex;

            SparseMatrix<double> J = Matrix.CreateSparse(matrixSize, matrixSize,
                Misc.ConcatArrays(rowIndexesSingle), Misc.ConcatArrays(columnIndexesSingle), Misc.ConcatArrays(valuesSingle));
            return J;
        }

        /// <summary>
        /// jacobi matrix called by GetJacobiVanRoosbroeckAsync
        /// </summary>
        /// <param name="phiPsinPsip"></param>
        /// <returns></returns>
        private (int[] rowIndexes, int[] columnIndexes, double[] values) GetJacobiVanRoosbroeckBlock(int startPointIndexIncluding, int endPointIndexExcluding)
        {
            int amountOfNonZeroElements = 0;
            for (int i = startPointIndexIncluding; i < endPointIndexExcluding; i++)
            {
                if (mesh.finiteElements[i].hasBoundaryCondition)
                {
                    if (useFarrellBoundaryConditions == true)
                        amountOfNonZeroElements += 9;
                    else
                        amountOfNonZeroElements += 7 + 4 * mesh.finiteElements[i].neighbors.Count;
                }
                else
                    amountOfNonZeroElements += 9 + 5 * mesh.finiteElements[i].neighbors.Count;
            }


            int[] rowIndex = new int[amountOfNonZeroElements];
            int[] columnIndex = new int[amountOfNonZeroElements];
            double[] values = new double[amountOfNonZeroElements];

            int arrayIndex = 0;
            // Definition of single blocks of the matrix:███████████████████████████████████████████████████
            for (int i = startPointIndexIncluding; i < endPointIndexExcluding; i++)
            {
                //Boundary Conditions
                if (mesh.finiteElements[i].hasBoundaryCondition)
                {
                    #region BC Farrell /  Dirichlet 
                    if (useFarrellBoundaryConditions == true)
                    {
                        //Block 1
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 1;
                        //Block 2
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 3
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 4
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 0;
                        //Block 5
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 1;
                        //Block 6
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 7
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 0;
                        //Block 8
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 0;
                        //Block 9
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = 1;
                    }
                    #endregion

                    #region BC Surface recombination
                    if (useFarrellBoundaryConditions == false)
                    {
                        // Block 1 
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = 1;

                        // Block 2 and 3 equal 0

                        //Block 4
                        double sumB4 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.neighbor) * t.edgeSize;
                            sumB4 += GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.point) * t.edgeSize;
                        }
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = sumB4 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                                        + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 5
                        double sumB5 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_n, DerivateAtIndex.neighbor) * t.edgeSize;
                            sumB5 += GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_n, DerivateAtIndex.point) * t.edgeSize;
                        }
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this)
                            + mesh.finiteElements[i].surfaceRecElectronsDerivation(this, i, 1) * mesh.finiteElements[i].borderEdgeSize;

                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 6
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 7
                        double sumB7 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index;
                            values[arrayIndex++] = GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.neighbor) * t.edgeSize;
                            sumB7 += GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.point) * t.edgeSize;
                        }
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i;
                        values[arrayIndex++] = sumB7 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this)
                            + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 0) * mesh.finiteElements[i].borderEdgeSize;

                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 8
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                        //Block 9
                        double sumB9 = 0;
                        foreach (var t in mesh.finiteElements[i].neighbors)
                        {
                            rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                            columnIndex[arrayIndex] = t.index + 2 * mesh.nextAvailableFiniteElementIndex;
                            values[arrayIndex++] = GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_p, DerivateAtIndex.neighbor) * t.edgeSize;
                            sumB9 += GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_p, DerivateAtIndex.point) * t.edgeSize;
                        }
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this)
                            + mesh.finiteElements[i].surfaceRecHolesDerivation(this, i, 2) * mesh.finiteElements[i].borderEdgeSize;

                        //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    }

                    #endregion
                }

                // Points without boundary conditions
                else
                {
                    //  (x - -) Block 1█████████████████████████████████████████████████████████████████████████████
                    //  (- - -) 
                    //  (- - -) dF1/dPhi
                    #region Block 1 Mitte
                    double sumB1 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        rowIndex[arrayIndex] = i;
                        columnIndex[arrayIndex] = t.index;
                        values[arrayIndex++] = t.edgeSize / mesh.finiteElements[i].position.DistanceTo(mesh.finiteElements[t.index].position);
                        sumB1 += values[arrayIndex - 1];
                    }

                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = -sumB1 + BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                        * (mesh.finiteElements[i].DopingDerivationPhi(this) + mesh.finiteElements[i].pDerivationPhi(this) - mesh.finiteElements[i].nDerivationPhi(this));

                    #endregion
                    //  (- x -) Block 2█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - -)  dF1/dPsin
                    #region Block 2 Mitte
                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                    * (mesh.finiteElements[i].DopingDerivationPsin(this) - mesh.finiteElements[i].nDerivationPhi_n(this)));
                    #endregion
                    //  (- - x) Block 3█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - -)  dF1/dPsip
                    #region Block 3 Mitte
                    //JvRB[i, i + 2 * mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (BasicLib.physConstants.e * Math.Abs(mesh.finiteElements[i].size) / (mesh.finiteElements[i].material.propertiesSemiconductor.epsR * BasicLib.physConstants.eps0)
                                                    * (mesh.finiteElements[i].DopingDerivationPsip(this) + mesh.finiteElements[i].pDerivationPhi_p(this)));
                    #endregion

                    //  (- - -) Block 4█████████████████████████████████████████████████████████████████████████████
                    //  (x - -)
                    //  (- - -)  dF2/dPhi
                    #region Block 4 Mitte
                    double sumB4 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = t.index;
                        values[arrayIndex++] = GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.neighbor) * t.edgeSize;
                        sumB4 += GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.point) * t.edgeSize;
                    }
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = sumB4 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);
                    #endregion

                    //  (- - -) Block 5█████████████████████████████████████████████████████████████████████████████
                    //  (- x -)
                    //  (- - -)  dF2/dPsin
                    #region Block 5 Mitte
                    double sumB5 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = t.index + mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_n, DerivateAtIndex.neighbor) * t.edgeSize;
                        sumB5 += GetElectronCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_n, DerivateAtIndex.point) * t.edgeSize;
                    }
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = sumB5 - Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this);
                    #endregion

                    //  (- - -) Block 6█████████████████████████████████████████████████████████████████████████████
                    //  (- - x)
                    //  (- - -)  dF2/dPsip
                    #region Block 6 Mitte
                    //JvRB[i + mesh.highestPointIndex, i + 2 * mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (-Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this));
                    #endregion

                    //  (- - -) Block 7█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (x - -)  dF3/dPhi
                    #region Block 7 Mitte
                    double sumB7 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = t.index;
                        values[arrayIndex++] = GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.neighbor) * t.edgeSize;
                        sumB7 += GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi, DerivateAtIndex.point) * t.edgeSize;
                    }
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i;
                    values[arrayIndex++] = sumB7 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Phi(this);
                    #endregion

                    //  (- - -) Block 8█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- x -)  dF3/dPsin
                    #region Block 8 Mitte
                    //JvRB[i + 2 * mesh.highestPointIndex, i + mesh.highestPointIndex]
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = (Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_n(this));
                    #endregion

                    //  (- - -) Block 9█████████████████████████████████████████████████████████████████████████████
                    //  (- - -)
                    //  (- - x)  dF3/dPsip
                    #region Block 9 Mitte
                    double sumB9 = 0;
                    foreach (var t in mesh.finiteElements[i].neighbors)
                    {
                        rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                        columnIndex[arrayIndex] = t.index + 2 * mesh.nextAvailableFiniteElementIndex;
                        values[arrayIndex++] = GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_p, DerivateAtIndex.neighbor) * t.edgeSize;
                        sumB9 += GetHoleCurrentDerivationToNeighbor(i, t.index, DerivationVariable.Phi_p, DerivateAtIndex.point) * t.edgeSize;
                    }
                    rowIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    columnIndex[arrayIndex] = i + 2 * mesh.nextAvailableFiniteElementIndex;
                    values[arrayIndex++] = sumB9 + Math.Abs(mesh.finiteElements[i].size) * BasicLib.physConstants.e * mesh.finiteElements[i].Diff_Rekomb_Psi_p(this);
                    #endregion
                }

            }


            return (rowIndex, columnIndex, values);
        }


        //Current Equations------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// returns the electron current between points with indeices "pointIndex" and "indexNeighborPoint". SG or TE currents is choosen automatically.
        /// </summary>
        /// <param name="pointIndex"></param>
        /// <param name="indexNeighborPoint"></param>
        /// <returns></returns>
        public double GetElectronCurrentToNeighbor(int pointIndex, int indexNeighborPoint)
        {
            int listIndexNeighbor = mesh.finiteElements[pointIndex].GetIndexInNeighborList(indexNeighborPoint);
            var currentType = mesh.finiteElements[pointIndex].typeOfElectronCurrents[listIndexNeighbor];

            if (currentType == TypeOfCurrent.ScharfetterGummel)
                return jElectron(pointIndex, indexNeighborPoint);
            if (currentType == TypeOfCurrent.ThermionicEmission)
            {
                double diffEc = mesh.finiteElements[pointIndex].GetDiffEcToNeighbor(mesh.finiteElements[indexNeighborPoint]);
                return ThermionicEmissionCurrentElectron(pointIndex, indexNeighborPoint, diffEc);
            }
            else
                throw new Exception("Missing Type of electron current at Point " + pointIndex);
        }

        /// <summary>
        /// returns the derivation of the electron current depending on type of current (SG or TE) 
        /// </summary>
        /// <param name="pointIndex">index of the point itself</param>
        /// <param name="indexNeighborPoint">index of the neighbor</param>
        /// <param name="derivationWithRespectTo">integer wether derivation is with respect to Phi (0), Phi_n (1) or Phi_p (2)</param>
        /// <param name="indexSelector">integer which selects the equation belonging to index1 or index2</param>
        /// <returns></returns>
        public double GetElectronCurrentDerivationToNeighbor(int pointIndex, int indexNeighborPoint, DerivationVariable derivationVariable, DerivateAtIndex derivateAtIndex)
        {
            int listIndexNeighbor = mesh.finiteElements[pointIndex].GetIndexInNeighborList(indexNeighborPoint);
            var currentType = mesh.finiteElements[pointIndex].typeOfElectronCurrents[listIndexNeighbor];

            switch (derivateAtIndex)
            {

                //point eqaution
                case DerivateAtIndex.point:
                    {
                        if (currentType == TypeOfCurrent.ScharfetterGummel)
                        {
                            if (derivationVariable == DerivationVariable.Phi)
                                return jnDerivationPhiIndex1(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return jnDerivationPsinIndex1(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return 0;
                            else
                                throw new Exception("Error: TE electron current at Point " + pointIndex);
                        }
                        else if (currentType == TypeOfCurrent.ThermionicEmission)
                        {
                            double diffEc = mesh.finiteElements[pointIndex].GetDiffEcToNeighbor(mesh.finiteElements[indexNeighborPoint]);

                            if (derivationVariable == DerivationVariable.Phi)
                                return JnTE_PhiDerivationPoint(pointIndex, indexNeighborPoint, diffEc);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return JnTE_PhinDerivationPoint(pointIndex, indexNeighborPoint, diffEc);
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return 0;
                            else
                                throw new Exception("Error: TE electron current at Point " + pointIndex);
                        }
                        else
                            throw new Exception("Missing Type of electron current at Point " + pointIndex);
                    }
                //Neighbor equation
                case DerivateAtIndex.neighbor:
                    {
                        if (currentType == TypeOfCurrent.ScharfetterGummel)
                        {
                            if (derivationVariable == DerivationVariable.Phi)
                                return jnDerivationPhiIndex2(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return jnDerivationPsinIndex2(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return 0;
                            else
                                throw new Exception("Error: TE electron current at Point " + pointIndex);
                        }
                        else if (currentType == TypeOfCurrent.ThermionicEmission)
                        {
                            double diffEc = mesh.finiteElements[pointIndex].GetDiffEcToNeighbor(mesh.finiteElements[indexNeighborPoint]);

                            if (derivationVariable == DerivationVariable.Phi)
                                return JnTE_PhiDerivationNeighbor(pointIndex, indexNeighborPoint, diffEc);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return JnTE_PhinDerivationNeighbor(pointIndex, indexNeighborPoint, diffEc);
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return 0;
                            else
                                throw new Exception("Error: TE electron current at Point " + pointIndex);

                        }
                        else
                            throw new Exception("Missing Type of electron current at Point " + pointIndex);
                    }
                default:
                    throw new Exception("Error: TE electron current at Point " + pointIndex);
            }
        }

        /// <summary>
        /// returns the hole current between points with indeices "pointIndex" and "indexNeighborPoint". SG or TE currents is choosen automatically.
        /// </summary>
        /// <param name="pointIndex"></param>
        /// <param name="indexNeighborPoint"></param>
        /// <returns></returns>
        public double GetHoleCurrentToNeighbor(int pointIndex, int indexNeighborPoint)
        {

            int listIndexNeighbor = mesh.finiteElements[pointIndex].GetIndexInNeighborList(indexNeighborPoint);
            var currentType = mesh.finiteElements[pointIndex].typeOfHoleCurrents[listIndexNeighbor];

            if (currentType == TypeOfCurrent.ScharfetterGummel)
                return jHole(pointIndex, indexNeighborPoint);
            if (currentType == TypeOfCurrent.ThermionicEmission)
            {
                double diffEv = mesh.finiteElements[pointIndex].GetDiffEvToNeighbor(mesh.finiteElements[indexNeighborPoint]);
                return ThermionicEmissionCurrentHole(pointIndex, indexNeighborPoint, diffEv);
            }
            else
                throw new Exception("Missing Type of hole current at Point " + pointIndex);
        }
        /// <summary>
        /// returns the derivation of the holr current depending on type of current (SG or TE) 
        /// </summary>
        /// <param name="pointIndex">index of the point itself</param>
        /// <param name="indexNeighborPoint">index of the neighbor</param>
        /// <param name="derivationWithRespectTo">integer wether derivation is with respect to Phi (0), Phi_n (1) or Phi_p (2)</param>
        /// <param name="indexSelector">integer which selects the equation belonging to index1 or index2</param>
        /// <returns></returns>
        public double GetHoleCurrentDerivationToNeighbor(int pointIndex, int indexNeighborPoint, DerivationVariable derivationVariable, DerivateAtIndex derivateAtIndex)
        {
            int listIndexNeighbor = mesh.finiteElements[pointIndex].GetIndexInNeighborList(indexNeighborPoint);
            var currentType = mesh.finiteElements[pointIndex].typeOfHoleCurrents[listIndexNeighbor];

            switch (derivateAtIndex)
            {
                //point eqaution
                case DerivateAtIndex.point:
                    {
                        if (currentType == TypeOfCurrent.ScharfetterGummel)
                        {
                            if (derivationVariable == DerivationVariable.Phi)
                                return jpDerivationPhiIndex1(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return 0;
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return jpDerivationPsipIndex1(pointIndex, indexNeighborPoint);
                            else
                                throw new Exception("Error: TE hole current at Point " + pointIndex);
                        }
                        else if (currentType == TypeOfCurrent.ThermionicEmission)
                        {
                            double diffEv = mesh.finiteElements[pointIndex].GetDiffEvToNeighbor(mesh.finiteElements[indexNeighborPoint]);

                            if (derivationVariable == DerivationVariable.Phi)
                                return JpTE_PhiDerivationPoint(pointIndex, indexNeighborPoint, diffEv);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return 0;
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return JpTE_PhipDerivationPoint(pointIndex, indexNeighborPoint, diffEv);
                            else
                                throw new Exception("Error: TE hole current at Point " + pointIndex);
                        }
                        else
                            throw new Exception("Missing Type of hole current at Point " + pointIndex);
                    }
                //Neighbor equation
                case DerivateAtIndex.neighbor:
                    {
                        if (currentType == TypeOfCurrent.ScharfetterGummel)
                        {
                            if (derivationVariable == DerivationVariable.Phi)
                                return jpDerivationPhiIndex2(pointIndex, indexNeighborPoint);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return 0;
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return jpDerivationPsipIndex2(pointIndex, indexNeighborPoint);
                            else
                                throw new Exception("Error: TE hole current at Point " + pointIndex);
                        }
                        else if (currentType == TypeOfCurrent.ThermionicEmission)
                        {
                            double diffEv = mesh.finiteElements[pointIndex].GetDiffEvToNeighbor(mesh.finiteElements[indexNeighborPoint]);

                            if (derivationVariable == DerivationVariable.Phi)
                                return JpTE_PhiDerivationNeighbor(pointIndex, indexNeighborPoint, diffEv);
                            else if (derivationVariable == DerivationVariable.Phi_n)
                                return 0;
                            else if (derivationVariable == DerivationVariable.Phi_p)
                                return JpTE_PhipDerivationNeighbor(pointIndex, indexNeighborPoint, diffEv);
                            else
                                throw new Exception("Error: TE hole current at Point " + pointIndex);

                        }
                        else
                            throw new Exception("Error: TE hole current at Point " + pointIndex);


                    }
                default:
                    throw new Exception("Missing Type of hole current at Point " + pointIndex);
            }
        }





        /// <summary>
        /// electron current, scharfetter gummel current modified by a factor depending on carrier statistik
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jElectron(int index1, int index2)
        {
            double boltzmannCurrent = -physConstants.e * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n)
                / 2 * U_T / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                * Misc.Bernoulli(-((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T)
                * mesh.finiteElements[index1].nDensity(this)
                    * (1 - ((mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                    * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T))));
            //Console.WriteLine("Index: " + index1 + " nach " + index2 + " Eta: " + mesh.points[index1].eta(MeshpointSemiconductor.EtaChargeType.eta_Electrons, this) + " , " + mesh.points[index2].eta(MeshpointSemiconductor.EtaChargeType.eta_Electrons, this));
            //Console.WriteLine("boltzmannCurrent: " + boltzmannCurrent);
            //Console.WriteLine("FermiDiracCurrent: " + jElectronFermiDirac(index1, index2));
            //return jElectronFermiDirac(index1, index2);
            return /*modificationFactorCurrent(index1, index2)*/  jElectronBoltzmann(index1, index2);

        }
        /// <summary>
        /// electron current in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jElectronBoltzmann(int index1, int index2)
        {
            double boltzmannCurrent = -physConstants.e * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n)
                    / 2 * U_T / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                    * Misc.Bernoulli(-((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T)
                    * mesh.finiteElements[index1].nDensity(this)
                        * (1 - ((mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                        * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T))));
            return boltzmannCurrent;
        }

        /// <summary>
        /// hole current in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jHole(int index1, int index2)
        {
            double SGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                * (
                 Misc.Bernoulli(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDensity(this)
                 - Misc.Bernoulli((-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index2].pDensity(this)

                );

            double modSGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                * Misc.Bernoulli(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDensity(this)
                        * (1 - mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                        * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T));

            //Console.WriteLine(mesh.finiteElements[index1].index +  " , " + mesh.finiteElements[index1].material.name + ": SG " + SGcurrent +" , Mod SG " + modSGcurrent);

            return modSGcurrent;
        }
        public double modificationFactorCurrent(int index1, int index2)
        {
            var point1 = mesh.finiteElements[index1];
            var point2 = mesh.finiteElements[index2];

            double etaIndex1 = point1.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);
            double etaIndex2 = point2.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);

            double distrFctIndex1 = point1.distributionFuction(this, etaIndex1);// Math.Exp(FermiDiracSpline.ValueAt(etaIndex1));
            double distrFctIndex2 = point2.distributionFuction(this, etaIndex2);// Math.Exp(FermiDiracSpline.ValueAt(etaIndex2));

            double modifcationFactor = Math.Sqrt(distrFctIndex1 * distrFctIndex2 / (Math.Exp(etaIndex1) * Math.Exp(etaIndex2)));

            return modifcationFactor;

        }

        public double modificationFactorDerivationIndex1(int index1, int index2, FiniteElementSemiconductor.EtaChargeType etaChargeType, BasicLib.DerivationVariable derivationVariable)
        {
            var point1 = mesh.finiteElements[index1];
            var point2 = mesh.finiteElements[index2];

            double etaIndex1 = point1.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);
            double etaIndex2 = point2.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);

            double distrFctIndex1 = point1.distributionFuctionDerivation(this, etaIndex1);
            double distrFctIndex2 = point2.distributionFuctionDerivation(this, etaIndex2);
            double distrFctDerivation1 = point1.distributionFuctionDerivation(this, etaIndex1);
            double distrFctDerivation2 = point2.distributionFuctionDerivation(this, etaIndex2);

            double outerDerivative = 0.5 * (Math.Pow(modificationFactorCurrent(index1, index2), -0.5));
            double innerDerivative = distrFctIndex2 / Math.Exp(etaIndex2) * (distrFctDerivation1 / Math.Exp(etaIndex1) - distrFctIndex1 / Math.Exp(etaIndex1));

            return outerDerivative * innerDerivative;
        }


        //Fermi-Dirac modified Scharfetter Gummel Equations ---------------------------------------------------------------------------------------------------------------------------------------
        public double jElectronFermiDirac(int index1, int index2)
        {
            var point1 = mesh.finiteElements[index1];
            var point2 = mesh.finiteElements[index2];
            double prefactor = BasicLib.physConstants.e * (point1.material.propertiesSemiconductor.mu_n + point2.material.propertiesSemiconductor.mu_n)
                / 2 * U_T / point1.position.DistanceTo(point2.position);

            double geomAverageNc = Math.Sqrt(point2.material.propertiesSemiconductor.Nc(T) * point1.material.propertiesSemiconductor.Nc(T));

            double etaIndex1 = point1.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);// (mesh.points[index1].phi + mesh.points[index1].phi_n - mesh.points[index1].material.propertiesSemiconductor.chemicalPotential) / U_T;
            double etaIndex2 = point2.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);//(mesh.points[index2].phi + mesh.points[index2].phi_n - mesh.points[index2].material.propertiesSemiconductor.chemicalPotential) / U_T;


            double deltaPhi = (point2.phi - point1.phi) / U_T;

            // alternative current modification via modification of the bernoulli equation:
            //double diffusionEnhancement = etaIndex2 - etaIndex1 / (Math.Log(distrFctIndex2) - Math.Log(distrFctIndex1));
            //double BernoulliEnhanced = Misc.Bernoulli(-((mesh.points[index2].phi - mesh.points[index1].phi)) / (U_T * diffusionEnhancement));


            /*
            Console.WriteLine("eta: " + etaIndex1);
            Console.WriteLine("n Boltzmann " + mesh.points[index1].nDensity(this));
            Console.WriteLine("n Fermi " + mesh.points[index1].material.propertiesSemiconductor.Nc * distrFctIndex1);
            Console.WriteLine("etaIndex1: " + etaIndex1);
            Console.WriteLine("etaIndex2: " + etaIndex2);
            Console.WriteLine("distrFctIndex1: " + distrFctIndex1);
            Console.WriteLine("distrFctIndex2: " + distrFctIndex2);
            Console.WriteLine("diffusionEnhancement: " + diffusionEnhancement);
            return -prefactor * diffusionEnhancement * Nc1 * distrFctIndex1 * BernoulliEnhanced * (1 - (Nc2*distrFctIndex2)/(Nc1 * distrFctIndex1) * Math.Exp(-deltaPhi/diffusionEnhancement));

            return Misc.Bernoulli(-((mesh.points[index2].phi - mesh.points[index1].phi)) / (U_T*diffusionEnhancement))  * mesh.points[index1].nDensity(this)
                    * (1 - ((mesh.points[index2].material.propertiesSemiconductor.Nc / mesh.points[index1].material.propertiesSemiconductor.Nc
                    * Math.Exp((mesh.points[index2].phi_n - mesh.points[index1].phi_n) / U_T))));
            */


            return -prefactor * geomAverageNc * modificationFactorCurrent(index1, index2) * Misc.Bernoulli(-deltaPhi) * Math.Exp(etaIndex1)
                * (1 - Math.Exp(-deltaPhi + etaIndex2 - etaIndex1));
            //*( 1 - Math.Exp(-deltaPhi) * Math.Exp(etaIndex2)/Math.Exp(etaIndex1));


        }

        public double jnFD_derivationPhi1(int index1, int index2)
        {
            /*
            var point1 = mesh.finiteElements[index1];
            var point2 = mesh.finiteElements[index2];
            double prefactor = BasicLib.physConstants.e * (point1.material.propertiesSemiconductor.mu_n + point2.material.propertiesSemiconductor.mu_n)
                / 2 * U_T / point1.position.DistanceTo(point2.position);

            double geomAverageNc = Math.Sqrt(point2.material.propertiesSemiconductor.Nc * point1.material.propertiesSemiconductor.Nc);

            double etaIndex1 = point1.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);// (mesh.points[index1].phi + mesh.points[index1].phi_n - mesh.points[index1].material.propertiesSemiconductor.chemicalPotential) / U_T;
            double etaIndex2 = point2.eta(FiniteElementSemiconductor.EtaChargeType.eta_Electrons, this);//(mesh.points[index2].phi + mesh.points[index2].phi_n - mesh.points[index2].material.propertiesSemiconductor.chemicalPotential) / U_T;

            var DerivVariable = BasicLib.DerivationVariable.Phi;

            double deltaPhi = (point2.phi - point1.phi) / U_T;

            return -prefactor * geomAverageNc * modificationFactorDerivation(index1, index2, FiniteElementSemiconductor.EtaChargeType.eta_Electrons, DerivVariable);
            */
            var etaChargeType = FiniteElementSemiconductor.EtaChargeType.eta_Electrons;
            var point1 = mesh.finiteElements[index1];
            var derivVariable = BasicLib.DerivationVariable.Phi;

            double etaDerivative = point1.etaDerivation(derivVariable, etaChargeType, this);

            return (modificationFactorDerivationIndex1(index1, index2, etaChargeType, derivVariable) * jElectronBoltzmann(index1, index2)
                + modificationFactorCurrent(index1, index2) * jnBoltzmannDerivationPhiIndex1(index1, index2)) + etaDerivative;
        }
        public double jnBoltzmannDerivationPhiIndex1(int index1, int index2)
        {
            return -BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                    * (1 - (mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T) * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T)))
                    * (Misc.BernoulliDerivation(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * mesh.finiteElements[index1].nDensity(this) / U_T
                    + Misc.Bernoulli(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * mesh.finiteElements[index1].nDerivationPhi(this));
        }


        /// <summary>
        /// derivation (Phi) of electron current at point in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jnDerivationPhiIndex1(int index1, int index2)
        {


            double boltzmannDerivation = -BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                    * (1 - (mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                    * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T)))
                    * (Misc.BernoulliDerivation(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * mesh.finiteElements[index1].nDensity(this) / U_T
                    + Misc.Bernoulli(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * mesh.finiteElements[index1].nDerivationPhi(this));
            //Console.WriteLine("Boltzmann" + boltzmannDerivation);
            //Console.WriteLine("Fermi-Dirac" + jnFD_derivationPhi1(index1, index2));
            //return jnFD_derivationPhi1(index1, index2);
            return boltzmannDerivation;
        }

        /// <summary>
        /// derivation (Phi) of electron current at neighbor in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jnDerivationPhiIndex2(int index1, int index2)
        {
            return -BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)) *
                    (1 - (mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                    * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T)))
                    * (Misc.BernoulliDerivation(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T)
                    * (-mesh.finiteElements[index1].nDensity(this) / U_T));

        }
        /// <summary>
        /// derivation (Phi_n) of electron current at point in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jnDerivationPsinIndex1(int index1, int index2)
        {
            return -BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * Misc.Bernoulli(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T)
                * ((mesh.finiteElements[index1].nDerivationPhi_n(this)) * (1 - (mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T)))
                    - (mesh.finiteElements[index1].nDensity(this) / U_T) * (-mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                    * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T)));

        }
        /// <summary>
        /// derivation (Phi_n) of electron current at neighbor in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jnDerivationPsinIndex2(int index1, int index2)
        {
            return -BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_n + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_n) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * Misc.Bernoulli(-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * (mesh.finiteElements[index1].nDensity(this) / U_T)
                    * (-mesh.finiteElements[index2].material.propertiesSemiconductor.Nc(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nc(T)
                    * Math.Exp((mesh.finiteElements[index2].phi_n - mesh.finiteElements[index1].phi_n) / U_T));
        }


        /// <summary>
        /// derivation (Phi) of hole current at point in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jpDerivationPhiIndex1(int index1, int index2)
        {
            double SGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
    / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
    * (
     Misc.Bernoulli(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDerivationPhi(this) + (-1 / U_T) * Misc.BernoulliDerivation(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDensity(this)
     - (1 / U_T) * Misc.BernoulliDerivation((-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index2].pDensity(this)

    );

            // return SGcurrent;

            return BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * (1 - mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T))
                    * (Misc.BernoulliDerivation((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T)
                    * (-mesh.finiteElements[index1].pDensity(this) / U_T)
                            + Misc.Bernoulli((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T)
                            * (mesh.finiteElements[index1].pDerivationPhi(this)));
        }
        /// <summary>
        /// derivation (Phi) of hole current at neighbor in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jpDerivationPhiIndex2(int index1, int index2)
        {
            double SGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
               / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
               * (
                 (1 / U_T) * Misc.BernoulliDerivation(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDensity(this)
                - (Misc.Bernoulli((-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index2].pDerivationPhi(this) + (-1 / U_T) * Misc.BernoulliDerivation((-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index2].pDensity(this))

               );

            //return SGcurrent;

            return BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * (1 - mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T))
                    * (mesh.finiteElements[index1].pDensity(this) / U_T)
                    * Misc.BernoulliDerivation((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T);
        }
        /// <summary>
        /// derivation (Phi_p) of hole current at point in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jpDerivationPsipIndex1(int index1, int index2)
        {

            double SGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                * (
                 Misc.Bernoulli(((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index1].pDerivationPhi_p(this)

                );

            //return SGcurrent;

            return BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * Misc.Bernoulli((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) *
                ((mesh.finiteElements[index1].pDerivationPhi_p(this)) * (1 - mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T))
                    + (mesh.finiteElements[index1].pDensity(this) / U_T) * (-mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                    * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T)));

        }
        /// <summary>
        /// derivation (Phi_p) of hole current at neighbor in Scharfetter Gummel formalism (modified)
        /// </summary>
        /// <param name="index1">point index</param>
        /// <param name="index2"> neighbor index</param>
        /// <returns></returns>
        public double jpDerivationPsipIndex2(int index1, int index2)
        {
            double SGcurrent = BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                    / mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position)
                    * (

                     -Misc.Bernoulli((-(mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi)) / U_T) * mesh.finiteElements[index2].pDerivationPhi_p(this)

                    );

            //return SGcurrent;

            return BasicLib.physConstants.e * U_T * (mesh.finiteElements[index1].material.propertiesSemiconductor.mu_p + mesh.finiteElements[index2].material.propertiesSemiconductor.mu_p) / 2
                / (mesh.finiteElements[index1].position.DistanceTo(mesh.finiteElements[index2].position))
                * Misc.Bernoulli((mesh.finiteElements[index2].phi - mesh.finiteElements[index1].phi) / U_T) * (-mesh.finiteElements[index1].pDensity(this) / U_T)
                * (-mesh.finiteElements[index2].material.propertiesSemiconductor.Nv(T) / mesh.finiteElements[index1].material.propertiesSemiconductor.Nv(T)
                * Math.Exp(-(mesh.finiteElements[index2].phi_p - mesh.finiteElements[index1].phi_p) / U_T));

        }


        //Equations for thermionic emission at abrupt heterojunctions-----------------------------------------------------------------------------------

        /// <summary>
        /// Returns the thermionic emission current of electrons at a heterointerface
        /// </summary>
        /// <param name="indexPoint">index of calculation point</param>
        /// <param name="indexNeighbor">index of the neigbor across the heterojunction</param>
        /// <param name="diffEc">barrier height of the hetero junction (difference between Ec)</param>
        /// <returns></returns>
        public double ThermionicEmissionCurrentElectron(int indexPoint, int indexNeighbor, double diffEc)
        {

            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;

            Console.ForegroundColor = ConsoleColor.Green;
            //Console.WriteLine("diffEc: " + diffEc);
            Console.ForegroundColor = ConsoleColor.Gray;

            if (diffEc > 0)
            {
                double a = (physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDensity(this)
                          - physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDensity(this) * Math.Exp(-Math.Abs(diffEc) / U_T));
                // double a = physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDensity(this) *
                //   (1 -  (TwoVNDevided(indexNeighbor, indexPoint) * Math.Exp(-Math.Abs(diffEc) / U_T)));
                /*
                Console.WriteLine("v neighb: " + mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity);
                Console.WriteLine("n neighb: " + mesh.points[indexNeighbor].nDensity(this));
                Console.WriteLine("v point: " + mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity);
                Console.WriteLine("n point: " + mesh.points[indexPoint].nDensity(this));
                Console.WriteLine("ELECTRON Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE Strom: " + a);
                Console.WriteLine("Strom alt " + (physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDensity(this)
                           - physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDensity(this) * Math.Exp(-Math.Abs(diffEc) / U_T)));
                */
                return a;
            }
            else if (diffEc < 0)
            {
                double b = physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDensity(this) * Math.Exp(-Math.Abs(diffEc) / U_T)
                   - (physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDensity(this));
                //double b = physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDensity(this)
                //   * ( TwoVNDevided(indexPoint, indexNeighbor) * Math.Exp(-Math.Abs(diffEc) / U_T)- 1 );
                /*
                Console.WriteLine("\n");
                Console.WriteLine("ELECTRON Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE Strom: " + b);
                Console.WriteLine("1: " + physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDensity(this));
                Console.WriteLine("exp: " + Math.Exp(-Math.Abs(diffEc) / U_T));
                Console.WriteLine("2: " + (physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDensity(this)));
                */
                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Blue;
                Console.WriteLine("Heterojunction überprüfen TE Electron!!! Bei Punkt " + indexPoint + " mit diffEc: " + diffEc);
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of electrons with respect to Phi at the calculation point
        /// </summary>
        /// <param name="indexPoint">index of calculation point</param>
        /// <param name="diffEc">barrier height of the hetero junction (difference between Ec)</param>
        /// <returns></returns>
        public double JnTE_PhiDerivationPoint(int indexPoint, int indexNeighbor, double diffEc)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEc > 0)
            {
                //return -physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDerivationPhi(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
                return -physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDerivationPhi(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
            }
            else if (diffEc < 0)
            {
                //return -(physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDerivationPhi(this));
                return -(physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDerivationPhi(this));
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Blue;
                Console.WriteLine("Heterojunction überprüfen TE Electron Phi Index!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of electrons with respect to Phi at the neighboring point
        /// </summary>
        /// <param name="indexNeighbor">index of the neigbor across the heterojunction</param>
        /// <param name="diffEc">barrier height of the hetero junction (difference between Ec)</param>
        /// <returns></returns>
        public double JnTE_PhiDerivationNeighbor(int indexPoint, int indexNeighbor, double diffEc)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEc > 0)
            {
                //return (physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDerivationPhi(this));
                return (physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDerivationPhi(this));
            }
            else if (diffEc < 0)
            {
                //return physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDerivationPhi(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
                return physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDerivationPhi(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Blue;
                Console.WriteLine("Heterojunction überprüfen  TE Electron Phi Neighbor!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of electrons with respect to Phi_n at the calculation point
        /// </summary>
        /// <param name="indexPoint">index of calculation point</param>
        /// <param name="diffEc">barrier height of the hetero junction (difference between Ec)</param>
        /// <returns></returns>
        public double JnTE_PhinDerivationPoint(int indexPoint, int indexNeighbor, double diffEc)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEc > 0)
            {
                //return -physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDerivationPhi_n(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
                return -physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDerivationPhi_n(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
            }
            else if (diffEc < 0)
            {
                //return -(physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexPoint].nDerivationPhi_n(this));
                return -(physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexPoint].nDerivationPhi_n(this));
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Blue;
                Console.WriteLine("Heterojunction überprüfen  TE Electron Phin Index!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of electrons with respect to Phi_n at the neighboring point
        /// </summary>
        /// <param name="indexNeighbor">index of the neigbor across the heterojunction</param>
        /// <param name="diffEc">barrier height of the hetero junction (difference between Ec)</param>
        /// <returns></returns>
        public double JnTE_PhinDerivationNeighbor(int indexPoint, int indexNeighbor, double diffEc)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEc > 0)
            {
                //return (physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDerivationPhi_n(this));
                return (physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDerivationPhi_n(this));
            }
            else if (diffEc < 0)
            {
                //return physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.points[indexNeighbor].nDerivationPhi_n(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
                return physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.electronThermalVelocity * mesh.finiteElements[indexNeighbor].nDerivationPhi_n(this) * Math.Exp(-Math.Abs(diffEc) / U_T);
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Blue;
                Console.WriteLine("Heterojunction überprüfen  TE Electron Phin Neighbor!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        double holeDirection1 = 1;
        double holeDirection2 = 1;
        /// <summary>
        /// Returns the thermionic emission current of holes at a heterointerface
        /// </summary>
        /// <param name="indexPoint">index of calculation point</param>
        /// <param name="indexNeighbor">index of the neigbor across the heterojunction</param>
        /// <param name="diffEv">barrier height of the hetero junction (difference between Ev)</param>
        /// <returns></returns>
        public double ThermionicEmissionCurrentHole(int indexPoint, int indexNeighbor, double diffEv)
        {

            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;

            Console.ForegroundColor = ConsoleColor.DarkGreen;
            //Console.WriteLine("\n diffEv: " + diffEv);
            Console.ForegroundColor = ConsoleColor.Gray;

            if (diffEv > 0)
            {
                double a = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDensity(this) * Math.Exp(-Math.Abs(diffEv) / U_T)
                        + physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDensity(this);

                //Console.WriteLine("Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE HOLE 1 Strom: ---------------------- " + a);

                a *= holeDirection1;
                /*
                Console.WriteLine(" > 0");
                Console.WriteLine("HOLE Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE Strom: " + a);
                Console.WriteLine("1: " + physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.points[indexNeighbor].pDensity(this));
                Console.WriteLine("exp: " + Math.Exp(-Math.Abs(diffEv) / U_T));
                Console.WriteLine("2: " + physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.points[indexPoint].pDensity(this));
                */
                return a;
            }
            else if (diffEv < 0)
            {
                double b = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDensity(this)
                        + physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDensity(this) * Math.Exp(-Math.Abs(diffEv) / U_T);
                //Console.WriteLine("Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE HOLE 2 Strom:  ---------------------- " + b);
                b *= holeDirection2;
                /*
                Console.WriteLine(" < 0");
                Console.WriteLine("HOLE Von Punkt " + indexPoint + " zu Punkt " + indexNeighbor + " TE Strom: " + b);
                Console.WriteLine("1: " + physConstants.e * mesh.points[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.points[indexNeighbor].pDensity(this));
                Console.WriteLine("exp: " + Math.Exp(-Math.Abs(diffEv) / U_T));
                Console.WriteLine("2: " + physConstants.e * mesh.points[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.points[indexPoint].pDensity(this));
                */
                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.DarkRed;
                Console.WriteLine("Heterojunction überprüfen TE Strom an Punkt " + indexPoint + " mit diffEv " + diffEv + "!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of holes with respect to Phi at the calculation point
        /// </summary>
        /// <param name="indexPoint">index of the calculation point</param>
        /// <param name="diffEv">barrier height of the hetero junction (difference between Ev)</param>
        /// <returns></returns>
        public double JpTE_PhiDerivationPoint(int indexPoint, int indexNeighbor, double diffEv)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEv > 0)
            {
                double a = physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDerivationPhi(this);
                a *= holeDirection1;

                return a;
            }
            else if (diffEv < 0)
            {
                double b = physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDerivationPhi(this) * Math.Exp(-Math.Abs(diffEv) / U_T);
                b *= holeDirection2;

                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine("Heterojunction überprüfen!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of holes with respect to Phi at the neighboring point
        /// </summary>
        /// <param name="indexNeighbor">index of the neighboring point</param>
        /// <param name="diffEv">barrier height of the hetero junction (difference between Ev)</param>
        /// <returns></returns>
        public double JpTE_PhiDerivationNeighbor(int indexPoint, int indexNeighbor, double diffEv)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexNeighbor
            if (diffEv > 0)
            {
                double a = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDerivationPhi(this) * Math.Exp(-Math.Abs(diffEv) / U_T);
                a *= holeDirection1;

                return a;
            }
            else if (diffEv < 0)
            {
                double b = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDerivationPhi(this);
                b *= holeDirection2;

                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine("Heterojunction überprüfen!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of holes with respect to Phi_p at the calculation point
        /// </summary>
        /// <param name="indexPoint">index of the calculation point</param>
        /// <param name="diffEv">barrier height of the hetero junction (difference between Ev)</param>
        /// <returns></returns>
        public double JpTE_PhipDerivationPoint(int indexPoint, int indexNeighbor, double diffEv)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexPoint
            if (diffEv > 0)
            {
                double a = physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDerivationPhi_p(this);
                a *= holeDirection1;

                return a;
            }
            else if (diffEv < 0)
            {
                double b = physConstants.e * mesh.finiteElements[indexPoint].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexPoint].pDerivationPhi_p(this) * Math.Exp(-Math.Abs(diffEv) / U_T);
                b *= holeDirection2;

                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine("Heterojunction überprüfen!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

        /// <summary>
        /// Returns the derivation of the TE current of holes with respect to Phi_p at the neighboring point
        /// </summary>
        /// <param name="indexNeighbor">index of the neighboring point</param>
        /// <param name="diffEv">barrier height of the hetero junction (difference between Ev)</param>
        /// <returns></returns>
        public double JpTE_PhipDerivationNeighbor(int indexPoint, int indexNeighbor, double diffEv)
        {
            if (mesh.finiteElements[indexPoint].neighbors[mesh.finiteElements[indexPoint].neighbors.FindIndex(x => x.index == indexNeighbor)].edgeSize < 1e-15)
                return 0;
            //Ableitungen nach indexNeighbor
            if (diffEv > 0)
            {
                double a = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDerivationPhi_p(this) * Math.Exp(-Math.Abs(diffEv) / U_T);
                a *= holeDirection1;

                return a;
            }
            else if (diffEv < 0)
            {
                double b = -physConstants.e * mesh.finiteElements[indexNeighbor].material.propertiesSemiconductor.holeThermalVelocity * mesh.finiteElements[indexNeighbor].pDerivationPhi_p(this);
                b *= holeDirection2;

                return b;
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine("Heterojunction überprüfen!!!");
                Console.ForegroundColor = ConsoleColor.Gray;
                return 0;
            }
        }

    }
}