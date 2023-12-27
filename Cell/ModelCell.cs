using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using Extreme.Mathematics;
using Extreme.Mathematics.LinearAlgebra;
using MoreLinq;
using Geometry;
using BasicLib;
using Extreme.Mathematics.Optimization;
using TransferMatrix;
using Database;
using Extreme.Mathematics.Algorithms;
using Extreme.Mathematics.EquationSolvers;
using Extreme.Mathematics.Curves;

namespace Cell
{
    /// <summary>
    /// Class, which represents a model of the Solarcell
    /// </summary>
    public class ModelCell
    {
        #region basics
        /// <summary>
        /// ID of the cell
        /// </summary>
        public int ID { get; private set; }
        /// <summary>
        /// Name of the cell
        /// </summary>
        public string name { get; private set; }
        /// <summary>
        /// determines, if this model is marked as a favorite model
        /// </summary>
        public bool isFavorite { get; set; } = false;
        /// <summary>
        /// temperature of this cell in Kelvin
        /// </summary>
        public double temperature { get; private set; }
        /// <summary>
        /// illumination intensity (1 means 1000W/m^2)
        /// </summary>
        public double illuminationInSuns { get; private set; }
        /// <summary>
        /// active are of this cell in m^2 (used for calculating efficiency, already includes edgeArea, must be multiplied with periodic factors)
        /// </summary>
        public double totalArea { get; private set; }
        #endregion

        #region preferences
        /// <summary>
        /// voltage potentials, at which the cell is operated
        /// </summary>
        public (double lower, double upper) operatingVoltage { get; private set; }
        /// <summary>
        /// determines, which side's potential is simulated or just set to the applied voltage or 0 for every mesh-cell
        /// </summary>
        public SimulationSelector simulationSelector { get; private set; }
        /// <summary>
        /// mode, how the voltage is changed
        /// </summary>
        public VoltageSweepMode voltageSweepMode { get; private set; }
        /// <summary>
        /// list of indexes of positions, segments and areas, which are connected to the external cell contact (selector = 1 => front contact, selector = 2 => back contact)
        /// </summary>
        public (List<(int index, int selector, List<double> resistanceAbsolute)> points,
            List<(int index, int selector, List<double> resistanceLinedensity)> segments,
            List<(int index, int selector, List<double> resistanceAreadensity)> regions) externalCellContacts
        { get; set; }
        #endregion

        #region characteristic curves and losses
        /// <summary>
        /// characteristic curve of the semiconductor material in the size of the whole cell area
        /// </summary>
        public CharacteristicCurve semiconductorCharacteristic { get; private set; }
        /// <summary>
        /// characteristic curve of the whole cell
        /// </summary>
        public CharacteristicCurve cellCharacteristic { get; set; }
        /// <summary>
        /// list of loss analyses in the same order as experimental points in cellCharacteristic
        /// </summary>
        public List<(double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms)> lossAnalyses { get; private set; }
        #endregion

        #region module preferences
        /// <summary>
        /// determines, whether this model is a module (true) or a cell (false)
        /// </summary>
        public bool isModule { get; private set; } = false;
        /// <summary>
        /// determines the module connections (bottom left points of the connection between both connection-edges)
        /// </summary>
        (Position P3, Position cell) pivotPositionsModuleConnection { get; set; }
        /// <summary>
        /// amount of cells in one single row, which form the module (voltage is multiplied by this factor)
        /// </summary>
        public int amountOfCellsInRow { get; private set; } = 1;
        /// <summary>
        /// factor, how many times broader the cell-strips are than in the given geometry-file (current is multiplied by this factor)
        /// </summary>
        public double strechAlongP1 { get; private set; } = 1;
        /// <summary>
        /// area of the edge regions in m^2 per simulated geometry (already divided by both periodic factors)
        /// </summary>
        public double edgeAreaPerSimulatedFraction { get; private set; } = 0;
        #endregion

        #region mesh
        /// <summary>
        /// meshing algorithm, which is used for this cell model
        /// </summary>
        public MeshingAlgorithm<FiniteElementCell, RegionCell> meshingAlgorithm { get; private set; }
        /// <summary>
        /// mesh with all meshpoints
        /// </summary>
        public Mesh<FiniteElementCell> mesh;
        /// <summary>
        /// geometry lines of this cell
        /// </summary>
        public GeometryFileData2D geometryFileData2D { get; private set; }
        #endregion

        #region optical TMM model
        /// <summary>
        /// power in W/m², which is use to calculate the efficiency of the cell
        /// </summary>
        public double totalPowerFromSpectrum { get; set; } = 1000;
        /// <summary>
        /// material before the stack for transfer matrix method
        /// </summary>
        #endregion

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Contructor, which sets the name, control panels, unit, cornerpoints, regions points and simulation holes from the preferences file
        /// </summary>
        /// <param name="name">Name of the cell</param>
        public ModelCell(string name, int ID, double temperature)
        {
            if (Misc.printInformation)
                Console.WriteLine();
            Misc.WriteDividingLine();
            Misc.WriteFormatedLine();
            Misc.WriteFormatedLine(">>> SIMULATION OF A CELL <<<");

            this.name = name;
            this.ID = ID;
            this.temperature = temperature;
            cellCharacteristic = new CharacteristicCurve(temperature);
        }

        // Create mesh (Delaunay and Voronoi) ███████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Create mesh (Delaunay and Voronoi)
        /// </summary>
        public void SetMesh(string[] geometryLines, int desiredAmountOfPoints, MeshingMethod meshingMethod, Mesh<FiniteElementCell> meshfromJson, int dimensions)
        {
            if (meshfromJson == null)
            {
                switch (meshingMethod)
                {
                    case MeshingMethod.quasiEquidistant_1D:
                        break;

                    case MeshingMethod.delaunayVoronoi_2D:
                        geometryFileData2D = new GeometryFileData2D(geometryLines);
                        isModule = CeckIfGeometryIsModule();
                        externalCellContacts = geometryFileData2D.boundaryConditions;
                        var delaunayVoronoi = new Meshing2D_DelaunayVoronoi<FiniteElementCell, RegionCell>(out mesh, geometryFileData2D); // Create four points, which are outside the simulation
                        delaunayVoronoi.ConstructContours(desiredAmountOfPoints); // Construct contour junctions and segments
                        mesh = delaunayVoronoi.AddRegularAndRandomlyShiftedPoints(mesh); // Add additional regular points
                        mesh = delaunayVoronoi.AddContours(mesh); // Add points for outer contours, simulation holes and regions
                        mesh = delaunayVoronoi.CreateVoronoiFromDelaunay(mesh); // Create Voronoi mesh from Delaunay triangulation
                        meshingAlgorithm = delaunayVoronoi;
                        Misc.WriteFormatedLine("Delaunay- und Voronoi-Triangulation done with " + mesh.nextAvailableFiniteElementIndex + " Elements."); // Output
                        break;

                    case MeshingMethod.quadtree_2D:
                        geometryFileData2D = new GeometryFileData2D(geometryLines);
                        isModule = CeckIfGeometryIsModule();
                        externalCellContacts = geometryFileData2D.boundaryConditions;
                        var quadtree = new Meshing2D_Quadtree<FiniteElementCell, RegionCell>(out mesh, geometryFileData2D, desiredAmountOfPoints); // Create four points, which are outside the simulation
                        meshingAlgorithm = quadtree;
                        Misc.WriteFormatedLine("Quadtree mesh done with " + mesh.nextAvailableFiniteElementIndex + " Elements."); // Output
                        break;

                    case MeshingMethod.delaunayVoronoi_3D:
                        break;
                }
            }
            else
            {
                isModule = CeckIfGeometryIsModule(geometryLines);
                mesh = meshfromJson;

                switch (dimensions)
                {
                    case 1:
                        Console.WriteLine("in case 1");
                        break;

                    case 2:
                        Console.WriteLine("in case 2");
                        GeometryFileData2D geometryFileData2D = new GeometryFileData2D(geometryLines);
                        externalCellContacts = geometryFileData2D.boundaryConditions;
                        meshingAlgorithm = new MeshingAlgorithm<FiniteElementCell, RegionCell>(geometryFileData2D, mesh.CalculatMeanElementDistance());
                        Misc.WriteFormatedLine("Mesh loaded with " + mesh.nextAvailableFiniteElementIndex + " Elements.");
                        break;
                }
            }
        }

        // set preferences ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// set global optical cell preferences
        /// </summary>
        /// <param name="opticMode">mode, how optics is calculated</param>
        /// <param name="spectrum">spectrum, which is used for TMM (can be null for other optic modes)</param>
        /// <param name="angleOfIncidence">incident angle in degree (0 = perpendicular incident > effective area is physical area, 90 = flat > no area in illumination)</param>
        /// <param name="sunElevation">height angle of the sun at the sky in degree (90 = sun in zenith)</param>
        public void SetOpticsOfGlobalCell(OpticMode opticMode, Spectrum spectrum, double illuminationInSuns, double angleOfIncidence = 0, double sunElevation = 90)
        {
            this.illuminationInSuns = illuminationInSuns;

            // change illumination
            spectrum = new Spectrum(spectrum.data.Select(d => (d.lambda, d.deltaLambda, d.spectralIntensityDensity * illuminationInSuns)).ToArray());
            foreach (var region in meshingAlgorithm.regions)
            {
                var pn = Data.GetPNjunctionFromID(region.pnJunction.ID);
                pn.characteristicCurve.currentPhoto *= illuminationInSuns;
                region.pnJunction = pn;
            }

            switch (opticMode)
            {
                case OpticMode.coefficient:
                    totalPowerFromSpectrum = 1000 * illuminationInSuns;
                    break;

                case OpticMode.lambertBeer:
                    totalPowerFromSpectrum = 1000 * illuminationInSuns;
                    break;

                case OpticMode.transferMatrixInCell:
                    totalPowerFromSpectrum = spectrum.totalPower;
                    break;
            }

            foreach (var region in meshingAlgorithm.regions.Where(r => r.type != pointType.hole))
            {
                region.CalculateOpticalCoefficients(opticMode, spectrum, angleOfIncidence, sunElevation);
                double opticalPrefactor = 1;
                opticalPrefactor *= 1 - region.opticFactor_additionalAtmosphereScattering;
                opticalPrefactor *= 1 - region.opticFactor_shading;
                opticalPrefactor *= region.opticFactor_effectiveArea;
                opticalPrefactor *= region.opticFactor_transmissionGrid;
                opticalPrefactor *= 1 - region.opticFactor_reflection - region.opticFactor_absorption.Where(m => !m.isAbsorber).Sum(m => m.factor) - region.opticFactor_transmission;
                Misc.WriteFormatedLine("  optical prefactor of region " + region.index + ": " + opticalPrefactor);
            }
            Misc.WriteFormatedLine();
        }
        /// <summary>
        /// set global electrical cell preferences
        /// </summary>
        /// <param name="simulationSelector">determines, wheter only front, back, both or no potential is simulated</param>
        /// <param name="lowerOperatingVoltage">voltage of the back contact (for single cell enter: 0)</param>
        public void SetElectricsOfGlobalCell(SimulationSelector simulationSelector, double lowerOperatingVoltage)
        {
            this.simulationSelector = simulationSelector;
            operatingVoltage = (lowerOperatingVoltage, lowerOperatingVoltage + 0.65); // only set preliminary till semiconductorCurve is set --> then upper voltage to MPP of SC

            //  ██╗ 
            //  ╚██╗ set external cell contacts
            //  ██╔╝
            //  ╚═╝
            //  ██╗ points
            //  ╚═╝
            foreach (var pointContact in externalCellContacts.points)
            {
                // select point, where the cell contact is in
                foreach (FiniteElementCell point in mesh.finiteElements.Values)
                {
                    if (meshingAlgorithm.contourJunctions[pointContact.index].position.InPolygon(point, true))
                    {
                        if (pointContact.selector == 1)
                        {
                            point.isExternalCellFrontContact = true;
                            point.contactResistanceExternalCellFrontContact = pointContact.resistanceAbsolute[0];
                        }
                        if (pointContact.selector == 2)
                        {
                            point.isExternalCellBackContact = true;
                            point.contactResistanceExternalCellBackContact = pointContact.resistanceAbsolute[0];
                        }
                        goto pointFound;
                    }
                }

                // if no point is found, take the nearest one
                int nearestPointIndex = mesh.finiteElements.Values.MinBy(p => meshingAlgorithm.contourJunctions[pointContact.index].position.DistanceTo(p.position)).First().index;
                if (pointContact.selector == 1)
                {
                    mesh.finiteElements[nearestPointIndex].isExternalCellFrontContact = true;
                    mesh.finiteElements[nearestPointIndex].contactResistanceExternalCellFrontContact = pointContact.resistanceAbsolute[0];
                }
                if (pointContact.selector == 2)
                {
                    mesh.finiteElements[nearestPointIndex].isExternalCellBackContact = true;
                    mesh.finiteElements[nearestPointIndex].contactResistanceExternalCellBackContact = pointContact.resistanceAbsolute[0];
                }
            pointFound:;
            }

            //  ██╗ segments
            //  ╚═╝
            foreach (var point in mesh.finiteElements.Values)
                foreach (var segment in this.externalCellContacts.segments)
                    if (point.indexOfBorderElementCreatedFrom != -1)
                        if (segment.index == point.indexOfBorderElementCreatedFrom)
                        {
                            if (segment.selector == 1)
                            {
                                point.isExternalCellFrontContact = true;
                                point.contactResistanceExternalCellFrontContact = segment.resistanceLinedensity[0]
                                    / point.corners.Where(c => c.isOnContour).First().position.DistanceTo(point.corners.Where(c => c.isOnContour).Last().position);
                            }
                            if (segment.selector == 2)
                            {
                                point.isExternalCellBackContact = true;
                                point.contactResistanceExternalCellBackContact = segment.resistanceLinedensity[0]
                                    / point.corners.Where(c => c.isOnContour).First().position.DistanceTo(point.corners.Where(c => c.isOnContour).Last().position);
                            }
                        }

            //  ██╗ regions
            //  ╚═╝
            foreach (var point in mesh.finiteElements.Values)
                foreach (var region in this.externalCellContacts.regions)
                    if (point.position.InPolygon(meshingAlgorithm.regions[region.index].orderedPoints.Select(p => p.position).ToList(), true))
                    {
                        if (region.selector == 1)
                        {
                            point.isExternalCellFrontContact = true;
                            point.contactResistanceExternalCellFrontContact = region.resistanceAreadensity[0] / point.size;
                        }
                        if (region.selector == 2)
                        {
                            point.isExternalCellBackContact = true;
                            point.contactResistanceExternalCellBackContact = region.resistanceAreadensity[0] / point.size;
                        }
                    }
        }
        /// <summary>
        /// sets the global variables to the single meshspoints
        /// </summary>
        /// <param name="TOdensityArray">array for topological optimization</param>
        public void SetPreferencesOfSingleMeshpoints(double[] TOdensityArray = null, bool linkTopAndBottomElementsWithPeriodicBoundaryConditions = false)
        {
            //  ██╗ 
            //  ╚██╗ set parameters of each single point
            //  ██╔╝
            //  ╚═╝
            for (int i = 0; i < mesh.finiteElements.Count; i++)
            {
                if (TOdensityArray == null)
                    mesh.finiteElements[i].RefreshParametersFirst(this);
                else
                    mesh.finiteElements[i].RefreshParametersFirst(this, TOdensityArray[i]);
            }
            for (int i = 0; i < mesh.finiteElements.Count; i++)
                mesh.finiteElements[i].RefreshParametersSecond(mesh);

            //  ██╗ 
            //  ╚██╗ calculate total active cell area
            //  ██╔╝
            //  ╚═╝
            totalArea = mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Sum(p => p.size);

            //  ██╗ 
            //  ╚██╗ use most used pn junction as semiconductor curve for reference (preliminary with the whole delaunay outer contour. later regions with countsAsActiveArea are subtracted, but no influence on MPP voltage)
            //  ██╔╝
            //  ╚═╝
            RegionCell mostUsedRegion = meshingAlgorithm.regions.Where(r => r.type != pointType.hole).GroupBy(r => r.pnJunction).OrderByDescending(grp => grp.Count()).First().First();
            semiconductorCharacteristic = new CharacteristicCurve(temperature,
                mostUsedRegion.pnJunction.characteristicCurve.currentPhoto * totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.currentSaturation * totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.diode1IdealityFactor,
                mostUsedRegion.pnJunction.characteristicCurve.Rseries == 0 ? double.NaN : mostUsedRegion.pnJunction.characteristicCurve.Rseries / totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.Rshunt / totalArea);

            if (linkTopAndBottomElementsWithPeriodicBoundaryConditions)
                LinkTopAndBottomElementsWithPeriodicBoundaryConditions();
        }
        /// <summary>
        /// links the top and the bottom elements with periodic boundary conditions
        /// </summary>
        private void LinkTopAndBottomElementsWithPeriodicBoundaryConditions()
        {
            int segmentIndexTop = meshingAlgorithm.contourSegments.MaxBy(s => s.firstAdjacentContourJunction.position.y + s.secondAdjacentContourJunction.position.y).First().index;
            int segmentIndexBottom = meshingAlgorithm.contourSegments.MinBy(s => s.firstAdjacentContourJunction.position.y + s.secondAdjacentContourJunction.position.y).First().index;
            foreach (var feTop in mesh.finiteElements.Values.Where(f => f.indexOfBorderElementCreatedFrom == segmentIndexTop))
            {
                // find matching element
                var feBottom = mesh.finiteElements.Values.Where(f => f.indexOfBorderElementCreatedFrom == segmentIndexBottom).MinBy(f => f.position.DistanceTo(feTop.position)).First();
                double distance = feTop.position.DistanceTo(meshingAlgorithm.contourSegments[segmentIndexTop].lineSegment) + feBottom.position.DistanceTo(meshingAlgorithm.contourSegments[segmentIndexBottom].lineSegment);
                double edgeLength = feTop.borderEdgeSize;

                if (feTop.neighbors.Any(n => n.index == feBottom.index))
                    continue;

                //  ██╗ top element
                //  ╚═╝

                // add to neighbors
                feTop.neighbors.Add((feBottom.index, distance, feTop.borderEdgeSize));

                // add to front resistors
                feTop.resistorsFrontToEdgeContact.Add(feTop.frontContact.propertiesContact.GetSpecificResistivity(feTop.thicknessFrontContact) * distance / (2 * feTop.thicknessFrontContact * edgeLength));
                if (feTop.frontGrid == null)
                    feTop.resistorsFrontToEdgeGrid.Add(double.NaN);
                else
                {
                    if (double.IsNaN(feTop.frontGridDensity)) // if in topology optimization
                        feTop.resistorsFrontToEdgeGrid.Add(feTop.frontGrid.propertiesContact.GetSpecificResistivity(feTop.thicknessFrontGrid) * distance / (2 * feTop.thicknessFrontGrid * edgeLength));
                    else
                    {
                        if (feTop.frontGridDensity == 0)
                            feTop.resistorsFrontToEdgeGrid.Add(double.NaN);
                        else
                            feTop.resistorsFrontToEdgeGrid.Add(feTop.frontGrid.propertiesContact.GetSpecificResistivity(feTop.thicknessFrontGrid) * distance / (2 * feTop.thicknessFrontGrid * edgeLength) * Misc.SIMPfunction_Conductivity(feTop.frontGridDensity));
                    }
                }
                if (double.IsNaN(feTop.resistorsFrontToEdgeGrid.Last()))
                    feTop.resistorsFrontToEdge.Add(feTop.resistorsFrontToEdgeContact.Last());
                else
                    feTop.resistorsFrontToEdge.Add(Misc.ParallelResistance(feTop.resistorsFrontToEdgeContact.Last(), feTop.resistorsFrontToEdgeGrid.Last()));

                // add to back resistors
                feTop.resistorsBackToEdgeContact.Add(feTop.backContact.propertiesContact.GetSpecificResistivity(feTop.thicknessBackContact) * distance / (2 * feTop.thicknessBackContact * edgeLength));
                if (feTop.backGrid == null)
                {
                    feTop.resistorsBackToEdgeGrid.Add(double.NaN);
                    feTop.resistorsBackToEdge.Add(feTop.resistorsBackToEdgeContact.Last());
                }
                else
                {
                    feTop.resistorsBackToEdgeGrid.Add(feTop.backGrid.propertiesContact.GetSpecificResistivity(feTop.thicknessBackGrid) * distance / (2 * feTop.thicknessBackGrid * edgeLength));
                    if (double.IsNaN(feTop.resistorsBackToEdgeGrid.Last()))
                        feTop.resistorsBackToEdge.Add(feTop.resistorsBackToEdgeContact.Last());
                    else
                        feTop.resistorsBackToEdge.Add(Misc.ParallelResistance(feTop.resistorsBackToEdgeContact.Last(), feTop.resistorsBackToEdgeGrid.Last()));
                }

                //  ██╗ bottom element
                //  ╚═╝

                // add to neighbors
                feBottom.neighbors.Add((feTop.index, distance, edgeLength));

                // add to front resistors
                feBottom.resistorsFrontToEdgeContact.Add(feBottom.frontContact.propertiesContact.GetSpecificResistivity(feBottom.thicknessFrontContact) * distance / (2 * feBottom.thicknessFrontContact * edgeLength));
                if (feBottom.frontGrid == null)
                    feBottom.resistorsFrontToEdgeGrid.Add(double.NaN);
                else
                {
                    if (double.IsNaN(feBottom.frontGridDensity)) // if in topology optimization
                        feBottom.resistorsFrontToEdgeGrid.Add(feBottom.frontGrid.propertiesContact.GetSpecificResistivity(feBottom.thicknessFrontGrid) * distance / (2 * feBottom.thicknessFrontGrid * edgeLength));
                    else
                    {
                        if (feBottom.frontGridDensity == 0)
                            feBottom.resistorsFrontToEdgeGrid.Add(double.NaN);
                        else
                            feBottom.resistorsFrontToEdgeGrid.Add(feBottom.frontGrid.propertiesContact.GetSpecificResistivity(feBottom.thicknessFrontGrid) * distance / (2 * feBottom.thicknessFrontGrid * edgeLength) * Misc.SIMPfunction_Conductivity(feBottom.frontGridDensity));
                    }
                }
                if (double.IsNaN(feBottom.resistorsFrontToEdgeGrid.Last()))
                    feBottom.resistorsFrontToEdge.Add(feBottom.resistorsFrontToEdgeContact.Last());
                else
                    feBottom.resistorsFrontToEdge.Add(Misc.ParallelResistance(feBottom.resistorsFrontToEdgeContact.Last(), feBottom.resistorsFrontToEdgeGrid.Last()));

                // add to back resistors
                feBottom.resistorsBackToEdgeContact.Add(feBottom.backContact.propertiesContact.GetSpecificResistivity(feBottom.thicknessBackContact) * distance / (2 * feBottom.thicknessBackContact * edgeLength));
                if (feBottom.backGrid == null)
                {
                    feBottom.resistorsBackToEdgeGrid.Add(double.NaN);
                    feBottom.resistorsBackToEdge.Add(feBottom.resistorsBackToEdgeContact.Last());
                }
                else
                {
                    feBottom.resistorsBackToEdgeGrid.Add(feBottom.backGrid.propertiesContact.GetSpecificResistivity(feBottom.thicknessBackGrid) * distance / (2 * feBottom.thicknessBackGrid * edgeLength));
                    if (double.IsNaN(feBottom.resistorsBackToEdgeGrid.Last()))
                        feBottom.resistorsBackToEdge.Add(feBottom.resistorsBackToEdgeContact.Last());
                    else
                        feBottom.resistorsBackToEdge.Add(Misc.ParallelResistance(feBottom.resistorsBackToEdgeContact.Last(), feBottom.resistorsBackToEdgeGrid.Last()));
                }
            }
            for (int i = 0; i < mesh.finiteElements.Count; i++)
                mesh.finiteElements[i].RefreshParametersSecond(mesh);
        }
        /// <summary>
        /// set initial guess for the module
        /// </summary>
        public void SetPreferencesOfModule(string[] geometryLines)
        {
            if (!isModule)
                return;

            this.simulationSelector = SimulationSelector.bothPotentials;

            // contact resistance
            double contactSeriesResistanceP2 = InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "contact series resistance P2:");
            foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.P2))
                point.contactSeriesResistanceP2 = contactSeriesResistanceP2 / point.size;

            //  ██╗ 
            //  ╚██╗ new additional area
            //  ██╔╝
            //  ╚═╝
            amountOfCellsInRow = Convert.ToInt32(InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "amount of cells:"));
            strechAlongP1 = InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "strech along P1:");
            edgeAreaPerSimulatedFraction = InputOutput.ReadSingleValueInArrayAfterStringAppearance(geometryLines, "edge area:") / (amountOfCellsInRow * strechAlongP1);

            // calculate total active cell area
            totalArea = mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Sum(p => p.size) + edgeAreaPerSimulatedFraction;

            // use most used pn junction as semiconductor curve for reference (preliminary with the whole delaunay outer contour. later regions with countsAsActiveArea are subtracted, but no influence on MPP voltage)
            RegionCell mostUsedRegion = meshingAlgorithm.regions.Where(r => r.type != pointType.hole).GroupBy(r => r.pnJunction).OrderByDescending(grp => grp.Count()).First().First();
            semiconductorCharacteristic = new CharacteristicCurve(temperature,
                mostUsedRegion.pnJunction.characteristicCurve.currentPhoto * totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.currentSaturation * totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.diode1IdealityFactor,
                mostUsedRegion.pnJunction.characteristicCurve.Rseries == 0 ? double.NaN : mostUsedRegion.pnJunction.characteristicCurve.Rseries / totalArea,
                mostUsedRegion.pnJunction.characteristicCurve.Rshunt / totalArea);

            // divide voltage
            semiconductorCharacteristic.factorToVoltage = amountOfCellsInRow;
            semiconductorCharacteristic.factorToCurrent = strechAlongP1;
            cellCharacteristic.factorToVoltage = amountOfCellsInRow;
            cellCharacteristic.factorToCurrent = strechAlongP1;

            // periodic boundary conditions
            int[] interconnectionSegmentsP3 = externalCellContacts.segments.Where(e => e.selector == 1).Select(e => e.index).ToArray();
            int[] interconnectionSegmentsCell = externalCellContacts.segments.Where(e => e.selector == 2).Select(e => e.index).ToArray();

            (int[] P3, int[] inCell) interconnectionSegments = (interconnectionSegmentsP3, interconnectionSegmentsCell);

            // get lower left point of all P3 noted segments
            List<ContourJunction> adjacentJunctionsP3 = new List<ContourJunction>();
            foreach (int segmentIndex in interconnectionSegments.P3)
            {
                adjacentJunctionsP3.Add(meshingAlgorithm.contourSegments[segmentIndex].firstAdjacentContourJunction);
                adjacentJunctionsP3.Add(meshingAlgorithm.contourSegments[segmentIndex].secondAdjacentContourJunction);
            }
            Position pivotPositionModuleConectionP3 = adjacentJunctionsP3.OrderBy(p => p.position.x).ThenBy(p => p.position.y).ToList().First().position;

            // get lower left point of all cell noted segments
            List<ContourJunction> adjacentJunctionsCell = new List<ContourJunction>();
            foreach (int segmentIndex in interconnectionSegments.inCell)
            {
                adjacentJunctionsCell.Add(meshingAlgorithm.contourSegments[segmentIndex].firstAdjacentContourJunction);
                adjacentJunctionsCell.Add(meshingAlgorithm.contourSegments[segmentIndex].secondAdjacentContourJunction);
            }
            Position pivotPositionModuleConectionCell = adjacentJunctionsCell.OrderBy(p => p.position.x).ThenBy(p => p.position.y).ToList().First().position;
            pivotPositionsModuleConnection = (pivotPositionModuleConectionP3, pivotPositionModuleConectionCell);

            // mark points as module connection
            foreach (var point in mesh.finiteElements.Values.Where(p => p.indexOfBorderElementCreatedFrom != -1 && interconnectionSegments.inCell.Any(e => e == p.indexOfBorderElementCreatedFrom)))
                point.hasModuleConnection_isInCell = true;
            foreach (var point in mesh.finiteElements.Values.Where(p => p.indexOfBorderElementCreatedFrom != -1 && interconnectionSegments.P3.Any(e => e == p.indexOfBorderElementCreatedFrom)))
            {
                point.hasModuleConnection_isInP3 = true;
                point.indexOfModuleConnectedPoint = mesh.finiteElements.Values.Where(p => p.hasModuleConnection_isInCell)
                    .MinBy(p => Math.Abs(p.position.DistanceTo(pivotPositionsModuleConnection.cell) - point.position.DistanceTo(pivotPositionsModuleConnection.P3))).First().index;
            }

            // external CELL contact not present in module (use P2 contact resistance herefore) -> set all to 0
            foreach (var point in mesh.finiteElements)
            {
                point.Value.contactResistanceExternalCellFrontContact = 0;
                point.Value.contactResistanceExternalCellBackContact = 0;
            }
        }
        /// <summary>
        /// Checks, whether the geometry file of this object is a module or a cell file (true = module, false = cell)
        /// </summary>
        public bool CeckIfGeometryIsModule(string[] geometryLines)
        {
            int materialsLine = InputOutput.GetLineOfStringInArray(geometryLines, "materials:") + 1;
            Dictionary<int, List<double>> materials = InputOutput.GetDoubleTupleOfStringArray(geometryLines, materialsLine);

            foreach (var m in materials.Values)
                if (m.Count > 11)
                    if (m[11] == 8002)
                        return true;
            return false;
        }
        /// <summary>
        /// Checks, whether the geometry file of this object is a module or a cell file (true = module, false = cell)
        /// </summary>
        public bool CeckIfGeometryIsModule()
        {
            foreach (var m in geometryFileData2D.materials.Values)
                if (m.Count > 11)
                    if (m[11] == 8002)
                        return true;
            return false;
        }

        // set initial Phi-vector ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// set initial guess
        /// </summary>
        public void SetInitialGuess(bool setOperatingVoltage = true)
        {
            if (setOperatingVoltage)
                operatingVoltage = (operatingVoltage.lower, operatingVoltage.lower + semiconductorCharacteristic.GetDataSetMaximumPowerPoint().voltage);

            for (int i = 0; i < mesh.finiteElements.Count; i++)
                mesh.finiteElements[i].SetInitialGuess(operatingVoltage, isModule, simulationSelector);
        }

        // solve differential equation system ███████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Solves the differential equation system and returns the produced power
        /// </summary>
        /// <param name="searchMPP">determines, how the equations will be solved (false = at certain voltage(s), true = search for MPP)</param>
        /// <param name="simulationVoltages">for search MPP == false: operating voltage(s), for searchMPP == true: array with only one value -> voltage precision of MPP</param>
        public void Solve(out (double voltage, double current, double power, double area, double efficiency) results,
            VoltageSweepMode voltageSweepMode, double[] simulationVoltages)
        {
            this.voltageSweepMode = voltageSweepMode;
            lossAnalyses = new List<(double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms)>();
            cellCharacteristic = new CharacteristicCurve(temperature);
            cellCharacteristic.factorToVoltage = amountOfCellsInRow;
            cellCharacteristic.factorToCurrent = strechAlongP1;

            //  ██╗ MPP data
            //  ╚═╝
            double maxPower = double.PositiveInfinity;
            double Vmpp = 0;
            double[] phiFrontVectorMPP = new double[mesh.nextAvailableFiniteElementIndex];
            double[] phiBackVectorMPP = new double[mesh.nextAvailableFiniteElementIndex];

            //  ██╗ calculate at certain voltages
            //  ╚═╝
            if (voltageSweepMode == VoltageSweepMode.singleVoltage || voltageSweepMode == VoltageSweepMode.multipleVoltages)
            {
                foreach (double voltage in simulationVoltages)
                {
                    double currentPower = Solve(new double[] { voltage + operatingVoltage.lower });
                    SaveMPPdataIfPowerIsMax(currentPower, voltage);
                }
            }

            //  ██╗ calculate till I = 2*Isc
            //  ╚═╝
            if (voltageSweepMode == VoltageSweepMode.autoIVcurve || voltageSweepMode == VoltageSweepMode.autoIVcurveAndMPP)
            {
                // calculate at 0V
                double currentPower = Solve(new double[] { operatingVoltage.lower });
                SaveMPPdataIfPowerIsMax(currentPower, 0);
                double Isc = -cellCharacteristic.experimentalData.First().current;
                double sweepVoltage = simulationVoltages[0] * amountOfCellsInRow;

                while (cellCharacteristic.experimentalData.Last().current < 2 * Isc)
                {
                    double currentSweepPower = Solve(new double[] { sweepVoltage + operatingVoltage.lower });
                    SaveMPPdataIfPowerIsMax(currentSweepPower, sweepVoltage);
                    sweepVoltage = Misc.RoundToSignificantDigits(sweepVoltage + simulationVoltages[0] * amountOfCellsInRow, 10);
                }
            }

            //  ██╗ search dynamically for MPP
            //  ╚═╝
            if (voltageSweepMode == VoltageSweepMode.searchMPP || voltageSweepMode == VoltageSweepMode.autoIVcurveAndMPP)
            {
                DownhillSimplex downhillSimplex;
                double initialVoltage;
                if (cellCharacteristic.experimentalData.Count < 2) // if no or only one volatage has been simulated yet -> semiconductor MPP as starting point
                    initialVoltage = semiconductorCharacteristic.GetDataSetMaximumPowerPoint().voltage + operatingVoltage.lower;
                else // best voltage from sweep as starting point
                {
                    var measuredDataOrdered = cellCharacteristic.experimentalData.OrderBy(d => d.power).ToArray();
                    initialVoltage = (measuredDataOrdered[0].voltage + measuredDataOrdered[1].voltage) / 2 * cellCharacteristic.factorToVoltage;
                }
                downhillSimplex = new DownhillSimplex(Solve, new double[] { initialVoltage }, null, 0.95);
                downhillSimplex.relativeDeltaParameterTolerance = double.PositiveInfinity;
                downhillSimplex.absoluteDeltaParameterTolerance = semiconductorCharacteristic.GetDataSetOpenCircuit().voltage * 0.003;
                double bestVoltage = downhillSimplex.Fit()[0];
                Solve(downhillSimplex.currentBestFitParameterSet);
                SaveMPPdataIfPowerIsMax(double.NegativeInfinity, bestVoltage);
            }

            //  ██╗ Set points to saved MPP
            //  ╚═╝
            operatingVoltage = (operatingVoltage.lower, operatingVoltage.lower + (Vmpp - operatingVoltage.lower) / amountOfCellsInRow);
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                mesh.finiteElements[i].phiFront = phiFrontVectorMPP[i];
                mesh.finiteElements[i].phiBack = phiBackVectorMPP[i];
            }
            foreach (FiniteElementCell finiteElement in mesh.finiteElements.Values)
                finiteElement.CalculateDerivedParameters(mesh);

            //  ██╗ calculate results
            //  ╚═╝
            results = CalculateOutputData();
            Misc.WriteFormatedLine();
            Misc.WriteDividingLine();

            //  ██╗ Save MPP data if the power is a new (negative) maximum power (power is negative)
            //  ╚═╝
            void SaveMPPdataIfPowerIsMax(double currentPower, double voltage)
            {
                if (currentPower < maxPower)
                {
                    Vmpp = voltage + operatingVoltage.lower;
                    maxPower = currentPower;
                    for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    {
                        phiFrontVectorMPP[i] = mesh.finiteElements[i].phiFront;
                        phiBackVectorMPP[i] = mesh.finiteElements[i].phiBack;
                    }
                }
            }
        }
        /// <summary>
        /// Solves the differential equation system
        /// </summary>
        double Solve(double[] parameters)
        {
            operatingVoltage = (operatingVoltage.lower, parameters[0] / amountOfCellsInRow);
            switch (simulationSelector)
            {
                case SimulationSelector.frontPotential:
                    for (int i = 0; i < mesh.finiteElements.Count; i++)
                    {
                        mesh.finiteElements[i].phiBackInit = operatingVoltage.lower;
                        mesh.finiteElements[i].phiBack = mesh.finiteElements[i].phiBackInit;
                    }
                    break;

                case SimulationSelector.backPotential:
                    for (int i = 0; i < mesh.finiteElements.Count; i++)
                    {
                        mesh.finiteElements[i].phiFrontInit = operatingVoltage.upper;
                        mesh.finiteElements[i].phiFront = mesh.finiteElements[i].phiFrontInit;
                    }
                    break;
            }

        //  ██╗ 
        //  ╚██╗ write initial guess from mesh to vector (Newtons method need an vector as input)
        //  ██╔╝
        //  ╚═╝
        beginningOfSolve:
            Vector<double> initialGuess;
            if (simulationSelector == SimulationSelector.bothPotentials) // if both potentials are simulated
            {
                initialGuess = Vector.Create<double>(2 * mesh.nextAvailableFiniteElementIndex);
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                {
                    initialGuess[2 * i] = mesh.finiteElements[i].phiFront;
                    initialGuess[2 * i + 1] = mesh.finiteElements[i].phiBack;
                }
            }
            else if (simulationSelector == SimulationSelector.frontPotential) // if only front potential is simulated
            {
                initialGuess = Vector.Create<double>(mesh.nextAvailableFiniteElementIndex);
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    initialGuess[i] = mesh.finiteElements[i].phiFront;
            }
            else // if only back potential is simulated
            {
                initialGuess = Vector.Create<double>(mesh.nextAvailableFiniteElementIndex);
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    initialGuess[i] = mesh.finiteElements[i].phiBack;
            }

            //  ██╗ 
            //  ╚██╗ Newtons method (solution is initial guess, 'getFunction' and 'getJacobi' are functions/instructions, which are given as parameter)
            //  ██╔╝
            //  ╚═╝
            double norm = 0;
            (initialGuess, norm) = NewtonMethod.Solve(initialGuess, GetFunctionAsync, GetJacobiAsync, 1e-10, 100);

            if (double.IsNaN(norm) || norm > 1e-3)
            {
                Console.WriteLine("Reset initial guess and recalculate.");
                SetInitialGuess(false);
                goto beginningOfSolve;
            }

            //  ██╗ 
            //  ╚██╗ write solution of Newtons method back to mesh
            //  ██╔╝
            //  ╚═╝
            if (simulationSelector == SimulationSelector.bothPotentials) // if both potentials are simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                {
                    mesh.finiteElements[i].phiFront = initialGuess[2 * i];
                    mesh.finiteElements[i].phiBack = initialGuess[2 * i + 1];
                }
            }
            else if (simulationSelector == SimulationSelector.frontPotential) // if only front potential is simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                {
                    mesh.finiteElements[i].phiFront = initialGuess[i];
                }
            }
            else // if only back potential is simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    mesh.finiteElements[i].phiBack = initialGuess[i];
            }

            //  ██╗ 
            //  ╚██╗ calculate currents and output power
            //  ██╔╝
            //  ╚═╝
            foreach (FiniteElementCell point in mesh.finiteElements.Values)
                point.CalculateDerivedParameters(mesh);

            var results = CalculateOutputData();
            cellCharacteristic.AddExperimentalPoint((results.voltage / amountOfCellsInRow, results.current / strechAlongP1, results.power, results.area, results.efficiency));
            lossAnalyses.Add(CalculateLosses());
            Misc.WriteFormatedLine("Cell efficiency at " + results.voltage + "V is " + results.efficiency + "%.");

            // calculate and return output power
            return results.power;
        }

        // Refine the mesh ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// refine mesh, calculate new voronoi mesh, solve differential equation again
        /// </summary>
        /// <param name="maxAmountOfRefinements">maximum amount of refinement iterations</param>
        public void RefineMesh(bool refineMesh, int maxAmountOfRefinements, double[] maximumDifferences,
            SimulationSelector simulationSelector, OpticMode opticMode, double lowerOperatingVoltage, double[] simulationVoltages)
        {
            if (!refineMesh)
                return;

            // run, till mesh is fine enough or maximum interations are done
            for (int i = 0; i < maxAmountOfRefinements; i++)
            {
                // Write to Console
                Misc.WriteFormatedLine();
                Console.WriteLine("Started refinement step " + (i + 1) + ".");

                // break when mesh is not refined anymore
                //mesh = delaunayVoronoi.RefineMesh(mesh, new Func<FiniteElementCell, double>[] { p => p.phiFront }, maximumDifferences, true,
                //out int amountBefore, out int amountAfter);
                //if (amountBefore == amountAfter)
                {
                    Console.WriteLine("Mesh didn't change.");
                    Misc.WriteFormatedLine();
                    Misc.WriteDividingLine();
                    break;
                }

                // Write timer to Console
                Console.WriteLine("Delaunay triangulation and Voronoi mesh of refinement step " + (i + 1) + " done.");

                //SetOpticsOfGlobalCell(opticMode, MiscTMM.spectrumAM15, illuminationIntensity);
                SetElectricsOfGlobalCell(simulationSelector, lowerOperatingVoltage);
                SetPreferencesOfSingleMeshpoints();
                SetInitialGuess();

                // solve differential equation via Newtons method
                Solve(out var simulationResults, voltageSweepMode, simulationVoltages);
            }
        }

        // Output to file ███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Write all points to file
        /// </summary>
        /// <param name="filepath">filepath of the outputfile</param>
        public void OutputToFileSpacialResolvedCell(string filepath)
        {
            using (StreamWriter file = new StreamWriter(filepath, false))
            {
                // Header
                file.WriteLine("x\ty\tindex\tPhi_{front}\tPhi_{back}\tfrontGrid?\tbackGrid?\tarea\tangle{I{out, lateral, front}}\t|I{out, lateral, front}|\tangle{I{out, lateral, back}}\t|I{out, lateral, back}|\tcreated heat-power");
                file.WriteLine("m\tm\t\tV\tV\t\t\tm^2\trad\tA\trad\tA\tW");

                // output points
                foreach (FiniteElementCell point in mesh.finiteElements.Values)
                    point.OutputToFile(file);

                file.Close();
            }
        }
        /// <summary>
        /// Write the characeristic curve to the given filepath
        /// </summary>
        /// <param name="filepath">filepath of the outputfile</param>
        /// <param name="paramNameList">list of all set parameter names</param>
        /// <param name="paramUnitList">list of all set parameter units</param>
        /// <param name="paramValueArray">list of all set parameter values</param>
        public void OutputToFileTotalCharacteristicCurve(string filepath, List<string> paramNameList = null, List<string> paramUnitList = null, double[] paramValueArray = null)
        {
            string filepathCharacteristics;
            if (paramNameList != null)
            {
                filepathCharacteristics = $"{ filepath.Remove(filepath.Length - 4) }";
                for (int i = 0; i < paramNameList.Count; i++)
                    filepathCharacteristics += "_" + paramNameList[i] + "=" + paramValueArray[i] + paramUnitList[i];
                filepathCharacteristics += ".dat";
            }
            else
                filepathCharacteristics = filepath;

            using (StreamWriter file = new StreamWriter(filepathCharacteristics, false))
            {
                // Header
                // properties
                file.Write("V\tI\tP\tPCE");

                file.Write("\t" + "PCE semiconductor");
                for (int i = lossAnalyses.First().lossMechanisms.Count - 1; i >= 0; i--)
                    file.Write("\t" + lossAnalyses.First().lossMechanisms[i].name);
                file.WriteLine();

                // units
                file.Write("V\tA\tW\t%");

                for (int i = lossAnalyses.First().lossMechanisms.Count; i >= 0; i--)
                    file.Write("\t" + "abs%");
                file.WriteLine();

                // output voltage points
                for (int i = 0; i < cellCharacteristic.experimentalData.Count; i++)
                {
                    // basic data
                    string datastring = InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[i].voltage * cellCharacteristic.factorToVoltage);
                    datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[i].current * cellCharacteristic.factorToCurrent);
                    datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[i].power);
                    datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[i].efficiency);

                    // losses
                    datastring += "\t" + InputOutput.ToStringWithSeparator(lossAnalyses[i].absoluteEfficiencyTheoretical);
                    for (int k = lossAnalyses.First().lossMechanisms.Count - 1; k >= 0; k--)
                        datastring += "\t" + InputOutput.ToStringWithSeparator(lossAnalyses[i].lossMechanisms[k].ratioAbsoluteLoss);

                    file.WriteLine(datastring);
                }
                file.Close();
            }
        }
        /// <summary>
        /// Write most efficient simulated voltage to sweep-file (append)
        /// </summary>
        /// <param name="filepath">filepath of the outputfile</param>
        /// <param name="paramValueArray">list of all set parameter values</param>
        public void OutputToFileOnlyMPP(string filepath, List<string> paramNameList, List<string> paramUnitList, double[] paramValueArray, bool createNewFile = false)
        {
            // open file
            StreamWriter fileSweep = new StreamWriter(filepath, !createNewFile);

            // Header
            if (createNewFile)
            {
                // properties
                foreach (var paramName in paramNameList)
                    fileSweep.Write(paramName + "\t");
                fileSweep.Write("V\tI\tP\tPCE\tVoc\tIsc\tVmpp\tImpp\tPmpp\tFF\tIph\tI0\tn\tRs\tRsh");

                fileSweep.Write("\t" + "PCE semiconductor");
                for (int i = lossAnalyses.First().lossMechanisms.Count - 1; i >= 0; i--)
                    fileSweep.Write("\t" + lossAnalyses.First().lossMechanisms[i].name);
                fileSweep.WriteLine();

                // units
                foreach (var paramUnit in paramUnitList)
                    fileSweep.Write(paramUnit + "\t");
                fileSweep.Write("V\tA\tW\t%\tV\tA\tV\tA\tW\t%\tA\tA\t\tOhm\tOhm");

                for (int i = lossAnalyses.First().lossMechanisms.Count; i >= 0; i--)
                    fileSweep.Write("\t" + "abs%");
                fileSweep.WriteLine();
            }

            // write parameters
            foreach (var paramValue in paramValueArray)
                fileSweep.Write(InputOutput.ToStringWithSeparator(paramValue) + "\t");

            // basic data
            int indexMPP = cellCharacteristic.experimentalData.Select((d, i) => new { Power = d.power, Index = i }).MinBy(x => x.Power).First().Index;
            string datastring = InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[indexMPP].voltage * cellCharacteristic.factorToVoltage);
            datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[indexMPP].current * cellCharacteristic.factorToCurrent);
            datastring += "\t" + InputOutput.ToStringWithSeparator(-cellCharacteristic.experimentalData[indexMPP].power);
            datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.experimentalData[indexMPP].efficiency);

            // fitting and more data
            if (voltageSweepMode == VoltageSweepMode.autoIVcurve || voltageSweepMode == VoltageSweepMode.autoIVcurveAndMPP
                || voltageSweepMode == VoltageSweepMode.multipleVoltages)
            {
                cellCharacteristic.ExecuteFit();

                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.GetDataSetOpenCircuit().voltage); // Voc in V
                datastring += "\t" + InputOutput.ToStringWithSeparator(-cellCharacteristic.GetDataSetShortCircuit().current); // Isc in A
                var MPP = cellCharacteristic.GetDataSetMaximumPowerPoint();
                datastring += "\t" + InputOutput.ToStringWithSeparator(MPP.voltage); // Vmpp in V
                datastring += "\t" + InputOutput.ToStringWithSeparator(MPP.current); // Impp in A
                datastring += "\t" + InputOutput.ToStringWithSeparator(MPP.power); // Pmpp in W
                datastring += "\t" + InputOutput.ToStringWithSeparator(MPP.fillfactor); // Fillfactor in %
                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.currentPhoto * cellCharacteristic.factorToCurrent); // Iph in A
                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.currentSaturation * cellCharacteristic.factorToCurrent); // I0 in A
                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.diode1IdealityFactor * cellCharacteristic.factorToVoltage); // n
                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.Rseries / cellCharacteristic.factorToCurrent * cellCharacteristic.factorToVoltage); // Rs in Ohm
                datastring += "\t" + InputOutput.ToStringWithSeparator(cellCharacteristic.Rshunt / cellCharacteristic.factorToCurrent * cellCharacteristic.factorToVoltage); // Rsh in Ohm
            }
            else
                for (int i = 0; i < 11; i++)
                    datastring += "\t";

            // losses
            datastring += "\t" + InputOutput.ToStringWithSeparator(lossAnalyses[indexMPP].absoluteEfficiencyTheoretical);
            for (int i = lossAnalyses[indexMPP].lossMechanisms.Count - 1; i >= 0; i--)
                datastring += "\t" + InputOutput.ToStringWithSeparator(lossAnalyses[indexMPP].lossMechanisms[i].ratioAbsoluteLoss);

            fileSweep.WriteLine(datastring);
            fileSweep.Close();
        }
        /// <summary>
        /// calculates the resulting characteristics of this cell
        /// </summary>
        /// <returns></returns>
        public (double voltage, double current, double power, double area, double efficiency) CalculateOutputData()
        {
            double current;
            if (simulationSelector == SimulationSelector.frontPotential)
                current = mesh.finiteElements.Values.Where(p => p.isExternalCellFrontContact).Select(p => p.GetCurrentGenerated()
                    + p.GetCurrentToNeighborsFrontSum(mesh)).Sum();
            else
                current = mesh.finiteElements.Values.Where(p => p.isExternalCellBackContact).Select(p => p.GetCurrentGenerated()
                    - p.GetCurrentToNeighborsBackSum(mesh)).Sum();

            return ((operatingVoltage.upper - operatingVoltage.lower) * amountOfCellsInRow, // voltage in [V]
                current * strechAlongP1, // current in [A]
                (operatingVoltage.upper - operatingVoltage.lower) * current * amountOfCellsInRow * strechAlongP1, // power in [W]
                totalArea * amountOfCellsInRow * strechAlongP1, // area [m²]
                -(operatingVoltage.upper - operatingVoltage.lower) * current / totalArea / totalPowerFromSpectrum * 100); // efficiency in [W / m² / W/m² * 100 = %]
        }

        // Calculate Losses █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// calculates all losses in this cell or module
        /// </summary>
        public (double powerTheoretical, double powerReal, double absoluteEfficiencyTheoretical, double absoluteEfficiencyReal, List<LossMechanism> lossMechanisms) CalculateLosses()
        {
            var cellData = CalculateOutputData();
            double powerCell_Vop = -cellData.power;
            double absoluteEfficiencyCell = cellData.efficiency;

            double powerSemiconductor_Vmpp = -semiconductorCharacteristic.GetDataSetMaximumPowerPoint().power;
            double powerDensitySemiconductor = powerSemiconductor_Vmpp / (totalArea * amountOfCellsInRow * strechAlongP1);
            double absoluteEfficiencySemiconductor = powerSemiconductor_Vmpp / (totalArea * amountOfCellsInRow * strechAlongP1 * totalPowerFromSpectrum) * 100;

            List<LossMechanism> lossMechanisms = new List<LossMechanism>();

            //  ██╗ 
            //  ╚██╗ edge area
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("edge area", powerSemiconductor_Vmpp,
                powerSemiconductor_Vmpp * edgeAreaPerSimulatedFraction / totalArea,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ additional atmospheric scattering
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("additional atmospheric scattering", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea)
                .Sum(p => p.size * p.opticFactor_additionalAtmosphereScattering) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗
            //  ╚██╗ shading
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("shading", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea)
                .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * p.opticFactor_shading) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ module area loss
            //  ██╔╝
            //  ╚═╝
            if (mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                lossMechanisms.Add(new LossMechanism("module P1-P2-P3-area loss", lossMechanisms.Last().powerAfterLoss,
                    powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type != pointType.cell)
                    .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * (1 - p.opticFactor_shading)) * amountOfCellsInRow * strechAlongP1,
                    powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ grid shading
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("grid shading", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type == pointType.cell)
                .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * (1 - p.opticFactor_shading) * (1 - p.opticFactor_transmissionGrid)) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ effective tilted area
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("effective tilted area", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type == pointType.cell)
                .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * (1 - p.opticFactor_shading) * p.opticFactor_transmissionGrid * (1 - p.opticFactor_effectiveArea)) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ reflection
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("reflection", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type == pointType.cell)
                .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * (1 - p.opticFactor_shading) * p.opticFactor_transmissionGrid * p.opticFactor_effectiveArea * p.opticFactor_reflection) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ parasitic absorption
            //  ██╔╝
            //  ╚═╝
            Dictionary<int, (string name, double loss)> parasiticLosses = new Dictionary<int, (string name, double loss)> ();
            foreach (var fe in mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type == pointType.cell))
                foreach (var pa in fe.opticFactor_absorption.Where(p => !p.isAbsorber))
                {
                    // if material ID is not yet assigned in dictionary, add it
                    if (!parasiticLosses.ContainsKey(pa.materialID))
                        parasiticLosses.Add(pa.materialID, (pa.materialName, 0));

                    // add loss value
                    parasiticLosses[pa.materialID] = (parasiticLosses[pa.materialID].name, parasiticLosses[pa.materialID].loss
                        + powerDensitySemiconductor * fe.size * (1 - fe.opticFactor_additionalAtmosphereScattering) * (1 - fe.opticFactor_shading)
                        * fe.opticFactor_transmissionGrid * fe.opticFactor_effectiveArea * pa.factor * amountOfCellsInRow * strechAlongP1);
                }

            // add each loss to analysis
            foreach (var pa in parasiticLosses)
                lossMechanisms.Add(new LossMechanism("parasitic absorption in " + pa.Value.name, lossMechanisms.Last().powerAfterLoss,
                        pa.Value.loss, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ incomplete absorption / transmission
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("transmission through device", lossMechanisms.Last().powerAfterLoss,
                powerDensitySemiconductor * mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.type == pointType.cell)
                .Sum(p => p.size * (1 - p.opticFactor_additionalAtmosphereScattering) * (1 - p.opticFactor_shading) * p.opticFactor_transmissionGrid * p.opticFactor_effectiveArea * p.opticFactor_transmission) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ local MPP missmatch
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("local MPP mismatch", lossMechanisms.Last().powerAfterLoss,
                lossMechanisms.Last().powerAfterLoss - mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.frontGrid == null)
                .Sum(p => -p.GetCurrentGenerated() * (p.phiFront - p.phiBack)) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ reverse current under grid
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("reverse current under grid", lossMechanisms.Last().powerAfterLoss,
                mesh.finiteElements.Values.Where(p => p.countsAsActiveArea).Where(p => p.frontGrid != null)
                .Sum(p => p.GetCurrentGenerated() * (p.phiFront - p.phiBack)) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ reverse current in non-active area
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("reverse current in non-active area", lossMechanisms.Last().powerAfterLoss,
                mesh.finiteElements.Values.Where(p => !p.countsAsActiveArea)
                .Sum(p => p.GetCurrentGenerated() * (p.phiFront - p.phiBack)) * amountOfCellsInRow * strechAlongP1,
                powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ ohmic losses in contact and grid layers in the cell
            //  ██╔╝
            //  ╚═╝
            double ohmicLossFrontContact = 0, ohmicLossFrontGrid = 0, ohmicLossBackContact = 0, ohmicLossBackGrid = 0;
            foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.cell))
                for (int n_inPoint = 0; n_inPoint < point.neighbors.Count; n_inPoint++) //  Iterate over each point-pair (in both directions)
                {
                    FiniteElementCell neighbor = mesh.finiteElements[point.neighbors[n_inPoint].index];
                    int p_inNeighbor = neighbor.neighbors.FindIndex(p => p.index == point.index);

                    // Front side
                    if (!double.IsNaN(point.resistorsFrontToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsFrontToEdge[p_inNeighbor]))
                    {
                        double potentialdifferenceOnOwnFrontSide = (neighbor.phiFront - point.phiFront) * point.resistorsFrontToEdge[n_inPoint]
                          / (point.resistorsFrontToEdge[n_inPoint] + neighbor.resistorsFrontToEdge[p_inNeighbor]);
                        if (!double.IsNaN(point.resistorsFrontToEdgeContact[n_inPoint])
                            && !double.IsInfinity(point.resistorsFrontToEdgeContact[n_inPoint]) && point.resistorsFrontToEdgeContact[n_inPoint] != 0)
                            ohmicLossFrontContact += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeContact[n_inPoint];
                        if (point.frontGrid != null)
                            if (!double.IsNaN(point.resistorsFrontToEdgeGrid[n_inPoint])
                                && !double.IsInfinity(point.resistorsFrontToEdgeGrid[n_inPoint]) && point.resistorsFrontToEdgeGrid[n_inPoint] != 0)
                                ohmicLossFrontGrid += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeGrid[n_inPoint];
                    }

                    // Back side
                    if (!double.IsNaN(point.resistorsBackToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsBackToEdge[p_inNeighbor]))
                    {
                        double potentialdifferenceOnOwnBackSide = (neighbor.phiBack - point.phiBack) * point.resistorsBackToEdge[n_inPoint]
                          / (point.resistorsBackToEdge[n_inPoint] + neighbor.resistorsBackToEdge[p_inNeighbor]);
                        if (!double.IsNaN(point.resistorsBackToEdgeContact[n_inPoint])
                            && !double.IsInfinity(point.resistorsBackToEdgeContact[n_inPoint]) && point.resistorsBackToEdgeContact[n_inPoint] != 0)
                            ohmicLossBackContact += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeContact[n_inPoint];
                        if (point.backGrid != null)
                            if (!double.IsNaN(point.resistorsBackToEdgeGrid[n_inPoint])
                                && !double.IsInfinity(point.resistorsBackToEdgeGrid[n_inPoint]) && point.resistorsBackToEdgeGrid[n_inPoint] != 0)
                                ohmicLossBackGrid += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeGrid[n_inPoint];
                    }
                }
            lossMechanisms.Add(new LossMechanism("ohmic front contact layer", lossMechanisms.Last().powerAfterLoss,
                ohmicLossFrontContact * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            lossMechanisms.Add(new LossMechanism("ohmic front grid layer", lossMechanisms.Last().powerAfterLoss,
                ohmicLossFrontGrid * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            lossMechanisms.Add(new LossMechanism("ohmic back contact layer", lossMechanisms.Last().powerAfterLoss,
                ohmicLossBackContact * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            lossMechanisms.Add(new LossMechanism("ohmic back grid layer", lossMechanisms.Last().powerAfterLoss,
                ohmicLossBackGrid * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            //  ██╗ 
            //  ╚██╗ ohmic losses in module transition
            //  ██╔╝
            //  ╚═╝
            if (mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                double ohmicLossModuleTransition = 0;

                // horizontal current in P1, gap12 and P2
                foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.P1 || p.type == pointType.gap12 || p.type == pointType.P2))
                    for (int n_inPoint = 0; n_inPoint < point.neighbors.Count; n_inPoint++) //  Iterate over each point-pair (in both directions)
                    {
                        FiniteElementCell neighbor = mesh.finiteElements[point.neighbors[n_inPoint].index];
                        int p_inNeighbor = neighbor.neighbors.FindIndex(p => p.index == point.index);

                        // Front side
                        if (!double.IsNaN(point.resistorsFrontToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsFrontToEdge[p_inNeighbor]))
                        {
                            double potentialdifferenceOnOwnFrontSide = (neighbor.phiFront - point.phiFront) * point.resistorsFrontToEdge[n_inPoint]
                              / (point.resistorsFrontToEdge[n_inPoint] + neighbor.resistorsFrontToEdge[p_inNeighbor]);
                            if (!double.IsNaN(point.resistorsFrontToEdgeContact[n_inPoint])
                                && !double.IsInfinity(point.resistorsFrontToEdgeContact[n_inPoint]) && point.resistorsFrontToEdgeContact[n_inPoint] != 0)
                                ohmicLossModuleTransition += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeContact[n_inPoint];
                            if (point.frontGrid != null)
                                if (!double.IsNaN(point.resistorsFrontToEdgeGrid[n_inPoint])
                                    && !double.IsInfinity(point.resistorsFrontToEdgeGrid[n_inPoint]) && point.resistorsFrontToEdgeGrid[n_inPoint] != 0)
                                    ohmicLossModuleTransition += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeGrid[n_inPoint];
                        }
                    }

                // horizontal current in P2, gap23 and P3
                foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.P2 || p.type == pointType.gap23 || p.type == pointType.P3))
                    for (int n_inPoint = 0; n_inPoint < point.neighbors.Count; n_inPoint++) //  Iterate over each point-pair (in both directions)
                    {
                        FiniteElementCell neighbor = mesh.finiteElements[point.neighbors[n_inPoint].index];
                        int p_inNeighbor = neighbor.neighbors.FindIndex(p => p.index == point.index);

                        // Back side
                        if (!double.IsNaN(point.resistorsBackToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsBackToEdge[p_inNeighbor]))
                        {
                            double potentialdifferenceOnOwnBackSide = (neighbor.phiBack - point.phiBack) * point.resistorsBackToEdge[n_inPoint]
                              / (point.resistorsBackToEdge[n_inPoint] + neighbor.resistorsBackToEdge[p_inNeighbor]);
                            if (!double.IsNaN(point.resistorsBackToEdgeContact[n_inPoint])
                                && !double.IsInfinity(point.resistorsBackToEdgeContact[n_inPoint]) && point.resistorsBackToEdgeContact[n_inPoint] != 0)
                                ohmicLossModuleTransition += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeContact[n_inPoint];
                            if (point.backGrid != null)
                                if (!double.IsNaN(point.resistorsBackToEdgeGrid[n_inPoint])
                                    && !double.IsInfinity(point.resistorsBackToEdgeGrid[n_inPoint]) && point.resistorsBackToEdgeGrid[n_inPoint] != 0)
                                    ohmicLossModuleTransition += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeGrid[n_inPoint];
                        }
                    }

                // vertical current in P2
                foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.P2))
                    ohmicLossModuleTransition += point.frontContact.propertiesContact.GetSpecificResistivity(point.thicknessFrontContact) * point.pnJunction.thicknessAbsorberLayer / point.size * Math.Pow(point.GetCurrentP2(), 2);

                lossMechanisms.Add(new LossMechanism("ohmic path-losses in module transition", lossMechanisms.Last().powerAfterLoss,
                    ohmicLossModuleTransition * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ losses in P2 contact resistance
            //  ██╔╝
            //  ╚═╝
            if (mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                double ohmicP2ContactResistance = 0;
                foreach (var P2point in mesh.finiteElements.Values.Where(p => p.type == pointType.P2))
                    ohmicP2ContactResistance += P2point.contactSeriesResistanceP2 * Math.Pow(P2point.GetCurrentP2(), 2);
                lossMechanisms.Add(new LossMechanism("P2 contact resistance", lossMechanisms.Last().powerAfterLoss,
                    ohmicP2ContactResistance * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ shunt losses in P1
            //  ██╔╝
            //  ╚═╝
            if (mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                double ohmicLossP1 = 0;
                // horizontal current in P1 and gap12
                foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.P1 || p.type == pointType.gap12))
                    for (int n_inPoint = 0; n_inPoint < point.neighbors.Count; n_inPoint++) //  Iterate over each point-pair (in both directions)
                    {
                        FiniteElementCell neighbor = mesh.finiteElements[point.neighbors[n_inPoint].index];
                        int p_inNeighbor = neighbor.neighbors.FindIndex(p => p.index == point.index);

                        // Back side
                        if (!double.IsNaN(point.resistorsBackToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsBackToEdge[p_inNeighbor]))
                        {
                            double potentialdifferenceOnOwnBackSide = (neighbor.phiBack - point.phiBack) * point.resistorsBackToEdge[n_inPoint]
                              / (point.resistorsBackToEdge[n_inPoint] + neighbor.resistorsBackToEdge[p_inNeighbor]);
                            if (!double.IsNaN(point.resistorsBackToEdgeContact[n_inPoint])
                                && !double.IsInfinity(point.resistorsBackToEdgeContact[n_inPoint]) && point.resistorsBackToEdgeContact[n_inPoint] != 0)
                                ohmicLossP1 += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeContact[n_inPoint];
                            if (point.backGrid != null)
                                if (!double.IsNaN(point.resistorsBackToEdgeGrid[n_inPoint])
                                    && !double.IsInfinity(point.resistorsBackToEdgeGrid[n_inPoint]) && point.resistorsBackToEdgeGrid[n_inPoint] != 0)
                                    ohmicLossP1 += Math.Pow(potentialdifferenceOnOwnBackSide, 2) / point.resistorsBackToEdgeGrid[n_inPoint];
                        }
                    }

                lossMechanisms.Add(new LossMechanism("shunt over P1", lossMechanisms.Last().powerAfterLoss,
                    ohmicLossP1 * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ shunt losses in P3
            //  ██╔╝
            //  ╚═╝
            if (mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                double ohmicLossP3 = 0;
                // horizontal current in gap23 and P3
                foreach (var point in mesh.finiteElements.Values.Where(p => p.type == pointType.gap23 || p.type == pointType.P3))
                    for (int n_inPoint = 0; n_inPoint < point.neighbors.Count; n_inPoint++) //  Iterate over each point-pair (in both directions)
                    {
                        FiniteElementCell neighbor = mesh.finiteElements[point.neighbors[n_inPoint].index];
                        int p_inNeighbor = neighbor.neighbors.FindIndex(p => p.index == point.index);

                        // Front side
                        if (!double.IsNaN(point.resistorsFrontToEdge[n_inPoint]) && !double.IsNaN(neighbor.resistorsFrontToEdge[p_inNeighbor]))
                        {
                            double potentialdifferenceOnOwnFrontSide = (neighbor.phiFront - point.phiFront) * point.resistorsFrontToEdge[n_inPoint]
                              / (point.resistorsFrontToEdge[n_inPoint] + neighbor.resistorsFrontToEdge[p_inNeighbor]);
                            if (!double.IsNaN(point.resistorsFrontToEdgeContact[n_inPoint])
                                && !double.IsInfinity(point.resistorsFrontToEdgeContact[n_inPoint]) && point.resistorsFrontToEdgeContact[n_inPoint] != 0)
                                ohmicLossP3 += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeContact[n_inPoint];
                            if (point.frontGrid != null)
                                if (!double.IsNaN(point.resistorsFrontToEdgeGrid[n_inPoint])
                                    && !double.IsInfinity(point.resistorsFrontToEdgeGrid[n_inPoint]) && point.resistorsFrontToEdgeGrid[n_inPoint] != 0)
                                    ohmicLossP3 += Math.Pow(potentialdifferenceOnOwnFrontSide, 2) / point.resistorsFrontToEdgeGrid[n_inPoint];
                        }
                    }

                lossMechanisms.Add(new LossMechanism("shunt over P3", lossMechanisms.Last().powerAfterLoss,
                    ohmicLossP3 * amountOfCellsInRow * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ ohmic losses in external cell contact resistances
            //  ██╔╝
            //  ╚═╝
            if (!mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                double ohmicContactResistanceLosses = 0;

                foreach (var externalFrontContactPoint in mesh.finiteElements.Values.Where(p => p.isExternalCellFrontContact))
                    ohmicContactResistanceLosses += externalFrontContactPoint.contactResistanceExternalCellFrontContact == 0 ? 0 :
                        (Math.Pow(externalFrontContactPoint.phiFront - operatingVoltage.upper, 2)
                        / externalFrontContactPoint.contactResistanceExternalCellFrontContact);

                foreach (var externalBackContactPoint in mesh.finiteElements.Values.Where(p => p.isExternalCellBackContact))
                    ohmicContactResistanceLosses += externalBackContactPoint.contactResistanceExternalCellBackContact == 0 ? 0 :
                        (Math.Pow(externalBackContactPoint.phiBack - operatingVoltage.lower, 2)
                        / externalBackContactPoint.contactResistanceExternalCellBackContact);

                lossMechanisms.Add(new LossMechanism("external cell contact resistances", lossMechanisms.Last().powerAfterLoss,
                    ohmicContactResistanceLosses * strechAlongP1, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));
            }

            //  ██╗ 
            //  ╚██╗ unadressable
            //  ██╔╝
            //  ╚═╝
            lossMechanisms.Add(new LossMechanism("unadressable", lossMechanisms.Last().powerAfterLoss,
                lossMechanisms.Last().powerAfterLoss - powerCell_Vop, powerSemiconductor_Vmpp, absoluteEfficiencySemiconductor));

            return (powerSemiconductor_Vmpp, powerCell_Vop, absoluteEfficiencySemiconductor, absoluteEfficiencyCell, lossMechanisms);
        }

        // get residuum of the function █████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the residuum of the function
        /// </summary>
        /// <param name="solution">array with all electrical potentials</param>
        /// <returns></returns>
        public Vector<double> GetFunctionAsync(Vector<double> solution)
        {
            // set solution vector of previous iteration to mesh
            if (simulationSelector == SimulationSelector.bothPotentials) // if both potentials are simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                {
                    mesh.finiteElements[i].phiFront = solution[2 * i];
                    mesh.finiteElements[i].phiBack = solution[2 * i + 1];
                }
            }
            else if (simulationSelector == SimulationSelector.frontPotential) // if only front potential is simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    mesh.finiteElements[i].phiFront = solution[i];
            }
            else // if only back potential is simulated
            {
                for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
                    mesh.finiteElements[i].phiBack = solution[i];
            }

            // Get the amount of logical processor units and create a list of as many tasks as processors
            int amountProcessors = Environment.ProcessorCount;
            Task<double[]>[] tasks = new Task<double[]>[amountProcessors];

            // Get the amount of rows, one single processor will calculate
            // (always rounded down and the last processor will have to do the modulo-rest on top)
            int lengthOfBlock = mesh.nextAvailableFiniteElementIndex / amountProcessors;

            // Create and start all tasks
            for (int p = 0; p < amountProcessors - 1; p++)
            {
                int processorNumber = p;
                tasks[p] = Task.Run(() => GetFunctionBlock(processorNumber * lengthOfBlock, (processorNumber + 1) * lengthOfBlock));
            }
            tasks[amountProcessors - 1] = Task.Run(() => GetFunctionBlock((amountProcessors - 1) * lengthOfBlock, mesh.nextAvailableFiniteElementIndex));

            // Wait for all tasks to be finished => all blocks are set
            double[][] resultListArray = Task.WhenAll(tasks).Result;

            // Create and return vector
            var F = Vector.Create(Misc.ConcatArrays(resultListArray));
            return F;
        }
        /// <summary>
        /// builds up one block of the Function vector (there are as many blocks as processors)
        /// </summary>
        /// <param name="startPointIndexIncluding">frist pointindex, which will be included in this block</param>
        /// <param name="endPointIndexExcluding">frist pointindex, which WON'T be included anymore in this block</param>
        /// <returns></returns>
        private double[] GetFunctionBlock(int startPointIndexIncluding, int endPointIndexExcluding)
        {
            //  ██╗ 
            //  ╚██╗ both potentials are simulated
            //  ██╔╝
            //  ╚═╝
            #region both potentials are simulated
            if (simulationSelector == SimulationSelector.bothPotentials)
            {
                // structure of F-vector:
                //
                //                  / F_1_front \
                //                 |  F_1_back   |
                //                 |  F_2_front  |
                //           F  =  |  F_2_back   |
                //                 |     ...     |
                //                 |  F_n_front  |
                //                  \ F_n_back  /

                // Initialize array, which will be used for creating the vector
                double[] values = new double[2 * (endPointIndexExcluding - startPointIndexIncluding)];

                // Fill array
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at front contact ————————————————————————————————————————————————————————————————————————————————————————————————————————
                    if (mesh.finiteElements[k].hasModuleConnection_isInP3)
                    {
                        // element, which is connected to the back side of the next cell
                        //
                        //           F_k  =  Φ_front - Φ_front_otherSeg - Vop
                        //
                        values[2 * (k - startPointIndexIncluding)]
                            = mesh.finiteElements[k].phiFront
                            - mesh.finiteElements[mesh.finiteElements[k].indexOfModuleConnectedPoint].phiFront
                            - (operatingVoltage.upper - operatingVoltage.lower);
                    }
                    else
                    {
                        // "normal" element
                        //
                        //           F_k  =  ∑ I_k,n  +  I_gen or I_P2
                        //
                        values[2 * (k - startPointIndexIncluding)]
                            = mesh.finiteElements[k].GetCurrentToNeighborsFrontSum(mesh)
                            + (mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2() : mesh.finiteElements[k].GetCurrentGenerated());

                        // external cell front contact element:
                        //
                        //           F_ecfc  =  F_k * Rcont_front  +  Φ_front - Vop_upper
                        //
                        if (mesh.finiteElements[k].isExternalCellFrontContact)
                            values[2 * (k - startPointIndexIncluding)] = values[2 * (k - startPointIndexIncluding)]
                                * mesh.finiteElements[k].contactResistanceExternalCellFrontContact
                                + mesh.finiteElements[k].phiFront - operatingVoltage.upper;
                    }

                    // residual at back contact —————————————————————————————————————————————————————————————————————————————————————————————————————————
                    if (mesh.finiteElements[k].hasModuleConnection_isInP3)
                    {
                        // element, which is connected to the back side of the next cell
                        //
                        //           F_k  =  Φ_back - Φ_back_otherSeg - Vop
                        //
                        values[2 * (k - startPointIndexIncluding) + 1]
                            = mesh.finiteElements[k].phiBack
                            - mesh.finiteElements[mesh.finiteElements[k].indexOfModuleConnectedPoint].phiBack
                            - (operatingVoltage.upper - operatingVoltage.lower);
                    }
                    else
                    {
                        // "normal" element
                        //
                        //           F_k  =  ∑ I_k,n  -  I_gen or I_P2
                        //
                        values[2 * (k - startPointIndexIncluding) + 1]
                            = mesh.finiteElements[k].GetCurrentToNeighborsBackSum(mesh)
                            - (mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2() : mesh.finiteElements[k].GetCurrentGenerated());

                        // external cell back contact element
                        //
                        //           F_ecbc  =  F_k * Rcont_back  +  Φ_back - Vop_lower
                        //
                        if (mesh.finiteElements[k].isExternalCellBackContact)
                            values[2 * (k - startPointIndexIncluding) + 1] = values[2 * (k - startPointIndexIncluding) + 1]
                                * mesh.finiteElements[k].contactResistanceExternalCellBackContact
                                + mesh.finiteElements[k].phiBack - operatingVoltage.lower;
                    }
                }

                return values;
            }
            #endregion

            //  ██╗ 
            //  ╚██╗ only front potential is simulated
            //  ██╔╝
            //  ╚═╝
            #region only front potential is simulated
            else if (simulationSelector == SimulationSelector.frontPotential)
            {
                // structure of F-vector:
                //
                //                  / F_1_front \
                //                 |  F_2_front  |
                //           F  =  |     ...     |
                //                 |     ...     |
                //                  \ F_n_front /

                // Initialize array, which will be used for creating the vector
                double[] values = new double[endPointIndexExcluding - startPointIndexIncluding];

                // Fill array
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at front contact ————————————————————————————————————————————————————————————————————————————————————————————————————————

                    // "normal" element
                    //
                    //           F_k  =  ∑ I_k,n  +  I_gen
                    //
                    values[k - startPointIndexIncluding]
                        = mesh.finiteElements[k].GetCurrentToNeighborsFrontSum(mesh) + mesh.finiteElements[k].GetCurrentGenerated();

                    // external cell front contact element:
                    //
                    //           F_ecfc  =  F_k * Rcont_front  +  Φ_front - Vop_upper
                    //
                    if (mesh.finiteElements[k].isExternalCellFrontContact)
                        values[k - startPointIndexIncluding] = values[k - startPointIndexIncluding]
                            * mesh.finiteElements[k].contactResistanceExternalCellFrontContact + mesh.finiteElements[k].phiFront - operatingVoltage.upper;
                }

                return values;
            }
            #endregion

            //  ██╗ 
            //  ╚██╗ only back potential is simulated
            //  ██╔╝
            //  ╚═╝
            #region only back potential is simulated
            else
            {
                // structure of F-vector:
                //
                //                  / F_1_back \
                //                 |  F_2_back  |
                //           F  =  |     ...    |
                //                 |     ...    |
                //                  \ F_n_back /

                // Initialize array, which will be used for creating the vector
                double[] values = new double[endPointIndexExcluding - startPointIndexIncluding];

                // Fill array
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at back contact —————————————————————————————————————————————————————————————————————————————————————————————————————————

                    // "normal" element
                    //
                    //           F_k  =  ∑ I_k,n  -  I_gen
                    //
                    values[k - startPointIndexIncluding]
                        = mesh.finiteElements[k].GetCurrentToNeighborsBackSum(mesh) - mesh.finiteElements[k].GetCurrentGenerated();

                    // external cell back contact element
                    //
                    //           F_ecbc  =  F_k * Rcont_back  +  Φ_back - Vop_lower
                    //
                    if (mesh.finiteElements[k].isExternalCellBackContact)
                        values[k - startPointIndexIncluding] = values[k - startPointIndexIncluding] *
                            mesh.finiteElements[k].contactResistanceExternalCellBackContact + mesh.finiteElements[k].phiBack - operatingVoltage.lower;
                }

                return values;
            }
            #endregion
        }

        // get Jacobi matrix ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the Jacobi matrix of the function
        /// </summary>
        /// <param name="solution">array with all electrical potentials</param>
        /// <returns></returns>
        public SparseMatrix<double> GetJacobiAsync(Vector<double> solution)
        {
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
                tasks[p] = Task.Run(() => GetJacobiBlock(processorNumber * lengthOfBlock, (processorNumber + 1) * lengthOfBlock));
            }
            tasks[amountProcessors - 1] = Task.Run(() => GetJacobiBlock((amountProcessors - 1) * lengthOfBlock, mesh.nextAvailableFiniteElementIndex));

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
            int matrixSize = mesh.nextAvailableFiniteElementIndex;
            if (simulationSelector == SimulationSelector.bothPotentials)
                matrixSize *= 2;
            var J = Matrix.CreateSparse(matrixSize, matrixSize,
                Misc.ConcatArrays(rowIndexesSingle), Misc.ConcatArrays(columnIndexesSingle), Misc.ConcatArrays(valuesSingle));
            return J;
        }
        /// <summary>
        /// builds up one block of the Jacobi matrix (there are as many blocks as processors)
        /// </summary>
        /// <param name="startPointIndexIncluding">frist pointindex, which will be included in this block</param>
        /// <param name="endPointIndexExcluding">frist pointindex, which WON'T be included anymore in this block</param>
        /// <returns></returns>
        private (int[] rowIndexes, int[] columnIndexes, double[] values) GetJacobiBlock(int startPointIndexIncluding, int endPointIndexExcluding)
        {
            //  ██╗ 
            //  ╚██╗ both potentials are simulated
            //  ██╔╝
            //  ╚═╝
            #region both potentials are simulated
            if (simulationSelector == SimulationSelector.bothPotentials)
            {
                // structure of Jacobimatrix:
                //
                //                  /  ∂F_1f/∂Φ_1f  ∂F_1f/∂Φ_1b  · · ·  ∂Igen_1f/∂Φ_1f     0     \
                //                 |   ∂F_1b/∂Φ_1f  ∂F_1b/∂Φ_1b  · · ·       0     ∂Igen_1b/∂Φ_1b |
                //                 |                                                              |
                //                 |        ·            ·       ·           ·            ·       |
                //           J  =  |        ·            ·         ·         ·            ·       |
                //                 |        ·            ·           ·       ·            ·       |
                //                 |                                                              |
                //                 | ∂Igen_Nf/∂Φ_1f     0        · · ·  ∂F_Nf/∂Φ_Nf  ∂F_Nf/∂Φ_Nb  |
                //                  \     0     ∂Igen_Nb/∂Φ_1b   · · ·  ∂F_Nb/∂Φ_Nf  ∂F_Nb/∂Φ_Nb /

                // Get amount of non-zero elements (amount of points/rows plus sum of neighbors of all points)
                int amountOfNonZeroElements = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    if (mesh.finiteElements[k].hasModuleConnection_isInP3)
                        amountOfNonZeroElements += 4;
                    else
                        amountOfNonZeroElements += 4 + 2 * mesh.finiteElements[k].neighbors.Count;
                }

                // Initialize arrays, which will be used for creating the sparse matrix
                int[] rowIndexes = new int[amountOfNonZeroElements];
                int[] columnIndexes = new int[amountOfNonZeroElements];
                double[] values = new double[amountOfNonZeroElements];

                // Fill arrays at non-zero positions
                int indexInSparseArrays = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at front contact (if external cell contact: J_ij = J_ij * Rcontact + 1) —————————————————————————————————————————————————
                    if (mesh.finiteElements[k].hasModuleConnection_isInP3)
                    {
                        rowIndexes[indexInSparseArrays] = 2 * k;
                        columnIndexes[indexInSparseArrays] = 2 * k;
                        values[indexInSparseArrays++] = 1;

                        rowIndexes[indexInSparseArrays] = 2 * k;
                        columnIndexes[indexInSparseArrays] = 2 * mesh.finiteElements[k].indexOfModuleConnectedPoint;
                        values[indexInSparseArrays++] = -1;
                    }
                    else
                    {
                        // diagonal element
                        //                        ___
                        //            ∂F_k        \    I_k,n       I_generated / I_P2
                        //           ——————  =     ⟩   ———————  +  ———————————————————
                        //            ∂Φ_k        /     ∂Φ_k              ∂Φ_k
                        //                        ‾‾‾
                        //                      n ∈ N(k)

                        // ∂F_f/∂Φ_f
                        rowIndexes[indexInSparseArrays] = 2 * k;
                        columnIndexes[indexInSparseArrays] = 2 * k;
                        values[indexInSparseArrays++] = mesh.finiteElements[k].neighbors.Select(n =>
                            mesh.finiteElements[k].GetCurrentToNeighborsDerivativeFrontSingle(mesh, n.index)).Sum()
                            + (mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2_DerivativePhiFront() : mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiFront());
                        if (mesh.finiteElements[k].isExternalCellFrontContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellFrontContact + 1;

                        // ∂F_f/∂Φ_b
                        rowIndexes[indexInSparseArrays] = 2 * k;
                        columnIndexes[indexInSparseArrays] = 2 * k + 1;
                        values[indexInSparseArrays++] = (mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2_DerivativePhiBack() : mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiBack());
                        if (mesh.finiteElements[k].isExternalCellFrontContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellFrontContact;

                        // offdiagonal element
                        //
                        //            ∂F_k       ∂I_k,n
                        //           ——————  =  ————————
                        //            ∂Φ_n        ∂Φ_n
                        //
                        for (int n = 0; n < mesh.finiteElements[k].neighbors.Count; n++)
                        {
                            // ∂F_f/∂Φ_f
                            rowIndexes[indexInSparseArrays] = 2 * k;
                            columnIndexes[indexInSparseArrays] = 2 * mesh.finiteElements[k].neighbors[n].index;
                            values[indexInSparseArrays++] = -mesh.finiteElements[k].GetCurrentToNeighborsDerivativeFrontSingle(mesh, mesh.finiteElements[k].neighbors[n].index);
                            if (mesh.finiteElements[k].isExternalCellFrontContact)
                                values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellFrontContact;
                        }
                    }

                    // residual at back contact (if external cell contact: J_ij = J_ij * Rcontact + 1) ——————————————————————————————————————————————————
                    if (mesh.finiteElements[k].hasModuleConnection_isInP3)
                    {
                        rowIndexes[indexInSparseArrays] = 2 * k + 1;
                        columnIndexes[indexInSparseArrays] = 2 * k + 1;
                        values[indexInSparseArrays++] = 1;

                        rowIndexes[indexInSparseArrays] = 2 * k + 1;
                        columnIndexes[indexInSparseArrays] = 2 * mesh.finiteElements[k].indexOfModuleConnectedPoint + 1;
                        values[indexInSparseArrays++] = -1;
                    }
                    else
                    {
                        // diagonal element
                        //                        ___
                        //            ∂F_k        \    I_k,n       I_generated / I_P2
                        //           ——————  =     ⟩   ———————  -  ———————————————————
                        //            ∂Φ_k        /     ∂Φ_k              ∂Φ_k
                        //                        ‾‾‾
                        //                      n ∈ N(k)

                        // ∂F_b/∂Φ_b
                        rowIndexes[indexInSparseArrays] = 2 * k + 1;
                        columnIndexes[indexInSparseArrays] = 2 * k + 1;
                        values[indexInSparseArrays++] = mesh.finiteElements[k].neighbors.Select(n =>
                            mesh.finiteElements[k].GetCurrentToNeighborsDerivativeBackSingle(mesh, n.index)).Sum()
                            - (mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2_DerivativePhiBack() : mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiBack());
                        if (mesh.finiteElements[k].isExternalCellBackContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellBackContact + 1;

                        // ∂F_b/∂Φ_f
                        rowIndexes[indexInSparseArrays] = 2 * k + 1;
                        columnIndexes[indexInSparseArrays] = 2 * k;
                        values[indexInSparseArrays++] = -(mesh.finiteElements[k].type == pointType.P2 ? mesh.finiteElements[k].GetCurrentP2_DerivativePhiFront() : mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiFront());
                        if (mesh.finiteElements[k].isExternalCellBackContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellBackContact;

                        // offdiagonal element
                        //
                        //            ∂F_k       ∂I_k,n
                        //           ——————  =  ————————
                        //            ∂Φ_n        ∂Φ_n
                        //
                        for (int n = 0; n < mesh.finiteElements[k].neighbors.Count; n++)
                        {
                            // ∂F_b/∂Φ_b
                            rowIndexes[indexInSparseArrays] = 2 * k + 1;
                            columnIndexes[indexInSparseArrays] = 2 * mesh.finiteElements[k].neighbors[n].index + 1;
                            values[indexInSparseArrays++] = -mesh.finiteElements[k].GetCurrentToNeighborsDerivativeBackSingle(mesh, mesh.finiteElements[k].neighbors[n].index);
                            if (mesh.finiteElements[k].isExternalCellBackContact)
                                values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellBackContact;
                        }
                    }
                }

                return (rowIndexes, columnIndexes, values);
            }
            #endregion

            //  ██╗ 
            //  ╚██╗ only front potential is simulated
            //  ██╔╝
            //  ╚═╝
            #region only front potential is simulated
            else if (simulationSelector == SimulationSelector.frontPotential)
            {
                // structure of Jacobimatrix:
                //
                //                  /  ∂F_1f/∂Φ_1f   · · ·   ∂Igen_1f/∂Φ_1f  \
                //                 |                                          |
                //                 |        ·         ·             ·         |
                //           J  =  |        ·           ·           ·         |
                //                 |        ·             ·         ·         |
                //                 |                                          |
                //                  \  ∂Igen_Nf/∂Φ_1f   · · ·   ∂F_Nf/∂Φ_Nf  /

                // Get amount of non-zero elements (amount of points/rows plus sum of neighbors of all points)
                int amountOfNonZeroElements = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                    amountOfNonZeroElements += 1 + mesh.finiteElements[k].neighbors.Count;

                // Initialize arrays, which will be used for creating the sparse matrix
                int[] rowIndexes = new int[amountOfNonZeroElements];
                int[] columnIndexes = new int[amountOfNonZeroElements];
                double[] values = new double[amountOfNonZeroElements];

                // Fill arrays at non-zero positions
                int indexInSparseArrays = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at front contact (if external cell contact: J_ij = J_ij * Rcontact + 1) —————————————————————————————————————————————————

                    // diagonal element
                    //                        ___
                    //            ∂F_k        \    I_k,n       I_generated
                    //           ——————  =     ⟩   ———————  +  —————————————
                    //            ∂Φ_k        /     ∂Φ_k           ∂Φ_k
                    //                        ‾‾‾
                    //                      n ∈ N(k)

                    // ∂F_f/∂Φ_f
                    rowIndexes[indexInSparseArrays] = k;
                    columnIndexes[indexInSparseArrays] = k;
                    values[indexInSparseArrays++] = mesh.finiteElements[k].neighbors.Select(n =>
                        mesh.finiteElements[k].GetCurrentToNeighborsDerivativeFrontSingle(mesh, n.index)).Sum()
                        + mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiFront();
                    if (mesh.finiteElements[k].isExternalCellFrontContact)
                        values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellFrontContact + 1;

                    // offdiagonal element
                    //
                    //            ∂F_k       ∂I_k,n
                    //           ——————  =  ————————
                    //            ∂Φ_n        ∂Φ_n
                    //
                    for (int n = 0; n < mesh.finiteElements[k].neighbors.Count; n++)
                    {
                        // ∂F_f/∂Φ_f
                        rowIndexes[indexInSparseArrays] = k;
                        columnIndexes[indexInSparseArrays] = mesh.finiteElements[k].neighbors[n].index;
                        values[indexInSparseArrays++] = -mesh.finiteElements[k].GetCurrentToNeighborsDerivativeFrontSingle(mesh, mesh.finiteElements[k].neighbors[n].index);
                        if (mesh.finiteElements[k].isExternalCellFrontContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellFrontContact;
                    }
                }

                return (rowIndexes, columnIndexes, values);
            }
            #endregion

            //  ██╗ 
            //  ╚██╗ only back potential is simulated
            //  ██╔╝
            //  ╚═╝
            #region only back potential is simulated
            else
            {
                // structure of Jacobimatrix:
                //
                //                  /  ∂F_1b/∂Φ_1b   · · ·   ∂Igen_1b/∂Φ_1b  \
                //                 |                                          |
                //                 |        ·         ·             ·         |
                //           J  =  |        ·           ·           ·         |
                //                 |        ·             ·         ·         |
                //                 |                                          |
                //                  \  ∂Igen_Nb/∂Φ_1b   · · ·   ∂F_Nb/∂Φ_Nb  /

                // Get amount of non-zero elements (amount of points/rows plus sum of neighbors of all points)
                int amountOfNonZeroElements = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                    amountOfNonZeroElements += 1 + mesh.finiteElements[k].neighbors.Count;

                // Initialize arrays, which will be used for creating the sparse matrix
                int[] rowIndexes = new int[amountOfNonZeroElements];
                int[] columnIndexes = new int[amountOfNonZeroElements];
                double[] values = new double[amountOfNonZeroElements];

                // Fill arrays at non-zero positions
                int indexInSparseArrays = 0;
                for (int k = startPointIndexIncluding; k < endPointIndexExcluding; k++)
                {
                    // residual at back contact (if external cell contact: J_ij = J_ij * Rcontact + 1) ——————————————————————————————————————————————————

                    // diagonal element
                    //                        ___
                    //            ∂F_k        \    I_k,n       I_generated
                    //           ——————  =     ⟩   ———————  -  —————————————
                    //            ∂Φ_k        /     ∂Φ_k           ∂Φ_k
                    //                        ‾‾‾
                    //                      n ∈ N(k)

                    // ∂F_b/∂Φ_b
                    rowIndexes[indexInSparseArrays] = k;
                    columnIndexes[indexInSparseArrays] = k;
                    values[indexInSparseArrays++] = mesh.finiteElements[k].neighbors.Select(n =>
                        mesh.finiteElements[k].GetCurrentToNeighborsDerivativeBackSingle(mesh, n.index)).Sum()
                        - mesh.finiteElements[k].GetCurrentGenerated_DerivativePhiBack();
                    if (mesh.finiteElements[k].isExternalCellBackContact)
                        values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellBackContact + 1;

                    // offdiagonal element
                    //
                    //            ∂F_k       ∂I_k,n
                    //           ——————  =  ————————
                    //            ∂Φ_n        ∂Φ_n
                    //
                    for (int n = 0; n < mesh.finiteElements[k].neighbors.Count; n++)
                    {
                        // ∂F_b/∂Φ_b
                        rowIndexes[indexInSparseArrays] = k;
                        columnIndexes[indexInSparseArrays] = mesh.finiteElements[k].neighbors[n].index;
                        values[indexInSparseArrays++] = -mesh.finiteElements[k].GetCurrentToNeighborsDerivativeBackSingle(mesh, mesh.finiteElements[k].neighbors[n].index);
                        if (mesh.finiteElements[k].isExternalCellBackContact)
                            values[indexInSparseArrays - 1] = values[indexInSparseArrays - 1] * mesh.finiteElements[k].contactResistanceExternalCellBackContact;
                    }
                }

                return (rowIndexes, columnIndexes, values);
            }
            #endregion
        }

        // Topological Optimization █████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// returns the vector with all front phi values
        /// </summary>
        public Vector<double> PhiFrontVector()
        {
            return Vector.Create(mesh.finiteElements.Select(p => p.Value.phiFront).ToArray());
        }
        /// <summary>
        /// returns the vector with all generated current values
        /// </summary>
        public Vector<double> IgenVector()
        {
            var Igen = Vector.Create<double>(mesh.nextAvailableFiniteElementIndex);
            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (mesh.finiteElements[i].isExternalCellFrontContact)
                    Igen[i] = 1;
                else
                    Igen[i] = mesh.finiteElements[i].GetCurrentGenerated();
            }

            return Igen;
        }
        /// <summary>
        /// returns the stiffness matrix G with all conductivities
        /// </summary>
        public SparseMatrix<double> ConductivityStiffnessMatrix()
        {
            int amountOfNonZeroElements = 0;
            for (int k = 0; k < mesh.nextAvailableFiniteElementIndex; k++)
            {
                if (mesh.finiteElements[k].isExternalCellFrontContact)
                    amountOfNonZeroElements += 1;
                else
                    amountOfNonZeroElements += 1 + mesh.finiteElements[k].neighbors.Count;
            }


            int[] rowIndexes = new int[amountOfNonZeroElements];
            int[] columnIndexes = new int[amountOfNonZeroElements];
            double[] values = new double[amountOfNonZeroElements];

            int indexInSparseArrays = 0;
            for (int k = 0; k < mesh.nextAvailableFiniteElementIndex; k++)
            {
                if (mesh.finiteElements[k].isExternalCellFrontContact)
                {
                    rowIndexes[indexInSparseArrays] = k;
                    columnIndexes[indexInSparseArrays] = k;
                    values[indexInSparseArrays++] = -1;
                }
                else
                {
                    double sigmaSum = 0;
                    for (int n_i = 0; n_i < mesh.finiteElements[k].neighbors.Count; n_i++)
                    {
                        int n = mesh.finiteElements[k].neighbors[n_i].index; // global index of neighbor
                        double matrixEntry = -1 / mesh.finiteElements[k].resistorsFront[n_i];
                        rowIndexes[indexInSparseArrays] = k;
                        columnIndexes[indexInSparseArrays] = n;
                        values[indexInSparseArrays++] = matrixEntry;
                        sigmaSum += matrixEntry;
                    }
                    rowIndexes[indexInSparseArrays] = k;
                    columnIndexes[indexInSparseArrays] = k;
                    values[indexInSparseArrays++] = -sigmaSum;
                }
            }
            SparseMatrix<double> G = Matrix.CreateSparse(mesh.nextAvailableFiniteElementIndex, mesh.nextAvailableFiniteElementIndex, rowIndexes, columnIndexes, values);
            return G;
        }
        /// <summary>
        /// returns the derivation of the stiffness matrix G with all conductivities
        /// </summary>
        public SparseMatrix<double> ConductivityStiffnessMatrix_DerivativeXe(SparseMatrix<double> ConductivityStiffnessMatrix, int indexOfXe, Vector<bool> TOdensityArray_isFixed)
        {
            List<int> rowIndexes = new List<int>();
            List<int> columnIndexes = new List<int>();
            List<double> values = new List<double>();

            double sigmaSum = 0;
            foreach (var n in mesh.finiteElements[indexOfXe].neighbors)
            {

                double delSigma_delXe = Math.Pow(ConductivityStiffnessMatrix[indexOfXe, n.index], 2) * mesh.finiteElements[indexOfXe].ResistorToEdge_DerivativeXe(n.index, TOdensityArray_isFixed);

                if (delSigma_delXe != 0)
                {
                    rowIndexes.Add(indexOfXe);
                    columnIndexes.Add(n.index);
                    values.Add(mesh.finiteElements[indexOfXe].isExternalCellFrontContact ? 0 : delSigma_delXe);

                    rowIndexes.Add(n.index);
                    columnIndexes.Add(indexOfXe);
                    values.Add(delSigma_delXe);

                    sigmaSum += delSigma_delXe;

                    rowIndexes.Add(n.index);
                    columnIndexes.Add(n.index);
                    values.Add(-delSigma_delXe);
                }
            }

            rowIndexes.Add(indexOfXe);
            columnIndexes.Add(indexOfXe);
            values.Add(mesh.finiteElements[indexOfXe].isExternalCellFrontContact ? 0 : -sigmaSum);

            SparseMatrix<double> delG_delXe = Matrix.CreateSparse(mesh.nextAvailableFiniteElementIndex, mesh.nextAvailableFiniteElementIndex, rowIndexes.ToArray(), columnIndexes.ToArray(), values.ToArray());
            return delG_delXe;
        }
        /// <summary>
        /// returns the derivation of the generated current Vector with respect to the given design parameter
        /// </summary>
        public Vector<double> GeneratedCurrent_DerivativeXe(int indexOfXe)
        {
            var delIgen_delXe = Vector.CreateSparse<double>(mesh.nextAvailableFiniteElementIndex);
            delIgen_delXe[indexOfXe] = mesh.finiteElements[indexOfXe].GetCurrentGenerated_DerivativeXe();
            return delIgen_delXe;
        }
        /// <summary>
        /// returns the derivation of the generated current Vector with respect to the design parameter vector
        /// </summary>
        public SparseMatrix<double> GeneratedCurrent_DerivativeX(Vector<bool> TOdensityArray_isFixed)
        {
            /*int[] rowsAndColumns = new int[mesh.nextAvailablePointIndex];
            double[] values = new double[mesh.nextAvailablePointIndex];

            for (int i = 0; i < mesh.nextAvailablePointIndex; i++)
            {
                if (!fixedGridDensities.Select(e => e.index).Any(e => e == i))
                {
                    rowsAndColumns[i] = i;
                    values[i] = mesh.points[i].GetCurrentGenerated_DerivativeXe();
                }
            }

            return Matrix.CreateSparse(mesh.nextAvailablePointIndex, mesh.nextAvailablePointIndex, rowsAndColumns, rowsAndColumns, values);*/

            List<int> rowsAndColumns = new List<int>();
            List<double> values = new List<double>();

            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (!TOdensityArray_isFixed[i])
                {
                    rowsAndColumns.Add(i);
                    values.Add(mesh.finiteElements[i].GetCurrentGenerated_DerivativeXe());
                }
            }

            var rowsAndColumns_array = rowsAndColumns.ToArray();
            return Matrix.CreateSparse(mesh.nextAvailableFiniteElementIndex, mesh.nextAvailableFiniteElementIndex, rowsAndColumns_array, rowsAndColumns_array, values.ToArray());
        }
        /// <summary>
        /// returns the derivation of the generated current Vector with respect to the design parameter vector
        /// </summary>
        public SparseMatrix<double> GeneratedCurrent_DerivativeX_Manipulated(Vector<bool> TOdensityArray_isFixed)
        {
            /*int[] rowsAndColumns = new int[mesh.nextAvailablePointIndex];
            double[] values = new double[mesh.nextAvailablePointIndex];

            for (int i = 0; i < mesh.nextAvailablePointIndex; i++)
            {
                if (!fixedGridDensities.Select(e => e.index).Any(e => e == i))
                {
                    rowsAndColumns[i] = i;
                    values[i] = mesh.points[i].isExternalCellFrontContact ? 0 : mesh.points[i].GetCurrentGenerated_DerivativeXe();
                }
            }

            return Matrix.CreateSparse(mesh.nextAvailablePointIndex, mesh.nextAvailablePointIndex, rowsAndColumns, rowsAndColumns, values);*/

            List<int> rowsAndColumns = new List<int>();
            List<double> values = new List<double>();

            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                if (!TOdensityArray_isFixed[i])
                {
                    rowsAndColumns.Add(i);
                    values.Add(mesh.finiteElements[i].isExternalCellFrontContact ? 0 : mesh.finiteElements[i].GetCurrentGenerated_DerivativeXe());
                }
            }

            var rowsAndColumns_array = rowsAndColumns.ToArray();
            return Matrix.CreateSparse(mesh.nextAvailableFiniteElementIndex, mesh.nextAvailableFiniteElementIndex, rowsAndColumns_array, rowsAndColumns_array, values.ToArray());
        }
        /// <summary>
        /// returns the derivation of the generated current Vector with respect to the voltage vector
        /// </summary>
        public SparseMatrix<double> GeneratedCurrent_DerivativePhiFront()
        {
            int[] rowsAndColumns = new int[mesh.nextAvailableFiniteElementIndex];
            double[] values = new double[mesh.nextAvailableFiniteElementIndex];

            for (int i = 0; i < mesh.nextAvailableFiniteElementIndex; i++)
            {
                rowsAndColumns[i] = i;
                values[i] = mesh.finiteElements[i].GetCurrentGenerated_DerivativePhiFront();
            }

            return Matrix.CreateSparse(mesh.nextAvailableFiniteElementIndex, mesh.nextAvailableFiniteElementIndex, rowsAndColumns, rowsAndColumns, values);
        }
    }
}