using BasicLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Geometry
{
    public class GeometryFileData
    {
        /// <summary>
        /// length scale of the geometry
        /// </summary>
        public string unit { get; private set; }

        /// <summary>
        /// materials with all properties in a double list
        /// </summary>
        public Dictionary<int, List<double>> materials { get; private set; }

        public Dictionary<int, List<double>> incoherentLayersBefore { get; private set; }
        public Dictionary<int, List<double>> coherentLayersBefore { get; private set; }
        public Dictionary<int, List<double>> coherentLayersBehind { get; private set; }
        public Dictionary<int, List<double>> incoherentLayersBehind { get; private set; }

        public Dictionary<int, List<double>> gradingInfo { get; set; }

        public int materialBefore { get; private set; }
        public (int ID, double roughnessOnTop) materialBehind { get; private set; }

        /// <summary>
        /// optical layerstacks in a list
        /// </summary>
        public Dictionary<int, (double roughnessFrontGrid, double roughnessFrontContact, double roughnessAbsorber, double roughnessBackContact, double roughnessBackGrid,
            int IDbefore, (int ID, double roughness) behind, List<(int ID, double thickness)> incoherent,
            List<(int ID, double thickness, double roughness)> aboveFrontGrid, List<(int ID, double thickness, double roughness)> aboveAbsorber,
            List<(int ID, double thickness, double roughness)> belowAbsorber, List<(int ID, double thickness, double roughness)> belowBackGrid)> opticalModels
        { get; private set; }


        /// <summary>
        /// additional optical layers in a semiconductor stack in a list
        /// </summary>
         /*
        public (int materialBeforeID, int materialBehindID, List<(int ID, double thickness, double roughness)> semiconductorStack,
            List<(int ID, double thickness)> incoherentBeforeStack, List<(int ID, double thickness, double roughness)> coherentBeforeStack,
            List<(int ID, double thickness, double roughness)> coherentBehindStack, List<(int ID, double thickness)> incoherentBehindStack) opticalmodelSemiconductor
        { get; private set; }
        */

        /// <summary>
        /// boundary conditions
        /// </summary>
        public (List<(int index, int selector, List<double> conditions)> boundaryPoints,
            List<(int index, int selector, List<double> conditions)> boundarySegments,
            List<(int index, int selector, List<double> conditions)> boundaryAreas) boundaryConditions;

        /// <summary>
        /// determines the dimension of this geometry file
        /// </summary>
        public int dimension { get; set; }

        /// <summary>
        /// constructor of parent class
        /// </summary>
        /// <param name="geometryLines">lines with all properties</param>
        public GeometryFileData(string[] geometryLines)
        {
            // initialize
            materials = new Dictionary<int, List<double>>();
            incoherentLayersBefore = new Dictionary<int, List<double>>();
            coherentLayersBefore = new Dictionary<int, List<double>>();
            coherentLayersBehind = new Dictionary<int, List<double>>();
            incoherentLayersBehind = new Dictionary<int, List<double>>();
            boundaryConditions.boundaryPoints = new List<(int index, int selector, List<double> conditions)>();
            boundaryConditions.boundarySegments = new List<(int index, int selector, List<double> conditions)>();
            boundaryConditions.boundaryAreas = new List<(int index, int selector, List<double> conditions)>();

            // Read unit
            int lengthScaleLine = InputOutput.GetLineOfStringInArray(geometryLines, "length_scale:") + 1;
            unit = (geometryLines[lengthScaleLine].Trim());

            // Read materials before and after layerstack
            int materialBeforeLine = InputOutput.GetLineOfStringInArray(geometryLines, "material_before:") + 1;
            if (materialBeforeLine != 0)
            {
                materialBefore = InputOutput.ToIntWithArbitrarySeparator(geometryLines[materialBeforeLine]);
                int materialBehindLine = InputOutput.GetLineOfStringInArray(geometryLines, "material_behind:") + 1;
                materialBehind = (InputOutput.ToIntWithArbitrarySeparator(geometryLines[materialBehindLine].Split('\t')[0]), InputOutput.ToDoubleWithArbitrarySeparator(geometryLines[materialBehindLine].Split('\t')[1]));
            }

            //Gradings
            int gradingLines = InputOutput.GetLineOfStringInArray(geometryLines, "gradings:") + 1;
            if (gradingLines != 0)
            {
                gradingInfo = InputOutput.GetDoubleTupleOfStringArray(geometryLines, gradingLines);

            }

            // Read materials
            int materialsLine = InputOutput.GetLineOfStringInArray(geometryLines, "materials:") + 1;
            materials = InputOutput.GetDoubleTupleOfStringArray(geometryLines, materialsLine);


            //Read additonal optical semiconductor layers
            int additionalIncoherentBeforeStart = InputOutput.GetLineOfStringInArray(geometryLines, "incoh_before:") + 1;
            int additionalCoherentBeforeStart = InputOutput.GetLineOfStringInArray(geometryLines, "coherent_before:") + 1;
            int additionalCoherentBehindStart = InputOutput.GetLineOfStringInArray(geometryLines, "coherent_behind:") + 1;
            int additionalIncoherentBehindStart = InputOutput.GetLineOfStringInArray(geometryLines, "incoh_behind:") + 1;

            if (additionalIncoherentBeforeStart != 0)
            {
                incoherentLayersBefore = InputOutput.GetDoubleTupleOfStringArray(geometryLines, additionalIncoherentBeforeStart);
            }
            if (additionalCoherentBeforeStart != 0)
            {
                coherentLayersBefore = InputOutput.GetDoubleTupleOfStringArray(geometryLines, additionalCoherentBeforeStart);
            }
            if (additionalCoherentBehindStart != 0)
            {
                coherentLayersBehind = InputOutput.GetDoubleTupleOfStringArray(geometryLines, additionalCoherentBehindStart);
            }
            if (additionalIncoherentBehindStart != 0)
            {
                incoherentLayersBehind = InputOutput.GetDoubleTupleOfStringArray(geometryLines, additionalIncoherentBehindStart);
            }


            // Read optical models
            int opticsLineStart = InputOutput.GetLineOfStringInArray(geometryLines, "opticalModels:") + 1;
            if (opticsLineStart != 0)
            {
                int opticsLineStop = opticsLineStart;
                for (int i = opticsLineStart; i < geometryLines.Length; i++)
                    if (!geometryLines[i].Contains("\t"))
                    {
                        opticsLineStop = i - 1;
                        break;
                    }
                var opticalModelsInString = InputOutput.ReadLinesToString2DJaggedArray(geometryLines.Skip(opticsLineStart).Take(opticsLineStop - opticsLineStart + 1).ToArray());

                opticalModels = new Dictionary<int, (double roughnessFrontGrid, double roughnessFrontContact, double roughnessAbsorber, double roughnessBackContact, double roughnessBackGrid,
                    int IDbefore, (int ID, double roughness) behind, List<(int ID, double thickness)> incoherent,
                    List<(int ID, double thickness, double roughness)> aboveFrontGrid, List<(int ID, double thickness, double roughness)> aboveAbsorber,
                    List<(int ID, double thickness, double roughness)> belowAbsorber, List<(int ID, double thickness, double roughness)> belowBackGrid)>();
                foreach (var model in opticalModelsInString)
                {
                    double roughnessFrontGrid = double.Parse(model[1]);
                    double roughnessFrontContact = double.Parse(model[2]);
                    double roughnessAbsorber = double.Parse(model[3]);
                    double roughnessBackContact = double.Parse(model[4]);
                    double roughnessBackGrid = double.Parse(model[5]);

                    int IDbefore = int.Parse(model[6]);

                    (int ID, double roughness) behind = (int.Parse(model[7]), double.Parse(model[8]));

                    int currentIndex = 10;
                    List<(int ID, double thickness)> incoherent = new List<(int ID, double thickness)>();
                    for (; currentIndex < model.Length && !model[currentIndex].Equals("aboveFrontGrid"); currentIndex += 2)
                        incoherent.Add((int.Parse(model[currentIndex]), double.Parse(model[currentIndex + 1])));

                    currentIndex++;
                    List<(int ID, double thickness, double roughness)> aboveFrontGrid = new List<(int ID, double thickness, double roughness)>();
                    for (; currentIndex < model.Length && !model[currentIndex].Equals("aboveAbsorber"); currentIndex += 3)
                        aboveFrontGrid.Add((int.Parse(model[currentIndex]), double.Parse(model[currentIndex + 1]), double.Parse(model[currentIndex + 2])));

                    currentIndex++;
                    List<(int ID, double thickness, double roughness)> aboveAbsorber = new List<(int ID, double thickness, double roughness)>();
                    for (; currentIndex < model.Length && !model[currentIndex].Equals("belowAbsorber"); currentIndex += 3)
                        aboveAbsorber.Add((int.Parse(model[currentIndex]), double.Parse(model[currentIndex + 1]), double.Parse(model[currentIndex + 2])));

                    currentIndex++;
                    List<(int ID, double thickness, double roughness)> belowAbsorber = new List<(int ID, double thickness, double roughness)>();
                    for (; currentIndex < model.Length && !model[currentIndex].Equals("belowBackGrid"); currentIndex += 3)
                        belowAbsorber.Add((int.Parse(model[currentIndex]), double.Parse(model[currentIndex + 1]), double.Parse(model[currentIndex + 2])));

                    currentIndex++;
                    List<(int ID, double thickness, double roughness)> belowBackGrid = new List<(int ID, double thickness, double roughness)>();
                    for (; currentIndex < model.Length; currentIndex += 3)
                        belowBackGrid.Add((int.Parse(model[currentIndex]), double.Parse(model[currentIndex + 1]), double.Parse(model[currentIndex + 2])));

                    opticalModels.Add(int.Parse(model[0]), (roughnessFrontGrid, roughnessFrontContact, roughnessAbsorber, roughnessBackContact, roughnessBackGrid, IDbefore, behind, incoherent, aboveFrontGrid, aboveAbsorber, belowAbsorber, belowBackGrid));
                }
            }

            // Read boundary conditions
            int lineIndexOfStringAppearance = InputOutput.GetLineOfStringInArray(geometryLines, "boundary_conditions:");
            for (int i = lineIndexOfStringAppearance + 1; true; i++)
            {
                try
                {
                    // assign zu points, segments or areas
                    string assignString = geometryLines[i].Split('\t')[0];
                    string indexString = Regex.Replace(assignString, "[^0-9.]", "");

                    // get selector
                    int selector = InputOutput.ToIntWithArbitrarySeparator(geometryLines[i].Split('\t')[1].Trim());

                    // get boundaryCondition
                    List<double> boundaryConditionList;
                    if (geometryLines[i].Split('\t').Length > 2)
                        boundaryConditionList = geometryLines[i].Split('\t').Skip(2).Select(s => InputOutput.ToDoubleWithArbitrarySeparator(s)).ToList();
                    else
                        boundaryConditionList = new List<double> { 0, 0, 0, 0 };

                    // set
                    if (assignString.Contains("point"))
                        boundaryConditions.boundaryPoints.Add((Convert.ToInt32(indexString.Trim()), selector, boundaryConditionList));
                    else if (assignString.Contains("segment"))
                        boundaryConditions.boundarySegments.Add((Convert.ToInt32(indexString.Trim()), selector, boundaryConditionList));
                    else if (assignString.Contains("area"))
                        boundaryConditions.boundaryAreas.Add((Convert.ToInt32(indexString.Trim()), selector, boundaryConditionList));
                }
                catch
                {
                    break;
                }
            }
        }

        /// <summary>
        /// Returns regions and additional points
        /// </summary>
        /// <param name="unitAsSeperateOuput">give back the unit as a string</param>
        public virtual (List<R> regions, List<ContourJunction> additionalContourJunctions, string unit) GetRegionsAndAdditionalPoints<R>(bool unitAsSeperateOuput = false)
            where R : Region, new()
        {
            Console.WriteLine("Override the function \"" + System.Reflection.MethodBase.GetCurrentMethod().Name
                + "()\" in the class \"" + this.GetType().Name + "\"!");
            return (null, null, null);
        }

        /// <summary>
        /// prints length scale in console
        /// </summary>
        public void PrintTopPart()
        {
            Console.WriteLine("length_scale:");
            Console.WriteLine(unit);
            Console.WriteLine();
        }
        /// <summary>
        /// prints materials and boundary conditions in console
        /// </summary>
        public void PrintBottomPart()
        {
            Console.WriteLine("materials:");
            foreach (var material in materials)
            {
                Console.Write(material.Key);
                foreach (var entry in material.Value)
                    Console.Write("\t" + entry);
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("boundary conditions:");
            foreach (var point in boundaryConditions.boundaryPoints)
                Console.WriteLine("point " + point.index + "\t" + point.selector + "\t" + point.conditions);
            foreach (var segment in boundaryConditions.boundarySegments)
                Console.WriteLine("segment " + segment.index + "\t" + segment.selector + "\t" + segment.conditions);
            foreach (var area in boundaryConditions.boundaryAreas)
                Console.WriteLine("area " + area.index + "\t" + area.selector + "\t" + area.conditions);
        }
    }

    public class GeometryFileData1D : GeometryFileData
    {
        public Dictionary<int, double> points { get; private set; }
        public Dictionary<int, (int pointIndex1, int pointIndex2)> segments { get; private set; }

        public GeometryFileData1D(string[] geometryLines) : base(geometryLines)
        {
            dimension = 1;

            // initialize
            points = new Dictionary<int, double>();
            segments = new Dictionary<int, (int pointIndex1, int pointIndex2)>();

            // points
            int pointsLine = InputOutput.GetLineOfStringInArray(geometryLines, "points:") + 1;
            foreach (var point in InputOutput.GetDoubleTupleOfStringArray(geometryLines, pointsLine))
                points.Add(point.Key, point.Value[0]);

            // segments
            int segmentsLine = InputOutput.GetLineOfStringInArray(geometryLines, "segments:") + 1;
            foreach (var segment in InputOutput.GetIntTupleOfStringArray(geometryLines, segmentsLine))
                segments.Add(segment.Key, (segment.Value[0], segment.Value[1]));
        }

        /// <summary>
        /// Returns regions and additional points
        /// </summary>
        /// <param name="unitAsSeperateOuput">give back the unit as a string</param>
        public override (List<R> regions, List<ContourJunction> additionalContourJunctions, string unit) GetRegionsAndAdditionalPoints<R>(bool unitAsSeperateOuput = false)
        {
            //  ██╗ 
            //  ╚██╗ points
            //  ██╔╝
            //  ╚═╝
            Dictionary<int, ContourJunction> contourJunctions = new Dictionary<int, ContourJunction>();
            foreach (var point in points)
            {
                if (unitAsSeperateOuput)
                    contourJunctions.Add(point.Key, new ContourJunction(point.Key, new Position(point.Value, 0)));
                else
                    contourJunctions.Add(point.Key, new ContourJunction(point.Key, new Position(point.Value * InputOutput.GetUnitFromString(unit), 0)));
            }

            //  ██╗ 
            //  ╚██╗ segments
            //  ██╔╝
            //  ╚═╝
            Dictionary<int, ContourSegment> contourSegments = new Dictionary<int, ContourSegment>();
            foreach (var segment in segments)
            {
                contourSegments.Add(segment.Key, new ContourSegment(segment.Key, contourJunctions[segment.Value.pointIndex1], contourJunctions[segment.Value.pointIndex2]));
                contourJunctions[segment.Value.pointIndex1].neighbors.Add(contourJunctions[segment.Value.pointIndex2]);
                contourJunctions[segment.Value.pointIndex1].adjacentContourSegments.Add(contourSegments.Last().Value);
                contourJunctions[segment.Value.pointIndex2].neighbors.Add(contourJunctions[segment.Value.pointIndex1]);
                contourJunctions[segment.Value.pointIndex2].adjacentContourSegments.Add(contourSegments.Last().Value);
            }

            //  ██╗ 
            //  ╚██╗ regions
            //  ██╔╝
            //  ╚═╝
            List<R> regions = new List<R>();
            foreach (var segment in contourSegments.Values)
            {
                R region = new R();
                if (materials[segment.index].Count == 1 && materials[segment.index][0] == 0)
                {
                    // simulation hole
                    region.Initialize(segment.index, new List<ContourSegment> { segment }, pointType.hole);
                    regions.Add(region);
                }
                else
                {
                    if (materials[segment.index].Count == 3)
                    {
                        // semiconductor point
                        region.Initialize(segment.index, new List<ContourSegment> { segment }, pointType.semiconductor);

                        //Set prefernces array

                        if (gradingInfo.ContainsKey(segment.index))
                        {
                            region.SetProperties(new double[] {
                            (int)Math.Round(materials[segment.index][0]), materials[segment.index][1], materials[segment.index][2],
                            gradingInfo[segment.index][0], gradingInfo[segment.index][1], gradingInfo[segment.index][2], gradingInfo[segment.index][3], gradingInfo[segment.index][4],
                            gradingInfo[segment.index][5], gradingInfo[segment.index][6], gradingInfo[segment.index][7], gradingInfo[segment.index][8], gradingInfo[segment.index][9]
                        }, pointType.semiconductor);
                        }
                        else
                        {
                            region.SetProperties(new double[] {
                            (int)Math.Round(materials[segment.index][0]), materials[segment.index][1], materials[segment.index][2],
                            double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN
                        }, pointType.semiconductor);
                        }


                        regions.Add(region);
                    }
                    else
                    {
                        // cell point
                        pointType regionType = pointType.cell;
                        if (materials[segment.index].Count > 11)
                        {
                            switch (materials[segment.index][11])
                            {
                                case 8001:
                                    regionType = pointType.P1;
                                    break;
                                case 8012:
                                    regionType = pointType.gap12;
                                    break;
                                case 8002:
                                    regionType = pointType.P2;
                                    break;
                                case 8023:
                                    regionType = pointType.gap23;
                                    break;
                                case 8003:
                                    regionType = pointType.P3;
                                    break;
                            }
                        }
                        region.Initialize(segment.index, new List<ContourSegment> { segment }, regionType);
                        region.SetProperties(materials[segment.index].ToArray(), regionType);
                        region.SetOpticalModel(opticalModels[segment.index]);
                        regions.Add(region);
                    }
                }
            }

            return (regions, contourJunctions.Values.Where(p => p.neighbors.Count == 0).ToList(), unit);
        }

        /// <summary>
        /// prints to console
        /// </summary>
        public void Print()
        {
            PrintTopPart();

            Console.WriteLine("points:");
            foreach (var point in points)
                Console.WriteLine(point.Key + "\t" + point.Value);
            Console.WriteLine();

            Console.WriteLine("segments:");
            foreach (var segment in segments)
                Console.WriteLine(segment.Key + "\t" + segment.Value.pointIndex1 + "\t" + segment.Value.pointIndex2);
            Console.WriteLine();

            PrintBottomPart();
        }
    }

    public class GeometryFileData2D : GeometryFileData
    {
        public Dictionary<int, (double x, double y)> points { get; private set; }
        public Dictionary<int, (int pointIndex1, int pointIndex2)> segments { get; private set; }
        public Dictionary<int, List<int>> areas { get; private set; }

        public GeometryFileData2D(string[] geometryLines) : base(geometryLines)
        {
            dimension = 2;

            // initialize
            points = new Dictionary<int, (double, double)>();
            segments = new Dictionary<int, (int pointIndex1, int pointIndex2)>();
            areas = new Dictionary<int, List<int>>();

            // points
            int pointsLine = InputOutput.GetLineOfStringInArray(geometryLines, "points:") + 1;
            foreach (var point in InputOutput.GetDoubleTupleOfStringArray(geometryLines, pointsLine))
                points.Add(point.Key, (point.Value[0], point.Value[1]));

            // segments
            int segmentsLine = InputOutput.GetLineOfStringInArray(geometryLines, "segments:") + 1;
            foreach (var segment in InputOutput.GetIntTupleOfStringArray(geometryLines, segmentsLine))
                segments.Add(segment.Key, (segment.Value[0], segment.Value[1]));

            // areas
            int areasLine = InputOutput.GetLineOfStringInArray(geometryLines, "areas:") + 1;
            areas = InputOutput.GetIntTupleOfStringArray(geometryLines, areasLine);
        }

        /// <summary>
        /// Returns regions and additional points
        /// </summary>
        /// <param name="unitAsSeperateOuput">give back the unit as a string</param>
        public override (List<R> regions, List<ContourJunction> additionalContourJunctions, string unit) GetRegionsAndAdditionalPoints<R>(bool unitAsSeperateOuput = false)
        {
            //  ██╗ 
            //  ╚██╗ points
            //  ██╔╝
            //  ╚═╝
            Dictionary<int, ContourJunction> contourJunctions = new Dictionary<int, ContourJunction>();
            foreach (var point in points)
            {
                if (unitAsSeperateOuput)
                    contourJunctions.Add(point.Key, new ContourJunction(point.Key, new Position(point.Value.x, point.Value.y)));
                else
                    contourJunctions.Add(point.Key, new ContourJunction(point.Key, new Position(point.Value.x * InputOutput.GetUnitFromString(unit), point.Value.y * InputOutput.GetUnitFromString(unit))));
            }

            //  ██╗ 
            //  ╚██╗ segments
            //  ██╔╝
            //  ╚═╝
            Dictionary<int, ContourSegment> contourSegments = new Dictionary<int, ContourSegment>();
            foreach (var segment in segments)
            {
                contourSegments.Add(segment.Key, new ContourSegment(segment.Key, contourJunctions[segment.Value.pointIndex1], contourJunctions[segment.Value.pointIndex2]));
                contourJunctions[segment.Value.pointIndex1].neighbors.Add(contourJunctions[segment.Value.pointIndex2]);
                contourJunctions[segment.Value.pointIndex1].adjacentContourSegments.Add(contourSegments.Last().Value);
                contourJunctions[segment.Value.pointIndex2].neighbors.Add(contourJunctions[segment.Value.pointIndex1]);
                contourJunctions[segment.Value.pointIndex2].adjacentContourSegments.Add(contourSegments.Last().Value);
            }

            //  ██╗ 
            //  ╚██╗ regions
            //  ██╔╝
            //  ╚═╝
            List<R> regions = new List<R>();

            foreach (var area in areas)
            {
                List<ContourSegment> orderedSegments = new List<ContourSegment>();
                foreach (var s in area.Value)
                    orderedSegments.Add(contourSegments[s]);

                R region = new R();
                if (materials[area.Key].Count == 1 && materials[area.Key][0] == 0)
                {
                    // simulation hole
                    region.Initialize(area.Key, orderedSegments, pointType.hole);
                    regions.Add(region);
                }
                else
                {
                    if (materials[area.Key].Count == 3)
                    {
                        // semiconductor point
                        region.Initialize(area.Key, area.Value.Select(i => contourSegments[i]).ToList(), pointType.semiconductor);
                        //region.SetProperties(new double[] { (int)Math.Round(materials[area.Key][0]), materials[area.Key][1], materials[area.Key][2] }, pointType.semiconductor);
                        if (gradingInfo.ContainsKey(area.Key))
                        {
                            region.SetProperties(new double[] {
                            (int)Math.Round(materials[area.Key][0]), materials[area.Key][1], materials[area.Key][2],
                            gradingInfo[area.Key][0], gradingInfo[area.Key][1], gradingInfo[area.Key][2], gradingInfo[area.Key][3], gradingInfo[area.Key][4],
                            gradingInfo[area.Key][5], gradingInfo[area.Key][6], gradingInfo[area.Key][7], gradingInfo[area.Key][8], gradingInfo[area.Key][9]
                        }, pointType.semiconductor);
                        }
                        else
                        {
                            region.SetProperties(new double[] {
                            (int)Math.Round(materials[area.Key][0]), materials[area.Key][1], materials[area.Key][2],
                            double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN
                        }, pointType.semiconductor);
                        }
                        regions.Add(region);
                    }
                    else
                    {
                        // cell point
                        pointType regionType = pointType.cell;
                        if (materials[area.Key].Count > 11)
                        {
                            switch (materials[area.Key][11])
                            {
                                case 8001:
                                    regionType = pointType.P1;
                                    break;
                                case 8012:
                                    regionType = pointType.gap12;
                                    break;
                                case 8002:
                                    regionType = pointType.P2;
                                    break;
                                case 8023:
                                    regionType = pointType.gap23;
                                    break;
                                case 8003:
                                    regionType = pointType.P3;
                                    break;
                            }
                        }
                        region.Initialize(area.Key, area.Value.Select(i => contourSegments[i]).ToList(), regionType);
                        region.SetProperties(materials[area.Key].ToArray(), regionType);
                        region.SetOpticalModel(opticalModels[area.Key]);
                        regions.Add(region);
                    }
                }

                // write region to adjacentRegions of points and segments
                foreach (var segment in orderedSegments)
                {
                    if (!segment.adjacentRegions.Any(a => a.index == regions.Last().index))
                        segment.adjacentRegions.Add(regions.Last());
                    if (!segment.firstAdjacentContourJunction.adjacentRegions.Any(a => a.index == regions.Last().index))
                        segment.firstAdjacentContourJunction.adjacentRegions.Add(regions.Last());
                    if (!segment.secondAdjacentContourJunction.adjacentRegions.Any(a => a.index == regions.Last().index))
                        segment.secondAdjacentContourJunction.adjacentRegions.Add(regions.Last());
                }
            }

            return (regions, contourJunctions.Values.Where(p => p.neighbors.Count == 0).ToList(), unit);
        }

        /// <summary>
        /// prints to console
        /// </summary>
        public void Print()
        {
            PrintTopPart();

            Console.WriteLine("points:");
            foreach (var point in points)
                Console.WriteLine(point.Key + "\t" + point.Value);
            Console.WriteLine();

            Console.WriteLine("segments:");
            foreach (var segment in segments)
                Console.WriteLine(segment.Key + "\t" + segment.Value.pointIndex1 + "\t" + segment.Value.pointIndex2);
            Console.WriteLine();

            Console.WriteLine("areas:");
            foreach (var area in areas)
            {
                Console.Write(area.Key + " >>");
                foreach (var i in area.Value)
                    Console.Write("\t" + i);
                Console.WriteLine();
            }
            Console.WriteLine();

            PrintBottomPart();
        }
    }
}
