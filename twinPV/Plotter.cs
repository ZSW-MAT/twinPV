using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.ComponentModel;
using AtomicusChart.Interface;
using AtomicusChart.Interface.Data;
using AtomicusChart.Interface.DataReaders;
using AtomicusChart.Interface.GeometryFactory;
using AtomicusChart.Interface.PresentationData;
using AtomicusChart.Interface.PresentationData.BaseTypes;
using AtomicusChart.Interface.PresentationData.Collections;
using AtomicusChart.Interface.PresentationData.Primitives;
using AtomicusChart.Interface.UtilityTypes;
using AtomicusChart.ValueData.DataReaders;
using AtomicusChart.ValueData.PresentationData;
using Geometry;
using AtomicusChart.Interface.Processing;
using Cell;
using Semiconductor;
using BasicLib;
using MoreLinq;

namespace twinPV
{
    public static class Plotter
    {
        /// <summary>
        /// colormap "Rainbow" from Mathematica
        /// </summary>
        public static ColorMap colormap_Rainbow_Mathematica { get; private set; } = new ColorMap();
        /// <summary>
        /// colormap "Sunsetcolors"
        /// </summary>
        public static ColorMap colormap_Sunsetcolors { get; private set; } = new ColorMap();
        /// <summary>
        /// colormap cut from the map "Sunsetcolors"
        /// </summary>
        public static ColorMap colormap_Sunsetcolors_cut { get; private set; } = new ColorMap();

        /// <summary>
        /// colormap "Maple" from Origin Lab
        /// </summary>
        public static ColorMap colormap_Origin_Maple { get; private set; } = new ColorMap();

        /// <summary>
        /// colormap fitting to GUI accent colors
        /// </summary>
        public static ColorMap colormap_GUI { get; private set; } = new ColorMap();

        // Constructor ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        static Plotter()
        {

            // Rainbow Mathematica
            colormap_Rainbow_Mathematica.OutOfRangeTop = new Color4(218, 33, 34, 150);
            colormap_Rainbow_Mathematica.Top = new Color4(218, 33, 34, 150);
            ColorStop c9 = new ColorStop(new Color4(228, 101, 46, 150), 0.9f);
            ColorStop c8 = new ColorStop(new Color4(228, 153, 56, 150), 0.8f);
            ColorStop c7 = new ColorStop(new Color4(206, 182, 65, 150), 0.7f);
            ColorStop c6 = new ColorStop(new Color4(170, 189, 82, 150), 0.6f);
            ColorStop c5 = new ColorStop(new Color4(130, 186, 113, 150), 0.5f);
            ColorStop c4 = new ColorStop(new Color4(99, 172, 154, 150), 0.4f);
            ColorStop c3 = new ColorStop(new Color4(76, 143, 192, 150), 0.3f);
            ColorStop c2 = new ColorStop(new Color4(63, 99, 207, 150), 0.2f);
            ColorStop c1 = new ColorStop(new Color4(71, 46, 185, 150), 0.1f);
            colormap_Rainbow_Mathematica.SetColorStops(new List<ColorStop>() { c1, c2, c3, c4, c5, c6, c7, c8, c9 });
            colormap_Rainbow_Mathematica.Bottom = new Color4(120, 28, 134, 150);
            colormap_Rainbow_Mathematica.OutOfRangeBottom = new Color4(120, 28, 134, 150);
            colormap_Rainbow_Mathematica.NanColor = new Color4(0, 0, 0, 0);

            // Sunsetcolors Mathematica
            colormap_Sunsetcolors.OutOfRangeTop = new Color4(255, 255, 255, 255);
            colormap_Sunsetcolors.Top = new Color4(255, 255, 255, 255);
            ColorStop sunc9 = new ColorStop(new Color4(255, 237, 179), 0.9f);
            ColorStop sunc8 = new ColorStop(new Color4(255, 215, 107), 0.8f);
            ColorStop sunc7 = new ColorStop(new Color4(255, 185, 53), 0.7f);
            ColorStop sunc6 = new ColorStop(new Color4(253, 151, 25), 0.6f);
            ColorStop sunc5 = new ColorStop(new Color4(250, 116, 13), 0.5f);
            ColorStop sunc4 = new ColorStop(new Color4(221, 86, 47), 0.4f);
            ColorStop sunc3 = new ColorStop(new Color4(181, 60, 80), 0.3f);
            ColorStop sunc2 = new ColorStop(new Color4(116, 41, 117), 0.2f);
            ColorStop sunc1 = new ColorStop(new Color4(58, 21, 79), 0.1f);
            colormap_Sunsetcolors.SetColorStops(new List<ColorStop>() { sunc1, sunc2, sunc3, sunc4, sunc5, sunc6, sunc7, sunc8, sunc9 });
            colormap_Sunsetcolors.Bottom = new Color4(0, 0, 0, 255);
            colormap_Sunsetcolors.OutOfRangeBottom = new Color4(0, 0, 0, 255);
            colormap_Sunsetcolors.NanColor = new Color4(0, 0, 0, 0);

            // Sunsetcolors_cut
            colormap_Sunsetcolors_cut.OutOfRangeTop = new Color4(255, 188, 59, 255);
            colormap_Sunsetcolors_cut.Top = new Color4(255, 188, 59, 255);
            ColorStop suncutc6 = new ColorStop(new Color4(253, 151, 25), 6f / 7f);
            ColorStop suncutc5 = new ColorStop(new Color4(250, 116, 13), 5f / 7f);
            ColorStop suncutc4 = new ColorStop(new Color4(221, 86, 47), 4f / 7f);
            ColorStop suncutc3 = new ColorStop(new Color4(181, 60, 80), 3f / 7f);
            ColorStop suncutc2 = new ColorStop(new Color4(116, 41, 117), 2f / 7f);
            ColorStop suncutc1 = new ColorStop(new Color4(58, 21, 79), 1f / 7f);
            colormap_Sunsetcolors_cut.SetColorStops(new List<ColorStop>() { suncutc1, suncutc2, suncutc3, suncutc4, suncutc5, suncutc6 });
            colormap_Sunsetcolors_cut.Bottom = new Color4(0, 0, 0, 255);
            colormap_Sunsetcolors_cut.OutOfRangeBottom = new Color4(0, 0, 0, 255);
            colormap_Sunsetcolors_cut.NanColor = new Color4(0, 0, 0, 0);

            // Origin Maple
            colormap_Origin_Maple.OutOfRangeTop = new Color4(1, 70, 45);
            colormap_Origin_Maple.Top = new Color4(1, 70, 45);
            ColorStop maple0 = new ColorStop(new Color4(21, 87, 44), 8f / 9f);
            ColorStop maple1 = new ColorStop(new Color4(106, 137, 43), 7f / 9f);
            ColorStop maple2 = new ColorStop(new Color4(183, 184, 42), 6f / 9f);
            ColorStop maple3 = new ColorStop(new Color4(234, 246, 41), 5f / 9f);
            ColorStop maple4 = new ColorStop(new Color4(239, 174, 40), 4f / 9f);
            ColorStop maple5 = new ColorStop(new Color4(244, 120, 39), 3f / 9f);
            ColorStop maple6 = new ColorStop(new Color4(250, 86, 73), 2f / 9f);
            ColorStop maple7 = new ColorStop(new Color4(255, 0, 37), 1f / 9f);
            colormap_Origin_Maple.SetColorStops(new List<ColorStop>() { maple0, maple1, maple2, maple3, maple4, maple5, maple6, maple7 });
            colormap_Origin_Maple.Bottom = new Color4(255, 0, 30);
            colormap_Origin_Maple.OutOfRangeBottom = new Color4(255, 0, 30);
            colormap_Origin_Maple.NanColor = new Color4(0, 0, 0, 0);

            // New Color Map
            colormap_GUI.OutOfRangeTop = new Color4(57, 120, 106);
            colormap_GUI.Top = new Color4(57, 120, 106);
            ColorStop gui1 = new ColorStop(new Color4(80, 167, 148), 6f / 7f);
            ColorStop gui2 = new ColorStop(new Color4(80, 164, 167), 5f / 7f);
            ColorStop gui3 = new ColorStop(new Color4(80, 142, 167), 4f / 7f);
            ColorStop gui4 = new ColorStop(new Color4(80, 121, 167), 3f / 7f);
            ColorStop gui5 = new ColorStop(new Color4(80, 99, 167), 2f / 7f);
            ColorStop gui6 = new ColorStop(new Color4(83, 80, 167), 1f / 7f);
            colormap_GUI.SetColorStops(new List<ColorStop>() { gui1, gui2, gui3, gui4, gui5, gui6 });
            colormap_GUI.Bottom = new Color4(59, 57, 120);
            colormap_GUI.OutOfRangeBottom = new Color4(59, 57, 120);
            colormap_GUI.NanColor = new Color4(0, 0, 0, 0);
        }

        // 3D plots █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// create surface plot
        /// </summary>
        /// <param name="valueRange">enter custom value range for the colormap</param>
        public static ValueSurface PlotSurface(string name, bool visible, Vector3F[] dataSurface, ColorMap colormap,
            ValueSurfacePresentationType valueSurfacePresentationType, OneAxisBounds? valueRange = null, float height = (float)double.NaN)
        {
            if (dataSurface.Length <= 0)
                return new ValueSurface { };

            if (valueRange == null)
                valueRange = new OneAxisBounds(dataSurface.Min(d => d.Z), dataSurface.Max(d => d.Z), 0, 0.01f);

            IrregularValueSurfaceDataReader reader;
            if (double.IsNaN(height))
                reader = new IrregularValueSurfaceDataReader(dataSurface, dataSurface.Select(d => d.Z).ToArray(), 2, valueRange);
            else
                reader = new IrregularValueSurfaceDataReader(dataSurface.Select(p => new Vector3F(p.X, p.Y, height)).ToArray(), dataSurface.Select(d => d.Z).ToArray(), 2, valueRange);

            ValueSurface Surface = new ValueSurface
            {
                Reader = reader,
                PresentationType = valueSurfacePresentationType,
                Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                Name = name,
                ColorMapContainer = new ColorMapContainer { ColorMap = colormap.Clone(), ValueRange = valueRange },
                IsVisible = visible
            };

            return Surface;
        }
        /// <summary>
        /// create surface plot
        /// </summary>
        /// <param name="valueRange">enter custom value range for the colormap</param>
        public static ValueSurface PlotSurface(string name, bool visible, Vector3F[] dataSurface, Color4 color,
            ValueSurfacePresentationType valueSurfacePresentationType, OneAxisBounds? valueRange = null, float height = (float)double.NaN)
        {
            ColorMap colormap = new ColorMap();
            colormap.OutOfRangeTop = color;
            colormap.Top = color;
            ColorStop c1 = new ColorStop(color, 0.1f);
            colormap.SetColorStops(new List<ColorStop>() { c1 });
            colormap.Bottom = color;
            colormap.OutOfRangeBottom = color;
            colormap.NanColor = color;

            return PlotSurface(name, visible, dataSurface, colormap, valueSurfacePresentationType, valueRange, height);
        }
        /// <summary>
        /// create contourline-plot
        /// </summary>
        public static TriangleContoursRenderData PlotContours(string name, bool visible, Vector3F[] dataSurface, int amountOfContourLines, OneAxisBounds? valueRange = null)
        {
            if (dataSurface.Length <= 0)
                return new TriangleContoursRenderData { };

            if (valueRange == null)
                valueRange = new OneAxisBounds(dataSurface.Min(d => d.Z), dataSurface.Max(d => d.Z), 0, 0.01f);

            // edit contour lines
            ObservableCollection<Contour> contours = new ObservableCollection<Contour>(Enumerable.Range(0, amountOfContourLines).Select(i =>
               new Contour(i / (float)amountOfContourLines * (valueRange.Value.Max - valueRange.Value.Min) + valueRange.Value.Min, Colors.White, 1f)));
            IrregularValueSurfaceDataReader reader = new IrregularValueSurfaceDataReader(dataSurface, dataSurface.Select(d => d.Z).ToArray(), 2,
                valueRange);

            // Create the data that is responsible for contour visualization.
            TriangleContoursRenderData lines = new TriangleContoursRenderData
            {
                // Set contours computer source.
                DataSource = reader,
                // Set contours collection source.
                ContoursSource = new CustomContoursOwner(new ObservableCollection<Contour>(contours)),
                // Set name.
                Name = name,
                IsBoundsVisible = false,
                IsVisible = visible,
            };

            return lines;
        }
        /// <summary>
        /// create logarithmic current arrows
        /// </summary>
        public static SingleColorPrimitiveCollection PlotArrowsLog(string name, bool visible, Vector3F[] dataSurface, (float angle, float magnitude)[] currentdensities,
            float scalingFactor, Color4 color, List<(float angle, float magnitude)[]> currentdensitiesToNormalize = null)
        {
            if (dataSurface.Length != currentdensities.Length)
                throw new Exception("dataSurface and currentdensities must be of the same length");

            // return no arrows if all currents are zero
            if (currentdensities.Length != dataSurface.Length || currentdensities.Where(p => p.magnitude != 0).ToList().Count == 0)
            {
                return new SingleColorPrimitiveCollection(currentdensities.Take(1).Select(p => Matrix4F.Scaling(1e-10f, 1e-10f, 1e-10f)).ToArray())
                {
                    Name = name,
                    Mesh = ArrowMeshFactory.GenerateArrowX(20, 0.5f, 0.1f),
                    IsVisible = false,
                };
            }

            // calculate minimum value
            if (currentdensitiesToNormalize == null)
                currentdensitiesToNormalize = new List<(float angle, float magnitude)[]>() { currentdensities };
            double minLog = Math.Log10(currentdensitiesToNormalize.Where(d => d.Length > 0).Min(d => d.Select(p => p.magnitude).Where(p => p != 0).DefaultIfEmpty((float)double.MaxValue).Min()));

            // link data
            (Vector3F coordinates, (float angle, float magnitude) currentdensities)[] data = new (Vector3F dataSurface, (float angle, float magnitude) currentdensities)[dataSurface.Length];
            for (int i = 0; i < dataSurface.Length; i++)
                data[i] = (dataSurface[i], currentdensities[i]);

            // create, rotate and move arrows
            SingleColorPrimitiveCollection currentArrows = new SingleColorPrimitiveCollection(data.Where(p => p.currentdensities.magnitude != 0).Select(p =>
                    Matrix4F.Scaling(1f, 0.5f, 1e-9f) // Make arrow proportions more interesting
                * Matrix4F.Scaling((float)scalingFactor / (float)Math.Sqrt(data.Length)) // reduce arrow size to feet our
                * Matrix4F.Scaling((float)(Math.Log10(p.currentdensities.magnitude) - minLog)) // Set scale to magnitude parameter
                * Matrix4F.RotationAxis(Vector3F.UnitZ, p.currentdensities.angle) //Rotate arrow
                * Matrix4F.Translation(new Vector3F(p.coordinates.X, p.coordinates.Y,
                p.coordinates.Z + (data.Max(d => d.coordinates.Z) - data.Min(d => d.coordinates.Z)) * 1e-6f)) // Move arrow to its position
            ).ToArray())
            {
                Name = name,
                Color = color,
                Mesh = ArrowMeshFactory.GenerateArrowX(20, 0.5f, 0.1f),
                IsVisible = visible
            };

            return currentArrows;
        }

        public static CompositeRenderData PlotArrowsLogSymbols(Vector3F[] positions)
        {
            var arrows = new ObservableCollection<RenderData>();
            for (int i = 0; i < positions.Length; i++)
            {
                arrows.Add(
                    new Marker(new Arrow
                    {
                        Direction = positions[i],
                        Radius = 0.2f,
                    })
                    {
                        Transform = Matrix4F.Translation(positions[i]),
                        PixelSize = 100,
                    });
            }

            return new CompositeRenderData(arrows) { Name = "arrows", IsVisible = true };
        }

        // invisible boundaries █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// invisible boundaries
        /// </summary>
        public static Series PlotBoundaries(double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
        {
            var reader = new DefaultPositionMaskDataReader(new Vector3F[] { new Vector3F((float)minX, (float)minY, (float)minZ), new Vector3F((float)maxX, (float)maxY, (float)maxZ) });

            Series boundaries = new Series
            {
                Reader = reader,
                IsLegendVisible = false,
                Thickness = 0,
                PatternStyle = PatternStyle.Solid,
                MarkerSize = 0,
            };

            return boundaries;
        }

        // points and lines █████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// create point / line plot
        /// </summary>
        public static Series PlotPoints(string name, bool visible, Vector3F[] dataXYZ, int linethickness = 2, Color4? lineColor = null, MarkerStyle markerStyle = MarkerStyle.Circle, int markerSize = 5, Color4? markerColor = null, bool IsLegendVisible = true, PatternStyle patternStyle = PatternStyle.Solid)
        {
            Color4 lineColorUsed = lineColor ?? Colors.Black;
            Color4 markerColorUsed = markerColor ?? Colors.Black;

            var reader = new DefaultPositionMaskDataReader(dataXYZ.Where(d => !float.IsInfinity(d.X) && !float.IsInfinity(d.Y) && !float.IsInfinity(d.Z)).ToArray());

            Series Points = new Series
            {
                Name = name,
                Color = lineColorUsed,
                Thickness = linethickness,
                PatternStyle = patternStyle,
                Reader = reader,
                MarkerStyle = markerStyle,
                MarkerSize = markerSize,
                MarkerColor = markerColorUsed,
                IsVisible = visible,
                IsLegendVisible = IsLegendVisible,
            };

            return Points;
        }
        public static Prism PlotAreaUnderCurve(string name, bool visible, Vector2F[] dataXY, float thicknessOfLayer, float height, Color4 areaColor)
        {
            dataXY = dataXY.Where(d => !float.IsInfinity(d.X) && !float.IsInfinity(d.Y)).Concat(new[] { new Vector2F(dataXY.Last().X, 0), new Vector2F(dataXY.First().X, 0) }).Distinct().ToArray();
            /*
            int lastOfLeadingZeros = 0;
            for (int i = 0; i < dataXY.Length; i++)
            {
                if (dataXY[i].Y < 1e-5f)
                    lastOfLeadingZeros = i;
                else
                    break;
            }
            dataXY = dataXY.Skip(lastOfLeadingZeros + 1).ToArray();

            int firstOfEndingZeros = dataXY.Length;
            for (int i = dataXY.Length - 1; i >= 0; i--)
            {
                if (dataXY[i].Y < 1e-5f)
                    firstOfEndingZeros = i;
                else
                    break;
            }
            dataXY = dataXY.Take(lastOfLeadingZeros + 1 - 1).ToArray();

            Console.WriteLine();
            foreach (var d in dataXY)
                Console.WriteLine(d.X + "\t" + d.Y);*/

            Prism area = new Prism
            {
                Name = name,
                IsVisible = visible,
                Side = dataXY,
                BottomToTopVector = new Vector3F(0, 0, thicknessOfLayer),
                Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                Color = areaColor,
                Transform = Matrix4F.Translation(0, 0, height),
            };

            return area;
        }
        /// <summary>
        /// create point plot with colormaps
        /// </summary>
        public static CompositeRenderData PlotPointsWithColormap(string name, bool visible, Vector3F[] dataXYZ, MarkerStyle markerStyle, int markerSize, ColorMap markerColormap, OneAxisBounds? valueRange = null, float height = (float)double.NaN)
        {
            if (valueRange == null)
                valueRange = new OneAxisBounds(dataXYZ.Min(d => d.Z), dataXYZ.Max(d => d.Z), 0, 0.01f);

            var points = new ObservableCollection<RenderData>();
            for (var i = 0; i < dataXYZ.Length; i++)
            {
                DefaultPositionMaskDataReader reader;
                if (double.IsNaN(height))
                    reader = new DefaultPositionMaskDataReader(new Vector3F[] { dataXYZ[i] });
                else
                    reader = new DefaultPositionMaskDataReader(new Vector3F[] { new Vector3F(dataXYZ[i].X, dataXYZ[i].Y, height) });

                points.Add(new Series
                {
                    Name = name,
                    Thickness = 0,
                    PatternStyle = PatternStyle.Solid,
                    Reader = reader,
                    MarkerStyle = markerStyle,
                    MarkerSize = markerSize,
                    MarkerColor = markerColormap.GetColor((dataXYZ[i].Z - valueRange.Value.Min) / (valueRange.Value.Max - valueRange.Value.Min)),
                    IsVisible = visible,
                });
            }

            return new CompositeRenderData(points) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create point index plot
        /// </summary>
        public static CompositeRenderData PlotPointIndexes(string name, bool visible, Vector3F[] dataXYZ, int[] pointIndizes, Color4 textColor, Color4 backgroundColor, Color4 markerColor)
        {
            if (dataXYZ.Length != pointIndizes.Length)
                return new CompositeRenderData(new ObservableCollection<RenderData>()) { };

            var labels = new ObservableCollection<RenderData>();

            for (var i = 0; i < dataXYZ.Length; i++)
            {
                labels.Add(new Label
                {
                    Text = Convert.ToString(pointIndizes[i]),
                    FontFamily = "Arial",
                    FontSize = 14,
                    Transform = Matrix4F.Translation(dataXYZ[i].X, dataXYZ[i].Y, dataXYZ[i].Z),
                    FontColor = textColor,
                    Background = backgroundColor,
                    MarkerColor = markerColor,
                });
            }

            return new CompositeRenderData(labels) { Name = name, IsVisible = visible };
        }

        // plot 3D-model ████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// create 3D-model of the semiconductor
        /// </summary>
        public static CompositeRenderData Plot3Dmodel(string name, bool visible, ModelSemiconductor semiconductor, float height, float thicknessOfLayers, float multiplicatorXYaxis)
        {
            var objects = new ObservableCollection<RenderData>();

            foreach (var r in semiconductor.meshingAlgorithm.regions.Where(r => r.type != pointType.hole))
            {
                Color4 color = new Color4(100, 100, 100, 255);
                if (r.material.propertiesSemiconductor.NAminus > r.material.propertiesSemiconductor.NDplus) // if positive doped
                    color = new Color4(100, 0, 0, 255);
                if (r.material.propertiesSemiconductor.NAminus < r.material.propertiesSemiconductor.NDplus) // if negative doped
                    color = new Color4(0, 0, 100, 255);
                if (r.material.name == "CIGS_GrainBoundary")
                    color = new Color4(200, 200, 200, 255);
                if (r.material.name == "CdS")
                    color = new Color4(205, 173, 0, 255);
                if (r.material.name == "CdS")
                    color = new Color4(205, 173, 0, 255);
                if (r.material.name == "CdS")
                    color = new Color4(205, 173, 0, 255);
                if (r.material.name == "CdS")
                    color = new Color4(205, 173, 0, 255);


                objects.Add(new Prism
                {
                    // 2D-contours of the grid regions
                    Side = r.orderedPoints.Select(p => p.position).Select(p => new Vector2F((float)p.x * multiplicatorXYaxis, (float)p.y * multiplicatorXYaxis)).ToArray(),
                    // Vector that define translate between top and bottom sides
                    BottomToTopVector = new Vector3F(0, 0, thicknessOfLayers),
                    Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                    Color = color,
                    Transform = Matrix4F.Translation(0, 0, height),
                });
            }

            return new CompositeRenderData(objects) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create 3D-model fo the cell
        /// </summary>
        public static CompositeRenderData Plot3Dmodel(string name, bool visible, ModelCell cell, float height, float thicknessOfLayers, float multiplicatorXYaxis, double[] TOdensityArray = null)
        {
            var objects = new ObservableCollection<RenderData>();

            objects.Add(new Prism
            {
                Side = cell.meshingAlgorithm.outerContour.orderedPoints.Select(p => p.position).Select(p => new Vector2F((float)p.x * multiplicatorXYaxis, (float)p.y * multiplicatorXYaxis)).ToArray(),
                BottomToTopVector = new Vector3F(0, 0, thicknessOfLayers),
                Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                Color = new Color4(10, 10, 10, 255),
                Name = "active Area",
                Transform = Matrix4F.Translation(0, 0, height),
            });

            // plot FINITE ELEMENTS with grid
            /*foreach (var p in cell.mesh.finiteElements.Values.Where(r => r.frontGrid != null))
            {
                List<Position> corners = new List<Position>();
                foreach (var c in p.corners)
                    if (!corners.Any(g => g.SameAs(c.position, cell.meshingAlgorithm.defaultTolerance)))
                        if (!corners.Any(cor => (float)cor.x == (float)c.position.x && (float)cor.y == (float)c.position.y)) // check if also not already present in float precision
                            corners.Add(c.position);

                float TOthickness = thicknessOfLayers;
                if (!double.IsNaN(p.frontGridDensity))
                    TOthickness *= (float)p.frontGridDensity;
                if (TOthickness > 0)
                {
                    objects.Add(new Prism
                    {
                        Side = corners.Select(c => new Vector2F((float)(c.x * multiplicatorXYaxis), (float)(c.y * multiplicatorXYaxis))).ToArray(),
                        BottomToTopVector = new Vector3F(0, 0, TOthickness),
                        Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                        Color = double.IsNaN(p.frontGridDensity) ? new Color4(200, 200, 200, 255) : new Color4((byte)(200 * p.frontGridDensity), (byte)(200 * p.frontGridDensity), (byte)(200 * p.frontGridDensity), 255),
                        Name = "Grid",
                        Transform = Matrix4F.Translation(0, 0, height + thicknessOfLayers),
                    });
                }
            }
            foreach (var p in cell.mesh.finiteElements.Values.Where(r => r.backGrid != null))
            {
                List<Position> corners = new List<Position>();
                foreach (var c in p.corners)
                    if (!corners.Any(g => g.SameAs(c.position, cell.meshingAlgorithm.defaultTolerance)))
                        if (!corners.Any(cor => (float)cor.x == (float)c.position.x && (float)cor.y == (float)c.position.y)) // check if also not already present in float precision
                            corners.Add(c.position);

                objects.Add(new Prism
                {
                    Side = corners.Select(c => new Vector2F((float)(c.x * multiplicatorXYaxis), (float)(c.y * multiplicatorXYaxis))).ToArray(),
                    BottomToTopVector = new Vector3F(0, 0, thicknessOfLayers),
                    Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                    Color = new Color4(200, 200, 200, 255),
                    Name = "Grid",
                    Transform = Matrix4F.Translation(0, 0, height - thicknessOfLayers),
                });
            }*/
            
            // plot REGIONS with grid
            foreach (var region in cell.meshingAlgorithm.regions.Where(r => r.frontGrid.ID != 990000000))
            {
                objects.Add(new Prism
                {
                    // 2D-contours of the grid regions
                    Side = region.orderedPoints.Select(p => new Vector2F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis)).ToArray(),
                    // Vector that define translate between top and bottom sides
                    BottomToTopVector = new Vector3F(0, 0, thicknessOfLayers),
                    Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                    Color = new Color4(200, 200, 200, 255),
                    Name = "Grid",
                    Transform = Matrix4F.Translation(0, 0, height + thicknessOfLayers),
                });
            }

            foreach (var region in cell.meshingAlgorithm.regions.Where(r => r.backGrid.ID != 990000000))
            {
                objects.Add(new Prism
                {
                    // 2D-contours of the grid regions
                    Side = region.orderedPoints.Select(p => new Vector2F((float)p.position.x * multiplicatorXYaxis, (float)p.position.y * multiplicatorXYaxis)).ToArray(),
                    // Vector that define translate between top and bottom sides
                    BottomToTopVector = new Vector3F(0, 0, thicknessOfLayers),
                    Material = new RenderMaterial(0.7f, 0.1f, 0f, 0.3f, 0f),
                    Color = new Color4(200, 200, 200, 255),
                    Name = "Grid",
                    Transform = Matrix4F.Translation(0, 0, height - thicknessOfLayers),
                });
            }

            return new CompositeRenderData(objects) { Name = name, IsVisible = visible };
        }

        // Special mesh plot functions ██████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// create voronoi shell plot
        /// </summary>
        public static CompositeRenderData PlotVoronoiEdges<T>(string name, bool visible, Mesh<T> mesh, float multiplicatorXYaxis)
            where T : FiniteElement
        {
            // create polygons
            ObservableCollection<RenderData> polygonLine = new ObservableCollection<RenderData>();
            foreach (var point in mesh.finiteElements.Values)
            {
                // read corners of voronoi shells
                Vector3F[] corners = new Vector3F[point.corners.Count + 1];
                if (point.corners.Count > 0)
                    for (int e = 0; e <= point.corners.Count; e++)
                    {
                        corners[e] = new Vector3F((float)(multiplicatorXYaxis * point.corners[mod(e, point.corners.Count)].position.x),
                            (float)(multiplicatorXYaxis * point.corners[mod(e, point.corners.Count)].position.y), 0);
                    }

                // edges of polygons
                polygonLine.Add(new AtomicusChart.Interface.PresentationData.Line
                {
                    Points = corners,
                    Thickness = 2,
                    Color = Colors.Black
                });
            }

            return new CompositeRenderData(polygonLine) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create voronoi area plot
        /// </summary>
        public static CompositeRenderData PlotVoronoiAreas<T>(string name, bool visible, Mesh<T> mesh, float multiplicatorXYaxis)
            where T : FiniteElement
        {
            // create new random class for coloring the polygons in different gray-colors
            var rand = new Random();
            // create polygons
            ObservableCollection<RenderData> polygonshape = new ObservableCollection<RenderData>();
            foreach (var point in mesh.finiteElements.Values)
            {
                // fill polygons
                byte color = Convert.ToByte(rand.Next(0, 255));
                polygonshape.Add(new PolygoneXYArea(point.corners.Select(p => (p.position.x * multiplicatorXYaxis, p.position.y * multiplicatorXYaxis)).ToArray(), 0)
                {
                    Color = new Color4(color, color, color, 100), // Color with transparency
                });
            }

            return new CompositeRenderData(polygonshape) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create triangle plot
        /// </summary>
        public static CompositeRenderData PlotNeighbors<T>(string name, bool visible, Mesh<T> mesh, float multiplicatorXYaxis)
            where T : FiniteElement, new()
        {
            if (mesh == null)
                return new CompositeRenderData(new ObservableCollection<RenderData>()) { };

            // create triangles
            ObservableCollection<RenderData> lineList = new ObservableCollection<RenderData>();
            foreach (var element in mesh.finiteElements.Values)
            {
                foreach (var n in element.neighbors.Where(n => n.index > element.index))
                {
                    Vector3F[] cornerpoints = new Vector3F[2];
                    cornerpoints[0] = new Vector3F((float)element.position.x * multiplicatorXYaxis, (float)element.position.y * multiplicatorXYaxis, 0);
                    cornerpoints[1] = new Vector3F((float)mesh.finiteElements[n.index].position.x * multiplicatorXYaxis, (float)mesh.finiteElements[n.index].position.y * multiplicatorXYaxis, 0);

                    lineList.Add(new AtomicusChart.Interface.PresentationData.Line
                    {
                        Points = cornerpoints,
                        Thickness = 1,
                        Color = Colors.DarkBlue
                    });
                }
            }

            return new CompositeRenderData(lineList) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create region edges plot
        /// </summary>
        public static CompositeRenderData PlotRegionLines<T, R>(string name, bool visible, MeshingAlgorithm<T, R> meshingAlgorithm, float multiplicatorXYaxis)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            ObservableCollection<RenderData> contours = new ObservableCollection<RenderData>();
            foreach (var segment in meshingAlgorithm.contourSegments)
            {
                // read contour points
                Vector3F[] corners = new Vector3F[2];
                corners[0] = new Vector3F((float)(multiplicatorXYaxis * segment.firstAdjacentContourJunction.position.x), (float)(multiplicatorXYaxis * segment.firstAdjacentContourJunction.position.y), 0);
                corners[1] = new Vector3F((float)(multiplicatorXYaxis * segment.secondAdjacentContourJunction.position.x), (float)(multiplicatorXYaxis * segment.secondAdjacentContourJunction.position.y), 0);
                // add each line
                contours.Add(new AtomicusChart.Interface.PresentationData.Line
                {
                    Points = corners,
                    Thickness = 2,
                    Color = Colors.DarkGreen
                });
            }

            var contourslines = new CompositeRenderData(contours) { Name = name, IsVisible = visible };
            return contourslines;
        }
        /// <summary>
        /// create region edge indexes plot
        /// </summary>
        public static CompositeRenderData PlotRegionLineIndexes<T, R>(string name, bool visible, MeshingAlgorithm<T, R> meshingAlgorithm, float multiplicatorXYaxis)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            var labels = new ObservableCollection<RenderData>();

            foreach (var segment in meshingAlgorithm.contourSegments)
            {
                Position pos = segment.firstAdjacentContourJunction.position.CenterWith(segment.secondAdjacentContourJunction.position) * multiplicatorXYaxis;
                labels.Add(new Label
                {
                    Text = Convert.ToString(segment.index),
                    FontFamily = "Arial",
                    FontSize = 14,
                    Transform = Matrix4F.Translation((float)pos.x, (float)pos.y, 0),
                    Background = Colors.DarkGreen,
                    MarkerColor = Colors.DarkGreen,
                    FontColor = Colors.White,
                });
            }

            return new CompositeRenderData(labels) { Name = name, IsVisible = visible };
        }
        /// <summary>
        /// create forbidden area plot for all contour junctions
        /// </summary>
        public static CompositeRenderData PlotForbiddenAreaJunctions<T, R>(string name, bool visible, MeshingAlgorithm<T, R> meshingAlgorithm, float multiplicatorXYaxis)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            if (meshingAlgorithm.contourJunctions.FirstOrDefault().radiusForbiddenArea < 0)
                return new CompositeRenderData(new ObservableCollection<RenderData>()) { Name = name, IsVisible = visible };

            ObservableCollection<RenderData> plotData = new ObservableCollection<RenderData>();
            foreach (var cs in meshingAlgorithm.contourSegments)
            {
                // circle around first contour junctions ────────────────────────────────────────────────────────────────────────────────────────────
                plotData.Add(new CircleXYLine(cs.firstAdjacentContourJunction.position.x * multiplicatorXYaxis, cs.firstAdjacentContourJunction.position.y * multiplicatorXYaxis,
                    meshingAlgorithm.dPointToJunctionGlobal * multiplicatorXYaxis, 0, 3600)
                {
                    Name = "circle with length dPoint", // Name

                    PatternStyle = PatternStyle.Solid, // Type of line
                    Thickness = 1f, // Thickness of line
                    Color = new Color4(0, 0, 0, 255), // Color of line

                    MarkerStyle = MarkerStyle.None, // Type of markers
                    MarkerSize = 5, // Size of markers
                    MarkerColor = new Color4(0, 0, 0, 255), // Color of markers
                });
                plotData.Add(new CircleXYArea(cs.firstAdjacentContourJunction.position.x * multiplicatorXYaxis, cs.firstAdjacentContourJunction.position.y * multiplicatorXYaxis,
                    meshingAlgorithm.dPointToJunctionGlobal * multiplicatorXYaxis, 0, 3600)
                {
                    Name = "circle with length dPoint", // Name
                    Color = new Color4(200, 0, 0, 100), // Color of filled area
                });

                // circle around second contour junctions ───────────────────────────────────────────────────────────────────────────────────────────
                plotData.Add(new CircleXYLine(cs.secondAdjacentContourJunction.position.x * multiplicatorXYaxis, cs.secondAdjacentContourJunction.position.y * multiplicatorXYaxis,
                    meshingAlgorithm.dPointToJunctionGlobal * multiplicatorXYaxis, 0, 3600)
                {
                    Name = "circle with length dPoint", // Name

                    PatternStyle = PatternStyle.Solid, // Type of line
                    Thickness = 1f, // Thickness of line
                    Color = new Color4(0, 0, 0, 255), // Color of line

                    MarkerStyle = MarkerStyle.None, // Type of markers
                    MarkerSize = 5, // Size of markers
                    MarkerColor = new Color4(0, 0, 0, 255), // Color of markers
                });
                plotData.Add(new CircleXYArea(cs.secondAdjacentContourJunction.position.x * multiplicatorXYaxis, cs.secondAdjacentContourJunction.position.y * multiplicatorXYaxis,
                    meshingAlgorithm.dPointToJunctionGlobal * multiplicatorXYaxis, 0, 3600)
                {
                    Name = "circle with length dPoint", // Name
                    Color = new Color4(200, 0, 0, 100), // Color of filled area
                });
            }

            var Hilfslinien = new CompositeRenderData(plotData) { Name = name, IsVisible = visible };
            return Hilfslinien;
        }
        /// <summary>
        /// create forbidden area for all contour segments
        /// </summary>
        public static CompositeRenderData PlotForbiddenAreaSegments<T, R>(string name, bool visible, MeshingAlgorithm<T, R> meshingAlgorithm, float multiplicatorXYaxis)
            where T : FiniteElement, new()
            where R : Region, new()
        {
            if (meshingAlgorithm.contourSegments.FirstOrDefault().forbiddenArea == null)
                return new CompositeRenderData(new ObservableCollection<RenderData>()) { Name = name, IsVisible = visible };

            ObservableCollection<RenderData> plotData = new ObservableCollection<RenderData>();
            foreach (var cs in meshingAlgorithm.contourSegments)
            {
                // forbidden area ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                plotData.Add(new PolygoneXYLine(cs.forbiddenArea.Select(p => (p.x * multiplicatorXYaxis, p.y * multiplicatorXYaxis)).ToArray(), 0)
                {
                    Name = "forbidden area", // Name

                    PatternStyle = PatternStyle.Solid, // Type of line
                    Thickness = 1f, // Thickness of line
                    Color = new Color4(0, 0, 0, 255), // Color of line

                    MarkerStyle = MarkerStyle.None, // Type of markers
                    MarkerSize = 5, // Size of markers
                    MarkerColor = new Color4(0, 0, 0, 255), // Color of markers
                });
                plotData.Add(new PolygoneXYArea(cs.forbiddenArea.Select(p => (p.x * multiplicatorXYaxis, p.y * multiplicatorXYaxis)).ToArray(), 0)
                {
                    Name = "forbidden area", // Name
                    Color = new Color4(200, 0, 0, 100), // Color of filled area
                });
            }

            var auxiliarLines = new CompositeRenderData(plotData) { Name = name, IsVisible = visible };
            return auxiliarLines;
        }

        // Special current plot functions ███████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// plot P2 current arrows
        /// </summary>
        public static CompositeRenderData PlotP2Arrows(string name, bool visible, ModelCell cell, float height, Color4 color, float multiplicatorXYaxis)
        {
            if (cell.mesh.finiteElements.Values.Any(p => p.type == pointType.P2))
            {
                float maxLengthP2 = (float)cell.mesh.finiteElements.Values.Where(p => p.type == pointType.P2).Max(p => p.GetCurrentP2() / p.size);

                ObservableCollection<RenderData> arrows = new ObservableCollection<RenderData>();
                foreach (var point in cell.mesh.finiteElements.Values.Where(p => p.type == pointType.P2))
                {
                    // add all triangles
                    arrows.Add(new SingleColorPrimitiveCollection(new Matrix4F[]{ Matrix4F.Scaling(0.1f)
                * Matrix4F.RotationAxis(Vector3F.UnitY, Math.PI / 2)
                * Matrix4F.Scaling(1f, 1f, height / 0.1f * (float)(point.GetCurrentP2() / (point.size * maxLengthP2)))
                * Matrix4F.Translation(new Vector3F((float)(point.position.x * multiplicatorXYaxis), (float)(point.position.y * multiplicatorXYaxis), height)) })
                    {
                        Color = color,
                        Mesh = ArrowMeshFactory.GenerateArrowX(20, 0.5f, 0.8f * height / 0.1f * (float)(point.GetCurrentP2() / (point.size * maxLengthP2))),
                    });
                }

                return new CompositeRenderData(arrows) { Name = name, IsVisible = visible };
            }
            else
                return new CompositeRenderData(new ObservableCollection<RenderData>());
        }
        /// <summary>
        /// plot generated current arrows
        /// </summary>
        public static CompositeRenderData PlotGeneratedArrows(string name, bool visible, ModelCell cell, float height, Color4 colorGenerated, Color4 colorShunt, float multiplicatorXYaxis)
        {
            float maxLengthGen = (float)cell.mesh.finiteElements.Values.Max(p => Math.Abs(p.GetCurrentGenerated()) / p.size);

            ObservableCollection<RenderData> arrows = new ObservableCollection<RenderData>();
            foreach (var point in cell.mesh.finiteElements.Values.Where(p => p.type != pointType.P2))
            {
                if (point.GetCurrentGenerated() <= 0)
                {
                    arrows.Add(new SingleColorPrimitiveCollection(new Matrix4F[]{ Matrix4F.Scaling(0.1f)
                        * Matrix4F.RotationAxis(Vector3F.UnitY, -Math.PI / 2)
                        * Matrix4F.Scaling(1f, 1f, height / 0.1f * (float)(-point.GetCurrentGenerated() / (point.size * maxLengthGen)))
                        * Matrix4F.Translation(new Vector3F((float)(point.position.x * multiplicatorXYaxis), (float)(point.position.y * multiplicatorXYaxis), 0)) })
                    {
                        Color = colorGenerated,
                        Mesh = ArrowMeshFactory.GenerateArrowX(20, 0.5f, 0.8f * height / 0.1f * (float)(-point.GetCurrentGenerated() / (point.size * maxLengthGen))),
                    });
                }
                else
                {
                    arrows.Add(new SingleColorPrimitiveCollection(new Matrix4F[]{ Matrix4F.Scaling(0.1f)
                        * Matrix4F.RotationAxis(Vector3F.UnitY, Math.PI / 2)
                        * Matrix4F.Scaling(1f, 1f, height / 0.1f * (float)(point.GetCurrentGenerated() / (point.size * maxLengthGen)))
                        * Matrix4F.Translation(new Vector3F((float)(point.position.x * multiplicatorXYaxis), (float)(point.position.y * multiplicatorXYaxis), height)) })
                    {
                        Color = colorShunt,
                        Mesh = ArrowMeshFactory.GenerateArrowX(20, 0.5f, 0.8f * height / 0.1f * (float)(point.GetCurrentGenerated() / (point.size * maxLengthGen))),
                    });
                }
            }

            return new CompositeRenderData(arrows) { Name = name, IsVisible = visible };
        }

        // Modulo function ██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████
        /// <summary>
        /// Modulo function, which also workes for negative numbers
        /// </summary>
        /// <param name="x">numerator</param>
        /// <param name="m">divisor</param>
        /// <returns>x mod m</returns>
        public static int mod(int x, int m)
        {
            int r = x % m;
            return r < 0 ? r + m : r;
        }
    }

    // Listener and Geometry classes ████████████████████████████████████████████████████████████████████████████████████████████████████████████████
    /// <summary>
    /// Class for contour plots
    /// </summary>
    class CustomContoursOwner : IContoursOwner
    {
        public CustomContoursOwner(ObservableCollection<Contour> contours) => Contours = contours;

        public event PropertyChangedEventHandler PropertyChanged;

        public event PropertyChangingEventHandler PropertyChanging;

        public ObservableCollection<Contour> Contours { get; }
    }
    /// <summary>
    /// Class for drawing contour lines of a circle in the xy-plane
    /// </summary>
    public class CircleXYLine : Series
    {
        public CircleXYLine(double centerPositionX, double centerPositionY, double radius, double height, int amountDrawnPoints)
        {
            Reader = new DefaultPositionMaskDataReader(CreateLinesStrip((float)radius,
                new Vector3F((float)centerPositionX, (float)centerPositionY, (float)height), amountDrawnPoints));
        }
        private static Vector3F[] CreateLinesStrip(float r, Vector3F pos, int resolution)
        {
            var res = new Vector3F[resolution];
            var start = r * Vector3F.UnitY;
            for (int i = 0; i < resolution; i++)
            {
                res[i] = pos + start * Matrix4F.RotationAxis(Vector3F.UnitZ, (double)i / (resolution - 1) * 2 * Math.PI);
            }
            return res;
        }
    }
    /// <summary>
    /// Class for drawing filled area of a circle in the xy-plane
    /// </summary>
    public class CircleXYArea : Surface
    {
        public CircleXYArea(double centerPositionX, double centerPositionY, double radius, double height, int amountDrawnPoints)
        {
            // read corners
            Vector3F[] corners = new Vector3F[amountDrawnPoints];
            for (int i = 0; i < amountDrawnPoints; i++)
                corners[i] = new Vector3F((float)centerPositionX, (float)centerPositionY, (float)height)
                    + (float)radius * Vector3F.UnitY * Matrix4F.RotationAxis(Vector3F.UnitZ, (double)i / (amountDrawnPoints - 1) * 2 * Math.PI);

            // create an array like { 0,1,2, 0,2,3, 0,3,4, 0,4,5, ... }
            int[] polygonTriangleIndexes = new int[(corners.Length - 3) * 3];
            for (int i = 0; i < polygonTriangleIndexes.Length; i++)
            {
                if (Plotter.mod(i, 3) == 0)
                    polygonTriangleIndexes[i] = 0;
                if (Plotter.mod(i, 3) == 1)
                    polygonTriangleIndexes[i] = (i + 2) / 3;
                if (Plotter.mod(i, 3) == 2)
                    polygonTriangleIndexes[i] = (i + 4) / 3;
            }

            var polygonLineIndexes = new int[] { 0 };

            //Alternate indexes are used for different presentations like SolidAndWireframe or Wireframe.
            SurfaceMesh = new Mesh(corners, polygonTriangleIndexes, polygonLineIndexes);
        }
    }
    /// <summary>
    /// Class for drawing contour lines of a polygon in the xy-plane
    /// </summary>
    public class PolygoneXYLine : Series
    {
        public PolygoneXYLine((double x, double y)[] corners, double height)
        {
            Vector3F[] outerPointsWithLast = new Vector3F[corners.Length + 1];
            for (int i = 0; i < corners.Length; i++)
                outerPointsWithLast[i] = new Vector3F((float)corners[i].x, (float)corners[i].y, (float)height);
            outerPointsWithLast[outerPointsWithLast.Length - 1] = new Vector3F((float)corners[0].x, (float)corners[0].y, (float)height);

            Reader = new DefaultPositionMaskDataReader(outerPointsWithLast);
        }
    }
    /// <summary>
    /// Class for drawing filled area of a polygon in the xy-plane
    /// </summary>
    public class PolygoneXYArea : Surface
    {
        public PolygoneXYArea((double x, double y)[] corners, double height)
        {
            // read corners
            Vector3F[] outerPoints = new Vector3F[corners.Length + 1];
            for (int e = 0; e <= corners.Length; e++)
            {
                outerPoints[e] = new Vector3F((float)corners[Plotter.mod(e, corners.Length)].x,
                    (float)corners[Plotter.mod(e, corners.Length)].y,
                    (float)height);
            }

            // create an array like { 0,1,2, 0,2,3, 0,3,4, 0,4,5, ... }
            int[] polygonTriangleIndexes = new int[(corners.Length - 2) * 3];
            for (int i = 0; i < polygonTriangleIndexes.Length; i++)
            {
                if (Plotter.mod(i, 3) == 0)
                    polygonTriangleIndexes[i] = 0;
                if (Plotter.mod(i, 3) == 1)
                    polygonTriangleIndexes[i] = (i + 2) / 3;
                if (Plotter.mod(i, 3) == 2)
                    polygonTriangleIndexes[i] = (i + 4) / 3;
            }

            var polygonLineIndexes = new int[] { 0 };

            //Alternate indexes are used for different presentations like SolidAndWireframe or Wireframe.
            SurfaceMesh = new Mesh(outerPoints, polygonTriangleIndexes, polygonLineIndexes);
        }
    }
}