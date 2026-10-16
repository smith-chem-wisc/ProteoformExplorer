using Chemistry;
using ProteoformExplorer.Deconvoluter;
using MassSpectrometry;
using MzLibUtil;
using ProteoformExplorer.Core;
using ScottPlot;
using ScottPlot.Drawing;
using ScottPlot.Plottable;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Drawing;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using MathNet.Numerics;

namespace ProteoformExplorer.GuiFunctions
{
    public class PlottingFunctions
    {
        public static void StylePlot(Plot plot)
        {
            plot.Palette = new Palette(GuiSettings.ChartColorPalette);
            plot.Grid(GuiSettings.ShowChartGrid);
            plot.XAxis.Ticks(major: true, minor: false);
            plot.YAxis.Ticks(major: true, minor: false);
            plot.Legend(location: GuiSettings.LegendLocation);
            plot.XAxis.TickMarkColor(Color.Transparent);
            plot.YAxis.TickMarkColor(Color.Transparent);

            plot.XAxis.Line(width: (float)GuiSettings.DpiScalingY);
            plot.XAxis.Label(size: (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingY)
                //    , bold: true
                );
            plot.XAxis.TickLabelStyle(fontSize: (float)(GuiSettings.ChartTickFontSize * GuiSettings.DpiScalingY)
                //    , fontBold: true
                );

            plot.YAxis.Line(width: (float)GuiSettings.DpiScalingX);
            plot.YAxis.Label(size: (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingX)
                //    , bold: true
                );
            plot.YAxis.TickLabelStyle(fontSize: (float)(GuiSettings.ChartTickFontSize * GuiSettings.DpiScalingX)
                //    , fontBold: true
                );

            var legend = plot.Legend(location: GuiSettings.LegendLocation);
            legend.FontSize = (float)(GuiSettings.ChartLegendFontSize * GuiSettings.DpiScalingX);

            // title
            plot.XAxis2.Label(size: (float)(GuiSettings.ChartHeaderFontSize * GuiSettings.DpiScalingX));

            foreach (var item in plot.GetPlottables())
            {
                if (item is ScatterPlot scatter)
                {
                    //TODO: this will affect all lines... spectra and xic..
                    scatter.LineWidth = GuiSettings.ChartLineWidth * GuiSettings.DpiScalingY;
                }
                else if (item is LollipopPlot lollipopPlot)
                {
                    //TODO
                    lollipopPlot.ErrorLineWidth = (float)(GuiSettings.ChartLineWidth * GuiSettings.DpiScalingX);
                    lollipopPlot.BarWidth = (float)(GuiSettings.ChartLineWidth * GuiSettings.DpiScalingX);
                }
            }

            plot.Benchmark(false);
        }

        public static void StyleDashboardPlot(Plot plot)
        {
            StylePlot(plot);

            plot.XAxis.TickLabelStyle(rotation: (float)GuiSettings.XLabelRotation);
            plot.Frame(top: false, right: false);
            plot.YAxis.Line(visible: false);
        }

        public static void StyleIonChromatogramPlot(Plot plot, VLine rtIndicator)
        {
            if (rtIndicator != null)
            {
                plot.Add(rtIndicator);
            }

            StylePlot(plot);

            plot.YAxis.TickLabelNotation(multiplier: true);
            plot.YAxis.Label("Intensity");
            plot.XAxis.Label("Retention Time");
            plot.Layout(padding: 10);
        }

        public static void StyleSpectrumPlot(Plot plot, MsDataScan scan, string dataFileName)
        {
            StylePlot(plot);

            try
            {
                plot.Title(dataFileName + "  Scan#" + scan.OneBasedScanNumber + "  RT: " + scan.RetentionTime + "  MS" + scan.MsnOrder + "\n" +
                    " " + scan.ScanFilter, bold: false);
            }
            catch (Exception e)
            {
                plot.Title(dataFileName + "  Scan#" + scan.OneBasedScanNumber + "  RT: " + scan.RetentionTime + "  MS" + scan.MsnOrder, bold: false);
            }

            plot.YAxis.TickLabelNotation(multiplier: true);
            plot.YAxis.Label("Intensity");
            plot.XAxis.Label("m/z");
            plot.SetViewLimits(scan.ScanWindowRange.Minimum, scan.ScanWindowRange.Maximum, 0, scan.MassSpectrum.YArray.Max() * 2.0);
        }

        public static void PlotSpeciesInSpectrum(List<AnnotatedSpecies> allSpeciesToPlot, int oneBasedScan, KeyValuePair<string, CachedSpectraFileData> data,
            Plot spectrumPlot, out MsDataScan scan, out List<string> errors, int? charge = null)
        {
            scan = null;
            errors = new List<string>();

            try
            {
                spectrumPlot.Clear();

                if (data.Value == null)
                {
                    return;
                }

                // get the scan
                scan = data.Value.GetOneBasedScan(oneBasedScan);

                if (scan == null)
                {
                    return;
                }

                // add non-annotated peaks
                List<Datum> spectrumData = new List<Datum>();

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    spectrumData.Add(new Datum(scan.MassSpectrum.XArray[i], scan.MassSpectrum.YArray[i]));
                }

                foreach (var item in spectrumData)
                {
                    spectrumPlot.AddLine(item.X, 0, item.X, item.Y.Value, GuiSettings.UnannotatedSpectrumColor, (float)GuiSettings.ChartLineWidth);
                }

                if (allSpeciesToPlot != null && scan.MsnOrder == 1)
                {
                    Dictionary<string, Color> speciesLabelToColor = new Dictionary<string, Color>();

                    // add annotated peaks
                    HashSet<double> claimedMzs = new HashSet<double>();
                    foreach (var species in allSpeciesToPlot)
                    {
                        if (!speciesLabelToColor.TryGetValue(species.SpeciesLabel, out var color))
                        {
                            speciesLabelToColor.Add(species.SpeciesLabel, spectrumPlot.GetNextColor());
                            color = speciesLabelToColor[species.SpeciesLabel];
                        }

                        List<Datum> annotatedData = new List<Datum>();
                        List<int> chargesToPlot = new List<int>();

                        if (species.DeconvolutionFeature != null)
                        {
                            foreach (var envelope in species.DeconvolutionFeature.AnnotatedEnvelopes)
                            {
                                int added = 0;

                                if (charge != null && charge.Value != envelope.Charge)
                                {
                                    continue;
                                }

                                foreach (var mz in envelope.PeakMzs)
                                {
                                    int index = scan.MassSpectrum.GetClosestPeakIndex(mz);
                                    double expMz = scan.MassSpectrum.XArray[index];
                                    double expIntensity = scan.MassSpectrum.YArray[index];

                                    if (!claimedMzs.Contains(expMz))
                                    {
                                        annotatedData.Add(new Datum(expMz, expIntensity));
                                        added++;
                                        claimedMzs.Add(expMz);
                                    }
                                }

                                if (added == 0)
                                {
                                    var firstMz = envelope.PeakMzs.First();

                                    int index = scan.MassSpectrum.GetClosestPeakIndex(firstMz);
                                    double expMz = scan.MassSpectrum.XArray[index];
                                    double expIntensity = scan.MassSpectrum.YArray[index];

                                    annotatedData.Add(new Datum(expMz, expIntensity));
                                }
                            }
                        }

                        foreach (var item in annotatedData)
                        {
                            spectrumPlot.AddLine(item.X, 0, item.X, item.Y.Value, color, (float)GuiSettings.AnnotatedEnvelopeLineWidth);
                        }
                    }
                }

                StyleSpectrumPlot(spectrumPlot, scan, data.Value.FileName);
                spectrumPlot.SetAxisLimits(scan.ScanWindowRange.Minimum, scan.ScanWindowRange.Maximum, 0, scan.MassSpectrum.YArray.Max() * 1.2);
            }
            catch (Exception e)
            {
                errors.Add("An error occurred while plotting the spectrum: " + e.Message);
            }
        }

        public static void PlotXic(double mz, int z, Tolerance tolerance, double rt, double rtWindow, KeyValuePair<string, CachedSpectraFileData> data, Plot xicPlot,
            bool clearOldPlot, VLine rtIndicator, out List<string> errors, double xOffset = 0, double yOffset = 0, string label = null)
        {
            errors = new List<string>();

            try
            {
                if (clearOldPlot)
                {
                    xicPlot.Clear();
                }

                var scans = data.Value.GetScansInRtWindow(rt, rtWindow);

                var xicData = PfmXplorerUtil.GetXicData(scans, mz, z, tolerance);

                var xs = xicData.Select(p => p.X).ToArray();
                var ys = xicData.Select(p => p.Y.Value).ToArray();
                var color = xicPlot.GetNextColor();

                xicPlot.AddScatterLines(xs, ys, color, label: label);
                xicPlot.Legend(label != null, location: GuiSettings.LegendLocation);

                bool isWaterfall = xOffset != 0 || yOffset != 0;

                if ((isWaterfall && GuiSettings.FillWaterfall) || (!isWaterfall && GuiSettings.FillSideView))
                {
                    var colorWithTransparency = Color.FromArgb(GuiSettings.FillAlpha, color.R, color.G, color.B);
                    xicPlot.AddFill(xs, ys, baseline: yOffset, color: colorWithTransparency);
                }

                if (clearOldPlot)
                {
                    xicPlot.SetAxisLimitsY(Math.Min(0, ys.Min()), ys.Max());
                    xicPlot.SetAxisLimitsX(xs.Min(), xs.Max());
                }
                else
                {
                    var axisLimits = xicPlot.GetAxisLimits(xicPlot.XAxis.AxisIndex, xicPlot.YAxis.AxisIndex);
                    double yMin = Math.Min(axisLimits.YMin, ys.Min());
                    double yMax = Math.Max(axisLimits.YMax, ys.Max());
                    double xMin = Math.Min(axisLimits.XMin, xs.Min());
                    double xMax = Math.Max(axisLimits.XMax, xs.Max());

                    xicPlot.SetAxisLimitsY(yMin, yMax);
                    xicPlot.SetAxisLimitsX(xMin, xMax);

                    //xicPlot.Plot.SetViewLimits(xMin, xMax, yMin, yMax * 2);
                }

                StyleIonChromatogramPlot(xicPlot, rtIndicator);
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void PlotSummedChargeStateXic(double modeMass, int z, double rt, double rtWindow, KeyValuePair<string, CachedSpectraFileData> data, Plot xicPlot,
            bool clearOldPlot, VLine rtIndicator, out List<string> errors, double xOffset = 0, double yOffset = 0, bool waterfall = false, string label = null)
        {
            errors = new List<string>();

            try
            {

                if (clearOldPlot)
                {
                    xicPlot.Clear();
                }

                var scans = data.Value.GetScansInRtWindow(rt, rtWindow);

                var xicData = PfmXplorerUtil.GetSummedChargeXics(scans, modeMass, z);

                var xs = xicData.Select(p => p.X + xOffset).ToArray();
                var ys = xicData.Select(p => p.Y.Value + yOffset).ToArray();
                var color = xicPlot.GetNextColor();

                xicPlot.AddScatterLines(xs, ys, color, lineWidth: (float)GuiSettings.ChartLineWidth,
                    label: waterfall ? null : label);

                if (waterfall && !string.IsNullOrWhiteSpace(label))
                {
                    var text = new Text
                    {
                        Label = label,
                        X = xs.Min(),
                        Y = yOffset,
                        Color = color,
                        FontSize = (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingX),
                        Alignment = Alignment.MiddleLeft
                    };

                    xicPlot.Add(text);
                }

                xicPlot.Legend(!waterfall && label != null, location: GuiSettings.LegendLocation);

                if ((waterfall && GuiSettings.FillWaterfall) || (!waterfall && GuiSettings.FillSideView))
                {
                    var colorWithTransparency = Color.FromArgb(GuiSettings.FillAlpha, color.R, color.G, color.B);
                    xicPlot.AddFill(xs, ys, baseline: yOffset, color: colorWithTransparency);
                }

                if (clearOldPlot)
                {
                    xicPlot.SetAxisLimitsY(Math.Min(0, ys.Min()), ys.Max());
                    xicPlot.SetAxisLimitsX(xs.Min(), xs.Max());
                }
                else
                {
                    var axisLimits = xicPlot.GetAxisLimits(xicPlot.XAxis.AxisIndex, xicPlot.YAxis.AxisIndex);
                    double yMin = Math.Min(axisLimits.YMin, ys.Min());
                    double yMax = Math.Max(axisLimits.YMax, ys.Max());
                    double xMin = Math.Min(axisLimits.XMin, xs.Min());
                    double xMax = Math.Max(axisLimits.XMax, xs.Max());

                    xicPlot.SetAxisLimitsY(yMin, yMax);
                    xicPlot.SetAxisLimitsX(xMin, xMax);

                    //xicPlot.Plot.SetViewLimits(xMin, xMax, yMin, yMax * 2);
                }

                StyleIonChromatogramPlot(xicPlot, rtIndicator);

                if (waterfall)
                {
                    xicPlot.Legend(false);
                }
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void PlotTotalIonChromatograms(Plot plot, VLine rtIndicator, out List<string> errors, double? minRt = null, double? maxRt = null)
        {
            errors = new List<string>();
            plot.Clear();

            try
            {
                if (DataManagement.CurrentlySelectedFile.Value == null)
                {
                    return;
                }

                // display TIC chromatogram
                var ticChromatogram = DataManagement.CurrentlySelectedFile.Value.GetTicChromatogram(GuiSettings.TicRollingAverage, minRt, maxRt);

                plot.AddScatterLines(
                    ticChromatogram.Select(p => p.X).ToArray(),
                    ticChromatogram.Select(p => p.Y.Value).ToArray(),
                    GuiSettings.TicColor, (float)GuiSettings.ChartLineWidth, label: "TIC");

                // display identified TIC chromatogram
                if (DataManagement.AllLoadedAnnotatedSpecies.Any())
                {
                    var deconvolutedTicChromatogram = DataManagement.CurrentlySelectedFile.Value.GetDeconvolutedTicChromatogram(GuiSettings.TicRollingAverage, minRt, maxRt);

                    if (deconvolutedTicChromatogram.Any())
                    {
                        plot.AddScatterLines(
                            deconvolutedTicChromatogram.Select(p => p.X).ToArray(),
                            deconvolutedTicChromatogram.Select(p => p.Y.Value).ToArray(),
                            GuiSettings.DeconvolutedColor, (float)GuiSettings.ChartLineWidth, label: "Deconvoluted TIC");
                    }

                    var identifiedTicChromatogram = DataManagement.CurrentlySelectedFile.Value.GetIdentifiedTicChromatogram(GuiSettings.TicRollingAverage, minRt, maxRt);

                    if (identifiedTicChromatogram.Any())
                    {
                        var nonIdentifiedValue = identifiedTicChromatogram.Where(p => p.Label is null).ToArray();
                        bool oneType = identifiedTicChromatogram.DistinctBy(p => p.Label).Count() == 1;
                        foreach (var labelSet in identifiedTicChromatogram.GroupBy(p => p.Label))
                        {
                            //if (labelSet.Key is null)
                            //    continue;

                            var toPlot = nonIdentifiedValue.Concat(labelSet).OrderBy(p => p.X);

                            plot.AddScatterLines(
                                toPlot.Select(p => p.X).ToArray(),
                                toPlot.Select(p => p.Y.Value).ToArray(),
                                labelSet.Key.ConvertStringToColor(), (float)GuiSettings.ChartLineWidth, label: oneType ? "Identified TIC" : labelSet.Key.ConvertName());
                        }
                    }
                }

                plot.SetViewLimits(xMin: ticChromatogram.Min(p => p.X), xMax: ticChromatogram.Max(p => p.X), yMin: 0, yMax: ticChromatogram.Max(p => p.Y.Value) * 2.0);

                StyleIonChromatogramPlot(plot, rtIndicator);
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void DrawPercentTicPerFileInfoDashboardPlot(Plot plot, out List<string> errors)
        {
            plot.Clear();
            errors = new List<string>();

            try
            {
                plot.Title(@"MS1 TIC");

                var ticValues = new List<(string file, double tic, double deconvolutedTic, double identifiedTic, string label)>();


                var identifiedTicDict = DataManagement.AllLoadedAnnotatedSpecies
                    .Where(p => p.Identification != null)
                    .Select(p => p.Identification.Dataset)
                    .Where(p => p != null)
                    .Distinct()
                    .ToDictionary(p => p, p => 0.0);

                bool multipleInputRuns = identifiedTicDict.Count > 1;
                foreach (var file in DataManagement.SpectraFiles.OrderBy(p => p.Key.ConvertName()))
                {
                    double tic = 0;
                    double deconvolutedTic = 0;
                    double identifiedTic = 0;
                    identifiedTicDict.ForEach(p => identifiedTicDict[p.Key] = 0);

                    var ticChromatogram = file.Value.GetTicChromatogram();
                    if (ticChromatogram != null)
                    {
                        tic = ticChromatogram.Sum(p => p.Y.Value);
                    }

                    var deconvolutedTicChromatogram = file.Value.GetDeconvolutedTicChromatogram();
                    if (deconvolutedTicChromatogram != null)
                    {
                        deconvolutedTic = deconvolutedTicChromatogram.Sum(p => p.Y.Value);
                    }

                    var identifiedTicChromatogram = file.Value.GetIdentifiedTicChromatogram();
                    if (identifiedTicChromatogram != null)
                    {
                        identifiedTicChromatogram.Where(p => p.Label != null).ForEach(id => identifiedTicDict[id.Label] += id.Y.Value);
                        foreach (var label in identifiedTicChromatogram.Select(m => m.Label).Distinct())
                        {
                            if (label is null)
                                continue;
                            ticValues.Add((PfmXplorerUtil.GetFileNameWithoutExtension(file.Key).ConvertName(), tic, deconvolutedTic, identifiedTicDict[label], label.ConvertName()));
                        }
                    }
                }

                double[] positions = Enumerable.Range(0, ticValues.DistinctBy(p => p.file).Count()).Select(p => (double)p).ToArray();
                string[] labels = ticValues.Select(p => p.file.ConvertName()).Distinct().ToArray();

                var ticPlot = new ScottPlot.Plottable.LollipopPlot(positions, ticValues.DistinctBy(p => p.file).Select(p => p.tic).ToArray());
                ticPlot.Label = @"TIC";
                ticPlot.LollipopColor = GuiSettings.TicColor;
                ticPlot.LollipopRadius = 10;
                plot.Add(ticPlot);

                var deconTicPlot = new ScottPlot.Plottable.LollipopPlot(positions, ticValues.DistinctBy(p => p.file).Select(p => p.deconvolutedTic).ToArray());
                deconTicPlot.Label = @"Deconvoluted TIC";
                deconTicPlot.LollipopColor = GuiSettings.DeconvolutedColor;
                deconTicPlot.LollipopRadius = 10;
                plot.Add(deconTicPlot);


                if (ticValues.Select(p => p.label).Distinct().Count() > 1)
                {
                    foreach (var labelGroup in ticValues.GroupBy(p => p.label))
                    {
                        var identTicPlot =
                            new LollipopPlot(positions, labelGroup.Select(p => p.identifiedTic).ToArray());
                        identTicPlot.Label = labelGroup.Key;
                        identTicPlot.LollipopColor = labelGroup.Key.ConvertStringToColor();
                        identTicPlot.LollipopRadius = 10;
                        plot.Add(identTicPlot);
                    }
                }
                else
                {
                    var identTicPlot = new ScottPlot.Plottable.LollipopPlot(positions, ticValues.Select(p => p.identifiedTic).ToArray());
                    identTicPlot.Label = @"Identified TIC";
                    identTicPlot.LollipopColor = GuiSettings.IdentifiedColor;
                    identTicPlot.LollipopRadius = 10;
                    plot.Add(identTicPlot);
                }
                

                plot.YAxis.Label("Intensity");
                plot.SetAxisLimitsY(ticValues.Max(p => p.tic) * -0.03, ticValues.Max(p => p.tic) * 1.4);
                plot.XTicks(positions, labels);
                plot.YAxis.TickLabelNotation(multiplier: true);

                StyleDashboardPlot(plot);
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void DrawNumEnvelopesDashboardPlot(Plot plot, out List<string> errors)
        {
            plot.Clear();
            errors = new List<string>();

            try
            {
                plot.Title(@"MS1 Envelope Counts");

                var numFilteredEnvelopesPerFile = new List<(string file, int numFilteredEnvs)>();

                foreach (var file in DataManagement.SpectraFiles.OrderBy(p => p.Key.ConvertName()))
                {
                    int envs = file.Value.GetDistinctEnvelopes().Count;
                    numFilteredEnvelopesPerFile.Add((PfmXplorerUtil.GetFileNameWithoutExtension(file.Key), envs));
                }

                double[] values = numFilteredEnvelopesPerFile.Select(p => (double)p.numFilteredEnvs).ToArray();
                double[] positions = Enumerable.Range(0, numFilteredEnvelopesPerFile.Count).Select(p => (double)p).ToArray();
                string[] labels = numFilteredEnvelopesPerFile.Select(p => p.file.ConvertName()).ToArray();

                var envsCountPlot = new ScottPlot.Plottable.LollipopPlot(positions, values);
                envsCountPlot.Label = @"Deconvoluted Envelopes";
                envsCountPlot.LollipopColor = GuiSettings.DeconvolutedColor;
                envsCountPlot.LollipopRadius = 10;
                plot.Add(envsCountPlot);

                plot.YAxis.Label("Count");
                plot.SetAxisLimitsY(0, numFilteredEnvelopesPerFile.Max(p => p.numFilteredEnvs) * 1.2);
                plot.XTicks(positions, labels);

                StyleDashboardPlot(plot);
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void DrawMassDistributionsDashboardPlot(Plot plot, out List<string> errors)
        {
            plot.Clear();
            errors = new List<string>();

            try
            {
                plot.Title(@"MS1 Envelope Mass Histograms");
                var color = GuiSettings.DeconvolutedColor;

                int fileNum = 0;
                double maxMass = 0;
                foreach (var file in DataManagement.SpectraFiles.OrderBy(p => p.Key.ConvertName()))
                {
                    //TODO: decon feature could be null
                    var envelopeMasses = file.Value
                        .GetDistinctEnvelopes()
                        .Select(v => v.PeakMzs.First().ToMass(v.Charge))
                        .ToArray();
                    maxMass = Math.Max(maxMass, envelopeMasses.Max());

                    // generate histogram bars
                    List<(double start, double end, int count)> histogramBins = new List<(double start, double end, int count)>();

                    for (int i = 0; i < int.MaxValue; i++)
                    {
                        double binMin = i * GuiSettings.MassHistogramBinWidth;
                        double binMax = (i + 1) * GuiSettings.MassHistogramBinWidth;

                        if (binMin > envelopeMasses.Max())
                        {
                            break;
                        }

                        int itemsInBin = envelopeMasses.Count(p => p >= binMin && p < binMax);
                        histogramBins.Add((binMin, binMax, itemsInBin));
                    }

                    double[] fractionalAbundance = histogramBins.Select(p => p.count / ((double)envelopeMasses.Length)).ToArray();
                    double[] locations = histogramBins.Select(p => (p.start + p.end) / 2.0).ToArray();

                    var bar = plot.AddBar(values: fractionalAbundance, positions: locations);
                    bar.BarWidth = GuiSettings.MassHistogramBinWidth;
                    bar.FillColor = Color.FromArgb(GuiSettings.FillAlpha, color);
                    bar.BorderLineWidth = 0;
                    bar.Orientation = ScottPlot.Orientation.Horizontal;
                    bar.ValueOffsets = locations.Select(p => (double)fileNum).ToArray();
                    bar.Label = PfmXplorerUtil.GetFileNameWithoutExtension(file.Key);

                    fileNum++;
                }

                plot.YAxis.Label("Daltons");
                plot.SetAxisLimitsY(0, maxMass * 1.2);

                double[] xPositions = Enumerable.Range(0, DataManagement.SpectraFiles.Count).Select(p => (double)p).ToArray();
                string[] xLabels = DataManagement.SpectraFiles.Select(p => PfmXplorerUtil.GetFileNameWithoutExtension(p.Key).ConvertName()).ToArray();
                plot.XTicks(xPositions, xLabels);

                StyleDashboardPlot(plot);
                plot.Legend(enable: false);
            }
            catch (Exception e)
            {
                errors.Add(e.Message);
            }
        }

        public static void OnSpeciesChanged()
        {

        }

        public static Text OnSpectrumPeakSelected(double xClickLocation, MsDataScan scan, Text textAnnotation)
        {
            if (textAnnotation == null)
            {
                textAnnotation = new Text();
            }

            var mzIndex = scan.MassSpectrum.GetClosestPeakIndex(xClickLocation);

            double mz = scan.MassSpectrum.XArray[mzIndex];
            double intensity = scan.MassSpectrum.YArray[mzIndex];

            textAnnotation.Label = mz.ToString("F3");
            textAnnotation.X = mz;
            textAnnotation.Y = intensity;
            textAnnotation.Alignment = Alignment.LowerCenter;
            textAnnotation.FontSize = (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingX);

            Tolerance t = new AbsoluteTolerance(0.0001);

            if (DataManagement.CurrentlySelectedFile.Value.OneBasedScanToAnnotatedEnvelopes.TryGetValue(scan.OneBasedScanNumber, out var envelopes))
            {
                var items = envelopes.Where(p => p.PeakMzs.Any(v => t.Within(mz, v))).ToList();

                if (items.Any())
                {
                    string label = mz.ToString("F3") + "\n" + "z=";

                    label += string.Join('|', items.Select(p => p.Charge));

                    textAnnotation.Label = label;
                }
            }

            return textAnnotation;
        }

        public static VLine UpdateRtIndicator(MsDataScan scan, VLine indicator, Plot plot)
        {
            if (scan == null)
            {
                return indicator;
            }

            if (indicator == null)
            {
                indicator = plot.AddVerticalLine(scan.RetentionTime, color: GuiSettings.RtIndicatorColor,
                    width: (float)(GuiSettings.ChartLineWidth * GuiSettings.DpiScalingX));
            }

            indicator.X = scan.RetentionTime;
            indicator.Color = GuiSettings.RtIndicatorColor;

            // this does not seem to update properly. need to call .Render() on the WpfPlot object, not the Plot
            //plot.Render();

            return indicator;
        }
    }
}
