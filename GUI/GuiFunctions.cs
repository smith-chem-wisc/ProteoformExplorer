﻿using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using mzPlot;
using ProteoformExplorer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GUI
{
    public class GuiFunctions
    {
        private static Dictionary<(string, int), MsDataScan> CachedScans = new Dictionary<(string, int), MsDataScan>();

        public static Plot PlotSpeciesInSpectrum(List<AnnotatedSpecies> allSpeciesToPlot, int oneBasedScan, KeyValuePair<string, DynamicDataConnection> data,
            Plot plot, OxyPlot.Wpf.PlotView plotView, bool clearOldPlot)
        {
            // get the scan
            var scan = data.Value.GetOneBasedScanFromDynamicConnection(oneBasedScan);

            // add non-annotated peaks
            List<Datum> spectrumData = new List<Datum>();

            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
            {
                spectrumData.Add(new Datum(scan.MassSpectrum.XArray[i], scan.MassSpectrum.YArray[i]));
            }

            plot = new SpectrumPlot(plotView, spectrumData, refreshAfterAddingData: false);

            // add annotated peaks
            HashSet<double> claimedMzs = new HashSet<double>();
            foreach (var species in allSpeciesToPlot)
            {
                List<Datum> annotatedData = new List<Datum>();
                List<int> chargesToPlot = new List<int>();

                if (species.DeconvolutionFeature != null)
                {
                    double mass = HomePage.DeconvolutionEngine.GetModeMassFromMonoisotopicMass(species.DeconvolutionFeature.MonoisotopicMass);
                    chargesToPlot.AddRange(species.DeconvolutionFeature.Charges);

                    foreach (var z in chargesToPlot)
                    {
                        int index = scan.MassSpectrum.GetClosestPeakIndex(mass.ToMz(z));
                        double expMz = scan.MassSpectrum.XArray[index];
                        double expIntensity = scan.MassSpectrum.YArray[index];

                        var envelope = HomePage.DeconvolutionEngine.GetIsotopicEnvelope(scan.MassSpectrum, index, z, new List<Deconvoluter.DeconvolutedPeak>(),
                            claimedMzs, new List<(double, double)>());

                        if (envelope != null)
                        {
                            annotatedData.AddRange(envelope.Peaks.Select(p => new Datum(p.ExperimentalMz, p.ExperimentalIntensity)));
                        }
                        else
                        {
                            annotatedData.Add(new Datum(expMz, expIntensity));
                        }
                    }

                    foreach (var item in annotatedData)
                    {
                        claimedMzs.Add(item.X);
                    }
                }
                else if (species.Identification != null)
                {
                    double mass = HomePage.DeconvolutionEngine.GetModeMassFromMonoisotopicMass(species.Identification.MonoisotopicMass);
                    int z = species.Identification.PrecursorChargeState;

                    int index = scan.MassSpectrum.GetClosestPeakIndex(mass.ToMz(z));
                    double expMz = scan.MassSpectrum.XArray[index];
                    double expIntensity = scan.MassSpectrum.YArray[index];

                    var envelope = HomePage.DeconvolutionEngine.GetIsotopicEnvelope(scan.MassSpectrum, index, z, new List<Deconvoluter.DeconvolutedPeak>(),
                        claimedMzs, new List<(double, double)>());
                }

                plot.AddSpectrumPlot(annotatedData, OxyPlot.OxyColors.Blue, 2.0);
                ZoomAxes(annotatedData, plot);
            }

            return plot;
        }

        public static Plot PlotSpeciesInXic(double mz, int z, Tolerance tolerance, double rt, double rtWindow, KeyValuePair<string, DynamicDataConnection> data,
            Plot plot, OxyPlot.Wpf.PlotView plotView, bool clearOldPlot)
        {
            double rtWindowHalfWidth = rtWindow / 2;
            var startScan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(data, rt - rtWindowHalfWidth);
            var endScan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(data, rt + rtWindowHalfWidth);

            List<MsDataScan> scans = new List<MsDataScan>();
            for (int i = startScan.OneBasedScanNumber; i <= endScan.OneBasedScanNumber; i++)
            {
                var theScan = data.Value.GetOneBasedScanFromDynamicConnection(i);

                if (theScan.RetentionTime > startScan.RetentionTime + rtWindowHalfWidth)
                {
                    break;
                }

                if (theScan.MsnOrder == 1)
                {
                    scans.Add(theScan);
                }
            }

            List<Datum> xicData = new List<Datum>();
            foreach (var scan in scans)
            {
                int ind = scan.MassSpectrum.GetClosestPeakIndex(mz);
                double expMz = scan.MassSpectrum.XArray[ind];
                double expIntensity = scan.MassSpectrum.YArray[ind];

                if (tolerance.Within(expMz.ToMass(z), mz.ToMass(z)))
                {
                    xicData.Add(new Datum(scan.RetentionTime, expIntensity));
                }
                else
                {
                    xicData.Add(new Datum(scan.RetentionTime, 0));
                }
            }

            if (clearOldPlot)
            {
                plot = new LinePlot(plotView, xicData);
            }
            else
            {
                plot.AddLinePlot(xicData);
            }

            return plot;
        }

        public static void ZoomAxes(List<Datum> annotatedIons, Plot plot, double yZoom = 1.2)
        {
            double highestAnnotatedIntensity = 0;
            double highestAnnotatedMz = double.MinValue;
            double lowestAnnotatedMz = double.MaxValue;

            foreach (var ion in annotatedIons)
            {
                double mz = ion.X;
                double intensity = ion.Y.Value;

                highestAnnotatedIntensity = Math.Max(highestAnnotatedIntensity, intensity);
                highestAnnotatedMz = Math.Max(highestAnnotatedMz, mz);
                lowestAnnotatedMz = Math.Min(lowestAnnotatedMz, mz);
            }

            if (highestAnnotatedIntensity > 0)
            {
                plot.Model.Axes[1].Zoom(0, highestAnnotatedIntensity * yZoom);
            }

            if (highestAnnotatedMz > double.MinValue && lowestAnnotatedMz < double.MaxValue)
            {
                plot.Model.Axes[0].Zoom(lowestAnnotatedMz - 100, highestAnnotatedMz + 100);
            }
        }

        //public static MsDataScan TryGetCachedScan(string fullFilePath, int oneBasedScanNumber, bool autoCache = true)
        //{
        //    if (CachedScans.TryGetValue((fullFilePath, oneBasedScanNumber), out var scan))
        //    {
        //        return scan;
        //    }


        //}
    }
}