﻿using MassSpectrometry;
using MzLibUtil;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ProteoformExplorer.Core
{
    public class CachedSpectraFileData
    {
        public KeyValuePair<string, DynamicDataConnection> DataFile { get; private set; }
        public Dictionary<int, HashSet<AnnotatedSpecies>> OneBasedScanToAnnotatedSpecies { get; private set; }
        public Dictionary<int, List<AnnotatedEnvelope>> OneBasedScanToAnnotatedEnvelopes { get; private set; }
        private List<Datum> TicData;
        private List<Datum> IdentifiedTicData;
        private List<Datum> DeconvolutedTicData;
        private static Dictionary<(string, int), MsDataScan> CachedScans;
        private static int NumScansToCache;
        private static Queue<(string, int)> CachedScanNumberQueue;

        public CachedSpectraFileData(KeyValuePair<string, DynamicDataConnection> loadedDataFile)
        {
            DataFile = loadedDataFile;
            TicData = new List<Datum>();
            IdentifiedTicData = new List<Datum>();
            DeconvolutedTicData = new List<Datum>();
            OneBasedScanToAnnotatedSpecies = new Dictionary<int, HashSet<AnnotatedSpecies>>();
            OneBasedScanToAnnotatedEnvelopes = new Dictionary<int, List<AnnotatedEnvelope>>();
            CachedScans = new Dictionary<(string, int), MsDataScan>();
            NumScansToCache = 1000;
            CachedScanNumberQueue = new Queue<(string, int)>();
        }

        public void BuildScanToSpeciesDictionary(List<AnnotatedSpecies> allAnnotatedSpecies)
        {
            OneBasedScanToAnnotatedSpecies.Clear();

            foreach (var species in allAnnotatedSpecies.Where(p => p.SpectraFileNameWithoutExtension == Path.GetFileNameWithoutExtension(DataFile.Key)))
            {
                if (species.DeconvolutionFeature == null ||
                    species.DeconvolutionFeature.AnnotatedEnvelopes == null ||
                    species.DeconvolutionFeature.AnnotatedEnvelopes.Count == 0)
                {
                    // the deconvoluted species is from a file type that does not specify the envelopes in the deconvolution feature
                    // therefore, we'll have to do some peakfinding and guess what envelopes are part of this feature
                    if (species.DeconvolutionFeature == null)
                    {
                        if (species.Identification != null)
                        { 
                            //TODO: figure this out
                            species.DeconvolutionFeature = new DeconvolutionFeature(species.Identification.MonoisotopicMass,
                                0, 0, 0, new List<int> { species.Identification.PrecursorChargeState }, DataFile.Key, null);
                        }
                        else
                        {
                            // TODO: some kind of error message? or just skip? this species doesn't have a deconvolution feature or an identification...
                        }
                    }

                    species.DeconvolutionFeature.FindAnnotatedEnvelopesInData(new KeyValuePair<string, CachedSpectraFileData>(DataFile.Key, this));
                }

                if (species.DeconvolutionFeature != null)
                {
                    foreach (var envelope in species.DeconvolutionFeature.AnnotatedEnvelopes)
                    {
                        int scanNum = envelope.OneBasedScanNumber;

                        if (!OneBasedScanToAnnotatedSpecies.ContainsKey(scanNum))
                        {
                            OneBasedScanToAnnotatedSpecies.Add(scanNum, new HashSet<AnnotatedSpecies>());
                        }
                        if (!OneBasedScanToAnnotatedEnvelopes.ContainsKey(scanNum))
                        {
                            OneBasedScanToAnnotatedEnvelopes.Add(scanNum, new List<AnnotatedEnvelope>());
                        }

                        OneBasedScanToAnnotatedSpecies[scanNum].Add(species);
                        OneBasedScanToAnnotatedEnvelopes[scanNum].Add(envelope);
                    }
                }

                if (species.Identification != null)
                {
                    int scanNum = species.Identification.OneBasedPrecursorScanNumber;

                    if (!OneBasedScanToAnnotatedSpecies.ContainsKey(scanNum))
                    {
                        OneBasedScanToAnnotatedSpecies.Add(scanNum, new HashSet<AnnotatedSpecies>());
                    }

                    OneBasedScanToAnnotatedSpecies[scanNum].Add(species);
                }
            }
        }

        public MsDataScan GetOneBasedScan(int oneBasedScanNum)
        {
            lock (CachedScans)
            {
                if (!CachedScans.TryGetValue((this.DataFile.Key, oneBasedScanNum), out var scan))
                {
                    scan = DataFile.Value.GetOneBasedScanFromDynamicConnection(oneBasedScanNum);

                    if (scan == null)
                    {
                        return scan;
                    }

                    while (CachedScans.Count > NumScansToCache)
                    {
                        var scanToRemove = CachedScanNumberQueue.Dequeue();
                        CachedScans.Remove(scanToRemove);
                    }

                    if (CachedScans.TryAdd((this.DataFile.Key, scan.OneBasedScanNumber), scan))
                    {
                        CachedScanNumberQueue.Enqueue((this.DataFile.Key, oneBasedScanNum));
                    }
                }

                return scan;
            }
        }

        public List<Datum> GetTicChromatogram()
        {
            HashSet<double> claimedMzs = new HashSet<double>();
            if (TicData.Count == 0)
            {
                int lastScanNum = PfmXplorerUtil.GetLastOneBasedScanNumber(new KeyValuePair<string, CachedSpectraFileData>(DataFile.Key, this));

                for (int i = 1; i <= lastScanNum; i++)
                {
                    claimedMzs.Clear();
                    var scan = GetOneBasedScan(i);
                    double deconvolutedTic = 0;

                    // tic
                    if (scan != null && scan.MsnOrder == 1)
                    {
                        TicData.Add(new Datum(scan.RetentionTime, scan.TotalIonCurrent, scan.OneBasedScanNumber));
                    }

                    // deconvoluted tic
                    if (OneBasedScanToAnnotatedEnvelopes.TryGetValue(i, out var annotatedEnvelopes) && scan.MsnOrder == 1)
                    {
                        foreach (var envelope in annotatedEnvelopes)
                        {
                            foreach (double mz in envelope.PeakMzs)
                            {
                                int index = scan.MassSpectrum.GetClosestPeakIndex(mz);
                                double actualMz = scan.MassSpectrum.XArray[index];

                                if (claimedMzs.Contains(actualMz))
                                {
                                    continue;
                                }

                                claimedMzs.Add(actualMz);
                                deconvolutedTic += scan.MassSpectrum.YArray[index];
                            }
                        }

                        DeconvolutedTicData.Add(new Datum(scan.RetentionTime, deconvolutedTic, scan.OneBasedScanNumber));
                    }

                    // TODO: identified tic
                }
            }

            return TicData;
        }

        public List<Datum> GetDeconvolutedTicChromatogram()
        {
            GetTicChromatogram();

            return DeconvolutedTicData;
        }

        public List<Datum> GetIdentifiedTicChromatogram()
        {
            GetTicChromatogram();

            return IdentifiedTicData;
        }

        public HashSet<AnnotatedSpecies> SpeciesInScan(int oneBasedScan)
        {
            if (OneBasedScanToAnnotatedSpecies.TryGetValue(oneBasedScan, out var species))
            {
                return species;
            }

            return null;
        }
    }
}