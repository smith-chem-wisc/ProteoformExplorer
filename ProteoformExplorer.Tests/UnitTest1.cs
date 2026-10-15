using Chemistry;
using ProteoformExplorer.Deconvoluter;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using ProteoformExplorer.Core;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using Readers;
using ThermoFisher.CommonCore.Data.Business;

namespace Test
{
    public class Tests
    {
        [OneTimeSetUp]
        public void Setup()
        {
            Loaders.LoadElements();
        }

        [Test]
        public void TestReadingMetaMorpheus()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\MetaMorpheus.tsv");
            var species = InputReaderParser.ReadSpeciesFromFile(filePath, out var errors);

            Assert.That(!errors.Any());
            Assert.That(species.Count == 72);
            Assert.That(species.All(p => p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep2"
                || p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep1"));
            Assert.That(species.All(p => p.Identification.PrecursorChargeState > 0));
            Assert.That(species.All(p => p.Identification.OneBasedPrecursorScanNumber > 0));
        }

        [Test]
        public void TestReadingTdPortal()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\td_portal.tsv");
            var species = InputReaderParser.ReadSpeciesFromFile(filePath, out var errors);

            Assert.That(!errors.Any());
            Assert.That(species.Count == 12);
            Assert.That(species.All(p => p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep2"
                || p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep1"));

            species = species.Where(p => p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep1" 
                && p.Identification.IdentificationScanNum == 9).Take(1).ToList();
            string rawFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\05-26-17_B7A_yeast_td_fract7_rep1.raw");
            var connection = new KeyValuePair<string, MsDataFile>(rawFile, MsDataFileReader.GetDataFile(rawFile));
            var cachedData = new KeyValuePair<string, CachedSpectraFileData>(rawFile, new CachedSpectraFileData(connection));
            DataManagement.SpectraFiles = new Dictionary<string, CachedSpectraFileData> { { cachedData.Key, cachedData.Value } };
            cachedData.Value.CreateAnnotatedDeconvolutionFeatures(species);

            Assert.That(species.All(p => p.Identification.PrecursorChargeState > 0));
            Assert.That(species.All(p => p.Identification.OneBasedPrecursorScanNumber > 0));
        }

        [Test]
        public void TestReadingFlashDecon()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\FlashDeconTest.tsv");
            var species = InputReaderParser.ReadSpeciesFromFile(filePath, out var errors);

            Assert.That(!errors.Any());
            Assert.That(species.Count == 10);
            Assert.That(species.All(p => p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep2"
                || p.SpectraFileNameWithoutExtension == "05-26-17_B7A_yeast_td_fract7_rep1"));

            Assert.That(species[0].SpeciesLabel == "1 (6983.009)");
            Assert.That(species[0].DeconvolutionFeature.Charges.SequenceEqual(new List<int> { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 }));

            Assert.That(species[3].SpeciesLabel == "4 (13553.956)");
            Assert.That(species[3].DeconvolutionFeature.Charges.SequenceEqual(new List<int> { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 }));

            species = species.Where(p => (int)p.DeconvolutionFeature.MonoisotopicMass == 15229).ToList();
            string rawFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\05-26-17_B7A_yeast_td_fract7_rep1.raw");
            var connection = new KeyValuePair<string, MsDataFile>(rawFile, MsDataFileReader.GetDataFile(rawFile));
            var cachedData = new KeyValuePair<string, CachedSpectraFileData>(rawFile, new CachedSpectraFileData(connection));
            cachedData.Value.CreateAnnotatedDeconvolutionFeatures(species);

            Assert.That(species.First().DeconvolutionFeature.AnnotatedEnvelopes.Count >= 1);
        }

        [Test]
        public void TestReadingThermoDecon()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\3-19-18_MCF7_IM_br1_f1_2_results+_calibrated.txt");
            var species = InputReaderParser.ReadSpeciesFromFile(filePath, out var errors);

            Assert.That(!errors.Any());
            Assert.That(species.Count == 20);
            Assert.That(species.First().DeconvolutionFeature.SpectraFileNameWithoutExtension == "3-19-18_MCF7_IM_br1_f1_2_results+_calibrated");
            Assert.That(species.First().SpectraFileNameWithoutExtension == "3-19-18_MCF7_IM_br1_f1_2_results+_calibrated");

            Assert.That(species[0].SpeciesLabel == "1 (10835.862)");
            Assert.That(species[0].DeconvolutionFeature.Charges.SequenceEqual(new List<int> { 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 }));

            Assert.That(species[3].SpeciesLabel == "4 (13992.938)");
            Assert.That(species[3].DeconvolutionFeature.Charges.SequenceEqual(new List<int> { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 }));
        }

        [Test]
        public void TestDynamicConnectionRtBinarySearch()
        {
            string rawFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mcf7_sliced\mcf7_sliced_td.raw");
            var connection = new KeyValuePair<string, MsDataFile>(rawFile, MsDataFileReader.GetDataFile(rawFile));
            var cachedData = new KeyValuePair<string, CachedSpectraFileData>(rawFile, new CachedSpectraFileData(connection));

            var lastScanNum = PfmXplorerUtil.GetLastOneBasedScanNumber(cachedData);
            Assert.That(lastScanNum == 52);

            var scan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(cachedData, 57.0);
            Assert.That(scan.OneBasedScanNumber == 1);

            scan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(cachedData, 58.21);
            Assert.That(scan.OneBasedScanNumber == 30);

            scan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(cachedData, 59.0);
            Assert.That(scan.OneBasedScanNumber == 52);

            connection.Value.CloseDynamicConnection();
        }

        [Test]
        public void TestAdjacentMs1ScanTraversal()
        {
            string rawFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mcf7_sliced\mcf7_sliced_td.raw");
            var connection = new KeyValuePair<string, MsDataFile>(rawFile, MsDataFileReader.GetDataFile(rawFile));
            var cachedData = new KeyValuePair<string, CachedSpectraFileData>(rawFile, new CachedSpectraFileData(connection));
            var scans = connection.Value.Scans.ToList();

            foreach (var currentScan in scans)
            {
                var expectedPrevious = scans
                    .Where(scan => scan.OneBasedScanNumber < currentScan.OneBasedScanNumber && scan.MsnOrder == 1)
                    .LastOrDefault();
                var expectedNext = scans
                    .Where(scan => scan.OneBasedScanNumber > currentScan.OneBasedScanNumber && scan.MsnOrder == 1)
                    .FirstOrDefault();

                var previous = PfmXplorerUtil.GetAdjacentScanByMsOrder(cachedData, currentScan.OneBasedScanNumber, -1, 1);
                var next = PfmXplorerUtil.GetAdjacentScanByMsOrder(cachedData, currentScan.OneBasedScanNumber, 1, 1);

                Assert.That(previous?.OneBasedScanNumber, Is.EqualTo(expectedPrevious?.OneBasedScanNumber));
                Assert.That(next?.OneBasedScanNumber, Is.EqualTo(expectedNext?.OneBasedScanNumber));
            }

            connection.Value.CloseDynamicConnection();
        }

        [Test]
        public static void TestDeconvolution()
        {
            // >sp|P49703|ARL4D_HUMAN ADP-ribosylation factor-like protein 4D OS=Homo sapiens OX=9606 GN=ARL4D PE=1 SV=2
            string sequence = "MGNHLTEMAPTASSFLPHFQALHVVVIGLDSAGKTSLLYRLKFKEFVQSVPTKGFNTEKIRVPLGGSRGITFQ" +
                              "VWDVGGQEKLRPLWRSYTRRTDGLVFVVDAAEAERLEEAKVELHRISRASDNQGVPVLVLANKQDQPGALSAA" +
                              "EVEKRLAVRELAAATLTHVQGCSAVDGLGLQQGLERLYEMILKRKKAARGGKKRR";
            int charge = 15;
            double intensityMultiplier = 1e6;

            Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(sequence);
            var formula = baseSequence.GetChemicalFormula();
            var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);

            double[] masses = isotopicDistribution.Masses.ToArray();
            double[] abundances = isotopicDistribution.Intensities.ToArray();
            double max = abundances.Max();

            List<(double, double)> peaks = new List<(double, double)>();

            for (int i = 0; i < masses.Length; i++)
            {
                abundances[i] /= max;

                if (abundances[i] >= 0.05)
                {
                    peaks.Add((masses[i].ToMz(charge), abundances[i] * intensityMultiplier));
                }
            }

            Random r = new Random(1);
            for (int i = 0; i < 1000; i++)
            {
                double mz = r.NextDouble() * 1200 + 400;
                double intensity = r.NextDouble() * 30000 + 30000;
                peaks.Add((mz, intensity));
            }

            peaks = peaks.OrderBy(p => p.Item1).ToList();

            var spectrum = new MzSpectrum(peaks.Select(p => p.Item1).ToArray(), peaks.Select(p => p.Item2).ToArray(), true);

            var engine = new DeconvolutionEngine(0, 0.4, 3, 0.4, 1.5, 5, 1, 60, 2);
            var envs = engine.Deconvolute(spectrum, spectrum.Range, MassSpectrometry.Polarity.Positive).ToList();

            List<string> output = new List<string> { DeconvolutedEnvelope.TabDelimitedHeader };

            foreach (var env in envs)
            {
                env.SpectraFileName = "temp";
                output.Add(env.ToOutputString());
            }

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPfmExplorerDeconOutput.tsv");
            File.WriteAllLines(path, output);

            var species = InputReaderParser.ReadSpeciesFromFile(path, out var errors);

            Assert.That(!errors.Any());
            Assert.That(species.Count > 0);
            Assert.That(species.All(p => p.DeconvolutionFeature != null));
            Assert.That(species.All(p => p.DeconvolutionFeature.AnnotatedEnvelopes != null && p.DeconvolutionFeature.AnnotatedEnvelopes.Count > 0));
            Assert.That(species.All(p => p.DeconvolutionFeature.AnnotatedEnvelopes.All(v => v.PeakMzs.Count > 0)));

            Assert.That(species.Count == 1);
            Assert.That(species.First().DeconvolutionFeature.MonoisotopicMass > formula.MonoisotopicMass - 10
                && species.First().DeconvolutionFeature.MonoisotopicMass < formula.MonoisotopicMass + 10);
        }

        [Test]
        public static void TestDeconvolutionRealFile()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mcf7_sliced\mcf7_sliced_td.raw");
            var data = IO.ThermoRawFileReader.ThermoRawFileReader.LoadAllStaticData(filePath);
            var deconEngine = new DeconvolutionEngine(2000, 0.4, 3, 0.4, 3, 5, 2, 60, 3);
            var envelopes = deconEngine.Deconvolute(data, filePath, 1).ToList();

            List<string> output = new List<string>() { DeconvolutedEnvelope.TabDelimitedHeader };

            foreach (var env in envelopes)
            {
                output.Add(env.ToOutputString());
            }

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPfmExplorerDeconOutput_mcf7.tsv");
            File.WriteAllLines(path, output);

            var species = InputReaderParser.ReadSpeciesFromFile(path, out var errors);

            // populate the decon feature with data from the spectra file
            var file = new CachedSpectraFileData(new KeyValuePair<string, MsDataFile>(filePath, MsDataFileReader.GetDataFile(filePath)));
            file.CreateAnnotatedDeconvolutionFeatures(species);

            Assert.That(!errors.Any());
            Assert.That(species.Count > 0);
            Assert.That(species.All(p => p.DeconvolutionFeature != null));
            var debug = species.Where(p => p.DeconvolutionFeature.AnnotatedEnvelopes == null || p.DeconvolutionFeature.AnnotatedEnvelopes.Count == 0).ToList();

            Assert.That(species.All(p => p.DeconvolutionFeature.AnnotatedEnvelopes != null && p.DeconvolutionFeature.AnnotatedEnvelopes.Count > 0));
            Assert.That(species.All(p => p.DeconvolutionFeature.AnnotatedEnvelopes.All(v => v.PeakMzs.Count > 0)));
            Assert.That(species.All(p => p.DeconvolutionFeature.SpectraFileNameWithoutExtension == PfmXplorerUtil.GetFileNameWithoutExtension(filePath)));
        }

        //[Test]
        //[TestCase(@"032421_MALAT1Capture_11903_30632", 1136.1993, 30632.139, 27)] // charge too high
        //[TestCase(@"032421_MALAT1Capture_3033_5648", 628.9741, 5648.693, 9)] // very noisy/difficult

        //[TestCase(@"032421_MALAT1Capture_3090_13765", 689.6843, 13765.519, 20)] // easy case
        //[TestCase(@"032421_MALAT1Capture_5081_7011", 638.7745, 7011.428, 11)] // charge too low
        //[TestCase(@"032421_MALAT1Capture_5790_7045", 784.3348, 7045.937, 9)] // charge too low

        //public static void TestDeconvolutionDifferentCases(string dataFilePath, double mz, double monoMass, int z)
        //{
        //    string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DeconCases\" + dataFilePath + ".mzML");
        //    var data = IO.MzML.Mzml.LoadAllStaticData(filePath);
        //    var deconEngine = new DeconvolutionEngine(2000, 0.3, 6, 0.4, 3, 5, 2, 60, 2);

        //    var scan = data.GetAllScansList().First();

        //    Tolerance t = new AbsoluteTolerance(0.001);

        //    var candidates = deconEngine.GetEnvelopeCandidates(scan.MassSpectrum, scan.ScanWindowRange);
        //    var candidatesWithMz = candidates.Where(p => p.Peaks.Any(v => t.Within(v.ExperimentalMz, mz))).OrderByDescending(p => p.Score).ToList();

        //    var parsimonyEnvelopes = deconEngine.RunEnvelopeParsimony(candidates, scan.MassSpectrum);
        //    var parsimonyEnvsWithMz = parsimonyEnvelopes.Where(p => p.Peaks.Any(v => t.Within(v.ExperimentalMz, mz))).ToList();
        //}
    }
}
