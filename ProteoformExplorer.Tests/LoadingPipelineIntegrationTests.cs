using MassSpectrometry;
using NUnit.Framework;
using ProteoformExplorer.Core;
using ProteoformExplorer.Core.IO;
using ProteoformExplorer.GuiFunctions;
using ProteoformExplorer.Wpf;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Reflection;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Task 6 loading-pipeline coverage: OSM identification rows and FlashLFQ quantified peak rows
    /// are appended as separate per-row records (no dedup/merge), normalized file names associate
    /// with loaded spectra keys, and FlashLFQ Peak intensity never becomes model or plot input.
    /// </summary>
    public class LoadingPipelineIntegrationTests
    {
        private static string TestDataRoot => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
        private static string OsmFilePath => Path.Combine(TestDataRoot, @"OSM\AllOSMs_mini.osmtsv");
        private static string QuantifiedPeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\TestQuantifiedPeaks.tsv");
        private static string EdgePeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\QuantifiedPeakMappingEdgePeaks.tsv");
        private static string IntensityVariantPeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\IntensityVariantPeaks.tsv");

        [OneTimeSetUp]
        public void Setup()
        {
            Loaders.LoadElements();
        }

        private static List<AnnotatedSpecies> ProcessResultFile(string path)
        {
            var processor = FileProcessor.GetProcessor(path);
            Assert.That(processor, Is.Not.Null, $"{Path.GetFileName(path)} must be routed to a processor");
            return processor.Process(path);
        }

        private static void WithFreshLoadedSpeciesCollection(Action<ObservableCollection<AnnotatedSpecies>> body)
        {
            var previous = DataManagement.AllLoadedAnnotatedSpecies;
            DataManagement.AllLoadedAnnotatedSpecies = new ObservableCollection<AnnotatedSpecies>();
            try
            {
                body(DataManagement.AllLoadedAnnotatedSpecies);
            }
            finally
            {
                DataManagement.AllLoadedAnnotatedSpecies = previous;
            }
        }

        [Test]
        public void TestOsmRowsRemainSeparateAndQuantifiedRowsCollapse()
        {
            var osmGridFile = new FileForDataGrid(OsmFilePath, out var osmErrors);
            var quantifiedGridFile = new FileForDataGrid(QuantifiedPeakPath, out var quantifiedErrors);

            Assert.That(osmErrors, Is.Empty);
            Assert.That(quantifiedErrors, Is.Empty);
            Assert.That(osmGridFile.FileType, Is.EqualTo(FileType.Identification));
            Assert.That(quantifiedGridFile.FileType, Is.EqualTo(FileType.Deconvolution));

            var osmSpecies = ProcessResultFile(OsmFilePath);
            var quantifiedSpecies = ProcessResultFile(QuantifiedPeakPath);

            Assert.That(osmSpecies, Has.Count.EqualTo(5), "fixture contains 5 OSM rows, one exact duplicate");
            Assert.That(quantifiedSpecies, Has.Count.EqualTo(2), "fixture contains 3 quantified rows, one exact duplicate");

            WithFreshLoadedSpeciesCollection(loaded =>
            {
                DataManagement.AddLoadedSpecies(osmSpecies);
                DataManagement.AddLoadedSpecies(quantifiedSpecies);

                Assert.That(loaded, Has.Count.EqualTo(7),
                    "every OSM row and every grouped quantified peak must yield one loaded record (5 OSM + 2 quantified)");

                for (int i = 0; i < osmSpecies.Count; i++)
                {
                    Assert.That(loaded[i], Is.SameAs(osmSpecies[i]), $"OSM row {i} must be appended as its own instance");
                }

                for (int i = 0; i < quantifiedSpecies.Count; i++)
                {
                    Assert.That(loaded[osmSpecies.Count + i], Is.SameAs(quantifiedSpecies[i]),
                        $"quantified row {i} must be appended as its own instance");
                }

                Assert.That(loaded.Take(5).All(p => p.Identification != null && p.DeconvolutionFeature == null), Is.True,
                    "OSM records stay identification-backed");
                Assert.That(loaded.Skip(5).All(p => p.Identification == null && p.DeconvolutionFeature != null), Is.True,
                    "quantified records stay deconvolution-only");

                var osmFirst = loaded[0];
                var osmDuplicate = loaded[3];
                Assert.That(osmDuplicate, Is.Not.SameAs(osmFirst));
                Assert.That(osmDuplicate.Identification, Is.Not.SameAs(osmFirst.Identification));
                Assert.That(osmDuplicate.Identification.IdentificationScanNum,
                    Is.EqualTo(osmFirst.Identification.IdentificationScanNum));

                var quantifiedFirst = loaded[5];
                var quantifiedSecond = loaded[6];
                Assert.That(quantifiedSecond, Is.Not.SameAs(quantifiedFirst));
                Assert.That(quantifiedSecond.DeconvolutionFeature, Is.Not.SameAs(quantifiedFirst.DeconvolutionFeature));
                Assert.That(quantifiedSecond.DeconvolutionFeature.MonoisotopicMass,
                    Is.Not.EqualTo(quantifiedFirst.DeconvolutionFeature.MonoisotopicMass));
            });
        }

        [Test]
        public void TestNormalizedFileNamesAssociateWithLoadedSpectraKeysAndUnmatchedRowsAreRetained()
        {
            var edgeSpecies = ProcessResultFile(EdgePeakPath);
            var quantifiedSpecies = ProcessResultFile(QuantifiedPeakPath);
            var osmSpecies = ProcessResultFile(OsmFilePath);

            Assert.That(edgeSpecies, Has.Count.EqualTo(9), "every mzLib row is mapped, with identical rows collapsed");

            string edgeRawKey = PfmXplorerUtil.GetFileNameWithoutExtension(@"loaded\EdgeRaw.raw");
            string testRawKey = PfmXplorerUtil.GetFileNameWithoutExtension(@"loaded\TestRawFile.raw");
            string osmRawKey = PfmXplorerUtil.GetFileNameWithoutExtension(@"loaded\20250612_RNA-Mix_10V.raw");

            Assert.That(edgeRawKey, Is.EqualTo("edgeraw"));
            Assert.That(edgeSpecies.Where(p => !string.IsNullOrEmpty(p.SpectraFileNameWithoutExtension))
                .All(p => p.SpectraFileNameWithoutExtension == edgeRawKey), Is.True);
            Assert.That(quantifiedSpecies.All(p => p.SpectraFileNameWithoutExtension == testRawKey), Is.True,
                "a File Name without the .raw suffix must match the key computed for the loaded raw file");
            Assert.That(osmSpecies.Select(p => p.SpectraFileNameWithoutExtension), Does.Contain(osmRawKey));

        }

        [Test]
        public void TestPeakIntensityDoesNotChangeMappedModelInputs()
        {
            var rows = FileReader.ReadFile<QuantifiedPeakFile>(IntensityVariantPeakPath).Results;
            Assert.That(rows, Has.Count.EqualTo(2));
            Assert.That(rows[0].PeakIntensity, Is.Not.EqualTo(rows[1].PeakIntensity),
                "fixture rows must differ in FlashLFQ Peak intensity");
            Assert.That(rows[1].PeptideMonoisotopicMass, Is.EqualTo(rows[0].PeptideMonoisotopicMass));
            Assert.That(rows[1].PeakRTApex, Is.EqualTo(rows[0].PeakRTApex));
            Assert.That(rows[1].PeakCharge, Is.EqualTo(rows[0].PeakCharge));

            var species = ProcessResultFile(IntensityVariantPeakPath);
            Assert.That(species, Has.Count.EqualTo(1));

            var first = species[0].DeconvolutionFeature;
            Assert.That(first.MonoisotopicMass, Is.EqualTo(rows[0].PeptideMonoisotopicMass));
            Assert.That(first.ApexRt, Is.EqualTo(rows[0].PeakRTApex));
            Assert.That(first.Charges, Is.EqualTo(new List<int> { rows[0].PeakCharge.Value }));
            Assert.That(species[0].SpeciesLabel, Is.EqualTo(rows[0].FullSequence));
            Assert.That(species[0].Identification, Is.Null);
            Assert.That(first.AnnotatedEnvelopes, Is.Null);
        }

        [Test]
        public void TestModelAndPlotInputsExposeNoFlashLfqIntensityChannel()
        {
            foreach (Type modelType in new[] { typeof(AnnotatedSpecies), typeof(ProteoformExplorer.Core.DeconvolutionFeature) })
            {
                var suspects = modelType.GetProperties(BindingFlags.Public | BindingFlags.Instance)
                    .Where(p => p.Name.Contains("intensity", StringComparison.OrdinalIgnoreCase)
                             || p.PropertyType.Name.Contains("QuantifiedPeak", StringComparison.OrdinalIgnoreCase))
                    .Select(p => $"{modelType.Name}.{p.Name}")
                    .ToList();

                Assert.That(suspects, Is.Empty, "mapped model types must not expose a FlashLFQ intensity channel");
            }

            var quantifiedRowTypes = new[] { typeof(QuantifiedPeak), typeof(QuantifiedPeakFile) };

            foreach (MethodInfo method in typeof(PlottingFunctions).GetMethods(BindingFlags.Public | BindingFlags.Static))
            {
                foreach (ParameterInfo parameter in method.GetParameters())
                {
                    Assert.That(quantifiedRowTypes, Does.Not.Contain(parameter.ParameterType),
                        $"plot entry point {method.Name} must not consume FlashLFQ peak data");
                }
            }

            foreach (string name in new[] { nameof(PfmXplorerUtil.GetXicData), nameof(PfmXplorerUtil.GetSummedChargeXics) })
            {
                var xicMethods = typeof(PfmXplorerUtil).GetMethods(BindingFlags.Public | BindingFlags.Static)
                    .Where(m => m.Name == name)
                    .ToList();
                Assert.That(xicMethods, Has.Count.EqualTo(1), $"{name} must exist exactly once");

                var parameters = xicMethods[0].GetParameters();
                Assert.That(parameters.Select(p => p.ParameterType), Does.Contain(typeof(List<MsDataScan>)),
                    $"{name} must produce XIC values from raw scans");

                foreach (ParameterInfo parameter in parameters)
                {
                    Assert.That(quantifiedRowTypes, Does.Not.Contain(parameter.ParameterType),
                        $"{name} must not consume FlashLFQ peak data");
                }
            }
        }
    }
}
