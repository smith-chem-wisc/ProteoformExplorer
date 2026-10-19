using NUnit.Framework;
using ProteoformExplorer.Core;
using ProteoformExplorer.Core.IO;
using ProteoformExplorer.Wpf;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    public class OsmLoadingTests
    {
        private static string OsmFilePath => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\OSM\AllOSMs_mini.osmtsv");

        [OneTimeSetUp]
        public void Setup()
        {
            Loaders.LoadElements();
        }

        [Test]
        public void TestOsmFileIsRoutedToAProcessor()
        {
            var processor = FileProcessor.GetProcessor(OsmFilePath);
            Assert.That(processor, Is.Not.Null);
        }

        [Test]
        public void TestOsmRowsLoadAsIndependentIdentificationBackedSpecies()
        {
            var species = FileProcessor.GetProcessor(OsmFilePath).Process(OsmFilePath);

            Assert.That(species.Count, Is.EqualTo(5));
            Assert.That(species.All(p => p.Identification != null), Is.True);
            Assert.That(species.All(p => p.DeconvolutionFeature == null), Is.True);
            Assert.That(species.All(p => p.SpeciesLabel == p.Identification.FullSequence), Is.True);
        }

        [Test]
        public void TestOsmRowMappingPreservesObservedMassAndNegativeCharge()
        {
            var species = FileProcessor.GetProcessor(OsmFilePath).Process(OsmFilePath);
            var id = species[0].Identification;

            Assert.That(id.BaseSequence, Is.EqualTo("ACGUACGUACGU"));
            Assert.That(id.FullSequence, Is.EqualTo("ACGUACGUACGU"));
            // the observed "Precursor Mass" value must be mapped, not the theoretical "Monoisotopic Mass" column
            Assert.That(id.MonoisotopicMass, Is.EqualTo(7000.12345).Within(1e-5));
            Assert.That(id.MonoisotopicMass, Is.Not.EqualTo(7000.15));
            Assert.That(id.PrecursorChargeState, Is.EqualTo(-5));
            Assert.That(id.OneBasedPrecursorScanNumber, Is.EqualTo(2580));
            Assert.That(id.IdentificationScanNum, Is.EqualTo(2590));
            Assert.That(id.SpectraFileNameWithoutExtension, Is.EqualTo("20250612_rna-mix_10v"));
            Assert.That(id.Dataset, Is.EqualTo("TestData"));
        }

        [Test]
        public void TestOsmDuplicateRowsRemainSeparateSpecies()
        {
            var species = FileProcessor.GetProcessor(OsmFilePath).Process(OsmFilePath);

            var first = species[0];
            var duplicate = species[3];

            Assert.That(duplicate.Identification.MonoisotopicMass, Is.EqualTo(first.Identification.MonoisotopicMass));
            Assert.That(duplicate.Identification.IdentificationScanNum, Is.EqualTo(first.Identification.IdentificationScanNum));
            Assert.That(duplicate.SpeciesLabel, Is.EqualTo(first.SpeciesLabel));
            Assert.That(duplicate, Is.Not.SameAs(first));
            Assert.That(duplicate.Identification, Is.Not.SameAs(first.Identification));
        }

        [Test]
        public void TestOsmChargesAndModifiedSequenceArePreserved()
        {
            var species = FileProcessor.GetProcessor(OsmFilePath).Process(OsmFilePath);

            Assert.That(species.Select(p => p.Identification.PrecursorChargeState), Is.EqualTo(new[] { -5, -5, -4, -5, -6 }));
            Assert.That(species[2].Identification.FullSequence, Is.EqualTo("ACGUACGU[Biological:Methylation on A]ACGU"));
            Assert.That(species[2].Identification.BaseSequence, Is.EqualTo("ACGUACGUACGU"));
        }

        [Test]
        public void TestOsmFileIsClassifiedAsIdentificationInput()
        {
            string header = File.ReadLines(OsmFilePath).First();

            Assert.That(InputReaderParser.AcceptedTextFileFormats.Contains(".osmtsv"), Is.True);
            Assert.That(InputReaderParser.GetFileTypeFromHeader(header), Is.EqualTo(InputReaderParser.InputSourceType.MetaMorpheus));

            var file = new FileForDataGrid(OsmFilePath, out var errors);
            Assert.That(errors, Is.Empty);
            Assert.That(file.FileType, Is.EqualTo(FileType.Identification));
        }

        [Test]
        public void TestPsmtsvProcessorRouteStillLoads()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mcf7_sliced\AllPSMs.psmtsv");
            var species = FileProcessor.GetProcessor(psmFilePath).Process(psmFilePath);

            Assert.That(species.Count, Is.EqualTo(12));
            Assert.That(species.All(p => p.Identification != null), Is.True);
        }
    }
}
