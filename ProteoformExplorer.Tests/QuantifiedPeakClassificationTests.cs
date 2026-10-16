using NUnit.Framework;
using ProteoformExplorer.Core.IO;
using ProteoformExplorer.Wpf;
using Readers;
using System.IO;

namespace Test
{
    public class QuantifiedPeakClassificationTests
    {
        private static string TestDataRoot => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
        private static string ValidPeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\TestQuantifiedPeaks.tsv");
        private static string MalformedPeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\MalformedPeaks.tsv");
        private static string OsmPath => Path.Combine(TestDataRoot, @"OSM\AllOSMs_mini.osmtsv");
        private static string PsmPath => Path.Combine(TestDataRoot, @"mcf7_sliced\AllPSMs.psmtsv");
        private static string FlashDeconPath => Path.Combine(TestDataRoot, "FlashDeconTest.tsv");
        private static string MetaMorpheusPath => Path.Combine(TestDataRoot, "MetaMorpheus.tsv");
        private static string TdPortalPath => Path.Combine(TestDataRoot, "td_portal.tsv");

        [Test]
        public void TestValidQuantifiedPeakFileIsClassifiedAsDeconvolution()
        {
            var file = new FileForDataGrid(ValidPeakPath, out var errors);

            Assert.That(errors, Is.Empty);
            Assert.That(file.FileType, Is.EqualTo(FileType.Deconvolution));
        }

        [Test]
        public void TestValidQuantifiedPeakFileIsRoutedToQuantifiedPeakProcessor()
        {
            Assert.That(ValidPeakPath.ParseFileType(), Is.EqualTo(SupportedFileType.FlashLFQQuantifiedPeak));
            Assert.That(FileProcessor.GetProcessor(ValidPeakPath), Is.Not.Null);
        }

        [Test]
        public void TestPeaksSuffixIsRecognizedAndRoutedByMzLibType()
        {
            Assert.That(MalformedPeakPath.ParseFileType(), Is.EqualTo(SupportedFileType.FlashLFQQuantifiedPeak));
            var file = new FileForDataGrid(MalformedPeakPath, out _);
            Assert.That(file.FileType, Is.EqualTo(FileType.Deconvolution));
            Assert.That(FileProcessor.GetProcessor(MalformedPeakPath), Is.Not.Null);
        }

        [Test]
        public void TestExistingLegacyFileClassificationsAreUnchanged()
        {
            var flashDecon = new FileForDataGrid(FlashDeconPath, out _);
            Assert.That(flashDecon.FileType, Is.EqualTo(FileType.Deconvolution));

            var metaMorpheus = new FileForDataGrid(MetaMorpheusPath, out _);
            Assert.That(metaMorpheus.FileType, Is.EqualTo(FileType.Identification));

            var tdPortal = new FileForDataGrid(TdPortalPath, out _);
            Assert.That(tdPortal.FileType, Is.EqualTo(FileType.Identification));
        }

        [Test]
        public void TestOsmAndPsmRoutesAreUnchanged()
        {
            Assert.That(FileProcessor.GetProcessor(OsmPath), Is.Not.Null);
            Assert.That(new FileForDataGrid(OsmPath, out _).FileType, Is.EqualTo(FileType.Identification));
            Assert.That(FileProcessor.GetProcessor(PsmPath), Is.Not.Null);
        }

        [Test]
        public void TestQuantifiedPeakFixtureIsReadableByMzLibReader()
        {
            var quantifiedPeaks = FileReader.ReadFile<QuantifiedPeakFile>(ValidPeakPath).Results;

            Assert.That(quantifiedPeaks.Count, Is.EqualTo(3));
            Assert.That(quantifiedPeaks[0].FileName, Is.EqualTo("TestRawFile"));
            Assert.That(quantifiedPeaks[0].PeptideMonoisotopicMass, Is.EqualTo(3200.12345).Within(1e-5));
            Assert.That(quantifiedPeaks[0].PeakRTApex, Is.EqualTo(25.4).Within(1e-5));
            Assert.That(quantifiedPeaks[0].PeakCharge, Is.EqualTo(-3));
            Assert.That(quantifiedPeaks[1].PeakRTApex, Is.EqualTo(30.0).Within(1e-5));
            Assert.That(quantifiedPeaks[2].PeakCharge, Is.EqualTo(quantifiedPeaks[0].PeakCharge));
        }
    }
}
