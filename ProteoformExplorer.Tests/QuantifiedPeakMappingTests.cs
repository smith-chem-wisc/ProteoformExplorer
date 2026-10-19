using MassSpectrometry;
using NUnit.Framework;
using ProteoformExplorer.Core;
using ProteoformExplorer.Core.IO;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    public class QuantifiedPeakMappingTests
    {
        private static string TestDataRoot => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
        private static string ValidPeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\TestQuantifiedPeaks.tsv");
        private static string EdgePeakPath => Path.Combine(TestDataRoot, @"QuantifiedPeaks\QuantifiedPeakMappingEdgePeaks.tsv");

        private static List<AnnotatedSpecies> ProcessValid()
        {
            var processor = FileProcessor.GetProcessor(ValidPeakPath);
            Assert.That(processor, Is.Not.Null, "valid quantified peak file must be routed to a processor");
            return processor.Process(ValidPeakPath);
        }

        private static List<AnnotatedSpecies> ProcessEdge()
        {
            var processor = FileProcessor.GetProcessor(EdgePeakPath);
            Assert.That(processor, Is.Not.Null, "edge quantified peak file must be routed to a processor");
            return processor.Process(EdgePeakPath);
        }

        [Test]
        public void TestEachValidRowBecomesOneDeconvolutionOnlySpecies()
        {
            var species = ProcessValid();
            Assert.That(species, Has.Count.EqualTo(2));
            foreach (var item in species)
            {
                Assert.That(item.Identification, Is.Null, "quantified rows must not create an Identification");
                Assert.That(item.DeconvolutionFeature, Is.Not.Null);
            }
        }

        [Test]
        public void TestBasicFieldMappingUsesRowValues()
        {
            var species = ProcessValid();
            var feature = species[0].DeconvolutionFeature;

            Assert.That(species[0].SpeciesLabel, Is.EqualTo("ACGUACGU"));
            Assert.That(species[0].SpectraFileNameWithoutExtension, Is.EqualTo("testrawfile"));
            Assert.That(feature.SpectraFileNameWithoutExtension, Is.EqualTo("testrawfile"));
            Assert.That(feature.MonoisotopicMass, Is.EqualTo(3200.12345).Within(1e-5));
            Assert.That(feature.Charges, Is.EqualTo(new List<int> { -3 }));
            // RT values are minutes in the source row and must not be converted.
            Assert.That(feature.ApexRt, Is.EqualTo(25.4).Within(1e-5));
            Assert.That(feature.RtElutionRange.Minimum, Is.EqualTo(24.9).Within(1e-5));
            Assert.That(feature.RtElutionRange.Maximum, Is.EqualTo(25.9).Within(1e-5));
        }

        [Test]
        public void TestNegativePeakChargeIsPreservedAsSingletonChargeList()
        {
            var species = ProcessValid();
            var feature = species[0].DeconvolutionFeature;

            Assert.That(feature.Charges, Is.EqualTo(new List<int> { -3 }));
            Assert.That(feature.Polarity, Is.EqualTo(Polarity.Negative));
        }

        [Test]
        public void TestIdenticalRowsCollapseIntoOneFeature()
        {
            var species = ProcessValid();

            Assert.That(species, Has.Count.EqualTo(2), "identical quantified rows should collapse");
            Assert.That(species.Select(p => p.SpeciesLabel), Is.EquivalentTo(new[] { "ACGUACGU", "AUGC" }));
        }

        [Test]
        public void TestPeakChargeFallsBackToPrecursorChargeWhenCellEmpty()
        {
            var species = ProcessEdge();
            var feature = species[1].DeconvolutionFeature;

            Assert.That(feature.Charges, Is.EqualTo(new List<int> { 3 }));
            Assert.That(feature.Polarity, Is.EqualTo(Polarity.Positive));
        }

        [Test]
        public void TestRtValuesAreMappedDirectly()
        {
            var species = ProcessEdge();

            var emptyBoundaries = species[1].DeconvolutionFeature;
            Assert.That(emptyBoundaries.ApexRt, Is.EqualTo(14.0).Within(1e-5));
            Assert.That(emptyBoundaries.RtElutionRange.Minimum, Is.EqualTo(0).Within(1e-5));
            Assert.That(emptyBoundaries.RtElutionRange.Maximum, Is.EqualTo(0).Within(1e-5));
        }

        [Test]
        public void TestFullSequenceIsMappedDirectly()
        {
            var species = ProcessEdge();

            Assert.That(species[1].SpeciesLabel, Is.Empty);
            Assert.That(species[1].Identification, Is.Null);
            Assert.That(species[2].SpeciesLabel, Is.Empty);
            Assert.That(species[2].Identification, Is.Null);
        }

        [Test]
        public void TestSourceFilenameNormalization()
        {
            var species = ProcessEdge();

            Assert.That(species[0].SpectraFileNameWithoutExtension, Is.EqualTo("edgeraw"));
            Assert.That(species[0].DeconvolutionFeature.SpectraFileNameWithoutExtension, Is.EqualTo("edgeraw"));
        }

        [Test]
        public void TestMzLibRowsAreMappedWithoutLocalValidation()
        {
            var species = ProcessEdge();

            Assert.That(species, Has.Count.EqualTo(9));
        }
    }
}
