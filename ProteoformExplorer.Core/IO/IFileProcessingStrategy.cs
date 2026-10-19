#nullable enable
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ProteoformExplorer.Core.IO;

public interface IFileProcessingStrategy
{
    List<AnnotatedSpecies> ProcessFile(string filePath, string dataset = "");
}

public class MetaMorpheusProcessingStrategy : IFileProcessingStrategy
{
    public List<AnnotatedSpecies> ProcessFile(string filePath, string dataset = "")
    {
        var psms = SpectrumMatchTsvReader.ReadTsv(filePath, out _);

        List<AnnotatedSpecies> species = new List<AnnotatedSpecies>(psms.Count);
        foreach (var spectralMatch in psms)
        {
            var id = new Identification(spectralMatch.BaseSequence, spectralMatch.FullSequence, spectralMatch.PrecursorMass, spectralMatch.ChargeState,
                spectralMatch.PrecursorScanNum, spectralMatch.Ms2ScanNumber, spectralMatch.FileNameWithoutExtension, dataset);

            species.Add(new AnnotatedSpecies(id));
        }
        return species;
    }
}

public class Ms1FeatureProcessingStrategy : IFileProcessingStrategy
{
    public List<AnnotatedSpecies> ProcessFile(string filePath, string dataset = "")
    {
        var features = FileReader.ReadFile<Ms1FeatureFile>(filePath).Results;
        List<AnnotatedSpecies> species = new List<AnnotatedSpecies>(features.Count);

        // TODO: File name may have issues here with feature alignment for plots. This parses the file name of the align file, not the mass spec file. 
        var fileName = PfmXplorerUtil.GetFileNameWithoutExtension(filePath.Replace("_ms1.feature", ""));
        foreach (var feature in features)
        {
            var decon = new DeconvolutionFeature(feature.Mass, feature.RetentionTimeApex / 60.0, feature.RetentionTimeBegin / 60.0, feature.RetentionTimeEnd / 60.0,
                Enumerable.Range(feature.ChargeStateMin, feature.ChargeStateMax - feature.ChargeStateMin + 1).ToList(), fileName);
            var specie = new AnnotatedSpecies(decon, feature.Id.ToString());
            species.Add(specie);
        }
        return species;
    }
}

public class FlashDeconvTsvProcessingStrategy : IFileProcessingStrategy
{
    public List<AnnotatedSpecies> ProcessFile(string filePath, string dataset = "")
    {
        var features = FileReader.ReadFile<FlashDeconvTsvFile>(filePath).Results;

        List<AnnotatedSpecies> species = new List<AnnotatedSpecies>(features.Count);
        foreach (var feature in features)
        {
            var decon = new DeconvolutionFeature(feature.MonoisotopicMass, feature.RetentionTimeApex / 60.0, feature.RetentionTimeBegin / 60.0, feature.RetentionTimeEnd / 60.0,
                Enumerable.Range(feature.ChargeStateMin, feature.ChargeStateMax - feature.ChargeStateMin + 1).ToList(), feature.FilePath);
            var specie = new AnnotatedSpecies(decon, feature.FeatureIndex.ToString());
            species.Add(specie);
        }
        return species;
    }
}

public class QuantifiedPeakProcessingStrategy : IFileProcessingStrategy
{
    public List<AnnotatedSpecies> ProcessFile(string filePath, string dataset = "")
    {
        var rows = FileReader.ReadFile<QuantifiedPeakFile>(filePath).Results;

        var groupedRows = rows.GroupBy(row => new
        {
            FileName = PfmXplorerUtil.GetFileNameWithoutExtension(row.FileName),
            row.PeptideMonoisotopicMass,
            row.FullSequence
        });

        List<AnnotatedSpecies> species = new List<AnnotatedSpecies>(groupedRows.Count());
        foreach (var group in groupedRows)
        {
            var first = group.First();
            var charges = group
                .Select(row => row.PeakCharge is int peakCharge && peakCharge != 0 ? peakCharge : row.PrecursorCharge)
                .Distinct()
                .ToList();
            double rtStart = group.Select(row => row.PeakRTStart.GetValueOrDefault()).Min();
            double rtEnd = group.Select(row => row.PeakRTEnd.GetValueOrDefault()).Max();

            var deconvolutionFeature = new DeconvolutionFeature(first.PeptideMonoisotopicMass,
                first.PeakRTApex.GetValueOrDefault(), rtStart, rtEnd, charges, group.Key.FileName);

            species.Add(new AnnotatedSpecies(deconvolutionFeature, group.Key.FullSequence));
        }

        return species;
    }
}
