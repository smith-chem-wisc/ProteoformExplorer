#nullable enable
using Readers;
using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
namespace ProteoformExplorer.Core.IO;

public class FileProcessor
{
    private readonly IFileProcessingStrategy _strategy;

    public FileProcessor(IFileProcessingStrategy strategy)
    {
        _strategy = strategy;
    }

    public List<AnnotatedSpecies> Process(string filePath)
    {
        string? dataset = null;
        if (_strategy is MetaMorpheusProcessingStrategy)
        {
            // assumes dataset name is a component of the file path: the name of the main search output in which the psmtsv file is located
            dataset = Path.GetFileNameWithoutExtension(Path.GetDirectoryName(Path.GetDirectoryName(filePath)));
        }
        return _strategy.ProcessFile(filePath, dataset ?? "");
    }

    public static FileProcessor? GetProcessor(string filePath)
    {
        var type = filePath.ParseFileType();
        switch (type)
        {
            case SupportedFileType.psmtsv:
            case SupportedFileType.osmtsv:
                return new FileProcessor(new MetaMorpheusProcessingStrategy());
            case SupportedFileType.FlashLFQQuantifiedPeak:
                return new FileProcessor(new QuantifiedPeakProcessingStrategy());
            case SupportedFileType.Ms1Feature:
                return new FileProcessor(new Ms1FeatureProcessingStrategy());
            case SupportedFileType.Tsv_FlashDeconv:
                return new FileProcessor(new FlashDeconvTsvProcessingStrategy());
            default:
                return null;
        }
    }
}
