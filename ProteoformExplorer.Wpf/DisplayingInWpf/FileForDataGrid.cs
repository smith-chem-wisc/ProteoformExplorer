using ProteoformExplorer.Core;
using ProteoformExplorer.Core.IO;
using ProteoformExplorer.GuiFunctions;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Windows.Media;

namespace ProteoformExplorer.Wpf
{
    public enum FileType { Identification, Deconvolution, Spectra, Unknown }

    public class FileForDataGrid
    {
        public bool IsLoaded = false;
        public string FullFilePath { get; set; }
        public string FileNameWithExtension { get; set; }
        public string LowerFileNameWithoutExtensions { get; set; }
        public FileType FileType { get; set; }
        public SolidColorBrush BackgroundColor { get; private set; }
        public SolidColorBrush ForegroundColor { get; private set; }

        public FileForDataGrid(string fullFilePath, out List<string> errors)
        {
            FullFilePath = fullFilePath;
            FileNameWithExtension = Path.GetFileName(fullFilePath);
            LowerFileNameWithoutExtensions = PfmXplorerUtil.GetFileNameWithoutExtension(fullFilePath);
            DetermineFileType(out errors);
        }

        private void DetermineFileType(out List<string> errors)
        {
            errors = new List<string>();
            var extension = Path.GetExtension(FullFilePath).ToLowerInvariant();

            if (InputReaderParser.AcceptedSpectraFileFormats.Contains(extension))
            {
                FileType = FileType.Spectra;
            }
            else if (InputReaderParser.AcceptedTextFileFormats.Contains(extension))
            {
                StreamReader reader;
                try
                {
                    using (reader = new StreamReader(FullFilePath))
                    {
                        string line = reader.ReadLine();

                        SupportedFileType? mzLibFileType = null;
                        try
                        {
                            mzLibFileType = FullFilePath.ParseFileType();
                        }
                        catch
                        {
                        }

                        if (mzLibFileType == SupportedFileType.FlashLFQQuantifiedPeak)
                        {
                            FileType = FileType.Deconvolution;
                        }
                        else
                        {
                            var fileType = InputReaderParser.GetFileTypeFromHeader(line);

                            if (fileType == InputReaderParser.InputSourceType.ThermoDecon || fileType == InputReaderParser.InputSourceType.Promex
                                || fileType == InputReaderParser.InputSourceType.FlashDeconv || fileType == InputReaderParser.InputSourceType.ProteoformExplorer
                                || fileType == InputReaderParser.InputSourceType.TopFD)
                            {
                                FileType = FileType.Deconvolution;
                            }
                            else if (fileType == InputReaderParser.InputSourceType.MetaMorpheus || fileType == InputReaderParser.InputSourceType.TDPortal)
                            {
                                FileType = FileType.Identification;
                            }
                            else
                            {
                                FileType = FileType.Unknown;
                            }
                        }
                    }
                }
                catch (Exception e)
                {
                    FileType = FileType.Unknown;
                    errors.Add("Error reading file: " + e.Message);
                }
            }
            else
            {
                FileType = FileType.Unknown;
                errors.Add("Error reading file: The extension " + extension + " is unknown");
            }

            SolidColorBrush brush = Brushes.White;

            switch (FileType)
            {
                case FileType.Deconvolution:
                    var color = Color.FromArgb(
                            GuiSettings.DeconvolutedColor.A,
                            GuiSettings.DeconvolutedColor.R,
                            GuiSettings.DeconvolutedColor.G,
                            GuiSettings.DeconvolutedColor.B);
                    brush = new SolidColorBrush(color);
                    break;
                case FileType.Identification:
                    color = Color.FromArgb(
                            GuiSettings.IdentifiedColor.A,
                            GuiSettings.IdentifiedColor.R,
                            GuiSettings.IdentifiedColor.G,
                            GuiSettings.IdentifiedColor.B);
                    brush = new SolidColorBrush(color);
                    break;
                case FileType.Spectra:
                    color = Color.FromArgb(
                            GuiSettings.TicColor.A,
                            GuiSettings.TicColor.R,
                            GuiSettings.TicColor.G,
                            GuiSettings.TicColor.B);
                    brush = new SolidColorBrush(color);
                    break;
                case FileType.Unknown:
                    color = Color.FromArgb(
                            Colors.Red.A,
                            Colors.Red.R,
                            Colors.Red.G,
                            Colors.Red.B);
                    brush = new SolidColorBrush(color);
                    break;
                default:
                    color = Color.FromArgb(
                            Colors.Red.A,
                            Colors.Red.R,
                            Colors.Red.G,
                            Colors.Red.B);
                    brush = new SolidColorBrush(color);
                    break;
            }

            BackgroundColor = brush;
            ForegroundColor = Brushes.White;
        }

        public override string ToString()
        {
            return FileNameWithExtension.ConvertName();
        }
    }
}
