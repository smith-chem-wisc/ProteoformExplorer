using MassSpectrometry;
using ProteoformExplorer.Core;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using ScottPlot.Plottable;
using System;
using System.Linq;
using System.Drawing;
using System.Collections.Generic;
using ProteoformExplorer.GuiFunctions;

namespace ProteoformExplorer.Wpf
{
    public class TicPageViewModel : BaseViewModel
    {
        private static double? _minRt;
        private static double? _maxRt;

        public double? MinRt
        {
            get => _minRt;
            set
            {
                if (_minRt != value)
                {
                    _minRt = value;
                    OnPropertyChanged(nameof(MinRt));
                }
            }
        }

        public double? MaxRt
        {
            get => _maxRt;
            set
            {
                if (_maxRt != value)
                {
                    _maxRt = value;
                    OnPropertyChanged(nameof(MaxRt));
                }
            }
        }

        public void InitializeMinMaxRt()
        {
            if (_minRt.HasValue && _maxRt.HasValue)
                return;

            MinRt = 0;
            MaxRt = DataManagement.SpectraFiles.Select(p => p.Value.DataFile.Value.Scans.Last().RetentionTime).Max();
        }
    }

    /// <summary>
    /// Interaction logic for Page1_QuantifiedTic.xaml
    /// </summary>
    public partial class Page1_QuantifiedTic : Page
    {
        private MsDataScan CurrentScan;
        private VLine IntegratedAreaStart;
        private VLine IntegratedAreaEnd;
        private VLine CurrentRtIndicator;
        private Text PercentDeconAnnotation;
        private Text PercentIdentifiedAnnotation;
        private Text SelectedSpectrumItemAnnotation;
        private readonly TicPageViewModel _viewModel = new();

        public Page1_QuantifiedTic()
        {
            InitializeComponent();
            DataContext = _viewModel;
            DataListView.ItemsSource = DataLoading.LoadedSpectraFiles;
            Loaded += Page1_QuantifiedTic_Loaded;
            _viewModel.InitializeMinMaxRt();

            // right click to zoom is replaced by right click to integrate for this chart
            topPlotView.Configuration.RightClickDragZoom = false;

            PercentDeconAnnotation = new Text();
            PercentDeconAnnotation.Color = GuiSettings.DeconvolutedColor;
            PercentDeconAnnotation.FontSize = 16;

            PercentIdentifiedAnnotation = new Text();
            PercentIdentifiedAnnotation.Color = GuiSettings.IdentifiedColor;
            PercentIdentifiedAnnotation.FontSize = 16;
        }

        private void Page1_QuantifiedTic_Loaded(object sender, RoutedEventArgs e)
        {
            if (DataLoading.LoadedSpectraFiles?.Any() != true)
            {
                return;
            }

            DataListView.SelectedItem = DataLoading.LoadedSpectraFiles.First();

            if (DataManagement.CurrentlySelectedFile.Value == null)
            {
                WpfFunctions.OnSpectraFileChanged(DataListView, null);
                DisplayTic();
            }
        }

        private void Home_Click(object sender, RoutedEventArgs e)
        {
            this.NavigationService.Navigate(new Dashboard());
        }

        private void DisplayTic()
        {
            PlottingFunctions.PlotTotalIonChromatograms(topPlotView.Plot, CurrentRtIndicator, out var errors, _viewModel.MinRt, _viewModel.MaxRt);

            if (errors.Any())
            {
                MessageBox.Show(errors.First());
            }
        }

        private void DisplayAnnotatedSpectrum(int scanNum)
        {
            var scan = DataManagement.CurrentlySelectedFile.Value.GetOneBasedScan(scanNum);

            if (scan == null)
            {
                return;
            }

            var speciesInScan = DataManagement.CurrentlySelectedFile.Value.GetSpeciesInScan(scanNum);
            CurrentScan = scan;
            CurrentRtIndicator = PlottingFunctions.UpdateRtIndicator(scan, CurrentRtIndicator, topPlotView.Plot);
            topPlotView.Render();

            PlottingFunctions.PlotSpeciesInSpectrum(speciesInScan, scanNum, DataManagement.CurrentlySelectedFile, bottomPlotView.Plot, out scan, out var errors);

            if (errors.Any())
            {
                MessageBox.Show(errors.First());
            }
        }

        private void Window_PreviewKeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Left)
            {
                if (CurrentScan != null)
                {
                    DisplayAnnotatedSpectrum(CurrentScan.OneBasedScanNumber - 1);
                }
                e.Handled = true;
            }
            else if (e.Key == Key.Right)
            {
                if (CurrentScan != null)
                {
                    DisplayAnnotatedSpectrum(CurrentScan.OneBasedScanNumber + 1);
                }
                e.Handled = true;
            }
        }

        private void topPlotView_PreviewMouseLeftOrRightButtonDown(object sender, MouseButtonEventArgs e)
        {
            if (DataManagement.CurrentlySelectedFile.Value == null)
            {
                return;
            }

            double rt = WpfFunctions.GetXPositionFromMouseClickOnChart(sender, e);
            var theScan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(DataManagement.CurrentlySelectedFile, rt, 1);
            DisplayAnnotatedSpectrum(theScan.OneBasedScanNumber);

            IntegratedAreaStart = new VLine();
            IntegratedAreaStart.X = theScan.RetentionTime;

            ClearIntegratedAreasAndTextAnnotations();
        }

        private void DataListView_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            WpfFunctions.OnSpectraFileChanged(sender, e);

            DisplayTic();
        }

        private void openFileListViewButton_Click(object sender, RoutedEventArgs e)
        {
            WpfFunctions.ShowOrHideSpectraFileList(DataListView, gridSplitter);
        }

        private void topPlotView_PreviewMouseRightButtonUp(object sender, MouseButtonEventArgs e)
        {
            if (IntegratedAreaStart == null)
            {
                return;
            }

            IntegratedAreaEnd = new VLine();
            var loc = WpfFunctions.GetXPositionFromMouseClickOnChart(sender, e);
            IntegratedAreaEnd.X = loc;

            // this is in case the user dragged right-to-left instead of left-to-right
            double integratedAreaStart = Math.Min(IntegratedAreaStart.X, IntegratedAreaEnd.X);
            double integratedAreaEnd = Math.Max(IntegratedAreaStart.X, IntegratedAreaEnd.X);

            // shade the area
            if (Math.Abs(integratedAreaEnd - integratedAreaStart) > 0.001)
            {
                List<(string label, double summedTic, double xLoc, double yLoc)> integrations = new List<(string, double, double, double)>();

                foreach (var item in topPlotView.Plot.GetPlottables())
                {
                    if (item is ScatterPlot scatter)
                    {
                        double summedTic = 0;

                        var pointsX = scatter.Xs;
                        var pointsY = scatter.Ys;

                        int startInd = Array.IndexOf(pointsX, pointsX.First(p => p >= integratedAreaStart));
                        if (startInd > 0 && Math.Abs(pointsX[startInd - 1] - integratedAreaStart) < Math.Abs(pointsX[startInd] - integratedAreaStart))
                        {
                            startInd--;
                        }

                        int endInd = integratedAreaEnd > pointsX.Last() ? pointsX.Length - 1 : Array.IndexOf(pointsX, pointsX.First(p => p >= integratedAreaEnd));
                        if (endInd > 0 && Math.Abs(pointsX[endInd - 1] - integratedAreaEnd) < Math.Abs(pointsX[endInd] - integratedAreaEnd))
                        {
                            endInd--;
                        }

                        if (startInd == endInd)
                        {
                            continue;
                        }

                        double[] xs = new double[endInd - startInd + 3];
                        double[] ys = new double[endInd - startInd + 3];

                        xs[0] = pointsX[startInd];
                        ys[0] = 0;

                        for (int i = startInd; i <= endInd; i++)
                        {
                            xs[i - startInd + 1] = pointsX[i];
                            ys[i - startInd + 1] = pointsY[i];
                        }

                        xs[xs.Length - 1] = pointsX[endInd];
                        ys[xs.Length - 1] = 0;

                        var color = Color.FromArgb(GuiSettings.IntegrationFillAlpha, scatter.Color);
                        var poly = topPlotView.Plot.AddPolygon(xs, ys, fillColor: color);
                        summedTic = ys.Sum();

                        integrations.Add((scatter.Label, summedTic, xs[1], ys[1]));
                    }
                }

                if (integrations.Any(p => p.label == "Deconvoluted TIC"))
                {
                    var tic = integrations.First(p => p.label == "TIC");
                    var deconTic = integrations.First(p => p.label == "Deconvoluted TIC");
                    double percentTicDeconvoluted = deconTic.summedTic / tic.summedTic;

                    var axisDims = topPlotView.Plot.GetAxisLimits();
                    PercentDeconAnnotation.Label = (percentTicDeconvoluted * 100).ToString("F1") + "%";
                    PercentDeconAnnotation.X = deconTic.xLoc;
                    PercentDeconAnnotation.Y = axisDims.YMax;
                    PercentDeconAnnotation.FontSize = (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingX);
                    PercentDeconAnnotation.FontBold = true;

                    if (!topPlotView.Plot.GetPlottables().Contains(PercentDeconAnnotation))
                    {
                        topPlotView.Plot.Add(PercentDeconAnnotation);
                    }
                }
                else if (integrations.Any(p => p.label == "Identified TIC"))
                {
                    var tic = integrations.First(p => p.label == "TIC");
                    var identTic = integrations.First(p => p.label == "Identified TIC");
                    double percentTicIdentified = identTic.summedTic / tic.summedTic;

                    var axisDims = topPlotView.Plot.GetAxisLimits();
                    PercentIdentifiedAnnotation.Label = (percentTicIdentified * 100).ToString("F1") + "%";
                    PercentIdentifiedAnnotation.X = identTic.xLoc;
                    PercentIdentifiedAnnotation.Y = axisDims.YCenter;
                    PercentIdentifiedAnnotation.FontSize = (float)(GuiSettings.ChartAxisLabelFontSize * GuiSettings.DpiScalingX);
                    PercentDeconAnnotation.FontBold = true;

                    if (!topPlotView.Plot.GetPlottables().Contains(PercentIdentifiedAnnotation))
                    {
                        topPlotView.Plot.Add(PercentIdentifiedAnnotation);
                    }
                }
            }
        }

        private void ClearIntegratedAreasAndTextAnnotations()
        {
            topPlotView.Plot.Remove(PercentDeconAnnotation);
            topPlotView.Plot.Remove(PercentIdentifiedAnnotation);

            List<IPlottable> itemsToRemove = topPlotView.Plot.GetPlottables().Where(p => p is Polygon).ToList();

            foreach (var item in itemsToRemove)
            {
                topPlotView.Plot.Remove(item);
            }
        }

        private void bottomPlotView_PreviewMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            var xval = WpfFunctions.GetXPositionFromMouseClickOnChart(sender, e);
            SelectedSpectrumItemAnnotation = PlottingFunctions.OnSpectrumPeakSelected(xval, CurrentScan, SelectedSpectrumItemAnnotation);

            if (!bottomPlotView.Plot.GetPlottables().Contains(SelectedSpectrumItemAnnotation))
            {
                bottomPlotView.Plot.Add(SelectedSpectrumItemAnnotation);
            }
        }
    }
}
