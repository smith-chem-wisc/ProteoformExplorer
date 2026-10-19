using MassSpectrometry;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using ProteoformExplorer.Core;
using System.Windows;
using System.Windows.Controls;
using System.Linq;
using System.IO;
using ProteoformExplorer.GuiFunctions;

namespace ProteoformExplorer.Wpf
{
    /// <summary>
    /// Interaction logic for Page3_StackedIons.xaml
    /// </summary>
    public partial class Page3_StackedIons : Page
    {
        private ObservableCollection<INode> SelectableAnnotatedSpecies;
        private const double mmPerInch = 25.4;

        public Page3_StackedIons()
        {
            InitializeComponent();
            DataListView.ItemsSource = DataLoading.LoadedSpectraFiles;

            SelectableAnnotatedSpecies = new ObservableCollection<INode>();
            SpeciesListView.ItemsSource = SelectableAnnotatedSpecies;
        }

        public void SelectItem(AnnotatedSpecies species)
        {
            if (species == null)
            {
                return;
            }

            // select the spectra file
            DataListView.SelectedItem = DataLoading.LoadedSpectraFiles.First(p =>
                p.LowerFileNameWithoutExtensions == species.SpectraFileNameWithoutExtension);

            // notify the data manager that the spectra file has changed
            WpfFunctions.OnSpectraFileChanged(DataListView, null);

            // this will populate the treeview for species found in the currently selected file
            WpfFunctions.PopulateTreeViewWithSpeciesAndCharges(SelectableAnnotatedSpecies);

            // select the item in the treeview
            var theNodeToSelect = SelectableAnnotatedSpecies.First(p => p is AnnotatedSpeciesNode node && node.AnnotatedSpecies == species);
            AnnotatedSpeciesNode node = (AnnotatedSpeciesNode)theNodeToSelect;
            AnnotatedSpeciesNodeSpecificCharge child = (AnnotatedSpeciesNodeSpecificCharge)node.Charges.First();

            node.IsExpanded = true;
            child.IsSelected = true;

            // make the plot
            PlotSpeciesStacked(child.AnnotatedSpecies, child.Charge);
            GuiFunctions.PlottingFunctions.OnSpeciesChanged();
        }

        private void Home_Click(object sender, RoutedEventArgs e)
        {
            this.NavigationService.Navigate(new Dashboard());
        }

        private void PlotSpeciesStacked(AnnotatedSpecies species, int? charge = null)
        {
            List<int> chargesToPlot = new List<int>();
            double modeMass;
            MsDataScan initialScan;

            // get apex or precursor scan
            if (species.DeconvolutionFeature != null)
            {
                initialScan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(DataManagement.CurrentlySelectedFile, species.DeconvolutionFeature.ApexRt, 1);
                modeMass = PfmXplorerUtil.DeconvolutionEngine.GetModeMassFromMonoisotopicMass(species.DeconvolutionFeature.MonoisotopicMass);
            }
            else
            {
                initialScan = DataManagement.CurrentlySelectedFile.Value.GetOneBasedScan(species.Identification.OneBasedPrecursorScanNumber);
                modeMass = PfmXplorerUtil.DeconvolutionEngine.GetModeMassFromMonoisotopicMass(species.Identification.MonoisotopicMass);
            }

            // decide on charges to plot
            if (charge == null)
            {
                if (species.DeconvolutionFeature != null)
                {
                    chargesToPlot.AddRange(species.DeconvolutionFeature.Charges);
                }
                else if (species.Identification != null)
                {
                    chargesToPlot.Add(species.Identification.PrecursorChargeState);
                }
            }
            else
            {
                chargesToPlot.Add(charge.Value);
            }

            double? globalMinRt = null;
            double? globalMaxRt = null;
            if (charge != null)
            {
                foreach (var file in DataManagement.SpectraFiles.Values)
                {
                    var scans = file.GetScansInRtWindow(initialScan.RetentionTime, GuiFunctions.GuiSettings.RtExtractionWindow);
                    if (!scans.Any())
                    {
                        continue;
                    }

                    globalMinRt = globalMinRt.HasValue
                        ? System.Math.Min(globalMinRt.Value, scans.First().RetentionTime)
                        : scans.First().RetentionTime;
                    globalMaxRt = globalMaxRt.HasValue
                        ? System.Math.Max(globalMaxRt.Value, scans.Last().RetentionTime)
                        : scans.Last().RetentionTime;
                }
            }

            int fileNum = 0;
            foreach (var file in DataManagement.SpectraFiles.OrderBy(p => p.Key.ConvertName()))
            {
                if (charge == null)
                {
                    // plot summed charge states, one line per file
                    //TODO
                    topPlotView.Plot.Clear();
                }
                else
                {
                    // plot the charge state envelope, one line per file
                    GuiFunctions.PlottingFunctions.PlotSummedChargeStateXic(modeMass, charge.Value, initialScan.RetentionTime, GuiFunctions.GuiSettings.RtExtractionWindow, file, topPlotView.Plot,
                        fileNum == 0, null, out var errors, 0, 0, true, file.Key.ConvertName(), globalMinRt, globalMaxRt);

                    if (errors.Any())
                    {
                        MessageBox.Show(errors.First());
                    }

                    fileNum++;
                }
            }

            double inchOffset = GuiFunctions.GuiSettings.WaterfallSpacing / mmPerInch;
            var axisLimits = topPlotView.Plot.GetAxisLimits(topPlotView.Plot.XAxis.AxisIndex, topPlotView.Plot.YAxis.AxisIndex);

            // y axis offset
            double yAxisRange = topPlotView.Plot.GetPixelY(axisLimits.YMax) - topPlotView.Plot.GetPixelY(axisLimits.YMin);
            double yUnitsPerDot = (axisLimits.YMax - axisLimits.YMin) / yAxisRange;
            double offsetUnitStepY = inchOffset * GuiFunctions.GuiSettings.DpiScalingY * 96 * yUnitsPerDot;

            // x axis offset
            double xAxisRange = topPlotView.Plot.GetPixelX(axisLimits.XMax) - topPlotView.Plot.GetPixelX(axisLimits.XMin);
            double xUnitsPerDot = (axisLimits.XMax - axisLimits.XMin) / xAxisRange;
            double offsetUnitStepX = inchOffset * GuiFunctions.GuiSettings.DpiScalingX * 96 * xUnitsPerDot;

            fileNum = 0;
            double xOffset = 0;
            double yOffset = 0;
            foreach (var file in DataManagement.SpectraFiles.OrderBy(p => p.Key.ConvertName()))
            {
                if (charge == null)
                {
                    // plot summed charge states, one line per file
                    //TODO
                    topPlotView.Plot.Clear();
                }
                else
                {
                    // plot the charge state envelope, one line per file
                    GuiFunctions.PlottingFunctions.PlotSummedChargeStateXic(modeMass, charge.Value, initialScan.RetentionTime, GuiFunctions.GuiSettings.RtExtractionWindow, file, topPlotView.Plot,
                        fileNum == 0, null, out var errors, xOffset, yOffset, true, file.Key.ConvertName(), globalMinRt, globalMaxRt);

                    if (errors.Any())
                    {
                        MessageBox.Show(errors.First());
                    }

                    fileNum++;
                    xOffset -= offsetUnitStepX;
                    yOffset -= offsetUnitStepY;
                }
            }

            // zoom axes
            axisLimits = topPlotView.Plot.GetAxisLimits(topPlotView.Plot.XAxis.AxisIndex, topPlotView.Plot.YAxis.AxisIndex);
            topPlotView.Plot.SetAxisLimits(
                axisLimits.XMin - axisLimits.XSpan / 2,
                axisLimits.XMax + axisLimits.XSpan / 2,
                axisLimits.YMin - axisLimits.YSpan / 2,
                axisLimits.YMax + axisLimits.YSpan / 2);
        }

        private void SpeciesListView_SelectedItemChanged(object sender, RoutedPropertyChangedEventArgs<object> e)
        {
            INode species = ((TreeView)sender).SelectedItem as INode;

            if (species is AnnotatedSpeciesNode annotatedSpeciesNode)
            {
                PlotSpeciesStacked(annotatedSpeciesNode.AnnotatedSpecies);
            }
            else if (species is AnnotatedSpeciesNodeSpecificCharge annotatedSpeciesNodeSpecificCharge)
            {
                PlotSpeciesStacked(annotatedSpeciesNodeSpecificCharge.AnnotatedSpecies, annotatedSpeciesNodeSpecificCharge.Charge);
            }

            GuiFunctions.PlottingFunctions.OnSpeciesChanged();
        }

        private void DataListView_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            WpfFunctions.OnSpectraFileChanged(sender, e);

            WpfFunctions.PopulateTreeViewWithSpeciesAndCharges(SelectableAnnotatedSpecies);
        }

        private void openFileListViewButton_Click(object sender, RoutedEventArgs e)
        {
            WpfFunctions.ShowOrHideSpectraFileList(DataListView, gridSplitter);
        }

        private void SpeciesListView_PreviewMouseDoubleClick(object sender, System.Windows.Input.MouseButtonEventArgs e)
        {
            // double click doesn't actually do anything. the program seems to crash when double-clicking, so 
            // this event handler is here just to do nothing if a double-click occurs
        }

        // https://stackoverflow.com/questions/3225940/prevent-automatic-horizontal-scroll-in-treeview/34269542
        private void TreeViewItem_RequestBringIntoView(object sender, RequestBringIntoViewEventArgs e)
        {
            e.Handled = true;
        }
    }
}
