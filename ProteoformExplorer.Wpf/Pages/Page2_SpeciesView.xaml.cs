using Chemistry;
using MassSpectrometry;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using ProteoformExplorer.Core;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using ScottPlot.Plottable;
using System.Windows.Threading;
using System;
using ProteoformExplorer.GuiFunctions;

namespace ProteoformExplorer.Wpf
{
    /// <summary>
    /// Interaction logic for Page2_SpeciesView.xaml
    /// </summary>
    public partial class Page2_SpeciesView : Page
    {
        private ObservableCollection<INode> SelectableAnnotatedSpecies;
        private MsDataScan CurrentScan;
        private AnnotatedSpecies CurrentlyDisplayedSpecies;
        private int? CurrentlyDisplayedCharge;
        private VLine CurrentRtIndicator;
        private Text SelectedSpectrumItemAnnotation;

        public Page2_SpeciesView()
        {
            InitializeComponent();
            DataListView.ItemsSource = DataLoading.LoadedSpectraFiles;

            SelectableAnnotatedSpecies = new ObservableCollection<INode>();
            SpeciesListView.ItemsSource = SelectableAnnotatedSpecies;
        }

        private void Home_Click(object sender, RoutedEventArgs e)
        {
            this.NavigationService.Navigate(new Dashboard());
        }

        private void PlotSpecies(AnnotatedSpecies species, int? charge = null)
        {
            PlotSpeciesIsotopeXics(species, out var apexScan, charge);
            PlotSpeciesInSpectrum(species, apexScan, charge);
            CurrentlyDisplayedSpecies = species;
            CurrentlyDisplayedCharge = charge;
        }

        private void PlotSpeciesIsotopeXics(AnnotatedSpecies species, out MsDataScan initialScan, int? charge = null)
        {
            List<int> chargesToPlot = new List<int>();
            double modeMass;

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
                    chargesToPlot.AddRange(species.DeconvolutionFeature.Charges.OrderBy(p => p));
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

            if (charge == null)
            {
                // plot summed isotopes, one line per charge
                for (int i = 0; i < chargesToPlot.Count; i++)
                {
                    int z = chargesToPlot[i];

                    GuiFunctions.PlottingFunctions.PlotSummedChargeStateXic(modeMass, z, initialScan.RetentionTime, GuiFunctions.GuiSettings.RtExtractionWindow,
                        DataManagement.CurrentlySelectedFile, topPlotView.Plot, i == 0, CurrentRtIndicator, out var errors, 0, 0, false, "z=" + z);

                    if (errors.Any())
                    {
                        MessageBox.Show(errors.First());
                    }
                }
            }
            else
            {
                // plot isotopes, one line per isotope

                // decide which isotopes to plot
                List<(double mz, int z)> peaksToMakeXicsFor = new List<(double mz, int z)>();
                foreach (var z in chargesToPlot)
                {
                    double mz = modeMass.ToMz(z);
                    int ind = initialScan.MassSpectrum.GetClosestPeakIndex(mz);

                    var env = PfmXplorerUtil.DeconvolutionEngine.GetIsotopicEnvelope(initialScan.MassSpectrum, ind, z,
                        new List<Deconvoluter.DeconvolutedPeak>(), new HashSet<double>(), new List<(double, double)>());

                    if (env != null)
                    {
                        peaksToMakeXicsFor.AddRange(env.Peaks.OrderBy(p => p.ExperimentalMz).Select(p => (p.ExperimentalMz, p.Charge)));
                    }
                    else
                    {
                        peaksToMakeXicsFor.Add((mz, z));
                    }
                }

                double xOffset = 0;
                double yOffset = 0;
                if (GuiFunctions.GuiSettings.WaterfallXics)
                {
                    //TODO
                }

                // make the plots
                for (int i = 0; i < peaksToMakeXicsFor.Count; i++)
                {
                    var peak = peaksToMakeXicsFor[i];

                    GuiFunctions.PlottingFunctions.PlotXic(peak.mz, peak.z, PfmXplorerUtil.DeconvolutionEngine.PpmTolerance, initialScan.RetentionTime, GuiFunctions.GuiSettings.RtExtractionWindow,
                        DataManagement.CurrentlySelectedFile, topPlotView.Plot, i == 0, CurrentRtIndicator, out var errors, xOffset, yOffset, peak.mz.ToMass(peak.z).ToString("F2"));

                    if (errors.Any())
                    {
                        MessageBox.Show(errors.First());
                    }
                }
            }
        }

        private void PlotSpeciesInSpectrum(AnnotatedSpecies species, MsDataScan scan, int? charge = null)
        {
            PlottingFunctions.PlotSpeciesInSpectrum(new List<AnnotatedSpecies> { species }, scan.OneBasedScanNumber, DataManagement.CurrentlySelectedFile,
                bottomPlotView.Plot, out var scan2, out var errors, charge);

            if (errors.Any())
            {
                MessageBox.Show(errors.First());
            }

            CurrentScan = scan;

            CurrentRtIndicator = PlottingFunctions.UpdateRtIndicator(scan, CurrentRtIndicator, topPlotView.Plot);
            topPlotView.Render();
        }

        private void Window_PreviewKeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Left)
            {
                if (CurrentScan != null && CurrentlyDisplayedSpecies != null)
                {
                    var previousScan = PfmXplorerUtil.GetAdjacentScanByMsOrder(
                        DataManagement.CurrentlySelectedFile,
                        CurrentScan.OneBasedScanNumber,
                        -1,
                        1);

                    if (previousScan != null)
                    {
                        PlotSpeciesInSpectrum(CurrentlyDisplayedSpecies, previousScan, CurrentlyDisplayedCharge);
                    }
                }
                e.Handled = true;
            }
            else if (e.Key == Key.Right)
            {
                if (CurrentScan != null && CurrentlyDisplayedSpecies != null)
                {
                    var nextScan = PfmXplorerUtil.GetAdjacentScanByMsOrder(
                        DataManagement.CurrentlySelectedFile,
                        CurrentScan.OneBasedScanNumber,
                        1,
                        1);

                    if (nextScan != null)
                    {
                        PlotSpeciesInSpectrum(CurrentlyDisplayedSpecies, nextScan, CurrentlyDisplayedCharge);
                    }
                }
                e.Handled = true;
            }
        }

        private void topPlotView_PreviewMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            if (DataManagement.CurrentlySelectedFile.Value == null)
            {
                return;
            }

            double rt = WpfFunctions.GetXPositionFromMouseClickOnChart(sender, e);
            var theScan = PfmXplorerUtil.GetClosestScanToRtFromDynamicConnection(DataManagement.CurrentlySelectedFile, rt);

            PlotSpeciesInSpectrum(CurrentlyDisplayedSpecies, theScan, CurrentlyDisplayedCharge);
        }

        private void SpeciesListView_SelectedItemChanged(object sender, RoutedPropertyChangedEventArgs<object> e)
        {
            INode species = ((TreeView)sender).SelectedItem as INode;

            if (species is AnnotatedSpeciesNode annotatedSpeciesNode)
            {
                PlotSpecies(annotatedSpeciesNode.AnnotatedSpecies);
            }
            else if (species is AnnotatedSpeciesNodeSpecificCharge annotatedSpeciesNodeSpecificCharge)
            {
                PlotSpecies(annotatedSpeciesNodeSpecificCharge.AnnotatedSpecies, annotatedSpeciesNodeSpecificCharge.Charge);
            }

            PlottingFunctions.OnSpeciesChanged();
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

        private void SpeciesListView_PreviewMouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            // double click doesn't actually do anything. the program seems to crash when double-clicking, so 
            // this event handler is here just to do nothing if a double-click occurs
        }

        // https://stackoverflow.com/questions/3225940/prevent-automatic-horizontal-scroll-in-treeview/34269542
        private void TreeViewItem_RequestBringIntoView(object sender, RequestBringIntoViewEventArgs e)
        {
            e.Handled = true;
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
