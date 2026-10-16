using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ProteoformExplorer.Core
{
    public static class DataManagement
    {
        public static Dictionary<string, CachedSpectraFileData> SpectraFiles;
        public static KeyValuePair<string, CachedSpectraFileData> CurrentlySelectedFile;
        public static ObservableCollection<AnnotatedSpecies> AllLoadedAnnotatedSpecies;

        /// <summary>
        /// Appends every processed species to <see cref="AllLoadedAnnotatedSpecies"/> as its own
        /// record, preserving input order and instance identity. No deduplication, grouping, or
        /// dictionary keying is applied: duplicate source rows remain separate records and species
        /// from different result files are never merged. Called by the loading pipeline
        /// (<c>DataLoading.LoadFile</c>) once per processed file.
        /// </summary>
        public static void AddLoadedSpecies(IEnumerable<AnnotatedSpecies> species)
        {
            lock (AllLoadedAnnotatedSpecies)
            {
                foreach (AnnotatedSpecies item in species)
                {
                    AllLoadedAnnotatedSpecies.Add(item);
                }
            }
        }
    }
}
