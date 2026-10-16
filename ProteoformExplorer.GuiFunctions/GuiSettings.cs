using Nett;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;

namespace ProteoformExplorer.GuiFunctions
{
    public static class GuiSettings
    {
        static GuiSettings()
        {
            NameMappings = [];
            if (File.Exists(SettingsFilePath))
                LoadSettingsFromFile();
            else
                SaveSettingsToFile();
        }

        // from: http://seaborn.pydata.org/tutorial/color_palettes.html (qualitative bright palette)
        public static string[] ChartColorPalette = new string[]
        {
            @"#00ca3d",  // green
            @"#e91405",  // red
            @"#9214e1",  // purple
            @"#9f4a00",  // brown
            @"#f44ac1",  // pink
            @"#fdc807",  // gold
            @"#17d4ff",  // teal
            @"#3729fe",  // blue
            @"#ff8000",  // orange
        };

        public static Color UnannotatedSpectrumColor = Color.Gray;
        public static Color TicColor = Color.FromArgb(53, 59, 72); // dark gray
        public static Color DeconvolutedColor = Color.FromArgb(0, 255, 255); // blue
        public static Color IdentifiedColor = Color.FromArgb(102, 0, 204); // purple
        public static Color RtIndicatorColor = Color.Red;
        public static int IntegrationFillAlpha = 50;
        public static int FillAlpha = 190;
        public static double ChartTickFontSize = 10;
        public static double ChartHeaderFontSize = 18;
        public static double ChartLegendFontSize = 14;
        public static double ChartAxisLabelFontSize = 13;
        public static double ChartLineWidth = 1;
        public static double AnnotatedEnvelopeLineWidth = 2;
        public static ScottPlot.Alignment LegendLocation = ScottPlot.Alignment.UpperLeft;
        public static bool ShowChartGrid = false;
        public static double XLabelRotation = 30;
        public static bool FillWaterfall = true;
        public static bool FillSideView = false;
        public static bool WaterfallXics = false;
        public static double RtExtractionWindow = 6.0; // in minutes
        public static double WaterfallSpacing = 10; // in millimeters
        public static double MassHistogramBinWidth = 500; // in daltons
        public static bool DpiScaling = true;
        public static double DpiScalingX = 1;
        public static double DpiScalingY = 1;
        public static int TicRollingAverage = 3;

        #region Coloring identified Tics

        private static Dictionary<string, Color> UsedColorsWithKey { get; } = [];
        private static Queue<Color> ColorQueue = new Queue<Color>(new[]
        {
            Color.Purple, 
            Color.Green,
            Color.OrangeRed,
            Color.Yellow, 
            Color.Pink,
            Color.DarkMagenta,
            Color.DarkCyan,
            Color.DarkGoldenrod,
            Color.DarkOliveGreen,
            Color.DarkRed,
            Color.DarkSlateBlue,
            Color.DarkTurquoise,
            Color.DarkViolet,
            Color.DeepPink,
            Color.DeepSkyBlue,
        });
        public static List<NameMapping> NameMappings;

        public static string ConvertName(this string input)
        {
            if (string.IsNullOrWhiteSpace(input))
                return input;

            // Helper to normalize names: lower case, remove extension if present
            static string Normalize(string name)
            {
                if (string.IsNullOrEmpty(name))
                    return string.Empty;
                var lower = name.ToLowerInvariant();
                var idx = lower.LastIndexOf('.');
                return idx > 0 ? lower[..idx] : lower;
            }

            var normalizedInput = Normalize(input);

            foreach (var mapping in NameMappings)
            {
                // Check normalized long names
                foreach (var longName in mapping.LongNames)
                {
                    if (Normalize(longName) == normalizedInput)
                        return mapping.ShortName;
                }
                // Also check if input is a normalized version of the short name
                if (Normalize(mapping.ShortName) == normalizedInput)
                    return mapping.ShortName;
            }

            return input;
        }

        public static Color ConvertStringToColor(this string input)
        {
            // Saved to file as short name
            var mapping = NameMappings.FirstOrDefault(m => m.ShortName == input);
            if (mapping is not null) 
                return mapping.Color;

            // Saved to file as long name
            mapping = NameMappings.FirstOrDefault(m => m.LongNames.Contains(input));
            if (mapping is not null) 
                return mapping.Color;

            // Saved during this instance
            if (UsedColorsWithKey.TryGetValue(input, out var color)) return color;

            if (ColorQueue.Count == 0)
            {
                ColorQueue = new Queue<Color>(new[]
                {
                    Color.Purple,
                    Color.Green,
                    Color.OrangeRed,
                    Color.Yellow,
                    Color.Pink,
                    Color.DarkMagenta,
                    Color.DarkCyan,
                    Color.DarkGoldenrod,
                    Color.DarkOliveGreen,
                    Color.DarkRed,
                    Color.DarkSlateBlue,
                    Color.DarkTurquoise,
                    Color.DarkViolet,
                    Color.DeepPink,
                    Color.DeepSkyBlue,
                });
            }

            var newColor = ColorQueue.Dequeue();
            UsedColorsWithKey[input] = newColor;
            return newColor;
        }

        #endregion

        #region Settings Loading and Writing

        private static readonly string SettingsDirectory = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.ApplicationData),
            "ProteoformExplorer"
        );

        private static readonly string SettingsFilePath = Path.Combine(SettingsDirectory, "ProteoformExplorerSettings.toml");

        public static void LoadSettingsFromFile()
        {
            if (!File.Exists(SettingsFilePath)) 
                return;

            var table = Toml.ReadFile(SettingsFilePath).ToDictionary();
            FromTomlTable(table);
        }

        public static void SaveSettingsToFile()
        {
            if (!Directory.Exists(SettingsDirectory))
                Directory.CreateDirectory(SettingsDirectory);

            var table = ToTomlTable();
            File.WriteAllText(SettingsFilePath, table.ToString());
        }

        private static TomlTable ToTomlTable()
        {
            TomlTable table = Toml.Create();

            Type type = typeof(GuiSettings);

            var fields = type.GetFields();

            foreach (var field in fields)
            {
                var name = field.Name;

                if (name.ToLower() == "dpiscalingx" || name.ToLower() == "dpiscalingy")
                {
                    continue;
                }

                var v = field.GetValue(null);

                if (v is int val)
                {
                    table.Add(name, val);
                }
                else if (v is double dval)
                {
                    table.Add(name, dval);
                }
                else if (v is bool b)
                {
                    table.Add(name, b);
                }
                else if (v is ScottPlot.Alignment alignment)
                {
                    table.Add(name, alignment);
                }
                else if (v is Color color)
                {
                    table.Add(name, color);
                }
                else if (v is string[] strarr)
                {
                    table.Add(name, strarr);
                }
                else if (v is Dictionary<string, string> dict)
                {
                    table.Add(name, dict);
                }
                else if (v is Dictionary<string, Color> colorDict)
                {
                    table.Add(name, colorDict);
                }
                else if (v is Dictionary<string, object> objDict)
                {
                    table.Add(name, objDict);
                }
                else if (v is object[] objarr)
                {
                    table.Add(name, objarr);
                }
                else if (v is List<NameMapping> maps)
                {
                    table.Add(name, maps);
                }
            }

            return table;
        }

        private static void FromTomlTable(Dictionary<string, object> tomlDictionary)
        {
            Type type = typeof(GuiSettings);
            var fields = type.GetFields();

            foreach (var field in fields)
            {
                var name = field.Name;

                if (tomlDictionary.TryGetValue(name, out var v))
                {
                    if (v is int val)
                    {
                        field.SetValue(null, val);
                    }
                    else if (v is double dval)
                    {
                        field.SetValue(null, dval);
                    }
                    else if (v is bool b)
                    {
                        field.SetValue(null, b);
                    }
                    else if (v is ScottPlot.Alignment alignment)
                    {
                        field.SetValue(null, alignment);
                    }
                    else if (v is Color color)
                    {
                        field.SetValue(null, color);
                    }
                    else if (v is object[] objarr)
                    {
                        if (objarr.Length == 0)
                            continue;
                            
                        if (objarr[0] is Dictionary<string, object> { Count: 3 })
                        {
                            foreach (var obj in objarr)
                            {
                                var dict = (Dictionary<string, object>)obj;
                                var colorDict = (Dictionary<string, object>)dict["Color"];
                                Color col = Color.Black;

                                string shortName = dict["ShortName"] as string;
                                var longNames = (dict["LongNames"] as object[]).Select(p => p.ToString()).ToList();

                                if (colorDict is not null)
                                {
                                    int a = int.Parse(colorDict["A"].ToString());
                                    int r = int.Parse(colorDict["R"].ToString());
                                    int g = int.Parse(colorDict["G"].ToString());
                                    int b2 = int.Parse(colorDict["B"].ToString());
                                    col = Color.FromArgb(a, r, g, b2);
                                }

                                var mapping = new NameMapping()
                                {
                                    ShortName = shortName,
                                    LongNames = longNames,
                                    Color = col
                                };
                                NameMappings.Add(mapping);
                            }
                        }
                        else if (objarr[1] is string str)
                        {
                            field.SetValue(null, objarr.Select(p => (string)p).ToArray());
                        }
                    }
                    else if (v is Dictionary<string, string> dict)
                    {
                        field.SetValue(null, dict);
                    }
                    else if (v is Dictionary<string, Color> colorDict)
                    {
                        field.SetValue(null, colorDict);
                    }
                    else if (v is List<NameMapping> maps)
                    {
                        field.SetValue(null, maps);
                    }
                }
            }
        }

        #endregion
    }
}
