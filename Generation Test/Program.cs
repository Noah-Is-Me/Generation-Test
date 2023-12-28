using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics.Distributions;

namespace Generation_Test
{
    internal static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            //Application.Run(new Form1());
            Random random = new Random();
            Stopwatch stopwatch = new Stopwatch();

            int targetFitness = 100;
            double mutationRate = 0.0000001;
            int attempts = 30;
            double increment = 0.0000001; // 0.0002000
            int totalIncrements = 1500;   // 0.0000005
            int individualLength = 25000; // average for humans is 3200000000
            int populationSize = 100; // Or 100?
            int roundingDigits = 6;

            double graphInterval = (mutationRate + increment*totalIncrements)/10;

            double mutationStdDev = 1;

            double germlineMutationMean = -mutationStdDev / 1;
            double somaticMutationMean = -mutationStdDev / 1;

            double probabilityOfPositiveGermlineMutation = (1 - Normal.CDF(germlineMutationMean, mutationStdDev, 0));
            Debug.WriteLine("Probability of positive mutation: " + probabilityOfPositiveGermlineMutation);

            Dictionary<double, double> averageGenerations = new Dictionary<double, double>();


            int generationCount;
            double[] fittestIndividual;
            double highestFitness;
            double germlineMutationRate;
            double somaticMutationRate;

            for (int z = 0; z < totalIncrements; z++)
            {
                averageGenerations.Add(mutationRate, 0);
                stopwatch.Restart();

                germlineMutationRate = mutationRate;
                somaticMutationRate = mutationRate;

                Poisson probabilityOfGermlineSuccess = new Poisson(germlineMutationRate * individualLength);
                Poisson probabilityOfSomaticSuccess = new Poisson(somaticMutationRate * individualLength);

                for (int n = 0; n < attempts; n++)
                {
                    generationCount = 0;

                    fittestIndividual = new double[3];

                    fittestIndividual[0] = germlineMutationRate;
                    fittestIndividual[1] = somaticMutationRate;


                    fittestIndividual[2] = 0;

                    while (fittestIndividual[2] < targetFitness)
                    {
                        generationCount++;

                        highestFitness = double.NegativeInfinity;

                        double[] tempFittestIndividual = new double[3];
                        fittestIndividual.CopyTo(tempFittestIndividual, 0);

                        for (int i = 0; i < populationSize; i++)
                        {
                            double[] currentIndividual = new double[3];
                            fittestIndividual.CopyTo(currentIndividual, 0);

                            // GERMLINE Fitness Mutate:
                            int rollSuccesses = probabilityOfGermlineSuccess.Sample();
                            double fitnessIncrease = normalDistribution(germlineMutationMean * rollSuccesses, mutationStdDev * Math.Sqrt(rollSuccesses));
                            currentIndividual[2] += fitnessIncrease;

                            // SOMATIC Fitness Mutate:
                            double[] currentSomaticIndividual = new double[3];
                            currentIndividual.CopyTo(currentSomaticIndividual, 0);

                            rollSuccesses = probabilityOfSomaticSuccess.Sample();
                            fitnessIncrease = normalDistribution(somaticMutationMean * rollSuccesses, mutationStdDev * Math.Sqrt(rollSuccesses));
                            currentSomaticIndividual[2] += fitnessIncrease;

                            if (currentSomaticIndividual[2] > highestFitness)
                            {
                                highestFitness = currentSomaticIndividual[2];
                                tempFittestIndividual = currentIndividual;
                            }


                        }

                        if (!double.IsNegativeInfinity(highestFitness))
                        {
                            fittestIndividual = tempFittestIndividual;
                        }

                    }
                    averageGenerations[mutationRate] += generationCount;
                }

                Debug.WriteLine("");
                Debug.WriteLine(z);
                Debug.WriteLine(stopwatch.ElapsedMilliseconds + "ms");

                averageGenerations[mutationRate] /= attempts;

                mutationRate += increment;
            }

            /*
            Debug.WriteLine("");
            foreach (KeyValuePair<double, double> dataPair in averageGenerations)
            {
                Debug.WriteLine(dataPair.Key + "; " + dataPair.Value);
            }
            */
            CreateGraph(averageGenerations);

            double normalDistribution(double mean, double stdDev)
            {
                double u1 = 1.0 - random.NextDouble(); //uniform(0,1] random doubles
                double u2 = 1.0 - random.NextDouble();
                double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
                double randNormal = mean + stdDev * randStdNormal;
                return randNormal;
            }

            void CreateGraph(Dictionary<double, double> data)
            {
                Chart chart = new Chart();

                ChartArea CA = chart.ChartAreas.Add("A1");
                Series generations = chart.Series.Add("Average Generations");
                generations.ChartType = SeriesChartType.FastLine;

                chart.BackColor = Color.White;
                CA.BackColor = Color.White;

                CA.AxisY.Title = "Average Generations";
                CA.AxisX.Title = "Mutation Rates (mutations per gene)";

                CA.AxisX.TitleAlignment = StringAlignment.Center;
                CA.AxisY.TitleAlignment = StringAlignment.Center;

                CA.AxisX.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
                CA.AxisY.TitleFont = new Font("Ariel", 15, FontStyle.Bold);

                CA.AxisX.Minimum = 0;
                CA.AxisX.Interval = graphInterval;

                chart.Titles.Add($"Average Generations to Mutation Rate (Target Fitness = {targetFitness})");
                chart.Titles.ElementAt(0).Font = new Font("Ariel", 15, FontStyle.Bold);
                chart.Size = new Size(1920, 1080);
                chart.Series["Average Generations"].BorderWidth = 4;
                chart.Series["Average Generations"].Color = Color.Blue;

                chart.AntiAliasing = AntiAliasingStyles.Graphics;
                chart.TextAntiAliasingQuality = TextAntiAliasingQuality.High;

                foreach (KeyValuePair<double, double> dataPoint in data)
                {
                    chart.Series["Average Generations"].Points.AddXY(Math.Round(dataPoint.Key, roundingDigits), dataPoint.Value*100);
                }

                string graphInfo = $"Generation Test {DateTime.Now.ToString("MMM-dd-yyyy hh-mm-ss-fff tt")}";
                string baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
                string projectDirectory = Directory.GetParent(Directory.GetParent(Directory.GetParent(baseDirectory).FullName).FullName).FullName;
                string graphsFolderPath = Path.Combine(projectDirectory, "Graphs");
                string imagePath = Path.Combine(graphsFolderPath, graphInfo + ".png");

                Debug.WriteLine("Image Path: " + imagePath);

                chart.SaveImage(imagePath, ChartImageFormat.Png);
            }
        }
    }
}
