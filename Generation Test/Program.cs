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
            double mutationRate = 0.0000001; //0.0000001
            int attempts = 5;
            double increment = 0.000001; // 0.0002000
            int totalIncrements = 330;   // 0.0000005
            int individualLength = 25000; // average for humans is 3200000000
            int populationSize = 1; // Or 100?
            int roundingDigits = 6;
            bool countPercentPositive = true;

            bool absoluteLethal = true;
            bool showExpectedPercent = true;

            bool doSomaticMutation = false;

            double graphInterval = (mutationRate + increment*totalIncrements)/10;

            double mutationStdDev = 1;

            double germlineMutationMean = -mutationStdDev / 1.5; // -mutationStdDev / 1
            double somaticMutationMean = -mutationStdDev / 1.5; // -mutationStdDev / 1

            double probabilityOfPositiveGermlineMutation = (1 - Normal.CDF(germlineMutationMean, mutationStdDev, 0));
            Debug.WriteLine("Probability of positive mutation: " + probabilityOfPositiveGermlineMutation);

            Dictionary<double, double> averageGenerations = new Dictionary<double, double>();


            int generationCount;
            double[] fittestIndividual;
            double highestFitness;
            double germlineMutationRate;
            double somaticMutationRate;

            double beneficialMutationCount;
            bool countedBeneficialMutation;

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
                    beneficialMutationCount = 0;

                    fittestIndividual = new double[3];

                    fittestIndividual[0] = germlineMutationRate;
                    fittestIndividual[1] = somaticMutationRate;


                    fittestIndividual[2] = 0;

                    while (fittestIndividual[2] < targetFitness)
                    //for (int j=0; j<100; j++)
                    {
                        countedBeneficialMutation = false;
                        generationCount++;

                        highestFitness = double.NegativeInfinity;

                        double[] tempFittestIndividual = new double[3];
                        fittestIndividual.CopyTo(tempFittestIndividual, 0);

                        //Debug.WriteLine("----------");

                        for (int i = 0; i < populationSize; i++)
                        {
                            double[] currentIndividual = new double[3];
                            fittestIndividual.CopyTo(currentIndividual, 0);

                            // GERMLINE Fitness Mutate:
                            int rollSuccesses = probabilityOfGermlineSuccess.Sample();

                            Normal normal = new Normal(germlineMutationMean * rollSuccesses, mutationStdDev * Math.Sqrt(rollSuccesses));
                            double fitnessIncrease = normal.Sample();

                            //Debug.WriteLine("Mutations: " + rollSuccesses);
                            //Debug.WriteLine("Fitnes Increase: " + fitnessIncrease);

                            double oldFitness = currentIndividual[2];

                            
                            if (fitnessIncrease < 0)
                            {
                                if (absoluteLethal) continue;
                            }
                            else if (!countedBeneficialMutation && fitnessIncrease > 0)
                            {
                                //Debug.WriteLine(fitnessIncrease);
                                countedBeneficialMutation = true;
                                beneficialMutationCount++;
                            }


                            currentIndividual[2] += fitnessIncrease;

                            
                            if (doSomaticMutation)
                            {

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

                                Debug.WriteLine("Triggered Somatic");
                            }
                            else
                            {
                            
                            
                                if (currentIndividual[2] > highestFitness)
                                {
                                    highestFitness = currentIndividual[2];
                                    tempFittestIndividual = currentIndividual;
                                }
                            
                                
                            }


                        }

                        
                        if (!double.IsNegativeInfinity(highestFitness))
                        {
                            fittestIndividual = tempFittestIndividual;
                        }
                        
                        

                    }

                    if (countPercentPositive)
                    {
                        averageGenerations[mutationRate] += (double)beneficialMutationCount / (double)generationCount;
                        //Debug.WriteLine(beneficialMutationCount + "; "+ generationCount);
                        //Debug.WriteLine((double)beneficialMutationCount / (double)generationCount);
                    }
                    else
                        averageGenerations[mutationRate] += generationCount;
                }

                Debug.WriteLine("");
                Debug.WriteLine(z);
                Debug.WriteLine(stopwatch.ElapsedMilliseconds + "ms");

                averageGenerations[mutationRate] /= attempts;
                //Debug.WriteLine(averageGenerations[mutationRate]);

                mutationRate += increment;
            }

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
                Series expected = chart.Series.Add("Expected");

                expected.ChartType = SeriesChartType.FastLine;
                expected.Color = Color.CadetBlue;
                expected.BorderWidth = 3;

                generations.ChartType = SeriesChartType.FastLine;

                chart.BackColor = Color.White;
                CA.BackColor = Color.White;

                if (countPercentPositive)
                    CA.AxisY.Title = "Percentage Positive Mutations";
                else
                    CA.AxisY.Title = "Average Generations";

                CA.AxisX.Title = "Mutation Rates (mutations per gene)";

                CA.AxisX.TitleAlignment = StringAlignment.Center;
                CA.AxisY.TitleAlignment = StringAlignment.Center;

                CA.AxisX.TitleFont = new Font("Ariel", 15, FontStyle.Bold);
                CA.AxisY.TitleFont = new Font("Ariel", 15, FontStyle.Bold);

                CA.AxisX.Minimum = 0;
                CA.AxisX.Interval = graphInterval;

                CA.AxisY.Maximum = 1.2;

                if (countPercentPositive)
                {
                    generations.Color = Color.Red;
                    chart.Titles.Add($"Percentage Positve Mutations to Mutation Rate");
                }
                else
                {
                    generations.Color = Color.Blue;
                    chart.Titles.Add($"Average Generations to Mutation Rate (Target Fitness = {targetFitness})");
                }

                chart.Titles.ElementAt(0).Font = new Font("Ariel", 15, FontStyle.Bold);
                chart.Size = new Size(1920, 1080);
                generations.BorderWidth = 4;
                

                chart.AntiAliasing = AntiAliasingStyles.Graphics;
                chart.TextAntiAliasingQuality = TextAntiAliasingQuality.High;

                foreach (KeyValuePair<double, double> dataPoint in data)
                {
                    generations.Points.AddXY(Math.Round(dataPoint.Key, roundingDigits), dataPoint.Value);

                    if (showExpectedPercent)
                    {
                        expected.Points.AddXY(Math.Round(dataPoint.Key, roundingDigits), expectedEquation(dataPoint.Key));
                    }
                }

                string graphInfo = $"Generation Test {DateTime.Now.ToString("MMM-dd-yyyy hh-mm-ss-fff tt")}";
                string baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
                string projectDirectory = Directory.GetParent(Directory.GetParent(Directory.GetParent(baseDirectory).FullName).FullName).FullName;
                string graphsFolderPath = Path.Combine(projectDirectory, "Graphs");
                string imagePath = Path.Combine(graphsFolderPath, graphInfo + ".png");

                Debug.WriteLine("Image Path: " + imagePath);

                chart.SaveImage(imagePath, ChartImageFormat.Png);
            }

            double expectedEquation(double m)
            {
                double pRootc = populationSize*Math.Sqrt(probabilityOfPositiveGermlineMutation);
                double divisor = Math.E + (Math.Pow(probabilityOfPositiveGermlineMutation, 2) * Math.Pow(populationSize, 1d / Math.E));

                double y = Math.Pow(probabilityOfPositiveGermlineMutation, individualLength * m / divisor);
                y *= individualLength * m / divisor;
                y = 1d - y;
                y = Math.Pow(y, pRootc);
                y = 1d - y;

                return y;

                // 1d - Math.Pow(1d - (m * Math.Pow(probabilityOfPositiveGermlineMutation, m * individualLength/ Math.E) / Math.E), populationSize);
            }
        }
    }
}
