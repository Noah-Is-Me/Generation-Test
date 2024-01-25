﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
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
            int totalIncrements = 660;   // 0.0000005
            int individualLength = 25000; // average for humans is 3200000000
            int populationSize = 1; // Or 100?
            int roundingDigits = 6;
            bool countPercentPositive = true;

            bool absoluteLethal = false;
            bool showExpectedPercent = false;
            bool getRawData = true;

            bool doSomaticMutation = false;

            double graphInterval = (mutationRate + increment*totalIncrements)/10;

            double mutationStdDev = 1;

            double germlineMutationMean = -mutationStdDev / 1; // -mutationStdDev / 1
            double somaticMutationMean = -mutationStdDev / 1; // -mutationStdDev / 1

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

                    //while (fittestIndividual[2] < targetFitness)
                    for (int j=0; j<100; j++)
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
                            //double rollSuccesses = mutationRate * individualLength;
                            
                            //double probOfOneNegative = 1d - Math.Pow(probabilityOfPositiveGermlineMutation, rollSuccesses);
                            //double probOfOnePositive = 1d - Math.Pow(1d-probabilityOfPositiveGermlineMutation, rollSuccesses);

                            Normal normal = new Normal(germlineMutationMean * rollSuccesses, mutationStdDev * Math.Sqrt(rollSuccesses));
                            //double fitnessIncrease = normal.Sample();
                            double fitnessIncrease = 1-(normal.CumulativeDistribution(0));


                            // SOMETHING IS WRONG WITH LINE 110 THE POISSON DISTRIBUTION
                            // It should match m*g but it does not when comparing cdf of normal distribution


                            fitnessIncrease = rollSuccesses;

                            /*
                            if (absoluteLethal && random.NextDouble() < probOfOneNegative)
                            {
                                continue;
                            } 
                            
                            if (!countedBeneficialMutation && fitnessIncrease >= 0)
                            {
                                //Debug.WriteLine(fitnessIncrease);
                                countedBeneficialMutation = true;
                                beneficialMutationCount++;
                            }
                            */

                            /*
                            if (!countedBeneficialMutation && random.NextDouble() < (1 - normal.CumulativeDistribution(0)))
                            {
                                countedBeneficialMutation = true;
                                beneficialMutationCount++;
                            }
                            */

                            if (!countedBeneficialMutation)
                            {
                                countedBeneficialMutation = true;
                                beneficialMutationCount+= fitnessIncrease;
                            }

                            //Debug.WriteLine("Mutations: " + rollSuccesses);
                            //Debug.WriteLine("Fitnes Increase: " + fitnessIncrease);

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

                //CA.AxisY.Maximum = 1.2;

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

                string csv = "";
                csv += $"{CA.AxisX.Title},{CA.AxisY.Title}";

                foreach (KeyValuePair<double, double> dataPoint in data)
                {
                    double x = Math.Round(dataPoint.Key, roundingDigits);

                    generations.Points.AddXY(x, dataPoint.Value);

                    if (showExpectedPercent)
                    {
                        expected.Points.AddXY(x, expectedEquation(dataPoint.Key));
                    }

                    if (getRawData)
                    {
                        csv += $"\n{x},{dataPoint.Value}";
                    }
                }

                string graphInfo = $"Generation Test {DateTime.Now.ToString("MMM-dd-yyyy hh-mm-ss-fff tt")}";
                string baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
                string projectDirectory = Directory.GetParent(Directory.GetParent(Directory.GetParent(baseDirectory).FullName).FullName).FullName;
                string graphsFolderPath = Path.Combine(projectDirectory, "Graphs");
                string imagePath = Path.Combine(graphsFolderPath, graphInfo + ".png");

                string csvPath = Path.Combine(graphsFolderPath, graphInfo + ".txt");


                if (getRawData)
                {
                    //Debug.WriteLine(csv);

                    File.WriteAllText(csvPath, csv);
                }

                Debug.WriteLine("Image Path: " + imagePath);

                chart.SaveImage(imagePath, ChartImageFormat.Png);
            }

            double expectedEquation(double m)
            {
                /*
                double pRootc = Math.Sqrt(probabilityOfPositiveGermlineMutation * populationSize);
                double divisor = 1; // Math.E + (Math.Pow(probabilityOfPositiveGermlineMutation, 2) * Math.Pow(populationSize, 1d / Math.E));

                double y = Math.Pow(probabilityOfPositiveGermlineMutation, individualLength * m / divisor);
                y *= individualLength * m / divisor;
                y = 1d - y;
                y = Math.Pow(y, pRootc);
                y = 1d - y;
                */

                /*
                double t = 1d;

                // PART 1!!!
                double y1 = 1d - probabilityOfPositiveGermlineMutation;
                y1 = Math.Pow(y1, individualLength * m * t);
                y1 = Math.Pow(y1, populationSize);
                y1 = 1d - y1;
                

                // PART 2!!!
                double y2 = t * m * (1d - probabilityOfPositiveGermlineMutation);
                y2 = 1 - y2;
                y2 = Math.Pow(y2, individualLength);
                y2 = 1 - y2;
                y2 = Math.Pow(y2, populationSize);
                y2 = 1 - y2;

                // COMBINE THEM!!!
                double y = y1 * y2;

                */

                // NEW MIGHT WORK

                double c = Math.Pow(probabilityOfPositiveGermlineMutation, 2);

                double y1 = 1d - m;
                y1 += c * m;
                y1 = Math.Pow(y1, individualLength);

                double y2 = 1d - (c * m);
                y2 = Math.Pow(y2, individualLength);
                y2 = 1d - y2;

                double y = y1 * y2;
                y = 1d - y;
                y = Math.Pow(y, populationSize);
                y = 1d - y;


                return y;

                // 1d - Math.Pow(1d - (m * Math.Pow(probabilityOfPositiveGermlineMutation, m * individualLength/ Math.E) / Math.E), populationSize);
            }
        }
    }
}
