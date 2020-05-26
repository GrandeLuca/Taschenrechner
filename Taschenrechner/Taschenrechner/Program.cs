using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;

namespace Taschenrechner
{
    class Program
    {
        public static void WriteResults(double[] results,string[] names) 
        {
            for (int i = 0; i < results.Length; i++)
            {
                Console.WriteLine("{0} = {1}", names[i], results[i]);
            }
        }


        static void Main(string[] args)
        {
            //Konstanten

            double e = 1.602176634 * Math.Pow(10,-19) ;//C
            double m_e = 9.1093837015 * Math.Pow(10, -31);//kg
            double h = 6.62607015 * Math.Pow(10, -34); //Js
            double h_bar = h / (2 *Math.PI);
            double c = 299792458; //m/s
            double k_b=1.380649 * Math.Pow(10, -23);
            double micro = Math.Pow(10, -6);
            double nano = Math.Pow(10, -9);
            double pico = Math.Pow(10, -12);
            double sigma = 2 * Math.Pow(Math.PI, 5) * Math.Pow(k_b, 4) / (15 * Math.Pow(h, 3) * c * c);

            
            


            //Aufgabe 3.2
            if (false)
            {
                double lambda_1 = 589 * Math.Pow(10, -9);
                double lambda_2 = 253.7 * Math.Pow(10, -9);
                double E_1 = 0.36 * e;
                double E_2 = 3.14 * e;
                double h_exp = (E_2 - E_1) / (c * (1 / lambda_2 - 1 / lambda_1));
                double W_A_1 = h_exp * c / lambda_1 - E_1;
                double W_A_2 = h_exp * c / lambda_2 - E_2;
                double lambda_max = h_exp * c / W_A_1;

                Console.WriteLine("h={0}", h_exp);
                Console.WriteLine("W_A1={0}", W_A_1);
                Console.WriteLine("W_A2={0}", W_A_2);
                Console.WriteLine("lambda_max={0}", lambda_max);
                
            }

            //Aufgabe 3.3
            if (false)
            {
                double f = 1000000;
                double lambda = 556 * nano;
                double E_gamma = h * f;
                double W = 1000;
                double eta=0.16 * nano;
                double A_Auge=0.0001;
                double E_gamma2 = h * c / lambda;
                double P_Auge = eta * A_Auge;
                double N_t = P_Auge / E_gamma2;

                var results = new double[] {P_Auge,N_t};
                var resultnames = new string[] { "P_Auge","N/t" };
                WriteResults(results, resultnames);
            }

            //Aufgabe 4.1
            if (false)
            {
                double cw = 4.184*1000;
                double d = 0.1;
                double h1 = 0.25;
                double A = Math.PI * d * h1;
                double T_0 = 350;
                double T_t = 310;
                double V = Math.PI * d / 2 * d / 2 * h1;
                double m = V * 1000;
                double n = m / 0.0018;
                double c_n = cw * m;
                double T_w = 300;
                double alpha = 1 / (4 * Math.Pow(T_w, 3)) * Math.Log((T_w + T_t) / (T_t - T_w)) + 1 / (2 * Math.Pow(T_w, 3)) * Math.Atan(T_t / T_w) - (1 / (4 * Math.Pow(T_w, 3)) * Math.Log((T_w + T_0) / (T_0 - T_w)) + 1 / (2 * Math.Pow(T_w, 3)) * Math.Atan(T_0 / T_w));
                double t = c_n*alpha/ (sigma * A);
                Console.Write(c_n/3/sigma/A*(Math.Pow(1/T_t,3)- Math.Pow(1 / T_0,3)));
            }

            //Aufgabe 4.2
            if (true)
            {
                double T = 300;
                double ny;
                double alpha;
                double lambda;
                var results = new List<double>();
                double nystep = Math.Pow(10, 8);
                double lambdastep = Math.Pow(10, -12);

                 //Sucht ny_max
                for (int i = 0; i < Math.Pow(10,8); i++)
                {
                    ny = i * nystep;
                    alpha = h * ny / k_b / T;
                    double res = 8 * Math.PI*h*ny*ny*ny/c/c/c/(Math.Exp(alpha)-1) ;
                    results.Add(res);
                }
                double v_max = results.IndexOf(results.Max()) * nystep;
                Console.WriteLine("ny_max = {0}",v_max );
                results.Clear();

                //Sucht lambda_max
                for (int i = 0; i < Math.Pow(10, 8); i++)
                {
                    lambda = i * lambdastep;
                    alpha = h * c/lambda / k_b / T;
                    double res = 8 * Math.PI * h*Math.Pow(c,1)/Math.Pow(lambda,5)/  (Math.Exp(alpha) - 1);
                    results.Add(res);
                }
                double l_max = results.IndexOf(results.Max()) * lambdastep;
                Console.WriteLine("lamdba_max = {0}", l_max );
                double c_max = l_max * v_max;
                Console.WriteLine("c_max = {0}", c_max);
                Console.WriteLine("c_max/c = {0}", c_max / c);
            }
        }
    }
}
